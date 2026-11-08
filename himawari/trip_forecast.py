#!/usr/bin/env python3
"""
Sunshine forecast for the car along the race route.

Input: start time and speed. The script works out where the car is at every moment
(driving only in the daily driving window, default 08:00–17:00 local time, standing still
overnight), fetches Open-Meteo irradiance forecasts for the route points and gives the
global horizontal irradiance (GHI, W/m²) the car gets on its way, plus the energy per day.

If the trip lies in the past (e.g. a replay of the race week 2025), the archived forecasts
(Historical Forecast API) are used and the actual satellite-derived values are added, so the
forecast can be checked against what really happened.

Times are local race time (ACST = UTC+9:30, valid for NT and SA in August) unless --utc.
At most the next 5 days are calculated (--max-days): further ahead the forecasts are too
unreliable. If the finish is not reached by then, the trip ends where the car is.

Output (folder trip_forecast/):
  trip_forecast.csv    every 5 minutes: time, km, position, driving/stopped, GHI per model
  trip_forecast.png    sunshine along the trip and position over time
  summary              printed: arrival, km and energy per day

Usage:
    pip install numpy requests matplotlib
    python trip_forecast.py --start "2027-08-22 08:30" --speed 80
    python trip_forecast.py --start now --speed 85 --start-km 850       # car is at km 850 now
    python trip_forecast.py --start "2025-08-24 08:00" --speed 80       # replay with actual values
    python trip_forecast.py --start now --speed 80 --models best_match,ecmwf_ifs025,gfs_global
    python trip_forecast.py --start now --speed 80 --panel-area 4 --efficiency 0.24   # + kWh for the car
    python trip_forecast.py --start now --speed 80 --no-driving-window  # drive non-stop
"""
import argparse
import ast
import concurrent.futures as cf
import csv
import datetime as dt
import json
import os
import sys
import threading

import numpy as np

FC_URL = "https://api.open-meteo.com/v1/forecast"
HIST_URL = "https://historical-forecast-api.open-meteo.com/v1/forecast"
SAT_URL = "https://satellite-api.open-meteo.com/v1/archive"
LOCAL = dt.timedelta(hours=9, minutes=30)        # ACST
STEP = dt.timedelta(minutes=5)
BATCH = 50
TIMEOUT = 30
UA = "himawari-nowcast-trip/1.0 (research)"


# ---------------------------------------------------------------- route and sun
def load_route(path):
    """Route points [(lat, lon, name)] from a .py file with a list of dicts (read, not run),
    a .json list or a .csv with latitude/longitude columns."""
    ext = os.path.splitext(path)[1].lower()
    data = None
    if ext == ".py":
        for node in ast.parse(open(path, encoding="utf-8").read()).body:
            if isinstance(node, ast.Assign):
                try:
                    val = ast.literal_eval(node.value)
                except ValueError:
                    continue
                if isinstance(val, list):
                    data = val
                    break
    elif ext == ".json":
        data = json.load(open(path, encoding="utf-8"))
        data = data.get("locations", []) if isinstance(data, dict) else data
    else:
        data = list(csv.DictReader(open(path, encoding="utf-8")))
    if not data or len(data) < 2:
        sys.exit(f"No route found in {path}")
    return [(float(d["latitude"]), float(d["longitude"]), str(d.get("name", f"P{i + 1}")))
            for i, d in enumerate(data)]


def haversine_km(lat1, lon1, lat2, lon2):
    p1, p2 = np.radians(lat1), np.radians(lat2)
    a = (np.sin((p2 - p1) / 2) ** 2
         + np.cos(p1) * np.cos(p2) * np.sin(np.radians(lon2 - lon1) / 2) ** 2)
    return 6371.0 * 2 * np.arcsin(np.sqrt(a))


def cos_zenith(lat, lon, t):
    la = np.radians(np.asarray(lat, float))
    hour = t.hour + t.minute / 60 + t.second / 3600
    g = 2 * np.pi / 365 * (t.timetuple().tm_yday - 1 + (hour - 12) / 24)
    eqt = 229.18 * (0.000075 + 0.001868 * np.cos(g) - 0.032077 * np.sin(g)
                    - 0.014615 * np.cos(2 * g) - 0.040849 * np.sin(2 * g))
    decl = (0.006918 - 0.399912 * np.cos(g) + 0.070257 * np.sin(g) - 0.006758 * np.cos(2 * g)
            + 0.000907 * np.sin(2 * g) - 0.002697 * np.cos(3 * g) + 0.00148 * np.sin(3 * g))
    ha = np.radians((hour * 60 + eqt + 4 * np.asarray(lon, float)) / 4 - 180)
    return np.sin(la) * np.sin(decl) + np.cos(la) * np.cos(decl) * np.cos(ha)


def clear_sky(lat, lon, t):
    """Clear-sky GHI in W/m² (Haurwitz), 0 at night."""
    cz = cos_zenith(lat, lon, t)
    czp = np.maximum(cz, 1e-3)
    return np.where(cz > 0, 1098.0 * czp * np.exp(-0.057 / czp), 0.0)


# ---------------------------------------------------------------- trip
def simulate(start_utc, speed, start_km, end_km, window, until_utc, step=STEP):
    """Car position every `step` until the finish or until_utc: [(utc, km, driving)]. Outside
    the daily driving window (local times) the car stands still; window=None: non-stop."""
    out, t, km = [], start_utc, start_km
    while km < end_km and t < until_utc:
        driving = window is None or window[0] <= (t + LOCAL).time() < window[1]
        out.append((t, km, driving))
        if driving:
            km = min(end_km, km + speed * step.total_seconds() / 3600)
        t += step
    out.append((t, km, False))
    return out


# ---------------------------------------------------------------- Open-Meteo
DOWN = set()


def fetch(url, pts, params, label, quiet=False):
    """All points (BATCH per request, in parallel); per point if lists are refused.
    Returns a list with the JSON answer per point (None if missing) or None."""
    import requests
    host = url.split("/")[2]
    if host in DOWN:
        return None
    local = threading.local()
    err = {"msg": ""}

    def get(idx):
        if not hasattr(local, "s"):
            local.s = requests.Session()
            local.s.headers["User-Agent"] = UA
        p = dict(params, latitude=",".join(f"{pts[i][0]:.4f}" for i in idx),
                 longitude=",".join(f"{pts[i][1]:.4f}" for i in idx))
        try:
            r = local.s.get(url, params=p, timeout=TIMEOUT)
        except requests.RequestException as e:
            err["msg"] = f"no answer from {host} ({type(e).__name__})"
            DOWN.add(host)
            return None
        if r.status_code != 200:
            err["msg"] = f"HTTP {r.status_code}: {r.text[:200]}"
            return None
        d = r.json()
        d = d if isinstance(d, list) else [d]
        return d if len(d) == len(idx) else None

    batches = [list(range(b, min(len(pts), b + BATCH))) for b in range(0, len(pts), BATCH)]
    with cf.ThreadPoolExecutor(4) as ex:
        res = list(ex.map(get, batches))
    if all(r is None for r in res) and host not in DOWN and len(batches[0]) > 1:
        if get(batches[0][:1]) is not None:
            with cf.ThreadPoolExecutor(4) as ex:
                single = list(ex.map(lambda i: get([i]), range(len(pts))))
            return [s[0] if s else None for s in single]
    if all(r is None for r in res):
        if not quiet or host in DOWN:
            print(f"  ! {label}: {err['msg'] or 'no data'}")
        return None
    out = [None] * len(pts)
    for idx, r in zip(batches, res):
        for i, d in zip(idx, r or []):
            out[i] = d
    return out


class Series:
    """GHI time series per route point, stored as clear-sky index between its time stamps;
    value(i, t) = clear sky at point i and time t × interpolated index. Smooth day curve,
    clouds from the forecast."""

    def __init__(self, answers, var, block_key, shift, pts):
        self.lat = np.array([p[0] for p in pts])
        self.lon = np.array([p[1] for p in pts])
        stamps = None
        vals = []
        for d in answers:
            blk = (d or {}).get(block_key) or {}
            if stamps is None and blk.get("time"):
                stamps = [dt.datetime.strptime(s[:16], "%Y-%m-%dT%H:%M") - shift for s in blk["time"]]
            vals.append(blk.get(var))
        if stamps is None:
            raise ValueError("empty")
        T = len(stamps)
        g = np.array([[np.nan if (v is None or v[j] is None) else v[j] for j in range(T)]
                      if v is not None and len(v) == T else [np.nan] * T for v in vals], float)
        self.sec = np.array([s.replace(tzinfo=dt.timezone.utc).timestamp() for s in stamps])
        self.k = np.full_like(g, np.nan)
        for j, s in enumerate(stamps):
            cs = clear_sky(self.lat, self.lon, s)
            ok = (cs > 30) & np.isfinite(g[:, j])
            self.k[ok, j] = np.clip(g[ok, j] / cs[ok], 0, 1.3)
        self.n_values = int(np.isfinite(g).sum())

    def value(self, i, t):
        s = t.replace(tzinfo=dt.timezone.utc).timestamp()
        j = int(np.searchsorted(self.sec, s))
        if j <= 0 or j >= len(self.sec):
            return np.nan
        a, b = self.k[i, j - 1], self.k[i, j]
        w = (s - self.sec[j - 1]) / (self.sec[j] - self.sec[j - 1])
        if np.isnan(a) and np.isnan(b):
            k = np.nan
        elif np.isnan(a):
            k = b
        elif np.isnan(b):
            k = a
        else:
            k = (1 - w) * a + w * b
        cs = float(clear_sky(self.lat[i], self.lon[i], t))
        if cs <= 0:
            return 0.0
        return cs * k if np.isfinite(k) else np.nan


def get_series(url, pts, base, label):
    """Try 15-minute values, then hourly. Open-Meteo labels a mean value with the end of its
    interval, so it is shifted to the middle."""
    for block, var, shift, extra in (
            ("minutely_15", "shortwave_radiation", dt.timedelta(minutes=7, seconds=30),
             {"minutely_15": "shortwave_radiation"}),
            ("hourly", "shortwave_radiation", dt.timedelta(minutes=30),
             {"hourly": "shortwave_radiation"})):
        # 15-minute values are tried quietly; hourly values are the fallback
        ans = fetch(url, pts, dict(base, **extra), f"{label} ({block})", quiet=block == "minutely_15")
        if ans is None:
            if url.split("/")[2] in DOWN:
                return None, None
            continue
        try:
            s = Series(ans, var, block, shift, pts)
        except ValueError:
            continue
        if s.n_values:
            return s, block
    return None, None


# ---------------------------------------------------------------- clouds
def cloud_flags(rows, name, thr):
    """Per row: True = cloud (less than thr of the clear-sky sunshine arrives), False =
    clear, None = no judgement possible (night, sun low: clear sky < 200 W/m², no forecast)."""
    out = []
    for r in rows:
        v, cs = r.get(f"v:{name}"), r["cs"]
        # with the sun very low the clear-sky estimate is too rough to judge clouds
        out.append(None if cs < 200 or v is None or not np.isfinite(v) else bool(v / cs < thr))
    return out


def cloud_share(rows, flags):
    """Share of the driven distance (in daylight) that lies under cloud, and that distance."""
    under = total = 0.0
    for i in range(len(rows) - 1):
        if rows[i]["driving"] and flags[i] is not None:
            d = rows[i + 1]["km"] - rows[i]["km"]
            total += d
            under += d if flags[i] else 0.0
    return (under / total if total > 0 else None), under, total


# ---------------------------------------------------------------- map
# simplified coastline near the route (for orientation only, not exact)
COAST = [(-38.06, 141.0), (-37.6, 140.1), (-36.9, 139.7), (-36.0, 139.6), (-35.6, 138.1),
         (-35.0, 138.5), (-34.6, 138.45), (-34.2, 138.15), (-34.6, 137.9), (-35.1, 137.75),
         (-35.25, 137.0), (-34.6, 137.45), (-33.9, 137.6), (-33.2, 137.9), (-32.5, 137.75),
         (-33.0, 137.55), (-33.9, 136.7), (-34.75, 135.85), (-34.0, 135.2), (-33.3, 134.3),
         (-32.1, 133.7), (-31.6, 131.0), (-31.7, 129.0), (-31.65, 127.5)]
COAST_N = [(-15.0, 127.5), (-14.1, 129.5), (-14.9, 129.6), (-13.5, 130.0), (-12.46, 130.84),
           (-12.2, 131.0), (-11.3, 132.2), (-11.9, 133.5), (-12.0, 135.0), (-12.3, 136.8),
           (-13.9, 136.0), (-14.9, 135.5), (-16.0, 136.5), (-16.5, 137.8), (-17.3, 139.3),
           (-17.6, 140.8), (-16.6, 141.4)]
TOWNS = [("Darwin", -12.46, 130.84), ("Katherine", -14.47, 132.26),
         ("Tennant Creek", -19.65, 134.19), ("Alice Springs", -23.70, 133.88),
         ("Coober Pedy", -29.01, 134.75), ("Port Augusta", -32.49, 137.77),
         ("Adelaide", -34.93, 138.60)]


def save_trip_map(rows, names, pts, route, km, a, end_km, path):
    """Route with the car's calculated positions: a dot every 30 min while driving, coloured by
    the forecast sunshine (first model); stars at the start of each driving day, squares at
    each overnight stop, labelled with time and km."""
    import matplotlib.pyplot as plt
    first = names[0] if names[0] != "actual (satellite)" else (names[1] if len(names) > 1 else names[0])
    fig, ax = plt.subplots(figsize=(8.5, 11))
    for line in (COAST, COAST_N):
        ax.plot([p[1] for p in line], [p[0] for p in line], color="0.6", lw=1)
    ax.plot([p[1] for p in route], [p[0] for p in route], color="0.75", lw=5, solid_capstyle="round",
            zorder=1, label="route")
    for nme, tla, tlo in TOWNS:
        ax.plot(tlo, tla, "s", color="0.3", ms=4, zorder=2)
        ax.annotate(nme, (tlo, tla), xytext=(6, -3), textcoords="offset points", fontsize=8,
                    color="0.3")
    # the driven route as a band: red = under cloud, green = clear (like route_map.png)
    flags = cloud_flags(rows, first, a.cloud_index)
    for i in range(len(rows) - 1):
        if rows[i]["driving"] and flags[i] is not None:
            ax.plot([rows[i]["lon"], rows[i + 1]["lon"]], [rows[i]["lat"], rows[i + 1]["lat"]],
                    color="#d62728" if flags[i] else "#2ca02c", lw=7, alpha=0.75,
                    solid_capstyle="butt", zorder=3)
    lines = []
    for d in sorted({(r["t"] + LOCAL).date() for r in rows}):
        idx = [i for i, r in enumerate(rows) if (r["t"] + LOCAL).date() == d]
        sub = [rows[i] for i in idx] + ([rows[idx[-1] + 1]] if idx[-1] + 1 < len(rows) else [])
        sh, under, tot = cloud_share(sub, [flags[i] for i in idx] + [None])
        if sh is not None:
            lines.append(f"{d:%a %d.%m.}: {100 * sh:3.0f} %  ({under:.0f} of {tot:.0f} km)")
    sh, under, tot = cloud_share(rows, flags)
    if sh is not None:
        txt = (f"Cloud on the driven route ({first})\n"
               f"total: {100 * sh:.0f} %  ({under:.0f} of {tot:.0f} km)\n" + "\n".join(lines)
               + f"\n(cloud = below {100 * a.cloud_index:.0f} % of clear-sky sunshine)")
        ax.text(0.98, 0.98, txt, transform=ax.transAxes, ha="right", va="top", fontsize=8,
                multialignment="left",
                family="monospace", bbox=dict(boxstyle="round", fc="white", ec="0.6", alpha=0.9),
                zorder=7)
    # dots every 30 min of the trip (counted from the start, whatever minute it starts at)
    every = max(1, int(round(1800 / STEP.total_seconds())))
    drv = [r for i, r in enumerate(rows) if r["driving"] and i % every == 0]
    v = np.array([r[f"v:{first}"] for r in drv], float)
    nan = ~np.isfinite(v)
    if nan.any():   # no forecast there (e.g. beyond the forecast range): grey, not invisible
        ax.scatter([r["lon"] for r, m in zip(drv, nan) if m], [r["lat"] for r, m in zip(drv, nan) if m],
                   s=28, facecolors="0.85", edgecolors="0.4", linewidths=0.5, zorder=4,
                   label="no forecast")
    sc = ax.scatter([r["lon"] for r, m in zip(drv, nan) if not m],
                    [r["lat"] for r, m in zip(drv, nan) if not m], c=v[~nan], cmap="inferno",
                    vmin=0, vmax=1100, s=28, edgecolors="k", linewidths=0.3, zorder=4)
    fig.colorbar(sc, ax=ax, shrink=0.5, pad=0.02, label=f"forecast sunshine ({first}) [W/m²]")
    # day starts and overnight stops
    events = []
    for i, r in enumerate(rows):
        prev = rows[i - 1]["driving"] if i else False
        if r["driving"] and not prev:
            events.append(("start", r))
        if not r["driving"] and prev and i < len(rows) - 1:
            events.append(("stop", r))
    events.append(("finish" if rows[-1]["km"] >= end_km - 0.1 else "stop", rows[-1]))
    for kind, r in events:
        mk = {"start": ("*", "tab:green", 15), "stop": ("s", "tab:blue", 9),
              "finish": ("P", "tab:red", 14)}[kind]
        ax.plot(r["lon"], r["lat"], mk[0], color=mk[1], ms=mk[2], mec="k", mew=0.6, zorder=5)
        txt = f"{(r['t'] + LOCAL):%a %H:%M}  km {r['km']:.0f}"
        dx = -10 if kind == "start" else 10
        dy = -14 if kind == "finish" else 0
        ax.annotate(txt, (r["lon"], r["lat"]), xytext=(dx, dy), textcoords="offset points",
                    ha="right" if kind == "start" else "left", va="center", fontsize=8,
                    bbox=dict(boxstyle="round,pad=0.2", fc="white", ec=mk[1], alpha=0.85), zorder=6)
    ax.plot([], [], "*", color="tab:green", ms=12, mec="k", label="driving starts")
    ax.plot([], [], "s", color="tab:blue", ms=8, mec="k", label="stop (overnight / end)")
    ax.plot([], [], "P", color="tab:red", ms=11, mec="k", label="finish")
    ax.scatter([], [], c="orange", s=28, edgecolors="k", linewidths=0.3,
               label="car every 30 min (colour = sunshine)")
    ax.plot([], [], color="#d62728", lw=6, alpha=0.75, label="driven under cloud")
    ax.plot([], [], color="#2ca02c", lw=6, alpha=0.75, label="driven in clear sky")
    lats = [p[0] for p in route]
    lons = [p[1] for p in route]
    ax.set_xlim(min(lons) - 3, max(lons) + 3)
    ax.set_ylim(min(lats) - 1.5, max(lats) + 1.5)
    ax.set_aspect(1 / np.cos(np.radians(np.mean(lats))))
    ax.set_xlabel("longitude [°E]")
    ax.set_ylabel("latitude [°]")
    ax.grid(alpha=0.25)
    ax.legend(loc="lower left", fontsize=8)
    t0 = rows[0]["t"] + LOCAL
    ax.set_title(f"Where the car is: start {t0:%a %d.%m. %H:%M} local, {a.speed:.0f} km/h\n"
                 f"(simplified coastline)", fontsize=11)
    fig.tight_layout()
    fig.savefig(path, dpi=110)
    plt.close(fig)


# ---------------------------------------------------------------- main
def main():
    here = os.path.dirname(os.path.abspath(__file__))
    p = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("--start", required=True,
                   help='start time "YYYY-MM-DD HH:MM" (local ACST unless --utc) or "now"')
    p.add_argument("--speed", type=float, required=True, help="driving speed in km/h")
    p.add_argument("--start-km", type=float, default=0.0, help="where the car is at the start (km)")
    p.add_argument("--end-km", type=float, default=None, help="where to stop (default: finish)")
    p.add_argument("--utc", action="store_true", help="--start is given in UTC")
    p.add_argument("--drive-from", default="08:00", help="daily driving window start, local")
    p.add_argument("--drive-until", default="17:00", help="daily driving window end, local")
    p.add_argument("--no-driving-window", action="store_true", help="drive day and night")
    p.add_argument("--models", default="best_match", help="Open-Meteo models, comma-separated")
    p.add_argument("--cloud-index", type=float, default=0.7,
                   help="'cloud' = less than this share of the clear-sky sunshine (default 0.7)")
    p.add_argument("--max-days", type=float, default=5,
                   help="calculate at most this many days (from now, or from the start of a "
                        "replay); default 5")
    p.add_argument("--panel-area", type=float, default=None, help="solar panel area in m²")
    p.add_argument("--efficiency", type=float, default=None, help="panel+system efficiency (0–1)")
    p.add_argument("--route", default=os.path.join(here, "route_locations.py"))
    p.add_argument("--out", default=os.path.join(here, "trip_forecast"))
    a = p.parse_args()

    route = load_route(a.route)
    lat = np.array([q[0] for q in route])
    lon = np.array([q[1] for q in route])
    km = np.concatenate([[0.0], np.cumsum(haversine_km(lat[:-1], lon[:-1], lat[1:], lon[1:]))])
    end_km = min(a.end_km if a.end_km is not None else km[-1], km[-1])
    if not 0 <= a.start_km < end_km:
        sys.exit(f"--start-km must be between 0 and {end_km:.0f}")

    now_utc = dt.datetime.now(dt.timezone.utc).replace(tzinfo=None)
    if a.start.strip().lower() == "now":
        start_utc = now_utc.replace(second=0, microsecond=0)
    else:
        start = dt.datetime.strptime(a.start.strip(), "%Y-%m-%d %H:%M")
        start_utc = start if a.utc else start - LOCAL
    window = None if a.no_driving_window else (
        dt.datetime.strptime(a.drive_from, "%H:%M").time(),
        dt.datetime.strptime(a.drive_until, "%H:%M").time())

    # on the 5-minute grid, so the car leaves exactly at the start of the driving window
    rem = (start_utc.minute % 5) * 60 + start_utc.second
    if rem:
        start_utc += dt.timedelta(seconds=300 - rem)
    # at most --max-days: counted from now for a trip that starts now or later,
    # from its start for a replay of the past
    if start_utc >= now_utc - dt.timedelta(minutes=10):
        until_utc = now_utc + dt.timedelta(days=a.max_days)
        if start_utc >= until_utc:
            sys.exit(f"The start lies more than {a.max_days:g} days ahead – forecasts are not "
                     f"reliable that far (change with --max-days).")
    else:
        until_utc = start_utc + dt.timedelta(days=a.max_days)
    trip = simulate(start_utc, a.speed, a.start_km, end_km, window, until_utc)
    t_first, t_last = trip[0][0], trip[-1][0]
    finish = next((t for t, k, d in trip if k >= end_km), None)
    print(f"Route {km[-1]:.0f} km; trip km {a.start_km:.0f} → {end_km:.0f} at {a.speed:.0f} km/h, "
          f"start {t_first + LOCAL:%a %d.%m. %H:%M} local")
    print(f"  arrival: {(finish + LOCAL):%a %d.%m. %H:%M} local" if finish else
          f"  finish not reached within {a.max_days:g} days – stops at km {trip[-1][1]:.0f}, "
          f"{trip[-1][0] + LOCAL:%a %d.%m. %H:%M} local")

    # route points the car passes (plus neighbours), forecasts for the whole trip time
    i0 = max(0, int(np.searchsorted(km, a.start_km)) - 1)
    i1 = min(len(route), int(np.searchsorted(km, end_km)) + 1)
    pts = route[i0:i1]
    pkm = km[i0:i1]
    d0 = (t_first - dt.timedelta(hours=1)).date()
    d1 = (t_last + dt.timedelta(hours=1)).date()
    past = d1 < now_utc.date()
    url = HIST_URL if past else FC_URL
    base = dict(start_date=f"{d0:%Y-%m-%d}", end_date=f"{d1:%Y-%m-%d}", timezone="GMT")
    print(f"Fetching {'archived ' if past else ''}forecasts for {len(pts)} route points, "
          f"{d0} … {d1} …")
    series = {}
    for m in [x.strip() for x in a.models.split(",") if x.strip()]:
        s, block = get_series(url, pts, dict(base, models=m), m)
        if s:
            series[m] = s
            print(f"  {m}: {block.replace('_', ' ')} values")
    if past:
        s, block = get_series(SAT_URL, pts, dict(base, models="satellite_radiation_seamless"),
                              "actual (satellite)")
        if s:
            series["actual (satellite)"] = s
            print(f"  actual (satellite): {block.replace('_', ' ')} values")
    if not any(k != "actual (satellite)" for k in series):
        sys.exit("No forecast could be fetched – see the messages above.")

    # sunshine at the car: linear between the two route points around its km
    names = list(series)
    rows = []
    for t, k, driving in trip:
        j = int(np.clip(np.searchsorted(pkm, k) - 1, 0, len(pkm) - 2))
        w = float(np.clip((k - pkm[j]) / max(pkm[j + 1] - pkm[j], 1e-9), 0, 1))
        clat = (1 - w) * pts[j][0] + w * pts[j + 1][0]
        clon = (1 - w) * pts[j][1] + w * pts[j + 1][1]
        vals = {}
        for nme, s in series.items():
            va, vb = s.value(j, t), s.value(j + 1, t)
            vals[nme] = (1 - w) * va + w * vb if np.isfinite(va) and np.isfinite(vb) else \
                (va if np.isfinite(va) else vb)
        near = pts[j][2] if w < 0.5 else pts[j + 1][2]
        rows.append(dict(t=t, km=k, driving=driving, lat=clat, lon=clon, near=near,
                         cs=float(clear_sky(clat, clon, t)), **{f"v:{n}": vals[n] for n in names}))

    os.makedirs(a.out, exist_ok=True)
    fmt = lambda x: "" if x is None or not np.isfinite(x) else round(float(x))
    fc_name = next((n for n in names if n != "actual (satellite)"), names[0])
    flags_main = cloud_flags(rows, fc_name, a.cloud_index)
    with open(os.path.join(a.out, "trip_forecast.csv"), "w", newline="") as f:
        wr = csv.writer(f)
        wr.writerow(["time_local", "time_utc", "km", "nearest_point", "latitude", "longitude",
                     "status", "clear_sky_wm2"] + [f"ghi_{n.replace(' ', '_')}_wm2" for n in names]
                    + [f"cloud_{fc_name.replace(' ', '_')}"])
        for r, fl in zip(rows, flags_main):
            wr.writerow([f"{r['t'] + LOCAL:%Y-%m-%d %H:%M}", f"{r['t']:%Y-%m-%d %H:%M}",
                         round(r["km"], 1), r["near"], round(r["lat"], 4), round(r["lon"], 4),
                         "driving" if r["driving"] else "stopped", fmt(r["cs"])]
                        + [fmt(r[f"v:{n}"]) for n in names] + ["" if fl is None else int(fl)])

    # summary per local day: km driven, energy (all day, and while driving)
    h = STEP.total_seconds() / 3600
    kwh = (a.panel_area * a.efficiency) if a.panel_area and a.efficiency else None
    print("\nPer day (local)            km driven   cloud   sunshine energy [kWh/m²]"
          + ("   car [kWh]" if kwh else ""))
    days = sorted({(r["t"] + LOCAL).date() for r in rows})
    for d in days:
        rs = [r for r in rows if (r["t"] + LOCAL).date() == d]
        kms = [r["km"] for r in rs]
        parts = []
        for n in names:
            v = np.array([r[f"v:{n}"] for r in rs], float)
            e = np.nansum(v) * h / 1000
            miss = np.isnan(v).mean()
            parts.append(f"{n} {e:.2f}" + (" (gaps)" if miss > 0.05 else ""))
        csd = sum(r["cs"] for r in rs) * h / 1000
        idx = [i for i, r in enumerate(rows) if (r["t"] + LOCAL).date() == d]
        sub = rs + ([rows[idx[-1] + 1]] if idx[-1] + 1 < len(rows) else [])
        sh, _, _ = cloud_share(sub, [flags_main[i] for i in idx] + [None])
        line = (f"  {d:%a %d.%m.}  km {kms[0]:5.0f}→{kms[-1]:5.0f} ({kms[-1] - kms[0]:4.0f})  "
                + (f"{100 * sh:4.0f} %  " if sh is not None else "    –   ")
                + ", ".join(parts) + f"; clear sky {csd:.2f}")
        if kwh:
            first = np.array([r[f"v:{names[0]}"] for r in rs], float)
            line += f"   {names[0]}: {np.nansum(first) * h / 1000 * kwh:.1f} kWh"
        print(line)
    sh, under, tot = cloud_share(rows, flags_main)
    if sh is not None:
        print(f"  Cloud on the driven route ({fc_name}): {100 * sh:.0f} % ({under:.0f} of {tot:.0f} km, "
              f"cloud = below {100 * a.cloud_index:.0f} % of clear-sky sunshine)")

    # figure
    try:
        import matplotlib
        matplotlib.use("Agg")
        import matplotlib.dates as mdates
        import matplotlib.pyplot as plt
        tl = [r["t"] + LOCAL for r in rows]
        fig, ax = plt.subplots(2, 1, figsize=(15, 7.5), sharex=True,
                               gridspec_kw=dict(height_ratios=[2, 1]))
        ax[0].fill_between(tl, [r["cs"] for r in rows], color="0.9", label="clear sky")
        styles = {"actual (satellite)": dict(color="k", lw=1.8)}
        for n in names:
            ax[0].plot(tl, [r[f"v:{n}"] for r in rows], label=n, **styles.get(n, dict(lw=1.4)))
        stopped = np.array([not r["driving"] for r in rows])
        for axx in ax:
            axx.fill_between(tl, 0, 1, where=stopped, transform=axx.get_xaxis_transform(),
                             color="tab:blue", alpha=0.06, lw=0)
        ax[0].set_ylabel("sunshine at the car, GHI [W/m²]")
        ax[0].legend(loc="upper right", fontsize=9)
        ax[0].set_title(f"Start {tl[0]:%d.%m. %H:%M} local, {a.speed:.0f} km/h, "
                        f"km {a.start_km:.0f}→{end_km:.0f} (shaded = car stopped)", pad=22)
        ax[1].plot(tl, [r["km"] for r in rows], color="tab:green", lw=2)
        ax[1].set_ylabel("km along route")
        ax[1].set_xlabel("local time (ACST)")
        # time axis: a label every 1–6 h depending on the trip length, the date at midnight,
        # a small tick every hour, light vertical lines for orientation
        span_h = (tl[-1] - tl[0]).total_seconds() / 3600
        step_h = 1 if span_h <= 12 else 2 if span_h <= 48 else 3 if span_h <= 120 else 6
        loc = mdates.HourLocator(byhour=range(0, 24, step_h))
        ax[1].xaxis.set_major_locator(loc)
        ax[1].xaxis.set_minor_locator(mdates.HourLocator())
        ax[1].xaxis.set_major_formatter(mdates.DateFormatter("%H:%M"))
        for axx in ax:
            axx.grid(axis="x", which="major", alpha=0.3)
        # each day: a clear line at midnight and the day's name centred above that day
        d = tl[0].date()
        while dt.datetime.combine(d, dt.time()) <= tl[-1]:
            d0 = max(dt.datetime.combine(d, dt.time()), tl[0])
            d1 = min(dt.datetime.combine(d + dt.timedelta(days=1), dt.time()), tl[-1])
            if d0 > tl[0]:
                for axx in ax:
                    axx.axvline(d0, color="0.35", lw=1.2, zorder=0)
            if (d1 - d0).total_seconds() > 3 * 3600:
                ax[0].text(d0 + (d1 - d0) / 2, 1.01, f"{d:%a %d.%m.}", ha="center", va="bottom",
                           fontsize=10, fontweight="bold", transform=ax[0].get_xaxis_transform())
            d += dt.timedelta(days=1)
        plt.setp(ax[1].get_xticklabels(), fontsize=8)
        fig.tight_layout()
        fig.savefig(os.path.join(a.out, "trip_forecast.png"), dpi=110)
        plt.close(fig)

        save_trip_map(rows, names, pts, route, km, a, end_km,
                      os.path.join(a.out, "trip_map.png"))
    except ImportError:
        pass
    print(f"\n→ {a.out}/trip_forecast.csv, trip_forecast.png and trip_map.png")


if __name__ == "__main__":
    main()