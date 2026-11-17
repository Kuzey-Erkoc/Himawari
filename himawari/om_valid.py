#!/usr/bin/env python3
"""
How good are the weather models' sunshine forecasts along our race route?

For a period (default: the race week 22–29 Aug 2026, local dates) the script fetches, for
every point of route_locations.py (--every N: every N-th point):

  forecasts   hourly global irradiance (GHI, W/m²) of several weather models via Open-Meteo
              (ECMWF, GFS, the Australian ACCESS-G, JMA …; ICON etc. with --models)
                * "short lead": Historical Forecast API – the newest forecast for each hour,
                  made a few hours before
                * "1 day ahead" / "2 days ahead": Previous Runs API – what the model said
                  24 h / 48 h before (what the team would have seen the evening before)
  actual      hourly GHI derived from the Himawari satellite (Open-Meteo Satellite Radiation
              API) – a satellite estimate, not a ground measurement

and compares them hour by hour, day by day (daily energy in kWh/m², what the race strategy
cares about) and along the route. Two simple references show whether a forecast is worth
anything:
  * "always clear sky"   – no clouds at all (computed from the sun position)
  * "same as yesterday"  – the actual values of the day before

A model that Open-Meteo has no archive for (or with gaps) is left out with a note.

Output (folder validation_openmeteo_<start date>/):
  data_hourly.csv          all hourly values (point, time UTC and local, actual, every forecast)
  metrics_hourly.csv       bias, MAE, RMSE, correlation per forecast (daytime hours)
  metrics_daily.csv        error of the daily energy per forecast (complete local days),
                           also split north / south / race week / change days / steady days
  metrics_by_region.csv    hourly metrics split north / south of 23.4°S, race week,
                           weather-change days and steady days
  metrics_by_point.csv     error per route point and forecast
  metrics_by_section.csv   error per route section (Darwin – Katherine – … – Adelaide)
  bias_factors.csv         correction factor per forecast, and the factor that would fit
  summary.md               the main results in tables
  fig_timeseries.png       actual vs. ECMWF at Darwin, Alice Springs, Coober Pedy, Adelaide
  fig_along_route.png      error and bias of every model along the route (km)
  fig_days_along_route.png daily energy along the route, one panel per day: actual vs. models
  fig_errors.png           mean error per forecast and lead time
  fig_daily_energy.png     daily energy: forecast vs. actual
  fig_hour_of_day.png      bias by hour of day (does a model miss mornings/afternoons?)
  fig_regions_changes.png  error north vs. south, and on steady vs. weather-change days

Bias correction: all models turned out systematically too dark. Every forecast is therefore
also shown multiplied by a factor ("… | corrected"): with one period the fixed factors of
main.py / trip_forecast.py (BIAS_FACTOR, ECMWF × 1.06); with two or more periods a factor
learned on the other period (fair test). Weather-change days: days whose daily energy differs by
more than 1 kWh/m² from the day before (--change-kwh).

Open-Meteo's free use allows about 600 location-requests per minute and 10 000 per day; all 150
points × 5 models need roughly 1 700. When the minute limit is hit the script waits and goes on.

Usage:
    pip install numpy requests matplotlib
    python om_valid.py                                    # race week 22–29 Aug 2026, all points
    python om_valid.py --period 2025-08-22 2025-08-29     # race week 2025
    python om_valid.py --period 2025-08-15 2025-09-05 --period 2026-08-15 2026-09-05 --every 5
    python om_valid.py --models ecmwf_ifs025,gfs_global,bom_access_global
    python om_valid.py --day-lead "short lead"            # route plots with the newest forecasts
"""
import argparse
import ast
import collections
import concurrent.futures as cf
import csv
import datetime as dt
import json
import os
import sys
import threading
import time

import numpy as np

HIST_URL = "https://historical-forecast-api.open-meteo.com/v1/forecast"
PREV_URL = "https://previous-runs-api.open-meteo.com/v1/forecast"
SAT_URL = "https://satellite-api.open-meteo.com/v1/archive"
# weather models (Open-Meteo names). A model Open-Meteo has no archive for is left out with a note.
# (UK Met Office, Météo-France and ICON were left out after the race-week-2026 check: worse
# than ECMWF/GFS, false cloud alarms; add them with --models if wanted)
DEFAULT_MODELS = "best_match,ecmwf_ifs025,gfs_global,bom_access_global,jma_gsm"
DEFAULT_PERIODS = [("2026-08-22", "2026-08-29")]       # race week 2026, local dates (ACST)
# fixed bias-correction factors – the same as in main.py / trip_forecast.py (from the 2025+2026
# seasons and the car's weather station). Used when only one period is analysed.
BIAS_FACTOR = {"ecmwf_ifs025": 1.06, "gfs_global": 1.03, "icon_global": 1.09, "best_match": 1.06}
# towns along the route: marks in the along-route plots and the route sections
TOWNS = [("Darwin", -12.46, 130.84), ("Katherine", -14.47, 132.26), ("Tennant Creek", -19.65, 134.19),
         ("Alice Springs", -23.70, 133.88), ("Coober Pedy", -29.01, 134.76),
         ("Port Augusta", -32.49, 137.77), ("Adelaide", -34.93, 138.60)]
MIN_COVERAGE = 0.9          # a forecast needs values for 90 % of the daytime hours, else left out
RACE_WEEKS = [("2025-08-22", "2025-08-29"), ("2026-08-22", "2026-08-29")]
NORTH_LAT = -23.44          # Tropic of Capricorn
LOCAL_OFFSET_H = 9.5        # ACST: NT and SA in August (no daylight saving yet)
BATCH = 50
TIMEOUT = 60
RATE_WAIT = 65              # Open-Meteo allows ~600 location-requests per minute: wait and retry
UA = "himawari-nowcast-validation/1.0 (research)"
DOWN = set()


# ---------------------------------------------------------------- helpers
def load_route(path):
    """Route points [(lat, lon, name)] from a .py file with a list of dicts (read, not run),
    a .json list or a .csv with latitude/longitude columns."""
    ext = os.path.splitext(path)[1].lower()
    if ext == ".py":
        data = None
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
    if not data:
        sys.exit(f"No route points found in {path}")
    return [(float(d["latitude"]), float(d["longitude"]), str(d.get("name", f"P{i + 1}")))
            for i, d in enumerate(data)]


def local_date(t):
    """Local (ACST) date of the hour that ends at UTC stamp t."""
    return (t - dt.timedelta(hours=1) + dt.timedelta(hours=LOCAL_OFFSET_H)).date()


def haversine_km(lat1, lon1, lat2, lon2):
    p1, p2 = np.radians(lat1), np.radians(lat2)
    a = (np.sin((p2 - p1) / 2) ** 2
         + np.cos(p1) * np.cos(p2) * np.sin(np.radians(lon2 - lon1) / 2) ** 2)
    return 6371.0 * 2 * np.arcsin(np.sqrt(a))


def cos_zenith(lat, lon, t):
    la = np.radians(np.asarray(lat, float))
    hour = t.hour + t.minute / 60
    g = 2 * np.pi / 365 * (t.timetuple().tm_yday - 1 + (hour - 12) / 24)
    eqt = 229.18 * (0.000075 + 0.001868 * np.cos(g) - 0.032077 * np.sin(g)
                    - 0.014615 * np.cos(2 * g) - 0.040849 * np.sin(2 * g))
    decl = (0.006918 - 0.399912 * np.cos(g) + 0.070257 * np.sin(g) - 0.006758 * np.cos(2 * g)
            + 0.000907 * np.sin(2 * g) - 0.002697 * np.cos(3 * g) + 0.00148 * np.sin(3 * g))
    ha = np.radians((hour * 60 + eqt + 4 * np.asarray(lon, float)) / 4 - 180)
    return np.sin(la) * np.sin(decl) + np.cos(la) * np.cos(decl) * np.cos(ha)


def clear_sky_hour(lat, lon, t_end):
    """Clear-sky GHI averaged over the hour before t_end (Haurwitz, 6 sub-steps), W/m²."""
    vals = []
    for m in range(5, 60, 10):
        cz = cos_zenith(lat, lon, t_end - dt.timedelta(minutes=60 - m))
        czp = np.maximum(cz, 1e-3)
        vals.append(np.where(cz > 0, 1098.0 * czp * np.exp(-0.057 / czp), 0.0))
    return np.mean(vals, axis=0)


def fetch(url, pts, params, label):
    """All points (BATCH per request, parallel), per-point fallback if lists are refused.
    Returns list of JSON per point (None where missing) or None + prints the reason."""
    import requests
    host = url.split("/")[2]
    if host in DOWN:
        print(f"  ! {label}: {host} did not answer before – skipped")
        return None
    local = threading.local()
    err = {"msg": ""}
    lock = threading.Lock()
    note = {"rate": False}

    def get(idx):
        if not hasattr(local, "s"):
            local.s = requests.Session()
            local.s.headers["User-Agent"] = UA
        p = dict(params, latitude=",".join(f"{pts[i][0]:.4f}" for i in idx),
                 longitude=",".join(f"{pts[i][1]:.4f}" for i in idx))
        for attempt in range(4):
            try:
                r = local.s.get(url, params=p, timeout=TIMEOUT)
            except requests.RequestException as e:
                err["msg"] = f"no answer from {host} ({type(e).__name__})"
                DOWN.add(host)
                return None
            if r.status_code != 429 or attempt == 3:
                break
            # too many requests per minute: Open-Meteo's free limit – wait and try again
            with lock:
                if not note["rate"]:
                    print(f"  … Open-Meteo request limit reached, waiting {RATE_WAIT} s")
                    note["rate"] = True
            time.sleep(RATE_WAIT)
        if r.status_code != 200:
            err["msg"] = f"HTTP {r.status_code}: {r.text[:200]}"
            return None
        d = r.json()
        d = d if isinstance(d, list) else [d]
        return d if len(d) == len(idx) else None

    batches = [list(range(b, min(len(pts), b + BATCH))) for b in range(0, len(pts), BATCH)]
    with cf.ThreadPoolExecutor(2) as ex:
        res = list(ex.map(get, batches))
    if all(r is None for r in res) and host not in DOWN and len(batches[0]) > 1:
        if get(batches[0][:1]) is not None:          # lists refused -> one point per request
            with cf.ThreadPoolExecutor(2) as ex:
                single = list(ex.map(lambda i: get([i]), range(len(pts))))
            return [s[0] if s else None for s in single]
    if all(r is None for r in res):
        print(f"  ! {label}: {err['msg'] or 'no data'}")
        return None
    out = [None] * len(pts)
    for idx, r in zip(batches, res):
        if r:
            for i, d in zip(idx, r):
                out[i] = d
    return out


def to_grid(answers, var, times, n):
    """(n points × len(times)) array of var, NaN where missing."""
    pos = {t: j for j, t in enumerate(times)}
    arr = np.full((n, len(times)), np.nan)
    if answers is None:
        return arr
    for i, d in enumerate(answers):
        if not d:
            continue
        block = next((d[k] for k in ("hourly", "minutely_15") if isinstance(d.get(k), dict)
                      and var in d[k]), None)
        if block is None:
            continue
        for ts, v in zip(block["time"], block[var]):
            j = pos.get(dt.datetime.strptime(ts[:16], "%Y-%m-%dT%H:%M"))
            if j is not None and v is not None:
                arr[i, j] = v
    return arr


def metrics(f, y):
    ok = np.isfinite(f) & np.isfinite(y)
    if ok.sum() < 10:
        return None
    e = f[ok] - y[ok]
    return dict(n=int(ok.sum()), mean_actual=float(y[ok].mean()), bias=float(e.mean()),
                mae=float(np.abs(e).mean()), rmse=float(np.sqrt((e ** 2).mean())),
                rel_mae=float(np.abs(e).mean() / max(y[ok].mean(), 1e-9)),
                corr=float(np.corrcoef(f[ok], y[ok])[0, 1]) if np.std(f[ok]) > 0 else float("nan"))


# ---------------------------------------------------------------- main
def main():
    here = os.path.dirname(os.path.abspath(__file__))
    p = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("--route", default=os.path.join(here, "route_locations.py"))
    p.add_argument("--period", nargs=2, action="append", metavar=("START", "END"),
                   help="local dates (ACST) YYYY-MM-DD YYYY-MM-DD, both included (repeatable); "
                        "default the race week 22–29 Aug 2026")
    p.add_argument("--models", default=DEFAULT_MODELS)
    p.add_argument("--every", type=int, default=1,
                   help="use every N-th route point (default 1 = all points)")
    p.add_argument("--day-lead", default="1 day ahead", choices=["short lead", "1 day ahead", "2 days ahead"],
                   help="which forecast the day-by-day route plot shows (default: the one from the "
                        "day before, what the team plans with)")
    p.add_argument("--change-kwh", type=float, default=1.0,
                   help="a day counts as weather change if its energy differs from the day before "
                        "by more than this (kWh/m², default 1.0)")
    p.add_argument("--out", default=None,
                   help="output folder (default validation_openmeteo_<first start date>)")
    a = p.parse_args()

    route = load_route(a.route)
    km = np.concatenate([[0.0], np.cumsum([haversine_km(route[i][0], route[i][1], route[i + 1][0],
                                                         route[i + 1][1]) for i in range(len(route) - 1)])])
    sel = list(range(0, len(route), max(1, a.every)))
    if sel[-1] != len(route) - 1:
        sel.append(len(route) - 1)
    pts = [route[i] for i in sel]
    pkm = km[sel]
    lat = np.array([q[0] for q in pts])
    lon = np.array([q[1] for q in pts])
    north = lat > NORTH_LAT
    models = [m.strip() for m in a.models.split(",") if m.strip()]
    periods = [tuple(x) for x in (a.period or DEFAULT_PERIODS)]
    pdates = [(dt.date.fromisoformat(s), dt.date.fromisoformat(e)) for s, e in periods]
    if a.out is None:
        a.out = os.path.join(here, f"validation_openmeteo_{periods[0][0]}")
    os.makedirs(a.out, exist_ok=True)
    print(f"{len(pts)} route points ({'all' if a.every == 1 else f'every {a.every}th'}, "
          f"{pkm[-1]:.0f} km), models: {', '.join(models)}")

    # hourly time axis (stamp = end of the averaging hour, as Open-Meteo labels it). The periods
    # are LOCAL days (ACST): a local day runs from 14:30 UTC the day before. One extra local day
    # before each period is fetched for the reference 'same as yesterday'.
    times = []
    for d0, d1 in pdates:
        t = dt.datetime.combine(d0 - dt.timedelta(days=2), dt.time())
        t_end = dt.datetime.combine(d1 + dt.timedelta(days=1), dt.time())
        while t < t_end:
            if d0 - dt.timedelta(days=1) <= local_date(t) <= d1:
                times.append(t)
            t += dt.timedelta(hours=1)
    times = sorted(set(times))
    n, T = len(pts), len(times)

    series = collections.OrderedDict()      # name -> (n × T) array
    for (s, e), (d0, d1) in zip(periods, pdates):
        print(f"\nPeriod {s} … {e} (local dates)")
        base = dict(start_date=f"{d0 - dt.timedelta(days=2)}", end_date=e, timezone="GMT")
        print("  actual (satellite) …")
        sat = fetch(SAT_URL, pts, dict(base, hourly="shortwave_radiation",
                                       models="satellite_radiation_seamless"), "satellite")
        arr = to_grid(sat, "shortwave_radiation", times, n)
        series["actual"] = np.where(np.isfinite(series["actual"]), series["actual"], arr) \
            if "actual" in series else arr
        for m in models:
            print(f"  {m}: short lead …")
            h = fetch(HIST_URL, pts, dict(base, hourly="shortwave_radiation", models=m), f"{m} short lead")
            key = f"{m} | short lead"
            arr = to_grid(h, "shortwave_radiation", times, n)
            series[key] = np.where(np.isfinite(series[key]), series[key], arr) if key in series else arr
            print(f"  {m}: 1 and 2 days ahead …")
            vars_ = "shortwave_radiation_previous_day1,shortwave_radiation_previous_day2"
            pr = fetch(PREV_URL, pts, dict(base, hourly=vars_, models=m), f"{m} previous runs")
            for d, lab in ((1, "1 day ahead"), (2, "2 days ahead")):
                key = f"{m} | {lab}"
                arr = to_grid(pr, f"shortwave_radiation_previous_day{d}", times, n)
                series[key] = np.where(np.isfinite(series[key]), series[key], arr) if key in series else arr

    y = series["actual"]
    if np.isfinite(y).sum() == 0:
        sys.exit("No satellite values – nothing to compare. Check the messages above.")
    empty = [k for k, v in series.items() if np.isfinite(v).sum() == 0]
    for k in empty:
        print(f"  (no data for {k} – left out)")
        del series[k]

    # references: clear sky, and 'same as yesterday' (for day-ahead comparisons)
    cs = np.stack([clear_sky_hour(lat, lon, t) for t in times], axis=1)
    series["always clear sky"] = cs
    yday = np.full_like(y, np.nan)
    idx = {t: j for j, t in enumerate(times)}
    for j, t in enumerate(times):
        jj = idx.get(t - dt.timedelta(days=1))
        if jj is not None:
            yday[:, j] = y[:, jj]
    series["same as yesterday"] = yday

    # keep only the hours of the requested local days (the extra day was for 'yesterday')
    keep = np.array([any(d0 <= local_date(t) <= d1 for d0, d1 in pdates) for t in times])
    times = [t for t, k_ in zip(times, keep) if k_]
    T = len(times)
    series = collections.OrderedDict((k, v[:, keep]) for k, v in series.items())
    cs = cs[:, keep]
    y = series["actual"]
    day = cs > 50                       # daytime hours (sun clearly up)

    # a forecast with gaps would shrink the common hours of all others: leave it out
    base_h = day & np.isfinite(y)
    for k in [k for k in series if " | " in k]:
        cov = float(np.isfinite(series[k])[base_h].mean()) if base_h.any() else 0.0
        if cov < MIN_COVERAGE:
            print(f"  ({k}: values for only {100 * cov:.0f} % of the daytime hours – left out)")
            del series[k]
    # implausible values (e.g. a broken archive: far more or less sunshine than the satellite
    # saw on average) would spoil every table: leave such a forecast out
    for k in [k for k in series if " | " in k]:
        ok = base_h & np.isfinite(series[k])
        ratio = float(series[k][ok].sum() / max(y[ok].sum(), 1e-9)) if ok.any() else np.nan
        if not 0.6 <= ratio <= 1.5:
            print(f"  ({k}: implausible values – on average {ratio:.2f} × the satellite value – left out)")
            del series[k]
    models = [m for m in models if any(k.startswith(f"{m} | ") for k in series)]

    # ---- bias correction: the models are systematically too dark. One scale factor per
    #      model and lead time, always learned on the OTHER period (fitted on 2025 → applied
    #      to 2026 and vice versa), so the corrected scores are a fair test, not hindsight.
    #      With only one period there is nothing to learn on fairly, so the fixed factors of
    #      main.py / trip_forecast.py (BIAS_FACTOR) are used – exactly what the team would use.
    fold = np.zeros(T, int)
    plabel = []
    for pi, (d0, d1) in enumerate(pdates):
        fold[np.array([d0 <= local_date(t) <= d1 for t in times])] = pi
        plabel.append(f"{d0} – {d1}")
    cross = len(periods) >= 2
    raw = [k for k in series if " | " in k]
    fac_rows = []
    for k in raw:
        f = series[k]
        ok_all = day & np.isfinite(y) & np.isfinite(f)
        fitted = float(y[ok_all].sum() / max(f[ok_all].sum(), 1e-9)) if ok_all.sum() >= 50 else np.nan
        out = np.full_like(f, np.nan)
        if cross:
            for fo in range(int(fold.max()) + 1):
                train = ok_all & (fold != fo)[None, :]
                if train.sum() < 50:
                    continue
                fac = float(y[train].sum() / max(f[train].sum(), 1e-9))
                test = np.broadcast_to((fold == fo)[None, :], f.shape)
                out[test] = f[test] * fac
                fac_rows.append(dict(forecast=k, applied_to=plabel[fo], factor=round(fac, 3),
                                     fitted_on_these_data=round(fitted, 3)))
        else:
            fac = BIAS_FACTOR.get(k.split(" | ")[0])
            if fac:
                out = f * fac
            fac_rows.append(dict(forecast=k, applied_to="fixed (BIAS_FACTOR)" if fac else "– (no fixed factor)",
                                 factor=fac or "", fitted_on_these_data=round(fitted, 3)))
        if np.isfinite(out).any():
            series[f"{k} | corrected"] = out
    order = (["actual"] + [x for k in raw for x in (k, f"{k} | corrected")]
             + ["always clear sky", "same as yesterday"])
    series = collections.OrderedDict((k, series[k]) for k in order if k in series)
    with open(os.path.join(a.out, "bias_factors.csv"), "w", newline="") as f:
        w = csv.DictWriter(f, fieldnames=["forecast", "applied_to", "factor", "fitted_on_these_data"])
        w.writeheader()
        w.writerows(fac_rows)
    if fac_rows:
        print("\nBias correction (forecast × factor; 'fits here' = factor that would fit these data):")
        for r in fac_rows:
            fs = f"× {r['factor']:.3f}" if r["factor"] != "" else "none"
            print(f"  {r['forecast']:<45} {r['applied_to']}: {fs}   (fits here × {r['fitted_on_these_data']})")

    # ---- hourly data file
    with open(os.path.join(a.out, "data_hourly.csv"), "w", newline="") as f:
        w = csv.writer(f)
        w.writerow(["time_utc_end_of_hour", "time_local_end_of_hour", "point", "km", "latitude",
                    "longitude", "region",
                    "clear_sky_wm2"] + list(series))
        for j, t in enumerate(times):
            if not day[:, j].any():
                continue
            for i in range(n):
                w.writerow([t.strftime("%Y-%m-%d %H:%M"),
                            (t + dt.timedelta(hours=LOCAL_OFFSET_H)).strftime("%Y-%m-%d %H:%M"),
                            pts[i][2], round(float(pkm[i]), 1),
                            pts[i][0], pts[i][1], "north" if north[i] else "south",
                            round(float(cs[i, j]), 1)]
                           + ["" if not np.isfinite(v[i, j]) else round(float(v[i, j]), 1)
                              for v in series.values()])

    # ---- local days (ACST) and daily energy (kWh/m²). Only days that lie completely inside
    #      the fetched time range count (all 24 hours present) – the half days at the edges of
    #      each period would otherwise show up as nonsense values of ~0.2 kWh/m².
    local_day = [local_date(t) for t in times]
    dpos = collections.defaultdict(list)
    for j, d in enumerate(local_day):
        dpos[d].append(j)
    days = sorted(d for d in dpos if len(dpos[d]) == 24)
    dindex = {d: di for di, d in enumerate(days)}
    jd = np.array([dindex.get(d, -1) for d in local_day])

    def daily(v):
        """Daily energy (kWh/m²) over the daytime hours (clear sky > 50 W/m²). Up to 2 missing
        hours per point and day are filled with the day's mean clear-sky index – single gaps
        at sunrise/sunset in one model would otherwise throw the whole day away."""
        out = np.full((n, len(days)), np.nan)
        for di, d in enumerate(days):
            js = dpos[d]
            c = cs[:, js]
            need = c > 50
            vv = v[:, js]
            have = need & np.isfinite(vv)
            with np.errstate(invalid="ignore", divide="ignore"):
                kmean = np.where(have.any(axis=1),
                                 np.nansum(np.where(have, vv, 0), axis=1)
                                 / np.maximum(np.sum(np.where(have, c, 0), axis=1), 1e-9), np.nan)
            miss = (need & ~np.isfinite(vv)).sum(axis=1)
            filled = np.where(need & ~np.isfinite(vv), c * kmean[:, None], vv)
            complete = (miss <= 2) & have.any(axis=1) & np.isfinite(kmean)
            out[:, di] = np.where(complete, np.sum(np.where(need, filled, 0), axis=1) / 1000.0, np.nan)
        return out
    dy = daily(y)
    fc_names = [k for k in series if k != "actual"]
    dser = {k: daily(series[k]) for k in fc_names}

    # ---- weather-change days: the daily energy at a point differs from the day before by
    #      more than --change-kwh. 'Same as yesterday' is wrong on exactly these days, so here
    #      a weather model has to prove itself.
    change = np.zeros_like(dy, bool)
    steady = np.zeros_like(dy, bool)
    for di in range(1, len(days)):
        if days[di] - days[di - 1] == dt.timedelta(days=1):
            jump = np.abs(dy[:, di] - dy[:, di - 1])
            ok = np.isfinite(jump)
            change[:, di] = ok & (jump > a.change_kwh)
            steady[:, di] = ok & (jump <= a.change_kwh)
    # fair comparison for the daily energy too: only point-days where every forecast exists
    # a forecast that has a daily value for fewer than 80 % of the point-days is left out of the
    # daily comparison (instead of throwing those days away for all others)
    base_d = np.isfinite(dy)
    for k in list(dser):
        cov = float(np.isfinite(dser[k])[base_d].mean()) if base_d.any() else 0.0
        if cov < 0.8:
            print(f"  ({k}: daily energy for only {100 * cov:.0f} % of the point-days – "
                  f"left out of the daily comparison)")
            dser[k] = np.full_like(dy, np.nan)
    kept_d = [v for v in dser.values() if np.isfinite(v).any()]
    d_all = (np.all([np.isfinite(v) for v in kept_d], axis=0) if kept_d else base_d) & base_d
    dy = np.where(d_all, dy, np.nan)
    dser = {k: np.where(d_all, v, np.nan) for k, v in dser.items()}
    change &= d_all
    steady &= d_all
    print(f"\n{len(days)} days, {int(d_all.sum())} point-days with a daily value for every forecast; "
          f"point-days with a weather change "
          f"(> {a.change_kwh:g} kWh/m² vs. day before): {int(change.sum())} "
          f"(north {int(change[north].sum())}, south {int(change[~north].sum())}), "
          f"steady: {int(steady.sum())}")

    def to_hours(dm):                   # (points × days) mask → (points × hours) mask
        m = np.zeros((n, T), bool)
        ok = jd >= 0
        m[:, ok] = dm[:, jd[ok]]
        return m

    race = np.zeros(T, bool)
    race_d = np.zeros(len(days), bool)
    for s, e in RACE_WEEKS:
        r0, r1 = dt.date.fromisoformat(s), dt.date.fromisoformat(e)
        race |= np.array([r0 <= d <= r1 for d in local_day])
        race_d |= np.array([r0 <= d <= r1 for d in days], bool)
    ones_h = np.ones((n, T), bool)
    ones_d = np.ones((n, len(days)), bool)
    subsets_h = [("all", ones_h), ("north", north[:, None] & ones_h), ("south", ~north[:, None] & ones_h),
                 ("race week", race[None, :] & ones_h),
                 ("change days", to_hours(change)), ("steady days", to_hours(steady))]
    subsets_d = [("all", ones_d), ("north", north[:, None] & ones_d), ("south", ~north[:, None] & ones_d),
                 ("race week", race_d[None, :] & ones_d), ("change days", change), ("steady days", steady)]
    SUBS = [s for s, _ in subsets_h]

    def rnd(r):
        return {x: round(v, 3) if isinstance(v, float) else v for x, v in r.items()}

    # ---- hourly metrics (daytime hours where the actual value and the forecast exist;
    #      n in the tables says how many hours that was)
    common = day & np.isfinite(y)
    # fair comparison: every forecast is scored on the same hours (where all of them exist)
    all_fc = np.all([np.isfinite(series[k]) for k in fc_names], axis=0)
    rows_h = []
    for sub, msk in subsets_h:
        for k in fc_names:
            mm = common & msk & all_fc
            r = metrics(np.where(mm, series[k], np.nan), np.where(mm, y, np.nan))
            if r:
                rows_h.append(dict(subset=sub, forecast=k, **rnd(r)))
    keys = ["subset", "forecast", "n", "mean_actual", "bias", "mae", "rmse", "rel_mae", "corr"]
    with open(os.path.join(a.out, "metrics_hourly.csv"), "w", newline="") as f:
        w = csv.DictWriter(f, fieldnames=keys)
        w.writeheader()
        w.writerows([r for r in rows_h if r["subset"] == "all"])
    with open(os.path.join(a.out, "metrics_by_region.csv"), "w", newline="") as f:
        w = csv.DictWriter(f, fieldnames=keys)
        w.writeheader()
        w.writerows([r for r in rows_h if r["subset"] != "all"])

    # ---- daily metrics, same subsets
    rows_d = []
    for sub, msk in subsets_d:
        for k in fc_names:
            r = metrics(np.where(msk, dser[k], np.nan), np.where(msk, dy, np.nan))
            if r:
                rows_d.append(dict(subset=sub, forecast=k, **rnd(r)))
    with open(os.path.join(a.out, "metrics_daily.csv"), "w", newline="") as f:
        w = csv.DictWriter(f, fieldnames=keys)
        w.writeheader()
        w.writerows(rows_d)
    rows_d_all = [r for r in rows_d if r["subset"] == "all"]

    # ---- along the route: every route point, and route sections between the towns
    hmask = common & all_fc
    town_km = []
    for tname, tla, tlo in TOWNS:
        i = int(np.argmin(haversine_km(lat, lon, tla, tlo)))
        if haversine_km(lat[i], lon[i], tla, tlo) < 60:
            town_km.append((tname, float(pkm[i])))
    sections = [(f"{a_} – {b_}", ka, kb) for (a_, ka), (b_, kb) in zip(town_km[:-1], town_km[1:])]
    rows_pt = []
    for i in range(n):
        for k in fc_names:
            e_h = (series[k][i] - y[i])[hmask[i]]
            e_d = (dser[k][i] - dy[i])[np.isfinite(dser[k][i]) & np.isfinite(dy[i])]
            if len(e_h) < 5:
                continue
            rows_pt.append(dict(point=pts[i][2], km=round(float(pkm[i]), 1), latitude=pts[i][0],
                                longitude=pts[i][1], forecast=k,
                                hourly_mae=round(float(np.abs(e_h).mean()), 1),
                                hourly_bias=round(float(e_h.mean()), 1),
                                daily_mae=round(float(np.abs(e_d).mean()), 3) if len(e_d) else "",
                                daily_actual_mean=round(float(np.nanmean(dy[i])), 3)
                                if np.isfinite(dy[i]).any() else "",
                                daily_forecast_mean=round(float(np.nanmean(dser[k][i])), 3)
                                if np.isfinite(dser[k][i]).any() else ""))
    with open(os.path.join(a.out, "metrics_by_point.csv"), "w", newline="") as f:
        w = csv.DictWriter(f, fieldnames=list(rows_pt[0]) if rows_pt else ["point"])
        w.writeheader()
        w.writerows(rows_pt)
    rows_sec = []
    for sname, ka, kb in sections:
        pm = (pkm >= ka) & (pkm <= kb)
        for k in fc_names:
            mm = hmask & pm[:, None]
            rh = metrics(np.where(mm, series[k], np.nan), np.where(mm, y, np.nan))
            dm = pm[:, None] & np.isfinite(dy)
            rd = metrics(np.where(dm, dser[k], np.nan), np.where(dm, dy, np.nan))
            if rh:
                rows_sec.append(dict(section=sname, km_from=round(ka), km_to=round(kb), forecast=k,
                                     hourly_mae=round(rh["mae"], 1), hourly_bias=round(rh["bias"], 1),
                                     daily_mae=round(rd["mae"], 3) if rd else "",
                                     daily_bias=round(rd["bias"], 3) if rd else ""))
    with open(os.path.join(a.out, "metrics_by_section.csv"), "w", newline="") as f:
        w = csv.DictWriter(f, fieldnames=list(rows_sec[0]) if rows_sec else ["section"])
        w.writeheader()
        w.writerows(rows_sec)

    # ---- figures
    try:
        import matplotlib
        matplotlib.use("Agg")
        import matplotlib.pyplot as plt
        import matplotlib.dates as mdates
        pal = ["#1f77b4", "#ff7f0e", "#2ca02c", "#d62728", "#9467bd", "#8c564b", "#e377c2",
               "#17becf", "#bcbd22"]                      # no grey: grey = references
        mcol = {m: pal[i % len(pal)] for i, m in enumerate(models)}
        lstyle = {"short lead": "-", "1 day ahead": "--", "2 days ahead": ":"}

        def shown(m, lab):
            """The series a plot shows for model m: bias-corrected if there is a factor."""
            k = f"{m} | {lab} | corrected"
            if k in series:
                return k, f"{m} × {BIAS_FACTOR.get(m, 0):.2f}" if not cross else f"{m}, corrected"
            k = f"{m} | {lab}"
            return (k, m) if k in series else (None, None)
        loc = lambda tt: [t + dt.timedelta(hours=LOCAL_OFFSET_H) - dt.timedelta(minutes=30) for t in tt]

        def mark_towns(ax, labels=True):
            for tname, tk in town_km:
                ax.axvline(tk, color="0.75", lw=0.8, zorder=0)
                if labels:
                    ax.text(tk, 1.01, tname, transform=ax.get_xaxis_transform(), rotation=30,
                            fontsize=7, ha="left", va="bottom")

        # time series at four places (local time; middle of each hour)
        places = [("Darwin", -12.46, 130.84), ("Alice Springs", -23.70, 133.88),
                  ("Coober Pedy", -29.01, 134.76), ("Adelaide", -34.93, 138.60)]
        best = "ecmwf_ifs025" if "ecmwf_ifs025" in models else models[0]
        nper = len(periods)
        fig, axs = plt.subplots(len(places), nper, figsize=(8 * nper, 2.6 * len(places)), squeeze=False)
        for r_i, (pname, pla, plo) in enumerate(places):
            i = int(np.argmin(haversine_km(lat, lon, pla, plo)))
            for c_i, (d0, d1) in enumerate(pdates):
                ax = axs[r_i, c_i]
                js = [j for j, t in enumerate(times) if d0 <= local_date(t) <= d1]
                tt = loc([times[j] for j in js])
                ax.fill_between(tt, cs[i, js], color="0.9", label="clear sky")
                ax.plot(tt, y[i, js], color="k", lw=1.6, label="actual (satellite)")
                for m in models:
                    if m not in (best, "gfs_global"):
                        continue
                    for lab in (("short lead", "1 day ahead") if m == best else ("1 day ahead",)):
                        k, name = shown(m, lab)
                        if k:
                            ax.plot(tt, series[k][i, js], lstyle[lab], color=mcol[m], lw=1.1,
                                    label=f"{name}, {lab}")
                ax.set_title(f"{pname} ({pts[i][2]}, km {pkm[i]:.0f})", fontsize=10)
                ax.set_ylabel("GHI [W/m²]")
                ax.xaxis.set_major_formatter(mdates.DateFormatter("%d.%m"))
                if r_i == 0 and c_i == 0:
                    ax.legend(fontsize=7, loc="upper right", ncol=2)
        fig.suptitle("Sunshine at four places: actual (satellite) vs. forecast (local time ACST)")
        fig.tight_layout()
        fig.savefig(os.path.join(a.out, "fig_timeseries.png"), dpi=110)
        plt.close(fig)

        # along the route: error per route point, all models (forecast of the day before)
        def smooth(v, w=5):
            v = np.asarray(v, float)
            out = np.full_like(v, np.nan)
            for j in range(len(v)):
                seg = v[max(0, j - w // 2): j + w // 2 + 1]
                out[j] = np.nanmean(seg) if np.isfinite(seg).any() else np.nan
            return out
        lead_r = a.day_lead
        fig, axs = plt.subplots(3, 1, figsize=(13, 11), sharex=True)
        refs = [("always clear sky", "0.55", ":"), ("same as yesterday", "0.35", "-.")]
        lines = [(shown(m, lead_r)[0], mcol[m], "-", shown(m, lead_r)[1]) for m in models
                 if shown(m, lead_r)[0]]
        lines += [(k, c_, ls, k) for k, c_, ls in refs if k in series]
        for k, c_, ls, lab in lines:
            e = np.where(hmask, series[k] - y, np.nan)
            cnt = np.isfinite(e).sum(axis=1)
            mae = np.where(cnt > 0, np.nanmean(np.abs(e), axis=1), np.nan)
            bia = np.where(cnt > 0, np.nanmean(e, axis=1), np.nan)
            de = np.abs(dser[k] - dy)
            dmae = np.where(np.isfinite(de).any(axis=1), np.nanmean(de, axis=1), np.nan)
            lw = 2.0 if k.startswith(best) else 1.1
            axs[0].plot(pkm, smooth(mae), ls, color=c_, lw=lw, label=lab)
            axs[1].plot(pkm, smooth(bia), ls, color=c_, lw=lw, label=lab)
            axs[2].plot(pkm, smooth(dmae), ls, color=c_, lw=lw, label=lab)
        axs[0].set_ylabel("mean error per\ndaytime hour [W/m²]")
        axs[1].set_ylabel("bias [W/m²]\n(< 0: too dark)")
        axs[1].axhline(0, color="k", lw=0.8)
        axs[2].set_ylabel("mean error of the\ndaily energy [kWh/m²]")
        axs[2].set_xlabel("km along the route (Darwin → Adelaide)")
        for ax_i, ax in enumerate(axs):
            mark_towns(ax, labels=ax_i == 0)
            ax.grid(alpha=0.3)
        axs[0].legend(fontsize=7, ncol=4, loc="upper left")
        fig.suptitle(f"Forecast error along the route – forecasts {lead_r}, bias-corrected where a "
                     f"factor exists (smoothed over 5 points)", y=0.995)
        fig.tight_layout()
        fig.savefig(os.path.join(a.out, "fig_along_route.png"), dpi=110)
        plt.close(fig)

        # day by day: daily energy along the route, actual vs. every model
        nd_ = len(days)
        if nd_:
            ncol = 2
            nrow = (nd_ + 1) // 2
            fig, axs = plt.subplots(nrow, ncol, figsize=(15, 3.2 * nrow), sharex=True, sharey=True,
                                    squeeze=False)
            for di, d in enumerate(days):
                ax = axs[di // ncol, di % ncol]
                ax.plot(pkm, dser["always clear sky"][:, di], ":", color="0.55", lw=1, label="clear sky")
                k = f"{best} | {lead_r}"
                if k in dser and f"{k} | corrected" in dser:      # raw ECMWF, to see the correction
                    ax.plot(pkm, dser[k][:, di], ":", color=mcol[best], lw=1.0,
                            label=f"{best} (raw, without factor)")
                for m in models:
                    k, name = shown(m, lead_r)
                    if k and k in dser:
                        ax.plot(pkm, dser[k][:, di], color=mcol[m], lw=2.0 if m == best else 1.1,
                                alpha=1 if m == best else 0.85, label=name)
                ax.plot(pkm, dy[:, di], color="k", lw=2.2, label="actual (satellite)")
                ax.text(0.99, 0.04, f"{d:%a %d.%m.%Y}", transform=ax.transAxes, ha="right",
                        va="bottom", fontsize=10, fontweight="bold",
                        bbox=dict(facecolor="white", edgecolor="0.7", boxstyle="round,pad=0.25"))
                mark_towns(ax, labels=di < ncol)
                ax.grid(alpha=0.3)
            for j in range(nd_, nrow * ncol):
                axs[j // ncol, j % ncol].axis("off")
            for r_ in range(nrow):
                axs[r_, 0].set_ylabel("daily energy [kWh/m²]")
            for c_ in range(ncol):
                axs[-1, c_].set_xlabel("km along the route")
            h_, l_ = axs[0, 0].get_legend_handles_labels()
            fig.legend(h_, l_, loc="lower center", ncol=6, fontsize=8)
            fig.suptitle(f"Daily sunshine energy along the route: actual vs. forecast {lead_r}", y=0.995)
            fig.tight_layout(rect=(0, 0.04, 1, 1))
            fig.savefig(os.path.join(a.out, "fig_days_along_route.png"), dpi=110)
            plt.close(fig)

        # errors per forecast
        def style(bars, names):
            for b, k in zip(bars, names):
                if k.endswith("corrected"):
                    b.set_hatch("///"); b.set_edgecolor("white")
        allr = [r for r in rows_h if r["subset"] == "all"]
        names = [r["forecast"] for r in allr]
        labels = [k.replace(" | ", ", ") for k in names]
        colors = [mcol.get(k.split(" | ")[0], "0.5") for k in names]
        fig, ax = plt.subplots(1, 2, figsize=(14, 0.3 * len(names) + 2))
        yv = np.arange(len(names))
        style(ax[0].barh(yv, [r["mae"] for r in allr], color=colors), names)
        ax[0].set_yticks(yv, labels, fontsize=8); ax[0].invert_yaxis()
        ax[0].set_xlabel("mean error per daytime hour [W/m²]  (lower = better)")
        ax[0].set_title("Hourly error  (hatched = bias-corrected)")
        dd = {r["forecast"]: r for r in rows_d_all}
        style(ax[1].barh(yv, [dd[k]["mae"] if k in dd else 0 for k in names], color=colors), names)
        ax[1].set_yticks(yv, labels, fontsize=8); ax[1].invert_yaxis()
        ax[1].set_xlabel("mean error of the daily energy [kWh/m²]  (lower = better)")
        ax[1].set_title("Daily energy error")
        fig.suptitle("Open-Meteo forecasts vs. satellite-derived actual values (all route points)")
        fig.tight_layout()
        fig.savefig(os.path.join(a.out, "fig_errors.png"), dpi=110)
        plt.close(fig)

        # daily energy scatter
        fig, axs = plt.subplots(1, 3, figsize=(15, 5), sharex=True, sharey=True)
        for ax, lab in zip(axs, ("short lead", "1 day ahead", "2 days ahead")):
            for m in models:
                k, name = shown(m, lab)
                if k and k in dser:
                    ax.scatter(dy.ravel(), dser[k].ravel(), s=8, alpha=0.5, color=mcol[m], label=name)
            lim = np.nanmax(dy) * 1.1 if np.isfinite(dy).any() else 10
            ax.plot([0, lim], [0, lim], "k-", lw=0.8)
            ax.set_xlabel("actual daily energy [kWh/m²]")
            ax.set_title(lab)
        axs[0].set_ylabel("forecast daily energy [kWh/m²]")
        axs[0].legend(fontsize=8)
        fig.suptitle("Daily energy per route point and day, bias-corrected where a factor exists "
                     "(on the line = perfect)")
        fig.tight_layout()
        fig.savefig(os.path.join(a.out, "fig_daily_energy.png"), dpi=110)
        plt.close(fig)

        # bias by local hour of day
        lh = np.array([((t - dt.timedelta(minutes=30)).hour + LOCAL_OFFSET_H + (t - dt.timedelta(minutes=30)).minute / 60) % 24
                       for t in times])
        bins = np.arange(6, 19)
        fig, ax = plt.subplots(figsize=(10, 5))
        for m in models:
            lab = "short lead"
            k, name = shown(m, lab)
            if not k:
                continue
            b = []
            for h0 in bins:
                msk = common & ((lh >= h0) & (lh < h0 + 1))[None, :] & np.isfinite(series[k])
                b.append(np.nanmean((series[k] - y)[msk]) if msk.any() else np.nan)
            ax.plot(bins + 0.5, b, lstyle[lab], color=mcol[m], marker="o", ms=3, label=name)
        ax.axhline(0, color="k", lw=0.8)
        ax.set_xlabel("local time (ACST)"); ax.set_ylabel("bias forecast − actual [W/m²]")
        ax.set_title("Systematic error by time of day, newest forecasts (above 0 = too sunny)")
        ax.legend(fontsize=7, ncol=2)
        fig.tight_layout()
        fig.savefig(os.path.join(a.out, "fig_hour_of_day.png"), dpi=110)
        plt.close(fig)

        # north vs south (hourly) and steady vs change days (daily energy)
        def grouped(ax, rows, subs, scols, xlabel, title):
            look = {(r["subset"], r["forecast"]): r["mae"] for r in rows}
            nm = [k for k in fc_names if any((s_, k) in look for s_ in subs)]
            yy = np.arange(len(nm))
            h = 0.8 / len(subs)
            for si, (s_, c_) in enumerate(zip(subs, scols)):
                ax.barh(yy + (si - (len(subs) - 1) / 2) * h,
                        [look.get((s_, k), np.nan) for k in nm], height=h, color=c_, label=s_)
            ax.set_yticks(yy, [k.replace(" | ", ", ") for k in nm], fontsize=8)
            ax.invert_yaxis()
            ax.set_xlabel(xlabel); ax.set_title(title); ax.legend(fontsize=8)
            ax.grid(axis="x", alpha=0.3)
        fig, ax = plt.subplots(1, 2, figsize=(15, 0.3 * len(fc_names) + 2.5))
        grouped(ax[0], rows_h, ["north", "south"], ["#E69F00", "#0072B2"],
                "mean error per daytime hour [W/m²]", "Hourly error: north vs. south of 23.4°S")
        grouped(ax[1], rows_d, ["steady days", "change days"], ["0.6", "#D55E00"],
                "mean error of the daily energy [kWh/m²]",
                f"Daily energy: steady vs. change days (jump > {a.change_kwh:g} kWh/m²)")
        fig.tight_layout()
        fig.savefig(os.path.join(a.out, "fig_regions_changes.png"), dpi=110)
        plt.close(fig)
    except ImportError:
        print("matplotlib not installed – no figures")

    # ---- summary
    def table(rows, unit, scale=1.0, nd=0):
        lines = [f"| Forecast | Mean error [{unit}] | Bias [{unit}] | Relative error | Correlation |",
                 "|---|---|---|---|---|"]
        for r in sorted(rows, key=lambda r: r["mae"]):
            lines.append(f"| {r['forecast'].replace(' | ', ', ')} | {r['mae'] * scale:.{nd}f} | {r['bias'] * scale:+.{nd}f} | "
                         f"{100 * r['rel_mae']:.0f} % | {r['corr']:.2f} |")
        return "\n".join(lines)
    md = [f"# Open-Meteo sunshine forecasts vs. satellite-derived actual values",
          "",
          f"Periods (local dates): {', '.join(f'{s} – {e}' for s, e in periods)}; {n} route points; "
          f"models: {', '.join(models)}; "
          f"daytime hours only (clear-sky GHI > 50 W/m²); every forecast is scored on the same hours. Actual = Himawari-derived GHI "
          f"(Open-Meteo Satellite Radiation API), not a ground measurement.",
          "",
          "## Hourly values",
          table([r for r in rows_h if r["subset"] == "all"], "W/m²"),
          "",
          "## Daily energy (complete days only)",
          table(rows_d_all, "kWh/m²", 1.0, 2),
          ""]

    def split_table(rows, title, fmt):
        by = collections.defaultdict(dict)
        for r in rows:
            by[r["forecast"]][r["subset"]] = r["mae"]
        out = [title, "| Forecast | " + " | ".join(SUBS[1:]) + " |", "|---" * len(SUBS) + "|"]
        for k in fc_names:
            if k in by:
                out.append(f"| {k.replace(' | ', ', ')} | " + " | ".join(
                    fmt.format(by[k][s_]) if s_ in by[k] else "–" for s_ in SUBS[1:]) + " |")
        return out
    md += split_table(rows_h, "## Mean error per subset – hourly [W/m²]", "{:.0f}") + [""]
    md += split_table(rows_d, "## Mean error per subset – daily energy [kWh/m²]", "{:.2f}") + [""]
    md += [f"Change days: daily energy differs by more than {a.change_kwh:g} kWh/m² from the day "
           f"before at that point ({int(change.sum())} point-days; steady: {int(steady.sum())}).", ""]
    if rows_sec:
        lead_r = a.day_lead
        snames = [x[0] for x in sections]
        look = {(r["section"], r["forecast"]): r for r in rows_sec}
        for field, unit, fmt in (("hourly_mae", "W/m², per daytime hour", "{:.0f}"),
                                 ("daily_mae", "kWh/m², daily energy", "{:.2f}")):
            md += [f"## Mean error per route section – forecasts {lead_r} [{unit}]",
                   "| Forecast | " + " | ".join(snames) + " |", "|---" * (len(snames) + 1) + "|"]
            for k in fc_names:
                if not (k.endswith(lead_r) or k.endswith(f"{lead_r} | corrected") or "|" not in k):
                    continue
                vals = [look.get((sn, k), {}).get(field, "") for sn in snames]
                md.append(f"| {k.replace(' | ', ', ')} | "
                          + " | ".join(fmt.format(v) if v != "" else "–" for v in vals) + " |")
            md.append("")
    if fac_rows:
        md += ["## Bias correction (forecast × factor)",
               "| Forecast | applied | factor | factor that would fit these data |", "|---|---|---|---|"]
        md += [f"| {r['forecast'].replace(' | ', ', ')} | {r['applied_to']} | "
               f"{r['factor'] if r['factor'] != '' else '–'} | {r['fitted_on_these_data']} |"
               for r in fac_rows]
        md.append("")
    md += ["Bias > 0: forecast too sunny. 'always clear sky' and 'same as yesterday' are the "
           "references a forecast has to beat. 'corrected' = forecast multiplied by "
           + ("a factor learned on the other period (fair test)." if cross else
              "the fixed factor of main.py / trip_forecast.py (BIAS_FACTOR). These factors were "
              "derived partly from the August 2026 satellite data, so for 2026 they are not a "
              "completely independent test.")]
    open(os.path.join(a.out, "summary.md"), "w", encoding="utf-8").write("\n".join(md))
    print("\n" + "\n".join(md))
    print(f"\n→ results in {a.out}/")


if __name__ == "__main__":
    main()
