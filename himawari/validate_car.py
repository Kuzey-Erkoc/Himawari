#!/usr/bin/env python3
"""
Which sunshine forecast is most reliable – checked against our own car's weather station?

The car's weather station measured global irradiance (GHI, W/m²) every ~1.6 s during the
race days (station_data/*.jsonl or *.jsonl.zip). That is an independent ground truth: unlike
the Himawari satellite estimate it is not derived from the same images our own nowcast uses.

The station logs no GPS position, so the car's position along the route is rebuilt:
  * car speed = true wind (absolute_wind_*) − apparent wind (wind_* + compass), as vectors;
    the sign convention is picked per file so that the resulting heading matches the compass
  * the distance driven is anchored to the WSC control stops (each a 30-minute stop at a
    known place): Darwin start, Katherine, Daly Waters, … , Port Augusta, Adelaide finish
  * check: the station's pressure altitude vs. the terrain height along the route
    (fig_track.png, summary.md)

Compared hour by hour along the car's track (each minute takes the forecast at the place
where the car was in that minute):
  satellite analysis      Himawari-derived GHI (Open-Meteo Satellite API) – what our own
                          analysis starts from (lead time 0)
  satellite +1/+2/+3 h    'the satellite clear-sky index stays': cloudiness seen by the
                          satellite 1–3 h earlier at the place where the car will be. Our
                          own nowcast scores about the same as this (skill_curve.csv).
  car sensor +1/+2/+3 h   'what our sensor measured stays' (clear-sky index of the car)
  weather models          short lead / 1 day ahead / 2 days ahead (Open-Meteo)
  always clear sky        no clouds at all

All methods are scored on the same hours (common sample), with 90 % bootstrap intervals.

Output (folder validation_car/):
  track_minutes.csv   per minute: time, km, lat/lon, speed, measured GHI, clear sky
  hourly.csv          per hour: measured GHI and every forecast along the track
  metrics.csv         error per method (common sample), with confidence intervals
  metrics_by_lead.csv satellite / car persistence vs. models per lead time
  daily_energy.csv    energy while driving per race day, measured vs. forecast
  summary.md          the results in words and tables
  fig_track.png       reconstructed position and the altitude check
  fig_days.png        measured vs. satellite vs. model, per race day
  fig_errors.png      error per method
  fig_lead.png        error vs. lead time

Usage:
    pip install numpy requests matplotlib
    python validate_car.py                         # station_data/ next to this script
    python validate_car.py --station a.jsonl b.zip # other files
    python validate_car.py --models ecmwf_ifs025,gfs_global
Open-Meteo answers are cached in validation_car/openmeteo_cache.json, so a rerun makes no
new requests (Open-Meteo limits free requests per day).
"""
import argparse
import collections
import csv
import datetime as dt
import glob
import hashlib
import io
import json
import os
import sys
import zipfile

import numpy as np

from om_valid import HIST_URL, PREV_URL, SAT_URL, fetch, to_grid, load_route, haversine_km, cos_zenith

LOCAL_OFFSET_H = 9.5            # ACST (NT and SA in August)
NORTH_LAT = -23.44
DEFAULT_MODELS = "best_match,ecmwf_ifs025,gfs_global,icon_global"
# WSC control stops in driving order (approximate positions, projected onto the route)
CONTROL_STOPS = [("Katherine", -14.465, 132.264), ("Daly Waters", -16.300, 133.380),
                 ("Tennant Creek", -19.648, 134.190), ("Barrow Creek", -21.530, 133.887),
                 ("Alice Springs", -23.700, 133.881), ("Kulgera", -25.839, 133.295),
                 ("Coober Pedy", -29.013, 134.754), ("Glendambo", -30.970, 135.750),
                 ("Port Augusta", -32.493, 137.765)]
MOVING_MS = 8.0                 # car counts as moving above this speed (m/s, 5-min median);
                                # parked, gusts make the reconstructed speed jump up to ~5-10 m/s
STOP_MIN, STOP_MAX = 15, 60     # a control stop lasts 30 min (+ slowing down / leaving)
STOP_GAP = 4                    # stop pieces less than this many minutes apart are one stop
LEADS_H = (1, 2, 3)
MAX_GAP_KM = 60                 # interpolate forecasts only between route points this close


# ---------------------------------------------------------------- station data
def read_jsonl(path):
    """{metric name: (time s, values)} from a Prometheus-style export (.jsonl or .zip)."""
    if path.lower().endswith(".zip"):
        with zipfile.ZipFile(path) as z:
            name = next(n for n in z.namelist() if n.endswith(".jsonl"))
            lines = io.TextIOWrapper(z.open(name), encoding="utf-8").readlines()
    else:
        lines = open(path, encoding="utf-8").readlines()
    out = {}
    for line in lines:
        if not line.strip():
            continue
        d = json.loads(line)
        k = d["metric"]["__name__"].replace("weather_station_", "")
        t = np.asarray(d["timestamps"], float) / 1000.0
        v = np.array([np.nan if x is None else x for x in d["values"]], float)
        out[k] = (t, v)
    return out


def car_velocity(m):
    """Car speed (m/s) at the wind_speed timestamps, from true minus apparent wind.
    Tries the possible direction conventions and keeps the one whose heading matches
    the compass best. Returns (t, speed, median heading error in degrees)."""
    t, rs = m["wind_speed"]
    rd = np.interp(t, *m["wind_direction"])
    cp = np.interp(t, m["compass"][0], np.degrees(np.unwrap(np.radians(m["compass"][1])))) % 360
    if "absolute_wind_speed" not in m:
        return t, np.where(rs > MOVING_MS, rs, 0.0), float("nan")
    ta, aws = m["absolute_wind_speed"]
    tb, awd = m["absolute_wind_direction"]
    awsi = np.interp(t, ta, aws, left=np.nan, right=np.nan)
    awdi = np.interp(t, tb, awd, left=np.nan, right=np.nan)
    best = None
    for sr in (1, -1):
        for fr in (0, 180):
            for ft in (0, 180):
                th = np.radians(cp + sr * rd + fr)
                tw = np.radians(awdi + ft)
                vx = -awsi * np.sin(tw) + rs * np.sin(th)
                vy = -awsi * np.cos(tw) + rs * np.cos(th)
                v = np.hypot(vx, vy)
                ok = np.isfinite(v) & (v > 5)
                if ok.sum() < 100:
                    continue
                hd = np.degrees(np.arctan2(vx, vy))
                err = float(np.median(np.abs((hd - cp + 180) % 360 - 180)[ok]))
                if best is None or err < best[0]:
                    best = (err, v)
    if best is None:
        return t, np.where(rs > MOVING_MS, rs, 0.0), float("nan")
    v = np.where(np.isfinite(best[1]), best[1], 0.0)     # no true wind -> parked
    return t, np.clip(v, 0, 45), best[0]


def per_minute(t, v, t0, n):
    """Mean of v per minute (index 0 = minute starting at t0), NaN where empty."""
    i = np.floor((t - t0) / 60).astype(int)
    ok = (i >= 0) & (i < n) & np.isfinite(v)
    s = np.bincount(i[ok], v[ok], n)
    c = np.bincount(i[ok], None, n)
    with np.errstate(invalid="ignore"):
        return np.where(c > 0, s / np.maximum(c, 1), np.nan), c


def load_station(paths):
    """Minute series over all files: dict of arrays + per-file info."""
    files = [read_jsonl(p) for p in paths]
    t0 = np.floor(min(f["global_radiation"][0].min() for f in files) / 60) * 60
    t1 = max(f["global_radiation"][0].max() for f in files)
    n = int((t1 - t0) // 60) + 1
    ghi = np.full(n, np.nan); spd = np.full(n, np.nan); prs = np.full(n, np.nan)
    tmp = np.full(n, np.nan); cnt = np.zeros(n)
    info = []
    for p, m in zip(paths, files):
        g, c = per_minute(*m["global_radiation"], t0, n)
        tv, v, herr = car_velocity(m)
        s, _ = per_minute(tv, v, t0, n)
        pr, _ = per_minute(*m["absolute_air_pressure"], t0, n)
        te, _ = per_minute(*m["air_temperature"], t0, n)
        for dst, src in ((ghi, g), (spd, s), (prs, pr), (tmp, te)):
            dst[:] = np.where(np.isfinite(src), src, dst)
        cnt += c
        info.append(dict(file=os.path.basename(p), heading_error_deg=round(herr, 1),
                         samples=int(len(m["global_radiation"][0]))))
    times = np.array([t0 + 60 * i for i in range(n)])
    return dict(t=times, ghi=ghi, speed=spd, pressure=prs, temp=tmp, n=cnt), info


# ---------------------------------------------------------------- track
def route_km(route):
    seg = [haversine_km(route[i][0], route[i][1], route[i + 1][0], route[i + 1][1])
           for i in range(len(route) - 1)]
    return np.concatenate([[0.0], np.cumsum(seg)])


def project_km(route, km, lat, lon):
    """km along the route of the point nearest to (lat, lon)."""
    best = (1e9, 0.0)
    for i in range(len(route) - 1):
        a, b = np.array(route[i][:2]), np.array(route[i + 1][:2])
        for f in np.linspace(0, 1, 21):
            q = a + f * (b - a)
            d = haversine_km(lat, lon, q[0], q[1])
            if d < best[0]:
                best = (d, km[i] + f * (km[i + 1] - km[i]))
    return best[1]


def find_stops(st):
    """Runs where the car stood still between driving, per local day.
    Returns list of dicts (start index, end index, minutes, local day)."""
    v = np.nan_to_num(st["speed"])
    pad = np.pad(v, 2, mode="edge")
    med = np.median(np.lib.stride_tricks.sliding_window_view(pad, 5), axis=1)
    moving = med > MOVING_MS
    # close short gaps inside a stop (a few gusty minutes)
    i = 0
    while i < len(moving):
        if moving[i]:
            j = i
            while j < len(moving) and moving[j]:
                j += 1
            if j - i < STOP_GAP and i > 0 and j < len(moving):
                moving[i:j] = False
            i = j
        else:
            i += 1
    local_day = np.array([(dt.datetime.fromtimestamp(x, dt.timezone.utc)
                           + dt.timedelta(hours=LOCAL_OFFSET_H)).date() for x in st["t"]])
    stops = []
    for d in sorted(set(local_day)):
        idx = np.where(local_day == d)[0]
        mv = idx[moving[idx]]
        if len(mv) < 30:
            continue
        a, b = mv[0], mv[-1]
        i = a
        while i <= b:
            if not moving[i]:
                j = i
                while j <= b and not moving[j]:
                    j += 1
                stops.append(dict(i0=int(i), i1=int(j - 1), minutes=int(j - i), day=d))
                i = j
            else:
                i += 1
    return stops, moving, local_day


def build_track(st, route, km):
    """km along the route for every minute, anchored to start, control stops and finish."""
    stops, moving, local_day = find_stops(st)
    v = np.nan_to_num(st["speed"]) * moving
    dist = np.concatenate([[0.0], np.cumsum(v[:-1] * 60 / 1000.0)])   # km driven
    cand = [s for s in stops if STOP_MIN <= s["minutes"] <= STOP_MAX]
    stop_km = [project_km(route, km, la, lo) for _, la, lo in CONTROL_STOPS]
    # choose len(CONTROL_STOPS) candidates in order whose driven distances fit the
    # expected route km best (dynamic programming over the ordered candidates)
    first = int(np.argmax(moving))
    last = int(len(moving) - 1 - np.argmax(moving[::-1]))
    scale = (km[-1] - 0.0) / max(dist[last] - dist[first], 1)
    K, C = len(CONTROL_STOPS), len(cand)
    cd = np.array([(dist[c["i0"]] - dist[first]) * scale for c in cand])
    INF = 1e18
    cost = np.full((K + 1, C + 1), INF)
    back = np.zeros((K + 1, C + 1), int)
    cost[0, :] = 0
    for k in range(1, K + 1):
        for c in range(k, C + 1):
            # stop k is candidate c-1; previous stop is some candidate < c-1
            prev = cost[k - 1, :c]
            j = int(np.argmin(prev))
            cost[k, c] = prev[j] + (cd[c - 1] - stop_km[k - 1]) ** 2 \
                + 5.0 * (cand[c - 1]["minutes"] - 30) ** 2
            back[k, c] = j
    if C < K:
        sys.exit(f"Only {C} control-stop candidates found, need {K} – check the station files.")
    c = int(np.argmin(cost[K, K:])) + K
    chosen = []
    for k in range(K, 0, -1):
        chosen.append(cand[c - 1])
        c = back[k, c]
    chosen = chosen[::-1]
    anchors_d = [dist[first]] + [dist[s["i0"]] for s in chosen] + [dist[last]]
    anchors_km = [0.0] + stop_km + [km[-1]]
    pos = np.interp(dist, anchors_d, anchors_km)
    lat = np.interp(pos, km, [p[0] for p in route])
    lon = np.interp(pos, km, [p[1] for p in route])
    seg_scale = [(anchors_km[i + 1] - anchors_km[i]) / max(anchors_d[i + 1] - anchors_d[i], 1e-6)
                 for i in range(len(anchors_d) - 1)]
    anchors = [dict(name="Darwin (start)", i=first, km=0.0)] + \
        [dict(name=n, i=s["i0"], km=k, minutes=s["minutes"]) for (n, _, _), s, k in
         zip(CONTROL_STOPS, chosen, stop_km)] + [dict(name="Adelaide (finish)", i=last, km=km[-1])]
    return dict(km=pos, lat=lat, lon=lon, moving=moving, local_day=local_day,
                anchors=anchors, seg_scale=seg_scale, stops=stops)


# ---------------------------------------------------------------- clear sky
def clear_sky_minute(lat, lon, t_epoch):
    """Clear-sky GHI (Haurwitz) at the middle of each minute, W/m²."""
    out = np.zeros(len(t_epoch))
    for i, (la, lo, x) in enumerate(zip(lat, lon, t_epoch)):
        tt = dt.datetime.fromtimestamp(x + 30, dt.timezone.utc).replace(tzinfo=None)
        cz = float(cos_zenith(la, lo, tt))
        out[i] = 1098.0 * cz * np.exp(-0.057 / cz) if cz > 1e-3 else 0.0
    return out


def clear_sky_grid(lat, lon, times):
    """(points × hours) clear-sky GHI averaged over the hour before each stamp."""
    from om_valid import clear_sky_hour
    return np.stack([clear_sky_hour(lat, lon, t) for t in times], axis=1)


# ---------------------------------------------------------------- open-meteo (cached)
class Cache:
    """Open-Meteo answers on disk, stored per batch of points, so an interrupted run
    (daily request limit, dropped connection) only fetches what is still missing."""
    BATCH = 25

    def __init__(self, path):
        self.path = path
        try:
            self.d = json.load(open(path, encoding="utf-8"))
        except (OSError, ValueError):
            self.d = {}

    def fetch(self, url, pts, params, label, tries=3):
        import time
        import om_valid
        out = []
        for b in range(0, len(pts), self.BATCH):
            sub = pts[b:b + self.BATCH]
            key = hashlib.sha1(json.dumps([url, [p[:2] for p in sub],
                                           sorted(params.items())]).encode()).hexdigest()
            r = self.d.get(key)
            for k in range(tries if r is None else 0):
                om_valid.DOWN.discard(url.split("/")[2])
                r = fetch(url, sub, params, f"{label} (points {b + 1}-{b + len(sub)})")
                if r is not None and all(x is not None for x in r):
                    self.d[key] = r
                    json.dump(self.d, open(self.path, "w", encoding="utf-8"))
                    break
                r = None
                time.sleep(5 * (k + 1))
            out += r if r is not None else [None] * len(sub)
        return out if any(x is not None for x in out) else None


# ---------------------------------------------------------------- statistics
def block_bootstrap(err_list, n_rep=2000, block=3, seed=1):
    """90 % interval of MAE for each error vector (same hours, paired), moving blocks
    of `block` consecutive hours. Returns (lo, hi) per vector and the replicate means."""
    rng = np.random.default_rng(seed)
    E = np.abs(np.vstack(err_list))           # methods × hours
    n = E.shape[1]
    nb = int(np.ceil(n / block))
    reps = np.empty((n_rep, E.shape[0]))
    for r in range(n_rep):
        starts = rng.integers(0, max(n - block + 1, 1), nb)
        idx = np.concatenate([np.arange(s, min(s + block, n)) for s in starts])[:n]
        reps[r] = E[:, idx].mean(axis=1)
    return reps


def stats(f, y):
    e = f - y
    return dict(n=int(len(e)), mean_measured=float(y.mean()), bias=float(e.mean()),
                mae=float(np.abs(e).mean()), rmse=float(np.sqrt((e ** 2).mean())),
                rel_mae=float(np.abs(e).mean() / max(y.mean(), 1e-9)),
                corr=float(np.corrcoef(f, y)[0, 1]) if np.std(f) > 0 else float("nan"))


# ---------------------------------------------------------------- main
def main():
    here = os.path.dirname(os.path.abspath(__file__))
    p = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("--station", nargs="+", help="station files (.jsonl or .jsonl.zip); "
                   "default: station_data/ next to this script")
    p.add_argument("--route", default=os.path.join(here, "route_locations.py"))
    p.add_argument("--models", default=DEFAULT_MODELS)
    p.add_argument("--every", type=int, default=2, help="use every N-th route point for the "
                   "Open-Meteo requests (default 2, ~40 km apart; values in between are interpolated)")
    p.add_argument("--out", default=os.path.join(here, "validation_car"))
    a = p.parse_args()

    paths = a.station or sorted(glob.glob(os.path.join(here, "station_data", "*.jsonl*")))
    paths = sorted(paths, key=lambda x: read_first_time(x))
    if not paths:
        sys.exit("No station files – put them into station_data/ or pass --station")
    os.makedirs(a.out, exist_ok=True)
    route = load_route(a.route)
    km = route_km(route)
    models = [m.strip() for m in a.models.split(",") if m.strip()]

    print(f"Reading {len(paths)} station files …")
    st, info = load_station(paths)
    for x in info:
        print(f"  {x['file']}: {x['samples']} samples, heading vs compass {x['heading_error_deg']}°")
    tr = build_track(st, route, km)
    print("Track anchors:")
    for an in tr["anchors"]:
        tl = dt.datetime.fromtimestamp(st["t"][an["i"]], dt.timezone.utc) + dt.timedelta(hours=LOCAL_OFFSET_H)
        print(f"  {an['name']:<18} km {an['km']:6.0f}  {tl:%d.%m. %H:%M} ACST"
              + (f"  ({an['minutes']} min stop)" if "minutes" in an else ""))
    print("  distance scale per leg (route km / reconstructed km): "
          + ", ".join(f"{s:.2f}" for s in tr["seg_scale"]))

    t = st["t"]
    cs_min = clear_sky_minute(tr["lat"], tr["lon"], t)

    # ---- Open-Meteo data for route points over the race days
    sel = list(range(0, len(route), max(1, a.every)))
    if sel[-1] != len(route) - 1:
        sel.append(len(route) - 1)
    pts = [route[i] for i in sel]
    pkm = km[sel]
    plat = np.array([q[0] for q in pts]); plon = np.array([q[1] for q in pts])
    d0 = dt.datetime.fromtimestamp(t[0], dt.timezone.utc).date()
    d1 = dt.datetime.fromtimestamp(t[-1], dt.timezone.utc).date()
    base = dict(start_date=f"{d0:%Y-%m-%d}", end_date=f"{d1:%Y-%m-%d}", timezone="GMT")
    times = []
    x = dt.datetime.combine(d0, dt.time())
    while x < dt.datetime.combine(d1 + dt.timedelta(days=1), dt.time()):
        times.append(x); x += dt.timedelta(hours=1)
    cache = Cache(os.path.join(a.out, "openmeteo_cache.json"))
    print(f"\nOpen-Meteo: {len(pts)} route points, {base['start_date']} … {base['end_date']}")
    grids = collections.OrderedDict()
    sat = cache.fetch(SAT_URL, pts, dict(base, hourly="shortwave_radiation",
                                          models="satellite_radiation_seamless"), "satellite")
    grids["satellite analysis"] = to_grid(sat, "shortwave_radiation", times, len(pts))
    elev = np.array([d.get("elevation", np.nan) if d else np.nan for d in (sat or [None] * len(pts))], float)
    for m in models:
        h = cache.fetch(HIST_URL, pts, dict(base, hourly="shortwave_radiation", models=m), f"{m} short lead")
        grids[f"{m} | short lead"] = to_grid(h, "shortwave_radiation", times, len(pts))
        pr = cache.fetch(PREV_URL, pts, dict(base, hourly="shortwave_radiation_previous_day1,"
                                             "shortwave_radiation_previous_day2", models=m), f"{m} previous runs")
        for dd, lab in ((1, "1 day ahead"), (2, "2 days ahead")):
            grids[f"{m} | {lab}"] = to_grid(pr, f"shortwave_radiation_previous_day{dd}", times, len(pts))
    for k in [k for k, g in grids.items() if not np.isfinite(g).any()]:
        print(f"  (no data for {k} – left out)")
        del grids[k]
    if "satellite analysis" not in grids:
        sys.exit("No satellite data – nothing to compare (see messages above).")
    csg = clear_sky_grid(plat, plon, times)

    # ---- along-track values per minute
    tpos = {int(x.replace(tzinfo=dt.timezone.utc).timestamp()): j for j, x in enumerate(times)}
    hour_stamp = (np.floor(t / 3600) * 3600 + 3600).astype(int)
    jh = np.array([tpos.get(int(h), -1) for h in hour_stamp])

    def along(grid, shift_h=0):
        """value of grid at hour stamp (H - shift) at the car's position of each minute."""
        out = np.full(len(t), np.nan)
        for i in range(len(t)):
            j = jh[i] - shift_h
            if jh[i] < 0 or j < 0:
                continue
            col = grid[:, j]
            ok = np.isfinite(col)
            if ok.sum() < 2:
                continue
            x = tr["km"][i]
            kk = pkm[ok]
            n = np.searchsorted(kk, x)
            gap = (kk[min(n, len(kk) - 1)] - kk[max(n - 1, 0)]) if 0 < n < len(kk) else np.inf
            if gap <= MAX_GAP_KM:          # only between nearby points with data
                out[i] = np.interp(x, kk, col[ok])
        return out

    cs_track_h = along(csg)                      # hourly clear sky along the track
    minute_series = collections.OrderedDict()
    for k, g in grids.items():
        minute_series[k] = along(g)
    ksat = np.where(csg > 50, np.clip(grids["satellite analysis"] / np.maximum(csg, 1), 0, 1.3), np.nan)
    for L in LEADS_H:
        minute_series[f"satellite +{L} h (cloudiness stays)"] = along(ksat, L) * cs_track_h
    minute_series["always clear sky"] = cs_track_h

    # ---- hourly means along the track (only hours with >= 45 measured minutes)
    valid_min = np.isfinite(st["ghi"])
    hours = sorted(set(hour_stamp[valid_min & (jh >= 0)]))
    H = []
    for h in hours:
        mm = (hour_stamp == h) & valid_min
        if mm.sum() < 45:
            continue
        row = dict(stamp=int(h), minutes=int(mm.sum()), measured=float(st["ghi"][mm].mean()),
                   clear_sky_min=float(cs_min[mm].mean()), km=float(tr["km"][mm].mean()),
                   lat=float(tr["lat"][mm].mean()), moving=float(tr["moving"][mm].mean()))
        for k, s in minute_series.items():
            v = s[mm]
            row[k] = float(np.nanmean(v)) if np.isfinite(v).mean() > 0.9 else np.nan
        H.append(row)
    stamp_idx = {r["stamp"]: i for i, r in enumerate(H)}
    # car-sensor persistence: measured clear-sky index L hours earlier x clear sky now
    for L in LEADS_H:
        k = f"car sensor +{L} h (measured stays)"
        for r in H:
            q = stamp_idx.get(r["stamp"] - 3600 * L)
            r[k] = np.nan
            if q is not None and H[q]["clear_sky_min"] > 50:
                kc = min(H[q]["measured"] / H[q]["clear_sky_min"], 1.3)
                r[k] = kc * r["clear_sky_min"]
    H = [r for r in H if r["clear_sky_min"] > 50]           # daytime only
    methods = [k for k in H[0] if k not in ("stamp", "minutes", "measured", "clear_sky_min",
                                             "km", "lat", "moving")]
    y_all = np.array([r["measured"] for r in H])

    # ---- write track and hourly data
    with open(os.path.join(a.out, "track_minutes.csv"), "w", newline="") as f:
        w = csv.writer(f)
        w.writerow(["time_utc", "time_acst", "km", "latitude", "longitude", "speed_kmh", "moving",
                    "ghi_measured_wm2", "clear_sky_wm2", "pressure_hpa"])
        for i in range(len(t)):
            if not np.isfinite(st["ghi"][i]):
                continue
            tu = dt.datetime.fromtimestamp(t[i], dt.timezone.utc)
            w.writerow([f"{tu:%Y-%m-%d %H:%M}", f"{tu + dt.timedelta(hours=LOCAL_OFFSET_H):%Y-%m-%d %H:%M}",
                        round(float(tr["km"][i]), 1), round(float(tr["lat"][i]), 4), round(float(tr["lon"][i]), 4),
                        round(float(np.nan_to_num(st["speed"][i]) * 3.6), 1), int(tr["moving"][i]),
                        round(float(st["ghi"][i]), 1), round(float(cs_min[i]), 1),
                        "" if not np.isfinite(st["pressure"][i]) else round(float(st["pressure"][i]), 1)])
    with open(os.path.join(a.out, "hourly.csv"), "w", newline="") as f:
        w = csv.writer(f)
        w.writerow(["time_utc_end_of_hour", "time_acst_end_of_hour", "km", "latitude", "minutes",
                    "share_moving", "measured", "clear_sky"] + methods)
        for r in H:
            tu = dt.datetime.fromtimestamp(r["stamp"], dt.timezone.utc)
            w.writerow([f"{tu:%Y-%m-%d %H:%M}", f"{tu + dt.timedelta(hours=LOCAL_OFFSET_H):%Y-%m-%d %H:%M}",
                        round(r["km"]), round(r["lat"], 3), r["minutes"], round(r["moving"], 2),
                        round(r["measured"], 1), round(r["clear_sky_min"], 1)]
                       + ["" if not np.isfinite(r[k]) else round(r[k], 1) for k in methods])

    # ---- metrics on the common sample (every method has a value)
    def sample(keys):
        ok = np.all([[np.isfinite(r[k]) for r in H] for k in keys], axis=0)
        return ok
    # main comparison: forecasts that need no earlier observation of the car (lead-time
    # persistence rows are compared separately per lead time below, so the main sample
    # keeps the first hours of each day)
    common_keys = [k for k in methods if " h (" not in k
                   and np.isfinite([r[k] for r in H]).mean() > 0.6]
    dropped = [k for k in methods if k not in common_keys]
    ok = sample(common_keys)
    y = y_all[ok]
    F = {k: np.array([r[k] for r in H])[ok] for k in common_keys}
    reps = block_bootstrap([F[k] - y for k in common_keys])
    ref = common_keys.index("satellite analysis")
    rows = []
    for i, k in enumerate(common_keys):
        s = stats(F[k], y)
        lo, hi = np.percentile(reps[:, i], [5, 95])
        dlo, dhi = np.percentile(reps[:, i] - reps[:, ref], [5, 95])
        rows.append(dict(method=k, **{x: round(v, 3) if isinstance(v, float) else v for x, v in s.items()},
                         mae_lo90=round(lo, 1), mae_hi90=round(hi, 1),
                         mae_minus_satellite_lo90=round(dlo, 1), mae_minus_satellite_hi90=round(dhi, 1)))
    rows.sort(key=lambda r: r["mae"])
    keys = list(rows[0])
    with open(os.path.join(a.out, "metrics.csv"), "w", newline="") as f:
        w = csv.DictWriter(f, fieldnames=keys); w.writeheader(); w.writerows(rows)

    # ---- per lead time: persistence vs. models, each row on its own common sample
    nwp_ref = [k for k in grids if k.endswith("short lead")]
    nwp_d1 = [k for k in grids if k.endswith("1 day ahead")]
    lead_rows = []
    for L in (0,) + LEADS_H:
        keys_L = (["satellite analysis"] if L == 0 else
                  [f"satellite +{L} h (cloudiness stays)", f"car sensor +{L} h (measured stays)"]) \
            + nwp_ref + nwp_d1 + ["always clear sky"]
        keys_L = [k for k in keys_L if k in methods]
        okL = sample(keys_L)
        for k in keys_L:
            s = stats(np.array([r[k] for r in H])[okL], y_all[okL])
            lead_rows.append(dict(lead_h=L, method=k, n=s["n"], mae=round(s["mae"], 1),
                                  bias=round(s["bias"], 1), rmse=round(s["rmse"], 1)))
    with open(os.path.join(a.out, "metrics_by_lead.csv"), "w", newline="") as f:
        w = csv.DictWriter(f, fieldnames=list(lead_rows[0])); w.writeheader(); w.writerows(lead_rows)

    # ---- energy per race day while driving (hours in the common sample)
    days = sorted(set((dt.datetime.fromtimestamp(r["stamp"] - 1800, dt.timezone.utc)
                       + dt.timedelta(hours=LOCAL_OFFSET_H)).date() for r in H))
    day_of = np.array([(dt.datetime.fromtimestamp(r["stamp"] - 1800, dt.timezone.utc)
                        + dt.timedelta(hours=LOCAL_OFFSET_H)).date() for r in H])[ok]
    drows = []
    for d in days:
        mk = day_of == d
        if not mk.any():
            continue
        row = dict(day=f"{d:%Y-%m-%d}", hours=int(mk.sum()), measured_kwh=round(y[mk].sum() / 1000, 2))
        for k in common_keys:
            row[k] = round(F[k][mk].sum() / 1000, 2)
        drows.append(row)
    with open(os.path.join(a.out, "daily_energy.csv"), "w", newline="") as f:
        w = csv.DictWriter(f, fieldnames=list(drows[0])); w.writeheader(); w.writerows(drows)
    dE = {k: np.array([r[k] - r["measured_kwh"] for r in drows]) for k in common_keys}

    # ---- track check: pressure altitude vs. terrain height along the route
    alt = 44330.0 * (1 - (st["pressure"] / 1013.25) ** (1 / 5.255))
    ter = np.interp(tr["km"], pkm[np.isfinite(elev)], elev[np.isfinite(elev)]) if np.isfinite(elev).any() \
        else np.full(len(t), np.nan)
    alt_corr = []
    for d in sorted(set(tr["local_day"])):
        mk = (tr["local_day"] == d) & tr["moving"] & np.isfinite(alt) & np.isfinite(ter)
        if mk.sum() > 60:
            off = np.median(alt[mk] - ter[mk])
            alt_corr.append((d, float(np.corrcoef(alt[mk], ter[mk])[0, 1]),
                             float(np.sqrt(np.mean((alt[mk] - off - ter[mk]) ** 2)))))

    # ---- sensor check: measured vs. satellite in clear hours
    clear = ok & np.array([r["satellite analysis"] / max(r["clear_sky_min"], 1) > 0.85
                           if np.isfinite(r["satellite analysis"]) else False for r in H])
    ratio = (np.median(y_all[clear] / np.array([r["satellite analysis"] for r in H])[clear])
             if clear.any() else float("nan"))

    # ---- figures
    figures(a.out, st, tr, km, H, ok, common_keys, rows, lead_rows, cs_min, alt, ter, grids)

    # ---- summary
    pretty = lambda k: k.replace(" | ", ", ")
    md = ["# Sunshine forecasts vs. our own car's weather station (WSC 2025)", "",
          f"Race days {days[0]:%d.%m.} – {days[-1]:%d.%m.%Y}; {len(y)} daytime hours along the car's track "
          f"(common sample: every method has a value). Truth = the car's own GHI sensor, hourly means. "
          f"Forecast values are taken at the place where the car was in each minute.", ""]
    if dropped:
        md += ["Left out of the common sample (too few values): " + ", ".join(pretty(k) for k in dropped), ""]
    md += ["## Error per method (hourly, W/m²)",
           "| Method | Mean error | 90 % interval | vs. satellite (90 %) | Bias | Correlation |",
           "|---|---|---|---|---|---|"]
    for r in rows:
        md.append(f"| {pretty(r['method'])} | {r['mae']:.0f} | {r['mae_lo90']:.0f}–{r['mae_hi90']:.0f} | "
                  f"{r['mae_minus_satellite_lo90']:+.0f} … {r['mae_minus_satellite_hi90']:+.0f} | "
                  f"{r['bias']:+.0f} | {r['corr']:.2f} |")
    md += ["", "'vs. satellite': how much larger the mean error is than that of the satellite analysis. "
           "If the interval does not include 0, the difference is clear; otherwise the data cannot "
           "separate the two.", "",
           "## Error vs. lead time (W/m²)", "",
           "Each column on its own common sample (hours where every method has a value). "
           "Persistence needs an observation L hours earlier, so later leads lose the first hours "
           "of each day.", "",
           "| Method | " + " | ".join(("0 h" if L == 0 else f"+{L} h") for L in (0,) + LEADS_H) + " |",
           "|---|" + "---|" * (len(LEADS_H) + 1)]
    lead_tab = collections.OrderedDict()
    for r in lead_rows:
        name = r["method"]
        for pre in ("satellite +", "car sensor +"):
            if name.startswith(pre):
                name = ("satellite (0 h: analysis, later: cloudiness stays)" if pre == "satellite +"
                        else "car sensor (measured stays)")
        if name == "satellite analysis":
            name = "satellite (0 h: analysis, later: cloudiness stays)"
        lead_tab.setdefault(name, {})[r["lead_h"]] = r
    for name, d in lead_tab.items():
        md.append(f"| {pretty(name)} | " + " | ".join(f"{d[L]['mae']:.0f}" if L in d else "–"
                                                     for L in (0,) + LEADS_H) + " |")
    md.append("| hours in sample | " + " | ".join(
        str(next((r["n"] for r in lead_rows if r["lead_h"] == L), 0)) for L in (0,) + LEADS_H) + " |")
    md += ["", "## Energy while driving, per race day (kWh/m²)", "",
           "| Method | " + " | ".join(r["day"][5:] for r in drows) + " | mean abs. error |",
           "|---|" + "---|" * (len(drows) + 1),
           "| **measured (car)** | " + " | ".join(f"{r['measured_kwh']:.2f}" for r in drows) + " | |"]
    for k in sorted(common_keys, key=lambda k: np.abs(dE[k]).mean()):
        md.append(f"| {pretty(k)} | " + " | ".join(f"{r[k]:.2f}" for r in drows)
                  + f" | {np.abs(dE[k]).mean():.2f} |")
    md.append("| hours | " + " | ".join(str(r["hours"]) for r in drows) + " | |")
    md += ["", "## Checks", "",
           f"* Sensor vs. satellite in clear hours (satellite ≥ 85 % of clear sky): measured / satellite = "
           f"{ratio:.2f} (median). Far from 1 means the sensor (tilt, calibration, shading) or the "
           f"satellite has a systematic offset – look at the bias column above.",
           "* Track: car speed from true − apparent wind, anchored to the control stops "
           "(route km / reconstructed km per leg: " + ", ".join(f"{s:.2f}" for s in tr["seg_scale"]) + ").",
           "* Pressure altitude vs. terrain height at the reconstructed position, per day "
           "(correlation, RMS m): " + "; ".join(f"{d:%d.%m.} {c:.2f} / {e:.0f} m" for d, c, e in alt_corr),
           "* Control stops found: " + ", ".join(
               f"{an['name']} {(dt.datetime.fromtimestamp(st['t'][an['i']], dt.timezone.utc) + dt.timedelta(hours=LOCAL_OFFSET_H)):%d.%m. %H:%M}"
               for an in tr["anchors"]),
           "", "## How to read this", "",
           "* 'satellite analysis' is the best our own satellite chain can do (lead time 0). Our nowcast "
           "for +1…3 h scores about like 'satellite +L h (cloudiness stays)' (see result_*/skill_curve.csv), "
           "so those rows stand for our own forecast at that lead time.",
           "* 'short lead' model values join the first hours of each model run and flatter the models; "
           "'1 day ahead' is what the team could see the evening before.",
           "* Only 5 race days, mostly clear weather: the intervals are wide. Add more days (e.g. test "
           "drives) to the station_data folder to narrow them."]
    open(os.path.join(a.out, "summary.md"), "w", encoding="utf-8").write("\n".join(md) + "\n")
    print("\n" + "\n".join(md))
    print(f"\n→ results in {a.out}/")


def read_first_time(path):
    try:
        m = read_jsonl(path)
        return float(m["global_radiation"][0][0])
    except Exception:
        return 0.0


# ---------------------------------------------------------------- figures
def figures(out, st, tr, km, H, ok, common_keys, rows, lead_rows, cs_min, alt, ter, grids):
    try:
        import matplotlib
        matplotlib.use("Agg")
        import matplotlib.pyplot as plt
        import matplotlib.dates as mdates
    except ImportError:
        print("matplotlib not installed – no figures")
        return
    INK, INK2, GRID = "#0b0b0b", "#52514e", "#e4e3df"
    SER = ["#2a78d6", "#eb6834", "#1baf7a", "#eda100", "#e87ba4", "#008300", "#6250d6", "#e34948"]
    plt.rcParams.update({"axes.edgecolor": INK2, "axes.labelcolor": INK, "xtick.color": INK2,
                         "ytick.color": INK2, "axes.grid": True, "grid.color": GRID, "grid.linewidth": 0.6,
                         "axes.spines.top": False, "axes.spines.right": False, "font.size": 9,
                         "lines.linewidth": 2})
    tl = np.array([dt.datetime.fromtimestamp(x, dt.timezone.utc).replace(tzinfo=None)
                   + dt.timedelta(hours=LOCAL_OFFSET_H) for x in st["t"]])
    days = [d for d in sorted(set(tr["local_day"])) if tr["moving"][tr["local_day"] == d].sum() > 30]

    # track and altitude check
    fig, axs = plt.subplots(2, len(days), figsize=(3.6 * len(days), 6.5), squeeze=False)
    for c, d in enumerate(days):
        mk = tr["local_day"] == d
        ax = axs[0, c]
        ax.plot(tl[mk], tr["km"][mk], color=SER[0])
        for an in tr["anchors"]:
            if tr["local_day"][an["i"]] == d:
                ax.plot(tl[an["i"]], an["km"], "o", ms=8, color=INK, mec="white", mew=2)
                ax.annotate(an["name"], (tl[an["i"]], an["km"]), xytext=(5, -12),
                            textcoords="offset points", fontsize=8, color=INK)
        ax.set_title(f"{d:%a %d.%m.}", color=INK)
        ax.xaxis.set_major_formatter(mdates.DateFormatter("%H"))
        if c == 0:
            ax.set_ylabel("km along the route")
        ax = axs[1, c]
        mv = mk & tr["moving"]
        off = np.nanmedian((alt - ter)[mv]) if mv.any() else 0
        ax.plot(tl[mk], ter[mk], color=INK2, lw=1.5, label="terrain at reconstructed position")
        ax.plot(tl[mk], alt[mk] - off, color=SER[1], lw=1.5, label="station pressure altitude (shifted)")
        ax.xaxis.set_major_formatter(mdates.DateFormatter("%H"))
        ax.set_xlabel("local time (ACST)")
        if c == 0:
            ax.set_ylabel("height [m]")
            ax.legend(fontsize=7, loc="upper left", frameon=False)
    fig.suptitle("Reconstructed car position (dots = control stops) and altitude check", color=INK)
    fig.tight_layout()
    fig.savefig(os.path.join(out, "fig_track.png"), dpi=110)
    plt.close(fig)

    # per day: measured vs. satellite vs. ECMWF
    model = next((k.split(" | ")[0] for k in common_keys if "ecmwf" in k), None) or \
        next((k.split(" | ")[0] for k in common_keys if " | " in k), None)
    show = [("satellite analysis", SER[0], "-"),
            (f"{model} | short lead", SER[1], "-"), (f"{model} | 1 day ahead", SER[1], "--")]
    hl = np.array([dt.datetime.fromtimestamp(r["stamp"] - 1800, dt.timezone.utc).replace(tzinfo=None)
                   + dt.timedelta(hours=LOCAL_OFFSET_H) for r in H])
    fig, axs = plt.subplots(1, len(days), figsize=(3.8 * len(days), 4.2), sharey=True, squeeze=False)
    for c, d in enumerate(days):
        ax = axs[0, c]
        mk = (tr["local_day"] == d) & np.isfinite(st["ghi"])
        ax.fill_between(tl[mk], cs_min[mk], color="#f0efec", label="clear sky", lw=0)
        ax.plot(tl[mk], st["ghi"][mk], color="#b9b8b2", lw=0.6, label="car sensor, 1 min")
        hm = np.array([x.date() == d for x in hl])
        ax.plot(hl[hm], [r["measured"] for r, m in zip(H, hm) if m], "o-", color=INK, ms=4,
                lw=2, label="car sensor, hourly mean")
        for k, col, ls in show:
            if k in H[0]:
                ax.plot(hl[hm], [r[k] for r, m in zip(H, hm) if m], ls, color=col, lw=2,
                        label=k.replace(" | ", ", "))
        ax.set_title(f"{d:%a %d.%m.}", color=INK)
        ax.xaxis.set_major_formatter(mdates.DateFormatter("%H"))
        ax.set_xlabel("local time (ACST)")
    axs[0, 0].set_ylabel("GHI [W/m²]")
    axs[0, 0].legend(fontsize=7, loc="upper left", frameon=False)
    fig.suptitle("Measured on the car vs. satellite and weather model along the car's track", color=INK)
    fig.tight_layout()
    fig.savefig(os.path.join(out, "fig_days.png"), dpi=110)
    plt.close(fig)

    # errors per method with 90 % interval
    def group_color(k):
        if k.startswith("satellite"):
            return SER[0]
        if k.startswith("car sensor"):
            return SER[6]
        if k.startswith("always"):
            return INK2
        return SER[1]
    fig, ax = plt.subplots(figsize=(9, 0.32 * len(rows) + 1.4))
    yv = np.arange(len(rows))
    ax.barh(yv, [r["mae"] for r in rows], height=0.6, color=[group_color(r["method"]) for r in rows])
    ax.errorbar([r["mae"] for r in rows], yv,
                xerr=[[r["mae"] - r["mae_lo90"] for r in rows], [r["mae_hi90"] - r["mae"] for r in rows]],
                fmt="none", ecolor=INK, elinewidth=1, capsize=3)
    for i, r in enumerate(rows):
        ax.text(r["mae_hi90"] + 2, i, f"{r['mae']:.0f}", va="center", fontsize=8, color=INK)
    ax.set_yticks(yv, [r["method"].replace(" | ", ", ") for r in rows])
    ax.invert_yaxis()
    ax.grid(axis="y", visible=False)
    ax.set_xlabel("mean error per daytime hour vs. car sensor [W/m²]  (lower = better, bar = 90 % interval)")
    ax.set_title("Which forecast is closest to what the car measured? (same hours for every method)",
                 color=INK, loc="left")
    fig.tight_layout()
    fig.savefig(os.path.join(out, "fig_errors.png"), dpi=110)
    plt.close(fig)

    # error vs. lead time
    fig, ax = plt.subplots(figsize=(8, 4.5))
    by = collections.defaultdict(dict)
    for r in lead_rows:
        by[r["method"]][r["lead_h"]] = r["mae"]
    sat = [by["satellite analysis"].get(0)] + [by.get(f"satellite +{L} h (cloudiness stays)", {}).get(L)
                                               for L in LEADS_H]
    car = [by.get(f"car sensor +{L} h (measured stays)", {}).get(L) for L in LEADS_H]
    ax.plot((0,) + LEADS_H, sat, "o-", color=SER[0], ms=8, label="satellite (our analysis / persistence)")
    ax.plot(LEADS_H, car, "s-", color=SER[6], ms=8, label="car sensor persistence")
    nm = [k for k in by if k.endswith("short lead") or k.endswith("1 day ahead")]
    for i, k in enumerate(nm):
        vals = [by[k].get(L) for L in (0,) + LEADS_H]
        ls = "--" if "1 day" in k else "-"
        col = SER[[1, 2, 3, 4][i // 2 % 4]] if len(nm) > 1 else SER[1]
        ax.plot((0,) + LEADS_H, vals, ls, color=col, lw=1.5, label=k.replace(" | ", ", "))
    ax.set_xticks((0,) + LEADS_H, ["0 h\n(analysis)"] + [f"+{L} h" for L in LEADS_H])
    ax.set_ylabel("mean error [W/m²]")
    ax.set_title("Error vs. lead time (models do not depend on it – lines show their level "
                 "on the same hours)", color=INK, loc="left", fontsize=9)
    ax.legend(fontsize=7, frameon=False, loc="upper left", ncol=2)
    ax.set_ylim(bottom=0)
    fig.tight_layout()
    fig.savefig(os.path.join(out, "fig_lead.png"), dpi=110)
    plt.close(fig)


if __name__ == "__main__":
    main()
