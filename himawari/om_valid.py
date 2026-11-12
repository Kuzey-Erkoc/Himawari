#!/usr/bin/env python3
"""
How good are Open-Meteo's sunshine forecasts along our race route?

For the race season (default: 15 Aug – 5 Sep 2025 and 2026, i.e. the race week
22–29 Aug ± one week) the script fetches, for points along route_locations.py:

  forecasts   hourly global irradiance (GHI, W/m²) of several weather models
                * "short lead": Historical Forecast API – Open-Meteo joins the first hours of
                  every model run, i.e. what the newest run said a few hours ahead
                * "1 day ahead" / "2 days ahead": Previous Runs API – what the model said
                  24 h / 48 h before (what the team would have seen the evening before)
  actual      hourly GHI derived from the Himawari satellite (Open-Meteo Satellite Radiation
              API) – a satellite estimate, not a ground measurement

and compares them hour by hour and day by day (daily energy in kWh/m², which is what the
race strategy cares about). Two simple references show whether a forecast is worth anything:
  * "always clear sky"   – no clouds at all (computed from the sun position)
  * "same as yesterday"  – the actual values of the day before (for the day-ahead forecasts)

Output (folder validation_openmeteo/):
  data_hourly.csv          all hourly values (point, time, actual, every forecast)
  metrics_hourly.csv       bias, MAE, RMSE, correlation per forecast (daytime hours)
  metrics_daily.csv        error of the daily energy per forecast (complete days only),
                           also split north / south / race week / change days / steady days
  metrics_by_region.csv    hourly metrics split north / south of 23.4°S, race week,
                           weather-change days and steady days
  bias_factors.csv         bias-correction factor per model and lead (learned on the other period)
  summary.md               the main results in words and tables
  fig_timeseries.png       actual vs. forecasts at Darwin, Alice Springs, Adelaide (race week)
  fig_errors.png           mean error and bias per forecast and lead time
  fig_daily_energy.png     daily energy: forecast vs. actual
  fig_hour_of_day.png      bias by hour of day (does the model miss mornings/afternoons?)
  fig_regions_changes.png  error north vs. south, and on steady vs. weather-change days

Bias correction: all models turned out systematically too dark. Every forecast is therefore
also evaluated multiplied by one factor per model and lead ("… | corrected"). The factor is
learned on the other period only (2025 → applied to 2026, 2026 → 2025), so it is a fair test.
Weather-change days: days whose daily energy differs by more than 1 kWh/m² from the day before
(--change-kwh). On these days "same as yesterday" is wrong by definition.

Usage:
    pip install numpy requests matplotlib
    python om_valid.py                       # race seasons 2025 and 2026
    python om_valid.py --period 2025-08-22 2025-08-29   # one period
    python om_valid.py --every 3             # more route points (default every 5th)
    python om_valid.py --models best_match,ecmwf_ifs025,gfs_global,icon_global
Send the folder validation_openmeteo/ (or at least the CSV files and summary.md) back for
the evaluation.
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

import numpy as np

HIST_URL = "https://historical-forecast-api.open-meteo.com/v1/forecast"
PREV_URL = "https://previous-runs-api.open-meteo.com/v1/forecast"
SAT_URL = "https://satellite-api.open-meteo.com/v1/archive"
DEFAULT_MODELS = "best_match,ecmwf_ifs025,gfs_global,icon_global"
DEFAULT_PERIODS = [("2025-08-15", "2025-09-05"), ("2026-08-15", "2026-09-05")]
RACE_WEEKS = [("2025-08-22", "2025-08-29"), ("2026-08-22", "2026-08-29")]
NORTH_LAT = -23.44          # Tropic of Capricorn
LOCAL_OFFSET_H = 9.5        # ACST: NT and SA in August (no daylight saving yet)
BATCH = 50
TIMEOUT = 60
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
        if get(batches[0][:1]) is not None:          # lists refused -> one point per request
            with cf.ThreadPoolExecutor(4) as ex:
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
                   help="date range YYYY-MM-DD YYYY-MM-DD (repeatable); default race seasons 2025+2026")
    p.add_argument("--models", default=DEFAULT_MODELS)
    p.add_argument("--every", type=int, default=5, help="use every N-th route point (default 5)")
    p.add_argument("--change-kwh", type=float, default=1.0,
                   help="a day counts as weather change if its energy differs from the day before "
                        "by more than this (kWh/m², default 1.0)")
    p.add_argument("--out", default=os.path.join(here, "validation_openmeteo"))
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
    periods = a.period or DEFAULT_PERIODS
    os.makedirs(a.out, exist_ok=True)
    print(f"{len(pts)} route points (every {a.every}th, {pkm[-1]:.0f} km), models: {', '.join(models)}")

    # hourly time axis (stamp = end of the averaging hour, as Open-Meteo labels it)
    times = []
    for s, e in periods:
        t = dt.datetime.strptime(s, "%Y-%m-%d")
        t_end = dt.datetime.strptime(e, "%Y-%m-%d") + dt.timedelta(days=1)
        while t < t_end:
            times.append(t)
            t += dt.timedelta(hours=1)
    times = sorted(set(times))
    n, T = len(pts), len(times)

    series = collections.OrderedDict()      # name -> (n × T) array
    for s, e in periods:
        print(f"\nPeriod {s} … {e}")
        base = dict(start_date=s, end_date=e, timezone="GMT")
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
    day = cs > 50                       # daytime hours (sun clearly up)

    # ---- bias correction: the models are systematically too dark. One scale factor per
    #      model and lead time, always learned on the OTHER period (fitted on 2025 → applied
    #      to 2026 and vice versa), so the corrected scores are a fair test, not hindsight.
    fold = np.zeros(T, int)
    plabel = []
    for pi, (s, e) in enumerate(periods):
        t0 = dt.datetime.strptime(s, "%Y-%m-%d")
        t1 = dt.datetime.strptime(e, "%Y-%m-%d") + dt.timedelta(days=1)
        fold[np.array([t0 <= t < t1 for t in times])] = pi
        plabel.append(f"{s} – {e}")
    if len(periods) == 1:               # only one period: learn on one half, test on the other
        fold = (np.arange(T) >= T // 2).astype(int)
        plabel = ["first half", "second half"]
    raw = [k for k in series if " | " in k]
    fac_rows = []
    for k in raw:
        f = series[k]
        out = np.full_like(f, np.nan)
        for fo in range(int(fold.max()) + 1):
            train = day & np.isfinite(y) & np.isfinite(f) & (fold != fo)[None, :]
            if train.sum() < 50:
                continue
            fac = float(y[train].sum() / max(f[train].sum(), 1e-9))
            test = np.broadcast_to((fold == fo)[None, :], f.shape)
            out[test] = f[test] * fac
            fac_rows.append(dict(forecast=k, applied_to=plabel[fo], factor=round(fac, 3),
                                 learned_from_hours=int(train.sum())))
        if np.isfinite(out).any():
            series[f"{k} | corrected"] = out
    order = (["actual"] + [x for k in raw for x in (k, f"{k} | corrected")]
             + ["always clear sky", "same as yesterday"])
    series = collections.OrderedDict((k, series[k]) for k in order if k in series)
    with open(os.path.join(a.out, "bias_factors.csv"), "w", newline="") as f:
        w = csv.DictWriter(f, fieldnames=["forecast", "applied_to", "factor", "learned_from_hours"])
        w.writeheader()
        w.writerows(fac_rows)
    if fac_rows:
        print("\nBias correction (forecast × factor, learned on the other period):")
        for r in fac_rows:
            print(f"  {r['forecast']:<32} → {r['applied_to']}: × {r['factor']:.3f}")

    # ---- hourly data file
    with open(os.path.join(a.out, "data_hourly.csv"), "w", newline="") as f:
        w = csv.writer(f)
        w.writerow(["time_utc_end_of_hour", "point", "km", "latitude", "longitude", "region",
                    "clear_sky_wm2"] + list(series))
        for j, t in enumerate(times):
            if not day[:, j].any():
                continue
            for i in range(n):
                w.writerow([t.strftime("%Y-%m-%d %H:%M"), pts[i][2], round(float(pkm[i]), 1),
                            pts[i][0], pts[i][1], "north" if north[i] else "south",
                            round(float(cs[i, j]), 1)]
                           + ["" if not np.isfinite(v[i, j]) else round(float(v[i, j]), 1)
                              for v in series.values()])

    # ---- local days (ACST) and daily energy (kWh/m²). Only days that lie completely inside
    #      the fetched time range count (all 24 hours present) – the half days at the edges of
    #      each period would otherwise show up as nonsense values of ~0.2 kWh/m².
    local_day = [(t - dt.timedelta(hours=1) + dt.timedelta(hours=LOCAL_OFFSET_H)).date() for t in times]
    dpos = collections.defaultdict(list)
    for j, d in enumerate(local_day):
        dpos[d].append(j)
    days = sorted(d for d in dpos if len(dpos[d]) == 24)
    dindex = {d: di for di, d in enumerate(days)}
    jd = np.array([dindex.get(d, -1) for d in local_day])

    def daily(v):
        out = np.full((n, len(days)), np.nan)
        for di, d in enumerate(days):
            js = dpos[d]
            need = cs[:, js] > 1
            vv = v[:, js]
            complete = np.all(np.isfinite(vv) | ~need, axis=1) & need.any(axis=1)
            out[:, di] = np.where(complete, np.nansum(np.where(need, vv, 0), axis=1) / 1000.0, np.nan)
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
    d_all = np.all([np.isfinite(v) for v in dser.values()], axis=0) & np.isfinite(dy)
    dy = np.where(d_all, dy, np.nan)
    dser = {k: np.where(d_all, v, np.nan) for k, v in dser.items()}
    change &= d_all
    steady &= d_all
    print(f"\n{len(days)} complete days; point-days with a weather change "
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
        t0 = dt.datetime.strptime(s, "%Y-%m-%d")
        t1 = dt.datetime.strptime(e, "%Y-%m-%d") + dt.timedelta(days=1)
        race |= np.array([t0 <= t < t1 for t in times])
        race_d |= np.array([t0.date() <= d < t1.date() for d in days], bool)
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

    # ---- figures
    try:
        import matplotlib
        matplotlib.use("Agg")
        import matplotlib.pyplot as plt
        import matplotlib.dates as mdates
        cols = plt.rcParams["axes.prop_cycle"].by_key()["color"]
        mcol = {m: cols[i % len(cols)] for i, m in enumerate(models)}
        lstyle = {"short lead": "-", "1 day ahead": "--", "2 days ahead": ":"}

        # time series at three places
        places = [("Darwin", -12.46, 130.84), ("Alice Springs", -23.70, 133.88), ("Adelaide", -34.93, 138.60)]
        best = models[0]
        nper = len(periods)
        fig, axs = plt.subplots(3, nper, figsize=(7 * nper, 9), squeeze=False)
        for r_i, (pname, pla, plo) in enumerate(places):
            i = int(np.argmin(haversine_km(lat, lon, pla, plo)))
            for c_i, (s, e) in enumerate(periods):
                ax = axs[r_i, c_i]
                t0 = dt.datetime.strptime(s, "%Y-%m-%d")
                t1 = dt.datetime.strptime(e, "%Y-%m-%d") + dt.timedelta(days=1)
                # show the race week if it lies in the period, else the first 8 days
                rw = [(dt.datetime.strptime(a_, "%Y-%m-%d"),
                       dt.datetime.strptime(b_, "%Y-%m-%d") + dt.timedelta(days=1))
                      for a_, b_ in RACE_WEEKS]
                rw = [(max(t0, x), min(t1, y_)) for x, y_ in rw if x < t1 and y_ > t0]
                t0, t1 = rw[0] if rw else (t0, min(t1, t0 + dt.timedelta(days=8)))
                js = [j for j, t in enumerate(times) if t0 <= t < t1]
                tt = [times[j] for j in js]
                ax.fill_between(tt, cs[i, js], color="0.9", label="clear sky")
                ax.plot(tt, y[i, js], color="k", lw=1.6, label="actual (satellite)")
                for lab in ("short lead", "1 day ahead"):
                    k = f"{best} | {lab}"
                    if k in series:
                        ax.plot(tt, series[k][i, js], lstyle[lab], color=mcol[best], lw=1.1,
                                label=f"{best}, {lab}")
                ax.set_title(f"{pname} ({pts[i][2]}, km {pkm[i]:.0f}) – {t0:%d.%m.}–"
                             f"{t1 - dt.timedelta(days=1):%d.%m.%Y}")
                ax.set_ylabel("GHI [W/m²]")
                ax.xaxis.set_major_formatter(mdates.DateFormatter("%d.%m"))
                if r_i == 0 and c_i == 0:
                    ax.legend(fontsize=8, loc="upper right")
        fig.suptitle("Sunshine along the route: actual (satellite) vs. weather-model forecast (UTC)")
        fig.tight_layout()
        fig.savefig(os.path.join(a.out, "fig_timeseries.png"), dpi=110)
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
        fig.suptitle("Open-Meteo forecasts vs. satellite-derived actual values (all periods, route points)")
        fig.tight_layout()
        fig.savefig(os.path.join(a.out, "fig_errors.png"), dpi=110)
        plt.close(fig)

        # daily energy scatter
        fig, axs = plt.subplots(1, 3, figsize=(15, 5), sharex=True, sharey=True)
        for ax, lab in zip(axs, ("short lead", "1 day ahead", "2 days ahead")):
            for m in models:
                k = f"{m} | {lab}"
                if k in dser:
                    ax.scatter(dy.ravel(), dser[k].ravel(), s=8, alpha=0.5, color=mcol[m], label=m)
            lim = np.nanmax(dy) * 1.1 if np.isfinite(dy).any() else 10
            ax.plot([0, lim], [0, lim], "k-", lw=0.8)
            ax.set_xlabel("actual daily energy [kWh/m²]")
            ax.set_title(lab)
        axs[0].set_ylabel("forecast daily energy [kWh/m²]")
        axs[0].legend(fontsize=8)
        fig.suptitle("Daily energy per route point and day (on the line = perfect)")
        fig.tight_layout()
        fig.savefig(os.path.join(a.out, "fig_daily_energy.png"), dpi=110)
        plt.close(fig)

        # bias by local hour of day
        lh = np.array([((t - dt.timedelta(minutes=30)).hour + LOCAL_OFFSET_H + (t - dt.timedelta(minutes=30)).minute / 60) % 24
                       for t in times])
        bins = np.arange(6, 19)
        fig, ax = plt.subplots(figsize=(10, 5))
        for k in fc_names:
            if "|" not in k or "2 days" in k or k.endswith("corrected"):
                continue
            m, lab = k.split(" | ")
            b = []
            for h0 in bins:
                msk = common & ((lh >= h0) & (lh < h0 + 1))[None, :] & np.isfinite(series[k])
                b.append(np.nanmean((series[k] - y)[msk]) if msk.any() else np.nan)
            ax.plot(bins + 0.5, b, lstyle[lab], color=mcol[m], marker="o", ms=3, label=k)
        ax.axhline(0, color="k", lw=0.8)
        ax.set_xlabel("local time (ACST)"); ax.set_ylabel("bias forecast − actual [W/m²]")
        ax.set_title("Systematic error by time of day (above 0 = forecast too sunny)")
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
          f"Periods: {', '.join(f'{s} – {e}' for s, e in periods)}; {n} route points; "
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
    if fac_rows:
        md += ["## Bias correction factors (forecast × factor; learned on the other period)",
               "| Forecast | applied to | factor |", "|---|---|---|"]
        md += [f"| {r['forecast'].replace(' | ', ', ')} | {r['applied_to']} | {r['factor']:.3f} |"
               for r in fac_rows]
        md.append("")
    md += ["Bias > 0: forecast too sunny. 'always clear sky' and 'same as yesterday' are the "
           "references a forecast has to beat. 'corrected' = forecast multiplied by a factor "
           "learned on the other period (fair test)."]
    open(os.path.join(a.out, "summary.md"), "w", encoding="utf-8").write("\n".join(md))
    print("\n" + "\n".join(md))
    print(f"\n→ results in {a.out}/")


if __name__ == "__main__":
    main()
