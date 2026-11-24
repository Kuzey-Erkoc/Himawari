#!/usr/bin/env python3
"""
Sunshine forecasts of the weather models vs. what actually arrived, along our race route.

For a period (default: the race week 22–29 Aug 2026, local dates) the script fetches, for
every point of route_locations.py (--every N: every N-th point):

  forecasts   hourly global irradiance (GHI, W/m²) of the weather models via Open-Meteo
              (default ECMWF and GFS; others with --models)
                * "short lead": the newest forecast for each hour, made a few hours before
                * "1 day ahead" … "5 days ahead": what the model said 1 … 5 days before
  actual      hourly GHI derived from the Himawari satellite (Open-Meteo Satellite Radiation
              API) – a satellite estimate, not a ground measurement

No averages are computed. Besides the values themselves the script counts cloud events: an
hour is cloudy if less than 60 % of the clear-sky sunshine arrives, a day if its energy is below
80 % of a clear day. For every forecast it counts how many cloudy hours/days it caught (hit),
missed, and how often it forecast clouds that did not come (false alarm).
The models forecast systematically too little sunshine, so each forecast is also written
multiplied by the fixed factor of main.py / trip_forecast.py ("… | corrected", ECMWF × 1.06,
GFS × 1.03; off with --no-bias-correction).

Output (folder validation_openmeteo_<start date>/):
  data_hourly.csv      all hourly values per route point (time UTC and local ACST, clear sky,
                       actual, every forecast) – daytime hours only
  cloud_events.csv     counts per model and lead time: cloud hours/days that were forecast (hit),
                       missed, false alarms – hourly, hourly with ±1 h tolerance, and per day
  fig_cloud_events.png the same for ECMWF as bar chart (counts, no averages)
  fig_timeseries.png   ECMWF × 1.06 only: newest forecast, 1 day and 2 days before in one figure,
                       actual vs. forecast hour by hour at the four places
  fig_timeseries_short_lead.png, fig_timeseries_1day.png … fig_timeseries_5days.png
                       one figure per lead time: actual vs. forecast hour by hour at Darwin,
                       Alice Springs, Coober Pedy and Adelaide

Open-Meteo gives hourly values as the average over the hour before the time stamp; that is how
the data come and not something this script computes.

Usage:
    pip install numpy requests matplotlib
    python om_valid.py                                    # race week 22–29 Aug 2026, all points
    python om_valid.py --period 2025-08-22 2025-08-29     # race week 2025
    python om_valid.py --models ecmwf_ifs025,gfs_global,icon_global
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
DEFAULT_MODELS = "ecmwf_ifs025,gfs_global"
CLOUDY_HOUR = 0.6          # an hour is "cloudy" if less than 60 % of the clear-sky sunshine arrives
CLOUDY_DAY = 0.8           # a day is "cloudy" if its energy is below 80 % of a clear day
MAX_LEAD_DAYS = 5           # forecasts from 1 … 5 days before (Open-Meteo Previous Runs API)
LEADS = [(0, "short lead")] + [(d, f"{d} day{'s' if d > 1 else ''} ahead") for d in range(1, MAX_LEAD_DAYS + 1)]
DEFAULT_PERIODS = [("2026-08-22", "2026-08-29")]       # race week 2026, local dates (ACST)
# fixed bias-correction factors – the same as in main.py / trip_forecast.py (from the 2025+2026
# seasons and the car's weather station).
BIAS_FACTOR = {"ecmwf_ifs025": 1.06, "gfs_global": 1.03, "icon_global": 1.09, "best_match": 1.06}
MIN_COVERAGE = 0.9          # a forecast needs values for 90 % of the daytime hours, else left out
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
    p.add_argument("--no-bias-correction", action="store_true",
                   help="do not add the corrected forecasts (× BIAS_FACTOR)")
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
    # are LOCAL days (ACST): a local day runs from 14:30 UTC the day before.
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
            print(f"  {m}: 1 to {MAX_LEAD_DAYS} days ahead …")
            vars_ = ",".join(f"shortwave_radiation_previous_day{d}" for d in range(1, MAX_LEAD_DAYS + 1))
            pr = fetch(PREV_URL, pts, dict(base, hourly=vars_, models=m), f"{m} previous runs")
            for d, lab in LEADS[1:]:
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

    # clear sky (upper limit), computed from the sun position
    cs = np.stack([clear_sky_hour(lat, lon, t) for t in times], axis=1)

    # keep only the hours of the requested local days
    keep = np.array([any(d0 <= local_date(t) <= d1 for d0, d1 in pdates) for t in times])
    times = [t for t, k_ in zip(times, keep) if k_]
    T = len(times)
    series = collections.OrderedDict((k, v[:, keep]) for k, v in series.items())
    cs = cs[:, keep]
    y = series["actual"]
    day = cs > 50                       # daytime hours (sun clearly up)


    # a forecast that is missing for many hours is left out (with a note)
    base_h = day & np.isfinite(y)
    for k in [k for k in series if " | " in k]:
        cov = float(np.isfinite(series[k])[base_h].sum() / max(base_h.sum(), 1))
        if cov < MIN_COVERAGE:
            print(f"  ({k}: values for only {100 * cov:.0f} % of the daytime hours – left out)")
            del series[k]
    # broken archive values: far above what is physically possible (clear sky) in many hours
    for k in [k for k in series if " | " in k]:
        f = series[k]
        ok = day & np.isfinite(f)
        too_high = (f > 1.3 * cs + 50) & ok
        if ok.any() and too_high.sum() > 0.2 * ok.sum():
            print(f"  ({k}: values far above clear sky in many hours – implausible, left out)")
            del series[k]
    models = [m for m in models if any(k.startswith(f"{m} | ") for k in series)]

    # bias correction with the fixed factors of main.py / trip_forecast.py
    raw = [k for k in series if " | " in k]
    if not a.no_bias_correction:
        for k in raw:
            fac = BIAS_FACTOR.get(k.split(" | ")[0])
            if fac:
                series[f"{k} | corrected"] = series[k] * fac
        used = sorted({k.split(" | ")[0] for k in raw if BIAS_FACTOR.get(k.split(" | ")[0])})
        if used:
            print("\nCorrected forecasts added: " + ", ".join(f"{m} × {BIAS_FACTOR[m]:.2f}" for m in used))
    order = ["actual"] + [x for k in raw for x in (k, f"{k} | corrected")]
    series = collections.OrderedDict((k, series[k]) for k in order if k in series)

    # ---- hourly data file (daytime hours)
    with open(os.path.join(a.out, "data_hourly.csv"), "w", newline="") as f:
        w = csv.writer(f)
        w.writerow(["time_utc_end_of_hour", "time_local_end_of_hour", "point", "km", "latitude",
                    "longitude", "region", "clear_sky_wm2"] + list(series))
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

    # ---- cloud events: counts, no averages. For every forecast and lead time:
    #      hit = clouds forecast and clouds came, miss = clouds came but not forecast,
    #      false alarm = clouds forecast but did not come, correct clear = neither.
    #      Hourly (exact hour, and with ±1 h tolerance for timing) and per whole day.
    strong = (cs > 200) & np.isfinite(y)          # hours with enough sun to judge clouds
    act_h = strong & (y < CLOUDY_HOUR * cs)

    def widen(m):                                  # cloudy within ±1 h at the same point
        out = m.copy()
        out[:, 1:] |= m[:, :-1]
        out[:, :-1] |= m[:, 1:]
        return out
    local_day = [local_date(t) for t in times]
    days = sorted(set(local_day))
    dpos = {d: [j for j, ld in enumerate(local_day) if ld == d] for d in days}

    def day_energy(v):                             # kWh/m² per point and local day (sum, daytime)
        out = np.full((n, len(days)), np.nan)
        for di, d in enumerate(days):
            js = dpos[d]
            need = day[:, js]
            vv = v[:, js]
            ok = np.all(np.isfinite(vv) | ~need, axis=1) & need.any(axis=1)
            out[:, di] = np.where(ok, np.where(need, vv, 0).sum(axis=1) / 1000.0, np.nan)
        return out
    e_act, e_cs = day_energy(y), day_energy(cs)
    act_d = np.isfinite(e_act) & (e_act < CLOUDY_DAY * e_cs)
    ok_d_act = np.isfinite(e_act)

    ev_rows = []
    for m in models:
        for d, lab in LEADS:
            k = f"{m} | {lab} | corrected" if f"{m} | {lab} | corrected" in series else f"{m} | {lab}"
            if k not in series:
                continue
            f = series[k]
            ok = strong & np.isfinite(f)
            pred_h = ok & (f < CLOUDY_HOUR * cs)
            a_h = act_h & ok
            for scale, A, P in (("hour", a_h, pred_h), ("hour ±1 h", widen(a_h) & ok, widen(pred_h) & ok)):
                # hit/miss: each actual cloud hour – was clouds forecast for it (within ±1 h)?
                # false alarm: each forecast cloud hour – no actual clouds (within ±1 h)?
                hits = int((a_h & P).sum())
                misses = int((a_h & ~P).sum())
                false_al = int((pred_h & ~A).sum())
                clear_ok = int((ok & ~a_h & ~pred_h).sum())
                ev_rows.append(dict(model=m, lead=lab, scale=scale, cloud_events=int(a_h.sum()),
                                    hits=hits, misses=misses, false_alarms=false_al, correct_clear=clear_ok))
            e_f = day_energy(f)
            ok_d = ok_d_act & np.isfinite(e_f)
            pred_d = ok_d & (e_f < CLOUDY_DAY * e_cs)
            a_d = act_d & ok_d
            ev_rows.append(dict(model=m, lead=lab, scale="day", cloud_events=int(a_d.sum()),
                                hits=int((pred_d & a_d).sum()), misses=int((a_d & ~pred_d).sum()),
                                false_alarms=int((pred_d & ~a_d).sum()),
                                correct_clear=int((ok_d & ~a_d & ~pred_d).sum())))
    with open(os.path.join(a.out, "cloud_events.csv"), "w", newline="") as fcsv:
        w = csv.DictWriter(fcsv, fieldnames=["model", "lead", "scale", "cloud_events", "hits", "misses",
                                             "false_alarms", "correct_clear"])
        w.writeheader()
        w.writerows(ev_rows)
    if ev_rows:
        print(f"\nCloud events (hour: < {CLOUDY_HOUR:.0%} of clear sky; day: < {CLOUDY_DAY:.0%} of a clear day):")
        for r in ev_rows:
            if r["scale"] in ("hour ±1 h", "day"):
                print(f"  {r['model']:<13} {r['lead']:<13} {r['scale']:<10} "
                      f"{r['hits']} of {r['cloud_events']} cloud {'days' if r['scale'] == 'day' else 'hours'} forecast, "
                      f"{r['misses']} missed, {r['false_alarms']} false alarms")

    # ---- figure: hour by hour at four places
    try:
        import matplotlib
        matplotlib.use("Agg")
        import matplotlib.pyplot as plt
        import matplotlib.dates as mdates
        pal = ["#1f77b4", "#ff7f0e", "#2ca02c", "#d62728", "#9467bd", "#8c564b", "#e377c2",
               "#17becf", "#bcbd22"]
        mcol = {m: pal[i % len(pal)] for i, m in enumerate(models)}
        mcol.update({"ecmwf_ifs025": "#ff7f0e", "gfs_global": "#2ca02c", "icon_global": "#d62728",
                     "best_match": "#1f77b4"})           # same colours as before

        def shown(m, lab):
            """The series the plot shows for model m: corrected if there is a factor."""
            k = f"{m} | {lab} | corrected"
            if k in series:
                return k, f"{m} × {BIAS_FACTOR.get(m, 0):.2f}"
            k = f"{m} | {lab}"
            return (k, m) if k in series else (None, None)
        loc = lambda tt: [t + dt.timedelta(hours=LOCAL_OFFSET_H) - dt.timedelta(minutes=30) for t in tt]

        # one figure per lead time: hour by hour at four places (local time; middle of each hour)
        places = [("Darwin", -12.46, 130.84), ("Alice Springs", -23.70, 133.88),
                  ("Coober Pedy", -29.01, 134.76), ("Adelaide", -34.93, 138.60)]
        nper = len(periods)
        saved = []
        for d, lab in LEADS:
            if not any(shown(m, lab)[0] for m in models):
                print(f"  (no forecasts '{lab}' – no figure)")
                continue
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
                        k, name = shown(m, lab)
                        if k:
                            ax.plot(tt, series[k][i, js], color=mcol[m], lw=1.2, label=name)
                    ax.set_title(f"{pname} ({pts[i][2]}, km {pkm[i]:.0f})", fontsize=10)
                    ax.set_ylabel("GHI [W/m²]")
                    ax.xaxis.set_major_formatter(mdates.DateFormatter("%d.%m"))
                    if r_i == 0 and c_i == 0:
                        ax.legend(fontsize=7, loc="upper right", ncol=2)
            when = "newest forecast (a few hours before)" if d == 0 else \
                f"forecast made {d} day{'s' if d > 1 else ''} before"
            fig.suptitle(f"Sunshine at four places: actual (satellite) vs. {when} – local time ACST")
            fig.tight_layout()
            name = "fig_timeseries_short_lead.png" if d == 0 else f"fig_timeseries_{d}day{'s' if d > 1 else ''}.png"
            fig.savefig(os.path.join(a.out, name), dpi=110)
            plt.close(fig)
            saved.append(name)

        # ECMWF only, corrected: newest, 1 day and 2 days before in one figure
        ec = "ecmwf_ifs025"
        ec_leads = [(lab, f"{ec} | {lab} | corrected", ls, col) for lab, ls, col in
                    (("short lead", "-", "#ff7f0e"), ("1 day ahead", "--", "#c0560a"), ("2 days ahead", ":", "#7a3300"))]
        ec_leads = [x for x in ec_leads if x[1] in series]
        if ec_leads:
            fig, axs = plt.subplots(len(places), nper, figsize=(8 * nper, 2.6 * len(places)), squeeze=False)
            for r_i, (pname, pla, plo) in enumerate(places):
                i = int(np.argmin(haversine_km(lat, lon, pla, plo)))
                for c_i, (d0, d1) in enumerate(pdates):
                    ax = axs[r_i, c_i]
                    js = [j for j, t in enumerate(times) if d0 <= local_date(t) <= d1]
                    tt = loc([times[j] for j in js])
                    ax.fill_between(tt, cs[i, js], color="0.9", label="clear sky")
                    ax.plot(tt, y[i, js], color="k", lw=1.6, label="actual (satellite)")
                    for lab, k, ls, col in ec_leads:
                        ax.plot(tt, series[k][i, js], ls, color=col, lw=1.2,
                                label=f"ECMWF × {BIAS_FACTOR[ec]:.2f}, {lab}")
                    ax.set_title(f"{pname} ({pts[i][2]}, km {pkm[i]:.0f})", fontsize=10)
                    ax.set_ylabel("GHI [W/m²]")
                    ax.xaxis.set_major_formatter(mdates.DateFormatter("%d.%m"))
                    if r_i == 0 and c_i == 0:
                        ax.legend(fontsize=7, loc="upper right", ncol=2)
            fig.suptitle("Sunshine at four places: actual (satellite) vs. ECMWF corrected – local time ACST")
            fig.tight_layout()
            fig.savefig(os.path.join(a.out, "fig_timeseries.png"), dpi=110)
            plt.close(fig)
            saved.append("fig_timeseries.png")

        # cloud events: hits / misses / false alarms per lead time (ECMWF corrected)
        ev = [r for r in ev_rows if r["model"] == "ecmwf_ifs025"] or \
             [r for r in ev_rows if r["model"] == (models[0] if models else "")]
        if ev:
            fig, axs = plt.subplots(1, 2, figsize=(14, 4.8))
            for ax, scale, unit in ((axs[0], "hour ±1 h", "hours"), (axs[1], "day", "point-days")):
                rows = [r for r in ev if r["scale"] == scale]
                labs = [r["lead"] for r in rows]
                x = np.arange(len(rows))
                wbar = 0.27
                for off, key, col, name in ((-wbar, "hits", "#2f9e44", "hit: clouds forecast and came"),
                                            (0, "misses", "#d62728", "missed: clouds came, not forecast"),
                                            (wbar, "false_alarms", "#f0a03c", "false alarm: forecast, did not come")):
                    vals = [r[key] for r in rows]
                    bars = ax.bar(x + off, vals, wbar, color=col, label=name)
                    for b_, v in zip(bars, vals):
                        ax.text(b_.get_x() + b_.get_width() / 2, b_.get_height(), str(v), ha="center",
                                va="bottom", fontsize=7)
                ax.set_xticks(x, [l_.replace(" ahead", "\nbefore").replace("short lead", "newest") for l_ in labs],
                              fontsize=8)
                ax.set_ylabel(f"number of {unit}")
                total = rows[0]["cloud_events"] if rows else 0
                ax.set_title(("Cloudy hours (timing ±1 h)" if scale != "day" else "Cloudy days") +
                             f" – about {total} actual cloud {unit}", fontsize=10)
                ax.grid(axis="y", alpha=0.3)
            h_, l_ = axs[0].get_legend_handles_labels()
            fig.legend(h_, l_, loc="lower center", ncol=3, fontsize=8)
            mname = ev[0]["model"]
            fig.suptitle(f"Did {mname}{' × %.2f' % BIAS_FACTOR[mname] if mname in BIAS_FACTOR else ''} "
                         f"forecast the clouds? Counts per lead time")
            fig.tight_layout(rect=(0, 0.07, 1, 1))
            fig.savefig(os.path.join(a.out, "fig_cloud_events.png"), dpi=110)
            plt.close(fig)

    except ImportError:
        print("matplotlib not installed – no figure (pip install matplotlib)")

    print(f"\n→ results in {a.out}/: data_hourly.csv, cloud_events.csv, fig_cloud_events.png, "
          f"fig_timeseries.png (ECMWF corrected: newest, "
          f"1 and 2 days before) and one figure per lead time "
          f"(fig_timeseries_short_lead.png, fig_timeseries_1day.png … fig_timeseries_{MAX_LEAD_DAYS}days.png)")


if __name__ == "__main__":
    main()
