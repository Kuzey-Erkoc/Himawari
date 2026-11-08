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
  metrics_daily.csv        error of the daily energy per forecast
  metrics_by_region.csv    the same split north / south of 23.4°S, and race week only
  summary.md               the main results in words and tables
  fig_timeseries.png       actual vs. forecasts at Darwin, Alice Springs, Adelaide (race week)
  fig_errors.png           mean error and bias per forecast and lead time
  fig_daily_energy.png     daily energy: forecast vs. actual
  fig_hour_of_day.png      bias by hour of day (does the model miss mornings/afternoons?)

Usage:
    pip install numpy requests matplotlib
    python validate_openmeteo.py                       # race seasons 2025 and 2026
    python validate_openmeteo.py --period 2025-08-22 2025-08-29   # one period
    python validate_openmeteo.py --every 3             # more route points (default every 5th)
    python validate_openmeteo.py --models best_match,ecmwf_ifs025,gfs_global,icon_global
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

    # ---- hourly metrics (daytime hours where the actual value and the forecast exist;
    #      n in the tables says how many hours that was)
    common = day & np.isfinite(y)
    fc_names = [k for k in series if k != "actual"]
    rows_h = []
    race = np.zeros(T, bool)
    for s, e in RACE_WEEKS:
        t0 = dt.datetime.strptime(s, "%Y-%m-%d")
        t1 = dt.datetime.strptime(e, "%Y-%m-%d") + dt.timedelta(days=1)
        race |= np.array([t0 <= t < t1 for t in times])
    subsets = [("all", np.ones((n, T), bool)), ("north", north[:, None] & np.ones(T, bool)),
               ("south", ~north[:, None] & np.ones(T, bool)), ("race week", race[None, :] & np.ones((n, 1), bool))]
    # fair comparison: every forecast is scored on the same hours (where all of them exist)
    all_fc = np.all([np.isfinite(series[k]) for k in fc_names], axis=0)
    for sub, msk in subsets:
        for k in fc_names:
            mm = common & msk & all_fc
            r = metrics(np.where(mm, series[k], np.nan), np.where(mm, y, np.nan))
            if r:
                rows_h.append(dict(subset=sub, forecast=k, **{x: round(v, 3) if isinstance(v, float) else v
                                                              for x, v in r.items()}))
    keys = ["subset", "forecast", "n", "mean_actual", "bias", "mae", "rmse", "rel_mae", "corr"]
    with open(os.path.join(a.out, "metrics_hourly.csv"), "w", newline="") as f:
        w = csv.DictWriter(f, fieldnames=keys)
        w.writeheader()
        w.writerows([r for r in rows_h if r["subset"] == "all"])
    with open(os.path.join(a.out, "metrics_by_region.csv"), "w", newline="") as f:
        w = csv.DictWriter(f, fieldnames=keys)
        w.writeheader()
        w.writerows([r for r in rows_h if r["subset"] != "all"])

    # ---- daily energy (local day, kWh/m²), only complete days
    local_day = [(t - dt.timedelta(hours=1) + dt.timedelta(hours=LOCAL_OFFSET_H)).date() for t in times]
    days = sorted(set(local_day))
    dpos = collections.defaultdict(list)
    for j, d in enumerate(local_day):
        dpos[d].append(j)

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
    dser = {k: daily(v) for k, v in series.items() if k != "actual"}
    d_all = np.all([np.isfinite(v) for v in dser.values()], axis=0) & np.isfinite(dy)
    dy = np.where(d_all, dy, np.nan)
    dser = {k: np.where(d_all, v, np.nan) for k, v in dser.items()}
    rows_d = []
    for k, v in dser.items():
        r = metrics(v, dy)
        if r:
            rows_d.append(dict(forecast=k, **{x: round(val, 3) if isinstance(val, float) else val
                                              for x, val in r.items()}))
    with open(os.path.join(a.out, "metrics_daily.csv"), "w", newline="") as f:
        w = csv.DictWriter(f, fieldnames=keys[1:])
        w.writeheader()
        w.writerows(rows_d)

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
        allr = [r for r in rows_h if r["subset"] == "all"]
        fig, ax = plt.subplots(1, 2, figsize=(14, 5.5))
        names = [r["forecast"] for r in allr]
        yv = np.arange(len(names))
        colors = []
        for k in names:
            m = k.split(" | ")[0]
            colors.append(mcol.get(m, "0.5"))
        ax[0].barh(yv, [r["mae"] for r in allr], color=colors)
        ax[0].set_yticks(yv, names, fontsize=8); ax[0].invert_yaxis()
        ax[0].set_xlabel("mean error per daytime hour [W/m²]  (lower = better)")
        ax[0].set_title("Hourly error")
        dd = {r["forecast"]: r for r in rows_d}
        ax[1].barh(yv, [dd[k]["mae"] if k in dd else 0 for k in names], color=colors)
        ax[1].set_yticks(yv, names, fontsize=8); ax[1].invert_yaxis()
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
            if "|" not in k or "2 days" in k:
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
          "## Daily energy",
          table(rows_d, "kWh/m²", 1.0, 2),
          "",
          "## North / south / race week (hourly, mean error W/m²)",
          "| Forecast | north | south | race week |", "|---|---|---|---|"]
    by = collections.defaultdict(dict)
    for r in rows_h:
        by[r["forecast"]][r["subset"]] = r["mae"]
    for k in fc_names:
        if k in by:
            md.append(f"| {k.replace(' | ', ', ')} | " + " | ".join(f"{by[k].get(s, float('nan')):.0f}"
                                                for s in ("north", "south", "race week")) + " |")
    md += ["", "Bias > 0: forecast too sunny. 'always clear sky' and 'same as yesterday' are the "
           "references a forecast has to beat."]
    open(os.path.join(a.out, "summary.md"), "w", encoding="utf-8").write("\n".join(md))
    print("\n" + "\n".join(md))
    print(f"\n→ results in {a.out}/")


if __name__ == "__main__":
    main()