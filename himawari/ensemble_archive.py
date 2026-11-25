"""
Save today's ECMWF ensemble forecast for the race route, so we can check the cloud
probabilities later.

Why: Open-Meteo keeps ensemble forecasts only for a few days. To find out whether
"70 % chance of clouds" really means clouds in 7 of 10 cases, we have to store the
forecasts ourselves and later compare them with the satellite values (which Open-Meteo
does keep). Run once a day (GitHub Actions: .github/workflows/ensemble_archive.yml).

What is stored (per route point, ECMWF ensemble with 51 versions, corrected × 1.06 like
the normal ECMWF forecast, lead time 0–5 days):
  <run date>_hours.csv.gz   per daylight hour (clear sky > 200 W/m²): how many versions
                            forecast a cloudy hour (< 60 % of clear sky), and the 10 / 50 / 90 %
                            values of the irradiance
  <run date>_days.csv.gz    per local day (ACST): how many versions forecast a cloudy day
                            (< 80 % of a clear day's energy), 10 / 50 / 90 % of the day energy
Cloud probability = cloudy versions / versions. Same thresholds as om_valid.py.

Usage:
    python ensemble_archive.py                       # writes to ensemble_archive/
    python ensemble_archive.py --out some/folder
"""
import argparse
import csv
import datetime as dt
import gzip
import os
import sys

import numpy as np

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import trip_forecast as tf  # noqa: E402  (route, clear sky, fetch, correction factor)

CLOUDY_HOUR = 0.6
CLOUDY_DAY = 0.8
MIN_CS = 200        # judge clouds only when clear sky > 200 W/m² (as om_valid.py)
DAYS = 6            # today + 5 days ahead


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--out", default="ensemble_archive")
    ap.add_argument("--route", default=os.path.join(os.path.dirname(os.path.abspath(__file__)),
                                                     "route_locations.py"))
    a = ap.parse_args()

    pts = tf.load_route(a.route)
    lat = np.array([p[0] for p in pts])
    lon = np.array([p[1] for p in pts])
    factor = tf.BIAS_FACTOR.get(tf.ENS_MODEL, 1.0)
    run_utc = dt.datetime.now(dt.timezone.utc).replace(tzinfo=None)
    run_local_day = (run_utc + tf.LOCAL).date()

    ans = tf.fetch(tf.ENS_URL, pts, dict(models=tf.ENS_MODEL, hourly="shortwave_radiation",
                                         timezone="GMT", forecast_days=DAYS + 1), "ECMWF ensemble")
    if not ans:
        sys.exit("No ensemble data received.")
    first = next(d for d in ans if d and d.get("hourly"))
    keys = sorted(k for k in first["hourly"] if k.startswith("shortwave_radiation"))
    stamps = first["hourly"]["time"]
    # value at time T = average of the hour before T → hour start = T − 1 h
    starts = [dt.datetime.fromisoformat(s) - dt.timedelta(hours=1) for s in stamps]
    print(f"{len(keys)} ensemble versions, {len(stamps)} hours, {len(pts)} route points")

    # g[point, member, hour] in W/m², corrected
    g = np.full((len(pts), len(keys), len(stamps)), np.nan)
    for i, d in enumerate(ans):
        h = (d or {}).get("hourly") or {}
        for k, key in enumerate(keys):
            v = h.get(key)
            if v and len(v) == len(stamps):
                g[i, k] = [np.nan if x is None else x for x in v]
    g *= factor
    cs = np.stack([tf.clear_sky(lat, lon, s + dt.timedelta(minutes=30)) for s in starts], axis=1)
    local_day = [(s + tf.LOCAL).date() for s in starts]
    days = [run_local_day + dt.timedelta(days=k) for k in range(DAYS)]

    os.makedirs(a.out, exist_ok=True)
    tag = f"{run_utc:%Y-%m-%d}"
    hp = os.path.join(a.out, f"{tag}_hours.csv.gz")
    with gzip.open(hp, "wt", newline="") as fh:
        w = csv.writer(fh)
        w.writerow(["run_utc", "point", "hour_start_utc", "lead_days", "clear_sky_wm2", "versions",
                    "cloudy_versions", "p10_wm2", "p50_wm2", "p90_wm2"])
        for j, s in enumerate(starts):
            if local_day[j] not in days:
                continue
            lead = (local_day[j] - run_local_day).days
            for i in range(len(pts)):
                if cs[i, j] <= MIN_CS:
                    continue
                v = g[i, :, j]
                v = v[np.isfinite(v)]
                if not len(v):
                    continue
                p10, p50, p90 = np.percentile(v, [10, 50, 90])
                w.writerow([f"{run_utc:%Y-%m-%dT%H:%M}", i + 1, f"{s:%Y-%m-%dT%H:%M}", lead, round(cs[i, j]),
                            len(v), int((v < CLOUDY_HOUR * cs[i, j]).sum()), round(p10), round(p50), round(p90)])

    dp = os.path.join(a.out, f"{tag}_days.csv.gz")
    with gzip.open(dp, "wt", newline="") as fh:
        w = csv.writer(fh)
        w.writerow(["run_utc", "point", "local_day", "lead_days", "clear_day_kwh_m2", "versions",
                    "cloudy_versions", "p10_kwh_m2", "p50_kwh_m2", "p90_kwh_m2"])
        for d in days:
            idx = [j for j, ld in enumerate(local_day) if ld == d]
            if len(idx) < 24:          # day not complete in the forecast
                continue
            for i in range(len(pts)):
                e_cs = cs[i, idx].sum() / 1000
                sub = g[i][:, idx]
                ok = np.isfinite(sub).all(axis=1)
                if not ok.any():
                    continue
                e = sub[ok].sum(axis=1) / 1000
                p10, p50, p90 = np.percentile(e, [10, 50, 90])
                w.writerow([f"{run_utc:%Y-%m-%dT%H:%M}", i + 1, d.isoformat(), (d - run_local_day).days,
                            round(e_cs, 2), int(ok.sum()), int((e < CLOUDY_DAY * e_cs).sum()),
                            round(p10, 2), round(p50, 2), round(p90, 2)])
    print(f"→ {hp}, {dp}")


if __name__ == "__main__":
    main()
