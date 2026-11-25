"""
Sunshine outlook for the whole race route for the next days.

For every one of the route points the script fetches the hourly irradiance forecast
(ECMWF × 1.06 and GFS × 1.03, same correction as trip_forecast.py) and shows
  - where along the route and when the sun shines (km × time),
  - where clouds are expected (share of the clear-sky sunshine that arrives),
  - the energy per day at every km, per model.
No averages: every value belongs to one route point and one hour or day.

Unlike trip_forecast.py this does not simulate the car: it shows the whole route at once,
so you can see where a cloud band lies and when it moves through.

Output (folder route_outlook/):
  route_outlook.png    3 panels: sunshine km × time, clouds km × time, energy per day along the route
  route_outlook.csv    per route point and local day: km, energy per model, clear-day energy
  printed              per day: cloudy stretches (km from–to) and the darkest point

Usage (in the folder with trip_forecast.py and route_locations.py):
    python route_outlook.py                  # next 5 days, ECMWF + GFS
    python route_outlook.py --days 3
    python route_outlook.py --models ecmwf_ifs025
    python route_outlook.py --no-bias-correction
"""
import argparse
import csv
import datetime as dt
import os
import sys

import numpy as np

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import trip_forecast as tf  # noqa: E402  (route, clear sky, fetch, correction factors)

CLOUDY_HOUR = 0.6   # hour cloudy if less than 60 % of clear-sky sunshine arrives (as om_valid.py)
CLOUDY_DAY = 0.8    # day cloudy if less than 80 % of a clear day's energy (as om_valid.py)
TOWNS = [("Katherine", -14.465, 132.264), ("Tennant Creek", -19.646, 134.191),
         ("Alice Springs", -23.698, 133.880), ("Coober Pedy", -29.013, 134.755),
         ("Port Augusta", -32.490, 137.770)]
NAMES = {"ecmwf_ifs025": "ECMWF", "gfs_global": "GFS", "icon_global": "ICON", "best_match": "best match"}


def route_km(pts):
    km = [0.0]
    for (a, b, _), (c, d, _) in zip(pts[:-1], pts[1:]):
        km.append(km[-1] + float(tf.haversine_km(a, b, c, d)))
    return np.array(km)


def get_model(pts, model, days):
    """Hourly GHI [points × hours] in W/m² and the UTC hour starts, or (None, None)."""
    ans = tf.fetch(tf.FC_URL, pts, dict(hourly="shortwave_radiation", models=model,
                                        timezone="GMT", forecast_days=days + 2), NAMES.get(model, model))
    if not ans:
        return None, None
    stamps = None
    for d in ans:
        if d and d.get("hourly"):
            stamps = d["hourly"]["time"]
            break
    if stamps is None:
        return None, None
    # Open-Meteo: value at time T = average over the hour before T → hour start = T − 1 h
    starts = [dt.datetime.fromisoformat(s) - dt.timedelta(hours=1) for s in stamps]
    m = np.full((len(pts), len(stamps)), np.nan)
    for i, d in enumerate(ans):
        v = ((d or {}).get("hourly") or {}).get("shortwave_radiation")
        if v and len(v) == len(stamps):
            m[i] = [np.nan if x is None else x for x in v]
    return m, starts


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--days", type=int, default=5, help="days ahead (default 5)")
    ap.add_argument("--models", default="ecmwf_ifs025,gfs_global")
    ap.add_argument("--no-bias-correction", action="store_true")
    ap.add_argument("--route", default=os.path.join(os.path.dirname(os.path.abspath(__file__)),
                                                     "route_locations.py"))
    ap.add_argument("--out", default="route_outlook")
    a = ap.parse_args()

    pts = tf.load_route(a.route)
    km = route_km(pts)
    lat = np.array([p[0] for p in pts])
    lon = np.array([p[1] for p in pts])
    models = [m.strip() for m in a.models.split(",") if m.strip()]

    now_local = dt.datetime.now(dt.timezone.utc).replace(tzinfo=None) + tf.LOCAL
    days = [(now_local + dt.timedelta(days=k)).date() for k in range(a.days)]
    print(f"Route: {len(pts)} points, {km[-1]:.0f} km. Days (local): {days[0]:%d.%m.} – {days[-1]:%d.%m.}")

    data = {}
    starts = None
    for m in models:
        print(f"  fetching {NAMES.get(m, m)} …")
        g, st = get_model(pts, m, a.days)
        if g is None:
            continue
        f = 1.0 if a.no_bias_correction else tf.BIAS_FACTOR.get(m, 1.0)
        data[m] = g * f
        starts = starts or st
    if not data:
        sys.exit("No forecast received.")

    # keep hours whose local date is one of the requested days
    loc = [s + tf.LOCAL for s in starts]
    keep = [j for j, t in enumerate(loc) if t.date() in days]
    loc = [loc[j] for j in keep]
    for m in data:
        data[m] = data[m][:, keep]
    # clear sky per point and hour (middle of the hour)
    cs = np.stack([tf.clear_sky(lat, lon, starts[j] + dt.timedelta(minutes=30)) for j in keep], axis=1)

    day_idx = {d: [j for j, t in enumerate(loc) if t.date() == d] for d in days}
    e_cs = {d: cs[:, day_idx[d]].sum(axis=1) / 1000 for d in days}   # kWh/m² per point
    energy = {m: {d: np.where(np.isfinite(data[m][:, day_idx[d]]).all(axis=1),
                              np.nansum(data[m][:, day_idx[d]], axis=1) / 1000, np.nan)
                  for d in days} for m in data}

    os.makedirs(a.out, exist_ok=True)
    with open(os.path.join(a.out, "route_outlook.csv"), "w", newline="") as fh:
        w = csv.writer(fh)
        w.writerow(["point", "km", "latitude", "longitude", "local_day", "clear_day_kwh_m2"]
                   + [f"{NAMES.get(m, m)}_kwh_m2" for m in data])
        for d in days:
            for i, p in enumerate(pts):
                w.writerow([p[2], round(km[i], 1), p[0], p[1], d.isoformat(), round(e_cs[d][i], 2)]
                           + [("" if not np.isfinite(energy[m][d][i]) else round(energy[m][d][i], 2))
                              for m in data])

    # ------------------------------------------------------------ printed summary (no averages)
    main_m = models[0] if models[0] in data else list(data)[0]
    corr = "" if a.no_bias_correction else f" × {tf.BIAS_FACTOR.get(main_m, 1.0):.2f}"
    print(f"\nCloudy stretches per day ({NAMES.get(main_m, main_m)}{corr}: day energy below "
          f"{CLOUDY_DAY:.0%} of a clear day):")
    for d in days:
        e = energy[main_m][d]
        cloudy = np.isfinite(e) & (e < CLOUDY_DAY * e_cs[d])
        runs, i = [], 0
        while i < len(pts):
            if cloudy[i]:
                j = i
                while j + 1 < len(pts) and cloudy[j + 1]:
                    j += 1
                runs.append(f"km {km[i]:.0f}–{km[j]:.0f}")
                i = j + 1
            else:
                i += 1
        k = int(np.nanargmin(e / e_cs[d])) if np.isfinite(e).any() else None
        dark = (f"darkest: km {km[k]:.0f} with {e[k]:.1f} of {e_cs[d][k]:.1f} kWh/m² "
                f"({e[k] / e_cs[d][k]:.0%})") if k is not None else "no data"
        print(f"  {d:%a %d.%m.}: {', '.join(runs) if runs else 'no cloudy stretch'}  |  {dark}")

    # ------------------------------------------------------------ plot
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    import matplotlib.dates as mdates

    g = data[main_m]
    ratio = np.where(cs > 50, g / np.maximum(cs, 1), np.nan)
    t_edges = mdates.date2num(loc + [loc[-1] + dt.timedelta(hours=1)])
    km_edges = np.concatenate([[km[0]], (km[:-1] + km[1:]) / 2, [km[-1]]])

    fig = plt.figure(figsize=(14, 13))
    gs = fig.add_gridspec(3, 2, width_ratios=[50, 1], height_ratios=[1, 1, 1.1], hspace=0.35, wspace=0.03)
    ax1, ax2, ax3 = fig.add_subplot(gs[0, 0]), fig.add_subplot(gs[1, 0]), fig.add_subplot(gs[2, 0])
    title_m = f"{NAMES.get(main_m, main_m)}{corr}"

    pc = ax1.pcolormesh(t_edges, km_edges, g, cmap="inferno", vmin=0, vmax=1100, shading="flat")
    fig.colorbar(pc, cax=fig.add_subplot(gs[0, 1]), label="W/m²")
    ax1.set_title(f"Sunshine forecast along the route – {title_m}")

    cmap = plt.get_cmap("RdYlBu").copy()
    cmap.set_bad("#2b2b2b")
    pc2 = ax2.pcolormesh(t_edges, km_edges, np.clip(ratio, 0, 1.1) * 100, cmap=cmap,
                         vmin=0, vmax=110, shading="flat")
    cb = fig.colorbar(pc2, cax=fig.add_subplot(gs[1, 1]), label="% of clear-sky sunshine")
    cb.ax.axhline(CLOUDY_HOUR * 100, color="k", lw=1.5)
    ax2.set_title(f"Clouds: share of clear-sky sunshine that arrives (below {CLOUDY_HOUR:.0%} = cloudy hour, "
                  f"red/orange) – night dark grey")

    for ax in (ax1, ax2):
        ax.invert_yaxis()
        ax.set_ylabel("route km (Darwin at top)")
        ax.xaxis.set_major_locator(mdates.DayLocator())
        ax.xaxis.set_major_formatter(mdates.DateFormatter("%a %d.%m."))
        ax.xaxis.set_minor_locator(mdates.HourLocator(byhour=[6, 12, 18]))
        ax.grid(axis="x", which="major", color="w", alpha=0.4)
        ax.set_xlabel("local time (ACST)")

    colors = plt.get_cmap("viridis")(np.linspace(0, 0.9, len(days)))
    styles = {m: s for m, s in zip(data, ["-", "--", ":", "-."])}
    for d, c in zip(days, colors):
        for m in data:
            ax3.plot(km, energy[m][d], styles[m], color=c, lw=1.6 if m == main_m else 1.1,
                     label=f"{d:%a %d.%m.}" if m == main_m else None)
        ax3.plot(km, e_cs[d], color=c, lw=0.6, alpha=0.5)
    ax3.plot([], [], "k-", label=title_m)
    for m in data:
        if m != main_m:
            f = "" if a.no_bias_correction else f" × {tf.BIAS_FACTOR.get(m, 1.0):.2f}"
            ax3.plot([], [], "k" + styles[m], label=f"{NAMES.get(m, m)}{f}")
    ax3.plot([], [], "k-", lw=0.6, alpha=0.5, label="clear day")
    ax3.set_xlim(km[0], km[-1])
    ax3.set_xlabel("route km (Darwin → Adelaide)")
    ax3.set_ylabel("energy per day, kWh/m²")
    ax3.set_title("Energy per day at every route point (dips = clouds)")
    ax3.grid(alpha=0.3)
    ax3.legend(ncol=4, fontsize=8, loc="lower left")
    for name, tlat, tlon in TOWNS:
        dist = tf.haversine_km(lat, lon, tlat, tlon)
        if np.min(dist) < 30:
            k_town = km[int(np.argmin(dist))]
            for ax in (ax1, ax2):
                ax.axhline(k_town, color="w", lw=0.5, ls=":")
                ax.text(t_edges[0], k_town, f" {name}", color="w", fontsize=7, va="bottom")
            ax3.axvline(k_town, color="grey", lw=0.6, ls=":")
            ax3.text(k_town, ax3.get_ylim()[1], name, fontsize=7, ha="center", va="bottom", color="grey")

    fig.suptitle(f"Route outlook from {now_local:%d.%m.%Y %H:%M} local time", fontsize=13)
    path = os.path.join(a.out, "route_outlook.png")
    fig.savefig(path, dpi=110, bbox_inches="tight")
    print(f"\n→ {path}, {os.path.join(a.out, 'route_outlook.csv')}")


if __name__ == "__main__":
    main()
