#!/usr/bin/env python3
"""
Himawari cloud nowcasting for Australia – B13 (infrared) and TRM (colour) images.

By default both bands are processed in one run. The satellite images are kept in
memory only and never saved; each run first deletes the previous run's results, so
result_b13/ and result_trm/ only contain what the current run produced
(exception: skill_history.csv, which collects forecast quality over many runs).

For each band the script:
1. Downloads the aus_<band> images of the last 24 h from JMA (10-minute steps,
   several requests in parallel).
2. Finds the coastline/state-border overlay automatically (overlay pixels stay constant
   for 24 h, clouds and the day/night cycle change everything else) and removes it
   (normalised blur fill), so the lines don't disturb the motion estimate.
3. Measures cloud motion (DIS optical flow, image pairs in parallel) over the full 24 h
   -> motion_24h.csv/.png; the forecast motion field averages the last 6 image pairs
   (1 h, --pairs).
4. Forecasts by advecting the newest image along the measured motion (each step continues
   from the previous one). The lead time is chosen automatically so the forecast reaches
   10 min past *now* (JMA images are published ~10-20 min late).
5. Skill curve: for many start times of the last 24 h (every 30 min, --skill-every) the
   forecast is run up to +120 min (--skill-max) and compared with the image that was
   actually observed, against persistence ("nothing changes"). With a route, also along
   the route corridor. Results of all runs are collected in skill_history.csv, so the
   curve becomes more reliable over days.
6. Route (route_locations.py next to the script, or --route): forecast along the race
   route up to +60 min (--route-lead), as a table, a km/time chart, and a map with the
   route drawn on the forecast image. --car-km/--speed add the car's expected position.
7. Sunshine (needs internet, --no-irradiance to switch off): satellite-derived irradiance
   (GHI, W/m², from Himawari) for the route points is fetched from Open-Meteo for the last
   24 h. From it the script learns which image brightness means how much sun
   (clear-sky index = GHI / clear-sky GHI), turns the cloud forecast into W/m² along the
   route, and rates it against "smart persistence" (today's clear-sky index stays),
   separately north and south of the Tropic of Capricorn (23.4°S): tropical north
   (dry season in August, thunderstorms from about October) and dry interior/south
   (winter fronts in August).
   Note: the reference is itself a satellite estimate, not a ground measurement, and the
   calibration is learned from the same 24 h it is tested on.
   Weather models: irradiance forecasts of ECMWF, GFS, … (Open-Meteo, --nwp-models) are
   fetched for the same points and scored in the same tests -> nwp_compare.png. For past
   times Open-Meteo joins the first hours of each model run, which flatters the models a
   little (in a race the newest available run is usually several hours old);
   --nwp-source previous_day1 shows the forecast made 24 h before instead. Model data are
   cached in openmeteo_nwp_cache.json for 60 min (--nwp-cache-min), because models update
   only every few hours and Open-Meteo limits free requests. All models forecast too little
   sunshine; their values are multiplied by a factor per model (BIAS_FACTOR, ECMWF × 1.06,
   from om_valid.py and the car's weather station), off with --no-bias-correction.
   Now + cloud change: the newest measured sunshine plus the change the moving clouds bring
   (difference of the calibrated brightness now and later), so the unsure brightness ->
   sunshine conversion only enters as a difference.
   Blend: the best estimate of sunshine combines "cloudiness stays", both satellite forecasts
   and the first weather model, with weights per lead time and route part that are learned
   from the retrospective tests (least squares, weights ≥ 0, sum 1). The test statistics are
   kept in blend_stats.json and grow with every run. The blend is only scored with weights
   learned from earlier days, so it is not judged on data it was fitted to.
   Season: the race is in late August (dry season in the north, winter fronts in the south).
   Data collected in other months describe other weather; --history-months 7,8,9 limits the
   history curves and the blend weights to race-season tests (all tests stay stored).
8. Saves (in result_<band>/):
     forecast_YYYYMMDD_HHMM.png   every 10-minute forecast step (up to the auto lead)
     final_forecast.png           the last forecast step
     motion_field.png             motion used for the forecast (arrows on the final forecast)
     change_HHMM_HHMM.png         what changes in the last forecast step
                                  (red = cloud arriving, blue = cloud leaving)
     animation.gif                the forecast steps (1.5 s per frame, --gif-delay)
     border_mask.png              detected overlay lines (for checking)
     skill_curve.csv/.png         forecast quality per lead time (this run + history)
     skill_history.csv            all skill tests of all runs (kept between runs)
     route_map.png                route on the forecast image, with borders and lat/lon grid
     route_forecast.csv/.png      brightness and sunshine (W/m²) per route point and time
     route_car.csv                what the car meets on its way (with --car-km)
     route_calibration.png        learned brightness -> sunshine relation
     nwp_compare.png              our sunshine forecast vs. weather models and blend
     blend_stats.json             test statistics for the blend weights (kept between runs)

Map position of the route: by default the image is assumed to cover 110-155°E, 10-45°S
(JMA's stated 'aus' area) in a plain lat/lon grid. Check route_map.png: the drawn grid
and route must match the coastline. If not, pass --georef georef.json with control
points read from the image (pixel x/y of known places, at least 3, better 6+):
    {"projection": "equirect",
     "points": [[-12.46, 130.84, 512, 101], [-34.93, 138.60, 678, 860], ...]}
         (each point: latitude, longitude, x, y in the full, uncropped image)
or a different extent: {"projection": "mercator", "extent": [110, 155, -45, -10],
                        "pixel_extent": [x0, y0, x1, y1]}

The 'cloud' value is image brightness (B13: bright = cold = high cloud). It is a proxy,
not solar irradiance: thin cirrus looks bright but lets much sunlight through, low warm
cloud can be hard to see in B13.

B13 (10.4 µm infrared) works day and night; TRM is a colour image that may be dark at night.

Usage:
    pip install opencv-python numpy requests matplotlib pillow
    python himawari_nowcast.py                     # B13 + TRM: download, analysis, forecast
    python himawari_nowcast.py --band b13          # only B13
    python himawari_nowcast.py --car-km 850 --speed 80   # + where the car will be
    python himawari_nowcast.py --route other_route.csv  # a different route file
    python himawari_nowcast.py --no-route          # without route
    python himawari_nowcast.py --history-months 7,8,9   # history/blend from race season only
  (route_locations.py next to this script is used automatically)
    python himawari_nowcast.py --no-download       # read images from himawari_<band>/ folders
    python himawari_nowcast.py --lead 30           # fixed +30 min from the newest image
    python himawari_nowcast.py --ahead 20          # reach 20 min past now (instead of 10)
    python himawari_nowcast.py --draw-borders      # draw coastline/state borders on outputs
    python himawari_nowcast.py --roi 400 300 600 500   # crop x y w h (pixels)
    python himawari_nowcast.py --workers 4         # CPU threads (default: all cores)
"""
import argparse
import ast
import collections
import concurrent.futures as cf
import csv
import datetime as dt
import json
import math
import os
import sys
import threading
import time

import cv2
import numpy as np

URL_TMPL = "https://www.data.jma.go.jp/mscweb/data/himawari/img/aus/aus_{band}_{hhmm}.jpg"
BANDS = {
    "b13": dict(cloud_v=150, can_be_dark=False),   # infrared: bright = cold = high cloud
    "trm": dict(cloud_v=170, can_be_dark=True),    # colour image
}
DEFAULT_ROUTE = "route_locations.py"        # used automatically if it lies next to the script
AUS_EXTENT = (110.0, 155.0, -45.0, -10.0)   # lon_min, lon_max, lat_min, lat_max (JMA 'aus')
STEP_MIN = 10
MAX_GAP_STEPS = 3          # pairs further apart than 3 × 10 min are not used for motion
NO_LM_MIN_AGE = 60         # without Last-Modified, only accept slots at least this old (min)
DL_WORKERS = 8             # parallel HTTP requests to JMA
DARK_GRAY = 30             # skill tests on frames darker than this are skipped (TRM at night)
OM_URL = "https://satellite-api.open-meteo.com/v1/archive"   # irradiance reference
OM_MODEL = "satellite_radiation_seamless"                    # picks Himawari for Australia
OM_BATCH = 50              # route points per request
OM_TIMEOUT = 20            # seconds to wait for an Open-Meteo answer
MIN_COSZ = 0.17            # irradiance only evaluated with the sun > ~10° above the horizon
NORTH_LAT = -23.44         # Tropic of Capricorn: route split into 'north' (tropics) / 'south'
OM_FC_URL = "https://api.open-meteo.com/v1/forecast"                 # weather-model forecasts
OM_PREV_URL = "https://previous-runs-api.open-meteo.com/v1/forecast"  # forecasts made 24 h ago
NWP_MODELS = ("best_match", "ecmwf_ifs025", "gfs_global", "icon_global")   # columns in history
NWP_DEFAULT = "ecmwf_ifs025,gfs_global,best_match"     # the first one is used in the blend
# Bias correction: all models forecast systematically too little sunshine. Factor per model,
# from om_valid.py (Aug 2025+2026 vs. satellite, learned on the other year) and checked against
# the car's own weather station in the WSC 2025 (validate_car.py; the satellite reads ~3 %
# too bright there, so the factors are taken ~3 % below the satellite-based ones).
BIAS_FACTOR = {"ecmwf_ifs025": 1.06, "gfs_global": 1.03, "icon_global": 1.09, "best_match": 1.06}
OUTPUT_PREFIXES = ("forecast_", "final_forecast", "change_", "motion_field", "motion_24h",
                   "animation", "border_mask", "skill_curve", "route_", "nwp_compare")


def chunked_map(ex, fn, items, chunk):
    """Like ex.map, but only `chunk` tasks in flight, so finished results (images) don't
    pile up in memory faster than they are consumed."""
    for i in range(0, len(items), chunk):
        yield from ex.map(fn, items[i:i + chunk])


# ---------------------------------------------------------------- download
def slot_times(now_utc):
    """All 144 slots of the last 24 h, oldest first. JMA keeps only HHMM (rolling)."""
    last = now_utc.replace(minute=now_utc.minute - now_utc.minute % STEP_MIN,
                           second=0, microsecond=0)
    return [last - dt.timedelta(minutes=STEP_MIN * i) for i in range(143, -1, -1)]


def download(band, now_utc):
    """Fetch the last 24 h into memory only – nothing is written to disk.
    Requests run in parallel; the checks run afterwards in time order.
    Returns (list of JPEG bytes, list of times)."""
    import requests
    local = threading.local()

    def fetch(t):
        if not hasattr(local, "sess"):  # one session per thread
            local.sess = requests.Session()
            local.sess.headers["User-Agent"] = "himawari-nowcast/2.4 (research)"
        url = URL_TMPL.format(band=band, hhmm=t.strftime("%H%M"))
        try:
            r = local.sess.get(url, timeout=30)
        except requests.RequestException as e:
            print(f"  ! {url}: {e}")
            return None
        return r.status_code, r.content, r.headers.get("Last-Modified")

    slots = slot_times(now_utc)
    with cf.ThreadPoolExecutor(DL_WORKERS) as ex:
        responses = list(ex.map(fetch, slots))

    srcs, times = [], []
    skipped_old = skipped_no_lm = skipped_dup = 0
    for t, resp in zip(slots, responses):
        if resp is None:
            continue
        status, content, lm = resp
        if status != 200 or len(content) < 1000:
            continue  # e.g. 02:40 UTC (maintenance slot) or not yet published
        lm_t = None
        if lm:  # file must belong to the assumed date, not to yesterday
            try:
                lm_t = dt.datetime.strptime(lm, "%a, %d %b %Y %H:%M:%S GMT")
            except ValueError:
                lm_t = None
            if lm_t is not None and abs((lm_t - t).total_seconds()) > 3 * 3600:
                skipped_old += 1
                continue
        if lm_t is None:
            # Can't verify the date. The newest slots may not be published yet, in which
            # case the file is still yesterday's – only accept slots that are surely updated.
            if (now_utc - t).total_seconds() < NO_LM_MIN_AGE * 60:
                skipped_no_lm += 1
                continue
        if srcs and content == srcs[-1]:
            skipped_dup += 1  # identical bytes as the previous slot -> stale copy
            continue
        srcs.append(content)
        times.append(t)
    print(f"Download: {len(srcs)} {band.upper()} images (kept in memory only)")
    if skipped_old or skipped_no_lm or skipped_dup:
        print(f"  skipped: {skipped_old} from yesterday, {skipped_no_lm} unverifiable "
              f"(no Last-Modified, < {NO_LM_MIN_AGE} min old), {skipped_dup} duplicates")
    return srcs, times


# ---------------------------------------------------------------- input
def list_frames(band, indir):
    prefix = band + "_"
    if not os.path.isdir(indir):
        return [], []
    files = sorted(f for f in os.listdir(indir) if f.startswith(prefix) and f.endswith(".jpg"))
    n = len(prefix)
    return ([os.path.join(indir, f) for f in files],
            [dt.datetime.strptime(f[n:n + 13], "%Y%m%d_%H%M") for f in files])


def read(src, roi=None):
    """src = JPEG bytes (downloaded) or a file path (--no-download)."""
    if isinstance(src, (bytes, bytearray)):
        img = cv2.imdecode(np.frombuffer(src, np.uint8), cv2.IMREAD_COLOR)
    else:
        img = cv2.imread(src, cv2.IMREAD_COLOR)
    if img is not None and roi:
        x, y, w, h = roi
        img = img[y:y + h, x:x + w]
    return img


# ---------------------------------------------------------------- route + map position
def load_route(path):
    """Route points as [(lat, lon, name), ...] from a .py file with a list of dicts
    (read with ast, never executed), a .json list/{"locations": [...]}, or a .csv with
    latitude,longitude[,name] columns."""
    ext = os.path.splitext(path)[1].lower()
    if ext == ".py":
        tree = ast.parse(open(path, encoding="utf-8").read())
        data = None
        for node in tree.body:
            if isinstance(node, ast.Assign):
                try:
                    val = ast.literal_eval(node.value)
                except ValueError:
                    continue
                if isinstance(val, list):
                    data = val
                    break
        if data is None:
            sys.exit(f"No list of locations found in {path}.")
    elif ext == ".json":
        data = json.load(open(path, encoding="utf-8"))
        if isinstance(data, dict):
            data = data.get("locations", [])
    else:
        data = list(csv.DictReader(open(path, encoding="utf-8")))
    pts = [(float(d["latitude"]), float(d["longitude"]), str(d.get("name", f"P{i + 1}")))
           for i, d in enumerate(data)]
    if len(pts) < 2:
        sys.exit(f"Route {path} needs at least 2 points.")
    return pts


def haversine_km(lat1, lon1, lat2, lon2):
    p1, p2 = np.radians(lat1), np.radians(lat2)
    a = (np.sin((p2 - p1) / 2) ** 2
         + np.cos(p1) * np.cos(p2) * np.sin(np.radians(lon2 - lon1) / 2) ** 2)
    return 6371.0 * 2 * np.arcsin(np.sqrt(a))


class Georef:
    """lat/lon -> pixel of the full (uncropped) image."""

    def __init__(self, cfg, width, height):
        self.proj = cfg.get("projection", "equirect")
        if self.proj not in ("equirect", "mercator"):
            sys.exit("georef projection must be 'equirect' or 'mercator'.")
        pts = cfg.get("points") or []
        if len(pts) >= 3:
            p = np.array(pts, float)
            X, Y = self._p(p[:, 0], p[:, 1])
            self.order = 2 if (len(pts) >= 6 and cfg.get("order", 1) == 2) else 1
            A = self._design(X, Y)
            self.cx, *_ = np.linalg.lstsq(A, p[:, 2], rcond=None)
            self.cy, *_ = np.linalg.lstsq(A, p[:, 3], rcond=None)
            res = np.hypot(A @ self.cx - p[:, 2], A @ self.cy - p[:, 3])
            self.desc = (f"{len(pts)} control points ({self.proj}, order {self.order}), "
                         f"max residual {res.max():.1f} px")
        else:
            lon0, lon1, lat0, lat1 = cfg.get("extent", AUS_EXTENT)
            x0, y0, x1, y1 = cfg.get("pixel_extent", [0, 0, width, height])
            Xa, Ya = self._p(np.array([lat0, lat1]), np.array([lon0, lon1]))
            # x = x0 + (X - X_west) * sx ; y = y0 + (Y_north - Y) * sy
            sx = (x1 - x0) / (Xa[1] - Xa[0])
            sy = (y1 - y0) / (Ya[1] - Ya[0])
            self.order = 1
            self.cx = np.array([sx, 0.0, x0 - Xa[0] * sx])
            self.cy = np.array([0.0, -sy, y0 + Ya[1] * sy])
            self.desc = (f"extent {lon0}-{lon1}°E, {lat0}-{lat1}° ({self.proj}) on pixels "
                         f"{x0},{y0}-{x1},{y1}" + ("" if "extent" in cfg else
                                                   " – default, check route_map.png"))

    def _p(self, lat, lon):
        lat, lon = np.asarray(lat, float), np.asarray(lon, float)
        if self.proj == "mercator":
            return lon, np.degrees(np.log(np.tan(np.pi / 4 + np.radians(lat) / 2)))
        return lon, lat

    def _design(self, X, Y):
        cols = [X, Y, np.ones_like(X)]
        if self.order == 2:
            cols += [X * X, X * Y, Y * Y]
        return np.stack(cols, axis=1)

    def to_px(self, lat, lon):
        X, Y = self._p(lat, lon)
        A = self._design(np.atleast_1d(X), np.atleast_1d(Y))
        return A @ self.cx, A @ self.cy


def build_route(pts, georef, roi, shape, radius):
    lat = np.array([p[0] for p in pts])
    lon = np.array([p[1] for p in pts])
    names = [p[2] for p in pts]
    km = np.concatenate([[0.0], np.cumsum(haversine_km(lat[:-1], lon[:-1], lat[1:], lon[1:]))])
    x, y = georef.to_px(lat, lon)
    if roi:
        x, y = x - roi[0], y - roi[1]
    h, w = shape
    inside = (x >= 0) & (x <= w - 1) & (y >= 0) & (y <= h - 1)
    ix = np.clip(np.round(x), 0, w - 1).astype(int)
    iy = np.clip(np.round(y), 0, h - 1).astype(int)
    corridor = np.zeros(shape, np.uint8)
    poly = np.stack([ix[inside], iy[inside]], axis=1).astype(np.int32)
    if len(poly) >= 2:
        cv2.polylines(corridor, [poly], False, 255, 2 * radius + 1)
    return dict(lat=lat, lon=lon, name=names, km=km, x=x, y=y, ix=ix, iy=iy, inside=inside,
                corridor=corridor > 0, k=2 * radius + 1)


def route_values(gray, route):
    """Mean brightness in a (2r+1)² box around each route point (NaN outside the image)."""
    b = cv2.blur(gray.astype(np.float32), (route["k"], route["k"]))
    v = b[route["iy"], route["ix"]]
    return np.where(route["inside"], v, np.nan)


# ---------------------------------------------------------------- irradiance
def cos_zenith(lat, lon, t):
    """Cosine of the solar zenith angle for UTC time t (NOAA approximation, ~0.5°)."""
    la = np.radians(np.asarray(lat, float))
    hour = t.hour + t.minute / 60 + t.second / 3600
    g = 2 * np.pi / 365 * (t.timetuple().tm_yday - 1 + (hour - 12) / 24)
    eqt = 229.18 * (0.000075 + 0.001868 * np.cos(g) - 0.032077 * np.sin(g)
                    - 0.014615 * np.cos(2 * g) - 0.040849 * np.sin(2 * g))
    decl = (0.006918 - 0.399912 * np.cos(g) + 0.070257 * np.sin(g) - 0.006758 * np.cos(2 * g)
            + 0.000907 * np.sin(2 * g) - 0.002697 * np.cos(3 * g) + 0.00148 * np.sin(3 * g))
    ha = np.radians((hour * 60 + eqt + 4 * np.asarray(lon, float)) / 4 - 180)
    return np.sin(la) * np.sin(decl) + np.cos(la) * np.cos(decl) * np.cos(ha)


def clear_sky_ghi(lat, lon, t):
    """Clear-sky global horizontal irradiance in W/m² (Haurwitz model; 0 at night).
    Returns (ghi_clear, cos_zenith)."""
    cz = cos_zenith(lat, lon, t)
    czp = np.maximum(cz, 1e-3)
    return np.where(cz > 0, 1098.0 * czp * np.exp(-0.057 / czp), 0.0), cz


_REF_CACHE = {}


_OM_DOWN = set()           # Open-Meteo hosts that did not answer in this run -> not retried


def om_points(url, pts, base_params, variables, label):
    """Query an Open-Meteo endpoint for many points: OM_BATCH points per request (after the
    first batch the others run in parallel), falling back to one request per point only if the
    service refuses lists. variables = candidate names for the 'hourly' parameter, tried in
    order. A host that does not answer (timeout, no connection) is not asked again in this run.
    Returns (variable used, {point index: JSON}) or (None, error text)."""
    import requests
    host = url.split("/")[2]
    if host in _OM_DOWN:
        return None, f"{host} did not answer earlier in this run – skipped"
    local = threading.local()
    state = dict(var=None, err="", down=False)

    def session():
        if not hasattr(local, "s"):
            local.s = requests.Session()
            local.s.headers["User-Agent"] = "himawari-nowcast/2.7 (research)"
        return local.s

    def query(idx):
        for var in ([state["var"]] if state["var"] else variables):
            if state["down"]:
                return None
            params = dict(base_params, hourly=var,
                          latitude=",".join(f"{pts[i][0]:.4f}" for i in idx),
                          longitude=",".join(f"{pts[i][1]:.4f}" for i in idx))
            try:
                r = session().get(url, params=params, timeout=OM_TIMEOUT)
            except requests.RequestException as e:
                state["err"], state["down"] = f"no answer from {host} ({type(e).__name__})", True
                _OM_DOWN.add(host)
                return None
            if r.status_code == 200:
                data = r.json()
                data = data if isinstance(data, list) else [data]
                if len(data) == len(idx):
                    state["var"] = var
                    return data
                state["err"] = f"{len(data)} answers for {len(idx)} points"
            else:
                state["err"] = f"HTTP {r.status_code}: {r.text[:150]}"
                if r.status_code == 429:          # rate limit: further requests won't help
                    state["down"] = True
                    return None
        return None

    def per_point(idx):
        with cf.ThreadPoolExecutor(4) as ex:
            single = list(ex.map(lambda i: query([i]), idx))
        return [s[0] if s else None for s in single]

    batches = [list(range(b, min(len(pts), b + OM_BATCH))) for b in range(0, len(pts), OM_BATCH)]
    first = query(batches[0])
    single_mode = False
    if first is None and len(batches[0]) > 1 and not state["down"]:
        if query(batches[0][:1]) is not None:   # only single points allowed -> one each
            single_mode = True
            first = per_point(batches[0])
    if first is None or all(r is None for r in first):
        return None, state["err"] or f"{label}: no data"
    answers = dict(zip(batches[0], first))
    if len(batches) > 1:
        fn = per_point if single_mode else query
        with cf.ThreadPoolExecutor(min(4, len(batches) - 1)) as ex:
            for idx, res in zip(batches[1:], ex.map(fn, batches[1:])):
                if res is None or all(r is None for r in res):
                    return None, state["err"] or f"{label}: no data"
                answers.update(zip(idx, res))
    return state["var"], answers


def om_series(answers, var, n, shift=dt.timedelta(0)):
    """{UTC datetime: array over points} from per-point Open-Meteo answers."""
    out = collections.defaultdict(lambda: np.full(n, np.nan))
    for i, d in answers.items():
        if not d:
            continue
        block = next((d[k] for k in ("minutely_10", "minutely_15", "hourly")
                      if isinstance(d.get(k), dict) and var in d[k]), None)
        if block is None:
            continue
        for ts, val in zip(block["time"], block[var]):
            if val is not None:
                out[dt.datetime.strptime(ts[:16], "%Y-%m-%dT%H:%M") - shift][i] = val
    return dict(out)


def fetch_reference(pts, t0, t1):
    """Satellite-derived global irradiance (GHI, W/m²) for the route points from Open-Meteo
    (from Himawari, 10-min steps, the past – no forecast). Returns {UTC datetime: array over
    points, NaN = missing}, or None if the service can't be used. Cached for the second band."""
    key = (len(pts), pts[0], pts[-1], t0.date(), t1.date())
    if key in _REF_CACHE:
        return _REF_CACHE[key]
    var, answers = om_points(
        OM_URL, pts, dict(models=OM_MODEL, temporal_resolution="native", timezone="GMT",
                          start_date=f"{t0:%Y-%m-%d}", end_date=f"{t1:%Y-%m-%d}"),
        ["shortwave_radiation_instant", "shortwave_radiation"], "reference")
    if var is None:
        print(f"! Irradiance reference (Open-Meteo) not available: {answers} "
              f"– continuing with brightness only.")
        _REF_CACHE[key] = None
        return None
    # a backward-averaged 10-min value stamped T covers T-10…T; the Australian part of the
    # image with slot time T is scanned a few minutes after T -> value stamped T+10 fits best
    shift = dt.timedelta(0) if var.endswith("_instant") else dt.timedelta(minutes=STEP_MIN)
    ref = om_series(answers, var, len(pts), shift)
    n_val = sum(int(np.isfinite(v).sum()) for v in ref.values())
    print(f"Irradiance reference (Open-Meteo, {var}): {len(ref)} time steps, {n_val} values")
    _REF_CACHE[key] = ref if n_val else None
    return _REF_CACHE[key]


# ---------------------------------------------------------------- weather-model forecasts
_NWP_MEM = {}
_NWP_FAIL = {}          # model fetches that failed in this run -> not retried for the 2nd band


class NWPForecast:
    """Weather-model irradiance forecast of one model at the route points, stored as clear-sky
    index k = GHI / GHI_clear at its time stamps. k(t) interpolates linearly between the
    stamps, so the day/night curve comes from our own clear-sky model (smooth), the clouds
    from the weather model. Points/times without daylight are NaN."""

    def __init__(self, model, stamps, ghi, lat, lon, factor=1.0):
        self.model = model
        self.factor = factor
        order = np.argsort(stamps)
        self.stamps = [stamps[i] for i in order]
        self.sec = np.array([s.timestamp() for s in self.stamps])
        g = np.asarray(ghi, float)[order]
        k = np.full_like(g, np.nan)
        for j, s in enumerate(self.stamps):
            gcs, cz = clear_sky_ghi(lat, lon, s)
            ok = (cz > MIN_COSZ) & np.isfinite(g[j])
            k[j, ok] = np.clip(g[j, ok] / gcs[ok] * factor, 0, 1.5)
        self.k = k

    def k_at(self, t):
        s = t.timestamp()
        j = int(np.searchsorted(self.sec, s))
        if j <= 0 or j >= len(self.sec):
            if len(self.sec) and abs(s - self.sec[min(max(j, 0), len(self.sec) - 1)]) < 1:
                return self.k[min(max(j, 0), len(self.sec) - 1)]
            return np.full(self.k.shape[1], np.nan)
        a, b = self.k[j - 1], self.k[j]
        w = (s - self.sec[j - 1]) / (self.sec[j] - self.sec[j - 1])
        out = (1 - w) * a + w * b
        # near sunrise/sunset one neighbour is NaN: use the other one if it is close in time
        out = np.where(np.isnan(out) & (w < 0.5), a, out)
        out = np.where(np.isnan(out) & (w >= 0.5), b, out)
        return out

    def ghi_at(self, t, lat, lon):
        return clear_sky_ghi(lat, lon, t)[0] * self.k_at(t)


def fetch_nwp(pts, t0, t1, models, source, cache_path, cache_min, bias=None):
    """Irradiance forecasts of weather models (ECMWF, GFS, …) from Open-Meteo for the route
    points, t0 … t1 (past and future). source 'forecast' = latest runs (for the past Open-Meteo
    joins the first hours of each run, i.e. short lead times); 'previous_day1' = what the model
    predicted 24 h earlier. Results are cached in cache_path for cache_min minutes, because
    weather models update only every 1–12 h and Open-Meteo limits free requests.
    bias: {model: factor} applied to the forecast (bias correction); the cache keeps the raw
    values. Returns {model: NWPForecast}."""
    lat = np.array([p[0] for p in pts])
    lon = np.array([p[1] for p in pts])
    key = json.dumps([source, len(pts), pts[0][:2], pts[-1][:2], f"{t0:%Y-%m-%d}",
                      f"{t1:%Y-%m-%d}"])
    now_s = dt.datetime.now(dt.timezone.utc).timestamp()
    cache = {}
    if cache_path and os.path.exists(cache_path):
        try:
            cache = json.load(open(cache_path))
        except (OSError, ValueError):
            cache = {}
    out, fresh, notes = {}, {}, []
    if source == "previous_day1":
        url = OM_PREV_URL
        variables = ["shortwave_radiation_instant_previous_day1", "shortwave_radiation_previous_day1"]
    else:
        url = OM_FC_URL
        variables = ["shortwave_radiation_instant", "shortwave_radiation"]

    def get(model):
        return om_points(url, pts, dict(models=model, timezone="GMT",
                                        start_date=f"{t0:%Y-%m-%d}", end_date=f"{t1:%Y-%m-%d}"),
                         variables, model)

    entries, todo = {}, []
    for model in models:
        mk = f"{key}|{model}"
        hit = _NWP_MEM.get(mk) or cache.get(mk)
        if mk in _NWP_FAIL:
            notes.append(f"{model}: {_NWP_FAIL[mk]} (not retried)")
        elif hit and now_s - hit["fetched"] < cache_min * 60:
            entries[model] = (hit, "cached")
        else:
            todo.append(model)
    if todo:   # all models at the same time
        with cf.ThreadPoolExecutor(len(todo)) as ex:
            results = dict(zip(todo, ex.map(get, todo)))
    for model in todo:
        var, answers = results[model]
        mk = f"{key}|{model}"
        if var is None:
            _NWP_FAIL[mk] = answers
            notes.append(f"{model}: {answers}")
            continue
        # hourly values without '_instant' are averages over the hour before the stamp
        shift = dt.timedelta(0) if "_instant" in var else dt.timedelta(minutes=30)
        series = om_series(answers, var, len(pts), shift)
        if not series or not any(np.isfinite(v).any() for v in series.values()):
            _NWP_FAIL[mk] = "no values"
            notes.append(f"{model}: no values")
            continue
        entry = dict(fetched=now_s, var=var,
                     stamps=[s.strftime("%Y-%m-%dT%H:%M") for s in series],
                     ghi=[[None if not np.isfinite(x) else round(float(x), 1) for x in v]
                          for v in series.values()])
        fresh[mk] = entry
        entries[model] = (entry, "fetched")
    for model in models:          # keep the order given on the command line
        if model not in entries:
            continue
        entry, how = entries[model]
        _NWP_MEM[f"{key}|{model}"] = entry
        stamps = [dt.datetime.strptime(s, "%Y-%m-%dT%H:%M") for s in entry["stamps"]]
        ghi = [[np.nan if x is None else x for x in row] for row in entry["ghi"]]
        fac = (bias or {}).get(model, 1.0)
        out[model] = NWPForecast(model, stamps, ghi, lat, lon, fac)
        notes.append(f"{model} {how}" + (f" (× {fac:.2f} bias-corrected)" if fac != 1.0 else ""))
    if fresh and cache_path:
        cache = {k: v for k, v in cache.items() if now_s - v.get("fetched", 0) < 86400}
        cache.update(fresh)
        try:
            with open(cache_path, "w") as f:
                json.dump(cache, f)
        except OSError:
            pass
    print(f"Weather-model forecasts (Open-Meteo, {source}): " + "; ".join(notes))
    return out


def fit_calibration(grays, tlist, route, ref, min_pairs=300):
    """Learn brightness -> clear-sky index k = GHI / GHI_clear from all daytime images of the
    last 24 h at the route points: median k per brightness bin, forced to fall with
    brightness (brighter = more cloud). Returns dict(k_of=function, centers, k, n) or None."""
    ok = route["inside"]
    bs, ks = [], []
    for g, t in zip(grays, tlist):
        G = ref.get(t)
        if G is None:
            continue
        gcs, cz = clear_sky_ghi(route["lat"], route["lon"], t)
        m = ok & (cz > MIN_COSZ) & np.isfinite(G)
        if m.any():
            bs.append(route_values(g, route)[m])
            ks.append(G[m] / gcs[m])
    if not bs:
        return None
    b, k = np.concatenate(bs), np.clip(np.concatenate(ks), 0, 1.5)
    if len(b) < min_pairs:
        print(f"! Calibration: only {len(b)} daytime pairs (need {min_pairs}) – brightness only.")
        return None
    edges = np.linspace(0, 256, 17)
    centers = (edges[:-1] + edges[1:]) / 2
    med = np.full(len(centers), np.nan)
    for i in range(len(centers)):
        sel = (b >= edges[i]) & (b < edges[i + 1])
        if sel.sum() >= 15:
            med[i] = np.median(k[sel])
    valid = np.isfinite(med)
    if valid.sum() < 2:
        return None
    med = np.minimum.accumulate(np.interp(centers, centers[valid], med[valid]))
    return dict(k_of=lambda x: np.clip(np.interp(x, centers, med), 0.0, 1.2),
                centers=centers, k=med, n=len(b), b=b, kk=k)


def save_calibration_plot(cal, band, cloud_v, path):
    try:
        import matplotlib
        matplotlib.use("Agg")
        import matplotlib.pyplot as plt
    except ImportError:
        return
    rng = np.random.default_rng(0)
    sel = rng.choice(len(cal["b"]), min(len(cal["b"]), 5000), replace=False)
    fig, ax = plt.subplots(figsize=(7, 4.5))
    ax.scatter(cal["b"][sel], cal["kk"][sel], s=3, alpha=0.25, label="route point × image")
    ax.plot(cal["centers"], cal["k"], "r-o", lw=2, label="calibration (median, monotone)")
    ax.axvline(cloud_v, color="grey", ls="--", lw=1, label=f"cloud threshold {cloud_v}")
    ax.set_xlabel(f"{band.upper()} brightness"); ax.set_ylabel("clear-sky index  GHI / GHI_clear")
    ax.set_xlim(0, 255); ax.set_ylim(0, 1.4); ax.legend(fontsize=8)
    ax.set_title(f"Brightness → sunshine, learned from {cal['n']} daytime values (last 24 h)")
    fig.tight_layout()
    fig.savefig(path, dpi=120)
    plt.close(fig)


# ---------------------------------------------------------------- borders
def detect_borders(paths, roi, static_thr, contrast_thr, color_rgb=None, tol=40):
    """uint8 mask (255 = coastline/border overlay pixel).
    Decodes the images in a separate pass on purpose: caching all decoded frames would cost
    ~300 MB RAM per band to save about a second of JPEG decoding."""
    if color_rgb is not None:
        img = read(paths[-1], roi)
        bgr = np.array(color_rgb[::-1], dtype=np.int16)
        return ((np.abs(img.astype(np.int16) - bgr).max(axis=2) <= tol)).astype(np.uint8) * 255

    s = ss = None
    n = 0
    for p in paths:
        img = read(p, roi)
        if img is None:
            continue
        f = img.astype(np.float64)
        if s is None:
            s, ss = np.zeros_like(f), np.zeros_like(f)
        if f.shape != s.shape:
            continue
        s += f
        ss += f * f
        n += 1
    if n < 6:
        sys.exit("Too few images to find borders automatically – use --border-color R G B.")
    mean = s / n
    std = np.sqrt(np.maximum(ss / n - mean ** 2, 0)).mean(axis=2)
    mean8 = np.clip(mean, 0, 255).astype(np.uint8)
    contrast = np.abs(mean8.astype(np.int16) - cv2.medianBlur(mean8, 9)).max(axis=2)
    mask = (std < static_thr) & (contrast > contrast_thr)

    # 2nd pass: learn the overlay colour and add static pixels of that colour
    # (catches low-contrast sections, e.g. yellow lines over desert)
    if mask.sum() > 50:
        line_col = np.median(mean[mask], axis=0)
        near = np.abs(mean - line_col).max(axis=2) < 35
        mask |= near & (std < static_thr * 1.5)
    mask = cv2.morphologyEx(mask.astype(np.uint8), cv2.MORPH_CLOSE, np.ones((3, 3), np.uint8))

    # drop isolated specks (JPEG noise)
    k, lab, stats, _ = cv2.connectedComponentsWithStats(mask, connectivity=8)
    keep = np.zeros(k, bool)
    keep[1:] = stats[1:, cv2.CC_STAT_AREA] >= 15
    return keep[lab].astype(np.uint8) * 255


def draw_borders(img, mask, color_bgr, thickness=1):
    m = mask if thickness <= 1 else cv2.dilate(mask, np.ones((thickness, thickness), np.uint8))
    out = img.copy()
    out[m > 0] = color_bgr
    return out


# ---------------------------------------------------------------- preprocessing
def make_cleaner(mask, sigma=3.0, sigma_wide=12.0):
    """Returns clean(img): paints over the overlay lines (incl. JPEG halo) with the
    weighted mean of the surrounding non-line pixels (normalised Gaussian blur).
    Everything that depends only on the mask is computed once here, so cleaning one
    frame costs two blurs instead of a full cv2.inpaint."""
    if not mask.any():
        return lambda img: img
    m = cv2.dilate(mask, np.ones((5, 5), np.uint8)) > 0
    h, wd = m.shape
    w = (~m).astype(np.float32)
    w3 = w[..., None]
    den = cv2.GaussianBlur(w, (0, 0), sigma)[m]
    wide = den < 0.05           # masked pixels with (almost) no valid neighbours nearby
    use_wide = bool(wide.any())
    den = np.maximum(den, 1e-6)[:, None]
    # wide fallback (line junctions, dense coastlines): a coarse fill is enough there,
    # so it is computed at 1/4 resolution, which makes the large kernel cheap
    if use_wide:
        small = (max(1, wd // 4), max(1, h // 4))
        w_s = cv2.resize(w, small, interpolation=cv2.INTER_AREA)
        den_s = np.maximum(cv2.GaussianBlur(w_s, (0, 0), sigma_wide / 4), 1e-6)[..., None]

    def clean(img):
        f = img.astype(np.float32) * w3
        val = cv2.GaussianBlur(f, (0, 0), sigma)[m] / den
        if use_wide:
            f_s = cv2.resize(f, small, interpolation=cv2.INTER_AREA)
            r = cv2.GaussianBlur(f_s, (0, 0), sigma_wide / 4) / den_s
            val[wide] = cv2.resize(r, (wd, h), interpolation=cv2.INTER_LINEAR)[m][wide]
        out = img.copy()
        out[m] = np.clip(val, 0, 255).astype(np.uint8)
        return out
    return clean


def cloud_mask(img, v_thr, s_thr):
    """Clouds are bright and grey/white (B13 is grey anyway; desert in TRM is orange)."""
    hsv = cv2.cvtColor(img, cv2.COLOR_BGR2HSV)
    return (hsv[..., 2] > v_thr) & (hsv[..., 1] < s_thr)


# ---------------------------------------------------------------- motion
_TLS = threading.local()


def flow(a, b):
    """Dense displacement a->b in pixels (DIS optical flow; one instance per thread)."""
    dis = getattr(_TLS, "dis", None)
    if dis is None:
        dis = _TLS.dis = cv2.DISOpticalFlow_create(cv2.DISOPTICAL_FLOW_PRESET_MEDIUM)
    return dis.calc(cv2.GaussianBlur(a, (3, 3), 0), cv2.GaussianBlur(b, (3, 3), 0), None)


def motion_field(recent, k, n_pairs, shape):
    """Motion at frame k: weighted mean of the flow fields of the last n_pairs frames
    (k-n_pairs+1 … k) in px per 10 min, newest weighs most. recent = [(frame index, flow)],
    flows already normalised to 10 min; missing pairs (gaps) are simply absent.
    Returns (field, number of pairs used)."""
    acc, wsum, used = 0, 0.0, 0
    for j, f in recent:
        if k - n_pairs < j <= k:
            wgt = j - (k - n_pairs)
            acc = acc + wgt * f
            wsum += wgt
            used += 1
    if wsum == 0:
        return np.zeros(shape + (2,), np.float32), 0
    return cv2.GaussianBlur((acc / wsum).astype(np.float32), (0, 0), 3), used


def advect_iter(img, v, substeps=10):
    """Semi-Lagrangian: I(x, t+dt) = I(x - v*dt, t); works for grey and colour images.
    Yields one uint8 image per 10-minute step; each step continues from the previous
    (float) state, so n steps cost n remaps instead of 1+2+…+n."""
    h, w = img.shape[:2]
    gx, gy = np.meshgrid(np.arange(w, dtype=np.float32), np.arange(h, dtype=np.float32))
    d = 1.0 / substeps
    mx, my = gx - v[..., 0] * d, gy - v[..., 1] * d
    out = img.astype(np.float32)
    while True:
        for _ in range(substeps):
            out = cv2.remap(out, mx, my, cv2.INTER_LINEAR, borderMode=cv2.BORDER_REPLICATE)
        yield np.clip(out, 0, 255).astype(np.uint8)


def auto_lead(last_obs, now_utc, ahead_min, max_min=120):
    """Lead time (multiple of 10) so the forecast reaches ahead_min past now."""
    latency = (now_utc - last_obs).total_seconds() / 60
    lead = int(np.ceil((latency + ahead_min) / STEP_MIN)) * STEP_MIN
    if lead > max_min:
        print(f"! Newest image is {latency:.0f} min old (stale data?) – capping lead at "
              f"{max_min} min; use --lead to set it yourself.")
        lead = max_min
    return max(STEP_MIN, lead), latency


# ---------------------------------------------------------------- skill
REGIONS = ("", "_north", "_south")
SKILL_FIELDS = (["band", "start_utc", "lead_min", "mae_adv", "mae_pers",
                 "mae_adv_route", "mae_pers_route", "acc_adv_route", "acc_pers_route"]
                + [f"mae_w_{m}{r}" for r in REGIONS for m in ("adv", "pers")]
                + [f"mae_w_nwp_{m}{r}" for m in NWP_MODELS for r in REGIONS]
                + [f"mae_w_blend{r}" for r in REGIONS]
                + [f"mae_w_anc{r}" for r in REGIONS])


def skill_test(band, grays, tlist, idx_by_time, k, v, max_steps, route, cloud_v, irr=None):
    """Forecast frame k up to max_steps × 10 min ahead and compare each step with the
    frame actually observed then (if there is one). With irr (reference + calibration),
    also the irradiance at the route points against the reference, compared with smart
    persistence (today's clear-sky index kept constant). Returns rows for skill_history.csv."""
    rows, grams = [], {}
    g0 = grays[k]
    p0 = g0.astype(np.float32)
    rv0 = route_values(g0, route) if route else None
    G0 = irr["ref"].get(tlist[k]) if irr else None
    if G0 is not None:
        gcs0, cz0 = clear_sky_ghi(route["lat"], route["lon"], tlist[k])
        k0 = np.clip(G0 / np.maximum(gcs0, 1.0), 0, 1.5)      # clear-sky index at the start
    it = advect_iter(g0, v)
    for st in range(1, max_steps + 1):
        pred = next(it)
        j = idx_by_time.get(tlist[k] + dt.timedelta(minutes=STEP_MIN * st))
        if j is None:
            continue
        truth = grays[j].astype(np.float32)
        ea = np.abs(pred.astype(np.float32) - truth)
        ep = np.abs(p0 - truth)
        row = {f: "" for f in SKILL_FIELDS}
        row.update(band=band, start_utc=tlist[k].strftime("%Y-%m-%d %H:%M"),
                   lead_min=st * STEP_MIN, mae_adv=round(float(ea.mean()), 3),
                   mae_pers=round(float(ep.mean()), 3))
        if route is not None and route["corridor"].any():
            c = route["corridor"]
            ok = route["inside"]
            rv_pred = route_values(pred, route)
            rt = route_values(grays[j], route)[ok] > cloud_v
            ra = rv_pred[ok] > cloud_v
            row.update(mae_adv_route=round(float(ea[c].mean()), 3),
                       mae_pers_route=round(float(ep[c].mean()), 3),
                       acc_adv_route=round(float((ra == rt).mean()), 4),
                       acc_pers_route=round(float(((rv0[ok] > cloud_v) == rt).mean()), 4))
            G = irr["ref"].get(tlist[j]) if G0 is not None else None
            if G is not None:
                gcs, cz = clear_sky_ghi(route["lat"], route["lon"], tlist[j])
                m = ok & (cz > MIN_COSZ) & (cz0 > MIN_COSZ) & np.isfinite(G) & np.isfinite(G0)
                # weather-model forecasts for the same target time (if fetched); all three
                # forecasts are scored on the same points, so they stay comparable
                w_nwp, p_nwp = {}, None
                for name, fc in irr.get("nwp", {}).items():
                    if name in NWP_MODELS:
                        kn = fc.k_at(tlist[j])
                        w_nwp[name] = np.abs(gcs * kn - G)
                        if name == irr.get("nwp_primary"):
                            p_nwp = gcs * kn
                        m = m & np.isfinite(kn)
                p_sat = gcs * irr["k_of"](rv_pred)                  # advected cloud -> W/m²
                p_per = gcs * k0                                    # smart persistence
                # anchored: today's measured value + the change the moving clouds bring
                p_anc = gcs * anchored_k(k0, irr["k_of"], rv0, rv_pred)
                w_adv, w_per = np.abs(p_sat - G), np.abs(p_per - G)
                w_anc = np.abs(p_anc - G)
                lead = st * STEP_MIN
                for tag, sel in zip(REGIONS, (m, m & irr["north"], m & ~irr["north"])):
                    if sel.any():
                        row[f"mae_w_adv{tag}"] = round(float(w_adv[sel].mean()), 2)
                        row[f"mae_w_pers{tag}"] = round(float(w_per[sel].mean()), 2)
                        row[f"mae_w_anc{tag}"] = round(float(w_anc[sel].mean()), 2)
                        for name, w in w_nwp.items():
                            row[f"mae_w_nwp_{name}{tag}"] = round(float(w[sel].mean()), 2)
                        if p_nwp is not None:
                            X = np.stack([p_per[sel], p_sat[sel], p_anc[sel], p_nwp[sel]])
                            grams.setdefault(f"{row['start_utc']}|{lead}", {})[tag] = \
                                gram_of(X, G[sel])
                            wb = pick_weights(irr.get("blend_w"), tag, lead)
                            if wb is not None:   # weights learned on earlier days only
                                pb = wb @ X
                                row[f"mae_w_blend{tag}"] = round(float(np.abs(pb - G[sel]).mean()), 2)
        rows.append(row)
    return rows, grams


# ---------------------------------------------------------------- blend
# Blend = w_now · (cloudiness stays) + w_sat · (our satellite forecast) + w_model · (weather model),
# weights ≥ 0 that sum to 1, one set per lead time and route part. They are the least-squares
# best weights over the retrospective tests, collected in blend_stats.json as small sums
# (Gram matrices), so they improve with every day of data.
BLEND_PREDICTORS = ("now", "satellite", "anchored", "model")
NB = len(BLEND_PREDICTORS)


def anchored_k(k_now, k_of, b_now, b_pred):
    """'Now + change': the measured clear-sky index now, plus the change of the calibrated
    index between the brightness now and the forecast brightness at the same place. The
    unsure brightness -> sunshine conversion then only enters as a difference."""
    return np.clip(k_now + k_of(b_pred) - k_of(b_now), 0.0, 1.2)
BLEND_MIN_STARTS = 10      # weights for a lead/region only with at least this many tests


def gram_of(X, y):
    """[n, Σy², Σx·y (NB), Σx·xᵀ (NB²)] for predictors X (NB × n) and truth y (n)."""
    return [int(X.shape[1]), float(y @ y)] + [float(v) for v in X @ y] + \
        [float(v) for v in (X @ X.T).ravel()]


def solve_blend(A, b, c):
    """Weights w ≥ 0 with Σw = 1 minimising Σ(w·x − y)² = wᵀAw − 2bᵀw + c (all subsets tried)."""
    k_all = len(b)
    A = A + np.eye(k_all) * (1e-6 * np.trace(A) / k_all + 1e-9)
    best = None
    for mask in range(1, 2 ** k_all):
        idx = [i for i in range(k_all) if mask >> i & 1]
        k = len(idx)
        K = np.zeros((k + 1, k + 1))
        K[:k, :k] = 2 * A[np.ix_(idx, idx)]
        K[:k, k] = K[k, :k] = 1
        sol = np.linalg.lstsq(K, np.r_[2 * b[idx], 1.0], rcond=None)[0]
        if (sol[:k] < -1e-9).any():
            continue
        w = np.zeros(k_all)
        w[idx] = np.clip(sol[:k], 0, None)
        w /= w.sum()
        f = w @ A @ w - 2 * b @ w + c
        if best is None or f < best[0]:
            best = (f, w)
    return best[1]


def history_filter(a, now):
    """Which test start times (string 'YYYY-MM-DD HH:MM') count for the history and the blend
    weights: --history-months (e.g. 8 or 7,8,9: only race-season weather) and/or
    --history-days (only the last N days). Returns (keep function, description or None)."""
    months = {int(m) for m in str(a.history_months or "").replace(" ", "").split(",") if m}
    oldest = (now - dt.timedelta(days=a.history_days)).strftime("%Y-%m-%d %H:%M") \
        if a.history_days else None
    if not months and not oldest:
        return (lambda start: True), None
    desc = ", ".join(x for x in (f"months {sorted(months)}" if months else "",
                                 f"last {a.history_days} days" if oldest else "") if x)
    return (lambda start: (not months or int(start[5:7]) in months)
            and (not oldest or start >= oldest)), desc


def weights_from_stats(tests, before=None, keep=None):
    """{region tag: {lead: weights (3,)}} from the collected tests; before = only tests that
    started before this time (so a run is scored with weights it did not learn from)."""
    acc = collections.defaultdict(lambda: [0, 0.0, np.zeros(NB), np.zeros((NB, NB)), set()])
    for key, per_tag in tests.items():
        start, lead = key.rsplit("|", 1)
        if keep is not None and not keep(start):
            continue
        if before is not None and start >= before.strftime("%Y-%m-%d %H:%M"):
            continue
        for tag, g in per_tag.items():
            a = acc[(tag, int(lead))]
            a[0] += g[0]
            a[1] += g[1]
            if len(g) != 2 + NB + NB * NB:
                continue
            a[2] += np.array(g[2:2 + NB])
            a[3] += np.array(g[2 + NB:]).reshape(NB, NB)
            a[4].add(start)
    out = collections.defaultdict(dict)
    for (tag, lead), (n, c, b, A, starts) in acc.items():
        if len(starts) >= BLEND_MIN_STARTS and n > 0:
            out[tag][lead] = solve_blend(A, b, c)
    return dict(out)


def pick_weights(W, tag, lead):
    """Weights for this route part and lead; falls back to the whole route, then to the
    nearest lead time that has weights."""
    if not W:
        return None
    for t in (tag, ""):
        w = W.get(t, {}).get(lead)
        if w is not None:
            return w
    leads = sorted(W.get("", {}))
    return W[""][min(leads, key=lambda L: abs(L - lead))] if leads else None


def latest_now(irr, route, grays, tlist, idx_by_time, t_obs):
    """The newest measured clear-sky index at the route points (reference at most 60 min
    before the newest image) and the image brightness at that same time.
    Returns (k_now, brightness) or (None, None)."""
    past = [t for t in irr["ref"] if t <= t_obs and t in idx_by_time]
    if not past or (t_obs - max(past)).total_seconds() > 3600:
        return None, None
    tr = max(past)
    gcs0, cz0 = clear_sky_ghi(route["lat"], route["lon"], tr)
    k0 = np.where(cz0 > MIN_COSZ, np.clip(irr["ref"][tr] / np.maximum(gcs0, 1.0), 0, 1.5), np.nan)
    return k0, route_values(grays[idx_by_time[tr]], route)


def route_blend(W, irr, rtimes, gclear, k_now, ghi_sat, ghi_anc, ghi_nwp):
    """Blended sunshine forecast per route point and time step (predictors in the order of
    BLEND_PREDICTORS). A predictor that is missing at a point (e.g. the model at dawn, no
    measured value) is dropped there and the others are re-weighted."""
    if k_now is None:
        k_now = np.full(len(irr["north"]), np.nan)
    out = []
    for s in range(len(rtimes)):
        lead = max(STEP_MIN, s * STEP_MIN)
        wn, ws = pick_weights(W, "_north", lead), pick_weights(W, "_south", lead)
        w = np.where(irr["north"][:, None], wn, ws)
        P = np.stack([gclear[s] * k_now, ghi_sat[s], ghi_anc[s], ghi_nwp[s]], axis=1)
        avail = np.isfinite(P)
        wv = np.where(avail, w, 0.0)
        sw = wv.sum(1)
        val = (np.where(avail, P, 0.0) * wv).sum(1) / np.maximum(sw, 1e-9)
        out.append(np.where(sw > 0, val, ghi_sat[s]))
    return out


def load_blend_stats(path, model):
    try:
        d = json.load(open(path))
        if d.get("model") == model and list(d.get("predictors", [])) == list(BLEND_PREDICTORS):
            return d.get("tests", {})
        print("  blend_stats.json was made for other blend parts – starting a new collection")
    except (OSError, ValueError):
        pass
    return {}


def save_blend_stats(path, model, tests):
    try:
        with open(path, "w") as f:
            json.dump({"model": model, "predictors": BLEND_PREDICTORS, "tests": tests}, f)
    except OSError:
        pass


def summarize_skill(rows):
    """Per lead: number of tests, mean MAE, skill = 1 - MAE_adv / MAE_pers, route values."""
    by = collections.defaultdict(list)
    for r in rows:
        by[int(r["lead_min"])].append(r)
    out = []
    for lead in sorted(by):
        rs = by[lead]
        adv = np.mean([float(r["mae_adv"]) for r in rs])
        pers = np.mean([float(r["mae_pers"]) for r in rs])
        rr = [r for r in rs if r["mae_adv_route"] not in ("", None)]
        s = dict(lead_min=lead, n_tests=len(rs), mae_adv=round(adv, 3), mae_pers=round(pers, 3),
                 skill=round(1 - adv / pers, 4) if pers > 0 else "")
        if rr:
            ar = np.mean([float(r["mae_adv_route"]) for r in rr])
            pr = np.mean([float(r["mae_pers_route"]) for r in rr])
            s.update(skill_route=round(1 - ar / pr, 4) if pr > 0 else "",
                     acc_adv_route=round(np.mean([float(r["acc_adv_route"]) for r in rr]), 4),
                     acc_pers_route=round(np.mean([float(r["acc_pers_route"]) for r in rr]), 4))
        for tag in REGIONS:
            rw = [r for r in rs if r.get(f"mae_w_adv{tag}") not in ("", None)]
            if rw:
                aw = np.mean([float(r[f"mae_w_adv{tag}"]) for r in rw])
                pw = np.mean([float(r[f"mae_w_pers{tag}"]) for r in rw])
                s.update({f"n_w{tag}": len(rw), f"mae_w_adv{tag}": round(aw, 1),
                          f"mae_w_pers{tag}": round(pw, 1),
                          f"skill_w{tag}": round(1 - aw / pw, 4) if pw > 0 else ""})
            for name in NWP_MODELS:
                rn = [r for r in rs if r.get(f"mae_w_nwp_{name}{tag}") not in ("", None)]
                if not rn:
                    continue
                nw = np.mean([float(r[f"mae_w_nwp_{name}{tag}"]) for r in rn])
                pw = np.mean([float(r[f"mae_w_pers{tag}"]) for r in rn])
                aw = np.mean([float(r[f"mae_w_adv{tag}"]) for r in rn])   # tool, same tests
                ra = [float(r[f"mae_w_anc{tag}"]) for r in rn
                      if r.get(f"mae_w_anc{tag}") not in ("", None)]
                if len(ra) == len(rn):
                    s[f"mae_w_anc_vs_{name}{tag}"] = round(float(np.mean(ra)), 1)
                s.update({f"n_nwp_{name}{tag}": len(rn),
                          f"mae_w_nwp_{name}{tag}": round(nw, 1),
                          f"mae_w_adv_vs_{name}{tag}": round(aw, 1),
                          f"mae_w_pers_vs_{name}{tag}": round(pw, 1),
                          f"skill_w_nwp_{name}{tag}": round(1 - nw / pw, 4) if pw > 0 else ""})
            ra = [r for r in rs if r.get(f"mae_w_anc{tag}") not in ("", None)]
            if ra:
                nw = np.mean([float(r[f"mae_w_anc{tag}"]) for r in ra])
                pw = np.mean([float(r[f"mae_w_pers{tag}"]) for r in ra])
                s.update({f"mae_w_anc{tag}": round(nw, 1),
                          f"skill_w_anc{tag}": round(1 - nw / pw, 4) if pw > 0 else ""})
            rb = [r for r in rs if r.get(f"mae_w_blend{tag}") not in ("", None)]
            if rb:
                bw = np.mean([float(r[f"mae_w_blend{tag}"]) for r in rb])
                pw = np.mean([float(r[f"mae_w_pers{tag}"]) for r in rb])
                aw = np.mean([float(r[f"mae_w_adv{tag}"]) for r in rb])
                s.update({f"n_blend{tag}": len(rb), f"mae_w_blend{tag}": round(bw, 1),
                          f"mae_w_adv_vs_blend{tag}": round(aw, 1),
                          f"mae_w_pers_vs_blend{tag}": round(pw, 1),
                          f"skill_w_blend{tag}": round(1 - bw / pw, 4) if pw > 0 else ""})
        out.append(s)
    return out


def save_nwp_plot(data, band, n_days, path):
    """Our forecast vs. weather models (Open-Meteo) vs. smart persistence: mean error in W/m²
    and skill, per lead time, whole route. data = summarize_skill output (history if present)."""
    try:
        import matplotlib
        matplotlib.use("Agg")
        import matplotlib.pyplot as plt
    except ImportError:
        return False
    names = [m for m in NWP_MODELS if any(s.get(f"mae_w_nwp_{m}") not in ("", None) for s in data)]
    if not names:
        return False
    fig, ax = plt.subplots(1, 2, figsize=(13, 4.8))
    ref = names[0]   # tool and persistence shown on the tests the first model was scored on
    pts = [s for s in data if s.get(f"mae_w_nwp_{ref}") not in ("", None)]
    L = [s["lead_min"] for s in pts]
    ax[0].plot(L, [s[f"mae_w_adv_vs_{ref}"] for s in pts], "-o", lw=2.5, color="tab:blue",
               label="our forecast (satellite)")
    ax[0].plot(L, [s[f"mae_w_pers_vs_{ref}"] for s in pts], "-o", color="tab:gray",
               label="'cloudiness stays' (smart persistence)")
    sk = [(s["lead_min"], 100 * s["skill_w"]) for s in data if s.get("skill_w") not in ("", None)]
    if sk:
        ax[1].plot(*zip(*sk), "-o", lw=2.5, color="tab:blue", label="our forecast (satellite)")
    pa = [s for s in pts if s.get(f"mae_w_anc_vs_{ref}") not in ("", None)]
    lab_a = "now + cloud change (satellite)"
    if pa:
        ax[0].plot([s["lead_min"] for s in pa], [s[f"mae_w_anc_vs_{ref}"] for s in pa], "-P",
                   ms=8, lw=2.5, color="magenta", label=lab_a)
    ska = [(s["lead_min"], 100 * s["skill_w_anc"]) for s in data
           if s.get("skill_w_anc") not in ("", None)]
    if ska:
        ax[1].plot(*zip(*ska), "-P", ms=8, lw=2.5, color="magenta", label=lab_a)
    styles = [("--s", "tab:orange"), ("--^", "tab:green"), ("--D", "tab:red"), ("--v", "tab:brown")]
    for name, (st, col) in zip(names, styles):
        pn = [s for s in data if s.get(f"mae_w_nwp_{name}") not in ("", None)]
        ax[0].plot([s["lead_min"] for s in pn], [s[f"mae_w_nwp_{name}"] for s in pn], st,
                   color=col, label=f"weather model: {name}")
        ax[1].plot([s["lead_min"] for s in pn], [100 * s[f"skill_w_nwp_{name}"] for s in pn], st,
                   color=col, label=f"weather model: {name}")
    pb = [s for s in data if s.get("mae_w_blend") not in ("", None)]
    if pb:
        lab = "blend (now + both satellite forecasts + model)"
        ax[0].plot([s["lead_min"] for s in pb], [s["mae_w_blend"] for s in pb], "-*", ms=10,
                   lw=2.5, color="black", label=lab)
        ax[1].plot([s["lead_min"] for s in pb], [100 * s["skill_w_blend"] for s in pb], "-*",
                   ms=10, lw=2.5, color="black", label=lab)
    ax[0].set_xlabel("lead time [min]"); ax[0].set_ylabel("mean error [W/m²]  (lower = better)")
    ax[0].set_title("Error of the sunshine forecast at the route points")
    ax[0].set_ylim(bottom=0)
    ax[1].axhline(0, color="k", lw=0.8)
    ax[1].set_xlabel("lead time [min]"); ax[1].set_ylabel("skill vs. smart persistence [%]")
    ax[1].set_title("Better than 'cloudiness stays'? (above 0 = yes)")
    fig.legend(*ax[0].get_legend_handles_labels(), loc="lower center", ncol=3, fontsize=9,
               frameon=False)
    fig.suptitle(f"Himawari {band.upper()} – our forecast vs. weather models "
                 f"({n_days} day{'s' if n_days != 1 else ''})")
    fig.tight_layout(rect=(0, 0.13, 1, 1))
    fig.savefig(path, dpi=120)
    plt.close(fig)
    return True


def update_history(path, rows):
    """Add new tests to skill_history.csv (one row per band/start/lead, no duplicates).
    The file is rewritten with the current columns, so older files gain new columns."""
    old = []
    if os.path.exists(path):
        with open(path, newline="") as f:
            old = list(csv.DictReader(f))
    key = lambda r: (r["band"], r["start_utc"], str(r["lead_min"]))
    filled = lambda r: sum(1 for f in SKILL_FIELDS if r.get(f) not in ("", None))
    merged = {key(r): r for r in old}
    for r in rows:
        # a test already in the history is replaced if the new one has more values
        # (e.g. the earlier run had no route or no irradiance reference)
        k = key(r)
        if k not in merged or filled(r) > filled(merged[k]):
            merged[k] = r
    out = sorted(merged.values(), key=lambda r: (r["band"], r["start_utc"], int(r["lead_min"])))
    with open(path, "w", newline="") as f:
        wr = csv.DictWriter(f, fieldnames=SKILL_FIELDS, restval="", extrasaction="ignore")
        wr.writeheader()
        wr.writerows(out)
    return out


def save_skill_plot(run, hist, band, n_days, outdir):
    try:
        import matplotlib
        matplotlib.use("Agg")
        import matplotlib.pyplot as plt
    except ImportError:
        return
    fig, ax = plt.subplots(1, 3, figsize=(17, 4.8))

    def series(data, key):
        """(leads, values in %) for the leads that have this value – the history can mix
        runs with and without route/irradiance, so a value may exist only for some leads."""
        pts = [(s["lead_min"], 100 * float(s[key])) for s in data
               if s.get(key) not in ("", None)]
        return tuple(zip(*pts)) if pts else None

    for data, style, name in ((run, "-o", "this run"), (hist, "--s", f"history ({n_days} days)")):
        if not data:
            continue
        # whole image drawn last and thicker, so it stays visible if the route lies on top
        for axi, key, kw in ((0, "skill_route", dict(ms=4, label=f"route corridor, {name}")),
                             (1, "acc_adv_route", dict(label=f"advection, {name}")),
                             (1, "acc_pers_route", dict(alpha=0.6, label=f"persistence, {name}")),
                             (0, "skill", dict(lw=2.5, mfc="none", ms=8, zorder=3,
                                               label=f"whole image, {name}"))):
            xy = series(data, key)
            if xy:
                ax[axi].plot(*xy, style, **kw)
        for tag, lab in zip(REGIONS, ("whole route", "north of 23.4°S", "south of 23.4°S")):
            xy = series(data, f"skill_w{tag}")
            if xy:
                ax[2].plot(*xy, style, lw=2.5 if not tag else 1.3, label=f"{lab}, {name}")
        xy = series(data, "skill_w_anc")
        if xy:
            ax[2].plot(*xy, style[:-1] + "P", ms=8, lw=2.5, color="magenta",
                       label=f"now + cloud change (whole route), {name}")
    ax[2].axhline(0, color="k", lw=0.8)
    ax[2].set_xlabel("lead time [min]"); ax[2].set_ylabel("skill vs. smart persistence [%]")
    ax[2].set_title("Sunshine (W/m²) at the route points")
    if ax[2].lines[1:]:
        ax[2].legend(fontsize=8)
    else:
        ax[2].text(0.5, 0.5, "no irradiance reference\n(route + Open-Meteo needed)",
                   ha="center", va="center", transform=ax[2].transAxes)
    ax[0].axhline(0, color="k", lw=0.8)
    ax[0].set_xlabel("lead time [min]"); ax[0].set_ylabel("skill vs. persistence [%]")
    ax[0].set_title("Error reduction vs. 'nothing changes'"); ax[0].legend(fontsize=8)
    ax[1].set_xlabel("lead time [min]"); ax[1].set_ylabel("correct cloud / no cloud [%]")
    ax[1].set_title("Route points classified correctly")
    if ax[1].lines:
        ax[1].legend(fontsize=8)
    else:
        ax[1].text(0.5, 0.5, "no route given (--route)", ha="center", transform=ax[1].transAxes)
    fig.suptitle(f"Himawari {band.upper()} – forecast quality by lead time")
    fig.tight_layout()
    fig.savefig(os.path.join(outdir, "skill_curve.png"), dpi=120)
    plt.close(fig)


# ---------------------------------------------------------------- output helpers
def overlay_vectors(img, v, step=24, scale=3):
    vis = img.copy()
    for y in range(step // 2, img.shape[0], step):
        for x in range(step // 2, img.shape[1], step):
            dx, dy = v[y, x] * scale
            if np.hypot(dx, dy) > 1:
                cv2.arrowedLine(vis, (x, y), (int(x + dx), int(y + dy)), (0, 140, 255), 1,
                                tipLength=0.3)
    return vis


def label(img, text, color, org=(10, 30)):
    out = img.copy()
    cv2.putText(out, text, org, cv2.FONT_HERSHEY_SIMPLEX, 0.7, (0, 0, 0), 4)
    cv2.putText(out, text, org, cv2.FONT_HERSHEY_SIMPLEX, 0.7, color, 2)
    return out


def change_map(prev_img, next_img, v):
    """Red = brighter (cloud arriving), blue = darker (cloud leaving), with motion arrows."""
    g0 = cv2.cvtColor(prev_img, cv2.COLOR_BGR2GRAY)
    g1 = cv2.cvtColor(next_img, cv2.COLOR_BGR2GRAY)
    diff = g1.astype(np.int16) - g0.astype(np.int16)
    vis = cv2.cvtColor(g1, cv2.COLOR_GRAY2BGR).astype(np.float32) * 0.5
    vis[..., 2] += np.clip(diff, 0, None) * 3
    vis[..., 0] += np.clip(-diff, 0, None) * 3
    return overlay_vectors(np.clip(vis, 0, 255).astype(np.uint8), v), float(np.abs(diff).mean())


def save_motion_plot(rows, band, outdir):
    try:
        import matplotlib
        matplotlib.use("Agg")
        import matplotlib.pyplot as plt
    except ImportError:
        return
    t = [dt.datetime.strptime(r["time_utc"], "%Y-%m-%d %H:%M") for r in rows]
    fig, ax = plt.subplots(3, 1, figsize=(10, 8), sharex=True)
    ax[0].plot(t, [r["cloud_fraction"] for r in rows]); ax[0].set_ylabel("Cloud fraction")
    ax[1].plot(t, [r["median_speed_px_per_10min"] for r in rows], label="median")
    ax[1].plot(t, [r["p90_speed_px_per_10min"] for r in rows], label="P90")
    ax[1].set_ylabel("px / 10 min"); ax[1].legend()
    ax[2].scatter(t, [r["heading_deg"] for r in rows], s=6)
    ax[2].set_ylabel("Heading [°]"); ax[2].set_ylim(0, 360); ax[2].set_xlabel("UTC")
    fig.suptitle(f"Himawari {band.upper()} – cloud motion over the last 24 h")
    fig.tight_layout()
    fig.savefig(os.path.join(outdir, "motion_24h.png"), dpi=120)
    plt.close(fig)


def car_km_at(t, a, now):
    """Expected car position (km along the route) at time t, or None without --car-km."""
    if a.car_km is None:
        return None
    return a.car_km + a.speed * (t - now).total_seconds() / 3600


def route_map(img, bmask, route, georef, roi, values, cloud_v, title, car_km):
    """Route on the image, with the detected borders and a 5° lat/lon grid, so the map
    position can be checked by eye. Route points: red = cloud, green = clear."""
    vis = draw_borders(img, bmask, (0, 255, 255), 1)
    h, w = vis.shape[:2]
    ox, oy = (roi[0], roi[1]) if roi else (0, 0)
    lon0, lon1, lat0, lat1 = AUS_EXTENT
    grid = (90, 90, 90)
    for lon in range(int(lon0), int(lon1) + 1, 5):
        lats = np.linspace(lat0 - 5, lat1 + 5, 200)
        x, y = georef.to_px(lats, np.full_like(lats, lon))
        pl = np.stack([x - ox, y - oy], 1).astype(np.int32)
        cv2.polylines(vis, [pl], False, grid, 1)
        xm, ym = georef.to_px(np.array([lat1]), np.array([float(lon)]))
        cv2.putText(vis, f"{lon}E", (int(xm[0] - ox) + 3, int(ym[0] - oy) + 15),
                    cv2.FONT_HERSHEY_SIMPLEX, 0.4, (200, 200, 200), 1)
    for lat in range(int(lat0), int(lat1) + 1, 5):
        lons = np.linspace(lon0 - 5, lon1 + 5, 200)
        x, y = georef.to_px(np.full_like(lons, lat), lons)
        pl = np.stack([x - ox, y - oy], 1).astype(np.int32)
        cv2.polylines(vis, [pl], False, grid, 1)
        xm, ym = georef.to_px(np.array([float(lat)]), np.array([lon0]))
        cv2.putText(vis, f"{-lat}S", (int(xm[0] - ox) + 3, int(ym[0] - oy) - 3),
                    cv2.FONT_HERSHEY_SIMPLEX, 0.4, (200, 200, 200), 1)

    ok = route["inside"]
    pts = np.stack([route["ix"], route["iy"]], 1)
    if ok.sum() >= 2:
        cv2.polylines(vis, [pts[ok].astype(np.int32)], False, (255, 255, 255), 3)
        cv2.polylines(vis, [pts[ok].astype(np.int32)], False, (60, 60, 60), 1)
    for i in np.flatnonzero(ok):
        cloudy = values[i] > cloud_v
        cv2.circle(vis, tuple(int(c) for c in pts[i]), 4,
                   (40, 40, 230) if cloudy else (60, 200, 60), -1)
    idx = np.flatnonzero(ok)
    if len(idx):
        marks = sorted(set([idx[0], idx[-1]] + [i for i in idx if i % 25 == 0]))
        for i in marks:
            txt = f"{route['name'][i]} ({route['km'][i]:.0f} km)"
            p = (int(pts[i][0]) + 9, int(pts[i][1]) + 5)
            (tw, th), base = cv2.getTextSize(txt, cv2.FONT_HERSHEY_SIMPLEX, 0.45, 1)
            cv2.rectangle(vis, (p[0] - 3, p[1] - th - 3), (p[0] + tw + 3, p[1] + base + 1),
                          (30, 30, 30), -1)
            cv2.putText(vis, txt, p, cv2.FONT_HERSHEY_SIMPLEX, 0.45, (255, 255, 255), 1,
                        cv2.LINE_AA)
    if car_km is not None:
        lat = np.interp(car_km, route["km"], route["lat"])
        lon = np.interp(car_km, route["km"], route["lon"])
        cx, cy = georef.to_px(np.array([lat]), np.array([lon]))
        c = (int(cx[0] - ox), int(cy[0] - oy))
        if 0 <= c[0] < w and 0 <= c[1] < h:
            cv2.drawMarker(vis, c, (255, 0, 255), cv2.MARKER_STAR, 22, 2)
            cv2.putText(vis, f"car now ({car_km:.0f} km)", (c[0] - 150, c[1] + 5),
                        cv2.FONT_HERSHEY_SIMPLEX, 0.55, (255, 0, 255), 2, cv2.LINE_AA)
    vis = label(vis, title, (255, 255, 255))
    return label(vis, "red = cloud, green = clear", (255, 255, 255), (10, 58))


def save_route_plot(route, times, values, cloud_v, now, a, band, path, ghi=None,
                    ghi_title="Expected sunshine: global horizontal irradiance"):
    """km × time chart of the forecast along the route: brightness, and (with a calibration)
    the expected global irradiance in W/m²."""
    try:
        import matplotlib
        matplotlib.use("Agg")
        import matplotlib.dates as mdates
        import matplotlib.pyplot as plt
    except ImportError:
        return
    ok = route["inside"]
    km = route["km"][ok]
    tnum = mdates.date2num(times)
    half = (tnum[1] - tnum[0]) / 2 if len(tnum) > 1 else 0.004
    ext = [km[0], km[-1], tnum[0] - half, tnum[-1] + half]
    panels = [("brightness", np.array(values)[:, ok], "Greys_r", 0, 255,
               f"{band.upper()} brightness (bright = cloud, red line = threshold {cloud_v})")]
    if ghi is not None:
        panels.append(("GHI [W/m²]", np.array(ghi)[:, ok], "inferno", 0, 1100, ghi_title))
    fig, axs = plt.subplots(len(panels), 1, figsize=(12, 4.6 * len(panels)), sharex=True,
                            squeeze=False)
    for ax, (cb, Z, cmap, vmin, vmax, title) in zip(axs[:, 0], panels):
        im = ax.imshow(Z, aspect="auto", origin="lower", cmap=cmap, vmin=vmin, vmax=vmax,
                       extent=ext, interpolation="nearest")
        if cb == "brightness" and Z.shape[0] > 1 and Z.shape[1] > 1:
            ax.contour(km, tnum, Z, levels=[cloud_v], colors="tab:red", linewidths=1)
        ax.axhline(mdates.date2num(now), color="tab:cyan", lw=1.5, label="now")
        if a.car_km is not None:
            tt = [t for t in times if t >= now] or [now]
            ax.plot([car_km_at(t, a, now) for t in tt], mdates.date2num(tt), color="magenta",
                    lw=2.5, label=f"car ({a.speed:.0f} km/h)")
            ax.set_xlim(km[0], km[-1])
        ax.yaxis_date()
        ax.yaxis.set_major_formatter(mdates.DateFormatter("%H:%M"))
        ax.set_ylabel("time [UTC]")
        ax.set_title(f"{title}; first row = newest image")
        fig.colorbar(im, ax=ax, label=cb)
        ax.legend(loc="upper right")
    axs[-1, 0].set_xlabel("km along route")
    fig.tight_layout()
    fig.savefig(path, dpi=120)
    plt.close(fig)


# ---------------------------------------------------------------- one band
def run_band(band, a, now, route_pts):
    cfg = BANDS[band]
    data_dir = os.path.join(a.data_root, f"himawari_{band}")
    out_dir = os.path.join(a.out_root, f"result_{band}")
    cloud_v = a.cloud_v if a.cloud_v is not None else cfg["cloud_v"]
    keep_hist, hist_desc = history_filter(a, now)
    print(f"\n===== {band.upper()} =====")
    timing, t_mark = collections.OrderedDict(), [time.perf_counter()]

    def tick(name):          # time spent since the last tick
        t = time.perf_counter()
        timing[name] = timing.get(name, 0) + t - t_mark[0]
        t_mark[0] = t

    if a.no_download:
        paths, times = list_frames(band, data_dir)   # read only, never written
    else:
        paths, times = download(band, now)           # in memory only
    tick("download")
    if len(paths) < a.pairs + 2:
        print(f"Too few {band} images – skipped.")
        return False

    # keep only what this run produces: remove results of earlier runs
    os.makedirs(out_dir, exist_ok=True)
    for f in os.listdir(out_dir):
        if f.endswith((".png", ".gif", ".csv")) and f.startswith(OUTPUT_PREFIXES):
            os.remove(os.path.join(out_dir, f))
    line_bgr = tuple(int(c) for c in a.line_color[::-1])

    # 1) borders: detect and remove before measuring motion
    bmask = detect_borders(paths, a.roi, a.static_thr, a.contrast_thr, a.border_color, a.border_tol)
    cv2.imwrite(os.path.join(out_dir, "border_mask.png"), bmask)
    print(f"Borders: {int((bmask > 0).sum())} pixels detected and removed")
    clean = make_cleaner(bmask)
    tick("borders")

    def B(img):
        return draw_borders(img, bmask, line_bgr, a.line_width) if a.draw_borders else img

    # route: map position (lat/lon -> pixel) in this band's images
    route = georef = None
    if route_pts:
        full = read(paths[-1], None).shape[:2]
        georef = Georef(json.load(open(a.georef)) if a.georef else {}, full[1], full[0])
        route = build_route(route_pts, georef, a.roi, bmask.shape, a.route_radius)
        n_in = int(route["inside"].sum())
        print(f"Route: {len(route_pts)} points, {route['km'][-1]:.0f} km, {n_in} inside the image "
              f"– map position: {georef.desc}")
        if n_in < 2:
            print("! Route is (almost) outside the image – check --georef / --roi. Route skipped.")
            route = None

    workers = max(1, a.workers)
    chunk = 2 * workers
    with cf.ThreadPoolExecutor(workers) as ex:
        # 2a) decode + clean every frame in parallel; keep grey image + packed cloud mask
        def prep(src):
            img = read(src, a.roi)
            if img is None or img.shape[:2] != bmask.shape:
                return None
            c = clean(img)
            return c, cv2.cvtColor(c, cv2.COLOR_BGR2GRAY), np.packbits(cloud_mask(c, cloud_v, a.cloud_s))

        grays, tlist, clouds = [], [], []
        last = None
        for res, t in zip(chunked_map(ex, prep, paths, chunk), times):
            if res is None:
                continue
            last, g, cl = res          # only the newest colour frame is kept
            grays.append(g)
            tlist.append(t)
            clouds.append(cl)
        if len(grays) < a.pairs + 2:
            print(f"Too few readable {band} images – skipped.")
            return False
        print(f"{len(grays)} images, {tlist[0]:%d.%m %H:%M} – {tlist[-1]:%d.%m %H:%M} UTC")
        idx_by_time = {t: i for i, t in enumerate(tlist)}
        tick("decode+clean")

        # irradiance: reference along the route, and brightness -> sunshine calibration
        irr = cal = None
        blend_tests, blend_path = {}, None
        if route is not None and not a.no_irradiance:
            ref = fetch_reference(route_pts, tlist[0], tlist[-1])
            if ref:
                cal = fit_calibration(grays, tlist, route, ref)
            if cal:
                irr = dict(ref=ref, k_of=cal["k_of"], north=route["lat"] > NORTH_LAT)
                kb = {int(b): float(cal["k_of"](b)) for b in (60, 120, cloud_v, 220)}
                print(f"Calibration from {cal['n']} daytime values: clear-sky index "
                      + ", ".join(f"{v:.2f} at brightness {b}" for b, v in kb.items()))
                save_calibration_plot(cal, band, cloud_v,
                                      os.path.join(out_dir, "route_calibration.png"))
                # weather-model forecasts for the same points: last 24 h (to compare) + ahead
                models = [m.strip() for m in a.nwp_models.split(",") if m.strip()]
                if models:
                    t_end = tlist[-1] + dt.timedelta(minutes=a.route_lead + 60)
                    irr["nwp"] = fetch_nwp(route_pts, tlist[0], t_end, models, a.nwp_source,
                                           os.path.join(a.out_root, "openmeteo_nwp_cache.json"),
                                           a.nwp_cache_min,
                                           None if a.no_bias_correction else BIAS_FACTOR)
                    prim = next((m for m in models if m in irr["nwp"] and m in NWP_MODELS), None)
                    if prim:
                        # blend: weights learned from tests that started before this 24-h
                        # window, so the blend is scored on data it has not seen
                        irr["nwp_primary"] = prim
                        # corrected and raw model values must not be mixed in the blend statistics
                        fac = irr["nwp"][prim].factor
                        irr["blend_key"] = prim if fac == 1.0 else f"{prim} x{fac:.2f}"
                        blend_path = os.path.join(out_dir, "blend_stats.json")
                        blend_tests = load_blend_stats(blend_path, irr["blend_key"])
                        irr["blend_w"] = weights_from_stats(blend_tests, before=tlist[0],
                                                            keep=keep_hist)

        tick("Open-Meteo")
        # lead time: fixed, or automatic so the forecast reaches past the current time
        if a.lead is not None:
            lead = max(STEP_MIN, a.lead // STEP_MIN * STEP_MIN)
            latency = None
        else:
            lead, latency = auto_lead(tlist[-1], now, a.ahead)
        steps = lead // STEP_MIN

        # 2b) optical flow for all pairs in parallel (each pair is independent). Flows are
        # consumed in time order: 24-h statistics, and skill tests started at regular times.
        # Only the last n_pairs flows are kept in memory.
        n = len(grays)
        shape = grays[-1].shape
        gaps = [None] + [(tlist[k] - tlist[k - 1]).total_seconds() / 60 / STEP_MIN
                         for k in range(1, n)]
        pairs = [k for k in range(1, n) if 0 < gaps[k] <= MAX_GAP_STEPS]
        max_skill = max(1, a.skill_max // STEP_MIN)

        def pair_flow(k):
            return flow(grays[k - 1], grays[k]) / gaps[k]   # normalised to px per 10 min

        recent = collections.deque(maxlen=a.pairs)
        rows, skill_rows, pending, new_grams = [], [], collections.deque(), {}
        for k, f10 in zip(pairs, chunked_map(ex, pair_flow, pairs, chunk)):
            if gaps[k] == 1:                                  # 24-h statistics: regular pairs
                cl = np.unpackbits(clouds[k], count=shape[0] * shape[1]).astype(bool).reshape(shape)
                sel = f10[cl] if cl.sum() > 100 else f10.reshape(-1, 2)
                spd = np.hypot(sel[:, 0], sel[:, 1])
                u, w = np.median(sel[:, 0]), np.median(sel[:, 1])
                row = dict(time_utc=tlist[k].strftime("%Y-%m-%d %H:%M"),
                           cloud_fraction=round(float(cl.mean()), 4),
                           median_speed_px_per_10min=round(float(np.median(spd)), 3),
                           p90_speed_px_per_10min=round(float(np.percentile(spd, 90)), 3),
                           heading_deg=round(float((np.degrees(np.arctan2(u, -w)) + 360) % 360), 1))
                if a.km_per_px:
                    row["median_speed_kmh"] = round(float(np.median(spd)) * a.km_per_px * 6, 1)
                rows.append(row)
            recent.append((k, f10))

            # skill test from frame k (at fixed clock times, so runs can be merged)
            slot = (tlist[k].hour * 60 + tlist[k].minute) // STEP_MIN
            if k < n - 1 and slot % a.skill_every == 0 and \
                    not (cfg["can_be_dark"] and grays[k].mean() < DARK_GRAY):
                v_k, used = motion_field(recent, k, a.pairs, shape)
                if used >= max(1, a.pairs // 2):
                    pending.append(ex.submit(skill_test, band, grays, tlist, idx_by_time, k,
                                             v_k, max_skill, route, cloud_v, irr))
                    while len(pending) > chunk:
                        r_rows, r_grams = pending.popleft().result()
                        skill_rows += r_rows
                        new_grams.update(r_grams)
        for fut in pending:
            r_rows, r_grams = fut.result()
            skill_rows += r_rows
            new_grams.update(r_grams)
        del clouds
        v, _ = motion_field(recent, n - 1, a.pairs, shape)

    tick("motion+tests")
    if rows:
        with open(os.path.join(out_dir, "motion_24h.csv"), "w", newline="") as f:
            wr = csv.DictWriter(f, fieldnames=rows[0].keys())
            wr.writeheader()
            wr.writerows(rows)
        save_motion_plot(rows, band, out_dir)
        r = rows[-1]
        print(f"Now: cloud fraction {r['cloud_fraction']:.0%}, "
              f"median {r['median_speed_px_per_10min']} px/10 min, heading {r['heading_deg']}°")

    if cfg["can_be_dark"] and cv2.cvtColor(last, cv2.COLOR_BGR2HSV)[..., 2].mean() < 40:
        print(f"! Latest {band.upper()} frame is very dark (night?) – motion estimate may be "
              f"unreliable; use the B13 result.")

    if latency is not None:
        print(f"Newest image {tlist[-1]:%H:%M} UTC is {latency:.0f} min old → lead +{lead} min "
              f"(reaches {tlist[-1] + dt.timedelta(minutes=lead):%H:%M} UTC, now {now:%H:%M} UTC)")

    # blend weights: add this run's tests to the statistics, then learn from all of them
    blend_W = None
    if blend_path and new_grams:
        blend_tests.update(new_grams)
        save_blend_stats(blend_path, irr.get("blend_key", irr["nwp_primary"]), blend_tests)
        blend_W = weights_from_stats(blend_tests, keep=keep_hist)
        whole = blend_W.get("", {})
        if whole:
            shown = [L for L in (10, 30, 60, 90, 120) if L in whole]
            print(f"Blend weights now / satellite / now+change / {irr['nwp_primary']} (whole route, from "
                  f"{len({k.split('|')[0] for k in blend_tests if keep_hist(k.split('|')[0])})}"
                  f" test starts{', ' + hist_desc if hist_desc else ''}): "
                  + ", ".join(f"+{L} " + "/".join(f"{x:.2f}" for x in whole[L]) for L in shown))
            if not irr.get("blend_w"):
                print("  (the blend is scored from the next day on: its weights must come "
                      "from earlier days than the tests)")
        else:
            blend_W = None
            print(f"Blend: collecting data – weights need ≥{BLEND_MIN_STARTS} test starts "
                  f"per lead time")

    # 3) skill curve: this run, and all runs so far
    run_sum = summarize_skill(skill_rows)
    hist_rows = update_history(os.path.join(out_dir, "skill_history.csv"), skill_rows)
    n_all = sum(1 for r in hist_rows if r["band"] == band)
    hist_rows = [r for r in hist_rows if r["band"] == band and keep_hist(r["start_utc"])]
    if hist_desc:
        print(f"History filter ({hist_desc}): {len(hist_rows)} of {n_all} tests used "
              f"for history curves and blend weights")
    hist_sum = summarize_skill(hist_rows)
    n_days = len({r["start_utc"][:10] for r in hist_rows})
    if run_sum:
        with open(os.path.join(out_dir, "skill_curve.csv"), "w", newline="") as f:
            fields = list(run_sum[0].keys())
            wr = csv.DictWriter(f, fieldnames=["source"] + fields)
            wr.writeheader()
            wr.writerows([dict(source="this_run", **s) for s in run_sum])
            wr.writerows([dict(source="history", **{k: s.get(k, "") for k in fields})
                          for s in hist_sum])
        save_skill_plot(run_sum, hist_sum if len(hist_rows) > len(skill_rows) else [],
                        band, n_days, out_dir)
        starts = len({r["start_utc"] for r in skill_rows})
        show = [s for s in run_sum if s["lead_min"] in (10, 30, 60, 90, 120, lead)]
        print(f"Skill vs. persistence ({starts} start times, every {a.skill_every * STEP_MIN} min):")
        for s in show:
            extra = ""
            if s.get("skill_route") not in ("", None):
                extra = (f" | route {s['skill_route']:+.0%}, cloud/clear correct "
                         f"{s['acc_adv_route']:.0%} (persistence {s['acc_pers_route']:.0%})")
            sk = f"{s['skill']:+.0%}" if s["skill"] != "" else "n/a"
            print(f"  +{s['lead_min']:>3} min: {sk} (MAE {s['mae_adv']:.2f} vs. "
                  f"{s['mae_pers']:.2f}){extra}")
            if s.get("skill_w") not in ("", None):
                ns = [f"{lab} {s[f'skill_w{t}']:+.0%}" for t, lab in
                      (("_north", "north"), ("_south", "south"))
                      if s.get(f"skill_w{t}") not in ("", None)]
                print(f"         sunshine: MAE {s['mae_w_adv']:.0f} vs. {s['mae_w_pers']:.0f} W/m²"
                      f" (smart persistence) → skill {s['skill_w']:+.0%}"
                      + (f" ({', '.join(ns)})" if ns else ""))
                if s.get("skill_w_anc") not in ("", None):
                    print(f"         now + cloud change: MAE {s['mae_w_anc']:.0f} W/m² → skill "
                          f"{s['skill_w_anc']:+.0%}")
                mods = [f"{m} {s[f'mae_w_nwp_{m}']:.0f} W/m² ({s[f'skill_w_nwp_{m}']:+.0%})"
                        for m in NWP_MODELS if s.get(f"skill_w_nwp_{m}") not in ("", None)]
                if mods:
                    print(f"         weather models: {', '.join(mods)}")
                if s.get("skill_w_blend") not in ("", None):
                    print(f"         blend: {s['mae_w_blend']:.0f} W/m² ({s['skill_w_blend']:+.0%})")
        if len(hist_rows) > len(skill_rows):
            print(f"  history: {len({r['start_utc'] for r in hist_rows})} start times over "
                  f"{n_days} days → skill_curve.png")
        # our forecast vs. weather models (history if there is one, else this run)
        cmp_sum = hist_sum if len(hist_rows) > len(skill_rows) else run_sum
        if save_nwp_plot(cmp_sum, band, n_days if cmp_sum is hist_sum else 1,
                         os.path.join(out_dir, "nwp_compare.png")):
            for name in NWP_MODELS:
                win = [s["lead_min"] for s in cmp_sum
                       if s.get(f"mae_w_nwp_{name}") not in ("", None)
                       and s[f"mae_w_adv_vs_{name}"] < s[f"mae_w_nwp_{name}"]]
                lose = [s["lead_min"] for s in cmp_sum
                        if s.get(f"mae_w_nwp_{name}") not in ("", None)
                        and s[f"mae_w_adv_vs_{name}"] >= s[f"mae_w_nwp_{name}"]]
                if win or lose:
                    print(f"  vs. {name}: our forecast better at "
                          + (f"+{min(win)}…+{max(win)} min" if win else "no lead time")
                          + (f", {name} better at +{min(lose)}…+{max(lose)} min" if lose else "")
                          + " → nwp_compare.png")
    else:
        print("Skill: no test possible (too few images or too many gaps).")

    tick("statistics+plots")
    # 4) forecast (incremental: each step continues from the previous one)
    route_steps = max(1, a.route_lead // STEP_MIN) if route else 0
    n_fc = max(steps, route_steps)
    gif = []  # forecast frames only (up to the auto lead)
    fcs = [last]
    steps_iter = advect_iter(last, v)
    for s in range(1, n_fc + 1):
        fc = next(steps_iter)
        fcs.append(fc)
        if s > steps:
            continue
        ft = tlist[-1] + dt.timedelta(minutes=STEP_MIN * s)
        cv2.imwrite(os.path.join(out_dir, f"forecast_{ft:%Y%m%d_%H%M}.png"), B(fc))
        gif.append(label(B(fc), f"FORECAST {ft:%H:%M} UTC", (80, 80, 255)))
    final_t = tlist[-1] + dt.timedelta(minutes=STEP_MIN * steps)
    prev_t = final_t - dt.timedelta(minutes=STEP_MIN)

    cv2.imwrite(os.path.join(out_dir, "final_forecast.png"),
                label(B(fcs[steps]), f"{band.upper()} FORECAST {final_t:%Y-%m-%d %H:%M} UTC",
                      (80, 80, 255)))

    t_from = tlist[max(0, len(tlist) - 1 - a.pairs)]
    cv2.imwrite(os.path.join(out_dir, "motion_field.png"),
                label(overlay_vectors(B(fcs[steps]), v),
                      f"motion from {t_from:%H:%M}-{tlist[-1]:%H:%M} UTC on forecast {final_t:%H:%M}",
                      (0, 140, 255)))

    vis, mean_change = change_map(fcs[steps - 1], fcs[steps], v)
    name = f"change_{prev_t:%H%M}_{final_t:%H%M}.png"
    cv2.imwrite(os.path.join(out_dir, name),
                label(B(vis), f"change {prev_t:%H:%M}->{final_t:%H:%M} UTC (red=arriving, blue=leaving)",
                      (255, 255, 255)))
    print(f"Change {prev_t:%H:%M}→{final_t:%H:%M} UTC: mean |Δ| {mean_change:.2f} → {name}")

    # 5) route forecast
    if route:
        rtimes = [tlist[-1] + dt.timedelta(minutes=STEP_MIN * s) for s in range(route_steps + 1)]
        rvals = [route_values(cv2.cvtColor(fcs[s], cv2.COLOR_BGR2GRAY), route)
                 for s in range(route_steps + 1)]
        # expected sunshine: clear-sky irradiance at that place and time × calibrated index
        gclear = [clear_sky_ghi(route["lat"], route["lon"], t)[0] for t in rtimes]
        ghi = [gc * cal["k_of"](rv) for gc, rv in zip(gclear, rvals)] if cal else None
        nwp = irr.get("nwp", {}) if irr else {}
        nwp_ghi = {m: [fc.ghi_at(t, route["lat"], route["lon"]) for t in rtimes]
                   for m, fc in nwp.items()}
        # now + cloud change: newest measured value + change predicted from the moving clouds
        ghi_anc, k_now = None, None
        if ghi:
            k_now, b_now = latest_now(irr, route, grays, tlist, idx_by_time, tlist[-1])
            if k_now is not None:
                ghi_anc = [np.where(np.isfinite(k_now),
                                    gc * anchored_k(k_now, cal["k_of"], b_now, rv), g)
                           for gc, rv, g in zip(gclear, rvals, ghi)]
        prim = irr.get("nwp_primary") if irr else None
        blend = (route_blend(blend_W, irr, rtimes, gclear, k_now, ghi,
                             ghi_anc if ghi_anc else ghi, nwp_ghi[prim])
                 if blend_W and ghi and prim in nwp_ghi else None)
        fmt = lambda x: round(float(x)) if np.isfinite(x) else ""
        with open(os.path.join(out_dir, "route_forecast.csv"), "w", newline="") as f:
            wr = csv.writer(f)
            wr.writerow(["time_utc", "kind", "minutes_from_now", "name", "km", "latitude",
                         "longitude", "brightness", "cloud", "ghi_forecast_wm2",
                         "ghi_clear_sky_wm2", "ghi_reference_wm2"]
                        + [f"ghi_{m}_wm2" for m in nwp_ghi] + ["ghi_now_plus_change_wm2",
                                                                "ghi_blend_wm2"])
            for s, (t, vals) in enumerate(zip(rtimes, rvals)):
                G = irr["ref"].get(t) if irr else None
                for i in np.flatnonzero(route["inside"]):
                    wr.writerow([t.strftime("%Y-%m-%d %H:%M"), "observed" if s == 0 else "forecast",
                                 round((t - now).total_seconds() / 60), route["name"][i],
                                 round(float(route["km"][i]), 1), route["lat"][i], route["lon"][i],
                                 round(float(vals[i]), 1), int(vals[i] > cloud_v),
                                 round(float(ghi[s][i])) if ghi else "",
                                 round(float(gclear[s][i])),
                                 round(float(G[i])) if G is not None and np.isfinite(G[i]) else ""]
                                + [fmt(g[s][i]) for g in nwp_ghi.values()]
                                + [fmt(ghi_anc[s][i]) if ghi_anc else "",
                                   fmt(blend[s][i]) if blend else ""])
        save_route_plot(route, rtimes, rvals, cloud_v, now, a, band,
                        os.path.join(out_dir, "route_forecast.png"), blend or ghi_anc or ghi,
                        "Expected sunshine (blend of now / satellite / weather model)" if blend
                        else "Expected sunshine (now + cloud change)" if ghi_anc
                        else "Expected sunshine (satellite forecast)")
        car_now = car_km_at(now, a, now)
        cv2.imwrite(os.path.join(out_dir, "route_map.png"),
                    route_map(fcs[route_steps], bmask, route, georef, a.roi, rvals[-1], cloud_v,
                              f"{band.upper()} forecast {rtimes[-1]:%H:%M} UTC with route",
                              car_now))
        ok = route["inside"]
        cloudy = [f"{t:%H:%M} {100 * np.mean(v[ok] > cloud_v):.0f}%"
                  for t, v in list(zip(rtimes, rvals))[::max(1, route_steps // 3)]]
        print(f"Route cloud cover (share of route points): {', '.join(cloudy)} "
              f"→ route_map.png, route_forecast.png")
        if ghi:
            mean_ghi = [f"{t:%H:%M} {np.mean(g[ok]):.0f}" for t, g in
                        list(zip(rtimes, ghi))[::max(1, route_steps // 3)]]
            print(f"Route mean sunshine [W/m²]: {', '.join(mean_ghi)}")
            for m, g in nwp_ghi.items():
                mg = [f"{t:%H:%M} {np.nanmean(x[ok]):.0f}" if np.isfinite(x[ok]).any()
                      else f"{t:%H:%M} –" for t, x in list(zip(rtimes, g))[::max(1, route_steps // 3)]]
                print(f"  weather model {m}: {', '.join(mg)}")
            if ghi_anc:
                ma = [f"{t:%H:%M} {np.nanmean(x[ok]):.0f}" for t, x in
                      list(zip(rtimes, ghi_anc))[::max(1, route_steps // 3)]]
                print(f"  now + cloud change: {', '.join(ma)}")
            if blend:
                mb = [f"{t:%H:%M} {np.nanmean(x[ok]):.0f}" for t, x in
                      list(zip(rtimes, blend))[::max(1, route_steps // 3)]]
                print(f"  blend (best estimate): {', '.join(mb)}")

        if a.car_km is not None:
            with open(os.path.join(out_dir, "route_car.csv"), "w", newline="") as f:
                wr = csv.writer(f)
                wr.writerow(["time_utc", "minutes_from_now", "car_km", "nearest_point",
                             "brightness", "cloud", "ghi_forecast_wm2", "ghi_clear_sky_wm2",
                             "ghi_now_plus_change_wm2", "ghi_blend_wm2"])
                seq = []
                for s, (t, vals) in enumerate(zip(rtimes, rvals)):
                    if t < now:
                        continue
                    ck = car_km_at(t, a, now)
                    i = int(np.argmin(np.abs(route["km"] - ck)))
                    if not ok[i] or ck > route["km"][-1]:
                        continue
                    wr.writerow([t.strftime("%Y-%m-%d %H:%M"), round((t - now).total_seconds() / 60),
                                 round(ck, 1), route["name"][i], round(float(vals[i]), 1),
                                 int(vals[i] > cloud_v),
                                 round(float(ghi[s][i])) if ghi else "",
                                 round(float(gclear[s][i])),
                                 fmt(ghi_anc[s][i]) if ghi_anc else "",
                                 fmt(blend[s][i]) if blend else ""])
                    best = blend or ghi_anc or ghi
                    seq.append(f"{t:%H:%M} km {ck:.0f} "
                               + (f"{best[s][i]:.0f} W/m²" if best else
                                  ("CLOUD" if vals[i] > cloud_v else "clear")))
            print(f"Car (km {a.car_km:.0f} now, {a.speed:.0f} km/h): " + (", ".join(seq) or
                  "no forecast time after now") + " → route_car.csv")

    # animation of the forecast steps (Pillow -> exact frame timing in ms)
    try:
        from PIL import Image
        frames = [Image.fromarray(cv2.cvtColor(g, cv2.COLOR_BGR2RGB)) for g in gif]
        ms = int(a.gif_delay * 1000)
        durations = [ms] * (len(frames) - 1) + [ms * 2]   # hold the final forecast longer
        frames[0].save(os.path.join(out_dir, "animation.gif"), save_all=True,
                       append_images=frames[1:], duration=durations, loop=0)
    except ImportError:
        pass
    print(f"Final forecast {final_t:%H:%M} UTC → {out_dir}/final_forecast.png")
    tick("forecast+outputs")
    print("Time: " + ", ".join(f"{k} {v:.0f} s" for k, v in timing.items())
          + f" → total {sum(timing.values()):.0f} s")
    return True


# ---------------------------------------------------------------- main
def main():
    p = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("--band", choices=["both", "b13", "trm"], default="both",
                   help="which images to process (default: both)")
    p.add_argument("--data-root", default=".", help="where the himawari_<band> image folders live")
    p.add_argument("--out-root", default=".", help="where the result_<band> folders are written")
    p.add_argument("--no-download", action="store_true",
                   help="read images from existing himawari_<band> folders instead of the web")
    p.add_argument("--lead", type=int, default=None,
                   help="fixed lead from the newest image in minutes (multiple of 10); "
                        "default: automatic so the forecast reaches --ahead minutes past now")
    p.add_argument("--ahead", type=int, default=10,
                   help="how far past the current time the forecast should reach (default 10)")
    p.add_argument("--pairs", type=int, default=6,
                   help="past image pairs used for the motion field (default 6 = last hour)")
    p.add_argument("--skill-max", type=int, default=120,
                   help="longest lead time tested in the skill curve, minutes (default 120)")
    p.add_argument("--history-months", default="",
                   help="use only tests from these months for history curves and blend weights, "
                        "e.g. 8 or 7,8,9 (race in late August); default: all")
    p.add_argument("--history-days", type=int, default=0,
                   help="use only tests from the last N days for history and blend (0 = all)")
    p.add_argument("--skill-every", type=int, default=3,
                   help="start a skill test every N frames = N×10 min (default 3 = 30 min)")
    p.add_argument("--route", default=None,
                   help="route points: .py (list of dicts), .json or .csv with latitude/longitude "
                        f"(default: {DEFAULT_ROUTE} next to this script, if it exists)")
    p.add_argument("--no-route", action="store_true", help="ignore the route file")
    p.add_argument("--georef", default=None,
                   help="JSON with the map position of the image (see top of this file)")
    p.add_argument("--route-lead", type=int, default=60,
                   help="route forecast horizon from the newest image, minutes (default 60)")
    p.add_argument("--route-radius", type=int, default=2,
                   help="pixels around each route point that are averaged (default 2 → 5×5)")
    p.add_argument("--car-km", type=float, default=None,
                   help="car position now, km along the route (enables route_car.csv)")
    p.add_argument("--speed", type=float, default=80.0, help="car speed for --car-km, km/h")
    p.add_argument("--nwp-models", default=NWP_DEFAULT,
                   help="weather models from Open-Meteo to compare with, comma-separated "
                        f"(default {NWP_DEFAULT}; also icon_global; '' = none)")
    p.add_argument("--nwp-source", choices=["forecast", "previous_day1"], default="forecast",
                   help="'forecast' = latest model runs (for the past: first hours of each run); "
                        "'previous_day1' = what the models predicted 24 h earlier")
    p.add_argument("--no-bias-correction", action="store_true",
                   help="use the weather-model values as they come (default: multiply by the "
                        "factors in BIAS_FACTOR, e.g. ECMWF × 1.06)")
    p.add_argument("--nwp-cache-min", type=float, default=60,
                   help="re-use fetched weather-model data for this many minutes (default 60)")
    p.add_argument("--no-irradiance", action="store_true",
                   help="don't fetch the Open-Meteo irradiance reference (brightness only)")
    p.add_argument("--cloud-v", type=int, default=None,
                   help="min brightness for 'cloud' (default 150 for b13, 170 for trm)")
    p.add_argument("--cloud-s", type=int, default=60, help="max saturation for 'cloud'")
    p.add_argument("--km-per-px", type=float, default=None, help="adds speeds in km/h")
    p.add_argument("--roi", type=int, nargs=4, metavar=("X", "Y", "W", "H"))
    p.add_argument("--border-color", type=int, nargs=3, metavar=("R", "G", "B"),
                   help="detect borders by overlay colour instead of automatically")
    p.add_argument("--border-tol", type=int, default=40)
    p.add_argument("--static-thr", type=float, default=12, help="max 24-h std of a border pixel")
    p.add_argument("--contrast-thr", type=int, default=25, help="min contrast of a border pixel")
    p.add_argument("--draw-borders", action="store_true",
                   help="draw the coastline/state borders on the output images (off by default)")
    p.add_argument("--line-color", type=int, nargs=3, default=(255, 255, 0), metavar=("R", "G", "B"))
    p.add_argument("--line-width", type=int, default=2)
    p.add_argument("--gif-delay", type=float, default=1.5,
                   help="seconds per frame in animation.gif (default 1.5; last frame held 2x)")
    p.add_argument("--workers", type=int, default=os.cpu_count() or 4,
                   help="CPU threads for cleaning, optical flow and skill tests (default: all cores)")
    a = p.parse_args()

    route_file = a.route
    if route_file is None and not a.no_route:
        # use route_locations.py automatically: next to the script, else in the current folder
        for d in (os.path.dirname(os.path.abspath(__file__)), os.getcwd()):
            cand = os.path.join(d, DEFAULT_ROUTE)
            if os.path.isfile(cand):
                route_file = cand
                break
    if a.no_route:
        route_file = None
    route_pts = load_route(route_file) if route_file else None
    if route_pts:
        print(f"Route file: {route_file}")
    now = dt.datetime.now(dt.timezone.utc).replace(tzinfo=None)
    bands = ["b13", "trm"] if a.band == "both" else [a.band]
    done = [b for b in bands if run_band(b, a, now, route_pts)]
    if not done:
        sys.exit("No band could be processed.")
    print(f"\nDone: {', '.join(b.upper() for b in done)}")


if __name__ == "__main__":
    main()