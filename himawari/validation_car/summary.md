# Sunshine forecasts vs. our own car's weather station (WSC 2025)

Race days 24.08. – 28.08.2025; 29 daytime hours along the car's track (common sample: every method has a value). Truth = the car's own GHI sensor, hourly means. Forecast values are taken at the place where the car was in each minute.

Left out of the common sample (too few values): ecmwf_ifs025, 1 day ahead, ecmwf_ifs025, 2 days ahead, satellite +1 h (cloudiness stays), satellite +2 h (cloudiness stays), satellite +3 h (cloudiness stays), car sensor +1 h (measured stays), car sensor +2 h (measured stays), car sensor +3 h (measured stays)

## Error per method (hourly, W/m²)
| Method | Mean error | 90 % interval | vs. satellite (90 %) | Bias | Correlation |
|---|---|---|---|---|---|
| satellite analysis | 27 | 16–38 | +0 … +0 | +15 | 0.99 |
| gfs_global, 1 day ahead | 31 | 20–36 | -9 … +12 | -10 | 0.99 |
| best_match, 2 days ahead | 39 | 33–50 | +2 … +27 | -33 | 0.99 |
| icon_global, 2 days ahead | 39 | 33–50 | +2 … +27 | -33 | 0.99 |
| best_match, short lead | 40 | 28–53 | -0 … +25 | -27 | 0.99 |
| gfs_global, 2 days ahead | 42 | 21–59 | -5 … +31 | -3 | 0.95 |
| ecmwf_ifs025, short lead | 43 | 31–57 | +2 … +30 | -18 | 0.98 |
| icon_global, short lead | 44 | 34–57 | +5 … +32 | -30 | 0.99 |
| gfs_global, short lead | 45 | 24–62 | -3 … +33 | -33 | 0.97 |
| best_match, 1 day ahead | 46 | 34–64 | +6 … +37 | -37 | 0.98 |
| icon_global, 1 day ahead | 46 | 34–64 | +6 … +37 | -37 | 0.98 |
| always clear sky | 53 | 24–74 | -4 … +45 | +15 | 0.92 |

'vs. satellite': how much larger the mean error is than that of the satellite analysis. If the interval does not include 0, the difference is clear; otherwise the data cannot separate the two.

## Error vs. lead time (W/m²)

Each column on its own common sample (hours where every method has a value). Persistence needs an observation L hours earlier, so later leads lose the first hours of each day.

| Method | 0 h | +1 h | +2 h | +3 h |
|---|---|---|---|---|
| satellite (0 h: analysis, later: cloudiness stays) | 19 | 21 | 19 | 16 |
| best_match, short lead | 27 | 28 | 29 | 29 |
| ecmwf_ifs025, short lead | 31 | 34 | 36 | 34 |
| gfs_global, short lead | 20 | 21 | 19 | 17 |
| icon_global, short lead | 32 | 33 | 33 | 33 |
| best_match, 1 day ahead | 30 | 31 | 31 | 30 |
| ecmwf_ifs025, 1 day ahead | 32 | 35 | 36 | 34 |
| gfs_global, 1 day ahead | 16 | 17 | 14 | 12 |
| icon_global, 1 day ahead | 30 | 31 | 31 | 30 |
| always clear sky | 26 | 28 | 27 | 27 |
| car sensor (measured stays) | – | 29 | 56 | 79 |
| hours in sample | 13 | 12 | 10 | 8 |

## Energy while driving, per race day (kWh/m²)

| Method | 08-24 | 08-25 | 08-26 | mean abs. error |
|---|---|---|---|---|
| **measured (car)** | 6.44 | 6.48 | 4.11 | |
| satellite analysis | 6.62 | 6.61 | 4.24 | 0.15 |
| gfs_global, 1 day ahead | 6.32 | 6.18 | 4.25 | 0.19 |
| ecmwf_ifs025, short lead | 6.18 | 6.18 | 4.16 | 0.20 |
| best_match, short lead | 6.20 | 6.19 | 3.87 | 0.26 |
| gfs_global, 2 days ahead | 6.31 | 6.18 | 4.47 | 0.26 |
| icon_global, short lead | 6.17 | 6.07 | 3.93 | 0.29 |
| best_match, 2 days ahead | 6.19 | 6.05 | 3.85 | 0.31 |
| gfs_global, short lead | 6.37 | 6.19 | 3.53 | 0.31 |
| icon_global, 2 days ahead | 6.19 | 6.05 | 3.85 | 0.31 |
| best_match, 1 day ahead | 6.19 | 6.06 | 3.71 | 0.36 |
| icon_global, 1 day ahead | 6.19 | 6.06 | 3.71 | 0.36 |
| always clear sky | 6.29 | 6.20 | 4.98 | 0.43 |
| hours | 10 | 11 | 8 | |

## Checks

* Sensor vs. satellite in clear hours (satellite ≥ 85 % of clear sky): measured / satellite = 0.98 (median). Far from 1 means the sensor (tilt, calibration, shading) or the satellite has a systematic offset – look at the bias column above.
* Track: car speed from true − apparent wind, anchored to the control stops (route km / reconstructed km per leg: 0.73, 0.68, 0.90, 0.98, 0.95, 1.14, 0.77, 0.86, 1.00, 0.98).
* Pressure altitude vs. terrain height at the reconstructed position, per day (correlation, RMS m): 24.08. 0.95 / 27 m; 25.08. 0.99 / 18 m; 26.08. 0.92 / 51 m; 27.08. 0.81 / 32 m; 28.08. 0.80 / 25 m
* Control stops found: Darwin (start) 24.08. 08:02, Katherine 24.08. 12:23, Daly Waters 24.08. 16:16, Tennant Creek 25.08. 11:34, Barrow Creek 25.08. 14:29, Alice Springs 26.08. 09:06, Kulgera 26.08. 11:59, Coober Pedy 27.08. 09:54, Glendambo 27.08. 13:29, Port Augusta 28.08. 08:24, Adelaide (finish) 28.08. 13:00

## How to read this

* 'satellite analysis' is the best our own satellite chain can do (lead time 0). Our nowcast for +1…3 h scores about like 'satellite +L h (cloudiness stays)' (see result_*/skill_curve.csv), so those rows stand for our own forecast at that lead time.
* 'short lead' model values join the first hours of each model run and flatter the models; '1 day ahead' is what the team could see the evening before.
* Only 5 race days, mostly clear weather: the intervals are wide. Add more days (e.g. test drives) to the station_data folder to narrow them.
