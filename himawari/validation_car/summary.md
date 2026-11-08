# Sunshine forecasts vs. our own car's weather station (WSC 2025)

Race days 24.08. – 28.08.2025; 46 daytime hours along the car's track (common sample: every method has a value). Truth = the car's own GHI sensor, hourly means. Forecast values are taken at the place where the car was in each minute.

Left out of the common sample (too few values): satellite +1 h (cloudiness stays), satellite +2 h (cloudiness stays), satellite +3 h (cloudiness stays), car sensor +1 h (measured stays), car sensor +2 h (measured stays), car sensor +3 h (measured stays)

## Error per method (hourly, W/m²)
| Method | Mean error | 90 % interval | vs. satellite (90 %) | Bias | Correlation |
|---|---|---|---|---|---|
| satellite analysis | 30 | 22–38 | +0 … +0 | +16 | 0.99 |
| best_match, short lead | 49 | 37–67 | +6 … +37 | -38 | 0.97 |
| ecmwf_ifs025, 1 day ahead | 49 | 37–66 | +6 … +37 | -31 | 0.97 |
| ecmwf_ifs025, 2 days ahead | 50 | 38–66 | +7 … +36 | -33 | 0.97 |
| ecmwf_ifs025, short lead | 51 | 38–68 | +7 … +38 | -34 | 0.97 |
| best_match, 2 days ahead | 51 | 42–62 | +9 … +34 | -45 | 0.98 |
| icon_global, 2 days ahead | 51 | 42–62 | +9 … +34 | -45 | 0.98 |
| best_match, 1 day ahead | 52 | 43–67 | +11 … +38 | -47 | 0.98 |
| icon_global, 1 day ahead | 52 | 43–67 | +11 … +38 | -47 | 0.98 |
| always clear sky | 54 | 32–80 | +3 … +49 | +22 | 0.93 |
| icon_global, short lead | 56 | 46–68 | +15 … +38 | -47 | 0.98 |
| gfs_global, short lead | 69 | 43–98 | +14 … +67 | -54 | 0.94 |
| gfs_global, 1 day ahead | 74 | 41–116 | +12 … +84 | -53 | 0.91 |
| gfs_global, 2 days ahead | 76 | 47–112 | +17 … +80 | -42 | 0.90 |

'vs. satellite': how much larger the mean error is than that of the satellite analysis. If the interval does not include 0, the difference is clear; otherwise the data cannot separate the two.

## Error vs. lead time (W/m²)

Each column on its own common sample (hours where every method has a value). Persistence needs an observation L hours earlier, so later leads lose the first hours of each day.

| Method | 0 h | +1 h | +2 h | +3 h |
|---|---|---|---|---|
| satellite (0 h: analysis, later: cloudiness stays) | 30 | 40 | 37 | 40 |
| best_match, short lead | 49 | 52 | 53 | 57 |
| ecmwf_ifs025, short lead | 51 | 54 | 57 | 61 |
| gfs_global, short lead | 69 | 73 | 77 | 83 |
| icon_global, short lead | 56 | 57 | 59 | 62 |
| best_match, 1 day ahead | 52 | 54 | 56 | 57 |
| ecmwf_ifs025, 1 day ahead | 49 | 53 | 55 | 58 |
| gfs_global, 1 day ahead | 74 | 79 | 83 | 90 |
| icon_global, 1 day ahead | 52 | 54 | 56 | 57 |
| always clear sky | 54 | 58 | 58 | 63 |
| car sensor (measured stays) | – | 47 | 68 | 75 |
| hours in sample | 46 | 42 | 37 | 32 |

## Energy while driving, per race day (kWh/m²)

| Method | 08-24 | 08-25 | 08-26 | 08-27 | 08-28 | mean abs. error |
|---|---|---|---|---|---|---|
| **measured (car)** | 6.44 | 6.48 | 4.62 | 5.11 | 2.25 | |
| satellite analysis | 6.62 | 6.61 | 4.80 | 5.16 | 2.43 | 0.14 |
| ecmwf_ifs025, 1 day ahead | 6.17 | 6.17 | 4.44 | 4.56 | 2.15 | 0.28 |
| ecmwf_ifs025, 2 days ahead | 6.19 | 6.15 | 4.42 | 4.53 | 2.08 | 0.31 |
| ecmwf_ifs025, short lead | 6.18 | 6.18 | 4.34 | 4.50 | 2.16 | 0.31 |
| best_match, short lead | 6.20 | 6.19 | 4.02 | 4.50 | 2.24 | 0.35 |
| always clear sky | 6.29 | 6.20 | 5.66 | 5.19 | 2.61 | 0.38 |
| gfs_global, 2 days ahead | 6.31 | 6.18 | 4.54 | 3.71 | 2.24 | 0.38 |
| best_match, 2 days ahead | 6.19 | 6.05 | 4.05 | 4.43 | 2.12 | 0.41 |
| icon_global, 2 days ahead | 6.19 | 6.05 | 4.05 | 4.43 | 2.12 | 0.41 |
| best_match, 1 day ahead | 6.19 | 6.06 | 3.95 | 4.55 | 2.02 | 0.43 |
| icon_global, 1 day ahead | 6.19 | 6.06 | 3.95 | 4.55 | 2.02 | 0.43 |
| icon_global, short lead | 6.17 | 6.07 | 4.18 | 4.37 | 1.96 | 0.43 |
| gfs_global, 1 day ahead | 6.32 | 6.18 | 4.37 | 3.38 | 2.23 | 0.48 |
| gfs_global, short lead | 6.37 | 6.19 | 3.63 | 4.03 | 2.21 | 0.49 |
| hours | 10 | 11 | 10 | 10 | 5 | |

## Checks

* Sensor vs. satellite in clear hours (satellite ≥ 85 % of clear sky): measured / satellite = 0.98 (median). Far from 1 means the sensor (tilt, calibration, shading) or the satellite has a systematic offset – look at the bias column above.
* Track: car speed from true − apparent wind, anchored to the control stops (route km / reconstructed km per leg: 0.73, 0.68, 0.90, 0.98, 0.95, 1.14, 0.77, 0.86, 1.00, 0.98).
* Pressure altitude vs. terrain height at the reconstructed position, per day (correlation, RMS m): 24.08. 0.95 / 27 m; 25.08. 0.99 / 18 m; 26.08. 0.92 / 51 m; 27.08. 0.81 / 32 m; 28.08. 0.80 / 25 m
* Control stops found: Darwin (start) 24.08. 08:02, Katherine 24.08. 12:23, Daly Waters 24.08. 16:16, Tennant Creek 25.08. 11:34, Barrow Creek 25.08. 14:29, Alice Springs 26.08. 09:06, Kulgera 26.08. 11:59, Coober Pedy 27.08. 09:54, Glendambo 27.08. 13:29, Port Augusta 28.08. 08:24, Adelaide (finish) 28.08. 13:00

## How to read this

* 'satellite analysis' is the best our own satellite chain can do (lead time 0). Our nowcast for +1…3 h scores about like 'satellite +L h (cloudiness stays)' (see result_*/skill_curve.csv), so those rows stand for our own forecast at that lead time.
* 'short lead' model values join the first hours of each model run and flatter the models; '1 day ahead' is what the team could see the evening before.
* Only 5 race days, mostly clear weather: the intervals are wide. Add more days (e.g. test drives) to the station_data folder to narrow them.
