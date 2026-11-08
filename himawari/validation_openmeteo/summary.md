# Open-Meteo sunshine forecasts vs. satellite-derived actual values

Periods: 2025-08-15 – 2025-09-05, 2026-08-15 – 2026-09-05; 31 route points; daytime hours only (clear-sky GHI > 50 W/m²). Actual = Himawari-derived GHI (Open-Meteo Satellite Radiation API), not a ground measurement.

## Hourly values
| Forecast | Mean error [W/m²] | Bias [W/m²] | Relative error | Correlation |
|---|---|---|---|---|
| gfs_global, short lead | 50 | -31 | 9 % | 0.96 |
| best_match, short lead | 53 | -43 | 9 % | 0.98 |
| ecmwf_ifs025, short lead | 55 | -45 | 10 % | 0.98 |
| same as yesterday | 65 | -4 | 12 % | 0.88 |
| always clear sky | 68 | +18 | 12 % | 0.90 |

## Daily energy
| Forecast | Mean error [kWh/m²] | Bias [kWh/m²] | Relative error | Correlation |
|---|---|---|---|---|
| same as yesterday | 0.41 | -0.02 | 7 % | 0.89 |
| gfs_global, short lead | 0.42 | -0.32 | 7 % | 0.96 |
| ecmwf_ifs025, short lead | 0.46 | -0.44 | 8 % | 0.98 |
| best_match, short lead | 0.46 | -0.43 | 8 % | 0.98 |
| always clear sky | 0.55 | +0.00 | 10 % | 0.89 |

## North / south / race week (hourly, mean error W/m²)
| Forecast | north | south | race week |
|---|---|---|---|
| best_match, short lead | 46 | 59 | 60 |
| ecmwf_ifs025, short lead | 51 | 60 | 61 |
| gfs_global, short lead | 37 | 63 | 60 |
| always clear sky | 49 | 88 | 69 |
| same as yesterday | 34 | 95 | 73 |

Bias > 0: forecast too sunny. 'always clear sky' and 'same as yesterday' are the references a forecast has to beat.