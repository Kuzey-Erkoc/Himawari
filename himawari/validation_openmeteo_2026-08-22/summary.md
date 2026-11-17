# Open-Meteo sunshine forecasts vs. satellite-derived actual values

Periods (local dates): 2026-08-22 – 2026-08-29; 150 route points; models: best_match, ecmwf_ifs025, gfs_global, icon_global; daytime hours only (clear-sky GHI > 50 W/m²); every forecast is scored on the same hours. Actual = Himawari-derived GHI (Open-Meteo Satellite Radiation API), not a ground measurement.

## Hourly values
| Forecast | Mean error [W/m²] | Bias [W/m²] | Relative error | Correlation |
|---|---|---|---|---|
| ecmwf_ifs025, short lead, corrected | 38 | -7 | 7 % | 0.97 |
| best_match, short lead, corrected | 40 | -4 | 7 % | 0.96 |
| ecmwf_ifs025, 1 day ahead, corrected | 44 | -8 | 8 % | 0.95 |
| ecmwf_ifs025, 2 days ahead, corrected | 45 | -1 | 8 % | 0.95 |
| best_match, 2 days ahead, corrected | 46 | +1 | 8 % | 0.94 |
| best_match, 1 day ahead, corrected | 46 | -6 | 8 % | 0.94 |
| icon_global, short lead, corrected | 50 | -20 | 9 % | 0.96 |
| gfs_global, short lead, corrected | 52 | -11 | 9 % | 0.94 |
| icon_global, 1 day ahead, corrected | 53 | -19 | 10 % | 0.94 |
| gfs_global, 1 day ahead, corrected | 54 | -6 | 10 % | 0.93 |
| icon_global, 2 days ahead, corrected | 56 | -18 | 10 % | 0.94 |
| gfs_global, 2 days ahead, corrected | 56 | -8 | 10 % | 0.92 |
| ecmwf_ifs025, short lead | 57 | -38 | 10 % | 0.97 |
| best_match, short lead | 57 | -35 | 10 % | 0.96 |
| gfs_global, short lead | 61 | -27 | 11 % | 0.94 |
| best_match, 2 days ahead | 61 | -31 | 11 % | 0.94 |
| ecmwf_ifs025, 1 day ahead | 62 | -40 | 11 % | 0.95 |
| ecmwf_ifs025, 2 days ahead | 62 | -33 | 11 % | 0.95 |
| gfs_global, 1 day ahead | 62 | -23 | 11 % | 0.93 |
| best_match, 1 day ahead | 62 | -38 | 11 % | 0.94 |
| gfs_global, 2 days ahead | 65 | -24 | 12 % | 0.92 |
| always clear sky | 69 | +25 | 12 % | 0.87 |
| same as yesterday | 76 | -16 | 13 % | 0.84 |
| icon_global, short lead | 79 | -65 | 14 % | 0.96 |
| icon_global, 1 day ahead | 80 | -64 | 14 % | 0.94 |
| icon_global, 2 days ahead | 82 | -63 | 14 % | 0.94 |

## Daily energy (complete days only)
| Forecast | Mean error [kWh/m²] | Bias [kWh/m²] | Relative error | Correlation |
|---|---|---|---|---|
| best_match, short lead, corrected | 0.21 | -0.04 | 4 % | 0.96 |
| ecmwf_ifs025, short lead, corrected | 0.21 | -0.06 | 4 % | 0.96 |
| ecmwf_ifs025, 1 day ahead, corrected | 0.28 | -0.07 | 5 % | 0.94 |
| best_match, 1 day ahead, corrected | 0.29 | -0.07 | 5 % | 0.94 |
| best_match, 2 days ahead, corrected | 0.30 | +0.01 | 5 % | 0.93 |
| ecmwf_ifs025, 2 days ahead, corrected | 0.31 | -0.01 | 6 % | 0.92 |
| icon_global, short lead, corrected | 0.31 | -0.21 | 6 % | 0.95 |
| icon_global, 1 day ahead, corrected | 0.34 | -0.20 | 6 % | 0.93 |
| gfs_global, short lead, corrected | 0.36 | -0.11 | 6 % | 0.94 |
| gfs_global, 1 day ahead, corrected | 0.38 | -0.07 | 7 % | 0.92 |
| icon_global, 2 days ahead, corrected | 0.38 | -0.19 | 7 % | 0.91 |
| gfs_global, 2 days ahead, corrected | 0.41 | -0.08 | 7 % | 0.90 |
| best_match, short lead | 0.43 | -0.36 | 8 % | 0.96 |
| ecmwf_ifs025, short lead | 0.44 | -0.37 | 8 % | 0.96 |
| gfs_global, short lead | 0.45 | -0.27 | 8 % | 0.94 |
| gfs_global, 1 day ahead | 0.46 | -0.23 | 8 % | 0.92 |
| best_match, 2 days ahead | 0.48 | -0.31 | 9 % | 0.93 |
| gfs_global, 2 days ahead | 0.48 | -0.24 | 9 % | 0.90 |
| ecmwf_ifs025, 1 day ahead | 0.49 | -0.38 | 9 % | 0.94 |
| best_match, 1 day ahead | 0.49 | -0.38 | 9 % | 0.94 |
| ecmwf_ifs025, 2 days ahead | 0.50 | -0.32 | 9 % | 0.92 |
| always clear sky | 0.64 | +0.25 | 11 % | 0.67 |
| same as yesterday | 0.67 | -0.16 | 12 % | 0.67 |
| icon_global, 1 day ahead | 0.69 | -0.65 | 12 % | 0.93 |
| icon_global, short lead | 0.69 | -0.66 | 12 % | 0.95 |
| icon_global, 2 days ahead | 0.70 | -0.64 | 13 % | 0.91 |

## Mean error per subset – hourly [W/m²]
| Forecast | north | south | race week | change days | steady days |
|---|---|---|---|---|---|
| best_match, short lead | 57 | 58 | 57 | 75 | 49 |
| best_match, short lead, corrected | 35 | 45 | 40 | 68 | 28 |
| best_match, 1 day ahead | 58 | 67 | 62 | 83 | 53 |
| best_match, 1 day ahead, corrected | 38 | 55 | 46 | 77 | 34 |
| best_match, 2 days ahead | 57 | 65 | 61 | 80 | 52 |
| best_match, 2 days ahead, corrected | 38 | 55 | 46 | 77 | 33 |
| ecmwf_ifs025, short lead | 57 | 56 | 57 | 71 | 50 |
| ecmwf_ifs025, short lead, corrected | 34 | 42 | 38 | 62 | 28 |
| ecmwf_ifs025, 1 day ahead | 59 | 65 | 62 | 76 | 55 |
| ecmwf_ifs025, 1 day ahead, corrected | 36 | 51 | 44 | 68 | 33 |
| ecmwf_ifs025, 2 days ahead | 58 | 66 | 62 | 77 | 53 |
| ecmwf_ifs025, 2 days ahead, corrected | 37 | 53 | 45 | 71 | 32 |
| gfs_global, short lead | 53 | 68 | 61 | 94 | 49 |
| gfs_global, short lead, corrected | 43 | 61 | 52 | 90 | 38 |
| gfs_global, 1 day ahead | 56 | 69 | 62 | 101 | 49 |
| gfs_global, 1 day ahead, corrected | 46 | 62 | 54 | 98 | 38 |
| gfs_global, 2 days ahead | 57 | 72 | 65 | 106 | 50 |
| gfs_global, 2 days ahead, corrected | 47 | 66 | 56 | 103 | 39 |
| icon_global, short lead | 71 | 87 | 79 | 97 | 72 |
| icon_global, short lead, corrected | 41 | 59 | 50 | 77 | 40 |
| icon_global, 1 day ahead | 74 | 86 | 80 | 103 | 71 |
| icon_global, 1 day ahead, corrected | 47 | 60 | 53 | 87 | 40 |
| icon_global, 2 days ahead | 71 | 93 | 82 | 85 | 82 |
| icon_global, 2 days ahead, corrected | 44 | 68 | 56 | 71 | 51 |
| always clear sky | 61 | 78 | 69 | 131 | 45 |
| same as yesterday | 66 | 86 | 76 | 226 | 35 |

## Mean error per subset – daily energy [kWh/m²]
| Forecast | north | south | race week | change days | steady days |
|---|---|---|---|---|---|
| best_match, short lead | 0.48 | 0.39 | 0.43 | 0.46 | 0.41 |
| best_match, short lead, corrected | 0.20 | 0.22 | 0.21 | 0.37 | 0.14 |
| best_match, 1 day ahead | 0.48 | 0.49 | 0.49 | 0.55 | 0.45 |
| best_match, 1 day ahead, corrected | 0.23 | 0.34 | 0.29 | 0.46 | 0.21 |
| best_match, 2 days ahead | 0.47 | 0.49 | 0.48 | 0.52 | 0.44 |
| best_match, 2 days ahead, corrected | 0.23 | 0.36 | 0.30 | 0.48 | 0.20 |
| ecmwf_ifs025, short lead | 0.49 | 0.40 | 0.44 | 0.47 | 0.41 |
| ecmwf_ifs025, short lead, corrected | 0.20 | 0.22 | 0.21 | 0.37 | 0.14 |
| ecmwf_ifs025, 1 day ahead | 0.48 | 0.49 | 0.49 | 0.51 | 0.46 |
| ecmwf_ifs025, 1 day ahead, corrected | 0.23 | 0.33 | 0.28 | 0.42 | 0.21 |
| ecmwf_ifs025, 2 days ahead | 0.49 | 0.52 | 0.50 | 0.55 | 0.45 |
| ecmwf_ifs025, 2 days ahead, corrected | 0.24 | 0.37 | 0.31 | 0.50 | 0.20 |
| gfs_global, short lead | 0.39 | 0.49 | 0.45 | 0.60 | 0.39 |
| gfs_global, short lead, corrected | 0.30 | 0.42 | 0.36 | 0.56 | 0.28 |
| gfs_global, 1 day ahead | 0.41 | 0.50 | 0.46 | 0.67 | 0.39 |
| gfs_global, 1 day ahead, corrected | 0.32 | 0.43 | 0.38 | 0.62 | 0.29 |
| gfs_global, 2 days ahead | 0.43 | 0.53 | 0.48 | 0.69 | 0.40 |
| gfs_global, 2 days ahead, corrected | 0.35 | 0.46 | 0.41 | 0.65 | 0.30 |
| icon_global, short lead | 0.64 | 0.74 | 0.69 | 0.74 | 0.67 |
| icon_global, short lead, corrected | 0.22 | 0.40 | 0.31 | 0.47 | 0.24 |
| icon_global, 1 day ahead | 0.65 | 0.72 | 0.69 | 0.79 | 0.65 |
| icon_global, 1 day ahead, corrected | 0.26 | 0.41 | 0.34 | 0.56 | 0.24 |
| icon_global, 2 days ahead | 0.63 | 0.78 | 0.70 | 0.59 | 0.75 |
| icon_global, 2 days ahead, corrected | 0.25 | 0.50 | 0.38 | 0.45 | 0.36 |
| always clear sky | 0.56 | 0.71 | 0.64 | 1.25 | 0.39 |
| same as yesterday | 0.60 | 0.73 | 0.67 | 2.19 | 0.25 |

Change days: daily energy differs by more than 1 kWh/m² from the day before at that point (237 point-days; steady: 813).

## Mean error per route section – forecasts 1 day ahead [W/m², per daytime hour]
| Forecast | Darwin – Katherine | Katherine – Tennant Creek | Tennant Creek – Alice Springs | Alice Springs – Coober Pedy | Coober Pedy – Port Augusta | Port Augusta – Adelaide |
|---|---|---|---|---|---|---|
| best_match, 1 day ahead | 54 | 55 | 67 | 68 | 66 | 63 |
| best_match, 1 day ahead, corrected | 34 | 30 | 54 | 56 | 51 | 56 |
| ecmwf_ifs025, 1 day ahead | 55 | 57 | 65 | 67 | 62 | 62 |
| ecmwf_ifs025, 1 day ahead, corrected | 32 | 30 | 50 | 54 | 47 | 53 |
| gfs_global, 1 day ahead | 35 | 47 | 82 | 69 | 67 | 71 |
| gfs_global, 1 day ahead, corrected | 30 | 35 | 71 | 60 | 59 | 68 |
| icon_global, 1 day ahead | 54 | 68 | 96 | 85 | 90 | 81 |
| icon_global, 1 day ahead, corrected | 37 | 36 | 68 | 56 | 62 | 64 |
| always clear sky | 35 | 39 | 109 | 70 | 69 | 109 |
| same as yesterday | 34 | 35 | 130 | 69 | 77 | 139 |

## Mean error per route section – forecasts 1 day ahead [kWh/m², daily energy]
| Forecast | Darwin – Katherine | Katherine – Tennant Creek | Tennant Creek – Alice Springs | Alice Springs – Coober Pedy | Coober Pedy – Port Augusta | Port Augusta – Adelaide |
|---|---|---|---|---|---|---|
| best_match, 1 day ahead | 0.42 | 0.51 | 0.51 | 0.54 | 0.50 | 0.36 |
| best_match, 1 day ahead, corrected | 0.21 | 0.19 | 0.33 | 0.38 | 0.32 | 0.27 |
| ecmwf_ifs025, 1 day ahead | 0.44 | 0.50 | 0.50 | 0.55 | 0.48 | 0.36 |
| ecmwf_ifs025, 1 day ahead, corrected | 0.20 | 0.19 | 0.32 | 0.38 | 0.30 | 0.27 |
| gfs_global, 1 day ahead | 0.23 | 0.38 | 0.58 | 0.51 | 0.48 | 0.52 |
| gfs_global, 1 day ahead, corrected | 0.24 | 0.26 | 0.47 | 0.41 | 0.39 | 0.51 |
| icon_global, 1 day ahead | 0.45 | 0.65 | 0.79 | 0.74 | 0.77 | 0.63 |
| icon_global, 1 day ahead, corrected | 0.17 | 0.20 | 0.41 | 0.38 | 0.46 | 0.42 |
| always clear sky | 0.27 | 0.33 | 1.06 | 0.64 | 0.63 | 0.98 |
| same as yesterday | 0.26 | 0.30 | 1.25 | 0.52 | 0.69 | 1.24 |

## Bias correction (forecast × factor)
| Forecast | applied | factor | factor that would fit these data |
|---|---|---|---|
| best_match, short lead | fixed (BIAS_FACTOR) | 1.06 | 1.067 |
| best_match, 1 day ahead | fixed (BIAS_FACTOR) | 1.06 | 1.072 |
| best_match, 2 days ahead | fixed (BIAS_FACTOR) | 1.06 | 1.058 |
| ecmwf_ifs025, short lead | fixed (BIAS_FACTOR) | 1.06 | 1.072 |
| ecmwf_ifs025, 1 day ahead | fixed (BIAS_FACTOR) | 1.06 | 1.075 |
| ecmwf_ifs025, 2 days ahead | fixed (BIAS_FACTOR) | 1.06 | 1.063 |
| gfs_global, short lead | fixed (BIAS_FACTOR) | 1.03 | 1.05 |
| gfs_global, 1 day ahead | fixed (BIAS_FACTOR) | 1.03 | 1.042 |
| gfs_global, 2 days ahead | fixed (BIAS_FACTOR) | 1.03 | 1.044 |
| icon_global, short lead | fixed (BIAS_FACTOR) | 1.09 | 1.131 |
| icon_global, 1 day ahead | fixed (BIAS_FACTOR) | 1.09 | 1.128 |
| icon_global, 2 days ahead | fixed (BIAS_FACTOR) | 1.09 | 1.126 |

Bias > 0: forecast too sunny. 'always clear sky' and 'same as yesterday' are the references a forecast has to beat. 'corrected' = forecast multiplied by the fixed factor of main.py / trip_forecast.py (BIAS_FACTOR). These factors were derived partly from the August 2026 satellite data, so for 2026 they are not a completely independent test.