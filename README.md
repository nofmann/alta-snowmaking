# Alta Snowmaking Outlook

Experimental dashboard showing hourly wet-bulb temperature for Alta, the main driver of snowmaking
conditions, from Open-Meteo HRRR (about 48 h) and GFS (7 days) forecasts. Wet-bulb is computed from
temperature, relative humidity and surface pressure, so it accounts for Alta's elevation.

## Run locally

```
Rscript snowmaking.R
```

Needs R with the `httr` and `jsonlite` packages. Writes `output/snowmaking.html`.

## Hosting

`.github/workflows/snowmaking.yml` runs the script every hour (17 min past) and publishes the
page with GitHub Pages. It can also be started by hand from the Actions tab.

Every run refreshes the observed wet-bulb, trend and GO/OFF signal (from the roundup scraper data). Model data
(HRRR/GFS) is only re-pulled when the cache in `data/forecast_cache.csv` is about 4 hours old (`REFRESH_H`);
the workflow commits the cache when it changes. Set `FORCE_FETCH=1` to pull regardless when running locally.

## Snowmaking rule

GO at wet-bulb 26 F or colder and falling; OFF at 25 F or warmer and rising; otherwise HOLD (WAIT while above
26). Trend is the least-squares slope of the observed wet-bulb over the last 3 hours, steady within +/-0.25 F/hr.
Thresholds are `RULE_GO`, `RULE_OFF`, `TREND_HOURS`, `TREND_FLAT` in `snowmaking.R`; "snowmaking hours" use
`SNOWMAKING_MAX` (26).

## Tuning

Settings are at the top of `snowmaking.R`: location, models, color scale, band cutoffs (`BANDS`) and
the "snowmaking hours" threshold (`SNOWMAKING_MAX`). Set `SNOW_ELEV_M` (meters) to downscale the
forecast to a specific elevation.
