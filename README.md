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

`.github/workflows/snowmaking.yml` runs the script six times a day (every 4 hours, UTC) and publishes the
page with GitHub Pages. It can also be started by hand from the Actions tab.

## Tuning

Settings are at the top of `snowmaking.R`: location, models, color scale, band cutoffs (`BANDS`) and
the "snowmaking hours" threshold (`SNOWMAKING_MAX`). Set `SNOW_ELEV_M` (meters) to downscale the
forecast to a specific elevation.
