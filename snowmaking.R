#!/usr/bin/env Rscript
# Alta snowmaking dashboard (physics-only version, no ML).
# Pull GFS + HRRR from Open-Meteo -> wet-bulb temperature at altitude -> quality score -> static HTML.

suppressPackageStartupMessages({ library(httr); library(jsonlite) })

# ---- config ---------------------------------------------------------------
LAT <- 40.58083; LON <- -111.63722
TZ  <- "America/Denver"
MODELS <- c(hrrr = "gfs_hrrr", gfs = "gfs_global")   # HRRR is preferred where it exists (~48 h); gfs_seamless is HRRR-blended, so not used
ELEV_M <- Sys.getenv("SNOW_ELEV_M", "")                 # optional: downscale to a specific elevation (m)

# Color scale and band names follow the Snow State wet-bulb chart (snow-state.com/snowmaking-charts).
# Colors are one per whole degree F, sampled from that chart; 11 and below share one color, 29 and above another.
# Band cutoffs (upper wet-bulb value, whole deg F) are inferred from where the chart places its labels; edit to suit.
SCALE <- c("11" = "#E463FF", "12" = "#D560FF", "13" = "#D35FFF", "14" = "#B658FF", "15" = "#B157FF",
           "16" = "#9750FF", "17" = "#8F4EFF", "18" = "#7849FF", "19" = "#6845FF", "20" = "#5941FF",
           "21" = "#493DFF", "22" = "#3939FF", "23" = "#525AFF", "24" = "#6B7BFF", "25" = "#859DFF",
           "26" = "#9EBEFF", "27" = "#B8E0FF", "28" = "#FFF2CC", "29" = "#FFB1B4")
BANDS <- c(great = 17, good = 22, marginal = 27, borderline = 28)   # anything warmer is "Too warm"
SNOWMAKING_MAX <- 27                                                 # "snowmaking hours" = wet-bulb at or below this
SCORE_X <- c(20, 24, 26, 28, 30); SCORE_Y <- c(100, 75, 50, 25, 0)

script_dir <- local({
  a <- grep("^--file=", commandArgs(FALSE), value = TRUE)
  if (length(a)) dirname(normalizePath(sub("^--file=", "", a[1]))) else getwd()
})
OUT_DIR <- file.path(script_dir, "output"); dir.create(OUT_DIR, showWarnings = FALSE, recursive = TRUE)

# ---- fetch ----------------------------------------------------------------
fetch_forecast <- function() {
  q <- list(latitude = LAT, longitude = LON,
            hourly = "temperature_2m,relative_humidity_2m,surface_pressure,wind_speed_10m",
            models = paste(MODELS, collapse = ","), timezone = TZ,
            temperature_unit = "fahrenheit", wind_speed_unit = "mph", forecast_days = 7)
  if (nzchar(ELEV_M)) q$elevation <- ELEV_M
  r <- GET("https://api.open-meteo.com/v1/forecast", query = q,
           user_agent("alta-snowmaking (R)"), timeout(30))
  stop_for_status(r)
  j <- fromJSON(content(r, "text", encoding = "UTF-8"))
  h <- j$hourly
  rows <- lapply(names(MODELS), function(nm) {
    m <- MODELS[[nm]]; col <- function(v) h[[paste0(v, "_", m)]]
    if (is.null(col("temperature_2m"))) return(NULL)
    data.frame(time = h$time, model = nm, temp_f = col("temperature_2m"),
               rh = col("relative_humidity_2m"), p_hpa = col("surface_pressure"),
               wind_mph = col("wind_speed_10m"), stringsAsFactors = FALSE)
  })
  list(data = do.call(rbind, rows), elev_m = j$elevation)
}

# ---- physics --------------------------------------------------------------
# Wet-bulb from T, RH and station pressure (psychrometric equation, solved numerically).
# Pressure matters: at ~9,000 ft the wet-bulb is noticeably lower than a sea-level formula gives.
wetbulb_f <- function(t_f, rh, p_hpa) {
  mapply(function(tf, rh, p) {
    if (anyNA(c(tf, rh, p))) return(NA_real_)
    t <- (tf - 32) * 5 / 9; rh <- min(max(rh, 1), 100)
    es <- function(x) 6.112 * exp(17.62 * x / (243.12 + x))
    e <- rh / 100 * es(t)
    f <- function(tw) es(tw) - 0.00066 * (1 + 0.00115 * tw) * p * (t - tw) - e
    tw <- uniroot(f, c(t - 40, t), tol = 1e-4)$root
    tw * 9 / 5 + 32
  }, t_f, rh, p_hpa)
}

quality <- function(wb) {
  v <- round(wb)
  ifelse(is.na(v), NA_character_,
  ifelse(v <= BANDS[["great"]], "Great",
  ifelse(v <= BANDS[["good"]], "Good",
  ifelse(v <= BANDS[["marginal"]], "Marginal",
  ifelse(v <= BANDS[["borderline"]], "Borderline", "Too warm")))))
}

scale_color <- function(wb) unname(SCALE[as.character(pmin(pmax(round(wb), 11), 29))])
# white or near-black text, whichever contrasts better with the cell color
text_color <- function(hex) {
  lin <- function(u) ifelse(u <= 0.03928, u / 12.92, ((u + 0.055) / 1.055)^2.4)
  lum <- function(h) { v <- col2rgb(h)[, 1] / 255; sum(c(0.2126, 0.7152, 0.0722) * lin(v)) }
  L <- lum(hex)
  if ((1.05) / (L + 0.05) >= (L + 0.05) / (lum("#1b2630") + 0.05)) "#ffffff" else "#1b2630"
}

# ---- build ----------------------------------------------------------------
fc <- fetch_forecast(); d <- fc$data
d$wb    <- round(wetbulb_f(d$temp_f, d$rh, d$p_hpa), 1)
d$score <- round(approx(SCORE_X, SCORE_Y, xout = d$wb, rule = 2)$y)
d$class <- quality(d$wb)

# best available: HRRR where it has data, else GFS
hr <- d[d$model == "hrrr" & !is.na(d$wb), ]
gf <- d[d$model == "gfs", ]
best <- rbind(hr, gf[!gf$time %in% hr$time, ])
best <- best[order(best$time), ]
best$src <- toupper(best$model)
best$date <- substr(best$time, 1, 10); best$hour <- as.integer(substr(best$time, 12, 13))
best$i <- seq_len(nrow(best)) - 1L   # index into the hover-popup data

now_str <- format(Sys.time(), "%Y-%m-%dT%H:00", tz = TZ)
fut <- best[best$time >= now_str, ]
good_plus <- function(x) sum(round(x$wb) <= SNOWMAKING_MAX, na.rm = TRUE)   # hours at Marginal or better
now_row <- fut[1, ]
w24 <- fut[seq_len(min(24, nrow(fut))), ]; w72 <- fut[seq_len(min(72, nrow(fut))), ]
bst <- w72[which.min(w72$wb), ]

message(sprintf("Elevation used by Open-Meteo: %.0f m (%.0f ft)", fc$elev_m, fc$elev_m * 3.281))
message(sprintf("Rows: %d (hrrr %d, gfs %d). Now: wb %.1f F (%s). Good+ hours next 24h: %d, 72h: %d",
                nrow(d), sum(d$model == "hrrr"), sum(d$model == "gfs"), now_row$wb, now_row$class,
                good_plus(w24), good_plus(w72)))

# ---- html -----------------------------------------------------------------
day_label <- function(ds) format(as.Date(ds), "%a %b %d")

strip_rows <- vapply(unique(best$date), function(dt) {
  x <- best[best$date == dt, ]
  cells <- vapply(0:23, function(hh) {
    r <- x[x$hour == hh, ]
    if (!nrow(r)) return("<td class='c-na'></td>")
    past <- if (r$time < now_str) " past" else ""
    bg <- scale_color(r$wb)
    sprintf("<td class='cell%s' data-i='%d' style='background:%s;color:%s'>%s</td>", past, r$i, bg, text_color(bg), round(r$wb))
  }, "")
  ok <- x[x$time >= now_str, ]
  summ <- if (nrow(ok)) sprintf("%d h", good_plus(ok)) else "&ndash;"
  sprintf("<tr><th>%s</th>%s<td class='sum'>%s</td></tr>", day_label(dt), paste(cells, collapse = ""), summ)
}, "")

gi <- match(paste(best$time, "gfs"), paste(d$time, d$model))
is_hrrr <- best$src == "HRRR"
chart <- list(
  labels = format(as.POSIXct(paste0(best$date, " ", sprintf("%02d", best$hour), ":00"), tz = "UTC"), "%a %H:00"),
  hrrr = ifelse(is_hrrr, best$wb, NA), hrrr_t = ifelse(is_hrrr, best$temp_f, NA),
  hrrr_rh = ifelse(is_hrrr, best$rh, NA), hrrr_w = ifelse(is_hrrr, best$wind_mph, NA),
  gfs = d$wb[gi], gfs_t = d$temp_f[gi], gfs_rh = d$rh[gi], gfs_w = d$wind_mph[gi]
)

# per-cell details for the hover popup
info <- list(l = sprintf("%s %02d:00", day_label(best$date), best$hour), s = best$src, c = best$class,
             wb = best$wb, t = best$temp_f, rh = best$rh, w = best$wind_mph)

# Reference chart (temp x RH -> wet-bulb), same layout as the Snow State chart but computed at Alta's pressure.
RC_T <- 14:38; RC_RH <- seq(20, 100, 5)
p_ref <- round(median(best$p_hpa, na.rm = TRUE))
g <- expand.grid(t = RC_T, rh = RC_RH)
g$wb <- round(wetbulb_f(g$t, g$rh, p_ref))
rc_rows <- vapply(RC_T, function(tt) {
  cells <- vapply(RC_RH, function(rr) {
    w <- g$wb[g$t == tt & g$rh == rr]; bg <- scale_color(w)
    sprintf("<td data-t='%d' data-rh='%d' style='background:%s;color:%s'>%d</td>", tt, rr, bg, text_color(bg), as.integer(w))
  }, "")
  sprintf("<tr><th>%d</th>%s</tr>", tt, paste(cells, collapse = ""))
}, "")
rc_head <- paste0("<tr><th>&deg;F</th>", paste0("<th>", RC_RH, "%</th>", collapse = ""), "</tr>")

tile <- function(label, value, sub) sprintf("<div class='tile'><div class='tl'>%s</div><div class='tv'>%s</div><div class='ts'>%s</div></div>", label, value, sub)
tiles <- paste0(
  tile("Wet-bulb now", sprintf("%.0f&deg;F", now_row$wb), sprintf("%s &middot; temp %.0f&deg;F, RH %.0f%%", now_row$class, now_row$temp_f, now_row$rh)),
  tile("Next 24 hours", sprintf("%d of %d", good_plus(w24), nrow(w24)), "hours with wet-bulb 27&deg;F or colder"),
  tile("Next 72 hours", sprintf("%d of %d", good_plus(w72), nrow(w72)), "hours with wet-bulb 27&deg;F or colder"),
  tile("Best hour, next 72 h", sprintf("%.0f&deg;F", bst$wb), sprintf("%s %02d:00 &middot; %s", day_label(bst$date), bst$hour, bst$src))
)

html <- r"---(<!doctype html><html lang="en"><head><meta charset="utf-8">
<meta name="viewport" content="width=device-width,initial-scale=1">
<title>Alta Snowmaking Outlook</title>
<style>
:root{--bg:#f6f8fa;--card:#fff;--ink:#1b2630;--dim:#5b6b79;--line:#d9e0e6;
 --ex:#0b3c8a;--go:#3b82c4;--ma:#bfd7ec;--no:#eceff2}
@media (prefers-color-scheme:dark){:root{--bg:#10161c;--card:#18212a;--ink:#e6edf3;--dim:#93a3b1;--line:#2c3945;
 --ex:#0b3c8a;--go:#3b82c4;--ma:#2b4a66;--no:#222d37}}
*{box-sizing:border-box}body{margin:0;background:var(--bg);color:var(--ink);font:14px/1.45 system-ui,Segoe UI,sans-serif}
main{max-width:1100px;margin:0 auto;padding:20px 16px 40px}
h1{margin:0 0 2px;font-size:22px}.sub{color:var(--dim);margin-bottom:16px}
.tiles{display:grid;grid-template-columns:repeat(auto-fit,minmax(200px,1fr));gap:10px;margin-bottom:18px}
.tile{background:var(--card);border:1px solid var(--line);border-radius:8px;padding:12px 14px}
.tl{color:var(--dim);font-size:12px}.tv{font-size:28px;font-weight:600}.ts{color:var(--dim);font-size:12px}
.card{background:var(--card);border:1px solid var(--line);border-radius:8px;padding:14px;margin-bottom:18px}
h2{font-size:15px;margin:0 0 10px}
.scroll{overflow-x:auto}
table{border-collapse:separate;border-spacing:2px;font-size:11px}
th{text-align:left;font-weight:600;white-space:nowrap;padding-right:8px}
td{min-width:28px;height:26px;text-align:center;border-radius:3px}
td.hd{background:none;color:var(--dim);height:auto}
.c-na{background:none}
.past{opacity:.45}.sum{min-width:48px;font-weight:600;background:none}
.lg{display:grid;grid-template-columns:repeat(14,1fr) 1.8fr 1.8fr;margin-top:12px;max-width:760px;font-size:11px;text-align:center}
.lg div{padding:5px 0;font-weight:600}
.lg .lab{font-weight:600;font-size:12px;padding:6px 2px 0;color:var(--ink)}
.note{color:var(--dim);font-size:12px;margin-top:8px}
canvas{max-height:320px}
td.cell{cursor:default}
#pop{position:fixed;display:none;z-index:10;pointer-events:none;width:min(640px,calc(100vw - 16px));
 background:var(--card);border:1px solid var(--line);border-radius:8px;box-shadow:0 8px 28px rgba(0,0,0,.3);padding:10px 12px}
#popHd{font-size:13px;line-height:1.5;margin-bottom:6px}
table.rc{width:100%;table-layout:fixed;border-spacing:1px;font-size:10.5px}
.rc th{text-align:center;padding:0;font-weight:600;color:var(--dim)}
.rc td{min-width:0;height:17px;border-radius:2px}
.rc td.hl{outline:3px solid var(--ink);outline-offset:-1px;box-shadow:0 0 0 2px var(--card);position:relative;z-index:1;font-weight:700}
#pop .lg{max-width:none;margin-top:8px}
</style></head><body><main>
<h1>Alta Snowmaking Outlook</h1>
<div class="sub">Generated __GEN__ &middot; Open-Meteo HRRR (to ~48 h) then GFS &middot; model elevation __ELEV__</div>
<div class="tiles">__TILES__</div>
<div class="card"><h2>Hourly wet-bulb temperature (&deg;F) &mdash; hour of day across, hours colored by snowmaking quality</h2>
<div class="scroll"><table id="strip"><tr><td class="hd"></td>__HDR__<td class="hd">Hours<br>&le;27&deg;F</td></tr>__ROWS__</table></div>
__LEGEND__
<div class="note">Hover (or tap) an hour to see its conditions and where it falls on the temperature / humidity chart. The last column counts upcoming hours at 27&deg;F wet-bulb or colder. Faded cells are already past. Colors and band names follow the Snow State wet-bulb chart (snow-state.com).</div></div>
<div class="card"><h2>Model comparison</h2><canvas id="ch"></canvas>
<div class="note">Wet-bulb computed from temperature, relative humidity and surface pressure. HRRR covers about 48 h; GFS covers the full period. Hover for temperature, humidity and wind (shown, not scored).</div></div>
<div id="pop"><div id="popHd"></div><table class="rc">__RCHEAD____RCROWS__</table>__LEGEND__
<div class="note">Temperature / humidity chart at Alta's pressure (__PREF__ hPa). Snow State's chart is at sea level, so its values run 1&ndash;2&deg;F warmer. The outlined cell is the hovered hour, with temp and RH rounded to the chart's grid.</div></div>
<script src="https://cdnjs.cloudflare.com/ajax/libs/Chart.js/4.4.1/chart.umd.min.js"></script>
<script>
const D=__DATA__, I=__INFO__;
const pop=document.getElementById('pop'), hd=document.getElementById('popHd');
let hl=null;
const BAND={'Great':'Great Snowmaking','Good':'Good Snowmaking','Marginal':'Marginal Snowmaking','Borderline':'Borderline','Too warm':'Too Warm'};
function hide(){pop.style.display='none';}
function show(td){
  const i=+td.dataset.i, t=I.t[i], rh=I.rh[i];
  const tr=Math.min(38,Math.max(14,Math.round(t))), rr=Math.min(100,Math.max(20,Math.round(rh/5)*5));
  if(hl) hl.classList.remove('hl');
  hl=pop.querySelector('.rc td[data-t="'+tr+'"][data-rh="'+rr+'"]'); if(hl) hl.classList.add('hl');
  const off=(t<13.5||t>38.5||rh<18)?' <span style="color:var(--dim)">(outside the chart range; nearest cell outlined)</span>':'';
  hd.innerHTML='<b>'+I.l[i]+'</b> &middot; '+I.s[i]+' &middot; '+BAND[I.c[i]]+'<br>Wet-bulb <b>'+I.wb[i].toFixed(1)+'&deg;F</b> &middot; temp '+Math.round(t)+'&deg;F &middot; RH '+Math.round(rh)+'% &middot; wind '+Math.round(I.w[i])+' mph'+off;
  pop.style.display='block';
  const r=td.getBoundingClientRect(), w=pop.offsetWidth, h=pop.offsetHeight;
  pop.style.left=Math.max(8,Math.min(r.left+r.width/2-w/2,innerWidth-w-8))+'px';
  let y=r.bottom+8; if(y+h>innerHeight-8) y=Math.max(8,r.top-h-8);
  pop.style.top=y+'px';
}
const strip=document.getElementById('strip');
strip.addEventListener('mouseover',e=>{const td=e.target.closest('td.cell'); if(td) show(td);});
strip.addEventListener('mouseleave',hide);
document.addEventListener('click',e=>{const td=e.target.closest('td.cell'); if(td) show(td); else hide();});
addEventListener('scroll',hide,{passive:true});
const dark=matchMedia('(prefers-color-scheme: dark)').matches, tc=dark?'#93a3b1':'#5b6b79', gc=dark?'#2c3945':'#e3e8ed';
const ref=v=>({data:D.labels.map(()=>v),borderColor:'#8a97a3',borderDash:[5,4],borderWidth:1,pointRadius:0,label:v+'°F',ref:true});
const tipLabel=c=>{const k=c.datasetIndex===0?'hrrr':'gfs',i=c.dataIndex;
 return [c.dataset.label+': '+c.parsed.y.toFixed(1)+'°F wet-bulb',
  '   temp '+Math.round(D[k+'_t'][i])+'°F · RH '+Math.round(D[k+'_rh'][i])+'% · wind '+Math.round(D[k+'_w'][i])+' mph'];};
new Chart(document.getElementById('ch'),{type:'line',data:{labels:D.labels,datasets:[
 {label:'HRRR',data:D.hrrr,borderColor:'#0b3c8a',backgroundColor:'#0b3c8a',borderWidth:2,pointRadius:0,spanGaps:false},
 {label:'GFS',data:D.gfs,borderColor:'#e07a1f',backgroundColor:'#e07a1f',borderWidth:2,pointRadius:0},
 ref(28),ref(22)]},
 options:{animation:false,interaction:{mode:'index',intersect:false},
  scales:{x:{ticks:{color:tc,maxTicksLimit:14,maxRotation:0},grid:{color:gc}},y:{title:{display:true,text:'Wet-bulb (°F)',color:tc},ticks:{color:tc},grid:{color:gc}}},
  plugins:{legend:{labels:{color:tc,filter:i=>!/°F$/.test(i.text)}},
   tooltip:{filter:i=>!i.dataset.ref,callbacks:{label:tipLabel}}}}});
</script></main></body></html>)---"

rep <- function(s, k, v) gsub(k, gsub("\\", "\\\\", v, fixed = TRUE), s, fixed = TRUE)
html <- rep(html, "__GEN__", format(Sys.time(), "%b %d %H:%M %Z", tz = TZ))
html <- rep(html, "__ELEV__", sprintf("%.0f ft", fc$elev_m * 3.281))
leg_vals <- c(11, 12, 14, 16, 18:27, 28, 29)
leg_txt  <- c("&le;11", leg_vals[2:14], "28", "29+")
leg_cells <- vapply(seq_along(leg_vals), function(i) {
  bg <- scale_color(leg_vals[i]); sprintf("<div style='background:%s;color:%s'>%s</div>", bg, text_color(bg), leg_txt[i])
}, "")
lab <- function(txt, from, to) sprintf("<div class='lab' style='grid-column:%d / %d'>%s</div>", from, to + 1, txt)
legend <- paste0("<div class='lg'>", paste(leg_cells, collapse = ""),
                 lab("Great Snowmaking", 1, 4), lab("Good Snowmaking", 5, 9), lab("Marginal Snowmaking", 10, 14),
                 lab("Borderline", 15, 15), lab("Too Warm", 16, 16), "</div>")
html <- rep(html, "__LEGEND__", legend)
html <- rep(html, "__TILES__", tiles)
html <- rep(html, "__HDR__", paste0("<td class='hd'>", 0:23, "</td>", collapse = ""))
html <- rep(html, "__ROWS__", paste(strip_rows, collapse = "\n"))
html <- rep(html, "__RCHEAD__", rc_head)
html <- rep(html, "__RCROWS__", paste(rc_rows, collapse = "\n"))
html <- rep(html, "__PREF__", as.character(p_ref))
html <- rep(html, "__INFO__", toJSON(info, auto_unbox = TRUE, na = "null", digits = NA))
html <- rep(html, "__DATA__", toJSON(chart, auto_unbox = TRUE, na = "null"))
out <- file.path(OUT_DIR, "snowmaking.html")
writeLines(html, out, useBytes = TRUE)
message("Report: ", out)
