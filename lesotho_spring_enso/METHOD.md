## What the tool does

For each spring it (1) delineates the topographic catchment upstream of the spring, (2) extracts monthly rainfall and NDVI averaged over that catchment, (3) labels every rainy season as El Niño, La Niña or neutral, (4) compares rainfall, vegetation and — where records exist — spring flow between those phases, and (5) combines the results into a 0–1 El Niño susceptibility score for prioritising communities.

### 1 · Catchments
* **DEM:** Copernicus GLO-30 (30 m) from Microsoft Planetary Computer, or any local DEM. Reprojected to UTM 35S (EPSG:32735).
* **Delineation:** pysheds D8 — pits and depressions filled, flats resolved, flow direction and accumulation computed in a window around the spring. If the catchment touches the window edge the window is doubled (up to the limit in the settings; beyond that the catchment is flagged *truncated*).
* **Snapping:** the spring is moved to the highest-accumulation cell within the snap radius, because springs sit on or beside drainage lines and GPS and DEM errors are tens of metres. `snap_dist_m` and `on_drainage_line` record what happened.
* **Caveat:** a topographic catchment is a *proxy* for the recharge area. In Lesotho's basalts and Clarens sandstone, fractures and dolerite dykes can carry groundwater across surface divides. Very small catchments (< 0.05 km²) mean the spring sits high on a slope, and the recharge area is poorly defined.

### 2 · Rainfall and NDVI
* **Rainfall:** CHIRPS monthly, 0.05° (~5 km), 1981–present. v3 is the default; CHC stops producing v2 after December 2026. Values are area-weighted by the fraction of each CHIRPS cell inside the catchment, so small catchments take the value of the cell(s) they sit in. SPI-1, 3, 6, 12 and 24 are computed by fitting a gamma distribution per calendar month.
* **NDVI:** MODIS MOD13Q1 (Terra; Aqua optional) v6.1, 250 m, 16-day composites from 2000, read from Planetary Computer. Only pixels flagged good or marginal are kept; snow and cloud are dropped. A composite is used if at least half of the catchment is good quality. Composites are averaged into calendar months and converted to z-scores per calendar month. MODIS is near the end of its life; VIIRS (VNP13) is the successor, and a local folder of monthly NDVI GeoTIFFs can be used instead (for example Landsat or Sentinel-2 NDVI prepared in QGIS).

### 3 · ENSO years
* **Index:** NOAA CPC Oceanic Niño Index (ONI). A month is El Niño (La Niña) when it falls in a run of at least 5 overlapping 3-month seasons ≥ +0.5 °C (≤ −0.5 °C).
* **Hydrological years** run October–September and are labelled like *2015/16*. A year takes the ENSO phase of its DJF season. In southern Africa, El Niño summers are usually drier and La Niña summers wetter, but the link is not guaranteed in every year or every valley. That is why the tool measures it for each catchment instead of assuming it.

### 4 · Per-year and composite statistics
| Measure | Definition |
|---|---|
| Rainy-season rain anomaly | Oct–Mar total as % departure from the mean of all complete seasons |
| Peak NDVI z | mean NDVI z-score over Jan–Apr |
| Flow impact z | mean flow anomaly in calendar year *Y*, the months after the rainy season that ended in *Y* |
| Composite | mean of each measure over El Niño, neutral and La Niña years; Mann–Whitney p-value for El Niño vs other years |

**Flow anomaly:** flow is log-transformed, and status is scored as dry 0, reduced 0.5, flowing 1. It is then standardised against the same calendar month when the record has at least 3 observations in at least 9 months. With fewer records it is standardised against the same half-year (wet Oct–Mar or dry Apr–Sep), and otherwise against the whole record.

**Rainfall memory:** the spring's flow anomaly is correlated (Spearman) with SPI-1 to SPI-24 at lags of 0–6 months. The best pair shows how quickly the spring follows rainfall. A spring that follows SPI-3 dries up within one failed season. A spring that follows SPI-24 is buffered by storage, but recovers slowly after a multi-year drought.

### 5 · Susceptibility score (0–1)
| Component | 0 | 1 | Default weight |
|---|---|---|---|
| **R** rain exposure: El Niño rainy-season anomaly | ≥ 0 % | ≤ −30 % | 0.35 |
| **V** vegetation: El Niño peak NDVI z | ≥ 0 | ≤ −1.5 | 0.15 |
| **H** spring response, from flow records | no El Niño drop, no rainfall link | large El Niño drop, follows short rainfall windows | 0.50 |
| **H** without flow records (proxy) | 0.25 at catchment ≥ 50 km² | 0.75 at catchment ≤ 0.1 km² | 0.50 |

The proxy is held between 0.25 and 0.75 on purpose, so a spring with no records cannot outrank one whose records show it failing in El Niño years.

The score is the weighted mean of the available components. **High** is ≥ 0.6, **Medium** is 0.4–0.6 and **Low** is < 0.4. **Confidence** is:
* *High* — flow records cover at least 2 El Niño years and at least 24 months.
* *Medium* — some flow records.
* *Low* — the catchment-size proxy only.

R and V change slowly across space, because CHIRPS cells are about 5 km across. Differences between neighbouring springs therefore come mostly from H, so flow records matter most. Where confidence is Low, treat the class as a screening result to verify in the field.

### Using the result
1. Sort by **priority rank**, then look at the confidence column.
2. For High/Low-confidence springs, schedule a flow check at the end of the next dry season (Aug–Oct).
3. When NOAA or SAWS issues an El Niño watch for the coming summer, the High-susceptibility springs, and the communities that depend on them, are the first candidates for pre-positioned water supply, tanks or spring protection work.
