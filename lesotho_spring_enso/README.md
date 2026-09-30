# Spring El Niño susceptibility tool — Lesotho

A Python application that, for any set of springs:

1. **delineates the upstream catchment** of each spring from a 30 m DEM,
2. extracts **monthly CHIRPS rainfall** (1981–present) and **MODIS NDVI** (2000–present) for each catchment,
3. labels every rainy season **El Niño / neutral / La Niña** from NOAA's Oceanic Niño Index,
4. compares rainfall, vegetation and your **spring flow records** between those phases, and finds how quickly each spring follows rainfall,
5. produces a **susceptibility map** (High / Medium / Low, with a confidence level) to prioritise communities for pre-positioned assistance.

The method is described in `METHOD.md`, which also appears in the app's *About the method* tab.

---

## Install

Conda is the easiest way on Windows, because it brings GDAL with it:

```bash
conda env create -f environment.yml
conda activate springs-enso
```

or with pip (Linux/macOS, or Windows with Python 3.10–3.12):

```bash
pip install -r requirements.txt
```

## Run the app

```bash
streamlit run app.py
```

A browser tab opens with the app.
1. **Springs:** click on the map to add springs, or upload a GeoPackage, GeoJSON, zipped shapefile, or a CSV/Excel file with `lon`,`lat` columns. Choose the ID column; you can edit IDs and names in the table.
2. **Flow records (optional):** upload a CSV/Excel file, one row per visit. Use the *Download template* button for the format. Columns are matched loosely:

   | spring_id | date | flow_lps | status |
   |---|---|---|---|
   | SP001 | 2016-02-15 | 0.42 | flowing |
   | SP002 | 2016-08-21 |  | dry |

   Measured flow, occasional readings and status-only records (`flowing` / `reduced` / `dry`) all work. The tool decides per spring which one to use.
3. Press **Run analysis**, then use the *Susceptibility map* and *Time series & ENSO* tabs. Download the GeoPackage (springs, catchments, snapped outlets), the Excel workbook (all tables) or the monthly CSV.

## Batch mode

```bash
python run_batch.py --points springs.gpkg --flows spring_flows.csv --out results
```

This writes `spring_enso_susceptibility.gpkg`, `spring_enso_results.xlsx`, three CSV files and `susceptibility_map.html`. Run `python run_batch.py -h` for all options (local data folders, snap radius, skip NDVI, and so on).

## Data sources and download volume

| Data | Source | Access | First-run volume (approx.) |
|---|---|---|---|
| DEM | Copernicus GLO-30 | Microsoft Planetary Computer, no account | a few MB per spring |
| Rainfall | CHIRPS v3 monthly, 0.05° | data.chc.ucsb.edu (only the Lesotho window is read) | small; ~540 months |
| NDVI | MODIS MOD13Q1 v6.1, 250 m, 16-day | Planetary Computer, no account | a few hundred MB for a cluster of nearby springs, up to a few GB for springs spread across the country |
| ENSO | NOAA CPC ONI (`oni.ascii.txt`) | cpc.ncep.noaa.gov | tiny |

Everything is **cached** in `~/.lesotho_spring_enso_cache`, set in the sidebar. Re-runs, and new springs near old ones, reuse what has already been downloaded.

**Slow or no internet:** every source can be replaced by local files, set in the sidebar or with command-line flags:
- `--dem` — any DEM GeoTIFF (for example a national 30 m DEM),
- `--chirps-dir` — a folder of monthly CHIRPS GeoTIFFs (names containing `YYYY.MM`),
- `--ndvi local --ndvi-dir` — monthly NDVI GeoTIFFs (names containing `YYYY_MM` or similar; integer files are treated as scaled by 10 000). Landsat or Sentinel-2 NDVI prepared in QGIS, as in the AWEP training guides, can go here.
- `--oni` — a saved copy of `oni.ascii.txt`.

Notes:
- CHC stops producing **CHIRPS v2** after December 2026, so v3 is the default.
- **MODIS** Terra and Aqua are near the end of their missions. When composites stop appearing, switch to a local NDVI folder (for example VIIRS VNP13 or Landsat/Sentinel-2).

## Check the data connections

```bash
python check_connection.py
```

This tests MODIS NDVI, CHIRPS, the Copernicus DEM and the NOAA ONI file one at a time, in under a minute, and prints the error for any that fails. If the app shows an NDVI warning, run this first.

## Test without internet

```bash
python -m tests.test_pipeline
```

This builds a synthetic DEM, rainfall, NDVI, ONI and flow dataset (clearly not real data) and runs the whole pipeline. It checks that a spring built to respond fast to El Niño is ranked above a buffered one.

## Project layout

```
app.py                 Streamlit app
run_batch.py           command-line batch runner
METHOD.md              method notes (also shown in the app)
springtool/
  config.py            settings
  dem.py               DEM download / reprojection / cache
  watershed.py         pysheds delineation, snapping, window expansion
  chirps.py            CHIRPS monthly rainfall
  ndvi.py              MODIS (or local) NDVI
  enso.py              ONI, El Niño / La Niña phases, hydrological years
  flows.py             flow-record reader (flow or status)
  analysis.py          SPI, anomalies, ENSO composites, flow response, scoring
  pipeline.py          end-to-end run and exports
  maps.py              folium map and plotly charts
tests/                 synthetic data generator and end-to-end test
```

## Important caveats

- A topographic catchment is only a **proxy for a spring's recharge area**. Fractured basalt and dolerite dykes can move groundwater across surface divides. Check `snap_dist_m`, `on_drainage_line` and very small catchments.
- CHIRPS cells are about 5 km across, so rainfall exposure varies slowly across space. Differences between nearby springs come mainly from **flow records**. Springs with *Low* confidence are screening results to confirm in the field.
- El Niño usually, but not always, means a dry summer in Lesotho. The tool measures the relationship for each catchment instead of assuming it.
