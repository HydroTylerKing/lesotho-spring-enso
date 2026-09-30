#!/usr/bin/env python
"""Batch mode: run the whole analysis from the command line.

Example
  python run_batch.py --points springs.gpkg --flows spring_flows.csv --out results

  python run_batch.py --points springs.gpkg --layer springs --id-col code \
      --start 1991-01 --chirps-dir D:/CHIRPS --ndvi none --out results_quick
"""
import argparse
from pathlib import Path

import geopandas as gpd
import pandas as pd

from springtool.config import Settings
from springtool.flows import normalise_flows, read_table
from springtool.maps import results_map
from springtool.pipeline import run, prepare_points


def read_points(path, layer=None):
    p = Path(path)
    if p.suffix.lower() in (".csv", ".xlsx", ".xls"):
        df = read_table(p)
        low = {c.lower(): c for c in df.columns}
        lon = next(low[c] for c in ("lon", "longitude", "x", "long") if c in low)
        lat = next(low[c] for c in ("lat", "latitude", "y") if c in low)
        return gpd.GeoDataFrame(df, geometry=gpd.points_from_xy(df[lon], df[lat]), crs="EPSG:4326")
    return gpd.read_file(p, layer=layer)


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--points", required=True, help="GeoPackage / shapefile / GeoJSON / CSV (lon, lat) of springs")
    ap.add_argument("--layer", help="layer name inside the GeoPackage")
    ap.add_argument("--id-col", help="column holding the spring ID (default: auto-detect)")
    ap.add_argument("--flows", help="CSV/XLSX of spring flow or status records")
    ap.add_argument("--out", default="results", help="output folder")
    ap.add_argument("--start", default="1981-01")
    ap.add_argument("--end", default=None)
    ap.add_argument("--dem", help="local DEM instead of Copernicus GLO-30")
    ap.add_argument("--chirps-dir", help="folder of local monthly CHIRPS GeoTIFFs")
    ap.add_argument("--chirps-version", default="v3", choices=["v3", "v2"])
    ap.add_argument("--ndvi", default="modis", choices=["modis", "local", "none"])
    ap.add_argument("--ndvi-dir", help="folder of local monthly NDVI GeoTIFFs (with --ndvi local)")
    ap.add_argument("--aqua", action="store_true", help="add MODIS Aqua composites (twice the download)")
    ap.add_argument("--oni", help="local copy of NOAA oni.ascii.txt")
    ap.add_argument("--snap", type=float, default=150.0, help="snap radius in metres")
    ap.add_argument("--min-accum", type=int, default=20, help="min upstream cells for a drainage line")
    ap.add_argument("--cache", help="cache folder")
    a = ap.parse_args()

    s = Settings(start=a.start, dem_path=a.dem, chirps_dir=a.chirps_dir, chirps_version=a.chirps_version,
                 ndvi_source=a.ndvi, ndvi_dir=a.ndvi_dir, oni_path=a.oni, snap_radius_m=a.snap,
                 min_accum_cells=a.min_accum, ndvi_platforms=("terra", "aqua") if a.aqua else ("terra",))
    if a.end:
        s.end = a.end
    if a.cache:
        s.cache_dir = a.cache

    pts = read_points(a.points, a.layer)
    if a.id_col:
        pts = prepare_points(pts, a.id_col)
    flows = None
    if a.flows:
        flows, note = normalise_flows(read_table(a.flows))
        print("Flow records:", note)

    last = {"msg": ""}

    def progress(msg, frac):
        if msg.split()[0] != last["msg"].split()[0] if last["msg"] else True:
            print()
        last["msg"] = msg
        print(f"\r[{frac * 100:5.1f}%] {msg:<60}", end="", flush=True)

    res = run(pts, s, flows, progress)
    print()
    out = Path(a.out)
    out.mkdir(parents=True, exist_ok=True)
    res.to_gpkg(out / "spring_enso_susceptibility.gpkg")
    res.to_excel(out / "spring_enso_results.xlsx")
    res.to_csvs(out)
    results_map(res).save(str(out / "susceptibility_map.html"))
    cols = ["spring_id", "name", "catchment_km2", "rain_elnino_mean", "flow_kind", "best_spi",
            "susceptibility", "susceptibility_class", "confidence", "priority_rank"]
    print(res.springs[[c for c in cols if c in res.springs]].sort_values("priority_rank").to_string(index=False))
    if res.info.get("delineation_errors"):
        print("Delineation failed for:", res.info["delineation_errors"])
    if res.info.get("ndvi", {}).get("problem"):
        print("NDVI:", res.info["ndvi"]["problem"])
    if res.info.get("flow_ids_not_matched"):
        print("Flow records with no matching spring point:", res.info["flow_ids_not_matched"])
    print("Outputs written to", out.resolve())


if __name__ == "__main__":
    main()
