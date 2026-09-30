#!/usr/bin/env python
"""Build a GeoPackage of spring locations in Lesotho from OpenStreetMap.

    python get_osm_springs.py                       -> lesotho_osm_springs.gpkg
    python get_osm_springs.py --out my_springs.gpkg

Layers written
  springs       natural=spring and man_made=spring_box (captured / protected springs)
  water_points  amenity=drinking_water, man_made=water_tap, man_made=water_well
                (village taps are often fed by protected springs; useful context)

Source: OpenStreetMap contributors (ODbL). The Overpass API is tried first
(three public servers); if all are busy, it falls back to Geofabrik's Lesotho
extract, which is updated daily but only carries springs, not spring boxes
or taps.

OSM coverage of springs in Lesotho is uneven: it depends on who has mapped
where. Treat this as a starting layer to merge with DRWS / district water
office inventories and your own GPS field data, not as a complete inventory.
The `springs` layer can be loaded straight into the spring El Niño tool
(spring_id and name columns are filled in).
"""
from __future__ import annotations

import argparse
import json
import sys
import tempfile
import urllib.parse
import urllib.request
from pathlib import Path

import geopandas as gpd
import pandas as pd
from shapely.geometry import Point

UA = "lesotho-spring-enso/1.0 (AWEP training; OSM spring extract)"
OVERPASS = ["https://overpass-api.de/api/interpreter",
            "https://overpass.kumi.systems/api/interpreter",
            "https://overpass.private.coffee/api/interpreter"]
QUERY = """
[out:json][timeout:240];
area["ISO3166-1"="LS"][admin_level=2]->.ls;
(
  nwr["natural"="spring"](area.ls);
  nwr["man_made"="spring_box"](area.ls);
  nwr["amenity"="drinking_water"](area.ls);
  nwr["man_made"~"^(water_tap|water_well)$"](area.ls);
);
out center tags;
"""
GEOFABRIK = "https://download.geofabrik.de/africa/lesotho-latest-free.shp.zip"
KEEP_TAGS = ["name", "name:en", "name:st", "drinking_water", "seasonal", "intermittent", "access",
             "operator", "pump", "capacity", "ele", "description", "note", "source", "check_date"]


def category(tags: dict) -> str:
    if tags.get("natural") == "spring":
        return "spring"
    if tags.get("man_made") == "spring_box":
        return "spring_box"
    if tags.get("amenity") == "drinking_water":
        return "drinking_water"
    return tags.get("man_made", "other")


def from_overpass() -> gpd.GeoDataFrame:
    data = urllib.parse.urlencode({"data": QUERY}).encode()
    last = None
    for url in OVERPASS:
        print(f"Querying {url} ...")
        try:
            req = urllib.request.Request(url, data=data, headers={"User-Agent": UA})
            with urllib.request.urlopen(req, timeout=300) as r:
                js = json.load(r)
            break
        except Exception as e:
            print(f"  failed: {type(e).__name__}: {e}")
            last = e
    else:
        raise RuntimeError(f"All Overpass servers failed ({last})")
    rows = []
    for el in js.get("elements", []):
        if el["type"] == "node":
            lon, lat = el["lon"], el["lat"]
        elif "center" in el:
            lon, lat = el["center"]["lon"], el["center"]["lat"]
        else:
            continue
        t = el.get("tags", {})
        row = dict(osm_type=el["type"], osm_id=el["id"], category=category(t),
                   **{k.replace(":", "_"): t.get(k) for k in KEEP_TAGS},
                   osm_tags=json.dumps(t, ensure_ascii=False), geometry=Point(lon, lat))
        rows.append(row)
    return gpd.GeoDataFrame(rows, geometry="geometry", crs="EPSG:4326") if rows else \
        gpd.GeoDataFrame(columns=["osm_type", "osm_id", "category", "geometry"], geometry="geometry", crs="EPSG:4326")


def from_geofabrik() -> gpd.GeoDataFrame:
    tmp = Path(tempfile.gettempdir()) / "lesotho-latest-free.shp.zip"
    print(f"Downloading {GEOFABRIK} (about 30-60 MB) ...")
    req = urllib.request.Request(GEOFABRIK, headers={"User-Agent": UA})
    with urllib.request.urlopen(req, timeout=600) as r, open(tmp, "wb") as fh:
        while chunk := r.read(1 << 20):
            fh.write(chunk)
    nat = gpd.read_file(f"zip://{tmp}!gis_osm_natural_free_1.shp")
    sp = nat[nat["fclass"] == "spring"].copy()
    parts = [gpd.GeoDataFrame(dict(osm_type="node", osm_id=sp["osm_id"].astype("int64"), category="spring",
                                   name=sp["name"]), geometry=sp.geometry, crs=nat.crs)]
    try:
        pois = gpd.read_file(f"zip://{tmp}!gis_osm_pois_free_1.shp")
        wp = pois[pois["fclass"].isin(["drinking_water", "water_well"])].copy()
        parts.append(gpd.GeoDataFrame(dict(osm_type="node", osm_id=wp["osm_id"].astype("int64"),
                                           category=wp["fclass"], name=wp["name"]), geometry=wp.geometry, crs=pois.crs))
    except Exception:
        pass
    return gpd.GeoDataFrame(pd.concat(parts, ignore_index=True), crs="EPSG:4326")


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--out", default="lesotho_osm_springs.gpkg")
    ap.add_argument("--geofabrik", action="store_true", help="skip Overpass and use the Geofabrik extract")
    a = ap.parse_args()

    if a.geofabrik:
        gdf, src = from_geofabrik(), "Geofabrik Lesotho extract"
    else:
        try:
            gdf, src = from_overpass(), "Overpass API"
        except Exception as e:
            print(e, "\nFalling back to Geofabrik.")
            gdf, src = from_geofabrik(), "Geofabrik Lesotho extract"

    if gdf.empty:
        sys.exit("No springs or water points found.")
    gdf = gdf.drop_duplicates(["osm_type", "osm_id"]).reset_index(drop=True)
    prefix = gdf["osm_type"].str[0].str.upper()
    gdf.insert(0, "spring_id", "OSM-" + prefix + gdf["osm_id"].astype(str))
    nm = gdf["name"] if "name" in gdf else pd.Series(None, index=gdf.index)
    gdf["name"] = nm.where(nm.notna() & (nm.astype(str) != ""), gdf["category"].str.replace("_", " ").str.title()
                           + " " + gdf["osm_id"].astype(str))
    gdf["lon"], gdf["lat"] = gdf.geometry.x.round(6), gdf.geometry.y.round(6)
    gdf["source"] = gdf.get("source", pd.Series(None, index=gdf.index)).fillna("") \
        .radd("OpenStreetMap contributors (ODbL) via " + src + "; osm source tag: ")

    springs = gdf[gdf["category"].isin(["spring", "spring_box"])]
    water = gdf[~gdf["category"].isin(["spring", "spring_box"])]
    out = Path(a.out)
    if out.exists():
        out.unlink()
    if len(springs):
        springs.to_file(out, layer="springs", driver="GPKG")
    if len(water):
        water.to_file(out, layer="water_points", driver="GPKG")
    print(f"\nSource: {src}")
    print(gdf["category"].value_counts().to_string())
    print(f"\nWrote {out.resolve()}")
    print(f"  layer 'springs'      : {len(springs)} points (load this one into the spring El Niño tool)")
    print(f"  layer 'water_points' : {len(water)} points")
    print("OSM coverage is uneven: merge with DRWS / district inventories and field GPS before relying on it.")


if __name__ == "__main__":
    main()
