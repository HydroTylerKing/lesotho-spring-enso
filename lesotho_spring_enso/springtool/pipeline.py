"""End-to-end run: points -> catchments -> rainfall + NDVI + ENSO -> flow response -> susceptibility."""
from __future__ import annotations

import json
from dataclasses import dataclass, field
from pathlib import Path

import geopandas as gpd
import numpy as np
import pandas as pd
from shapely.geometry import Point

from . import analysis as A
from .chirps import catchment_precip
from .config import Settings, WGS84
from .enso import load_oni, hydro_year_phase
from .flows import monthly_flows
from .ndvi import catchment_ndvi
from .watershed import delineate_point

ID_CANDIDATES = ["spring_id", "springid", "id", "spring", "site_id", "site", "name", "code", "fid"]


def prepare_points(gdf: gpd.GeoDataFrame, id_col: str | None = None) -> gpd.GeoDataFrame:
    """Points in EPSG:4326 with a unique string 'spring_id' and optional 'name'."""
    gdf = gdf.copy()
    if gdf.crs is None:
        gdf = gdf.set_crs(WGS84)
    gdf = gdf.to_crs(WGS84)
    gdf = gdf[gdf.geometry.notna() & ~gdf.geometry.is_empty]
    gdf["geometry"] = gdf.geometry.apply(lambda g: g if g.geom_type == "Point" else g.representative_point())
    if id_col is None:
        low = {c.lower(): c for c in gdf.columns}
        id_col = next((low[c] for c in ID_CANDIDATES if c in low), None)
    if id_col is not None:
        gdf["spring_id"] = gdf[id_col].astype(str).str.strip()
    else:
        gdf["spring_id"] = [f"SP{i + 1:03d}" for i in range(len(gdf))]
    if gdf["spring_id"].duplicated().any():
        dup = gdf["spring_id"].duplicated(keep=False)
        gdf.loc[dup, "spring_id"] = gdf.loc[dup, "spring_id"] + "_" + gdf.loc[dup].groupby("spring_id").cumcount().add(1).astype(str)
    if "name" not in gdf.columns:
        gdf["name"] = gdf["spring_id"]
    return gdf.reset_index(drop=True)


@dataclass
class Results:
    springs: gpd.GeoDataFrame
    catchments: gpd.GeoDataFrame
    outlets: gpd.GeoDataFrame
    monthly: pd.DataFrame        # long: spring_id, month, ...
    hydro: pd.DataFrame          # long: spring_id, hydro_year, ...
    enso_years: pd.DataFrame
    info: dict = field(default_factory=dict)

    # ---------------------------------------------------------- exports
    def to_gpkg(self, path):
        path = Path(path)
        if path.exists():
            path.unlink()
        self.springs.to_file(path, layer="springs_susceptibility", driver="GPKG")
        self.catchments.to_file(path, layer="catchments", driver="GPKG")
        self.outlets.to_file(path, layer="snapped_outlets", driver="GPKG")
        return path

    def to_excel(self, path):
        with pd.ExcelWriter(path, engine="openpyxl") as xw:
            readme = pd.DataFrame({"item": list(self.info), "value": [json.dumps(v, default=str) if not isinstance(v, str) else v
                                                                      for v in self.info.values()]})
            readme.to_excel(xw, sheet_name="run_info", index=False)
            self.springs.drop(columns="geometry").assign(
                lon=self.springs.geometry.x, lat=self.springs.geometry.y).to_excel(xw, sheet_name="susceptibility", index=False)
            self.hydro.to_excel(xw, sheet_name="hydrological_years", index=False)
            self.monthly.to_excel(xw, sheet_name="monthly", index=False)
            self.enso_years.reset_index().to_excel(xw, sheet_name="enso_years", index=False)
        return path

    def to_csvs(self, folder):
        folder = Path(folder)
        folder.mkdir(parents=True, exist_ok=True)
        self.monthly.to_csv(folder / "monthly_timeseries.csv", index=False)
        self.hydro.to_csv(folder / "hydrological_years.csv", index=False)
        self.springs.drop(columns="geometry").assign(
            lon=self.springs.geometry.x, lat=self.springs.geometry.y).to_csv(folder / "susceptibility.csv", index=False)
        return folder


def run(points: gpd.GeoDataFrame, settings: Settings, flows: pd.DataFrame | None = None,
        progress=None) -> Results:
    def say(msg, frac=None, stage=(0.0, 1.0)):
        if progress:
            lo, hi = stage
            progress(msg, lo + (hi - lo) * (frac if frac is not None else 0.0))

    pts = prepare_points(points)
    info = {"settings": settings.to_dict(), "n_springs": len(pts)}

    # 1 ── catchments
    cats, errors = {}, {}
    for i, r in pts.iterrows():
        say(f"Delineating catchment {i + 1}/{len(pts)} ({r.spring_id})", i / len(pts), (0.0, 0.25))
        try:
            cats[r.spring_id] = delineate_point(r.spring_id, r.geometry.x, r.geometry.y, settings)
        except Exception as e:
            errors[r.spring_id] = str(e)
    if not cats:
        raise RuntimeError(f"No catchments could be delineated: {errors}")
    info["delineation_errors"] = errors
    geoms = {sid: c.geometry for sid, c in cats.items()}

    # 2 ── rainfall
    say("Reading CHIRPS rainfall", 0, (0.25, 0.45))
    precip, pinfo = catchment_precip(geoms, settings, progress=lambda m, f: say(m, f, (0.25, 0.45)))
    info["precip"] = pinfo

    # 3 ── NDVI
    say("Reading NDVI", 0, (0.45, 0.85))
    ndvi, ninfo = catchment_ndvi(geoms, settings, progress=lambda m, f: say(m, f, (0.45, 0.85)))
    info["ndvi"] = ninfo

    # 4 ── ENSO
    say("Reading ONI (NOAA)", 0, (0.85, 0.9))
    oni = load_oni(settings)
    hy_phase = hydro_year_phase(oni)
    info["oni_last_month"] = str(oni.index.max().date())

    # 5 ── flows and analysis
    fmon = monthly_flows(flows) if flows is not None and len(flows) else {}
    fmon_ci = {k.lower(): v for k, v in fmon.items()}
    unmatched = sorted(set(fmon) - {s for s in cats} - {s.lower() for s in cats})
    info["flow_ids_not_matched"] = unmatched

    weights = dict(rain=settings.w_rain, ndvi=settings.w_ndvi, hydro=settings.w_hydro)
    rows, mons, hys = [], [], []
    for j, (sid, c) in enumerate(cats.items()):
        say(f"Analysing {sid}", j / len(cats), (0.9, 1.0))
        flow = fmon.get(sid) or fmon_ci.get(sid.lower())
        nd = ndvi[sid] if sid in getattr(ndvi, "columns", []) else None
        mon = A.build_monthly(precip[sid], nd, oni, flow)
        hy = A.hydro_table(mon, hy_phase)
        comp = A.enso_composite(hy)
        resp = A.flow_response(mon)
        sc = A.score(comp, resp, c.area_km2, flow[1] if flow else "", weights)
        rows.append(dict(spring_id=sid, catchment_km2=round(c.area_km2, 4), snap_dist_m=round(c.snap_dist_m, 1),
                         on_drainage_line=c.on_drainage_line, catchment_truncated=c.truncated,
                         outlet_elev_m=round(c.outlet_elev_m, 1), mean_elev_m=round(c.mean_elev_m, 1),
                         mean_rain_season_mm=round(float(hy["rain_season_mm"].mean()), 1),
                         flow_kind=flow[1] if flow else "none", **comp, **resp, **sc))
        mons.append(mon.reset_index().assign(spring_id=sid))
        hys.append(hy.reset_index().assign(spring_id=sid))

    res = pd.DataFrame(rows)
    res["priority_rank"] = res["susceptibility"].rank(ascending=False, method="min").astype("Int64")
    springs = pts.merge(res, on="spring_id", how="left")
    springs = gpd.GeoDataFrame(springs, geometry="geometry", crs=WGS84)
    for c in springs.columns:
        if springs[c].dtype == float:
            springs[c] = springs[c].round(4)

    catch_gdf = gpd.GeoDataFrame(
        [dict(spring_id=sid, area_km2=c.area_km2, truncated=c.truncated, geometry=c.geometry) for sid, c in cats.items()],
        geometry="geometry", crs=WGS84).merge(res[["spring_id", "susceptibility", "susceptibility_class", "confidence"]],
                                              on="spring_id", how="left")
    outlets = gpd.GeoDataFrame(
        [dict(spring_id=sid, snap_dist_m=c.snap_dist_m, accum_cells=c.outlet_accum_cells, geometry=c.outlet)
         for sid, c in cats.items()], geometry="geometry", crs=WGS84)

    monthly = pd.concat(mons, ignore_index=True)
    monthly = monthly[["spring_id"] + [c for c in monthly.columns if c != "spring_id"]]
    hydro = pd.concat(hys, ignore_index=True)
    hydro = hydro[["spring_id"] + [c for c in hydro.columns if c != "spring_id"]]
    say("Done", 1.0, (0.0, 1.0))
    return Results(springs, catch_gdf, outlets, monthly, hydro, hy_phase, info)
