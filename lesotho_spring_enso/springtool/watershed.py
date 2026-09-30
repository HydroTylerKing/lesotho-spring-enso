"""Upstream catchment delineation for spring points with pysheds (D8).

For each point:
  1. take a DEM window around it (starting at `window_start_km` half-width),
  2. fill pits/depressions, resolve flats, compute D8 flow direction and accumulation,
  3. snap the point to the highest-accumulation cell within `snap_radius_m`
     (springs sit on or next to drainage lines; GPS and DEM errors are tens of metres),
  4. trace the catchment; if it touches the window edge, double the window and repeat.

Topographic catchments are a proxy for spring recharge areas; they can be
wrong where geology (dolerite dykes, fractured basalt) routes groundwater
across topographic divides. Treat `area_km2` and `snap_dist_m` as quality flags.
"""
from __future__ import annotations

import warnings
from dataclasses import dataclass

import numpy as np
from pyproj import Transformer
from rasterio.features import shapes
from shapely.geometry import shape, Point
from shapely.ops import unary_union, transform as shp_transform

from .config import Settings, UTM_CRS, WGS84
from .dem import dem_window

warnings.filterwarnings("ignore", category=RuntimeWarning, module="pysheds")

# pysheds 0.4 still calls np.in1d, removed in NumPy 2.4
if not hasattr(np, "in1d"):
    np.in1d = lambda ar1, ar2, assume_unique=False, invert=False: np.isin(
        np.ravel(ar1), ar2, assume_unique=assume_unique, invert=invert)


@dataclass
class Catchment:
    spring_id: str
    geometry: object          # shapely polygon, EPSG:4326
    geometry_utm: object      # shapely polygon, EPSG:32735
    outlet: Point             # snapped outlet, EPSG:4326
    area_km2: float
    snap_dist_m: float
    outlet_accum_cells: int
    on_drainage_line: bool    # False if no cell within the snap radius met min_accum_cells
    truncated: bool           # True if catchment still touched the DEM edge at the largest window
    window_km: float
    outlet_elev_m: float
    mean_elev_m: float
    max_elev_m: float


def _condition(dem_path):
    from pysheds.grid import Grid
    grid = Grid.from_raster(str(dem_path))
    dem = grid.read_raster(str(dem_path))
    pit = grid.fill_pits(dem)
    flooded = grid.fill_depressions(pit)
    inflated = grid.resolve_flats(flooded)
    fdir = grid.flowdir(inflated)
    acc = grid.accumulation(fdir)
    return grid, dem, fdir, acc


def _snap(acc: np.ndarray, affine, x: float, y: float, radius_m: float, min_acc: int):
    res = abs(affine.a)
    col, row = ~affine * (x, y)
    r0, c0 = int(np.floor(row)), int(np.floor(col))
    k = int(np.ceil(radius_m / res))
    rs = slice(max(r0 - k, 0), min(r0 + k + 1, acc.shape[0]))
    cs = slice(max(c0 - k, 0), min(c0 + k + 1, acc.shape[1]))
    sub = np.asarray(acc[rs, cs], dtype="float64")
    rr, cc = np.mgrid[rs, cs]
    cx = affine.c + (cc + 0.5) * affine.a
    cy = affine.f + (rr + 0.5) * affine.e
    dist = np.hypot(cx - x, cy - y)
    inside = dist <= max(radius_m, res * 0.75)
    sub = np.where(inside & np.isfinite(sub), sub, -1)
    ok = sub >= min_acc
    pick = ok if ok.any() else sub >= 0
    # highest accumulation; ties broken by distance
    score = np.where(pick, sub - dist * 1e-6, -np.inf)
    i = np.unravel_index(np.argmax(score), score.shape)
    return int(rr[i]), int(cc[i]), float(dist[i]), int(sub[i]), bool(ok.any())


def delineate_point(spring_id: str, lon: float, lat: float, settings: Settings) -> Catchment:
    to_utm = Transformer.from_crs(WGS84, UTM_CRS, always_xy=True)
    to_ll = Transformer.from_crs(UTM_CRS, WGS84, always_xy=True)
    x, y = to_utm.transform(lon, lat)

    half = settings.window_start_km
    while True:
        dem_path = dem_window(lon, lat, half, settings)
        grid, dem, fdir, acc = _condition(dem_path)
        aff = grid.affine
        row, col, dist, acc_cells, on_line = _snap(acc, aff, x, y, settings.snap_radius_m,
                                                   settings.min_accum_cells)
        catch = grid.catchment(x=col, y=row, fdir=fdir, xytype="index")
        cmask = np.asarray(catch, dtype=bool)
        touches = cmask[0, :].any() or cmask[-1, :].any() or cmask[:, 0].any() or cmask[:, -1].any()
        if not touches or half * 2 > settings.window_max_km:
            break
        half *= 2

    polys = [shape(g) for g, v in shapes(cmask.astype("uint8"), mask=cmask, transform=aff,
                                         connectivity=8) if v == 1]
    if not polys:  # single-cell outlet: use the cell itself
        from shapely.geometry import box
        cx0 = aff.c + col * aff.a
        cy0 = aff.f + row * aff.e
        polys = [box(cx0, cy0 + aff.e, cx0 + aff.a, cy0)]
    geom_utm = unary_union(polys)
    geom_utm = geom_utm.buffer(0)
    geom_ll = shp_transform(lambda a, b, z=None: to_ll.transform(a, b), geom_utm)

    ox = aff.c + (col + 0.5) * aff.a
    oy = aff.f + (row + 0.5) * aff.e
    olon, olat = to_ll.transform(ox, oy)

    dem_a = np.asarray(dem, dtype="float64")
    dem_a = np.where(dem_a <= -9000, np.nan, dem_a)
    elev = dem_a[cmask]
    cell_km2 = abs(aff.a * aff.e) / 1e6
    return Catchment(
        spring_id=str(spring_id),
        geometry=geom_ll,
        geometry_utm=geom_utm,
        outlet=Point(olon, olat),
        area_km2=float(cmask.sum() * cell_km2),
        snap_dist_m=dist,
        outlet_accum_cells=acc_cells,
        on_drainage_line=on_line,
        truncated=bool(touches),
        window_km=half,
        outlet_elev_m=float(dem_a[row, col]),
        mean_elev_m=float(np.nanmean(elev)) if elev.size else np.nan,
        max_elev_m=float(np.nanmax(elev)) if elev.size else np.nan,
    )
