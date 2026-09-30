"""DEM access: Copernicus GLO-30 from Microsoft Planetary Computer, or a local DEM.

Every window is reprojected to UTM 35S at `dem_resolution_m` and cached as a
GeoTIFF so repeated runs on the same springs do not download again.
"""
from __future__ import annotations

import hashlib
from pathlib import Path

import numpy as np
import rasterio
from affine import Affine
from pyproj import Transformer
from rasterio.warp import reproject, Resampling, transform_bounds
from rasterio.windows import from_bounds

from .config import Settings, UTM_CRS, WGS84, PC_STAC

NODATA = -9999.0


def _target_grid(x: float, y: float, half_m: float, res: float):
    x0 = np.floor((x - half_m) / res) * res
    y1 = np.ceil((y + half_m) / res) * res
    n = int(round(2 * half_m / res))
    transform = Affine(res, 0, x0, 0, -res, y1)
    return transform, (n, n)


def _dem_sources(lonlat_bounds, settings: Settings) -> list[str]:
    if settings.dem_path:
        return [settings.dem_path]
    import planetary_computer
    import pystac_client
    cat = pystac_client.Client.open(PC_STAC, modifier=planetary_computer.sign_inplace)
    items = list(cat.search(collections=["cop-dem-glo-30"], bbox=list(lonlat_bounds)).items())
    if not items:
        raise RuntimeError("No Copernicus DEM tiles found for this area.")
    return [it.assets["data"].href for it in items]


def dem_window(lon: float, lat: float, half_km: float, settings: Settings) -> Path:
    """Return path to a cached UTM DEM GeoTIFF centred on (lon, lat)."""
    to_utm = Transformer.from_crs(WGS84, UTM_CRS, always_xy=True)
    x, y = to_utm.transform(lon, lat)
    res = settings.dem_resolution_m
    transform, shape = _target_grid(x, y, half_km * 1000.0, res)

    src_tag = settings.dem_path or "copdem30"
    key = hashlib.md5(f"{src_tag}|{transform.c:.0f}|{transform.f:.0f}|{shape}|{res}".encode()).hexdigest()[:16]
    out = settings.cache("dem") / f"dem_{key}.tif"
    if out.exists():
        return out

    left, top = transform.c, transform.f
    right, bottom = left + shape[1] * res, top - shape[0] * res
    lonlat_bounds = transform_bounds(UTM_CRS, WGS84, left, bottom, right, top, densify_pts=21)

    dst = np.full(shape, np.nan, dtype="float32")
    for href in _dem_sources(lonlat_bounds, settings):
        with rasterio.open(href) as src:
            sb = transform_bounds(UTM_CRS, src.crs, left, bottom, right, top, densify_pts=21)
            # pad by a few source pixels so bilinear resampling has neighbours at the edge
            px = abs(src.transform.a) * 3
            win = from_bounds(sb[0] - px, sb[1] - px, sb[2] + px, sb[3] + px, transform=src.transform)
            win = win.round_offsets().round_lengths()
            win = win.intersection(rasterio.windows.Window(0, 0, src.width, src.height))
            if win.width <= 0 or win.height <= 0:
                continue
            arr = src.read(1, window=win).astype("float32")
            nod = src.nodata
            if nod is not None:
                arr[arr == nod] = np.nan
            tmp = np.full(shape, np.nan, dtype="float32")
            reproject(arr, tmp, src_transform=src.window_transform(win), src_crs=src.crs,
                      src_nodata=np.nan, dst_transform=transform, dst_crs=UTM_CRS,
                      dst_nodata=np.nan, resampling=Resampling.bilinear)
            fill = np.isnan(dst) & ~np.isnan(tmp)
            dst[fill] = tmp[fill]

    if np.isnan(dst).all():
        raise RuntimeError(f"DEM has no data around lon={lon:.4f}, lat={lat:.4f}.")
    dst = np.where(np.isnan(dst), NODATA, dst).astype("float32")
    profile = dict(driver="GTiff", height=shape[0], width=shape[1], count=1, dtype="float32",
                   crs=UTM_CRS, transform=transform, nodata=NODATA, compress="deflate", tiled=True)
    with rasterio.open(out, "w", **profile) as ds:
        ds.write(dst, 1)
    return out
