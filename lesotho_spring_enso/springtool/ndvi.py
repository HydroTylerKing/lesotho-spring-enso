"""Monthly NDVI for catchments.

MODIS source: MOD13Q1 (Terra) / MYD13Q1 (Aqua) v6.1, 250 m, 16-day composites,
read as Cloud-Optimised GeoTIFFs from Microsoft Planetary Computer (no account
needed). Pixels are kept when pixel_reliability is 0 (good) or 1 (marginal);
snow (2) and cloud (3) are dropped. Each composite gets an area-weighted mean
per catchment; composites are averaged into calendar months by mid-date.

Local source: a folder of monthly NDVI GeoTIFFs with YYYY.MM / YYYY-MM /
YYYY_MM in their names (any CRS; integer files are assumed scaled by 10 000).

Results per catchment and composite are cached, so re-running (or adding
springs) only fetches what is new.
"""
from __future__ import annotations

import hashlib
import re
from concurrent.futures import ThreadPoolExecutor
from pathlib import Path

import numpy as np
import pandas as pd
import rasterio
from pyproj import Transformer
from shapely.ops import transform as shp_transform, unary_union

from .config import Settings, PC_STAC
from .zonal import coverage_weights, window_for, weighted_mean

COLLECTION = "modis-13Q1-061"
NDVI_ASSET = "250m_16_days_NDVI"
REL_ASSET = "250m_16_days_pixel_reliability"
GDAL_ENV = dict(GDAL_DISABLE_READDIR_ON_OPEN="EMPTY_DIR",
                GDAL_HTTP_MAX_RETRY="4", GDAL_HTTP_RETRY_DELAY="2", GDAL_CACHEMAX=256)


def geom_key(geom) -> str:
    return hashlib.md5(geom.simplify(1e-5).wkb).hexdigest()[:16]


def _to_crs(geom, crs):
    t = Transformer.from_crs("EPSG:4326", crs, always_xy=True)
    return shp_transform(lambda x, y, z=None: t.transform(x, y), geom)


def _zonal_from_dataset(ds_ndvi, ds_rel, geoms: dict):
    """Area-weighted NDVI and good-quality fraction for each geometry (4326) in one composite."""
    out = {}
    shape = (ds_ndvi.height, ds_ndvi.width)
    for sid, g in geoms.items():
        gp = _to_crs(g, ds_ndvi.crs)
        win = window_for(gp, ds_ndvi.transform, shape)
        if win is None:
            continue
        nd = ds_ndvi.read(1, window=win).astype("float64")
        wt = ds_ndvi.window_transform(win)
        w = coverage_weights(gp, wt, (win.height, win.width), factor=6)
        if w.sum() == 0:
            continue
        valid = (nd >= -2000) & (nd <= 10000)
        if ds_rel is not None:
            rel = ds_rel.read(1, window=win)
            valid &= (rel == 0) | (rel == 1)
        elif ds_ndvi.nodata is not None:
            valid &= nd != ds_ndvi.nodata
        m, frac = weighted_mean(nd * 1e-4, w, valid)
        out[sid] = (float(m), float(frac), float(w.sum()))
    return out


# ------------------------------------------------------------------ MODIS
def _modis_items(bbox, start, end, platforms):
    import planetary_computer
    import pystac_client
    cat = pystac_client.Client.open(PC_STAC)
    items = cat.search(collections=[COLLECTION], bbox=list(bbox), datetime=f"{start}/{end}").items()
    keep = []
    for it in items:
        plat = str(it.properties.get("platform", "")).lower()
        if not platforms or any(p in plat for p in platforms):
            keep.append(it)
    return keep


def _item_mid(it) -> pd.Timestamp:
    p = it.properties
    if p.get("start_datetime") and p.get("end_datetime"):
        a = pd.Timestamp(p["start_datetime"]).tz_localize(None)
        b = pd.Timestamp(p["end_datetime"]).tz_localize(None)
        return a + (b - a) / 2
    return pd.Timestamp(it.datetime).tz_localize(None) + pd.Timedelta(days=8)


def _read_item(it, need: dict):
    """Area-weighted NDVI for the catchments in `need` from one MODIS composite (one retry)."""
    import planetary_computer
    last = None
    for attempt in range(2):
        try:
            with rasterio.Env(**GDAL_ENV):
                h_nd = planetary_computer.sign(it.assets[NDVI_ASSET].href)
                h_rel = planetary_computer.sign(it.assets[REL_ASSET].href)
                with rasterio.open(h_nd) as dn, rasterio.open(h_rel) as dr:
                    return _zonal_from_dataset(dn, dr, need), None
        except Exception as e:
            last = f"{type(e).__name__}: {e}"
    return None, last


def _modis_composites(geoms: dict, settings: Settings, progress=None):
    cdir = settings.cache("ndvi", "modis")
    keys = {sid: geom_key(g) for sid, g in geoms.items()}
    cached = {}
    for sid, k in keys.items():
        f = cdir / f"{k}.csv"
        cached[sid] = pd.read_csv(f) if f.exists() else pd.DataFrame(
            columns=["item", "date", "platform", "ndvi", "valid_frac", "weight"])

    bbox = unary_union(list(geoms.values())).bounds
    start = max(settings.start, settings.ndvi_start) + "-01"
    end = pd.Timestamp(settings.end + "-01") + pd.offsets.MonthEnd(0)
    if progress:
        progress("Searching MODIS NDVI catalogue", 0.0)
    stats = dict(items_found=0, items_read=0, items_failed=0, items_cached=0, first_error=None)
    try:
        items = _modis_items(bbox, start, f"{end:%Y-%m-%d}", settings.ndvi_platforms)
    except Exception as e:
        stats["first_error"] = f"Catalogue search failed: {type(e).__name__}: {e}"
        items = []
    stats["items_found"] = len(items)

    todo = []
    for it in items:
        need = {sid: g for sid, g in geoms.items() if it.id not in set(cached[sid]["item"])}
        if need:
            todo.append((it, need))
    stats["items_cached"] = len(items) - len(todo)

    new_rows = {sid: [] for sid in geoms}
    with ThreadPoolExecutor(max_workers=settings.workers) as ex:
        futures = ex.map(lambda a: (a[0], a[1]) + _read_item(*a), todo)
        for i, (it, need, res, err) in enumerate(futures):
            if res is None:
                stats["items_failed"] += 1
                stats["first_error"] = stats["first_error"] or err
            else:
                stats["items_read"] += 1
                mid = _item_mid(it)
                plat = it.properties.get("platform", "") or "unknown"
                for sid in need:
                    if sid in res:
                        m, frac, w = res[sid]
                        new_rows[sid].append(dict(item=it.id, date=mid, platform=plat,
                                                  ndvi=m, valid_frac=frac, weight=w))
                    else:  # catchment not on this MODIS tile: remember so we skip it next time
                        new_rows[sid].append(dict(item=it.id, date=mid, platform=plat,
                                                  ndvi=np.nan, valid_frac=0.0, weight=0.0))
            if progress:
                progress(f"MODIS NDVI composite {i + 1}/{len(todo)}"
                         + (f" ({stats['items_failed']} failed)" if stats["items_failed"] else ""),
                         (i + 1) / max(len(todo), 1))

    frames = []
    for sid in geoms:
        df = pd.concat([cached[sid], pd.DataFrame(new_rows[sid])], ignore_index=True) \
            if new_rows[sid] else cached[sid]
        df.to_csv(cdir / f"{keys[sid]}.csv", index=False)
        frames.append(df.assign(spring_id=sid))
    out = pd.concat(frames, ignore_index=True)
    out["platform"] = out["platform"].fillna("unknown")
    return out, stats


# ------------------------------------------------------------------ local
def _local_composites(geoms: dict, settings: Settings, progress=None) -> pd.DataFrame:
    pat = re.compile(r"(19|20)(\d{2})[._-](\d{2})(?!\d)")
    files = sorted(p for p in Path(settings.ndvi_dir).rglob("*") if p.suffix.lower() in (".tif", ".tiff")
                   and pat.search(p.name))
    rows = []
    for i, p in enumerate(files):
        m = pat.search(p.name)
        date = pd.Timestamp(int(m.group(1) + m.group(2)), int(m.group(3)), 15)
        with rasterio.open(p) as ds:
            scale_int = np.issubdtype(np.dtype(ds.dtypes[0]), np.integer)
            shape = (ds.height, ds.width)
            for sid, g in geoms.items():
                gp = _to_crs(g, ds.crs)
                win = window_for(gp, ds.transform, shape)
                if win is None:
                    continue
                arr = ds.read(1, window=win).astype("float64")
                valid = np.isfinite(arr)
                if ds.nodata is not None:
                    valid &= arr != ds.nodata
                if scale_int:
                    arr = arr * 1e-4
                valid &= (arr >= -0.2) & (arr <= 1.0)
                w = coverage_weights(gp, ds.window_transform(win), (win.height, win.width), factor=6)
                mean, frac = weighted_mean(arr, w, valid)
                rows.append(dict(spring_id=sid, item=p.name, date=date, platform="local",
                                 ndvi=float(mean), valid_frac=float(frac), weight=float(w.sum())))
        if progress:
            progress(f"NDVI {p.name}", (i + 1) / len(files))
    return pd.DataFrame(rows)


def catchment_ndvi(catchments: dict, settings: Settings, progress=None):
    """Monthly NDVI (months x springs) and info dict."""
    if settings.ndvi_source == "none":
        return pd.DataFrame(), {"source": "none"}
    if settings.ndvi_source == "local":
        comp = _local_composites(catchments, settings, progress)
        src = "local folder"
        stats = {}
    else:
        comp, stats = _modis_composites(catchments, settings, progress)
        src = "MODIS " + "+".join(settings.ndvi_platforms) + " 250 m (MOD13Q1/MYD13Q1 v6.1)"
    info = {"source": src, **stats}
    if comp.empty:
        info.update(composites=0, problem=_explain(info, 0, 0))
        return pd.DataFrame(), info
    comp["date"] = pd.to_datetime(comp["date"])
    # merge pieces of a catchment that fall on two MODIS tiles (same composite date & platform)
    comp["good_w"] = comp["valid_frac"].fillna(0) * comp["weight"].fillna(0)
    comp["nd_w"] = comp["ndvi"].fillna(0) * comp["good_w"]
    g = comp.groupby(["spring_id", "date", "platform"], as_index=False)[["good_w", "nd_w", "weight"]].sum()
    with np.errstate(invalid="ignore", divide="ignore"):
        g["ndvi"] = np.where(g["good_w"] > 0, g["nd_w"] / g["good_w"], np.nan)
        g["valid_frac"] = np.where(g["weight"] > 0, g["good_w"] / g["weight"], 0.0)
    on_tile = g[g["weight"] > 0]
    comp = g[(g["valid_frac"] >= settings.ndvi_min_valid) & g["ndvi"].notna()].copy()
    info.update(composites_on_catchments=int(len(on_tile)), composites=int(len(comp)))
    if comp.empty:
        info["problem"] = _explain(info, len(on_tile), 0)
        return pd.DataFrame(), info
    comp["month"] = comp["date"].dt.to_period("M").dt.to_timestamp()
    monthly = comp.groupby(["month", "spring_id"])["ndvi"].mean().unstack("spring_id")
    monthly = monthly.reindex(columns=list(catchments))
    info["months_with_ndvi"] = int(monthly.notna().any(axis=1).sum())
    if info.get("items_failed"):
        info["problem"] = _explain(info, len(on_tile), len(comp))
    return monthly, info


def _explain(info: dict, on_tile: int, kept: int) -> str:
    if info.get("first_error") and not info.get("items_read") and not info.get("items_cached"):
        return ("No MODIS images could be read. First error: " + str(info["first_error"]) +
                " — run `python check_connection.py` to test each data source.")
    if info.get("items_found") == 0 and info["source"].startswith("MODIS"):
        return "The MODIS catalogue returned no images for these catchments and dates."
    if on_tile == 0:
        return "MODIS images were read but none covered the catchments."
    if kept == 0:
        return ("Every MODIS composite was rejected by the quality filter (cloud/snow). "
                "Lower 'Min. good-quality share per composite' in the NDVI settings.")
    return (f"{info.get('items_failed')} of {info.get('items_found')} MODIS images failed to download "
            f"(first error: {info.get('first_error')}). Run again to retry them; the rest are cached.")
