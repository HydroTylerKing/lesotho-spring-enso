"""Monthly CHIRPS rainfall for catchments.

Downloads only the Lesotho window of each global monthly GeoTIFF (HTTP range
reads through GDAL /vsicurl/), caches it as a small .npy, then computes
area-weighted catchment totals (mm/month).

Sources
  v3: https://data.chc.ucsb.edu/products/CHIRPS/v3.0/monthly/global/tifs/chirps-v3.0.YYYY.MM.tif
  v2: https://data.chc.ucsb.edu/products/CHIRPS-2.0/global_monthly/tifs/chirps-v2.0.YYYY.MM.tif.gz
      (CHC stops producing v2 after December 2026)
Local folder: any GeoTIFFs whose names contain YYYY.MM, YYYY-MM or YYYY_MM.
"""
from __future__ import annotations

import json
import re
from concurrent.futures import ThreadPoolExecutor
from pathlib import Path

import numpy as np
import pandas as pd
import rasterio
from affine import Affine
from rasterio.windows import from_bounds, Window

from .config import Settings
from .zonal import coverage_weights, window_for, weighted_mean

URLS = {
    "v3": "https://data.chc.ucsb.edu/products/CHIRPS/v3.0/monthly/global/tifs/chirps-v3.0.{y}.{m:02d}.tif",
    "v2": "https://data.chc.ucsb.edu/products/CHIRPS-2.0/global_monthly/tifs/chirps-v2.0.{y}.{m:02d}.tif.gz",
}
RES = 0.05
GDAL_ENV = dict(GDAL_DISABLE_READDIR_ON_OPEN="EMPTY_DIR", CPL_VSIL_CURL_ALLOWED_EXTENSIONS=".tif,.gz",
                GDAL_HTTP_MAX_RETRY="4", GDAL_HTTP_RETRY_DELAY="2", VSI_CACHE="TRUE")


def snap_region(bounds, pad=0.1):
    """Snap lon/lat bounds outward to the 0.05° CHIRPS grid."""
    x0 = np.floor((bounds[0] - pad) / RES) * RES
    y0 = np.floor((bounds[1] - pad) / RES) * RES
    x1 = np.ceil((bounds[2] + pad) / RES) * RES
    y1 = np.ceil((bounds[3] + pad) / RES) * RES
    return (round(x0, 4), round(y0, 4), round(x1, 4), round(y1, 4))


def _local_index(folder: str) -> dict:
    pat = re.compile(r"(19|20)(\d{2})[._-](\d{2})(?!\d)")
    idx = {}
    for p in Path(folder).rglob("*"):
        if p.suffix.lower() in (".tif", ".tiff", ".gz") and (m := pat.search(p.name)):
            idx[(int(m.group(1) + m.group(2)), int(m.group(3)))] = p
    return idx


def _read_month(src_path: str, region, shape):
    x0, y0, x1, y1 = region
    with rasterio.Env(**GDAL_ENV):
        with rasterio.open(src_path) as src:
            win = from_bounds(x0, y0, x1, y1, transform=src.transform).round_offsets().round_lengths()
            arr = src.read(1, window=win, boundless=True, fill_value=-9999).astype("float32")
    arr = arr[: shape[0], : shape[1]]
    arr[(arr < 0) | ~np.isfinite(arr)] = np.nan
    return arr


def _download_and_read(url: str, region, shape, tmpdir: Path):
    """Fallback when the server does not allow range reads: fetch the file, read the window, delete it."""
    import urllib.request
    import uuid
    tmp = Path(tmpdir) / f"_dl_{uuid.uuid4().hex}{'.tif.gz' if url.endswith('.gz') else '.tif'}"
    try:
        with urllib.request.urlopen(url, timeout=120) as r, open(tmp, "wb") as fh:
            ctype = r.headers.get("Content-Type", "")
            if "html" in ctype.lower():
                raise RuntimeError(f"server returned a web page ({ctype}), not a GeoTIFF")
            while chunk := r.read(1 << 20):
                fh.write(chunk)
        src = f"/vsigzip/{tmp}" if url.endswith(".gz") else str(tmp)
        return _read_month(src, region, shape)
    finally:
        tmp.unlink(missing_ok=True)


def load_cube(region, months: pd.DatetimeIndex, settings: Settings, progress=None):
    """Return (cube[t,y,x] mm, transform, months_available)."""
    x0, y0, x1, y1 = region
    shape = (int(round((y1 - y0) / RES)), int(round((x1 - x0) / RES)))
    transform = Affine(RES, 0, x0, 0, -RES, y1)
    tag = "local" if settings.chirps_dir else settings.chirps_version
    cdir = settings.cache("chirps", tag, "_".join(f"{v:.2f}" for v in region))
    local = _local_index(settings.chirps_dir) if settings.chirps_dir else None

    def job(ts):
        f = cdir / f"{ts.year}.{ts.month:02d}.npy"
        if f.exists():
            return ts, np.load(f), None
        if local is not None:
            p = local.get((ts.year, ts.month))
            if p is None:
                return ts, None, None
            src = f"/vsigzip/{p}" if p.suffix == ".gz" else str(p)
            try:
                arr = _read_month(src, region, shape)
            except Exception as e:
                return ts, None, f"{p.name}: {type(e).__name__}: {e}"
        else:
            url = URLS[settings.chirps_version].format(y=ts.year, m=ts.month)
            try:  # fast path: read only the Lesotho window over HTTP
                src = ("/vsigzip//vsicurl/" if url.endswith(".gz") else "/vsicurl/") + url
                arr = _read_month(src, region, shape)
            except Exception:
                try:  # slow path: download the whole file, read the window, delete it
                    arr = _download_and_read(url, region, shape, cdir)
                except Exception as e:
                    return ts, None, f"{url}: {type(e).__name__}: {e}"
        np.save(f, arr)
        return ts, arr, None

    out, errors = {}, []
    with ThreadPoolExecutor(max_workers=settings.workers) as ex:
        for i, (ts, arr, err) in enumerate(ex.map(job, months)):
            if arr is not None:
                out[ts] = arr
            elif err:
                errors.append(err)
            if progress:
                progress(f"CHIRPS {ts:%Y-%m}", (i + 1) / len(months))
    if not out:
        raise RuntimeError("No CHIRPS months could be read. First error: "
                           + (errors[0] if errors else "no matching files in the local CHIRPS folder")
                           + ". Run `python check_connection.py`, or set a local CHIRPS folder.")
    got = pd.DatetimeIndex(sorted(out))
    missing = months.difference(got)
    # trailing months not yet published are normal; gaps in the middle are reported
    gaps = missing[missing < got.max()]
    cube = np.stack([out[t] for t in got])
    return cube, transform, got, list(gaps)


def catchment_precip(catchments: dict, settings: Settings, progress=None):
    """catchments: {spring_id: shapely polygon EPSG:4326}. Returns (DataFrame months x springs, info)."""
    from shapely.ops import unary_union
    months = pd.date_range(settings.start + "-01", settings.end + "-01", freq="MS")
    region = snap_region(unary_union(list(catchments.values())).bounds)
    cube, transform, got, gaps = load_cube(region, months, settings, progress)
    shape = cube.shape[1:]
    series = {}
    for sid, geom in catchments.items():
        win = window_for(geom, transform, shape)
        if win is None:
            series[sid] = pd.Series(np.nan, index=got)
            continue
        wt = rasterio.windows.transform(win, transform)
        sub = cube[:, win.row_off: win.row_off + win.height, win.col_off: win.col_off + win.width]
        w = coverage_weights(geom, wt, (win.height, win.width), factor=10)
        mean, frac = weighted_mean(sub, w)
        s = pd.Series(mean, index=got)
        s[frac < 0.5] = np.nan
        series[sid] = s
    df = pd.DataFrame(series)
    df.index.name = "month"
    return df, {"region": region, "months": len(got), "gaps": [str(g.date()) for g in gaps],
                "source": "local folder" if settings.chirps_dir else f"CHIRPS {settings.chirps_version}"}
