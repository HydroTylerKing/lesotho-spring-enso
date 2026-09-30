"""Area-weighted zonal means for catchments that may be smaller than one pixel."""
from __future__ import annotations

import numpy as np
from affine import Affine
from rasterio.features import rasterize
from rasterio.windows import from_bounds, Window


def coverage_weights(geom, transform: Affine, shape: tuple[int, int], factor: int = 8) -> np.ndarray:
    """Fraction of each pixel covered by `geom` (same CRS as `transform`).

    Rasterises on a grid `factor` times finer and block-averages, so a catchment
    that covers a third of a 5 km CHIRPS cell gets weight ~0.33 there.
    If the geometry is too small to hit any fine cell, the pixel under its
    representative point gets weight 1.
    """
    h, w = shape
    fine_t = transform * Affine.scale(1.0 / factor)
    fine = rasterize([(geom, 1)], out_shape=(h * factor, w * factor), transform=fine_t,
                     fill=0, all_touched=False, dtype="uint8")
    wts = fine.reshape(h, factor, w, factor).mean(axis=(1, 3)).astype("float64")
    if wts.sum() == 0:
        p = geom.representative_point()
        col, row = ~transform * (p.x, p.y)
        r, c = int(np.floor(row)), int(np.floor(col))
        if 0 <= r < h and 0 <= c < w:
            wts[r, c] = 1.0
    return wts


def window_for(geom, transform: Affine, shape: tuple[int, int], pad: int = 1) -> Window | None:
    """Pixel window covering the geometry bounds (+pad pixels), clipped to the raster."""
    minx, miny, maxx, maxy = geom.bounds
    win = from_bounds(minx, miny, maxx, maxy, transform=transform)
    c0 = int(np.floor(win.col_off)) - pad
    r0 = int(np.floor(win.row_off)) - pad
    c1 = int(np.ceil(win.col_off + win.width)) + pad
    r1 = int(np.ceil(win.row_off + win.height)) + pad
    c0, r0 = max(c0, 0), max(r0, 0)
    c1, r1 = min(c1, shape[1]), min(r1, shape[0])
    if c1 <= c0 or r1 <= r0:
        return None
    return Window(c0, r0, c1 - c0, r1 - r0)


def weighted_mean(values: np.ndarray, weights: np.ndarray, valid: np.ndarray | None = None):
    """Return (mean, valid_fraction). `values` may be 2-D or 3-D (time first)."""
    if valid is None:
        valid = np.isfinite(values)
    w = weights * valid
    tot = weights.sum()
    if values.ndim == 3:
        num = np.nansum(np.where(valid, values, 0.0) * w, axis=(1, 2))
        den = w.sum(axis=(1, 2))
    else:
        num = np.nansum(np.where(valid, values, 0.0) * w)
        den = w.sum()
    with np.errstate(invalid="ignore", divide="ignore"):
        mean = np.where(den > 0, num / np.where(den > 0, den, 1), np.nan)
        frac = den / tot if tot > 0 else np.zeros_like(den)
    return mean, frac
