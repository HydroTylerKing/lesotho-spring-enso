"""Generate a complete synthetic dataset so the tool can be tested offline.

Creates, in the given folder:
  dem.tif            synthetic terrain, UTM 35S, 60 m
  chirps/            monthly rainfall GeoTIFFs 1981-2025, 0.05°, with an El Niño signal
  ndvi/              monthly NDVI GeoTIFFs 2000-2025, UTM 250 m, lagging rainfall
  oni.ascii.txt      NOAA-format ONI table (synthetic values, NOT real)
  springs.gpkg       4 test springs
  flows.csv          flow records: fast spring, slow spring, status-only spring, one with none

Nothing here is real data. It exists only to check that the code runs and that
a spring built to respond to El Niño is ranked above one built not to.
"""
from __future__ import annotations

import sys
from pathlib import Path

import geopandas as gpd
import numpy as np
import pandas as pd
import rasterio
from affine import Affine
from shapely.geometry import Point

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
from tests.synthetic_dem import make_dem  # noqa: E402

SEAS = ["DJF", "JFM", "FMA", "MAM", "AMJ", "MJJ", "JJA", "JAS", "ASO", "SON", "OND", "NDJ"]
EL_NINO = {1983: 2.1, 1987: 1.2, 1992: 1.6, 1995: 1.0, 1998: 2.4, 2003: 1.1, 2007: 0.8,
           2010: 1.5, 2016: 2.6, 2019: 0.8, 2024: 2.0}
LA_NINA = {1985: -1.0, 1989: -1.7, 1996: -0.9, 1999: -1.5, 2000: -1.7, 2008: -1.6, 2011: -1.6,
           2018: -0.9, 2021: -1.2, 2022: -1.0, 2023: -1.0}


def synthetic_oni(years=range(1980, 2027)):
    rng = np.random.default_rng(3)
    months = pd.date_range(f"{years.start}-01-01", f"{years.stop - 1}-06-01", freq="MS")
    v = rng.normal(0, 0.15, len(months))
    for i, t in enumerate(months):
        for tab in (EL_NINO, LA_NINA):
            for hy, peak in tab.items():
                # event builds from Jun (hy-1) to peak in Dec/Jan and decays by May hy
                d = (t - pd.Timestamp(hy - 1, 12, 15)).days / 30.4
                if -7 <= d <= 5:
                    v[i] += peak * np.cos(np.pi * d / 14) ** 2
    return pd.Series(v, index=months)


def write_oni(path, oni: pd.Series):
    lines = ["SEAS  YR   TOTAL   ANOM"]
    for t, a in oni.items():
        lines.append(f"{SEAS[t.month - 1]:>4} {t.year:4d} {26.5 + a:7.2f} {a:7.2f}")
    Path(path).write_text("\n".join(lines) + "\n")


def make(folder):
    folder = Path(folder)
    folder.mkdir(parents=True, exist_ok=True)
    rng = np.random.default_rng(7)
    make_dem(folder / "dem.tif", half_km=70, res=60)

    oni = synthetic_oni()
    write_oni(folder / "oni.ascii.txt", oni)

    # ---- rainfall: seasonal climatology, west-east gradient, El Niño deficit, noise
    clim = np.array([120, 100, 85, 45, 20, 8, 8, 12, 30, 70, 90, 110], float)  # mm, Jan..Dec
    x0, y1, res = 26.9, -28.4, 0.05
    nx, ny = 54, 48
    t = Affine(res, 0, x0, 0, -res, y1)
    lon = x0 + (np.arange(nx) + 0.5) * res
    grad = 0.8 + 0.4 * (lon - lon.min()) / (lon.max() - lon.min())
    months = pd.date_range("1981-01-01", "2025-12-01", freq="MS")
    (folder / "chirps").mkdir(exist_ok=True)
    rain_series = []
    for m in months:
        o = oni.get(m, 0.0)
        wet = m.month in (10, 11, 12, 1, 2, 3)
        factor = np.clip(1 - (0.32 * o if wet else 0.05 * o) + rng.normal(0, 0.25), 0.05, None)
        field = clim[m.month - 1] * factor * grad[None, :].repeat(ny, 0)
        field *= np.exp(rng.normal(0, 0.08, field.shape))
        rain_series.append(field.mean())
        with rasterio.open(folder / "chirps" / f"chirps-v3.0.{m.year}.{m.month:02d}.tif", "w", driver="GTiff",
                           height=ny, width=nx, count=1, dtype="float32", crs="EPSG:4326",
                           transform=t, nodata=-9999) as ds:
            ds.write(field.astype("float32"), 1)
    rain = pd.Series(rain_series, index=months)

    # ---- NDVI: follows 3-month rainfall anomaly, UTM 250 m, one file per month
    from pyproj import Transformer
    cx, cy = Transformer.from_crs("EPSG:4326", "EPSG:32735", always_xy=True).transform(28.2, -29.5)
    nres, n = 250, 600
    tn = Affine(nres, 0, cx - n * nres / 2, 0, -nres, cy + n * nres / 2)
    r3 = rain.rolling(3).sum()
    r3z = (r3 - r3.groupby(r3.index.month).transform("mean")) / r3.groupby(r3.index.month).transform("std")
    seas = 0.38 + 0.14 * np.cos(2 * np.pi * (np.arange(12) - 1.5) / 12)   # peaks Feb
    (folder / "ndvi").mkdir(exist_ok=True)
    base = 0.9 + 0.2 * rng.random((n, n))
    for m in pd.date_range("2000-03-01", "2025-12-01", freq="MS"):
        val = seas[m.month - 1] + 0.05 * np.nan_to_num(r3z.get(m, 0.0))
        arr = (np.clip(val * base + rng.normal(0, 0.02, (n, n)), -0.1, 1) * 10000).astype("int16")
        with rasterio.open(folder / "ndvi" / f"ndvi_{m.year}_{m.month:02d}.tif", "w", driver="GTiff",
                           height=n, width=n, count=1, dtype="int16", crs="EPSG:32735",
                           transform=tn, nodata=-3000) as ds:
            ds.write(arr, 1)

    # ---- springs
    pts = gpd.GeoDataFrame({"spring_id": ["HA-01", "HA-02", "MO-07", "TS-11"],
                            "name": ["Ha Fast (responsive)", "Ha Slow (buffered)", "Mohale status-only", "No records"],
                            "village": ["Ha Synthetic A", "Ha Synthetic B", "Synthetic C", "Synthetic D"]},
                           geometry=[Point(28.20, -29.50), Point(28.40, -29.60), Point(28.05, -29.45), Point(28.30, -29.35)],
                           crs="EPSG:4326")
    pts.to_file(folder / "springs.gpkg", layer="springs", driver="GPKG")

    # ---- flows
    spi_like = lambda k: (rain.rolling(k).sum() - rain.rolling(k).sum().groupby(rain.index.month).transform("mean")) \
        / rain.rolling(k).sum().groupby(rain.index.month).transform("std")
    fast = np.exp(0.6 * spi_like(3).shift(1)) * (0.5 + 0.2 * np.isin(months.month, [1, 2, 3, 4]))
    slow = np.exp(0.25 * spi_like(24)) * 0.8
    rows = []
    for d in pd.date_range("2005-01-01", "2025-12-01", freq="MS"):
        if rng.random() < 0.8:
            rows.append(("HA-01", d + pd.Timedelta(days=int(rng.integers(0, 27))),
                         round(float(fast[d] * np.exp(rng.normal(0, 0.1))), 3), ""))
        if rng.random() < 0.5:
            rows.append(("HA-02", d + pd.Timedelta(days=int(rng.integers(0, 27))),
                         round(float(slow[d] * np.exp(rng.normal(0, 0.1))), 3), ""))
    for d in pd.date_range("2008-02-01", "2025-11-01", freq="3MS"):
        z = spi_like(6)[d]
        st = "dry" if z < -0.8 else "reduced" if z < -0.1 else "flowing"
        rows.append(("MO-07", d + pd.Timedelta(days=5), np.nan, st))
    pd.DataFrame(rows, columns=["spring_id", "date", "flow_lps", "status"]).to_csv(folder / "flows.csv", index=False)
    return folder


if __name__ == "__main__":
    make(sys.argv[1] if len(sys.argv) > 1 else "synthetic_data")
