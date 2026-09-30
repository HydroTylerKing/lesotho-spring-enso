"""Make a synthetic 'Lesotho-like' DEM in UTM 35S for offline testing."""
import numpy as np, rasterio
from affine import Affine
from pyproj import Transformer

def make_dem(path, lon=28.2, lat=-29.5, half_km=70, res=60, seed=1):
    x, y = Transformer.from_crs("EPSG:4326", "EPSG:32735", always_xy=True).transform(lon, lat)
    n = int(2 * half_km * 1000 / res)
    t = Affine(res, 0, x - half_km * 1000, 0, -res, y + half_km * 1000)
    yy, xx = np.mgrid[0:n, 0:n] * res / 1000.0
    rng = np.random.default_rng(seed)
    # plateau sloping to the west, with sinuous valleys and ridges
    z = 3000 - 12 * xx + 8 * np.abs(yy - half_km - 10 * np.sin(xx / 9)) \
        + 60 * np.sin(yy / 3.1) * np.cos(xx / 4.3) + 25 * np.sin(xx / 1.7 + yy / 2.3)
    from scipy.ndimage import gaussian_filter
    z += gaussian_filter(rng.normal(0, 40, (n, n)), 3)
    with rasterio.open(path, "w", driver="GTiff", height=n, width=n, count=1, dtype="float32",
                       crs="EPSG:32735", transform=t, nodata=-9999) as ds:
        ds.write(z.astype("float32"), 1)
    return path

if __name__ == "__main__":
    make_dem("/tmp/claude-0/-home-claude/f17e18f2-ad14-534d-96e2-1b341e030a45/scratchpad/dem.tif")
