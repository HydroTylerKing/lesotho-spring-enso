#!/usr/bin/env python
"""Test each online data source the tool uses, one at a time, and print what fails.

    python check_connection.py

Takes under a minute. Nothing is cached or changed.
"""
import traceback

LON, LAT = 27.95, -29.35   # near Maseru


def step(name, fn):
    print(f"\n── {name}")
    try:
        msg = fn()
        print(f"   OK  {msg}")
        return True
    except Exception as e:
        print(f"   FAILED  {type(e).__name__}: {e}")
        traceback.print_exc(limit=3)
        return False


def modis():
    import planetary_computer
    import pystac_client
    import rasterio
    from pyproj import Transformer
    from springtool.ndvi import COLLECTION, NDVI_ASSET, REL_ASSET, GDAL_ENV
    from springtool.config import PC_STAC
    cat = pystac_client.Client.open(PC_STAC)
    items = list(cat.search(collections=[COLLECTION], bbox=[LON - .01, LAT - .01, LON + .01, LAT + .01],
                            datetime="2024-01-01/2024-03-31").items())
    terra = [i for i in items if i.properties.get("platform") == "terra"]
    print(f"   catalogue: {len(items)} composites Jan–Mar 2024 ({len(terra)} Terra)")
    it = (terra or items)[0]
    href = planetary_computer.sign(it.assets[NDVI_ASSET].href)
    print("   signed URL ok:", href.split("?")[0][-60:])
    with rasterio.Env(**GDAL_ENV):
        with rasterio.open(href) as ds:
            x, y = Transformer.from_crs("EPSG:4326", ds.crs, always_xy=True).transform(LON, LAT)
            r, c = ds.index(x, y)
            v = ds.read(1, window=((r, r + 1), (c, c + 1)))[0, 0]
        with rasterio.open(planetary_computer.sign(it.assets[REL_ASSET].href)) as dr:
            q = dr.read(1, window=((r, r + 1), (c, c + 1)))[0, 0]
    return f"{it.id}: NDVI at Maseru = {v * 1e-4:.3f}, reliability = {q}"


def probe(url):
    """Ask for the first 16 bytes and report what the server sends back."""
    import urllib.request
    req = urllib.request.Request(url, headers={"Range": "bytes=0-15"})
    try:
        with urllib.request.urlopen(req, timeout=30) as r:
            head = r.read(16)
            kind = "GeoTIFF" if head[:2] in (b"II", b"MM") else "gzip" if head[:2] == b"\x1f\x8b" else "NOT a GeoTIFF"
            print(f"   {url}\n     HTTP {r.status} | type {r.headers.get('Content-Type')} | "
                  f"range {'yes' if r.status == 206 else 'NO'} | starts {head[:8]!r} -> {kind}")
    except Exception as e:
        print(f"   {url}\n     {type(e).__name__}: {e}")


def chirps():
    import rasterio
    from rasterio.windows import from_bounds
    from springtool.chirps import URLS, GDAL_ENV, _download_and_read
    print("   What the CHC server returns:")
    for u in (URLS["v3"].format(y=2024, m=1), URLS["v2"].format(y=2024, m=1)):
        probe(u)
    url = URLS["v3"].format(y=2024, m=1)
    try:
        with rasterio.Env(**GDAL_ENV), rasterio.open("/vsicurl/" + url) as ds:
            w = from_bounds(LON, LAT - 0.05, LON + 0.05, LAT, ds.transform).round_offsets().round_lengths()
            v = ds.read(1, window=w)[0, 0]
        return f"CHIRPS v3 Jan 2024 rainfall at Maseru = {v:.0f} mm (fast windowed read works)"
    except Exception as e:
        print(f"   fast windowed read failed ({type(e).__name__}: {e}); trying full-file download…")
    import tempfile
    region = (LON, LAT - 0.05, LON + 0.05, LAT)
    arr = _download_and_read(url, region, (1, 1), tempfile.gettempdir())
    return (f"CHIRPS v3 Jan 2024 rainfall at Maseru = {arr[0, 0]:.0f} mm "
            "(full-file download works, but it is slow: consider a local CHIRPS folder)")


def dem():
    import planetary_computer
    import pystac_client
    import rasterio
    from springtool.config import PC_STAC
    cat = pystac_client.Client.open(PC_STAC, modifier=planetary_computer.sign_inplace)
    it = next(cat.search(collections=["cop-dem-glo-30"], bbox=[LON, LAT, LON + .001, LAT + .001]).items())
    with rasterio.open(it.assets["data"].href) as ds:
        r, c = ds.index(LON, LAT)
        v = ds.read(1, window=((r, r + 1), (c, c + 1)))[0, 0]
    return f"Copernicus DEM elevation at Maseru = {v:.0f} m"


def oni():
    import urllib.request
    from springtool.enso import ONI_URL
    with urllib.request.urlopen(ONI_URL, timeout=30) as r:
        last = r.read().decode().strip().splitlines()[-1]
    return f"latest ONI row: {last}"


if __name__ == "__main__":
    import rasterio
    print("rasterio", rasterio.__version__, "| GDAL", rasterio.__gdal_version__)
    ok = [step("MODIS NDVI (Planetary Computer)", modis), step("CHIRPS rainfall (CHC)", chirps),
          step("Copernicus DEM (Planetary Computer)", dem), step("ONI (NOAA CPC)", oni)]
    print("\nAll sources reachable." if all(ok) else
          "\nSome sources failed — send the FAILED lines above to whoever supports the tool.")
