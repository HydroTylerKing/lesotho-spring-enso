"""Settings shared by the app, the batch runner and the processing modules."""
from __future__ import annotations

import datetime as _dt
from dataclasses import dataclass, field, asdict
from pathlib import Path

# Lesotho bounding box (lon_min, lat_min, lon_max, lat_max) with a small margin
LESOTHO_BBOX = (26.95, -30.72, 29.50, -28.54)
LESOTHO_CENTER = (-29.6, 28.25)  # lat, lon
UTM_CRS = "EPSG:32735"          # WGS 84 / UTM zone 35S covers all of Lesotho
WGS84 = "EPSG:4326"

PC_STAC = "https://planetarycomputer.microsoft.com/api/stac/v1"


def _last_complete_month() -> str:
    today = _dt.date.today().replace(day=1)
    last = today - _dt.timedelta(days=1)
    return f"{last.year}-{last.month:02d}"


@dataclass
class Settings:
    # ---- period -------------------------------------------------------
    start: str = "1981-01"                 # CHIRPS starts Jan 1981
    end: str = field(default_factory=_last_complete_month)

    # ---- watershed delineation ---------------------------------------
    dem_path: str | None = None            # local DEM (any CRS). None = Copernicus GLO-30 from Planetary Computer
    dem_resolution_m: float = 30.0
    snap_radius_m: float = 150.0           # search radius for the drainage line (GPS error + DEM error)
    min_accum_cells: int = 20              # smallest upstream area (cells) to count as a drainage line
    window_start_km: float = 5.0           # half-width of first DEM window
    window_max_km: float = 60.0            # half-width limit before a catchment is flagged as truncated

    # ---- precipitation ------------------------------------------------
    chirps_version: str = "v3"             # "v3" (current) or "v2" (production ends Dec 2026)
    chirps_dir: str | None = None          # folder of local monthly CHIRPS GeoTIFFs; None = download from CHC

    # ---- NDVI ----------------------------------------------------------
    ndvi_source: str = "modis"             # "modis" (Planetary Computer MOD13Q1/MYD13Q1), "local", or "none"
    ndvi_dir: str | None = None            # folder of local monthly NDVI GeoTIFFs (ndvi_source="local")
    ndvi_platforms: tuple = ("terra",)     # ("terra",) or ("terra", "aqua"); aqua doubles download volume
    ndvi_min_valid: float = 0.5            # drop a composite if < this share of the catchment is good quality
    ndvi_start: str = "2000-03"

    # ---- ENSO ---------------------------------------------------------
    oni_path: str | None = None            # local copy of NOAA oni.ascii.txt; None = download

    # ---- susceptibility weights ---------------------------------------
    w_rain: float = 0.35                   # El Nino rainy-season rainfall deficit in the catchment
    w_ndvi: float = 0.15                   # El Nino vegetation (NDVI) response in the catchment
    w_hydro: float = 0.50                  # spring's own flow response (or catchment-size proxy)

    # ---- housekeeping ---------------------------------------------------
    cache_dir: str = str(Path.home() / ".lesotho_spring_enso_cache")
    workers: int = 8

    def cache(self, *parts: str) -> Path:
        p = Path(self.cache_dir).joinpath(*parts)
        p.mkdir(parents=True, exist_ok=True)
        return p

    def to_dict(self) -> dict:
        return asdict(self)
