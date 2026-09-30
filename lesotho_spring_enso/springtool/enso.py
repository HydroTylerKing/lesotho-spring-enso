"""ENSO phases from NOAA CPC's Oceanic Niño Index (ONI, ERSSTv5).

NOAA definition: El Niño (La Niña) conditions when the 3-month running mean
SST anomaly in Niño 3.4 is >= +0.5 °C (<= -0.5 °C) for at least 5 consecutive
overlapping seasons.

Lesotho's rain falls mostly October-March, so ENSO is assigned to *hydrological
years* (October-September, labelled e.g. "2015/16") using the ONI of the
DJF season (centred on January), which is when ENSO events peak and when
the southern African rainfall teleconnection is strongest.
"""
from __future__ import annotations

from pathlib import Path

import numpy as np
import pandas as pd

from .config import Settings

ONI_URL = "https://www.cpc.ncep.noaa.gov/data/indices/oni.ascii.txt"
SEASON_CENTRE = {"DJF": 1, "JFM": 2, "FMA": 3, "MAM": 4, "AMJ": 5, "MJJ": 6,
                 "JJA": 7, "JAS": 8, "ASO": 9, "SON": 10, "OND": 11, "NDJ": 12}


def load_oni(settings: Settings) -> pd.Series:
    """Monthly ONI series (index = first of centre month)."""
    if settings.oni_path:
        text = Path(settings.oni_path).read_text()
    else:
        cache = settings.cache("enso") / "oni.ascii.txt"
        try:
            import urllib.request
            with urllib.request.urlopen(ONI_URL, timeout=30) as r:
                text = r.read().decode()
            cache.write_text(text)
        except Exception:
            if not cache.exists():
                raise RuntimeError("Could not download the ONI table from NOAA and no cached copy exists. "
                                   f"Download {ONI_URL} and set 'ONI file' in the settings.")
            text = cache.read_text()
    rows = []
    for line in text.splitlines():
        parts = line.split()
        if len(parts) >= 4 and parts[0] in SEASON_CENTRE:
            try:
                rows.append((pd.Timestamp(int(parts[1]), SEASON_CENTRE[parts[0]], 1), float(parts[3])))
            except ValueError:
                continue
    s = pd.Series(dict(rows)).sort_index()
    s.index.name = "month"
    return s.rename("oni")


def _runs(flag: np.ndarray, min_len: int = 5) -> np.ndarray:
    out = np.zeros_like(flag, dtype=bool)
    i = 0
    while i < len(flag):
        if flag[i]:
            j = i
            while j < len(flag) and flag[j]:
                j += 1
            if j - i >= min_len:
                out[i:j] = True
            i = j
        else:
            i += 1
    return out


def monthly_phase(oni: pd.Series) -> pd.Series:
    v = oni.values
    el = _runs(v >= 0.5)
    la = _runs(v <= -0.5)
    ph = np.where(el, "El Niño", np.where(la, "La Niña", "Neutral"))
    return pd.Series(ph, index=oni.index, name="enso_phase")


def hydro_year(ts) -> int:
    """Hydrological year Oct-Sep, identified by its end year (Oct 2015-Sep 2016 -> 2016)."""
    ts = pd.Timestamp(ts)
    return ts.year + 1 if ts.month >= 10 else ts.year


def hy_label(hy: int) -> str:
    return f"{hy - 1}/{str(hy)[-2:]}"


def strength(peak: float) -> str:
    a = abs(peak)
    if a >= 2.0:
        return "very strong"
    if a >= 1.5:
        return "strong"
    if a >= 1.0:
        return "moderate"
    if a >= 0.5:
        return "weak"
    return ""


def hydro_year_phase(oni: pd.Series) -> pd.DataFrame:
    """One row per hydrological year: ENSO phase (from DJF), peak ONI Oct-Mar, strength."""
    ph = monthly_phase(oni)
    years = sorted({hydro_year(t) for t in oni.index})
    rows = []
    for hy in years:
        jan = pd.Timestamp(hy, 1, 1)
        season = oni[(oni.index >= pd.Timestamp(hy - 1, 10, 1)) & (oni.index <= pd.Timestamp(hy, 3, 1))]
        if jan not in ph.index:
            continue
        phase = ph[jan]
        # a run needs 5 seasons: the last few months of the record are provisional
        provisional = jan > oni.index.max() - pd.DateOffset(months=4)
        peak = season.max() if phase == "El Niño" else season.min() if phase == "La Niña" else \
            season.iloc[int(np.argmax(np.abs(season.values)))]
        rows.append(dict(hydro_year=hy, season=hy_label(hy), enso_phase=phase,
                         oni_djf=float(oni[jan]), oni_peak=float(peak),
                         strength=strength(peak) if phase != "Neutral" else "",
                         provisional=bool(provisional)))
    return pd.DataFrame(rows).set_index("hydro_year")
