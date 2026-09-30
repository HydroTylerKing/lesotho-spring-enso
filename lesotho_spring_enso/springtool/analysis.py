"""ENSO composites, spring-flow response and the susceptibility score.

Per spring the tool builds:

1. A monthly table: rainfall (mm), SPI-1/3/6/12/24, NDVI and its calendar-month
   z-score, ONI and ENSO phase, and the spring's flow (if any).
2. A hydrological-year table (Oct-Sep): rainy-season (Oct-Mar) rainfall and
   % anomaly, peak-season (Jan-Apr) NDVI z-score, flow z-score in the impact
   window, ENSO phase and strength.
3. ENSO composites: mean rainfall anomaly and NDVI z in El Niño, neutral and
   La Niña years, with a Mann-Whitney test of El Niño vs other years.
4. Flow response (where records exist): the rainfall accumulation period and
   lag that best explains flow (Spearman r against SPI-k), and the flow
   anomaly in El Niño impact windows vs other years.
5. A 0-1 susceptibility score = weighted mean of
      R  rain exposure   : El Niño rainy-season deficit, 0 at >=0 %, 1 at <=-30 %
      V  vegetation      : El Niño peak-season NDVI z, 0 at >=0, 1 at <=-1.5
      H  hydrological    : from flow records (El Niño flow anomaly and how fast
                           the spring follows rainfall), or, without records,
                           a catchment-size proxy (small catchments = small storage)
   classed High (>=0.6), Medium (0.4-0.6), Low (<0.4), with a confidence
   level that says how much of the score rests on real flow observations.

"Impact window": flow in calendar year Y is attributed to the rainy season
that ended in Y (hydrological year Y), i.e. Jan-Dec after an El Niño summer.
This covers both springs that fail during the drought and springs that fail
in the following dry season.
"""
from __future__ import annotations

import numpy as np
import pandas as pd
from scipy import stats

from .enso import hydro_year, monthly_phase

SPI_SCALES = (1, 3, 6, 12, 24)
RAINY = (10, 11, 12, 1, 2, 3)
PEAK_NDVI = (1, 2, 3, 4)


# ---------------------------------------------------------------- indices
def spi(precip: pd.Series, k: int) -> pd.Series:
    """Standardised Precipitation Index (gamma fit per calendar month, with zero-rain probability)."""
    acc = precip.rolling(k, min_periods=k).sum()
    out = pd.Series(np.nan, index=precip.index)
    for m in range(1, 13):
        sel = acc.index.month == m
        v = acc[sel].dropna()
        if len(v) < 10:
            continue
        q0 = (v <= 0).mean()
        pos = v[v > 0]
        if len(pos) < 8 or pos.std() == 0:
            continue
        a, _, scale = stats.gamma.fit(pos.values, floc=0)
        x = acc[sel]
        cdf = q0 + (1 - q0) * stats.gamma.cdf(x.clip(lower=0), a, loc=0, scale=scale)
        cdf = np.where(x <= 0, q0 / 2 if q0 > 0 else 1e-6, cdf)
        out[sel] = stats.norm.ppf(np.clip(cdf, 1e-6, 1 - 1e-6))
        out[sel & acc.isna().values] = np.nan
    return out


def month_z(s: pd.Series, min_n: int = 5) -> pd.Series:
    """z-score relative to the same calendar month in all years."""
    g = s.groupby(s.index.month)
    mu, sd, n = g.transform("mean"), g.transform("std"), g.transform("count")
    z = (s - mu) / sd.replace(0, np.nan)
    return z.where(n >= min_n)


# ---------------------------------------------------------------- tables
def build_monthly(precip: pd.Series, ndvi: pd.Series | None, oni: pd.Series,
                  flow: tuple | None) -> pd.DataFrame:
    idx = precip.index
    if ndvi is not None and len(ndvi.dropna()):
        idx = idx.union(ndvi.dropna().index)
    df = pd.DataFrame(index=pd.DatetimeIndex(idx, name="month"))
    df["precip_mm"] = precip.reindex(df.index)
    clim = df["precip_mm"].groupby(df.index.month).transform("mean")
    df["precip_anom_mm"] = df["precip_mm"] - clim
    for k in SPI_SCALES:
        df[f"spi{k}"] = spi(df["precip_mm"], k)
    if ndvi is not None:
        df["ndvi"] = ndvi.reindex(df.index)
        df["ndvi_z"] = month_z(df["ndvi"])
    else:
        df["ndvi"] = np.nan
        df["ndvi_z"] = np.nan
    df["oni"] = oni.reindex(df.index)
    df["enso_phase"] = monthly_phase(oni).reindex(df.index)
    df["hydro_year"] = [hydro_year(t) for t in df.index]
    if flow is not None:
        s, kind = flow
        df["flow"] = s.reindex(df.index)
        df["flow_kind"] = kind
        df["flow_z"] = flow_anomaly(df["flow"], kind)
    else:
        df["flow"] = np.nan
        df["flow_kind"] = ""
        df["flow_z"] = np.nan
    return df


def flow_anomaly(s: pd.Series, kind: str) -> pd.Series:
    """Seasonally adjusted flow anomaly (z). Uses calendar months if the record allows,
    otherwise wet (Oct-Mar) / dry (Apr-Sep) halves, otherwise the whole record."""
    x = np.log1p(s.clip(lower=0)) if kind == "flow" else s.astype(float)
    obs = x.dropna()
    if len(obs) < 4:
        return pd.Series(np.nan, index=s.index)
    per_month = obs.groupby(obs.index.month).count()
    if (per_month >= 3).sum() >= 9:
        return month_z(x, min_n=3)
    half = pd.Series(np.where(np.isin(x.index.month, RAINY), "wet", "dry"), index=x.index)
    cnt = obs.groupby(half[obs.index]).count()
    if len(cnt) == 2 and (cnt >= 3).all():
        mu = x.groupby(half).transform("mean")
        sd = x.groupby(half).transform("std").replace(0, np.nan)
        return (x - mu) / sd
    sd = obs.std()
    return (x - obs.mean()) / sd if sd and sd > 0 else pd.Series(np.nan, index=s.index)


def hydro_table(mon: pd.DataFrame, hy_phase: pd.DataFrame) -> pd.DataFrame:
    rows = []
    for hy, g in mon.groupby("hydro_year"):
        rain = g[g.index.month.isin(RAINY)]["precip_mm"]
        peak = g[g.index.month.isin(PEAK_NDVI)]
        imp = mon[mon.index.year == hy]          # impact window: calendar year hy
        rows.append(dict(
            hydro_year=hy,
            rain_season_mm=rain.sum() if rain.notna().sum() == 6 else np.nan,
            annual_mm=g["precip_mm"].sum() if g["precip_mm"].notna().sum() == 12 else np.nan,
            ndvi_peak=peak["ndvi"].mean(),
            ndvi_peak_z=peak["ndvi_z"].mean(),
            flow_impact_z=imp["flow_z"].mean(),
            flow_impact_n=int(imp["flow_z"].notna().sum()),
            flow_impact_min=imp["flow"].min(),
        ))
    t = pd.DataFrame(rows).set_index("hydro_year")
    mean_rain = t["rain_season_mm"].mean()
    t["rain_anom_pct"] = 100 * (t["rain_season_mm"] / mean_rain - 1)
    t = t.join(hy_phase[["season", "enso_phase", "oni_djf", "oni_peak", "strength", "provisional"]], how="left")
    t["season"] = t["season"].fillna(pd.Series({h: f"{h - 1}/{str(h)[-2:]}" for h in t.index}))
    t["enso_phase"] = t["enso_phase"].fillna("Unknown")
    return t


def _mwu(a, b):
    a, b = pd.Series(a).dropna(), pd.Series(b).dropna()
    if len(a) >= 3 and len(b) >= 3:
        return float(stats.mannwhitneyu(a, b, alternative="two-sided").pvalue)
    return np.nan


def enso_composite(hy: pd.DataFrame) -> dict:
    ok = hy[~hy["provisional"].fillna(False).astype(bool)]
    out = {}
    for col, tag in (("rain_anom_pct", "rain"), ("ndvi_peak_z", "ndvi"), ("flow_impact_z", "flow")):
        for ph, short in (("El Niño", "elnino"), ("Neutral", "neutral"), ("La Niña", "lanina")):
            v = ok.loc[ok["enso_phase"] == ph, col].dropna()
            out[f"{tag}_{short}_mean"] = float(v.mean()) if len(v) else np.nan
            out[f"{tag}_{short}_n"] = int(len(v))
        el = ok.loc[ok["enso_phase"] == "El Niño", col]
        other = ok.loc[ok["enso_phase"] != "El Niño", col]
        out[f"{tag}_elnino_p"] = _mwu(el, other)
    el_rain = ok.loc[ok["enso_phase"] == "El Niño", "rain_anom_pct"].dropna()
    out["rain_elnino_worst"] = float(el_rain.min()) if len(el_rain) else np.nan
    out["rain_elnino_drought_share"] = float((el_rain <= -20).mean()) if len(el_rain) else np.nan
    return out


def flow_response(mon: pd.DataFrame) -> dict:
    """Which rainfall timescale and lag the spring follows (Spearman r with SPI-k)."""
    z = mon["flow_z"]
    best = dict(best_spi=np.nan, best_lag=np.nan, best_r=np.nan, best_p=np.nan, n_flow_months=int(z.notna().sum()))
    if z.notna().sum() < 8:
        return best
    br = -np.inf
    for k in SPI_SCALES:
        for lag in range(0, 7):
            x = mon[f"spi{k}"].shift(lag)
            ok = z.notna() & x.notna()
            if ok.sum() < 8:
                continue
            r, p = stats.spearmanr(z[ok], x[ok])
            if np.isfinite(r) and r > br:
                br = r
                best.update(best_spi=k, best_lag=lag, best_r=float(r), best_p=float(p))
    return best


# ---------------------------------------------------------------- scoring
def _lin(x, zero, one):
    """Linear 0-1 ramp: 0 at `zero`, 1 at `one` (either direction)."""
    if x is None or not np.isfinite(x):
        return np.nan
    return float(np.clip((x - zero) / (one - zero), 0, 1))


def score(comp: dict, resp: dict, area_km2: float, flow_kind: str, weights: dict) -> dict:
    R = _lin(comp.get("rain_elnino_mean"), 0, -30)
    V = _lin(comp.get("ndvi_elnino_mean"), 0, -1.5)

    has_flow = resp.get("n_flow_months", 0) >= 4
    H_parts = []
    if has_flow:
        el, oth = comp.get("flow_elnino_mean"), np.nanmean([comp.get("flow_neutral_mean", np.nan),
                                                             comp.get("flow_lanina_mean", np.nan)])
        if comp.get("flow_elnino_n", 0) >= 1 and np.isfinite(el):
            diff = el - (oth if np.isfinite(oth) else 0.0)
            H_parts.append(_lin(diff, 0, -1.0))
        if np.isfinite(resp.get("best_r", np.nan)) and resp["best_r"] >= 0.3:
            fast = _lin(resp["best_spi"] + resp["best_lag"], 24, 3)   # follows short rainfall windows
            strength_ = _lin(resp["best_r"], 0.3, 0.7)
            H_parts.append(fast * (0.5 + 0.5 * strength_))
        elif np.isfinite(resp.get("best_r", np.nan)):
            H_parts.append(0.0)                                         # no rainfall signal: buffered
    if H_parts:
        H, H_basis = float(np.nanmean(H_parts)), "flow records"
    else:
        # Without records H is a guess: small catchments store less water, but keep the proxy
        # between 0.25 and 0.75 so an unmeasured spring cannot outrank one with observed El Niño failure.
        a = max(area_km2, 1e-3)
        H = 0.25 + 0.5 * _lin(np.log10(a), np.log10(50), np.log10(0.1))
        H_basis = "catchment-size proxy"

    comps = {"R": (R, weights["rain"]), "V": (V, weights["ndvi"]), "H": (H, weights["hydro"])}
    num = sum(v * w for v, w in comps.values() if np.isfinite(v))
    den = sum(w for v, w in comps.values() if np.isfinite(v))
    total = num / den if den else np.nan
    cls = "" if not np.isfinite(total) else "High" if total >= 0.6 else "Medium" if total >= 0.4 else "Low"

    n_el = comp.get("flow_elnino_n", 0)
    if H_basis == "flow records" and n_el >= 2 and resp["n_flow_months"] >= 24:
        conf = "High"
    elif H_basis == "flow records":
        conf = "Medium"
    else:
        conf = "Low"
    return dict(score_rain=R, score_ndvi=V, score_hydro=H, hydro_basis=H_basis,
                susceptibility=total, susceptibility_class=cls, confidence=conf)
