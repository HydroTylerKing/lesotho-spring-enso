"""Reading spring flow records in whatever form they exist.

Accepted: CSV or Excel, long format, one row per observation. Column names
are matched loosely (case-insensitive):

  spring id : spring_id, id, spring, site, site_id, name, code
  date      : date, datetime, obs_date, sample_date, visit_date, measured
  flow      : flow_lps, flow, discharge, q, yield, flow_l_s, lps, l_s, yield_lps
  status    : status, condition, state, flowing

Either a numeric flow column or a status column is needed (both is fine;
numeric flow wins). Status text is mapped to a 0-1 score:
  dry / no flow / not flowing / 0                   -> 0
  low / reduced / weak / trickle / seasonal          -> 0.5
  flowing / normal / good / strong / high / yes      -> 1
"""
from __future__ import annotations

import re
from pathlib import Path

import numpy as np
import pandas as pd

ID_COLS = ["spring_id", "springid", "id", "spring", "site", "site_id", "siteid", "name", "code", "fid"]
DATE_COLS = ["date", "datetime", "obs_date", "sample_date", "visit_date", "measured", "time", "month"]
FLOW_COLS = ["flow_lps", "flow", "discharge", "q", "yield", "flow_l_s", "lps", "l_s", "yield_lps", "flow_ls"]
STATUS_COLS = ["status", "condition", "state", "flowing", "flow_status"]

STATUS_MAP = [
    (r"^(dry|no[\s_-]*flow|not[\s_-]*flowing|none|0|nil|zero|dried)", 0.0),
    (r"^(low|reduced|weak|trickle|seasonal|intermittent|drip|little|0\.5)", 0.5),
    (r"^(flow|normal|good|strong|high|yes|running|perennial|1|2)", 1.0),
]


def _norm(c: str) -> str:
    return re.sub(r"[^a-z0-9]+", "_", str(c).strip().lower()).strip("_")


def _pick(cols: dict, options):
    for o in options:
        if o in cols:
            return cols[o]
    return None


def status_score(v) -> float:
    if pd.isna(v):
        return np.nan
    if isinstance(v, (int, float, np.number)):
        return float(v)
    s = str(v).strip().lower()
    for pat, val in STATUS_MAP:
        if re.match(pat, s):
            return val
    return np.nan


def read_table(path_or_buffer, name: str | None = None) -> pd.DataFrame:
    name = (name or str(path_or_buffer)).lower()
    if name.endswith((".xlsx", ".xls")):
        return pd.read_excel(path_or_buffer)
    return pd.read_csv(path_or_buffer)


def normalise_flows(raw: pd.DataFrame) -> tuple[pd.DataFrame, str]:
    """Return long table [spring_id, date, value, kind] and a note describing what was used."""
    cols = {_norm(c): c for c in raw.columns}
    cid, cdate = _pick(cols, ID_COLS), _pick(cols, DATE_COLS)
    cflow, cstat = _pick(cols, FLOW_COLS), _pick(cols, STATUS_COLS)
    if cid is None or cdate is None or (cflow is None and cstat is None):
        raise ValueError("Flow table needs a spring ID column, a date column and a flow or status column. "
                         f"Found columns: {list(raw.columns)}")
    dates = pd.to_datetime(raw[cdate], errors="coerce", format="ISO8601") \
        if raw[cdate].dtype == object else pd.to_datetime(raw[cdate], errors="coerce")
    if dates.isna().mean() > 0.2:   # e.g. 15/08/2016 as written in field books
        dates = pd.to_datetime(raw[cdate], errors="coerce", dayfirst=True, format="mixed")
    df = pd.DataFrame({"spring_id": raw[cid].astype(str).str.strip(), "date": dates})
    note = []
    flow = pd.to_numeric(raw[cflow], errors="coerce") if cflow is not None else pd.Series(np.nan, index=raw.index)
    stat = raw[cstat].map(status_score) if cstat is not None else pd.Series(np.nan, index=raw.index)
    if stat.max(skipna=True) and stat.max(skipna=True) > 1:  # numeric status codes like 0/1/2
        stat = stat / stat.max()
    df["flow"], df["status"] = flow, stat
    parts = []
    for sid, g in df.groupby("spring_id"):
        n_flow, n_stat = g["flow"].notna().sum(), g["status"].notna().sum()
        if n_flow >= max(n_stat * 0.5, 1):
            v = g["flow"].copy()
            v[v.isna() & (g["status"] == 0)] = 0.0      # 'dry' visits count as zero flow
            parts.append(g.assign(value=v, kind="flow"))
        else:
            parts.append(g.assign(value=g["status"], kind="status"))
    df = pd.concat(parts)[["spring_id", "date", "value", "kind"]]
    kinds = df.groupby("spring_id")["kind"].first().value_counts().to_dict()
    if kinds.get("flow"):
        note.append(f"{kinds['flow']} spring(s) with measured flow ('{cflow}')")
    if kinds.get("status"):
        note.append(f"{kinds['status']} spring(s) with status only ('{cstat}', scored dry 0 / reduced 0.5 / flowing 1)")
    df = df.dropna(subset=["date", "value"])
    return df, "; ".join(note)


def monthly_flows(df: pd.DataFrame) -> dict:
    """{spring_id: (monthly Series, kind)} — mean of all observations within each month."""
    out = {}
    for sid, g in df.groupby("spring_id"):
        s = g.set_index("date")["value"].sort_index()
        m = s.groupby(s.index.to_period("M")).mean()
        m.index = m.index.to_timestamp()
        out[str(sid)] = (m, g["kind"].iloc[0])
    return out


def template() -> pd.DataFrame:
    return pd.DataFrame({
        "spring_id": ["SP001", "SP001", "SP002", "SP002"],
        "date": ["2016-02-15", "2016-08-20", "2016-02-16", "2016-08-21"],
        "flow_lps": [0.42, 0.10, np.nan, np.nan],
        "status": ["flowing", "reduced", "flowing", "dry"],
    })
