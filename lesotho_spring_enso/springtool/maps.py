"""Folium map and Plotly charts shared by the app and batch mode."""
from __future__ import annotations

import folium
import numpy as np
import pandas as pd
import plotly.graph_objects as go
from plotly.subplots import make_subplots

from .config import LESOTHO_CENTER

CLASS_COLOURS = {"High": "#c62828", "Medium": "#ef8f00", "Low": "#2e7d6b", "": "#8a8a8a"}
PHASE_COLOURS = {"El Niño": "rgba(198,40,40,0.13)", "La Niña": "rgba(30,90,180,0.12)"}


def base_map(zoom=8):
    m = folium.Map(location=LESOTHO_CENTER, zoom_start=zoom, tiles=None, control_scale=True)
    folium.TileLayer("OpenStreetMap", name="Streets").add_to(m)
    folium.TileLayer("https://server.arcgisonline.com/ArcGIS/rest/services/World_Imagery/MapServer/tile/{z}/{y}/{x}",
                     attr="Esri World Imagery", name="Imagery").add_to(m)
    folium.TileLayer("https://{s}.tile.opentopomap.org/{z}/{x}/{y}.png", attr="OpenTopoMap", name="Topographic").add_to(m)
    return m


def _fmt(v, nd=2, suffix=""):
    return "–" if v is None or (isinstance(v, float) and not np.isfinite(v)) else f"{v:.{nd}f}{suffix}"


def results_map(res, zoom=8):
    m = base_map(zoom)
    cat = folium.FeatureGroup(name="Catchments").add_to(m)
    for _, r in res.catchments.iterrows():
        col = CLASS_COLOURS.get(r.get("susceptibility_class") or "", "#8a8a8a")
        folium.GeoJson(r.geometry.__geo_interface__, style_function=lambda _, c=col: dict(
            color=c, weight=1.5, fillColor=c, fillOpacity=0.18),
            tooltip=f"{r.spring_id}: {r.area_km2:.2f} km²").add_to(cat)
    pts = folium.FeatureGroup(name="Springs").add_to(m)
    for _, r in res.springs.iterrows():
        if r.geometry is None:
            continue
        cls = r.get("susceptibility_class") or ""
        html = (f"<b>{r.get('name', r.spring_id)}</b> ({r.spring_id})<br>"
                f"Susceptibility: <b>{cls or 'n/a'}</b> ({_fmt(r.get('susceptibility'))}), confidence {r.get('confidence', '')}<br>"
                f"Priority rank: {r.get('priority_rank', '')}<br>"
                f"El Niño rainy-season rain: {_fmt(r.get('rain_elnino_mean'), 0, ' %')}<br>"
                f"La Niña rainy-season rain: {_fmt(r.get('rain_lanina_mean'), 0, ' %')}<br>"
                f"El Niño peak NDVI z: {_fmt(r.get('ndvi_elnino_mean'))}<br>"
                f"Flow data: {r.get('flow_kind', 'none')}; El Niño flow z: {_fmt(r.get('flow_elnino_mean'))}<br>"
                f"Follows SPI-{_fmt(r.get('best_spi'), 0)} (lag {_fmt(r.get('best_lag'), 0)} mo, r={_fmt(r.get('best_r'))})<br>"
                f"Catchment: {_fmt(r.get('catchment_km2'), 3, ' km²')}, snapped {_fmt(r.get('snap_dist_m'), 0, ' m')}")
        folium.CircleMarker([r.geometry.y, r.geometry.x], radius=8, color="white", weight=1.5,
                            fill=True, fill_color=CLASS_COLOURS.get(cls, "#8a8a8a"), fill_opacity=0.95,
                            tooltip=f"{r.spring_id} – {cls or 'n/a'}", popup=folium.Popup(html, max_width=330)).add_to(pts)
    legend = "".join(f'<div><span style="display:inline-block;width:12px;height:12px;border-radius:6px;'
                     f'background:{c};margin-right:6px"></span>{k}</div>' for k, c in CLASS_COLOURS.items() if k)
    m.get_root().html.add_child(folium.Element(
        f'<div style="position:fixed;bottom:24px;left:12px;z-index:999;background:white;padding:8px 10px;'
        f'border-radius:6px;font:12px sans-serif;box-shadow:0 1px 4px rgba(0,0,0,.3)">'
        f'<b>El Niño susceptibility</b>{legend}</div>'))
    folium.LayerControl(collapsed=True).add_to(m)
    if len(res.catchments):
        b = res.catchments.total_bounds
        m.fit_bounds([[b[1], b[0]], [b[3], b[2]]], padding=(20, 20))
    return m


def _enso_bands(fig, hy: pd.DataFrame, rows):
    for h, r in hy.iterrows():
        col = PHASE_COLOURS.get(r.get("enso_phase"))
        if col:
            for row in rows:
                fig.add_vrect(x0=pd.Timestamp(h - 1, 10, 1), x1=pd.Timestamp(h, 4, 1), fillcolor=col,
                              line_width=0, layer="below", row=row, col=1)


def timeseries_figure(mon: pd.DataFrame, hy: pd.DataFrame, title=""):
    has_flow = mon["flow"].notna().any()
    nrows = 3 if has_flow else 2
    fig = make_subplots(rows=nrows, cols=1, shared_xaxes=True, vertical_spacing=0.04,
                        row_heights=[0.4, 0.3, 0.3][:nrows],
                        subplot_titles=["Monthly rainfall (CHIRPS, mm) and SPI-12", "NDVI (catchment mean)",
                                        "Spring flow" if has_flow else ""][:nrows],
                        specs=[[{"secondary_y": True}]] + [[{}]] * (nrows - 1))
    x = mon.index if "month" not in mon.columns else mon["month"]
    fig.add_bar(x=x, y=mon["precip_mm"], name="Rain (mm)", marker_color="#3b6fb6", row=1, col=1)
    fig.add_scatter(x=x, y=mon["spi12"], name="SPI-12", line=dict(color="#222", width=1.3),
                    row=1, col=1, secondary_y=True)
    fig.add_scatter(x=x, y=mon["ndvi"], name="NDVI", line=dict(color="#2e7d32", width=1.4),
                    connectgaps=False, row=2, col=1)
    if has_flow:
        kind = mon["flow_kind"].iloc[0]
        fig.add_scatter(x=x, y=mon["flow"], mode="lines+markers", name="Flow (L/s)" if kind == "flow" else "Status (0 dry – 1 flowing)",
                        marker=dict(size=4), line=dict(color="#6a3d9a", width=1), connectgaps=False, row=3, col=1)
    _enso_bands(fig, hy.set_index("hydro_year") if "hydro_year" in hy.columns else hy, range(1, nrows + 1))
    fig.update_layout(height=620 if has_flow else 470, title=title, margin=dict(l=40, r=20, t=60, b=30),
                      legend=dict(orientation="h", y=-0.08), bargap=0)
    fig.update_yaxes(title_text="SPI", secondary_y=True, row=1, col=1)
    return fig


def composite_figure(hy: pd.DataFrame):
    order = ["El Niño", "Neutral", "La Niña"]
    colours = {"El Niño": "#c62828", "Neutral": "#8a8a8a", "La Niña": "#1e5ab4"}
    hy = hy[hy["enso_phase"].isin(order) & ~hy["provisional"].fillna(False).astype(bool)]
    panels = [("rain_anom_pct", "Rainy-season rain anomaly (%)"), ("ndvi_peak_z", "Peak-season NDVI (z)")]
    if hy["flow_impact_z"].notna().any():
        panels.append(("flow_impact_z", "Spring flow in impact window (z)"))
    fig = make_subplots(rows=1, cols=len(panels), subplot_titles=[p[1] for p in panels])
    for i, (col, _) in enumerate(panels, start=1):
        for ph in order:
            v = hy.loc[hy["enso_phase"] == ph]
            fig.add_box(y=v[col], x=[ph] * len(v), name=ph, marker_color=colours[ph], boxpoints="all",
                        jitter=0.35, pointpos=0, text=v["season"], hovertemplate="%{text}: %{y:.2f}<extra></extra>",
                        showlegend=False, row=1, col=i)
        fig.add_hline(y=0, line_dash="dot", line_color="#666", row=1, col=i)
    fig.update_layout(height=380, margin=dict(l=40, r=20, t=50, b=30))
    return fig
