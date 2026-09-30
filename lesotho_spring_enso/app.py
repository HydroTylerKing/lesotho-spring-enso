"""Lesotho spring El Niño susceptibility tool — Streamlit app.

Run:  streamlit run app.py
"""
from __future__ import annotations

import io
import tempfile
import zipfile
from pathlib import Path

import geopandas as gpd
import numpy as np
import pandas as pd
import streamlit as st
from shapely.geometry import Point
from streamlit_folium import st_folium

from springtool.config import Settings
from springtool.flows import normalise_flows, read_table, template
from springtool.maps import base_map, results_map, timeseries_figure, composite_figure, CLASS_COLOURS
from springtool.pipeline import run, prepare_points

st.set_page_config(page_title="Spring El Niño susceptibility – Lesotho", page_icon="💧", layout="wide")

ss = st.session_state
ss.setdefault("points", pd.DataFrame(columns=["spring_id", "name", "lon", "lat"]))
ss.setdefault("last_click", None)
ss.setdefault("results", None)
ss.setdefault("ver", 0)

# ------------------------------------------------------------------ sidebar
with st.sidebar:
    st.header("Settings")
    c1, c2 = st.columns(2)
    start = c1.text_input("Start (YYYY-MM)", "1981-01")
    end = c2.text_input("End (YYYY-MM)", Settings().end)

    with st.expander("Watershed delineation", expanded=False):
        dem_path = st.text_input("Local DEM (optional)", "", help="Leave blank to use Copernicus GLO-30 (30 m) "
                                 "from Microsoft Planetary Computer.")
        snap = st.slider("Snap radius (m)", 30, 500, 150, 10,
                         help="The spring is moved to the strongest drainage line within this distance.")
        min_acc = st.slider("Min. upstream cells for a drainage line", 1, 500, 20,
                            help="20 cells at 30 m ≈ 1.8 ha.")
        wmax = st.select_slider("Max. DEM window half-width (km)", [20, 40, 60, 80, 120], 60)

    with st.expander("Rainfall (CHIRPS)", expanded=False):
        chirps_version = st.radio("Version", ["v3", "v2"], horizontal=True,
                                  help="CHC stops producing CHIRPS v2 after December 2026.")
        chirps_dir = st.text_input("Local CHIRPS folder (optional)", "",
                                   help="Monthly GeoTIFFs with YYYY.MM in the name. Blank = download.")

    with st.expander("NDVI", expanded=False):
        ndvi_label = st.radio("Source", ["MODIS 250 m (download)", "Local folder", "Skip NDVI"])
        ndvi_source = {"MODIS 250 m (download)": "modis", "Local folder": "local", "Skip NDVI": "none"}[ndvi_label]
        ndvi_dir = st.text_input("Local NDVI folder", "") if ndvi_source == "local" else None
        aqua = st.checkbox("Add MODIS Aqua (twice the download)", False) if ndvi_source == "modis" else False
        min_valid = st.slider("Min. good-quality share per composite", 0.1, 1.0, 0.5, 0.05)

    with st.expander("ENSO", expanded=False):
        oni_path = st.text_input("Local ONI file (optional)", "",
                                 help="Copy of NOAA oni.ascii.txt. Blank = download from NOAA CPC.")

    with st.expander("Susceptibility weights", expanded=False):
        w_rain = st.slider("Rain exposure (R)", 0.0, 1.0, 0.35, 0.05)
        w_ndvi = st.slider("Vegetation response (V)", 0.0, 1.0, 0.15, 0.05)
        w_hydro = st.slider("Spring response (H)", 0.0, 1.0, 0.50, 0.05)

    cache_dir = st.text_input("Cache folder", Settings().cache_dir)

settings = Settings(start=start, end=end, dem_path=dem_path or None, snap_radius_m=float(snap),
                    min_accum_cells=int(min_acc), window_max_km=float(wmax), chirps_version=chirps_version,
                    chirps_dir=chirps_dir or None, ndvi_source=ndvi_source, ndvi_dir=ndvi_dir or None,
                    ndvi_platforms=("terra", "aqua") if aqua else ("terra",), ndvi_min_valid=min_valid,
                    oni_path=oni_path or None, w_rain=w_rain, w_ndvi=w_ndvi, w_hydro=w_hydro, cache_dir=cache_dir)

st.title("💧 Spring El Niño susceptibility — Lesotho")
st.caption("Delineate spring catchments, pull monthly CHIRPS rainfall and NDVI, compare El Niño and La Niña years "
           "with spring flow records, and map which springs to prioritise for pre-positioned assistance.")

tab_in, tab_res, tab_ts, tab_help = st.tabs(["1 · Springs & flow records", "2 · Susceptibility map",
                                             "3 · Time series & ENSO", "About the method"])


# ------------------------------------------------------------------ inputs
def _read_upload(up) -> gpd.GeoDataFrame | None:
    suffix = Path(up.name).suffix.lower()
    tmp = Path(tempfile.mkdtemp()) / up.name
    tmp.write_bytes(up.getvalue())
    if suffix in (".csv", ".xlsx", ".xls"):
        df = read_table(tmp)
        low = {c.lower(): c for c in df.columns}
        lon = next((low[c] for c in ("lon", "longitude", "x", "long") if c in low), None)
        lat = next((low[c] for c in ("lat", "latitude", "y") if c in low), None)
        if not lon or not lat:
            st.error("The table needs lon/lat (or longitude/latitude) columns.")
            return None
        return gpd.GeoDataFrame(df, geometry=gpd.points_from_xy(df[lon], df[lat]), crs="EPSG:4326")
    if suffix == ".zip":
        return gpd.read_file(f"zip://{tmp}")
    if suffix == ".gpkg":
        import pyogrio
        layers = [l[0] for l in pyogrio.list_layers(tmp)]
        layer = st.selectbox("Layer", layers) if len(layers) > 1 else layers[0]
        return gpd.read_file(tmp, layer=layer)
    return gpd.read_file(tmp)


with tab_in:
    left, right = st.columns([3, 2])
    with left:
        st.subheader("Add springs")
        how = st.radio("Input", ["Click on the map", "Upload a GeoPackage / file"], horizontal=True,
                       label_visibility="collapsed")
        if how == "Click on the map":
            st.caption("Click to add a spring. Zoom in first: the point is snapped to the nearest drainage line "
                       f"within {snap} m.")
            m = base_map(8)
            for _, r in ss.points.iterrows():
                import folium
                folium.CircleMarker([r.lat, r.lon], radius=6, color="#1565c0", fill=True, fill_opacity=0.9,
                                    tooltip=f"{r.spring_id} {r['name']}").add_to(m)
            out = st_folium(m, height=520, use_container_width=True, returned_objects=["last_clicked"], key="pick")
            click = out.get("last_clicked") if out else None
            if click and click != ss.last_click:
                ss.last_click = click
                n = len(ss.points) + 1
                new = pd.DataFrame([dict(spring_id=f"SP{n:03d}", name=f"Spring {n}",
                                         lon=round(click["lng"], 6), lat=round(click["lat"], 6))])
                ss.points = pd.concat([ss.points, new], ignore_index=True)
                ss.ver += 1
                st.rerun()
        else:
            up = st.file_uploader("Springs file", type=["gpkg", "geojson", "json", "zip", "csv", "xlsx"],
                                  help="GeoPackage, GeoJSON, zipped shapefile, or a table with lon/lat columns.")
            if up is not None:
                gdf = _read_upload(up)
                if gdf is not None and len(gdf):
                    cols = [c for c in gdf.columns if c != "geometry"]
                    guess = next((c for c in cols if c.lower() in ("spring_id", "id", "code", "site_id", "name")), None)
                    id_col = st.selectbox("Spring ID column", ["(number automatically)"] + cols,
                                          index=(cols.index(guess) + 1) if guess else 0)
                    name_col = st.selectbox("Name column (optional)", ["(same as ID)"] + cols)
                    if st.button(f"Use these {len(gdf)} points", type="primary"):
                        g = prepare_points(gdf, None if id_col.startswith("(") else id_col)
                        names = g[name_col].astype(str) if not name_col.startswith("(") else g["spring_id"]
                        ss.points = pd.DataFrame(dict(spring_id=g["spring_id"], name=names,
                                                      lon=g.geometry.x.round(6), lat=g.geometry.y.round(6)))
                        ss.ver += 1
                        st.rerun()

    with right:
        st.subheader(f"Springs ({len(ss.points)})")
        ss.points = st.data_editor(ss.points, num_rows="dynamic", width="stretch", key=f"pts_editor_{ss.ver}",
                                   column_config={"lon": st.column_config.NumberColumn(format="%.5f"),
                                                  "lat": st.column_config.NumberColumn(format="%.5f")})
        if len(ss.points) and st.button("Clear all springs"):
            ss.points = ss.points.iloc[0:0]
            ss.ver += 1
            st.rerun()

        st.subheader("Spring flow records (optional)")
        st.caption("One row per visit: spring ID, date, and measured flow (L/s) and/or status "
                   "(flowing / reduced / dry). IDs must match the springs above.")
        st.download_button("Download template CSV", template().to_csv(index=False), "spring_flow_template.csv")
        fup = st.file_uploader("Flow records", type=["csv", "xlsx", "xls"], key="flows")
        flows = None
        if fup is not None:
            try:
                flows, note = normalise_flows(read_table(io.BytesIO(fup.getvalue()), fup.name))
                ids = set(ss.points["spring_id"].astype(str).str.lower())
                matched = flows["spring_id"].str.lower().isin(ids)
                st.success(f"{len(flows)} observations, {flows['spring_id'].nunique()} springs — {note}")
                if (~matched).any():
                    st.warning("No spring point for IDs: " + ", ".join(sorted(flows.loc[~matched, "spring_id"].unique())[:20]))
            except Exception as e:
                st.error(str(e))

    st.divider()
    ready = len(ss.points) > 0 and ss.points[["lon", "lat"]].notna().all().all()
    if st.button("▶ Run analysis", type="primary", disabled=not ready, width="stretch"):
        pts = gpd.GeoDataFrame(ss.points.copy(), geometry=[Point(xy) for xy in zip(ss.points.lon, ss.points.lat)],
                               crs="EPSG:4326")
        bar = st.progress(0.0, "Starting")
        try:
            ss.results = run(pts, settings, flows, progress=lambda msg, f: bar.progress(min(max(f, 0.0), 1.0), msg))
            bar.progress(1.0, "Finished — see the Susceptibility map tab")
        except Exception as e:
            bar.empty()
            st.exception(e)


# ------------------------------------------------------------------ results
res = ss.results

with tab_res:
    if res is None:
        st.info("Add springs and run the analysis first.")
    else:
        errs = res.info.get("delineation_errors") or {}
        if errs:
            st.warning("Catchment could not be delineated for: " + ", ".join(f"{k} ({v})" for k, v in errs.items()))
        if res.info.get("ndvi", {}).get("problem"):
            st.warning("NDVI: " + res.info["ndvi"]["problem"])
        c1, c2, c3, c4 = st.columns(4)
        vc = res.springs["susceptibility_class"].value_counts()
        c1.metric("Springs", len(res.springs))
        c2.metric("High", int(vc.get("High", 0)))
        c3.metric("Medium", int(vc.get("Medium", 0)))
        c4.metric("Low", int(vc.get("Low", 0)))
        st_folium(results_map(res), height=560, use_container_width=True, returned_objects=[], key="resmap")

        cols = ["priority_rank", "spring_id", "name", "susceptibility_class", "susceptibility", "confidence",
                "catchment_km2", "rain_elnino_mean", "rain_lanina_mean", "ndvi_elnino_mean", "flow_kind",
                "flow_elnino_mean", "best_spi", "best_lag", "best_r", "hydro_basis", "snap_dist_m", "catchment_truncated"]
        tbl = res.springs[[c for c in cols if c in res.springs]].sort_values("priority_rank")
        st.dataframe(tbl, width="stretch", hide_index=True, column_config={
            "susceptibility": st.column_config.ProgressColumn("score", min_value=0, max_value=1, format="%.2f"),
            "rain_elnino_mean": st.column_config.NumberColumn("El Niño rain %", format="%.0f"),
            "rain_lanina_mean": st.column_config.NumberColumn("La Niña rain %", format="%.0f"),
            "ndvi_elnino_mean": st.column_config.NumberColumn("El Niño NDVI z", format="%.2f"),
            "flow_elnino_mean": st.column_config.NumberColumn("El Niño flow z", format="%.2f"),
        })

        d1, d2, d3 = st.columns(3)
        with tempfile.TemporaryDirectory() as td:
            gp = res.to_gpkg(Path(td) / "spring_enso_susceptibility.gpkg")
            d1.download_button("⬇ GeoPackage (springs, catchments, outlets)", gp.read_bytes(),
                               "spring_enso_susceptibility.gpkg", width="stretch")
        xb = io.BytesIO()
        res.to_excel(xb)
        d2.download_button("⬇ Excel workbook (all tables)", xb.getvalue(), "spring_enso_results.xlsx",
                           width="stretch")
        d3.download_button("⬇ Monthly time series (CSV)", res.monthly.to_csv(index=False), "monthly_timeseries.csv",
                           width="stretch")
        ni = res.info["ndvi"]
        st.caption(f"Rainfall: {res.info['precip']['source']} · NDVI: {ni.get('source')}"
                   + (f" — {ni.get('composites', 0)} composites used, {ni.get('months_with_ndvi', 0)} months with NDVI"
                      if ni.get("source") != "none" else "")
                   + f" · ONI to {res.info.get('oni_last_month')}")

with tab_ts:
    if res is None:
        st.info("Run the analysis to see time series.")
    else:
        labels = {r.spring_id: f"{r.spring_id} – {r['name']} ({r.susceptibility_class or 'n/a'})"
                  for _, r in res.springs.sort_values("priority_rank").iterrows() if r.spring_id in set(res.monthly.spring_id)}
        sid = st.selectbox("Spring", list(labels), format_func=labels.get)
        mon = res.monthly[res.monthly.spring_id == sid].set_index("month")
        if res.info.get("ndvi", {}).get("source") != "none" and mon["ndvi"].notna().sum() == 0:
            st.warning("No NDVI for this spring. " + (res.info["ndvi"].get("problem") or
                       "The catchment may be too cloudy or snowy in every composite; lower the quality threshold."))
        hy = res.hydro[res.hydro.spring_id == sid]
        yrs = st.slider("Years", int(mon.index.year.min()), int(mon.index.year.max()),
                        (max(int(mon.index.year.min()), 2000), int(mon.index.year.max())))
        view = mon[(mon.index.year >= yrs[0]) & (mon.index.year <= yrs[1])]
        st.plotly_chart(timeseries_figure(view, hy[(hy.hydro_year >= yrs[0]) & (hy.hydro_year <= yrs[1] + 1)],
                                          "Shaded: Oct–Mar rainy seasons of El Niño (red) and La Niña (blue) years"),
                        width="stretch")
        st.subheader("El Niño vs neutral vs La Niña years")
        st.plotly_chart(composite_figure(hy), width="stretch")
        r = res.springs.set_index("spring_id").loc[sid]
        st.markdown(
            f"**Score {r.susceptibility:.2f} ({r.susceptibility_class}, confidence {r.confidence})** = "
            f"rain exposure R {r.score_rain:.2f} · vegetation V "
            f"{'–' if pd.isna(r.score_ndvi) else f'{r.score_ndvi:.2f}'} · spring response H {r.score_hydro:.2f} "
            f"(from {r.hydro_basis}).")
        show = ["season", "enso_phase", "strength", "oni_djf", "rain_season_mm", "rain_anom_pct", "ndvi_peak",
                "ndvi_peak_z", "flow_impact_z", "flow_impact_n", "flow_impact_min"]
        st.dataframe(hy.set_index("hydro_year")[show].sort_index(ascending=False), width="stretch")

with tab_help:
    st.markdown(Path(__file__).with_name("METHOD.md").read_text(encoding="utf-8"))
