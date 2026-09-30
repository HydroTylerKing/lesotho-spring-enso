"""Offline end-to-end test on synthetic data:  python -m tests.test_pipeline"""
import sys
import tempfile
import time
from pathlib import Path

import geopandas as gpd

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
from springtool.config import Settings          # noqa: E402
from springtool.flows import normalise_flows, read_table  # noqa: E402
from springtool.pipeline import run              # noqa: E402
from tests.make_synthetic import make            # noqa: E402


def main(folder=None):
    folder = Path(folder or tempfile.mkdtemp(prefix="springs_synth_"))
    if not (folder / "flows.csv").exists():
        print("building synthetic data in", folder)
        make(folder)
    s = Settings(start="1981-01", end="2025-12", dem_path=str(folder / "dem.tif"), dem_resolution_m=60,
                 chirps_dir=str(folder / "chirps"), ndvi_source="local", ndvi_dir=str(folder / "ndvi"),
                 oni_path=str(folder / "oni.ascii.txt"), cache_dir=str(folder / "cache"))
    pts = gpd.read_file(folder / "springs.gpkg")
    flows, note = normalise_flows(read_table(folder / "flows.csv"))
    print("flows:", note)
    t = time.time()
    res = run(pts, s, flows, progress=lambda m, f: None)
    print(f"run took {time.time() - t:.1f}s")
    cols = ["spring_id", "catchment_km2", "flow_kind", "rain_elnino_mean", "rain_lanina_mean", "ndvi_elnino_mean",
            "flow_elnino_mean", "best_spi", "best_lag", "best_r", "score_rain", "score_ndvi", "score_hydro",
            "hydro_basis", "susceptibility", "susceptibility_class", "confidence"]
    print(res.springs[cols].to_string())
    out = folder / "out"
    out.mkdir(exist_ok=True)
    res.to_gpkg(out / "results.gpkg")
    res.to_excel(out / "results.xlsx")
    res.to_csvs(out)

    sp = res.springs.set_index("spring_id")
    assert sp.loc["HA-01", "rain_elnino_mean"] < -10, "El Niño rainfall deficit not detected"
    assert sp.loc["HA-01", "best_spi"] <= 6, "fast spring should follow short rainfall windows"
    assert sp.loc["HA-02", "best_spi"] >= 12, "slow spring should follow long rainfall windows"
    assert sp.loc["HA-01", "susceptibility"] > sp.loc["HA-02", "susceptibility"], "fast spring should rank higher"
    assert sp.loc["TS-11", "hydro_basis"] == "catchment-size proxy"
    print("ALL CHECKS PASSED — outputs in", out)
    return res


if __name__ == "__main__":
    main(sys.argv[1] if len(sys.argv) > 1 else None)
