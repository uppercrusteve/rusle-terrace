from __future__ import annotations

import csv
import json
import struct
from pathlib import Path

import numpy as np
import pytest
import rasterio

from conftest import BASE, PARAMS_TEXT, south_plane, write_dem
from rusle_terrace import cli
from rusle_terrace.factors import params_from_mapping
from rusle_terrace.runoff import DemError, analyze, read_dem, write_erosion


def test_missing_crs_raises(tmp_path: Path) -> None:
    dem = write_dem(tmp_path / "d.tif", south_plane(4, 4), crs=None)
    with pytest.raises(DemError, match="no CRS"):
        read_dem(dem)


def test_geographic_crs_raises(tmp_path: Path) -> None:
    dem = write_dem(tmp_path / "d.tif", south_plane(4, 4), crs="EPSG:4326")
    with pytest.raises(DemError, match="projected"):
        read_dem(dem)


def test_non_metric_projected_crs_raises(tmp_path: Path) -> None:
    # EPSG:2225 uses US survey feet.
    dem = write_dem(tmp_path / "d.tif", south_plane(4, 4), crs="EPSG:2225")
    with pytest.raises(DemError, match="metres"):
        read_dem(dem)


def test_missing_dem_raises(tmp_path: Path) -> None:
    with pytest.raises(DemError, match="not found"):
        read_dem(tmp_path / "missing.tif")


def test_multiband_dem_rejected(tmp_path: Path) -> None:
    path = tmp_path / "two.tif"
    with rasterio.open(
        path,
        "w",
        driver="GTiff",
        height=4,
        width=4,
        count=2,
        dtype="float32",
        crs="EPSG:32647",
        transform=rasterio.Affine(30, 0, 500000, 0, -30, 2800000),
    ) as dst:
        dst.write(np.ones((2, 4, 4), dtype="float32"))
    with pytest.raises(DemError, match="one band"):
        read_dem(path)


def test_nodata_is_masked_in_outputs(tmp_path: Path) -> None:
    z = south_plane(8, 8)
    z[5:7, 5:7] = -9999.0
    dem = write_dem(tmp_path / "d.tif", z, nodata=-9999.0)
    params = params_from_mapping(BASE)
    result = analyze(dem, params)
    assert np.isnan(result.erosion[5, 5])
    assert result.stats[1]["cell_count"] < 64
    out = write_erosion(tmp_path / "e.tif", result)
    with rasterio.open(out) as src:
        masked = src.read(1, masked=True)
        assert masked[5, 5] is np.ma.masked
        assert masked[3, 3] is not np.ma.masked


def test_demo_end_to_end(tmp_path: Path, capsys: pytest.CaptureFixture[str]) -> None:
    out = tmp_path / "demo"
    assert cli.main(["demo", "-o", str(out)]) == 0
    for name in (
        "dem.tif",
        "params.yaml",
        "dem_truth.json",
        "erosion.tif",
        "erosion-map.png",
        "slope-class-stats.csv",
        "summary.json",
        "report.md",
    ):
        assert (out / name).is_file(), name
    payload = json.loads((out / "summary.json").read_text(encoding="utf-8"))
    assert payload["data"]["synthetic"] is True
    assert "合成" in payload["data"]["label"]
    assert payload["results_t_ha_yr"]["mean"] > 0
    assert "SYNTHETIC" in capsys.readouterr().out

    png = (out / "erosion-map.png").read_bytes()
    assert png[:8] == b"\x89PNG\r\n\x1a\n"
    width, height = struct.unpack(">II", png[16:24])
    assert width >= 1200 and height >= 1200


def test_demo_truth_file_records_analytic_slope(tmp_path: Path) -> None:
    out = tmp_path / "demo"
    cli.main(["demo", "-o", str(out)])
    truth = json.loads((out / "dem_truth.json").read_text(encoding="utf-8"))
    assert truth["synthetic"] is True
    assert truth["analytic_slope_deg"]["max"] > 15.0
    assert truth["analytic_slope_deg"]["min"] < 5.0
    assert truth["design"]["terraces"] == 3


def test_stats_csv_has_five_classes(tmp_path: Path) -> None:
    out = tmp_path / "demo"
    cli.main(["demo", "-o", str(out)])
    with (out / "slope-class-stats.csv").open(encoding="utf-8-sig", newline="") as handle:
        rows = list(csv.DictReader(handle))
    assert len(rows) == 5
    assert float(rows[0]["class_lo_deg"]) == 0.0
    assert float(rows[0]["mean_erosion_t_ha_yr"]) > 0
    assert float(rows[3]["mean_erosion_t_ha_yr"]) > float(rows[0]["mean_erosion_t_ha_yr"])
    assert any(row["severity_by_mean"] == "无数据 (no data)" for row in rows)


def test_report_markdown_marks_synthetic(tmp_path: Path) -> None:
    out = tmp_path / "demo"
    cli.main(["demo", "-o", str(out)])
    text = (out / "report.md").read_text(encoding="utf-8")
    assert "合成演示" in text
    assert "docs/methods.md" in text
    assert "规划级" in text


def test_cli_help_exits_zero(capsys: pytest.CaptureFixture[str]) -> None:
    with pytest.raises(SystemExit) as excinfo:
        cli.main(["--help"])
    assert excinfo.value.code == 0
    assert "run" in capsys.readouterr().out


def test_cli_version(capsys: pytest.CaptureFixture[str]) -> None:
    with pytest.raises(SystemExit) as excinfo:
        cli.main(["--version"])
    assert excinfo.value.code == 0
    assert "0.1.0" in capsys.readouterr().out


def test_cli_run_reports_missing_crs(tmp_path: Path, capsys: pytest.CaptureFixture[str]) -> None:
    dem = write_dem(tmp_path / "d.tif", south_plane(4, 4), crs=None)
    params = tmp_path / "p.yaml"
    params.write_text(PARAMS_TEXT, encoding="utf-8")
    code = cli.main(["run", str(dem), str(params), "-o", str(tmp_path / "out")])
    assert code == 2
    assert "no CRS" in capsys.readouterr().err


def test_cli_run_rejects_bad_p_override(
    tmp_path: Path, dem_file: Path, params_file: Path, capsys: pytest.CaptureFixture[str]
) -> None:
    code = cli.main(
        ["run", str(dem_file), str(params_file), "-o", str(tmp_path / "out"), "--p-factor", "2"]
    )
    assert code == 2
    assert "0-1" in capsys.readouterr().err


def test_cli_run_success(tmp_path: Path, dem_file: Path, params_file: Path) -> None:
    out = tmp_path / "out"
    assert cli.main(["run", str(dem_file), str(params_file), "-o", str(out)]) == 0
    summary = json.loads((out / "summary.json").read_text(encoding="utf-8"))
    assert summary["data"]["synthetic"] is False
    assert summary["results_t_ha_yr"]["valid_cells"] == 64
    assert (out / "erosion-map.png").stat().st_size > 20_000


def test_cli_factors_emits_json(capsys: pytest.CaptureFixture[str]) -> None:
    code = cli.main(
        ["factors", "--precip", "620", "--sand", "45", "--silt", "35",
         "--clay", "20", "--oc", "1.2"]
    )
    assert code == 0
    payload = json.loads(capsys.readouterr().out)
    assert 0.01 <= payload["K_t_ha_h_ha_MJ_mm"] <= 0.7
    assert payload["P_by_slope_deg"]["20.0"] == 1.0


def test_cli_factors_bad_texture(capsys: pytest.CaptureFixture[str]) -> None:
    code = cli.main(
        ["factors", "--precip", "620", "--sand", "60", "--silt", "35",
         "--clay", "20", "--oc", "1.2"]
    )
    assert code == 0  # factors does not enforce the texture sum; K still computes
    payload = json.loads(capsys.readouterr().out)
    assert payload["K_t_ha_h_ha_MJ_mm"] > 0


def test_cli_missing_dem(
    tmp_path: Path, params_file: Path, capsys: pytest.CaptureFixture[str]
) -> None:
    code = cli.main(["run", str(tmp_path / "nope.tif"), str(params_file), "-o", str(tmp_path)])
    assert code == 2
    assert "not found" in capsys.readouterr().err
