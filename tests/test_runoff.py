from __future__ import annotations

from pathlib import Path

import numpy as np
import pytest
import rasterio
from rasterio.transform import from_origin

from conftest import BASE, PARAMS_TEXT, south_plane, write_dem
from rusle_terrace.factors import params_from_mapping
from rusle_terrace.runoff import (
    LS_SLOPE_ZERO,
    analyze,
    ls_factor,
    severity_codes,
    slope_class_stats,
    write_erosion,
)

PARAMS = params_from_mapping(BASE)


def analytic_ls(slope_length_m: float, gradient: float, m: float = 0.5, n: float = 1.3) -> float:
    beta = np.arctan(gradient)
    return (m + 1.0) * (slope_length_m / 22.13) ** m * (np.sin(beta) / 0.0896) ** n


def test_ls_matches_analytic_5x5_relative_error_below_5pct() -> None:
    cellsize = 30.0
    x = np.arange(5)[None, :] * cellsize
    z = np.broadcast_to(100.0 - 0.09 * x, (5, 5)).astype("float64")
    valid = np.ones((5, 5), dtype=bool)
    ls = ls_factor(z, cellsize, valid)
    for col in range(5):
        expected = analytic_ls((col + 1) * cellsize, 0.09)
        assert abs(ls[2, col] - expected) / expected < 0.05


def test_ls_matches_analytic_long_plane_half_cell_convention() -> None:
    cellsize = 30.0
    rows = 60
    y = np.arange(rows)[:, None] * cellsize
    z = np.broadcast_to(100.0 - 0.09 * y, (rows, 3)).astype("float64")
    ls = ls_factor(z, cellsize, np.ones_like(z, dtype=bool))
    errors = []
    for row in range(10, 50):
        expected = analytic_ls((row + 0.5) * cellsize, 0.09)
        errors.append(abs(ls[row, 1] - expected) / expected)
    assert np.mean(errors) < 0.05
    assert max(errors) < 0.05


def test_ls_flat_slope_uses_lower_bound_not_inf() -> None:
    ls = ls_factor(np.full((5, 5), 800.0), 30.0, np.ones((5, 5), dtype=bool))
    assert np.all(ls == LS_SLOPE_ZERO)
    assert np.all(np.isfinite(ls))


def test_ls_increases_with_contributing_area() -> None:
    y = np.arange(8)[:, None] * 30.0
    z = np.broadcast_to(100.0 - 0.09 * y, (8, 3)).astype("float64")
    ls = ls_factor(z, 30.0, np.ones_like(z, dtype=bool))
    assert ls[-1, 1] > ls[2, 1] > ls[0, 1]


def test_ls_masks_invalid_cells() -> None:
    y = np.arange(5)[:, None] * 30.0
    z = np.broadcast_to(100.0 - 0.09 * y, (5, 5)).astype("float64")
    valid = np.ones((5, 5), dtype=bool)
    valid[2, 2] = False
    assert np.isnan(ls_factor(z, 30.0, valid)[2, 2])


def test_erosion_equals_factor_product(tmp_path: Path) -> None:
    dem = write_dem(tmp_path / "d.tif", south_plane(6, 6))
    result = analyze(dem, PARAMS)
    expected = result.r * result.k * result.ls * result.c_factor * result.p_factor
    body = slice(1, 5)
    assert np.allclose(result.erosion[body, body], expected[body, body], rtol=1e-6)


def test_p_factor_override(tmp_path: Path) -> None:
    dem = write_dem(tmp_path / "d.tif", south_plane(6, 6))
    base = analyze(dem, PARAMS)
    forced = analyze(dem, PARAMS, p_override=0.9)
    assert np.all(forced.p_factor == 0.9)
    assert base.p_factor[3, 3] == 0.20
    assert forced.erosion[3, 3] == pytest.approx(base.erosion[3, 3] * (0.9 / 0.20))


def test_c_factor_raster_overrides_scalar(tmp_path: Path) -> None:
    dem_path = write_dem(tmp_path / "d.tif", south_plane(6, 6))
    c_path = tmp_path / "c.tif"
    with rasterio.open(
        c_path,
        "w",
        driver="GTiff",
        height=6,
        width=6,
        count=1,
        dtype="float32",
        crs="EPSG:32647",
        transform=from_origin(500_000.0, 2_800_000.0, 30.0, 30.0),
    ) as dst:
        dst.write(np.full((6, 6), 0.5, dtype="float32"), 1)
    base = analyze(dem_path, PARAMS)
    overridden = analyze(dem_path, PARAMS, c_raster=c_path)
    assert np.all(overridden.c_factor == 0.5)
    assert overridden.erosion[3, 3] == pytest.approx(base.erosion[3, 3] * 2.0)


def test_c_factor_nan_cells_are_reported(tmp_path: Path) -> None:
    dem_path = write_dem(tmp_path / "d.tif", south_plane(6, 6))
    data = np.full((6, 6), 0.4, dtype="float32")
    data[2, 3] = np.nan
    c_path = tmp_path / "c_nan.tif"
    with rasterio.open(
        c_path, "w", driver="GTiff", height=6, width=6, count=1, dtype="float32",
        crs="EPSG:32647",
        transform=from_origin(500_000.0, 2_800_000.0, 30.0, 30.0),
        nodata=np.nan,
    ) as dst:
        dst.write(data, 1)
    result = analyze(dem_path, PARAMS, c_raster=c_path)
    assert any("1 个无效像元" in note for note in result.notes)
    assert result.c_factor[2, 3] == pytest.approx(PARAMS.cover_factor)


def test_c_factor_raster_grid_mismatch(tmp_path: Path) -> None:
    dem_path = write_dem(tmp_path / "d.tif", south_plane(6, 6))
    c_path = write_dem(tmp_path / "c.tif", np.full((5, 7), 0.3))
    from rusle_terrace.runoff import DemError

    with pytest.raises(DemError, match="shape"):
        analyze(dem_path, PARAMS, c_raster=c_path)


def test_severity_codes_boundaries() -> None:
    values = np.array([0.0, 1.99, 2.0, 4.99, 5.0, 9.99, 10.0, 24.99, 25.0, 100.0])
    assert severity_codes(values).tolist() == [0, 0, 1, 1, 2, 2, 3, 3, 4, 4]


def test_severity_codes_nan_is_minus_one() -> None:
    assert severity_codes(np.array([np.nan, 3.0])).tolist() == [-1, 1]


def test_slope_class_areas_sum_to_total() -> None:
    cellsize = 30.0
    rows = cols = 20
    slope = np.degrees(np.arctan(0.09)) * np.ones((rows, cols))
    erosion = np.ones((rows, cols))
    valid = np.ones((rows, cols), dtype=bool)
    stats = slope_class_stats(slope, erosion, valid, PARAMS.slope_classes_deg, cellsize)
    total_area = rows * cols * cellsize * cellsize / 10_000.0
    assert len(stats) == len(PARAMS.slope_classes_deg) - 1
    assert sum(row["area_ha"] for row in stats) == pytest.approx(total_area)
    assert sum(row["area_pct"] for row in stats) == pytest.approx(100.0)
    assert stats[1]["cell_count"] == rows * cols  # 5.14 deg falls in 5-10


def test_slope_class_empty_bin_reports_nan() -> None:
    slope = np.full((3, 3), 2.0)
    stats = slope_class_stats(
        slope, np.ones((3, 3)), np.ones((3, 3), dtype=bool), (0.0, 5.0, 90.0), 30.0
    )
    assert np.isnan(stats[1]["mean_erosion_t_ha_yr"])
    assert stats[1]["cell_count"] == 0


def test_write_erosion_uses_nodata_for_masked_cells(tmp_path: Path) -> None:
    dem = write_dem(tmp_path / "d.tif", south_plane(5, 5))
    result = analyze(dem, PARAMS)
    result.erosion[1, 1] = np.nan
    out = write_erosion(tmp_path / "e.tif", result)
    with rasterio.open(out) as src:
        assert src.nodata == -9999.0
        raw = src.read(1)
        assert raw[1, 1] == -9999.0
        assert np.isfinite(raw[3, 3])


def test_analyze_reports_params_and_stats(tmp_path: Path) -> None:
    dem = write_dem(tmp_path / "d.tif", south_plane(5, 5))
    result = analyze(dem, PARAMS)
    assert result.params.support_practice == "terraces"
    assert result.r > 0 and 0.01 <= result.k <= 0.7
    assert len(result.stats) == 5


def test_params_cover_factor_in_erosion(tmp_path: Path) -> None:
    dem = write_dem(tmp_path / "d.tif", south_plane(5, 5))
    low = analyze(dem, params_from_mapping(dict(BASE, cover_factor=0.1)))
    high = analyze(dem, params_from_mapping(dict(BASE, cover_factor=0.4)))
    assert high.erosion[2, 2] == pytest.approx(low.erosion[2, 2] * 4.0)


def test_p_params_file_used(tmp_path: Path) -> None:
    path = tmp_path / "p.yaml"
    path.write_text(PARAMS_TEXT.replace("terraces", "none"), encoding="utf-8")
    from rusle_terrace.factors import load_params

    params = load_params(path)
    assert params.support_practice == "none"
