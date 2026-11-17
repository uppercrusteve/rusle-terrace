from __future__ import annotations

from pathlib import Path

import numpy as np
import pytest

from conftest import BASE
from rusle_terrace.factors import (
    ParameterError,
    k_factor,
    load_params,
    p_factor,
    params_from_mapping,
    r_factor,
)


def test_r_factor_value_and_units() -> None:
    assert r_factor(620) == pytest.approx(0.0483 * 620**1.61)


def test_r_factor_monotonic_in_precip() -> None:
    values = [r_factor(p) for p in (300, 450, 620, 900, 1400)]
    assert all(b > a for a, b in zip(values, values[1:], strict=False))


@pytest.mark.parametrize("bad", [0, -5, float("nan")])
def test_r_factor_rejects_non_positive(bad: float) -> None:
    with pytest.raises(ParameterError):
        r_factor(bad)


def test_k_factor_in_expected_range() -> None:
    textures = [
        (45, 35, 20, 1.2),
        (70, 15, 15, 0.5),
        (20, 55, 25, 2.5),
        (35, 40, 25, 1.0),
    ]
    for sand, silt, clay, oc in textures:
        assert 0.01 <= k_factor(sand, silt, clay, oc) <= 0.7


def test_k_factor_reference_value() -> None:
    # Hand-computed from the EPIC equation with C = 0.012.
    assert k_factor(45, 35, 20, 1.2) == pytest.approx(0.04258, abs=2e-4)


def test_k_factor_decreases_with_organic_carbon() -> None:
    assert k_factor(45, 35, 20, 0.4) > k_factor(45, 35, 20, 3.0)


def test_k_factor_coarse_sand_is_less_erodible() -> None:
    assert k_factor(80, 10, 10, 1.0) < k_factor(35, 45, 20, 1.0)


def test_k_factor_rejects_zero_texture_denominator() -> None:
    with pytest.raises(ParameterError):
        k_factor(45, 0, 0, 1.0)


def test_p_factor_none_is_unity() -> None:
    assert np.all(p_factor(np.array([1.0, 8.0, 30.0]), "none") == 1.0)


def test_p_factor_table_terraces() -> None:
    slopes = np.array([0.0, 4.99, 5.0, 9.99, 10.0, 14.99, 15.0, 40.0])
    got = p_factor(slopes, "terraces")
    assert got.tolist() == pytest.approx([0.15, 0.15, 0.20, 0.20, 0.25, 0.25, 0.30, 0.30])


def test_p_factor_table_contour() -> None:
    slopes = np.array([2.0, 7.0, 12.0, 20.0])
    assert p_factor(slopes, "contour").tolist() == pytest.approx([0.55, 0.60, 0.65, 0.70])


def test_p_factor_rejects_unknown_practice() -> None:
    with pytest.raises(ParameterError):
        p_factor(10.0, "agroforestry")


def test_params_parse_mapping_and_carbon_fraction() -> None:
    params = params_from_mapping(BASE)
    assert params.organic_carbon_fraction == pytest.approx(0.012)
    assert params.support_practice == "terraces"
    assert params.slope_classes_deg[-1] == 90


def test_params_missing_required_key() -> None:
    data = dict(BASE)
    del data["clay_pct"]
    with pytest.raises(ParameterError, match="clay_pct"):
        params_from_mapping(data)


def test_params_rejects_texture_sum_outside_tolerance() -> None:
    data = dict(BASE, sand_pct=70.0)
    with pytest.raises(ParameterError, match="sum to 100"):
        params_from_mapping(data)


def test_params_rejects_bad_cover_factor() -> None:
    with pytest.raises(ParameterError, match="cover_factor"):
        params_from_mapping(dict(BASE, cover_factor=1.4))


def test_params_rejects_bad_slope_classes() -> None:
    with pytest.raises(ParameterError, match="start at 0"):
        params_from_mapping(dict(BASE, slope_classes_deg=[5, 10, 20]))
    with pytest.raises(ParameterError, match="increasing"):
        params_from_mapping(dict(BASE, slope_classes_deg=[0, 10, 5, 90]))


def test_params_rejects_non_numeric() -> None:
    with pytest.raises(ParameterError, match="annual_precip_mm"):
        params_from_mapping(dict(BASE, annual_precip_mm="wet"))


def test_load_params_missing_file(tmp_path: Path) -> None:
    with pytest.raises(ParameterError, match="not found"):
        load_params(tmp_path / "nope.yaml")


def test_load_params_invalid_yaml(tmp_path: Path) -> None:
    path = tmp_path / "bad.yaml"
    path.write_text("annual_precip_mm: [unclosed\n", encoding="utf-8")
    with pytest.raises(ParameterError, match="invalid YAML"):
        load_params(path)


def test_load_params_success(tmp_path: Path) -> None:
    path = tmp_path / "ok.yaml"
    path.write_text(
        "\n".join(f"{k}: {v}" for k, v in BASE.items()) + "\n", encoding="utf-8"
    )
    params = load_params(path)
    assert params.annual_precip_mm == 620.0
