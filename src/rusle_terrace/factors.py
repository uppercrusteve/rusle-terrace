"""RUSLE factor equations and parameter-file handling.

All factors follow the brief exactly; units are documented in ``docs/methods.md``.

References
----------
Renard & Freimund (1994) for R; Williams (1995) EPIC formulation for K;
Wischmeier & Smith (1978) and Renard et al. (1997) for C/P practice values.
"""

from __future__ import annotations

from dataclasses import dataclass, field
from pathlib import Path
from typing import Any

import numpy as np
import yaml

PRACTICES = ("none", "contour", "terraces")

# P-factor tables, keyed by practice. Bins are [0,5), [5,10), [10,15), [15, inf) degrees.
P_TABLE: dict[str, tuple[float, float, float, float]] = {
    "none": (1.0, 1.0, 1.0, 1.0),
    "contour": (0.55, 0.60, 0.65, 0.70),
    "terraces": (0.15, 0.20, 0.25, 0.30),
}
P_BREAKS_DEG = (5.0, 10.0, 15.0)


class ParameterError(ValueError):
    """Raised when the parameter file or a supplied value is out of range."""


def r_factor(annual_precip_mm: float) -> float:
    """Rainfall erosivity R in MJ mm ha-1 h-1 yr-1 (Renard & Freimund 1994, eq. 2)."""
    p = float(annual_precip_mm)
    if not np.isfinite(p) or p <= 0:
        raise ParameterError("annual_precip_mm must be a positive finite number")
    return 0.0483 * p**1.61


def k_factor(sand_pct: float, silt_pct: float, clay_pct: float, organic_carbon_pct: float) -> float:
    """Soil erodibility K in t ha h ha-1 MJ-1 mm-1 (EPIC, Williams 1995).

    Inputs are percentages (0-100); organic carbon is the percentage but enters the
    equation as the fraction ``C = organic_carbon_pct / 100``.
    """
    san = float(sand_pct)
    sil = float(silt_pct)
    cla = float(clay_pct)
    oc = float(organic_carbon_pct)
    for name, value in (("sand_pct", san), ("silt_pct", sil), ("clay_pct", cla)):
        if not np.isfinite(value) or not 0.0 <= value <= 100.0:
            raise ParameterError(f"{name} must be within 0-100")
    if not np.isfinite(oc) or not 0.0 <= oc <= 100.0:
        raise ParameterError("organic_carbon_pct must be within 0-100")
    if sil + cla <= 0:
        raise ParameterError("silt_pct + clay_pct must be greater than 0")
    c = oc / 100.0
    sni = 1.0 - san / 100.0
    k_us = (
        0.2
        + 0.3
        * np.exp(-0.0256 * san * (1.0 - sil / 100.0))
        * (sil / (cla + sil)) ** 0.3
        * (1.0 - 0.25 * c / (c + np.exp(3.72 - 2.95 * c)))
        * (1.0 - 0.7 * sni / (sni + np.exp(-5.51 + 22.9 * sni)))
    )
    return float(k_us * 0.1317)


def p_factor(slope_deg: np.ndarray | float, practice: str) -> np.ndarray:
    """Support-practice factor P from slope in degrees.

    ``none`` -> 1.0; ``contour`` -> 0.55/0.60/0.65/0.70; ``terraces`` -> 0.15/0.20/0.25/0.30
    for [0,5), [5,10), [10,15), >=15 degrees respectively.
    """
    if practice not in P_TABLE:
        raise ParameterError(f"support_practice must be one of {PRACTICES}, got {practice!r}")
    slope = np.asarray(slope_deg, dtype="float64")
    b1, b2, b3 = P_BREAKS_DEG
    table = P_TABLE[practice]
    return np.select(
        [slope < b1, slope < b2, slope < b3],
        [table[0], table[1], table[2]],
        default=table[3],
    ).astype("float64")


@dataclass(frozen=True)
class Params:
    """Validated RUSLE model parameters."""

    annual_precip_mm: float
    sand_pct: float
    silt_pct: float
    clay_pct: float
    organic_carbon_pct: float
    cover_factor: float
    support_practice: str = "none"
    slope_classes_deg: tuple[float, ...] = field(
        default=(0.0, 5.0, 10.0, 15.0, 25.0, 90.0)
    )

    def __post_init__(self) -> None:
        if self.support_practice not in PRACTICES:
            raise ParameterError(
                f"support_practice must be one of {PRACTICES}, got {self.support_practice!r}"
            )
        if not 0.0 <= self.cover_factor <= 1.0:
            raise ParameterError("cover_factor must be within 0-1")
        classes = tuple(float(v) for v in self.slope_classes_deg)
        if len(classes) < 2:
            raise ParameterError("slope_classes_deg needs at least two bounds")
        if classes[0] != 0.0:
            raise ParameterError("slope_classes_deg must start at 0")
        if any(b <= a for a, b in zip(classes, classes[1:], strict=False)):
            raise ParameterError("slope_classes_deg must be strictly increasing")
        if classes[-1] < 90.0:
            raise ParameterError("slope_classes_deg must end at 90 or above")

    @property
    def organic_carbon_fraction(self) -> float:
        return self.organic_carbon_pct / 100.0

    def to_dict(self) -> dict[str, Any]:
        return {
            "annual_precip_mm": self.annual_precip_mm,
            "sand_pct": self.sand_pct,
            "silt_pct": self.silt_pct,
            "clay_pct": self.clay_pct,
            "organic_carbon_pct": self.organic_carbon_pct,
            "cover_factor": self.cover_factor,
            "support_practice": self.support_practice,
            "slope_classes_deg": list(self.slope_classes_deg),
        }


_REQUIRED = (
    "annual_precip_mm",
    "sand_pct",
    "silt_pct",
    "clay_pct",
    "organic_carbon_pct",
    "cover_factor",
)


def params_from_mapping(data: dict[str, Any]) -> Params:
    """Build :class:`Params` from a parsed YAML mapping with clear errors."""
    if not isinstance(data, dict):
        raise ParameterError("parameter file must contain a YAML mapping")
    missing = [key for key in _REQUIRED if key not in data]
    if missing:
        raise ParameterError(f"missing required parameter(s): {', '.join(missing)}")
    numeric: dict[str, float] = {}
    for key in _REQUIRED:
        try:
            numeric[key] = float(data[key])
        except (TypeError, ValueError) as exc:
            raise ParameterError(f"{key} must be numeric") from exc
        if not np.isfinite(numeric[key]):
            raise ParameterError(f"{key} must be finite")
    texture = numeric["sand_pct"] + numeric["silt_pct"] + numeric["clay_pct"]
    if abs(texture - 100.0) > 1.0:
        raise ParameterError(
            f"sand_pct + silt_pct + clay_pct must sum to 100 (+/-1), got {texture:.1f}"
        )
    classes = data.get("slope_classes_deg", [0, 5, 10, 15, 25, 90])
    if not isinstance(classes, (list, tuple)):
        raise ParameterError("slope_classes_deg must be a list of numbers")
    try:
        classes_t = tuple(float(v) for v in classes)
    except (TypeError, ValueError) as exc:
        raise ParameterError("slope_classes_deg must contain numbers") from exc
    return Params(
        annual_precip_mm=numeric["annual_precip_mm"],
        sand_pct=numeric["sand_pct"],
        silt_pct=numeric["silt_pct"],
        clay_pct=numeric["clay_pct"],
        organic_carbon_pct=numeric["organic_carbon_pct"],
        cover_factor=numeric["cover_factor"],
        support_practice=str(data.get("support_practice", "none")),
        slope_classes_deg=classes_t,
    )


def load_params(path: str | Path) -> Params:
    """Read and validate a YAML parameter file."""
    p = Path(path)
    if not p.is_file():
        raise ParameterError(f"parameter file not found: {p}")
    try:
        data = yaml.safe_load(p.read_text(encoding="utf-8"))
    except yaml.YAMLError as exc:
        raise ParameterError(f"invalid YAML in {p}: {exc}") from exc
    return params_from_mapping(data)
