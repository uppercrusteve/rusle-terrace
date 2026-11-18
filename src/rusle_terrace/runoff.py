"""DEM loading, LS/A computation and slope-class statistics."""

from __future__ import annotations

from dataclasses import dataclass
from pathlib import Path

import numpy as np
import rasterio
from scipy import ndimage

from .factors import Params, k_factor, p_factor, r_factor
from .flow import d8_directions, flow_accumulation, slope_degrees

NODATA_VALUE = -9999.0
LS_SLOPE_ZERO = 0.01
M_EXP = 0.5
N_EXP = 1.3
_METRE_UNITS = {"metre", "meter", "m"}

# Simplified SL190-2007 erosion severity, t/ha/yr (upper bound exclusive).
SEVERITY_BANDS: tuple[tuple[str, float, float], ...] = (
    ("微度 (very slight)", 0.0, 2.0),
    ("轻度 (slight)", 2.0, 5.0),
    ("中度 (moderate)", 5.0, 10.0),
    ("强烈 (strong)", 10.0, 25.0),
    ("极强烈 (very strong)", 25.0, np.inf),
)


class DemError(ValueError):
    """Raised for unusable DEM inputs (missing CRS, wrong CRS or grid)."""


@dataclass(frozen=True)
class Dem:
    """A validated single-band DEM on a metric projected grid."""

    data: np.ndarray
    valid: np.ndarray
    transform: rasterio.Affine
    crs: rasterio.crs.CRS
    path: Path

    @property
    def cellsize(self) -> float:
        return float(self.transform.a)

    @property
    def shape(self) -> tuple[int, int]:
        return self.data.shape


def read_dem(path: str | Path) -> Dem:
    """Read and validate a DEM.

    Rejects missing CRS, geographic CRS, non-metre units and non-square pixels.
    """
    p = Path(path)
    if not p.is_file():
        raise DemError(f"DEM not found: {p}")
    with rasterio.open(p) as src:
        if src.crs is None:
            raise DemError("DEM has no CRS; assign a projected CRS in metres before running")
        if not src.crs.is_projected:
            raise DemError(f"DEM CRS must be projected, got geographic CRS {src.crs}")
        units = (src.crs.linear_units or "").lower()
        factor = src.crs.linear_units_factor[1] if src.crs.linear_units_factor else 1.0
        if units not in _METRE_UNITS or abs(factor - 1.0) > 1e-9:
            raise DemError(f"DEM CRS must use metres, got linear unit {units!r}")
        if src.count != 1:
            raise DemError(f"DEM must have exactly one band, got {src.count}")
        transform = src.transform
        if transform.a <= 0 or not np.isclose(abs(transform.a), abs(transform.e)):
            raise DemError("DEM must have square, positive-resolution pixels")
        if src.dtypes[0] not in ("float32", "float64", "int16", "int32", "uint16", "uint8"):
            raise DemError(f"unsupported DEM dtype {src.dtypes[0]}")
        arr = src.read(1, masked=True)
        crs = src.crs
    data = np.asarray(arr.filled(np.nan), dtype="float64")
    valid = np.isfinite(data)
    if not valid.any():
        raise DemError("DEM contains no valid elevation cells")
    return Dem(data=data, valid=valid, transform=transform, crs=crs, path=p)


def read_c_factor(path: str | Path, dem: Dem, fallback: float) -> tuple[np.ndarray, int]:
    """Read a C-factor raster aligned to the DEM; invalid cells fall back to ``fallback``.

    Returns the array and the number of cells that were replaced, so the
    caller can report the substitution instead of hiding it.
    """
    p = Path(path)
    if not p.is_file():
        raise DemError(f"C-factor raster not found: {p}")
    with rasterio.open(p) as src:
        if src.shape != dem.shape:
            raise DemError(f"C-factor raster shape {src.shape} differs from DEM {dem.shape}")
        if not np.allclose(
            [src.transform.a, src.transform.e, src.transform.c, src.transform.f],
            [dem.transform.a, dem.transform.e, dem.transform.c, dem.transform.f],
        ):
            raise DemError("C-factor raster grid does not match the DEM")
        values = np.asarray(src.read(1, masked=True).filled(np.nan), dtype="float64")
    invalid = ~np.isfinite(values)
    replaced = int(invalid.sum())
    c = np.where(invalid, fallback, values)
    if np.nanmin(c) < 0 or np.nanmax(c) > 1:
        raise DemError("C-factor values must be within 0-1")
    return c, replaced


def _fill_nearest(data: np.ndarray, valid: np.ndarray) -> np.ndarray:
    filled = data.copy()
    if valid.all():
        return filled
    indices = ndimage.distance_transform_edt(~valid, return_distances=False, return_indices=True)
    filled[~valid] = data[tuple(indices[:, ~valid])]
    return filled


def ls_factor(dem: np.ndarray, cellsize: float, valid: np.ndarray) -> np.ndarray:
    """Slope-length/steepness factor LS (dimensionless) from D8 unit contributing area.

    ``LS = (m+1) * (As/22.13)^m * (sin(beta)/0.0896)^n`` with m=0.5, n=1.3,
    ``As`` in m2/m. Flat cells get the 0.01 lower bound. Invalid cells become NaN.
    """
    filled = _fill_nearest(np.asarray(dem, dtype="float64"), valid)
    slope = slope_degrees(filled, cellsize)
    directions = d8_directions(filled, cellsize)
    acc = flow_accumulation(directions)
    contributing_area = acc * cellsize  # cell area / cell width = cell width
    beta = np.radians(slope)
    ls = (
        (M_EXP + 1.0)
        * (contributing_area / 22.13) ** M_EXP
        * (np.sin(beta) / 0.0896) ** N_EXP
    )
    ls = np.where(slope <= 0.0, LS_SLOPE_ZERO, ls)
    ls[~valid] = np.nan
    return ls


def severity_codes(erosion: np.ndarray) -> np.ndarray:
    """Integer severity codes 0-4; ``-1`` where erosion is NaN."""
    a = np.asarray(erosion, dtype="float64")
    codes = np.full(a.shape, -1, dtype="int8")
    good = np.isfinite(a)
    bounds = [band[2] for band in SEVERITY_BANDS[:-1]]
    codes[good] = np.digitize(a[good], bounds)
    return codes


def severity_shares(erosion: np.ndarray) -> list[dict[str, float | str]]:
    """Area share (percent of valid cells) per severity band."""
    a = np.asarray(erosion, dtype="float64")
    good = np.isfinite(a)
    total = int(good.sum())
    rows: list[dict[str, float | str]] = []
    lower = 0.0
    for label, _, upper in SEVERITY_BANDS:
        mask = good & (a >= lower) & (a < upper)
        pct = 100.0 * mask.sum() / total if total else 0.0
        rows.append(
            {
                "severity": label,
                "lower_t_ha_yr": lower,
                "upper_t_ha_yr": upper,
                "area_pct": float(pct),
            }
        )
        lower = upper
    return rows


def slope_class_stats(
    slope: np.ndarray,
    erosion: np.ndarray,
    valid: np.ndarray,
    classes: tuple[float, ...],
    cellsize: float,
) -> list[dict[str, float | int | str]]:
    """Per-slope-class area share, mean and P90 erosion, and severity by class mean."""
    cell_area_ha = cellsize * cellsize / 10_000.0
    total_cells = int(valid.sum())
    total_area = total_cells * cell_area_ha
    idx = np.digitize(slope, np.asarray(classes[1:-1]), right=False)
    rows: list[dict[str, float | int | str]] = []
    for k, (lo, hi) in enumerate(zip(classes, classes[1:], strict=False)):
        mask = valid & (idx == k)
        count = int(mask.sum())
        area = count * cell_area_ha
        values = erosion[mask]
        values = values[np.isfinite(values)]
        if values.size:
            mean = float(values.mean())
            p90 = float(np.percentile(values, 90))
            breaks = [band[2] for band in SEVERITY_BANDS[:-1]]
            sev = SEVERITY_BANDS[int(np.digitize([mean], breaks)[0])][0]
        else:
            mean = float("nan")
            p90 = float("nan")
            sev = "无数据 (no data)"
        rows.append(
            {
                "class_index": k,
                "class_lo_deg": float(lo),
                "class_hi_deg": float(hi),
                "cell_count": count,
                "area_ha": float(area),
                "area_pct": float(100.0 * area / total_area) if total_area else 0.0,
                "mean_erosion_t_ha_yr": mean,
                "p90_erosion_t_ha_yr": p90,
                "severity_by_mean": sev,
            }
        )
    return rows


@dataclass
class Analysis:
    """Computed erosion surface and statistics for one DEM/parameter pair."""

    dem: Dem
    params: Params
    slope_deg: np.ndarray
    ls: np.ndarray
    c_factor: np.ndarray
    p_factor: np.ndarray
    erosion: np.ndarray
    r: float
    k: float
    stats: list[dict[str, float | int | str]]
    severity: np.ndarray
    notes: tuple[str, ...] = ()

    @property
    def valid_erosion(self) -> np.ndarray:
        return self.erosion[np.isfinite(self.erosion)]


def analyze(
    dem_path: str | Path,
    params: Params,
    c_raster: str | Path | None = None,
    p_override: float | None = None,
) -> Analysis:
    """Run the full RUSLE chain for one DEM."""
    dem = read_dem(dem_path)
    filled = _fill_nearest(dem.data, dem.valid)
    slope = slope_degrees(filled, dem.cellsize)
    ls = ls_factor(dem.data, dem.cellsize, dem.valid)
    notes: list[str] = []
    if c_raster:
        c, replaced = read_c_factor(c_raster, dem, params.cover_factor)
        if replaced:
            notes.append(
                f"C 因子栅格有 {replaced} 个无效像元，已用默认值 {params.cover_factor:g} 填充。"
            )
    else:
        c = np.full(dem.shape, params.cover_factor, dtype="float64")
    if p_override is not None:
        if not 0.0 <= p_override <= 1.0:
            raise ValueError("--p-factor must be within 0-1")
        p = np.full(dem.shape, float(p_override), dtype="float64")
    else:
        p = p_factor(slope, params.support_practice)
    r = r_factor(params.annual_precip_mm)
    k = k_factor(params.sand_pct, params.silt_pct, params.clay_pct, params.organic_carbon_pct)
    erosion = r * k * ls * c * p
    erosion[~dem.valid] = np.nan
    stats = slope_class_stats(slope, erosion, dem.valid, params.slope_classes_deg, dem.cellsize)
    return Analysis(
        dem=dem,
        params=params,
        slope_deg=slope,
        ls=ls,
        c_factor=c,
        p_factor=p,
        erosion=erosion,
        r=r,
        k=k,
        stats=stats,
        severity=severity_codes(erosion),
        notes=tuple(notes),
    )


def write_erosion(path: str | Path, analysis: Analysis) -> Path:
    """Write the erosion raster as float32 GeoTIFF with -9999 NoData."""
    out = Path(path)
    profile = {
        "driver": "GTiff",
        "height": analysis.dem.shape[0],
        "width": analysis.dem.shape[1],
        "count": 1,
        "dtype": "float32",
        "crs": analysis.dem.crs,
        "transform": analysis.dem.transform,
        "nodata": NODATA_VALUE,
        "compress": "deflate",
    }
    data = np.where(np.isfinite(analysis.erosion), analysis.erosion, NODATA_VALUE).astype("float32")
    with rasterio.open(out, "w", **profile) as dst:
        dst.write(data, 1)
        dst.set_band_description(1, "RUSLE A (t/ha/yr), planning screening")
    return out
