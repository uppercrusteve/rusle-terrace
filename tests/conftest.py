from __future__ import annotations

from pathlib import Path

import numpy as np
import pytest
import rasterio
from rasterio.transform import from_origin

DEFAULT_CRS = "EPSG:32647"


def write_dem(
    path: Path,
    z: np.ndarray,
    cellsize: float = 30.0,
    crs: str | None = DEFAULT_CRS,
    nodata: float | None = None,
) -> Path:
    rows, cols = z.shape
    profile = {
        "driver": "GTiff",
        "height": rows,
        "width": cols,
        "count": 1,
        "dtype": "float32",
        "crs": crs,
        "transform": from_origin(500_000.0, 2_800_000.0, cellsize, cellsize),
        "nodata": nodata,
    }
    with rasterio.open(path, "w", **profile) as dst:
        dst.write(z.astype("float32"), 1)
    return path


def plane(
    rows: int = 5, cols: int = 5, gradient: float = 0.09, cellsize: float = 30.0
) -> np.ndarray:
    """West-to-east uniform descending plane."""
    x = np.arange(cols)[None, :] * cellsize
    return np.broadcast_to(100.0 - gradient * x, (rows, cols)).astype("float64")


def south_plane(
    rows: int = 5, cols: int = 5, gradient: float = 0.09, cellsize: float = 30.0
) -> np.ndarray:
    y = np.arange(rows)[:, None] * cellsize
    return np.broadcast_to(100.0 - gradient * y, (rows, cols)).astype("float64")


BASE = {
    "annual_precip_mm": 620,
    "sand_pct": 45,
    "silt_pct": 35,
    "clay_pct": 20,
    "organic_carbon_pct": 1.2,
    "cover_factor": 0.25,
    "support_practice": "terraces",
}

PARAMS_TEXT = """\
annual_precip_mm: 620
sand_pct: 45
silt_pct: 35
clay_pct: 20
organic_carbon_pct: 1.2
cover_factor: 0.25
support_practice: terraces
slope_classes_deg: [0, 5, 10, 15, 25, 90]
"""


@pytest.fixture
def params_file(tmp_path: Path) -> Path:
    path = tmp_path / "params.yaml"
    path.write_text(PARAMS_TEXT, encoding="utf-8")
    return path


@pytest.fixture
def dem_file(tmp_path: Path) -> Path:
    return write_dem(tmp_path / "dem.tif", south_plane(8, 8))
