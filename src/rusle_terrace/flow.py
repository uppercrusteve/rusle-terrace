"""D8 flow routing, flow accumulation and slope from a DEM."""

from __future__ import annotations

from collections import deque

import numpy as np

# Clockwise from east: E, SE, S, SW, W, NW, N, NE.
DIRECTIONS: tuple[tuple[int, int], ...] = (
    (0, 1),
    (1, 1),
    (1, 0),
    (1, -1),
    (0, -1),
    (-1, -1),
    (-1, 0),
    (-1, 1),
)
DIRECTION_NAMES = ("E", "SE", "S", "SW", "W", "NW", "N", "NE")
_DISTANCE = np.array([1.0, np.sqrt(2.0), 1.0, np.sqrt(2.0), 1.0, np.sqrt(2.0), 1.0, np.sqrt(2.0)])
NODATA_DIRECTION = np.int8(-1)


def d8_directions(dem: np.ndarray, cellsize: float) -> np.ndarray:
    """Steepest-descent D8 direction codes (int8), ``-1`` on flats and outlets.

    Drop is ``(z - z_neighbour) / distance`` with a diagonal distance of sqrt(2)
    cell widths, so planar slopes are routed orthogonally as expected. Ties are
    resolved by direction order (east first). No depression filling is applied.
    """
    z = np.asarray(dem, dtype="float64")
    if z.ndim != 2:
        raise ValueError("DEM must be a 2-D array")
    rows, cols = z.shape
    padded = np.full((rows + 2, cols + 2), np.nan)
    padded[1:-1, 1:-1] = z
    best = np.zeros((rows, cols))
    out = np.full((rows, cols), NODATA_DIRECTION, dtype="int8")
    for code, ((dr, dc), dist) in enumerate(zip(DIRECTIONS, _DISTANCE, strict=True)):
        neighbour = padded[1 + dr : 1 + dr + rows, 1 + dc : 1 + dc + cols]
        drop = np.where(np.isnan(neighbour), -np.inf, (z - neighbour) / (dist * cellsize))
        better = drop > best
        best[better] = drop[better]
        out[better] = code
    return out


def flow_accumulation(directions: np.ndarray) -> np.ndarray:
    """Number of cells draining through each cell (including itself).

    Uses a topological queue over the D8 graph. Cells inside a routing cycle
    (possible without depression filling) keep their partial accumulation and
    are not propagated further.
    """
    dirs = np.asarray(directions, dtype="int8")
    if dirs.ndim != 2:
        raise ValueError("direction grid must be a 2-D array")
    rows, cols = dirs.shape
    flat = dirs.ravel()
    n = flat.size
    acc = np.ones(n, dtype="float64")
    target = np.full(n, -1, dtype="int64")
    for code, (dr, dc) in enumerate(DIRECTIONS):
        idx = np.flatnonzero(flat == code)
        if idx.size == 0:
            continue
        r, c = np.divmod(idx, cols)
        rr, cc = r + dr, c + dc
        inside = (rr >= 0) & (rr < rows) & (cc >= 0) & (cc < cols)
        target[idx[inside]] = rr[inside] * cols + cc[inside]
    indegree = np.zeros(n, dtype="int64")
    routed = target[target >= 0]
    np.add.at(indegree, routed, 1)
    queue = deque(np.flatnonzero((indegree == 0) & (flat != NODATA_DIRECTION)).tolist())
    # Outlets (no downslope neighbour) must start the queue even though indegree > 0.
    queue.extend(np.flatnonzero(flat == NODATA_DIRECTION).tolist())
    seen = np.zeros(n, dtype=bool)
    seen[list(queue)] = True
    while queue:
        i = queue.popleft()
        j = target[i]
        if j < 0:
            continue
        acc[j] += acc[i]
        indegree[j] -= 1
        if indegree[j] == 0 and not seen[j]:
            seen[j] = True
            queue.append(j)
    # Cycle fallback: propagate once in reverse index order so no mass is lost.
    remaining = np.flatnonzero(~seen)
    for i in remaining[::-1]:
        j = target[i]
        if j >= 0:
            acc[j] += acc[i]
    return acc.reshape(rows, cols)


def slope_degrees(dem: np.ndarray, cellsize: float) -> np.ndarray:
    """Slope angle in degrees from central differences (metric grid)."""
    z = np.asarray(dem, dtype="float64")
    if z.ndim != 2:
        raise ValueError("DEM must be a 2-D array")
    if not np.isfinite(cellsize) or cellsize <= 0:
        raise ValueError("cellsize must be positive")
    dz_dy, dz_dx = np.gradient(z, cellsize)
    return np.degrees(np.arctan(np.hypot(dz_dx, dz_dy)))
