from __future__ import annotations

import numpy as np
import pytest

from rusle_terrace.flow import d8_directions, flow_accumulation, slope_degrees

# Codes: 0=E, 1=SE, 2=S, 3=SW, 4=W, 5=NW, 6=N, 7=NE.
SE, E, S = 1, 0, 2


def test_d8_picks_diagonal_on_diagonal_plane() -> None:
    z = np.array(
        [
            [5.0, 4.0, 3.0],
            [4.0, 3.0, 2.0],
            [3.0, 2.0, 1.0],
        ]
    )
    dirs = d8_directions(z, 30.0)
    assert dirs[1, 1] == SE


def test_d8_picks_east_on_east_plane_and_leaves_outlet() -> None:
    z = np.array(
        [
            [3.0, 2.0, 1.0],
            [3.0, 2.0, 1.0],
            [3.0, 2.0, 1.0],
        ]
    )
    dirs = d8_directions(z, 30.0)
    assert dirs[1, 0] == E
    assert dirs[1, 1] == E
    assert dirs[1, 2] == -1


def test_d8_flat_surface_has_no_direction() -> None:
    dirs = d8_directions(np.full((4, 4), 12.0), 30.0)
    assert np.all(dirs == -1)


def test_d8_south_plane_routes_south() -> None:
    y = np.arange(4)[:, None] * 30.0
    z = np.broadcast_to(100.0 - 0.09 * y, (4, 4))
    dirs = d8_directions(z, 30.0)
    assert np.all(dirs[0:3, :] == S)
    assert np.all(dirs[3, :] == -1)


def test_flow_accumulation_conserves_cell_count() -> None:
    y = np.arange(5)[:, None] * 30.0
    z = np.broadcast_to(100.0 - 0.09 * y, (5, 6))
    acc = flow_accumulation(d8_directions(z, 30.0))
    assert acc.sum() == pytest.approx(90.0)  # 6 columns x 15 (5 rows, 1+2+3+4+5)
    assert acc.min() == pytest.approx(1.0)


def test_flow_accumulation_increases_downslope() -> None:
    y = np.arange(7)[:, None] * 30.0
    z = np.broadcast_to(100.0 - 0.09 * y, (7, 3))
    acc = flow_accumulation(d8_directions(z, 30.0))
    column = acc[:, 1]
    assert np.all(np.diff(column) > 0)
    assert column[-1] == pytest.approx(7.0)


def test_flow_accumulation_on_flat_is_all_ones() -> None:
    acc = flow_accumulation(np.full((3, 3), -1, dtype="int8"))
    assert np.all(acc == 1.0)
    assert acc.sum() == 9.0


def test_slope_degrees_uniform_9_percent() -> None:
    x = np.arange(5)[None, :] * 30.0
    z = np.broadcast_to(100.0 - 0.09 * x, (5, 5))
    got = slope_degrees(z, 30.0)
    expected = np.degrees(np.arctan(0.09))
    assert np.allclose(got, expected, atol=1e-9)
    assert got[2, 2] == pytest.approx(5.1427, abs=1e-3)


def test_slope_zero_on_flat() -> None:
    assert np.all(slope_degrees(np.full((3, 3), 5.0), 30.0) == 0.0)


def test_slope_rejects_bad_cellsize() -> None:
    with pytest.raises(ValueError, match="cellsize"):
        slope_degrees(np.zeros((3, 3)), 0.0)


def test_d8_rejects_one_dimensional_input() -> None:
    with pytest.raises(ValueError, match="2-D"):
        d8_directions(np.zeros(5), 30.0)
