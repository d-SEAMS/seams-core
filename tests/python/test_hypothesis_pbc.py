"""Property-based tests for periodic boundary condition distance calculations.

Uses Hypothesis to generate random inputs and verify mathematical properties
of the minimum image convention distance formula.
"""

import math

import numpy as np
from hypothesis import given, settings, assume
from hypothesis import strategies as st


def periodic_dist_1d(dr, L):
    """Minimum image convention in one dimension."""
    return abs(dr - L * round(dr / L))


def periodic_dist_sq_3d(dx, dy, dz, bx, by, bz):
    """Squared periodic distance in 3D (reference implementation)."""
    ddx = abs(dx)
    ddy = abs(dy)
    ddz = abs(dz)
    ddx -= bx * round(ddx / bx)
    ddy -= by * round(ddy / by)
    ddz -= bz * round(ddz / bz)
    return ddx * ddx + ddy * ddy + ddz * ddz


# Coordinate and box strategies
coord = st.floats(min_value=-50, max_value=50, allow_nan=False, allow_infinity=False)
box_len = st.floats(min_value=1.0, max_value=100.0, allow_nan=False, allow_infinity=False)


# -- 1D PBC properties --


@given(a=coord, b=coord, L=box_len)
def test_pbc_symmetry(a, b, L):
    """periodic_dist(a-b, L) == periodic_dist(b-a, L)."""
    d1 = periodic_dist_1d(a - b, L)
    d2 = periodic_dist_1d(b - a, L)
    assert abs(d1 - d2) < 1e-10


@given(a=coord, L=box_len)
def test_pbc_identity(a, L):
    """periodic_dist(a, a, L) == 0."""
    d = periodic_dist_1d(a - a, L)
    assert abs(d) < 1e-15


@given(a=coord, b=coord, L=box_len)
def test_pbc_bound(a, b, L):
    """Result is always <= L/2."""
    d = periodic_dist_1d(a - b, L)
    assert d <= L / 2 + 1e-10


@given(a=coord, b=coord, L=box_len)
def test_pbc_non_negative(a, b, L):
    """Distance is always non-negative."""
    d = periodic_dist_1d(a - b, L)
    assert d >= -1e-15


@given(
    a=coord, b=coord, c=coord, L=box_len,
)
@settings(max_examples=500)
def test_pbc_triangle_inequality(a, b, c, L):
    """Triangle inequality: d(a,c) <= d(a,b) + d(b,c)."""
    d_ab = periodic_dist_1d(a - b, L)
    d_bc = periodic_dist_1d(b - c, L)
    d_ac = periodic_dist_1d(a - c, L)
    # Triangle inequality holds for the minimum image convention when
    # all pairwise distances are < L/2 (which our bound test guarantees).
    assert d_ac <= d_ab + d_bc + 1e-10


# -- 3D PBC properties --


@given(
    dx=coord, dy=coord, dz=coord,
    bx=box_len, by=box_len, bz=box_len,
)
def test_pbc_3d_non_negative(dx, dy, dz, bx, by, bz):
    """Squared distance is always non-negative."""
    r2 = periodic_dist_sq_3d(dx, dy, dz, bx, by, bz)
    assert r2 >= -1e-15


@given(
    dx=coord, dy=coord, dz=coord,
    bx=box_len, by=box_len, bz=box_len,
)
def test_pbc_3d_sign_invariance(dx, dy, dz, bx, by, bz):
    """Negating all deltas gives the same squared distance."""
    r2_pos = periodic_dist_sq_3d(dx, dy, dz, bx, by, bz)
    r2_neg = periodic_dist_sq_3d(-dx, -dy, -dz, bx, by, bz)
    assert abs(r2_pos - r2_neg) < 1e-10


@given(bx=box_len, by=box_len, bz=box_len)
def test_pbc_3d_zero_displacement(bx, by, bz):
    """Zero displacement gives zero distance."""
    r2 = periodic_dist_sq_3d(0.0, 0.0, 0.0, bx, by, bz)
    assert abs(r2) < 1e-15


@given(
    dx=coord, dy=coord, dz=coord,
    bx=box_len, by=box_len, bz=box_len,
)
def test_pbc_3d_upper_bound(dx, dy, dz, bx, by, bz):
    """Squared distance bounded by (bx/2)^2 + (by/2)^2 + (bz/2)^2."""
    r2 = periodic_dist_sq_3d(dx, dy, dz, bx, by, bz)
    max_r2 = (bx / 2) ** 2 + (by / 2) ** 2 + (bz / 2) ** 2
    assert r2 <= max_r2 + 1e-8


@given(
    dx=coord, dy=coord, dz=coord,
    bx=box_len, by=box_len, bz=box_len,
)
def test_pbc_3d_periodicity(dx, dy, dz, bx, by, bz):
    """Shifting by a full box length does not change the distance."""
    r2_orig = periodic_dist_sq_3d(dx, dy, dz, bx, by, bz)
    r2_shift = periodic_dist_sq_3d(dx + bx, dy + by, dz + bz, bx, by, bz)
    assert abs(r2_orig - r2_shift) < 1e-8
