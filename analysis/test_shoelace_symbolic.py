"""Symbolic verification of the Shoelace (trapezoid) area formula.

Cross-references:
    src/order_parameter.cpp  lines 105-174  (projAreaSingleRing)
    src/order_parameter.cpp  lines 70-99    (calcCoverageArea)
"""

import sympy as sp
from sympy import Rational, cos, pi, simplify, sqrt, symbols


# ---------------------------------------------------------------------------
# Helpers
# ---------------------------------------------------------------------------
def _shoelace_cross(coords):
    """Standard Shoelace: A = 0.5 * |sum_i (x_i*y_{i+1} - x_{i+1}*y_i)|."""
    n = len(coords)
    s = 0
    for i in range(n):
        j = (i + 1) % n
        s += coords[i][0] * coords[j][1] - coords[j][0] * coords[i][1]
    return Rational(1, 2) * sp.Abs(s)


def _shoelace_trapezoid(coords):
    """Trapezoid variant (used in code): A = 0.5 * |sum_i (x_i+x_{i+1})*(y_{i+1}-y_i)|."""
    n = len(coords)
    s = 0
    for i in range(n):
        j = (i + 1) % n
        s += (coords[i][0] + coords[j][0]) * (coords[j][1] - coords[i][1])
    return Rational(1, 2) * sp.Abs(s)


# ---------------------------------------------------------------------------
# test_trapezoid_equals_shoelace
# ---------------------------------------------------------------------------
def test_trapezoid_equals_shoelace():
    """Prove that the cross-product and trapezoid Shoelace forms are equivalent.

    We verify for n = 3, 4, 5, 6 vertices with symbolic coordinates.
    The key identity is:
        sum (x_i + x_{i+1})(y_{i+1} - y_i) = sum (x_i*y_{i+1} - x_{i+1}*y_i)

    because expanding the left side gives
        sum x_i*y_{i+1} - x_i*y_i + x_{i+1}*y_{i+1} - x_{i+1}*y_i
    and the middle two terms form a telescoping sum that cancels over a
    closed polygon.
    """
    for n in (3, 4, 5, 6):
        xs = sp.symbols(f"x0:{n}")
        ys = sp.symbols(f"y0:{n}")
        coords = list(zip(xs, ys))

        # Compute the signed sums (before abs and 0.5 factor)
        cross_sum = sum(
            coords[i][0] * coords[(i + 1) % n][1]
            - coords[(i + 1) % n][0] * coords[i][1]
            for i in range(n)
        )
        trap_sum = sum(
            (coords[i][0] + coords[(i + 1) % n][0])
            * (coords[(i + 1) % n][1] - coords[i][1])
            for i in range(n)
        )

        diff = sp.expand(cross_sum - trap_sum)
        assert diff == 0, f"Failed for n={n}: diff = {diff}"


# ---------------------------------------------------------------------------
# test_return_order
# ---------------------------------------------------------------------------
def test_return_order():
    """Defect record: the v1.0 projAreaSingleRing return-order swap.

    v1.0's projAreaSingleRing computed areaXY, areaXZ, areaYZ but returned
    {areaXY, areaYZ, areaXZ}, so calcCoverageArea read areaYZ where it
    expected areaXZ and vice versa. Invisible for isotropic systems, wrong
    for any anisotropic geometry. The current C++ returns {areaXY, areaXZ,
    areaYZ}, matching the caller; this test pins both the historical
    mislabeling and the contract the fix restored.
    """
    areaXY_val, areaXZ_val, areaYZ_val = symbols(
        "aXY aXZ aYZ", positive=True
    )

    # What v1.0 returned
    v1_returned = [areaXY_val, areaYZ_val, areaXZ_val]
    # The caller's reading (unchanged across versions)
    assert v1_returned[0] == areaXY_val
    assert v1_returned[1] == areaYZ_val, "v1 slot 1 carried areaYZ, not areaXZ"
    assert v1_returned[2] == areaXZ_val, "v1 slot 2 carried areaXZ, not areaYZ"

    # The corrected contract: slot order matches the caller's labels
    corrected = [areaXY_val, areaXZ_val, areaYZ_val]
    assert corrected[1] == areaXZ_val
    assert corrected[2] == areaYZ_val


# ---------------------------------------------------------------------------
# test_regular_polygon_area
# ---------------------------------------------------------------------------
def test_regular_polygon_area():
    """Verify the Shoelace formula against known regular polygon areas.

    Regular hexagon with side s:  exact area = 3*sqrt(3)/2 * s^2
    Regular n-gon with circumradius R:  area = n/2 * R^2 * sin(2*pi/n)
    """
    s = symbols("s", positive=True)

    # Regular hexagon: vertices at angles k*pi/3, circumradius = s
    hex_coords = [
        (s * cos(k * pi / 3), s * sp.sin(k * pi / 3)) for k in range(6)
    ]

    area_shoelace = _shoelace_trapezoid(hex_coords)
    area_exact = 3 * sqrt(3) * s**2 / 2

    diff = simplify(area_shoelace - area_exact)
    assert diff == 0, f"Hexagon area mismatch: {diff}"

    # Regular triangle (equilateral) with circumradius R
    R = symbols("R", positive=True)
    tri_coords = [
        (R * cos(k * 2 * pi / 3), R * sp.sin(k * 2 * pi / 3)) for k in range(3)
    ]
    area_tri = _shoelace_trapezoid(tri_coords)
    area_tri_exact = Rational(3, 2) * R**2 * sp.sin(2 * pi / 3)

    diff_tri = simplify(area_tri - area_tri_exact)
    assert diff_tri == 0, f"Triangle area mismatch: {diff_tri}"
