"""Symbolic and numerical verification of spherical harmonics lookup tables.

Verifies the hard-coded Y_lm lookup tables in src/bop.cpp against
SymPy's reference implementation of spherical harmonics.

The C++ code uses:
  phi   = atan2(x, y)   -- azimuthal angle (note: swapped args vs standard)
  theta = acos(z / r)   -- polar angle
  Index mapping: table index m_idx in [0, 2l] corresponds to m_actual = m_idx - l

The lookup tables implement real spherical harmonics with Condon-Shortley phase.
"""

import sympy
from sympy import (
    Abs,
    I,
    Rational,
    conjugate,
    cos,
    exp,
    pi,
    sin,
    sqrt,
    symbols,
    re,
    im,
)
from sympy.functions.special.spherical_harmonics import Ynm
import numpy as np


# ---------------------------------------------------------------------------
# C++ lookup table reimplementation in Python
# ---------------------------------------------------------------------------

def _cpp_q3(m_idx, theta, phi):
    """Reimplementation of sph::lookupTableQ3 from bop.cpp.
    m_idx: 0..6, mapping to m_actual = m_idx - 3.
    theta: polar angle, phi: azimuthal angle.
    """
    PI = np.pi
    if m_idx == 0:  # m = -3
        c = 0.125 * np.sqrt(35 / PI)
        return c * np.exp(-3j * phi) * np.sin(theta)**3
    elif m_idx == 1:  # m = -2
        c = 0.25 * np.sqrt(105 / (2 * PI))
        return c * np.exp(-2j * phi) * np.sin(theta)**2 * np.cos(theta)
    elif m_idx == 2:  # m = -1
        c = 0.125 * np.sqrt(21 / PI)
        return c * np.exp(-1j * phi) * np.sin(theta) * (5 * np.cos(theta)**2 - 1)
    elif m_idx == 3:  # m = 0
        c = 0.25 * np.sqrt(7 / PI)
        return c * (5 * np.cos(theta)**3 - 3 * np.cos(theta))
    elif m_idx == 4:  # m = +1
        c = -0.125 * np.sqrt(21 / PI)
        return c * np.exp(1j * phi) * np.sin(theta) * (5 * np.cos(theta)**2 - 1)
    elif m_idx == 5:  # m = +2
        c = 0.25 * np.sqrt(105 / (2 * PI))
        return c * np.exp(2j * phi) * np.sin(theta)**2 * np.cos(theta)
    elif m_idx == 6:  # m = +3
        c = -0.125 * np.sqrt(35 / PI)
        return c * np.exp(3j * phi) * np.sin(theta)**3
    return 0j


def _cpp_q6(m_idx, theta, phi):
    """Reimplementation of sph::lookupTableQ6 from bop.cpp.
    m_idx: 0..12, mapping to m_actual = m_idx - 6.
    """
    PI = np.pi
    ct = np.cos(theta)
    st = np.sin(theta)
    if m_idx == 0:  # m = -6
        c = 0.015625 * np.sqrt(3003 / PI)
        return c * np.exp(-6j * phi) * st**6
    elif m_idx == 1:  # m = -5
        c = 0.09375 * np.sqrt(1001 / PI)
        return c * np.exp(-5j * phi) * st**5 * ct
    elif m_idx == 2:  # m = -4
        c = 0.09375 * np.sqrt(91 / (2 * PI))
        return c * np.exp(-4j * phi) * st**4 * (11 * ct**2 - 1)
    elif m_idx == 3:  # m = -3
        c = 0.03125 * np.sqrt(1365 / PI)
        return c * np.exp(-3j * phi) * st**3 * (11 * ct**3 - 3 * ct)
    elif m_idx == 4:  # m = -2
        c = 0.015625 * np.sqrt(1365 / PI)
        return c * np.exp(-2j * phi) * st**2 * (33 * ct**4 - 18 * ct**2 + 1)
    elif m_idx == 5:  # m = -1
        c = 0.0625 * np.sqrt(273 / (2 * PI))
        return c * np.exp(-1j * phi) * st * (33 * ct**5 - 30 * ct**3 + 5 * ct)
    elif m_idx == 6:  # m = 0
        c = 0.03125 * np.sqrt(13 / PI)
        return c * (231 * ct**6 - 315 * ct**4 + 105 * ct**2 - 5)
    elif m_idx == 7:  # m = +1
        c = -0.0625 * np.sqrt(273 / (2 * PI))
        return c * np.exp(1j * phi) * st * (33 * ct**5 - 30 * ct**3 + 5 * ct)
    elif m_idx == 8:  # m = +2
        c = 0.015625 * np.sqrt(1365 / PI)
        return c * np.exp(2j * phi) * st**2 * (33 * ct**4 - 18 * ct**2 + 1)
    elif m_idx == 9:  # m = +3
        c = -0.03125 * np.sqrt(1365 / PI)
        return c * np.exp(3j * phi) * st**3 * (11 * ct**3 - 3 * ct)
    elif m_idx == 10:  # m = +4
        c = 0.09375 * np.sqrt(91 / (2 * PI))
        return c * np.exp(4j * phi) * st**4 * (11 * ct**2 - 1)
    elif m_idx == 11:  # m = +5
        c = -0.09375 * np.sqrt(1001 / PI)
        return c * np.exp(5j * phi) * st**5 * ct
    elif m_idx == 12:  # m = +6
        c = 0.015625 * np.sqrt(3003 / PI)
        return c * np.exp(6j * phi) * st**6
    return 0j


# ---------------------------------------------------------------------------
# Reference: SymPy spherical harmonics
# ---------------------------------------------------------------------------

def _sympy_ynm_numerical(l_val, m_val, theta_val, phi_val):
    """Evaluate SymPy Ynm at numerical points.

    SymPy Ynm(n, m, theta, phi) uses physics convention:
      theta = polar angle from z-axis [0, pi]
      phi = azimuthal angle [0, 2*pi)
    and includes the Condon-Shortley phase.
    """
    theta_s, phi_s = symbols("theta phi", real=True)
    expr = Ynm(l_val, m_val, theta_s, phi_s).expand(func=True)
    f = sympy.lambdify((theta_s, phi_s), expr, modules=["numpy"])
    return f(theta_val, phi_val)


# ---------------------------------------------------------------------------
# Test 1: Q3 lookup table vs SymPy
# ---------------------------------------------------------------------------

def test_q3_lookup_table():
    """Compare each of the 7 Q3 entries (l=3, m=-3..+3) against SymPy Ynm."""
    rng = np.random.default_rng(42)
    n_pts = 100
    theta_vals = rng.uniform(0.1, np.pi - 0.1, size=n_pts)
    phi_vals = rng.uniform(0, 2 * np.pi, size=n_pts)

    l_val = 3
    for m_idx in range(7):
        m_actual = m_idx - 3

        cpp_vals = np.array([
            _cpp_q3(m_idx, theta_vals[j], phi_vals[j])
            for j in range(n_pts)
        ])

        sympy_vals = _sympy_ynm_numerical(l_val, m_actual, theta_vals, phi_vals)

        np.testing.assert_allclose(
            cpp_vals.real, sympy_vals.real, atol=1e-10,
            err_msg=f"Q3 real part mismatch at m_idx={m_idx} (m={m_actual})"
        )
        np.testing.assert_allclose(
            cpp_vals.imag, sympy_vals.imag, atol=1e-10,
            err_msg=f"Q3 imag part mismatch at m_idx={m_idx} (m={m_actual})"
        )


# ---------------------------------------------------------------------------
# Test 2: Q6 lookup table vs SymPy
# ---------------------------------------------------------------------------

def test_q6_lookup_table():
    """Compare each of the 13 Q6 entries (l=6, m=-6..+6) against SymPy Ynm."""
    rng = np.random.default_rng(43)
    n_pts = 100
    theta_vals = rng.uniform(0.1, np.pi - 0.1, size=n_pts)
    phi_vals = rng.uniform(0, 2 * np.pi, size=n_pts)

    l_val = 6
    for m_idx in range(13):
        m_actual = m_idx - 6

        cpp_vals = np.array([
            _cpp_q6(m_idx, theta_vals[j], phi_vals[j])
            for j in range(n_pts)
        ])

        sympy_vals = _sympy_ynm_numerical(l_val, m_actual, theta_vals, phi_vals)

        np.testing.assert_allclose(
            cpp_vals.real, sympy_vals.real, atol=1e-10,
            err_msg=f"Q6 real part mismatch at m_idx={m_idx} (m={m_actual})"
        )
        np.testing.assert_allclose(
            cpp_vals.imag, sympy_vals.imag, atol=1e-10,
            err_msg=f"Q6 imag part mismatch at m_idx={m_idx} (m={m_actual})"
        )


# ---------------------------------------------------------------------------
# Test 3: Condon-Shortley consistency
# ---------------------------------------------------------------------------

def test_condon_shortley_consistency():
    """Verify Y_{l,-m} = (-1)^m * conj(Y_{l,m}) for the lookup tables.

    This is a fundamental identity of spherical harmonics with
    Condon-Shortley phase convention.
    """
    rng = np.random.default_rng(44)
    n_pts = 200
    theta_vals = rng.uniform(0.1, np.pi - 0.1, size=n_pts)
    phi_vals = rng.uniform(0, 2 * np.pi, size=n_pts)

    # Q3: l=3
    for m in range(1, 4):  # m = 1, 2, 3
        m_neg_idx = 3 - m   # index for -m
        m_pos_idx = 3 + m   # index for +m

        for j in range(n_pts):
            y_neg = _cpp_q3(m_neg_idx, theta_vals[j], phi_vals[j])
            y_pos = _cpp_q3(m_pos_idx, theta_vals[j], phi_vals[j])
            expected = (-1)**m * np.conj(y_neg)
            np.testing.assert_allclose(
                [y_pos.real, y_pos.imag],
                [expected.real, expected.imag],
                atol=1e-14,
                err_msg=f"Q3 Condon-Shortley failed for m={m}, point {j}"
            )

    # Q6: l=6
    for m in range(1, 7):  # m = 1..6
        m_neg_idx = 6 - m
        m_pos_idx = 6 + m

        for j in range(n_pts):
            y_neg = _cpp_q6(m_neg_idx, theta_vals[j], phi_vals[j])
            y_pos = _cpp_q6(m_pos_idx, theta_vals[j], phi_vals[j])
            expected = (-1)**m * np.conj(y_neg)
            np.testing.assert_allclose(
                [y_pos.real, y_pos.imag],
                [expected.real, expected.imag],
                atol=1e-14,
                err_msg=f"Q6 Condon-Shortley failed for m={m}, point {j}"
            )


# ---------------------------------------------------------------------------
# Test 4: normalization (addition theorem)
# ---------------------------------------------------------------------------

def test_normalization():
    """Verify sum_m |Y_lm(theta, phi)|^2 = (2l+1)/(4*pi) for l=3 and l=6.

    This is the spherical harmonics addition theorem evaluated at zero
    angular separation.
    """
    rng = np.random.default_rng(45)
    n_pts = 500
    theta_vals = rng.uniform(0.01, np.pi - 0.01, size=n_pts)
    phi_vals = rng.uniform(0, 2 * np.pi, size=n_pts)

    # l = 3: sum over 7 terms should equal 7/(4*pi)
    expected_q3 = 7.0 / (4.0 * np.pi)
    for j in range(n_pts):
        total = sum(
            abs(_cpp_q3(m_idx, theta_vals[j], phi_vals[j]))**2
            for m_idx in range(7)
        )
        np.testing.assert_allclose(
            total, expected_q3, atol=1e-12,
            err_msg=f"Q3 normalization failed at point {j}"
        )

    # l = 6: sum over 13 terms should equal 13/(4*pi)
    expected_q6 = 13.0 / (4.0 * np.pi)
    for j in range(n_pts):
        total = sum(
            abs(_cpp_q6(m_idx, theta_vals[j], phi_vals[j]))**2
            for m_idx in range(13)
        )
        np.testing.assert_allclose(
            total, expected_q6, atol=1e-12,
            err_msg=f"Q6 normalization failed at point {j}"
        )
