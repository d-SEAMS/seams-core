"""Symbolic verification of RDF normalization from rdf2d.cpp.

Cross-references:
    src/rdf2d.cpp  lines 162-203  (normalizeRDF)
"""

import sympy as sp
from sympy import Rational, oo, pi, simplify, symbols


# ---------------------------------------------------------------------------
# test_bin_volume_formula
# ---------------------------------------------------------------------------
def test_bin_volume_formula():
    """Verify the spherical shell volume used at rdf2d.cpp:194-195.

    The code computes:
        binVolume = 4*pi * dR^3 * ((i+1)^3 - i^3) / 3

    This must equal the exact shell volume  (4*pi/3)*(R_outer^3 - R_inner^3)
    where R_inner = i*dR and R_outer = (i+1)*dR.
    """
    i, dR = symbols("i dR", positive=True)

    R_inner = i * dR
    R_outer = (i + 1) * dR

    shell_exact = Rational(4, 3) * pi * (R_outer**3 - R_inner**3)
    shell_code = 4 * pi * dR**3 * ((i + 1) ** 3 - i**3) / 3

    assert simplify(shell_exact - shell_code) == 0


# ---------------------------------------------------------------------------
# test_slab_correction_factor
# ---------------------------------------------------------------------------
def test_slab_correction_factor():
    r"""Verify that the slab geometry correction integrates correctly.

    The slab correction factor in rdf2d.cpp:183-191 is:
        f(r) = 1 - r/(2h)     for r <= h
        f(r) = h/(2r)          for r > h

    For a uniform density rho in a slab of height h and plane area A,
    the number of neighbours within distance R of a central particle is:

        N(R) = rho * integral_0^R  4*pi*r^2 * f(r) dr

    For R <= h the integral should give the volume of the intersection of a
    sphere of radius R with a slab of thickness h centred on the particle.
    The exact intersection volume for R <= h is:

        V_cap(R) = 4*pi*R^3/3 - pi*R^4/(2h)

    (sphere minus two spherical-cap corrections at top and bottom of slab).
    """
    r, h, R = symbols("r h R", positive=True)

    factor_low = 1 - r / (2 * h)

    integrand = 4 * pi * r**2 * factor_low
    integral_low = sp.integrate(integrand, (r, 0, R))

    # Expected: full sphere volume minus slab-cap correction
    expected = Rational(4, 3) * pi * R**3 - pi * R**4 / (2 * h)

    assert simplify(integral_low - expected) == 0


# ---------------------------------------------------------------------------
# test_normalization_ideal_gas
# ---------------------------------------------------------------------------
def test_normalization_ideal_gas():
    """For an ideal gas (g(r)=1), the RDF normalization must return 1.

    The code (rdf2d.cpp:197-198) computes for each bin:
        rdfValues[ibin] += histogram[ibin] / (nIter * binVolume * nopA * rho * factor)

    For a single frame (nIter=1) of an ideal gas with N particles in a box
    of volume V = Lx*Ly*Lz, the expected histogram count in a thin shell is:

        histogram[ibin] = N * rho * binVolume * factor

    (the factor of N comes from summing over all reference particles, and
    rho * binVolume * factor is the expected count per reference particle).
    Dividing by (1 * binVolume * N * rho * factor) should yield exactly 1.
    """
    N, rho, bV, f = symbols("N rho bV f", positive=True)

    histogram_ideal = N * rho * bV * f
    g_r = histogram_ideal / (1 * bV * N * rho * f)

    assert simplify(g_r - 1) == 0


# ---------------------------------------------------------------------------
# Extra: verify that the high-r factor smoothly joins the low-r factor at r=h
# ---------------------------------------------------------------------------
def test_factor_continuity_at_boundary():
    """The two branches of the slab factor must match at r = h."""
    r, h = symbols("r h", positive=True)

    factor_low = 1 - r / (2 * h)
    factor_high = h / (2 * r)

    # Evaluate both at r = h
    val_low = factor_low.subs(r, h)
    val_high = factor_high.subs(r, h)

    assert simplify(val_low - val_high) == 0  # both equal 1/2
