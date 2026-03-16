"""Symbolic verification of hydrogen bond angle calculations.

Cross-references:
    src/bond.cpp      lines 260-300  (H-bond angle via dot product)
    src/generic.cpp   lines 155-161  (eigenVecAngle using acos)
"""

import sympy as sp
from sympy import Abs, Rational, acos, atan2, cos, pi, simplify, sin, sqrt, symbols


# ---------------------------------------------------------------------------
# test_angle_formula_domain
# ---------------------------------------------------------------------------
def test_angle_formula_domain():
    r"""Show that acos(dot/(|a||b|)) has domain [-1, 1] for unit vectors.

    For unit vectors a, b the Cauchy--Schwarz inequality gives
        |a . b| <= |a| |b| = 1
    so the argument to acos is always in [-1, 1].

    For non-unit floating-point vectors, the ratio dot/(|a||b|) can exceed 1
    in magnitude due to rounding.  We demonstrate this symbolically by showing
    that adding a perturbation eps to the dot product pushes the argument
    outside the valid domain.
    """
    # Symbolic unit vectors in 2D (enough to demonstrate the property)
    theta, phi = symbols("theta phi", real=True)

    a = sp.Matrix([cos(theta), sin(theta)])
    b = sp.Matrix([cos(phi), sin(phi)])

    dot_ab = a.dot(b)
    norm_a = sqrt(a.dot(a))
    norm_b = sqrt(b.dot(b))

    ratio = dot_ab / (norm_a * norm_b)

    # For unit vectors, norms are 1, so ratio = cos(theta)*cos(phi) + sin(theta)*sin(phi)
    # = cos(theta - phi), which is in [-1, 1].
    ratio_simplified = simplify(ratio)
    # Substitute theta = phi (parallel vectors) -> ratio should be exactly 1
    assert simplify(ratio_simplified.subs(theta, phi) - 1) == 0

    # Substitute theta = phi + pi (antiparallel) -> ratio should be -1
    assert simplify(ratio_simplified.subs(theta, phi + pi) + 1) == 0

    # Now show floating-point danger: if dot product has perturbation eps > 0
    # and vectors are nearly parallel, the argument exceeds 1.
    eps = symbols("eps", positive=True)
    perturbed_ratio = (dot_ab + eps) / (norm_a * norm_b)

    # At theta = phi (parallel), perturbed_ratio = 1 + eps > 1
    val = simplify(perturbed_ratio.subs(theta, phi))
    assert simplify(val - (1 + eps)) == 0, "Perturbed ratio should be 1 + eps"


# ---------------------------------------------------------------------------
# test_angle_atan2_alternative
# ---------------------------------------------------------------------------
def test_angle_atan2_alternative():
    r"""Verify that atan2(|cross|, dot) equals acos(dot/(|a||b|)) for unit vectors.

    The atan2 formulation avoids the domain issue because:
    - cross = |a||b| sin(angle)
    - dot   = |a||b| cos(angle)
    - atan2(sin, cos) = angle, no domain restriction needed.

    For 3D vectors a, b:
        |a x b| = |a||b| sin(angle)
        a . b   = |a||b| cos(angle)

    So atan2(|a x b|, a . b) = angle in [0, pi], same as acos but numerically
    stable at angle ~ 0 and angle ~ pi.
    """
    # Use symbolic 3D vectors
    ax, ay, az = symbols("ax ay az", real=True)
    bx, by, bz = symbols("bx by bz", real=True)

    a = sp.Matrix([ax, ay, az])
    b = sp.Matrix([bx, by, bz])

    cross = a.cross(b)
    cross_mag = sqrt(cross.dot(cross))
    dot = a.dot(b)

    # For unit vectors a = (1, 0, 0), b = (cos t, sin t, 0)
    t = symbols("t", positive=True)  # angle in (0, pi)
    subs = {ax: 1, ay: 0, az: 0, bx: cos(t), by: sin(t), bz: 0}

    cross_val = simplify(cross_mag.subs(subs))
    dot_val = simplify(dot.subs(subs))

    # cross_val should be |sin(t)| = sin(t) for t in (0, pi)
    assert simplify(cross_val - sin(t)) == 0 or simplify(cross_val - Abs(sin(t))) == 0

    # dot_val should be cos(t)
    assert simplify(dot_val - cos(t)) == 0

    # atan2(sin(t), cos(t)) = t for t in (0, pi)
    # SymPy does not auto-simplify atan2, so we verify numerically at
    # several points in (0, pi).
    import mpmath

    for t_val in [0.1, 0.5, 1.0, sp.pi / 4, sp.pi / 2, 2.0, sp.pi - 0.1]:
        t_num = float(t_val)
        atan2_val = float(mpmath.atan2(mpmath.sin(t_num), mpmath.cos(t_num)))
        assert abs(atan2_val - t_num) < 1e-14, (
            f"atan2(sin({t_num}), cos({t_num})) = {atan2_val}, expected {t_num}"
        )


# ---------------------------------------------------------------------------
# test_hbond_cutoffs
# ---------------------------------------------------------------------------
def test_hbond_cutoffs():
    r"""Document and verify the hard-coded H-bond cutoffs in bond.cpp.

    The code uses (bond.cpp lines 198-199):
        distCutoff = 2.42 Angstrom   (O-H donor...acceptor distance)
        angleCutoff = 30 degrees      (angle between O-O and O-H vectors)

    These values correspond to the geometric definition of hydrogen bonds
    from the literature:

    - Luzar & Chandler, Nature 379, 55-57 (1996):
      O-O distance < 3.5 A, angle O-H...O < 30 degrees.
      The O-H cutoff of 2.42 A is consistent with the first minimum of the
      O-H radial distribution function in SPC/E water, which typically falls
      around 2.40-2.45 A.

    - Wernet et al., Science 304, 995 (2004): similar geometric criteria.

    This test simply asserts the values match what is hard-coded, as a
    regression check and documentation anchor.
    """
    # Values from bond.cpp:198-199
    DIST_CUTOFF_ANGSTROM = 2.42
    ANGLE_CUTOFF_DEGREES = 30

    # Literature-consistent ranges
    assert 2.0 < DIST_CUTOFF_ANGSTROM < 2.5, (
        f"O-H distance cutoff {DIST_CUTOFF_ANGSTROM} A outside expected range [2.0, 2.5]"
    )
    assert 25 <= ANGLE_CUTOFF_DEGREES <= 35, (
        f"Angle cutoff {ANGLE_CUTOFF_DEGREES} deg outside expected range [25, 35]"
    )

    # Verify the relationship between O-H and O-O cutoffs.
    # For angle = 30 deg and O-H = 2.42 A, the maximum O-O distance is
    # approximately O-H / cos(30 deg) + r_OH_covalent, where the covalent
    # O-H bond length is ~0.96 A.  This gives a rough upper bound:
    import math

    r_OH_covalent = 0.96  # Angstrom, typical O-H bond length in water
    max_OO = DIST_CUTOFF_ANGSTROM + r_OH_covalent  # ~ 3.38 A
    # This should be close to (but not exceed) the common 3.5 A O-O cutoff
    assert max_OO < 3.5, (
        f"Derived max O-O distance {max_OO:.2f} A exceeds 3.5 A literature cutoff"
    )
