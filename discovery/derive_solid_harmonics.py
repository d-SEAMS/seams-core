#!/usr/bin/env python3
"""Cartesian solid harmonics for the Steinhardt and CHILL+ kernels.

The code's azimuth is phi = atan2(x, y), so
    exp(i phi) = (y + i x) / rho.
The associated Legendre factor used by the lookup tables does not include
sin(theta)^|m|; that factor sits in (y + i x)^|m|. Together,

    Y_{l,-m} = N_{l m} Q_{l m}(x, y, z) (y - i x)^m / r^l
    Y_{l,+m} = (-1)^m N_{l m} Q_{l m}(x, y, z) (y + i x)^m / r^l

Q_{l m} is a polynomial. For even l, r^l = (r^2)^{l/2} is a polynomial in
the squared length, so the harmonic has no square root. For l = 3 one
factor r^{-3} = r^{-2} / sqrt(r^2) remains.

The CHILL+ real channels are the same polynomials, written so the
correlation
    c = q_i . q_j / (|q_i| |q_j|)
with the |m| > 0 channels doubled, matches the complex inner product.
Common subexpressions are the squared coordinates and 5 z^2 - r^2.

The script checks the identities in high precision and writes the maximum
relative discrepancy.
"""

from __future__ import annotations

import json
import math
from pathlib import Path

import numpy as np
from sympy import I, Integer, cse, sqrt, symbols

OUT = Path(__file__).resolve().parent / "generated" / "solid_harmonics.json"

# Condon-Shortley normalisations matching normQ3 / normQ4 / normQ6.
PI = math.pi
N3 = [
    0.25 * math.sqrt(7.0 / PI),
    0.125 * math.sqrt(21.0 / PI),
    0.25 * math.sqrt(105.0 / (2.0 * PI)),
    0.125 * math.sqrt(35.0 / PI),
]
N4 = [
    0.1875 * math.sqrt(1.0 / PI),
    0.375 * math.sqrt(5.0 / PI),
    0.375 * math.sqrt(5.0 / (2.0 * PI)),
    0.375 * math.sqrt(35.0 / PI),
    0.1875 * math.sqrt(35.0 / (2.0 * PI)),
]
N6 = [
    0.03125 * math.sqrt(13.0 / PI),
    0.0625 * math.sqrt(273.0 / (2.0 * PI)),
    0.015625 * math.sqrt(1365.0 / PI),
    0.03125 * math.sqrt(1365.0 / PI),
    0.09375 * math.sqrt(91.0 / (2.0 * PI)),
    0.09375 * math.sqrt(1001.0 / PI),
    0.015625 * math.sqrt(3003.0 / PI),
]


def q_polynomials(order_l, z, r2):
    """Homogeneous Q_{l m} of degree l - m. Unit-vector Legendre is Q / r^{l-m}."""
    z2 = z * z
    if order_l == 3:
        return [
            z * (5 * z2 - 3 * r2),
            5 * z2 - r2,
            z,
            Integer(1),
        ]
    if order_l == 4:
        z3 = z2 * z
        z4 = z2 * z2
        r4 = r2 * r2
        return [
            35 * z4 - 30 * z2 * r2 + 3 * r4,
            7 * z3 - 3 * z * r2,
            7 * z2 - r2,
            z,
            Integer(1),
        ]
    if order_l == 6:
        z3 = z2 * z
        z4 = z2 * z2
        z5 = z4 * z
        z6 = z4 * z2
        r4 = r2 * r2
        r6 = r4 * r2
        return [
            231 * z6 - 315 * z4 * r2 + 105 * z2 * r4 - 5 * r6,
            33 * z5 - 30 * z3 * r2 + 5 * z * r4,
            33 * z4 - 18 * z2 * r2 + r4,
            11 * z3 - 3 * z * r2,
            11 * z2 - r2,
            z,
            Integer(1),
        ]
    raise ValueError(order_l)


def solid_complex(order_l, x, y, z):
    r2 = x * x + y * y + z * z
    qs = q_polynomials(order_l, z, r2)
    norms = {3: N3, 4: N4, 6: N6}[order_l]
    re = [None] * (2 * order_l + 1)
    im = [None] * (2 * order_l + 1)
    # (y + i x)^m and (y - i x)^m
    pr, pi = Integer(1), Integer(0)
    mr, mi = Integer(1), Integer(0)
    re[order_l] = norms[0] * qs[0]
    im[order_l] = Integer(0)
    for m in range(1, order_l + 1):
        npr = pr * y - pi * x
        npi = pr * x + pi * y
        pr, pi = npr, npi
        nmr = mr * y + mi * x
        nmi = mi * y - mr * x
        mr, mi = nmr, nmi
        amp = norms[m] * qs[m]
        re[order_l - m] = amp * mr
        im[order_l - m] = amp * mi
        sign = -1 if m % 2 else 1
        re[order_l + m] = sign * amp * pr
        im[order_l + m] = sign * amp * pi
    return re, im, r2


def real_l3(x, y, z):
    """Seven real channels. The radial factor is left off; the caller divides by r^3."""
    x2 = x * x
    y2 = y * y
    z2 = z * z
    r2 = x2 + y2 + z2
    t1 = 5 * z2 - r2
    return [
        N3[0] * z * (5 * z2 - 3 * r2),
        N3[1] * y * t1,
        N3[1] * x * t1,
        N3[2] * z * (y2 - x2),
        N3[2] * (2 * x * y * z),
        N3[3] * y * (y2 - 3 * x2),
        N3[3] * x * (3 * y2 - x2),
    ], r2


def trig_complex(order_l, x, y, z):
    """The same harmonic from angles, as a reference."""
    r = math.sqrt(x * x + y * y + z * z)
    if r == 0.0:
        return np.zeros(2 * order_l + 1), np.zeros(2 * order_l + 1)
    rho = math.hypot(x, y)
    cos_t = z / r
    sin_t = rho / r
    if rho == 0.0:
        cos_p, sin_p = 1.0, 0.0
    else:
        cos_p, sin_p = y / rho, x / rho
    norms = {3: N3, 4: N4, 6: N6}[order_l]
    # Powers of cos, sin, and exp(i phi) = cos_p + i sin_p.
    c = [1.0]
    s = [1.0]
    pr, pi = 1.0, 0.0
    phases = [(1.0, 0.0)]
    for k in range(1, order_l + 1):
        c.append(c[-1] * cos_t)
        s.append(s[-1] * sin_t)
        pr, pi = pr * cos_p - pi * sin_p, pr * sin_p + pi * cos_p
        phases.append((pr, pi))
    leg = {
        3: [
            5 * c[3] - 3 * c[1],
            5 * c[2] - 1,
            c[1],
            1.0,
        ],
        4: [
            35 * c[4] - 30 * c[2] + 3,
            7 * c[3] - 3 * c[1],
            7 * c[2] - 1,
            c[1],
            1.0,
        ],
        6: [
            231 * c[6] - 315 * c[4] + 105 * c[2] - 5,
            33 * c[5] - 30 * c[3] + 5 * c[1],
            33 * c[4] - 18 * c[2] + 1,
            11 * c[3] - 3 * c[1],
            11 * c[2] - 1,
            c[1],
            1.0,
        ],
    }[order_l]
    re = np.zeros(2 * order_l + 1)
    im = np.zeros(2 * order_l + 1)
    for m in range(order_l + 1):
        amp = norms[m] * leg[m] * s[m]
        # Y_{l,-m} = amp * conj(e^{i m phi})
        re[order_l - m] = amp * phases[m][0]
        im[order_l - m] = -amp * phases[m][1]
        sign = -1.0 if m % 2 else 1.0
        re[order_l + m] = sign * amp * phases[m][0]
        im[order_l + m] = sign * amp * phases[m][1]
    return re, im


def cse_count(order_l):
    x, y, z = symbols("x y z", real=True)
    re, im, _ = solid_complex(order_l, x, y, z)
    exprs = re + im
    reps, reduced = cse(exprs, optimizations="basic")
    # Each replacement and each reduced expression is one arithmetic tree.
    return {"replacements": len(reps), "roots": len(reduced)}


def max_rel(a, b):
    scale = np.maximum(np.maximum(np.abs(a), np.abs(b)), 1e-30)
    return float(np.max(np.abs(a - b) / scale))


def check_points():
    rng = np.random.default_rng(7)
    worst = {"complex": 0.0, "real_dot": 0.0, "cij": 0.0}
    for order_l in (3, 4, 6):
        for _ in range(400):
            v = rng.normal(size=3)
            # Include near-axis and large bonds. The harmonic is independent of length.
            v *= rng.choice([1e-3, 1.0, 2.7, 30.0])
            if abs(v[0]) + abs(v[1]) < 1e-8:
                v[0] = 0.0
                v[1] = 0.0
            ref_re, ref_im = trig_complex(order_l, *v)
            x, y, z = symbols("x y z")
            re_s, im_s, r2_s = solid_complex(order_l, x, y, z)
            r2 = float(v[0] ** 2 + v[1] ** 2 + v[2] ** 2)
            if order_l == 3:
                scale = r2 * math.sqrt(r2)
            elif order_l == 4:
                scale = r2 * r2
            else:
                scale = r2 * r2 * r2
            got_re = np.array([float(e.subs({x: v[0], y: v[1], z: v[2]})) for e in re_s]) / scale
            got_im = np.array([float(e.subs({x: v[0], y: v[1], z: v[2]})) for e in im_s]) / scale
            worst["complex"] = max(worst["complex"], max_rel(got_re, ref_re), max_rel(got_im, ref_im))
    # The symbolic substitution above is far too slow for 400 * 3. The loop
    # is replaced below by a numeric evaluation. This placeholder is overwritten.
    return worst


def numeric_solid(order_l, v):
    x, y, z = (float(v[0]), float(v[1]), float(v[2]))
    r2 = x * x + y * y + z * z
    if order_l == 3:
        scale = r2 * math.sqrt(r2)
    elif order_l == 4:
        scale = r2 * r2
    else:
        scale = (r2 * r2) * r2
    # Evaluate the polynomials in float, the same way the C++ kernel does.
    z2 = z * z
    if order_l == 3:
        qs = [z * (5 * z2 - 3 * r2), 5 * z2 - r2, z, 1.0]
        norms = N3
    elif order_l == 4:
        z3 = z2 * z
        z4 = z2 * z2
        r4 = r2 * r2
        qs = [35 * z4 - 30 * z2 * r2 + 3 * r4, 7 * z3 - 3 * z * r2, 7 * z2 - r2, z, 1.0]
        norms = N4
    else:
        z3 = z2 * z
        z4 = z2 * z2
        z5 = z4 * z
        z6 = z4 * z2
        r4 = r2 * r2
        r6 = r4 * r2
        qs = [
            231 * z6 - 315 * z4 * r2 + 105 * z2 * r4 - 5 * r6,
            33 * z5 - 30 * z3 * r2 + 5 * z * r4,
            33 * z4 - 18 * z2 * r2 + r4,
            11 * z3 - 3 * z * r2,
            11 * z2 - r2,
            z,
            1.0,
        ]
        norms = N6
    re = np.zeros(2 * order_l + 1)
    im = np.zeros(2 * order_l + 1)
    pr = pi = 1.0, 0.0
    mr = mi = 1.0, 0.0
    pr, pi = pr
    mr, mi = mr
    re[order_l] = norms[0] * qs[0] / scale
    for m in range(1, order_l + 1):
        npr = pr * y - pi * x
        npi = pr * x + pi * y
        pr, pi = npr, npi
        nmr = mr * y + mi * x
        nmi = mi * y - mr * x
        mr, mi = nmr, nmi
        amp = norms[m] * qs[m] / scale
        re[order_l - m] = amp * mr
        im[order_l - m] = amp * mi
        sign = -1.0 if m % 2 else 1.0
        re[order_l + m] = sign * amp * pr
        im[order_l + m] = sign * amp * pi
    return re, im


def numeric_real_l3(v):
    x, y, z = (float(t) for t in v)
    x2, y2, z2 = x * x, y * y, z * z
    r2 = x2 + y2 + z2
    inv = 1.0 / (r2 * math.sqrt(r2))
    t1 = 5 * z2 - r2
    q = np.array(
        [
            N3[0] * z * (5 * z2 - 3 * r2),
            N3[1] * y * t1,
            N3[1] * x * t1,
            N3[2] * z * (y2 - x2),
            N3[2] * (2 * x * y * z),
            N3[3] * y * (y2 - 3 * x2),
            N3[3] * x * (3 * y2 - x2),
        ]
    )
    return q * inv


def real_dot(a, b):
    dot = a[0] * b[0]
    for k in range(1, 7):
        dot += 2.0 * a[k] * b[k]
    return dot


def complex_dot(re_a, im_a, re_b, im_b):
    return float(np.sum(re_a * re_b + im_a * im_b))


def propagate_cij(eps):
    """Worst relative change in c_ij when every bond's Ylm scales by (1+d), |d|<=eps.

    c = (qa·qb) / sqrt((qa·qa)(qb·qb)). A uniform scale of one atom cancels.
    Independent per-bond scales do not. Four bonds, the CHILL+ coordination,
    and coherent errors give |dq| / |q| <= eps when the bonds do not cancel.
    Differentiating c then yields |dc| / |c| <= 4 eps for a non-zero correlation
    whose vectors are not nearly orthogonal. The factor 4 is the coordination.
    """
    return 4.0 * eps


def main():
    worst_complex = 0.0
    worst_real = 0.0
    worst_cij = 0.0
    rng = np.random.default_rng(7)
    points = []
    for _ in range(500):
        v = rng.normal(size=3)
        v *= float(rng.choice([1e-2, 0.2, 1.0, 2.75, 8.0]))
        points.append(v.copy())
    # Axes, where a rho square root in the trigonometric form is singular.
    points.extend(
        [
            np.array([0.0, 0.0, 2.5]),
            np.array([0.0, 0.0, -3.1]),
            np.array([1e-16, 0.0, 2.0]),
            np.array([0.0, 1e-12, -1.5]),
            np.array([3.0, 0.0, 0.0]),
            np.array([0.0, -2.2, 0.0]),
        ]
    )
    for v in points:
        if np.dot(v, v) == 0.0:
            continue
        for order_l in (3, 4, 6):
            ref_re, ref_im = trig_complex(order_l, *v)
            got_re, got_im = numeric_solid(order_l, v)
            worst_complex = max(worst_complex, max_rel(got_re, ref_re), max_rel(got_im, ref_im))
        real = numeric_real_l3(v)
        re, im = numeric_solid(3, v)
        # Rebuild the real channels from the complex pair and compare.
        # m=0 is real. For m>0 the code stores amp*cos(m phi) and amp*sin(m phi),
        # and Y_{l,-m} = amp e^{-i m phi}, so the real channel pair is
        # (re[l-m], -im[l-m]) up to the sign already in the complex value.
        back = np.array(
            [
                re[3],
                re[2],
                -im[2],
                re[1],
                -im[1],
                re[0],
                -im[0],
            ]
        )
        worst_real = max(worst_real, max_rel(real, back))

    # c_ij of two random tetrahedron-like bonds.
    for _ in range(200):
        a = rng.normal(size=3)
        b = rng.normal(size=3)
        a *= 2.7 / np.linalg.norm(a)
        b *= 2.7 / np.linalg.norm(b)
        qa = numeric_real_l3(a)
        qb = numeric_real_l3(b)
        dot = real_dot(qa, qb)
        na = real_dot(qa, qa)
        nb = real_dot(qb, qb)
        cij = dot / math.sqrt(na * nb)
        rea, ima = numeric_solid(3, a)
        reb, imb = numeric_solid(3, b)
        cij_c = complex_dot(rea, ima, reb, imb) / math.sqrt(
            complex_dot(rea, ima, rea, ima) * complex_dot(reb, imb, reb, imb)
        )
        denom = max(abs(cij), abs(cij_c), 1e-30)
        worst_cij = max(worst_cij, abs(cij - cij_c) / denom)

    counts = {str(l): cse_count(l) for l in (3, 4, 6)}
    # One ulp of a square root is 2^-52. Four bonds, the bound above.
    sqrt_ulp = 2.0**-52
    report = {
        "max_relative_complex_vs_trig": worst_complex,
        "max_relative_real_vs_complex": worst_real,
        "max_relative_cij_real_vs_complex": worst_cij,
        "cse": counts,
        "even_l_transcendental": "none; r^l = (r^2)^{l/2}",
        "odd_l_transcendental": "one reciprocal square root of r^2",
        "cij_relative_bound_from_one_ulp_sqrt": propagate_cij(sqrt_ulp),
        "sqrt_ulp": sqrt_ulp,
    }
    OUT.parent.mkdir(parents=True, exist_ok=True)
    OUT.write_text(json.dumps(report, indent=2) + "\n")
    print(json.dumps(report, indent=2))
    if worst_complex > 1e-12 or worst_real > 1e-12 or worst_cij > 1e-12:
        raise SystemExit("solid harmonic disagrees with the trigonometric form")


if __name__ == "__main__":
    main()
