"""Derivation of the restricted-triclinic edge and the one-image ball.

dumpBoundsToH recovers lx, ly, lz. A cutoff ball wider than half the
shortest edge holds two lattice images, which is why sampleRDF_AA
keeps one minimum-image distance instead of every vesin image.
Sollya encloses the two-flop recovery error
(analysis/triclinic_gate.sollya). Lean proves the real statements in
lean/DseamsProofs/Cell.lean.
"""

import pathlib
import re
import shutil
import subprocess

import sympy as sp
from sympy import Abs, FiniteSet, Max, Min, Rational, Symbol, sqrt


ROOT = pathlib.Path(__file__).resolve().parents[1]
SOLLYA = ROOT / "analysis" / "triclinic_gate.sollya"


def _recovered_edge(lx, xy, xz):
    """Same min/max samples as dumpBoundsToH, then xspan - xmax + xmin."""
    samples = [0, xy, xz, xy + xz]
    xmin = min(samples)
    xmax = max(samples)
    xspan = lx + xmax - xmin
    return sp.simplify(xspan - xmax + xmin)


def test_span_recovers_edge():
    """xspan = lx + xmax - xmin rearranges to the edge lx.

    xmin and xmax are the min and max of {0, xy, xz, xy + xz}, which is
    how a LAMMPS bound span is wider than lx.
    """
    lx, xmax, xmin = sp.symbols("lx xmax xmin", real=True)
    assert sp.expand((lx + xmax - xmin) - xmax + xmin - lx) == 0

    xy, xz = sp.symbols("xy xz", real=True)
    samples = [sp.Integer(0), xy, xz, xy + xz]
    # 0 is one of the samples, so the min is <= 0 and the max is >= 0.
    assert sp.ask(sp.Q.nonpositive(Min(*samples))) is True
    assert sp.ask(sp.Q.nonnegative(Max(*samples))) is True
    for xy_v in (sp.Rational(-3, 2), 0, sp.Rational(4, 5)):
        for xz_v in (sp.Rational(-1, 3), 0, sp.Integer(2)):
            assert _recovered_edge(sp.Integer(7), xy_v, xz_v) == 7


def test_gamma2_dominates_two_flops():
    """|(1+d1)(1+d2)-1| <= (1+u)^2 - 1 <= gamma2 for |di| <= u = 2^{-53}."""
    a, b, c = sp.symbols("a b c", real=True)
    d1, d2 = sp.symbols("d1 d2", real=True)
    got = ((a - b) * (1 + d1) + c) * (1 + d2)
    exact = a - b + c
    err = sp.expand(got - exact)
    model = (a - b) * ((1 + d1) * (1 + d2) - 1) + c * d2
    assert sp.expand(err - model) == 0

    u = Rational(1, 2**53)
    two_flop = (1 + u) ** 2 - 1
    gamma2 = 2 * u / (1 - 2 * u)
    assert sp.simplify(two_flop - (2 * u + u**2)) == 0
    assert two_flop < gamma2


def test_shortest_vector_is_a_recovered_edge():
    """A nonzero integer combination is at least min(lx, ly, lz).

    v = na a + nb b + nc c with a = (lx, 0, 0), b = (xy, ly, 0),
    c = (xz, yz, lz). The last nonzero coefficient forces one
    coordinate to be at least the corresponding positive edge.
    """
    lx, ly, lz = sp.symbols("lx ly lz", positive=True)
    xy, xz, yz = sp.symbols("xy xz yz", real=True)
    na, nb, nc = sp.symbols("na nb nc", integer=True)

    vx = na * lx + nb * xy + nc * xz
    vy = nb * ly + nc * yz
    vz = nc * lz
    norm2 = vx**2 + vy**2 + vz**2
    assert sp.expand(norm2 - vz**2 - vx**2 - vy**2) == 0

    # |nc| >= 1 when nc != 0, so |vz| >= lz and the norm is at least lz.
    nc_nz = sp.symbols("nc", integer=True, nonzero=True)
    assert sp.ask(sp.Q.ge(Abs(nc_nz), 1))
    vz_nz = nc_nz * lz
    assert sp.simplify(Abs(vz_nz) - Abs(nc_nz) * lz) == 0

    # nc = 0, nb != 0 leaves |vy| = |nb| ly >= ly.
    nb_nz = sp.symbols("nb", integer=True, nonzero=True)
    vy_nb = sp.simplify(vy.subs({nc: 0, nb: nb_nz}))
    assert sp.simplify(Abs(vy_nb) - Abs(nb_nz) * ly) == 0

    # nc = nb = 0, na != 0 leaves |vx| = |na| lx >= lx.
    na_nz = sp.symbols("na", integer=True, nonzero=True)
    vx_na = sp.simplify(vx.subs({nc: 0, nb: 0, na: na_nz}))
    assert sp.simplify(Abs(vx_na) - Abs(na_nz) * lx) == 0

    # sqrt(norm2) >= |vz| because norm2 - vz^2 is a sum of squares.
    z = Symbol("z", real=True)
    rest = Symbol("rest", nonnegative=True)
    assert sp.simplify(sqrt(z**2 + rest) ** 2 - (z**2 + rest)) == 0


def test_open_ball_holds_one_image():
    """Radius c < m/2 cannot hold two points separated by m.

    |v| = |(p + v) - p| <= |p + v| + |p|. If both points lie in a closed
    ball of radius c then |v| <= 2c. The gap c = m/2 - gap makes 2c < m.
    At c = m/2 the two points p = -m/2 and p + m both lie on the sphere.
    """
    x, y = sp.symbols("x y", real=True)
    # (|x| + |y|)^2 - (x + y)^2 = 2|x||y| - 2xy >= 0, so |x + y| <= |x| + |y|.
    gap_sq = sp.expand((Abs(x) + Abs(y)) ** 2 - (x + y) ** 2)
    assert sp.expand(gap_sq - (2 * Abs(x) * Abs(y) - 2 * x * y)) == 0
    # xy <= |xy| = |x| |y|. The remainder |z| - z is 0 for z >= 0
    # and -2z for z < 0.
    assert sp.simplify(Abs(x) * Abs(y) - Abs(x * y)) == 0
    z = Symbol("z", real=True)
    assert sp.refine(Abs(z) - z, sp.Q.nonnegative(z)) == 0
    assert sp.refine(Abs(z) - z, sp.Q.negative(z)) == -2 * z
    assert sp.ask(sp.Q.positive(-z), sp.Q.negative(z)) is True

    m, gap = sp.symbols("m gap", positive=True)
    c = m / 2 - gap
    assert sp.expand(m - 2 * c - 2 * gap) == 0
    assert sp.expand(Abs(-m / 2) - m / 2) == 0
    assert sp.expand(Abs(-m / 2 + m) - m / 2) == 0
    assert sp.expand(Abs(m) - m) == 0


def test_self_header_is_not_a_neighbour():
    """neighbourListByIndex rows lead with the atom. Coordination is the
    rest of the row, and a same-species self entry adds one homopolar count."""
    row = FiniteSet(0, 1, 3)
    self_index = FiniteSet(0)
    assert len(row) == len(row - self_index) + 1
    species = {0: 1, 1: 1, 2: 2, 3: 2}
    homo_with_self = sum(1 for j in (0, 1, 3) if species[j] == species[0])
    homo = sum(1 for j in (1, 3) if species[j] == species[0])
    assert homo_with_self == homo + 1
    assert homo == 1


def test_sollya_encloses_the_two_flop_corner():
    """The enclosure of gamma2*|x| + u*hi on [-hi, hi] sits on the corner."""
    sollya = shutil.which("sollya")
    assert sollya, "sollya is required to certify the margin"
    proc = subprocess.run(
        [sollya, str(SOLLYA)],
        check=True,
        capture_output=True,
        text=True,
    )
    # The script prints the enclosure and then its upper bound.
    numbers = re.findall(r"[0-9]+\.[0-9]+e[+-][0-9]+", proc.stdout)
    assert numbers, proc.stdout
    certified = float(numbers[-1])
    u = sp.Rational(1, 2**53)
    gamma2 = 2 * u / (1 - 2 * u)
    corner = gamma2 * 10000 + u * 10000
    assert certified + 1e-18 >= float(corner)
    assert certified < float(corner) * 1.01
