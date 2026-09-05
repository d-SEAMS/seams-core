"""Symbolic comparison of CHILL+ molecule bins vs TUM ring-plane stacking.

Literature I_sd bins CHILL+ cubic/hexagonal *molecules* along z.
TUM stacking bins the centroid of each HC-basal or DDC-equatorial
six-ring. This file proves the geometric identities that make the
ring vote a single layer coordinate, and that a chair hexagon's atoms
have a nonzero z-spread the molecule bin can see.
"""

import sympy as sp
from sympy import Rational, simplify


def test_chair_centroid_is_the_layer_coordinate():
    """Six chair vertices at z = +/- h have centroid 0, independent of h."""
    h = sp.symbols("h", positive=True)
    zs = [h, -h, h, -h, h, -h]
    centroid = simplify(sum(zs) / 6)
    assert centroid == 0
    atom_var = simplify(sum((z - centroid) ** 2 for z in zs) / 6)
    assert atom_var == h ** 2


def test_molecule_bin_can_split_a_chair_tum_cannot():
    """If the CHILL+ origin sits at the chair midplane, +h and -h can
    fall in different bins of width w when h is large enough. The ring
    centroid is identically 0, so TUM casts one vote."""
    h, w, origin = sp.symbols("h w origin", real=True)
    # floor((z - origin)/w) for z = +/- h
    # Choose origin = 0, w = 2*h: +h is in bin 0, -h wraps to bin -1
    plus_bin = sp.floor((h - 0) / (2 * h))
    minus_bin = sp.floor((-h - 0) / (2 * h))
    assert plus_bin.subs({}).doit() == 0
    assert minus_bin.doit() == -1
    centroid_bin = sp.floor(sp.Integer(0) / (2 * h))
    assert centroid_bin == 0


def test_mixed_chill_ring_is_still_one_tum_plane():
    """Three cubic and three hexagonal CHILL+ labels on one basal ring
    give Phi_c = 1/2 in the literature molecule count. TUM sees one
    basal plane (nHex=1, nCubic=0) so the layer letter is H."""
    n_ic, n_ih = 3, 3
    phi_c_molecules = Rational(n_ic, n_ic + n_ih)
    assert phi_c_molecules == Rational(1, 2)
    n_basal, n_equatorial = 1, 0
    phi_c_tum = (
        Rational(n_equatorial, n_basal + n_equatorial)
        if (n_basal + n_equatorial)
        else 0
    )
    assert phi_c_tum == 0
    assert phi_c_tum != phi_c_molecules


def test_clathrate_five_ring_is_not_a_plane():
    """A 5^12 face has no basal/equatorial six-ring vote. TUM Phi_c is
    undefined-as-zero; the literature molecule bin can still invent a
    letter from CHILL+ clathrate labels mapped into cubic/hexagonal."""
    n_basal = 0
    n_equatorial = 0
    n_five = 12
    assert n_five > 0
    phi_c_tum = 0 if (n_basal + n_equatorial) == 0 else Rational(
        n_equatorial, n_basal + n_equatorial
    )
    assert phi_c_tum == 0
