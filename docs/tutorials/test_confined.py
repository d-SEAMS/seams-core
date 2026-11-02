"""Confined water on the in-tree mixed TIP4P dump.

The 2020 nanotube and monolayer trajectories are figshare project
73545. They are not in this repository. This dump is the confined
TIP4P frame the engine ships.
"""

import math

import pytest

pydseams = pytest.importorskip("pydseams")


def test_in_plane_rdf_is_finite(tip4p, tmp_path):
    frame = pydseams.read(tip4p, atom_type=2)
    r, g = frame.rdf_2d(str(tmp_path), cutoff=6.0, binwidth=0.5)
    assert len(r) == len(g)
    assert r
    assert all(math.isfinite(value) for value in g)


def test_sheet_rings_and_prisms_return(tip4p, tmp_path):
    frame = pydseams.read(tip4p, atom_type=2)
    frame.monolayer_rings(str(tmp_path / "sheet"), sheet_area=1.0, max_depth=6)
    frame.find_prisms(output_dir=str(tmp_path / "prisms"), max_depth=4)
