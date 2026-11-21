"""Numbers behind the 2020 figures, on the dumps in this tree.

The nucleation, nanotube, monolayer, and in-plane RDF trajectories
are figshare project 73545. This checkout has the cubic mW lattice,
the sI unit, and the mixed TIP4P example.

Python is ``pip install pydseamslib``. On the cubic lattice CHILL+
is 4096 cubic and the cage score is 4096 double-diamond. The same
cubic count is what the R and Julia fronts report.
"""

import pytest

pydseams = pytest.importorskip("pydseams")


def test_cubic_ice_figure(cubic):
    frame = pydseams.read(cubic, frame=1, atom_type=1)
    assert frame.n_atoms == 4096
    assert frame.chill_plus().cubic == 4096
    assert frame.cages().n_ic == 4096


def test_sI_cage_figure(si):
    frame = pydseams.read(si, atom_type=1)
    found = frame.cages_by_signature("512")
    atoms = set()
    for cage in found:
        atoms.update(cage["vertices"])
    assert len(found) == 2
    assert len(atoms) == 40


def test_mixed_frame_disagreement(tip4p):
    frame = pydseams.read(tip4p, atom_type=2)
    ice = frame.chill_plus()
    cages = frame.cages()
    assert ice.water == 241
    assert ice.interClathrate == 9
    assert cages.n_water == 250
