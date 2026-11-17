"""One task across the bindings: the Python leg.

CHILL+ on frame 1 of ``mW_cubic.lammpstrj``, LAMMPS type 1, with the
cutoff neighbour list, is 4096 cubic. The R and Julia legs of the same
check run in d-SEAMS/dseams-r and d-SEAMS/Dseams.jl.
"""

import pytest

pydseams = pytest.importorskip("pydseams")


def test_python_cubic_count(cubic):
    frame = pydseams.read(cubic, atom_type=1)
    assert frame.chill_plus().cubic == 4096
