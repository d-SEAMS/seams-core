"""Primitive rings of the hydrogen-bond list.

``pydseams.read`` of a type-2 oxygen dump follows pydseamslib 2.11.1:
``ring_census`` is the Franzblau census of the hydrogen bonds, not
every oxygen pair inside the 3.5 A cutoff.
"""

import pytest

pydseams = pytest.importorskip("pydseams")


def test_primitive_ring_census(tip4p):
    from pydseams.band import ring_census

    frame = pydseams.read(tip4p, atom_type=2)
    assert ring_census(frame, max_size=8) == {3: 1, 4: 230, 5: 4, 6: 192}
