"""Primitive rings on the cutoff neighbour list.

``pydseams.band.ring_census`` counts Franzblau rings on
``Frame.bonds_by_index``. For this oxygen-only frame that list is
the cutoff list, not the mutual four-nearest list the CLI uses for
CHILL+.
"""

import pytest

pydseams = pytest.importorskip("pydseams")


def test_primitive_ring_census(tip4p):
    from pydseams.band import ring_census

    frame = pydseams.read(tip4p, atom_type=2)
    assert ring_census(frame, max_size=8) == {3: 1, 4: 230, 5: 4, 6: 192}
