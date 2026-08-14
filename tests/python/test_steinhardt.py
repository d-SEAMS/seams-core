"""Steinhardt order parameter tests using a synthetic FCC lattice.

Mirrors the "Steinhardt parameters reproduce the FCC reference values"
case in tests/test_bop.cpp: a perfect, periodic FCC lattice gives every
atom an identical local environment, so the neighbour-averaged qlBar
must coincide with the local ql.
"""

from pydseamslib import _core


def _fcc_cloud(reps, lattice):
    """Build a periodic FCC point cloud of `reps` unit cells per side."""
    basis = [
        (0.0, 0.0, 0.0),
        (0.5, 0.5, 0.0),
        (0.5, 0.0, 0.5),
        (0.0, 0.5, 0.5),
    ]

    cloud = _core.PointCloudDouble()
    pts = []
    id_index_map = {}
    atom_id = 1
    for i in range(reps):
        for j in range(reps):
            for k in range(reps):
                for bx, by, bz in basis:
                    pt = _core.PointDouble()
                    pt.c_type = 1
                    pt.atomID = atom_id
                    pt.molID = atom_id
                    pt.x = (i + bx) * lattice
                    pt.y = (j + by) * lattice
                    pt.z = (k + bz) * lattice
                    pts.append(pt)
                    id_index_map[atom_id] = atom_id - 1
                    atom_id += 1

    cloud.pts = pts
    cloud.nop = len(pts)
    cloud.currentFrame = 1
    box_length = reps * lattice
    cloud.box = [box_length, box_length, box_length]
    cloud.boxLow = [0.0, 0.0, 0.0]
    cloud.idIndexMap = id_index_map
    return cloud


def test_steinhardt_qbar_equals_ql_for_uniform_environment():
    lattice = 4.0
    cloud = _fcc_cloud(reps=4, lattice=lattice)

    # First FCC shell sits at a/sqrt(2); cut between that and the second
    # shell at a.
    cutoff = 0.85 * lattice
    nList = _core.neighListO(rcutoff=cutoff, yCloud=cloud, typeI=1)

    result = _core.steinhardtQl(yCloud=cloud, nList=nList, orderL=6)

    assert len(result.ql) == cloud.nop
    for ql, ql_bar in zip(result.ql, result.qlBar):
        assert abs(ql - 0.574524) <= 1e-5
        assert abs(ql_bar - ql) <= 1e-9, (
            f"qlBar ({ql_bar}) != ql ({ql}) in a uniform FCC environment"
        )


def test_steinhardt_q4_matches_fcc_reference():
    lattice = 4.0
    cloud = _fcc_cloud(reps=4, lattice=lattice)
    cutoff = 0.85 * lattice
    nList = _core.neighListO(rcutoff=cutoff, yCloud=cloud, typeI=1)

    result = _core.steinhardtQl(yCloud=cloud, nList=nList, orderL=4)

    assert len(result.ql) == cloud.nop
    for ql, ql_bar in zip(result.ql, result.qlBar):
        assert abs(ql - 0.190941) <= 1e-5
        assert abs(ql_bar - ql) <= 1e-9
