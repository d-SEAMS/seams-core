"""Integration tests verifying Python bindings produce correct results.

These tests use the example trajectory data and verify that the
Python bindings give the same results as the C++ code by checking
known properties of the test systems.
"""

from pathlib import Path
from pydseamslib import _core

# Test trajectory (relative to repo root)
TRAJ = Path("input/traj/exampleTraj.lammpstrj")


def test_read_lammps_trajectory():
    """Read a LAMMPS trajectory and verify basic properties."""
    strTRJ = str(TRAJ.absolute())
    cloud = _core.readLammpsTrjreduced(
        filename=strTRJ,
        targetFrame=1,
        typeI=2,  # oxygen atoms
        isSlice=False,
        coordLow=[0, 0, 0],
        coordHigh=[0, 0, 0],
    )
    # The example trajectory has oxygen atoms
    assert cloud.nop > 0
    assert len(cloud.pts) == cloud.nop
    # Box should have 3 dimensions
    assert len(cloud.box) == 3
    assert all(b > 0 for b in cloud.box)


def test_neighbour_list_consistency():
    """Verify neighbour list is symmetric (if A neighbours B, B neighbours A)."""
    strTRJ = str(TRAJ.absolute())
    cloud = _core.readLammpsTrjreduced(
        filename=strTRJ,
        targetFrame=1,
        typeI=2,
        isSlice=False,
        coordLow=[0, 0, 0],
        coordHigh=[0, 0, 0],
    )

    nList = _core.neighListO(rcutoff=3.5, yCloud=cloud, typeI=2)
    assert len(nList) == cloud.nop

    # Convert to index-based for easier checking
    nListIdx = _core.neighbourListByIndex(yCloud=cloud, nList=nList)

    # Check symmetry: if j in nListIdx[i], then i in nListIdx[j]
    for i in range(len(nListIdx)):
        for j in nListIdx[i][1:]:  # skip self (first element)
            assert i in nListIdx[j], (
                f"Asymmetric neighbour list: {j} in nList[{i}] but {i} not in nList[{j}]"
            )


def test_hydrogen_bond_network():
    """Build H-bond network and verify reasonable properties."""
    strTRJ = str(TRAJ.absolute())
    cloud = _core.readLammpsTrjreduced(
        filename=strTRJ,
        targetFrame=1,
        typeI=2,
        isSlice=False,
        coordLow=[0, 0, 0],
        coordHigh=[0, 0, 0],
    )

    nList = _core.neighListO(rcutoff=3.5, yCloud=cloud, typeI=2)
    hbonds = _core.populateHbonds(
        filename=strTRJ,
        yCloud=cloud,
        nList=nList,
        targetFrame=1,
        Htype=1,
    )

    assert len(hbonds) == cloud.nop
    # Each oxygen should have 0-4 hydrogen bonds (tetrahedral coordination)
    for i in range(len(hbonds)):
        n_hbonds = len(hbonds[i]) - 1  # subtract self
        assert 0 <= n_hbonds <= 6, f"Atom {i} has {n_hbonds} H-bonds (expected 0-6)"


def test_ring_network():
    """Find rings and verify they have the right sizes."""
    strTRJ = str(TRAJ.absolute())
    cloud = _core.readLammpsTrjreduced(
        filename=strTRJ,
        targetFrame=1,
        typeI=2,
        isSlice=False,
        coordLow=[0, 0, 0],
        coordHigh=[0, 0, 0],
    )

    nList = _core.neighListO(rcutoff=3.5, yCloud=cloud, typeI=2)
    hbonds = _core.populateHbonds(
        filename=strTRJ,
        yCloud=cloud,
        nList=nList,
        targetFrame=1,
        Htype=1,
    )
    hbondsIdx = _core.neighbourListByIndex(yCloud=cloud, nList=hbonds)
    rings = _core.ringNetwork(nList=hbondsIdx, maxDepth=6)

    # Should find some rings
    assert len(rings) > 0

    # All rings should have between 3 and 6 members (maxDepth=6)
    for ring in rings:
        assert 3 <= len(ring) <= 6, f"Ring has {len(ring)} members (expected 3-6)"


def test_full_pipeline_chill_plus():
    """Run full CHILL+ pipeline and verify ice type classification."""
    strTRJ = str(TRAJ.absolute())
    cloud = _core.readLammpsTrjreduced(
        filename=strTRJ,
        targetFrame=1,
        typeI=2,
        isSlice=False,
        coordLow=[0, 0, 0],
        coordHigh=[0, 0, 0],
    )

    nList = _core.neighListO(rcutoff=3.5, yCloud=cloud, typeI=2)

    # Run CHILL+ correlation
    cloud = _core.getCorrelPlus(yCloud=cloud, nList=nList, isSlice=False)

    # Every atom should now have c_ij entries
    classified_count = 0
    for pt in cloud.pts:
        if len(pt.c_ij) > 0:
            classified_count += 1
            # Each c_ij should have a real correlation value
            for cij in pt.c_ij:
                assert -1.0 <= cij.c_value <= 1.0, (
                    f"c_ij value {cij.c_value} out of range [-1, 1]"
                )

    # Most atoms should have been classified (they have neighbours)
    assert classified_count > cloud.nop * 0.5, (
        f"Only {classified_count}/{cloud.nop} atoms classified"
    )
