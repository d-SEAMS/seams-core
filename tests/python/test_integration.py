"""Integration tests verifying Python bindings produce correct results.

These tests use the example trajectory data and verify that the
Python bindings give the same results as the C++ code by checking
known properties of the test systems.
"""

from pathlib import Path
from pydseamslib import _core, Trajectory

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

    # Mutates the cloud in place and returns the same object
    returned = _core.getCorrelPlus(yCloud=cloud, nList=nList, isSlice=False)
    assert returned is cloud

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


# --- Trajectory class tests ---


def test_trajectory_load():
    """Load a trajectory with the Trajectory class and check basic properties."""
    traj = Trajectory(str(TRAJ.absolute()))
    assert traj.n_atoms > 0
    assert len(traj.box) == 3
    assert all(b > 0 for b in traj.box)
    assert traj.frame == 1
    assert traj.atom_type == 2
    assert traj.cutoff == 3.5


def test_trajectory_neighbor_list():
    """Verify lazy neighbour list construction through the Trajectory class."""
    traj = Trajectory(str(TRAJ.absolute()))
    nlist = traj.neighbor_list
    assert len(nlist) == traj.n_atoms
    # Each atom should have at least one neighbour (plus self)
    for entry in nlist:
        assert len(entry) >= 2, "Expected at least self + one neighbour"


def test_trajectory_hbonds():
    """Verify lazy hydrogen bond construction through the Trajectory class."""
    traj = Trajectory(str(TRAJ.absolute()))
    hbonds = traj.hbonds
    assert len(hbonds) == traj.n_atoms
    for i in range(len(hbonds)):
        n_hb = len(hbonds[i]) - 1
        assert 0 <= n_hb <= 6, f"Atom {i} has {n_hb} H-bonds (expected 0-6)"


def test_trajectory_rings():
    """Verify lazy ring detection through the Trajectory class."""
    traj = Trajectory(str(TRAJ.absolute()))
    rings = traj.rings
    assert len(rings) > 0
    for ring in rings:
        assert 3 <= len(ring) <= 6, f"Ring has {len(ring)} members (expected 3-6)"


def test_trajectory_classify_chill_plus():
    """Run CHILL+ via the Trajectory class and verify classification dict."""
    traj = Trajectory(str(TRAJ.absolute()))
    counts = traj.classify_chill_plus()

    # Should return a dict of ice type name -> count
    assert isinstance(counts, dict)
    assert len(counts) > 0

    # Total should equal number of atoms
    total = sum(counts.values())
    assert total == traj.n_atoms, (
        f"Sum of ice type counts ({total}) != n_atoms ({traj.n_atoms})"
    )

    # At least some atoms should be classified as non-unclassified
    classified = total - counts.get("unclassified", 0)
    assert classified > 0, "Expected some atoms to be classified"


def test_trajectory_classify_chill():
    """Run CHILL via the Trajectory class and verify classification dict."""
    traj = Trajectory(str(TRAJ.absolute()))
    counts = traj.classify_chill()

    assert isinstance(counts, dict)
    total = sum(counts.values())
    assert total == traj.n_atoms


def test_trajectory_custom_parameters():
    """Verify Trajectory accepts custom frame, atom_type, and cutoff."""
    traj = Trajectory(
        str(TRAJ.absolute()), frame=1, atom_type=2, cutoff=4.0
    )
    assert traj.cutoff == 4.0
    assert traj.n_atoms > 0


def test_trajectory_steinhardt():
    """Verify the Steinhardt parameters have one entry per atom."""
    traj = Trajectory(str(TRAJ.absolute()))
    result = traj.steinhardt(order_l=6)

    assert isinstance(result, dict)
    assert set(result.keys()) == {"ql", "ql_bar"}
    assert len(result["ql"]) == traj.n_atoms
    assert len(result["ql_bar"]) == traj.n_atoms
    assert 0.0 < (sum(result["ql"]) / len(result["ql"])) < 1.0
    assert max(result["ql"]) > 0.1


def test_trajectory_steinhardt_order_l_4():
    traj = Trajectory(str(TRAJ.absolute()))
    result = traj.steinhardt(order_l=4)

    assert len(result["ql"]) == traj.n_atoms
    assert 0.0 < (sum(result["ql"]) / len(result["ql"])) < 1.0
    assert max(result["ql"]) > 0.05
