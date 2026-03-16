"""Property-based tests for hydrogen bond network properties.

Loads the example trajectory via the nanobind bindings and verifies
structural invariants of the H-bond network.
"""

from pathlib import Path

from pydseamslib import _core

TRAJ = Path("input/traj/exampleTraj.lammpstrj")


def _load_hbond_network():
    """Load trajectory and build H-bond network."""
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
    return cloud, nList, hbonds


class TestHbondProperties:
    """Invariant checks on the hydrogen bond network."""

    @classmethod
    def setup_class(cls):
        cls.cloud, cls.nList, cls.hbonds = _load_hbond_network()

    def test_hbond_list_length(self):
        """H-bond list should have one entry per atom."""
        assert len(self.hbonds) == self.cloud.nop

    def test_hbond_count_per_oxygen(self):
        """Each oxygen should have 0-6 hydrogen bonds.

        In ice and liquid water, tetrahedral coordination gives up to 4,
        but defects and the neighbor-list cutoff can yield up to 6.
        """
        for i in range(len(self.hbonds)):
            # First element is self, rest are bonded neighbors
            n_hbonds = len(self.hbonds[i]) - 1
            assert 0 <= n_hbonds <= 6, (
                f"Atom {i} has {n_hbonds} H-bonds (expected 0-6)"
            )

    def test_no_self_bonds(self):
        """No atom should have itself as a bonded neighbor (beyond the self entry).

        The H-bond list convention stores the atom's own ID as the first
        element.  Beyond that, no entry should equal the atom's own ID.
        """
        for i in range(len(self.hbonds)):
            own_id = self.hbonds[i][0]
            neighbors = self.hbonds[i][1:]
            assert own_id not in neighbors, (
                f"Atom {i} (ID {own_id}) has self-bond in neighbor portion"
            )

    def test_hbond_symmetry(self):
        """If atom i lists j as bonded, j should list i.

        Convert to index-based lists first for consistent comparison.
        """
        hbondsIdx = _core.neighbourListByIndex(
            yCloud=self.cloud, nList=self.hbonds
        )
        for i in range(len(hbondsIdx)):
            for j in hbondsIdx[i][1:]:  # skip self
                assert i in hbondsIdx[j], (
                    f"H-bond asymmetry: {j} in hbonds[{i}] but {i} not in hbonds[{j}]"
                )

    def test_bonded_atoms_are_valid_indices(self):
        """All atom IDs in the H-bond list should be valid."""
        hbondsIdx = _core.neighbourListByIndex(
            yCloud=self.cloud, nList=self.hbonds
        )
        n = self.cloud.nop
        for i in range(len(hbondsIdx)):
            for j in hbondsIdx[i]:
                assert 0 <= j < n, (
                    f"H-bond index {j} out of range [0, {n}) for atom {i}"
                )

    def test_hbond_subset_of_neighbor_list(self):
        """H-bonded neighbors should be a subset of the spatial neighbor list.

        Every H-bond pair must also be spatial neighbors (within cutoff).
        """
        nListIdx = _core.neighbourListByIndex(
            yCloud=self.cloud, nList=self.nList
        )
        hbondsIdx = _core.neighbourListByIndex(
            yCloud=self.cloud, nList=self.hbonds
        )
        for i in range(len(hbondsIdx)):
            hb_set = set(hbondsIdx[i][1:])  # skip self
            nl_set = set(nListIdx[i])
            assert hb_set.issubset(nl_set), (
                f"Atom {i}: H-bond neighbors {hb_set - nl_set} not in spatial neighbor list"
            )


class TestHbondWithDifferentCutoffs:
    """Verify H-bond network changes sensibly with different spatial cutoffs."""

    def test_tighter_cutoff_fewer_neighbors(self):
        """A smaller cutoff should yield fewer or equal spatial neighbors."""
        strTRJ = str(TRAJ.absolute())
        cloud = _core.readLammpsTrjreduced(
            filename=strTRJ,
            targetFrame=1,
            typeI=2,
            isSlice=False,
            coordLow=[0, 0, 0],
            coordHigh=[0, 0, 0],
        )
        nList_tight = _core.neighListO(rcutoff=3.0, yCloud=cloud, typeI=2)
        nList_loose = _core.neighListO(rcutoff=3.5, yCloud=cloud, typeI=2)

        tight_idx = _core.neighbourListByIndex(yCloud=cloud, nList=nList_tight)
        loose_idx = _core.neighbourListByIndex(yCloud=cloud, nList=nList_loose)

        total_tight = sum(len(nl) - 1 for nl in tight_idx)
        total_loose = sum(len(nl) - 1 for nl in loose_idx)
        assert total_tight <= total_loose, (
            f"Tighter cutoff gave more neighbors: {total_tight} > {total_loose}"
        )
