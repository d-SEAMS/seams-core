"""Property-based tests for ring detection via the nanobind bindings.

Loads the example trajectory, builds the neighbor/H-bond network, finds
rings, and verifies structural invariants that must hold for any valid
primitive ring set.
"""

from pathlib import Path

from pydseamslib import _core

TRAJ = Path("input/traj/exampleTraj.lammpstrj")


def _load_rings(max_depth=6):
    """Helper: load trajectory, build H-bond graph, find rings."""
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
    rings = _core.ringNetwork(nList=hbondsIdx, maxDepth=max_depth)
    return rings, hbondsIdx, cloud


class TestRingProperties:
    """Invariant checks on the ring network."""

    @classmethod
    def setup_class(cls):
        cls.rings, cls.nListIdx, cls.cloud = _load_rings(max_depth=6)

    def test_rings_found(self):
        """At least one ring must exist in the test trajectory."""
        assert len(self.rings) > 0

    def test_ring_sizes_in_range(self):
        """Every ring size must be in [3, maxDepth]."""
        for i, ring in enumerate(self.rings):
            assert 3 <= len(ring) <= 6, (
                f"Ring {i} has size {len(ring)}, expected 3-6"
            )

    def test_no_duplicate_rings(self):
        """No two rings should contain the same set of atoms."""
        seen = set()
        for ring in self.rings:
            key = frozenset(ring)
            assert key not in seen, f"Duplicate ring found: {sorted(ring)}"
            seen.add(key)

    def test_rings_are_valid_cycles(self):
        """Each ring must be a valid cycle in the neighbor graph.

        For every consecutive pair (ring[i], ring[i+1]) and the closing
        edge (ring[-1], ring[0]), at least one must be a neighbor of the
        other in the H-bond index list.
        """
        for ri, ring in enumerate(self.rings):
            n = len(ring)
            for k in range(n):
                a = ring[k]
                b = ring[(k + 1) % n]
                a_neighbors = set(self.nListIdx[a]) if a < len(self.nListIdx) else set()
                b_neighbors = set(self.nListIdx[b]) if b < len(self.nListIdx) else set()
                assert b in a_neighbors or a in b_neighbors, (
                    f"Ring {ri}: edge ({a},{b}) not in H-bond neighbor list"
                )

    def test_ring_indices_in_bounds(self):
        """All atom indices in rings must be valid point cloud indices."""
        n_atoms = self.cloud.nop
        for ring in self.rings:
            for idx in ring:
                assert 0 <= idx < n_atoms, (
                    f"Ring index {idx} out of range [0, {n_atoms})"
                )

    def test_no_self_loops_in_rings(self):
        """No ring should contain the same atom index twice."""
        for i, ring in enumerate(self.rings):
            assert len(ring) == len(set(ring)), (
                f"Ring {i} has duplicate indices: {ring}"
            )


class TestRingDepthVariation:
    """Verify ring counts change sensibly with maxDepth."""

    def test_deeper_search_finds_at_least_as_many(self):
        """maxDepth=6 should find >= rings than maxDepth=4."""
        rings_4, _, _ = _load_rings(max_depth=4)
        rings_6, _, _ = _load_rings(max_depth=6)
        # Rings of size 3-4 should be a subset; overall count should not decrease
        sizes_4 = {frozenset(r) for r in rings_4}
        sizes_6 = {frozenset(r) for r in rings_6}
        assert sizes_4.issubset(sizes_6), (
            "Rings found at maxDepth=4 should also be found at maxDepth=6"
        )
