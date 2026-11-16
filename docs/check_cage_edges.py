#!/usr/bin/env python3
# Copyright (c) d-SEAMS developers
"""A cage face may not keep an edge longer than 1.5 times that cage's median.

The synthetic case is the periodic-image spike: one edge of a cell length
among hydrogen bonds. When libseams is on the library path, the sH and sII
cages are unwrapped the same way as the figure and must pass the same check.
"""

import os
import sys
from pathlib import Path

import numpy as np

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "scripts"))
sys.path.insert(0, str(ROOT / "docs" / "examples"))
from shared_library import built_library  # noqa: E402

from render_scenes import assert_cage_edges  # noqa: E402


def check_predicate():
    assert_cage_edges([np.array([2.80, 2.81, 2.79, 2.82])])
    try:
        assert_cage_edges([np.array([2.80, 2.81, 12.4])])
    except RuntimeError as exc:
        if "12.400" not in str(exc):
            raise SystemExit(f"the long edge was not named: {exc}")
    else:
        raise SystemExit("an edge 1.5 times the median was accepted")
    print("a cell-length cage edge is rejected")


def _libseams():
    if os.environ.get("SEAMS_LIB"):
        named = built_library(".", "seams", "SEAMS_LIB")
        if named is not None:
            return named
    for entry in os.environ.get("LD_LIBRARY_PATH", "").split(os.pathsep):
        if not entry:
            continue
        found = built_library(entry, "seams")
        if found is not None:
            return found
    return None


def check_structures():
    found = _libseams()
    if found is None:
        print("libseams is absent; the hydrate cages were not unwrapped")
        return
    import science_tutorials as tutorials
    from render_scenes import _box_cell, _positions, _rings, cage_meshes, find_cages

    library = tutorials._lib(found.parent)
    hosts = {
        "sh": (
            ROOT / "input" / "traj" / "genice_sH.lammpstrj",
            [
                ("512", {5: 12}, 1, 12),
                ("sh", {4: 3, 5: 6, 6: 3}, 3, 8),
                ("51268", {5: 12, 6: 8}, 4, 4),
            ],
        ),
        "sii": (
            ROOT / "input" / "traj" / "genice_sII.lammpstrj",
            [
                ("512", {5: 12}, 1, 16),
                ("51264", {5: 12, 6: 4}, 0, 8),
            ],
        ),
    }
    for name, (dump, kinds) in hosts.items():
        cycles = _rings(library, str(dump))
        positions = _positions(dump)
        _origin, vectors = _box_cell(dump)
        groups = []
        for _label, need, color, limit in kinds:
            for faces in find_cages(cycles, need, limit):
                groups.append((faces, color))
        _meshes, _stats, edge_sets = cage_meshes(cycles, groups, positions, vectors)
        longest = max(float(np.max(edges)) for edges in edge_sets)
        print(f"{name} cages {len(edge_sets)} longest edge {longest:.3f} A")


def main():
    check_predicate()
    check_structures()
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
