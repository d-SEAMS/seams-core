# Copyright (c) d-SEAMS developers
"""The Python symmetry front matches the C API on an ideal ring and on ice Ic."""

import os
import sys

try:
    import numpy  # noqa: F401
except ImportError:
    raise SystemExit(77)

sys.path.insert(0, os.path.join(os.path.dirname(__file__), "..", "bindings", "python"))

import symmetry


def main():
    ring = symmetry.template("ring-6")
    # Swap two atoms. The match still has to name the hexagon.
    ring[[0, 1]] = ring[[1, 0]]
    hit = symmetry.match(ring)
    if hit["name"] != "ring-6":
        raise SystemExit(f"ring matched {hit['name']}")
    if hit["rmsd"] > 1e-6:
        raise SystemExit(f"ring rmsd {hit['rmsd']}")
    if hit["permutation"].shape != (6,):
        raise SystemExit("permutation length")
    group = symmetry.point_group(symmetry.template("ring-6"))
    if symmetry.backend() == "ira":
        if group is None or group[0] != "D6h":
            raise SystemExit(f"point group {group}")
    lattice, xyz = symmetry.crystal("ic")
    labels = symmetry.classify(xyz, lattice)
    if not (labels == 1).all():
        raise SystemExit(f"ice Ic labels {labels}")
    tight = symmetry.classify(xyz, lattice, shell_rmsd=0.05)
    if not (tight == 1).all():
        raise SystemExit(f"tight ice Ic labels {tight}")
    space = symmetry.space_group(lattice, xyz)
    if space is not None and space[0] != 227:
        raise SystemExit(f"space group {space}")
    if space is None and symmetry.backend():
        # spglib is independent of the matcher. A missing module returns None.
        pass
    print("ok", hit["name"], symmetry.backend(), space)


if __name__ == "__main__":
    main()
