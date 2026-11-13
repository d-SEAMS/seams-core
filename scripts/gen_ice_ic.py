#!/usr/bin/env python3
"""Oxygen sites of ideal ice Ic.

The conventional cubic cell holds eight molecules. The lattice constant is
``4 * bond / sqrt(3)`` with ``bond = 2.75`` angstrom, so the nearest oxygen
neighbour sits at 2.75 angstrom and the second shell sits beyond 3.5.
``11**3`` cells are 10648 molecules and ``24**3`` cells are 110592.
The timing harness builds the same coordinates in memory. This script writes
an XYZ file only when an output path is given, so the frames are not stored
in the repository.

    python3 scripts/gen_ice_ic.py 11
    python3 scripts/gen_ice_ic.py 24 ice-ic-24.xyz
"""

from __future__ import annotations

import math
import sys

BOND_A = 2.75


def ice_ic(repetitions: int, bond: float = BOND_A):
    lattice = 4.0 * bond / math.sqrt(3.0)
    basis = (
        (0.0, 0.0, 0.0),
        (0.5, 0.5, 0.0),
        (0.5, 0.0, 0.5),
        (0.0, 0.5, 0.5),
        (0.25, 0.25, 0.25),
        (0.75, 0.75, 0.25),
        (0.75, 0.25, 0.75),
        (0.25, 0.75, 0.75),
    )
    points = []
    for ix in range(repetitions):
        for iy in range(repetitions):
            for iz in range(repetitions):
                for bx, by, bz in basis:
                    points.append(
                        ((ix + bx) * lattice, (iy + by) * lattice, (iz + bz) * lattice)
                    )
    return points, repetitions * lattice


def main(argv: list[str]) -> int:
    repetitions = int(argv[1]) if len(argv) > 1 else 11
    points, side = ice_ic(repetitions)
    print(f"molecules {len(points)} side_A {side:.6f} repetitions {repetitions}")
    if len(argv) > 2:
        path = argv[2]
        with open(path, "w", encoding="utf-8") as handle:
            handle.write(f"{len(points)}\n")
            handle.write(f"ideal ice Ic side {side:.8f} A\n")
            for x, y, z in points:
                handle.write(f"O {x:.8f} {y:.8f} {z:.8f}\n")
        print(path)
    return 0


if __name__ == "__main__":
    raise SystemExit(main(sys.argv))
