"""Ideal ice Ic frames for the d-SEAMS 1 against 2 size ladder.

Writes one LAMMPS dump per repetition count with the cubic-diamond oxygen
lattice of worked_example.py (O-O distance 2.75 angstrom), 8 m^3 sites for
m conventional cells along each edge.

Usage: python3 ladder_frames.py OUTDIR M [M ...]
"""

import pathlib
import sys

from worked_example import cubic_diamond, write_lammps


def main():
    out = pathlib.Path(sys.argv[1])
    out.mkdir(parents=True, exist_ok=True)
    for m in map(int, sys.argv[2:]):
        points, box = cubic_diamond(m)
        write_lammps(out / f"ic-{len(points)}.lammpstrj", points, box)
        sys.stdout.write(f"ic-{len(points)}.lammpstrj {len(points)} sites\n")


if __name__ == "__main__":
    main()
