#!/usr/bin/env python3
"""Write small GenIce frameworks as LAMMPS dumps for cage-signature tests."""

from __future__ import annotations

import argparse
import subprocess
from pathlib import Path

# name -> (genice type, reps)
FRAMEWORKS = {
    "sod": ("SOD", (3, 3, 3)),
    "fau": ("FAU", (1, 1, 1)),
    "sI": ("CS1", (1, 1, 1)),
    "sII": ("CS2", (1, 1, 1)),
}


def genice(exe: str, kind: str, rep) -> tuple[list[tuple[float, float, float]], list[float]]:
    cmd = [exe, kind, "--rep", *map(str, rep), "--format", "gromacs", "--seed", "1"]
    out = subprocess.run(cmd, check=True, capture_output=True, text=True).stdout
    lines = out.splitlines()
    n = int(lines[1])
    pos = []
    for line in lines[2 : 2 + n]:
        if line[10:15].strip().startswith("O"):
            pos.append((float(line[20:28]), float(line[28:36]), float(line[36:44])))
    cell = [float(x) for x in lines[2 + n].split()]
    if len(cell) != 3:
        raise SystemExit(f"{kind}: non-orthogonal cell {cell}")
    box = [c * 10.0 for c in cell]
    atoms = [(p[0] * 10.0 % box[0], p[1] * 10.0 % box[1], p[2] * 10.0 % box[2]) for p in pos]
    return atoms, box


def write_lammps(path: Path, atoms, box) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", encoding="ascii") as fh:
        fh.write("ITEM: TIMESTEP\n0\n")
        fh.write(f"ITEM: NUMBER OF ATOMS\n{len(atoms)}\n")
        fh.write("ITEM: BOX BOUNDS pp pp pp\n")
        for L in box:
            fh.write(f"0.0 {L:.10f}\n")
        fh.write("ITEM: ATOMS id type x y z\n")
        for i, (x, y, z) in enumerate(atoms, start=1):
            fh.write(f"{i} 1 {x:.8f} {y:.8f} {z:.8f}\n")


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--genice", default="genice2")
    ap.add_argument("--out", type=Path, required=True)
    args = ap.parse_args()
    for name, (kind, rep) in FRAMEWORKS.items():
        atoms, box = genice(args.genice, kind, rep)
        dest = args.out / f"genice_{name}.lammpstrj"
        write_lammps(dest, atoms, box)
        print(f"{name} n={len(atoms)} box={box} -> {dest}")


if __name__ == "__main__":
    main()
