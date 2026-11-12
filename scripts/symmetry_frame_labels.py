#!/usr/bin/env python3
# Copyright (c) d-SEAMS developers
# SPDX-License-Identifier: MIT
"""Write one CSV of per-atom ice labels for each nucleation or interface frame.

Inputs already in the tree, or the 200 ps seed frames used by the symmetry
notes. kind is ice Ih (0), ice Ic (1), ice Isd (2), or -1.
"""

from __future__ import annotations

import sys
from pathlib import Path

import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "bindings" / "python"))
import symmetry


def read_lammpstrj(path: Path):
    lines = path.read_text().splitlines()
    n = int(lines[3])
    box = np.diag(
        [float(lines[5 + i].split()[1]) - float(lines[5 + i].split()[0]) for i in range(3)]
    )
    ids = np.empty(n, dtype=np.int64)
    xyz = np.empty((n, 3))
    for i, line in enumerate(lines[9 : 9 + n]):
        parts = line.split()
        ids[i] = int(parts[0])
        xyz[i] = [float(parts[2]), float(parts[3]), float(parts[4])]
    return ids, np.ascontiguousarray(xyz), box


def read_seed(path: Path):
    lines = path.read_text().splitlines()
    side = float(lines[1].split()[1])
    xyz = np.array([[float(v) for v in ln.split()[1:4]] for ln in lines[2:]])
    ids = np.arange(1, len(xyz) + 1)
    return ids, np.ascontiguousarray(xyz), np.diag([side, side, side])


def write_csv(path: Path, ids, kind):
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w") as out:
        out.write("index,id,kind\n")
        for i, (atom, label) in enumerate(zip(ids, kind)):
            out.write(f"{i},{int(atom)},{int(label)}\n")


def main():
    root = Path(__file__).resolve().parents[1]
    out = root / "benchmarks" / "labels"
    jobs = [
        (root / "input" / "traj" / "nucleation_f1.lammpstrj", read_lammpstrj, out / "nucleation_f1.csv"),
        (root / "input" / "traj" / "interface_f0.lammpstrj", read_lammpstrj, out / "interface_f0.csv"),
    ]
    for seed in (91001, 92002, 93003):
        src = Path(f"/tmp/nuc200/seed-{seed}-200ps.xyz")
        if src.is_file():
            jobs.append((src, read_seed, out / f"seed-{seed}-200ps.csv"))
    for src, reader, dest in jobs:
        ids, xyz, box = reader(src)
        kind = symmetry.classify(xyz, box)
        write_csv(dest, ids, kind)
        ice = int((kind >= 0).sum())
        print(f"{dest.name} n={len(kind)} ice={ice} from {src}")


if __name__ == "__main__":
    main()
