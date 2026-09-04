#!/usr/bin/env python3
"""Write small GenIce frameworks as LAMMPS dumps for cage-signature tests."""

from __future__ import annotations

import argparse
import math
import subprocess
from pathlib import Path

# name -> (genice type, reps)
# sH is GenIce DOH / HS3 (clathrate H), not HS1 (ice Ih / sIV).
FRAMEWORKS = {
    "sod": ("SOD", (3, 3, 3)),
    "fau": ("FAU", (1, 1, 1)),
    "sI": ("CS1", (1, 1, 1)),
    "sII": ("CS2", (1, 1, 1)),
    # 2x1x1 so the 5^12 6^8 cage has 36 distinct vertices (1x1x1 wraps).
    "sH": ("sH", (2, 1, 1)),
}


def wrap_frac(sx: float, sy: float, sz: float) -> tuple[float, float, float]:
    return sx - math.floor(sx), sy - math.floor(sy), sz - math.floor(sz)


def wrap_cartesian(
    x: float, y: float, z: float, lx: float, ly: float, lz: float, tilt
) -> tuple[float, float, float]:
    xy = xz = yz = 0.0
    if tilt is not None:
        xy, xz, yz = tilt
    sz = z / lz
    sy = (y - yz * sz) / ly
    sx = (x - xy * sy - xz * sz) / lx
    sx, sy, sz = wrap_frac(sx, sy, sz)
    return (
        lx * sx + xy * sy + xz * sz,
        ly * sy + yz * sz,
        lz * sz,
    )


def parse_gro_cell(line: str):
    """Gromacs cell in nm: 3 lengths, or 9-vector v1x v2y v3z v1y v1z v2x v2z v3x v3y."""
    vals = [float(x) for x in line.split()]
    if len(vals) == 3:
        return [10.0 * v for v in vals], None
    if len(vals) == 9:
        v1x, v2y, v3z, v1y, v1z, v2x, v2z, v3x, v3y = vals
        if abs(v1y) > 1e-8 or abs(v1z) > 1e-8 or abs(v2z) > 1e-8:
            raise SystemExit(f"cell is not restricted-triclinic: {vals}")
        box = [10.0 * v1x, 10.0 * v2y, 10.0 * v3z]
        tilt = [10.0 * v2x, 10.0 * v3x, 10.0 * v3y]
        if all(abs(t) < 1e-12 for t in tilt):
            return box, None
        return box, tilt
    raise SystemExit(f"unrecognized gro cell {vals}")


def genice(exe: str, kind: str, rep):
    cmd = [exe, kind, "--rep", *map(str, rep), "--format", "gromacs", "--seed", "1"]
    out = subprocess.run(cmd, check=True, capture_output=True, text=True).stdout
    lines = out.splitlines()
    n = int(lines[1])
    pos = []
    for line in lines[2 : 2 + n]:
        if line[10:15].strip().startswith("O"):
            pos.append((float(line[20:28]), float(line[28:36]), float(line[36:44])))
    box, tilt = parse_gro_cell(lines[2 + n])
    atoms = [
        wrap_cartesian(p[0] * 10.0, p[1] * 10.0, p[2] * 10.0, box[0], box[1], box[2], tilt)
        for p in pos
    ]
    return atoms, box, tilt


def write_lammps(path: Path, atoms, box, tilt=None) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", encoding="ascii") as fh:
        fh.write("ITEM: TIMESTEP\n0\n")
        fh.write(f"ITEM: NUMBER OF ATOMS\n{len(atoms)}\n")
        if tilt is None:
            fh.write("ITEM: BOX BOUNDS pp pp pp\n")
            for L in box:
                fh.write(f"0.0 {L:.10f}\n")
        else:
            xy, xz, yz = tilt
            lx, ly, lz = box
            xmin = min(0.0, xy, xz, xy + xz)
            xmax = max(0.0, xy, xz, xy + xz)
            ymin = min(0.0, yz)
            ymax = max(0.0, yz)
            fh.write("ITEM: BOX BOUNDS xy xz yz pp pp pp\n")
            fh.write(f"{xmin:.10f} {lx + xmax:.10f} {xy:.10f}\n")
            fh.write(f"{ymin:.10f} {ly + ymax:.10f} {xz:.10f}\n")
            fh.write(f"0.0 {lz:.10f} {yz:.10f}\n")
        fh.write("ITEM: ATOMS id type x y z\n")
        for i, (x, y, z) in enumerate(atoms, start=1):
            fh.write(f"{i} 1 {x:.8f} {y:.8f} {z:.8f}\n")


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--genice", default="genice2")
    ap.add_argument("--out", type=Path, required=True)
    ap.add_argument("--only", default="", help="comma names to emit (default: all)")
    args = ap.parse_args()
    wanted = {s.strip() for s in args.only.split(",") if s.strip()}
    for name, (kind, rep) in FRAMEWORKS.items():
        if wanted and name not in wanted:
            continue
        atoms, box, tilt = genice(args.genice, kind, rep)
        dest = args.out / f"genice_{name}.lammpstrj"
        write_lammps(dest, atoms, box, tilt)
        print(f"{name} n={len(atoms)} box={box} tilt={tilt} -> {dest}")


if __name__ == "__main__":
    main()
