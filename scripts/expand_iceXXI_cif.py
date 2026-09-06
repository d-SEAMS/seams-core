#!/usr/bin/env python3
"""Expand Lee et al. Nat. Mater. 25, 302 (2026) ice XXI CIF to a LAMMPS dump.

The SI CIF is I-42d, Z=152 oxygens. Hydrogens were not located in the XRD.
Reads input/traj/iceXXI_lee2026.cif and writes input/traj/iceXXI_lee2026.lammpstrj.
"""
from __future__ import annotations

import math
import re
from pathlib import Path

OPS = [
    "x, y, z",
    "-x, -y, z",
    "y, -x, -z",
    "-y, x, -z",
    "-x+1/2, y, -z+3/4",
    "x+1/2, -y, -z+3/4",
    "-y+1/2, -x, z+3/4",
    "y+1/2, x, z+3/4",
    "x+1/2, y+1/2, z+1/2",
    "-x+1/2, -y+1/2, z+1/2",
    "y+1/2, -x+1/2, -z+1/2",
    "-y+1/2, x+1/2, -z+1/2",
    "-x+1, y+1/2, -z+5/4",
    "x+1, -y+1/2, -z+5/4",
    "-y+1, -x+1/2, z+5/4",
    "y+1, x+1/2, z+5/4",
]


def parse_frac(tok: str) -> float:
    return float(re.sub(r"\([^)]*\)", "", tok))


def apply_op(op: str, x: float, y: float, z: float) -> tuple[float, float, float]:
    env = {"x": x, "y": y, "z": z}
    out = []
    for part in op.split(","):
        expr = re.sub(
            r"(\d+)/(\d+)",
            lambda m: str(int(m.group(1)) / int(m.group(2))),
            part.strip(),
        )
        out.append(eval(expr, {"__builtins__": {}}, env))
    return out[0], out[1], out[2]


def wrap(u: float) -> float:
    return u - math.floor(u)


def near(a: tuple[float, float, float], b: tuple[float, float, float]) -> bool:
    return all(min(abs(a[i] - b[i]), 1.0 - abs(a[i] - b[i])) < 1e-4 for i in range(3))


def main() -> None:
    root = Path(__file__).resolve().parents[1]
    cif_path = root / "input" / "traj" / "iceXXI_lee2026.cif"
    dump_path = root / "input" / "traj" / "iceXXI_lee2026.lammpstrj"
    cif = cif_path.read_text(encoding="utf-8", errors="replace")
    sites = []
    for line in cif.splitlines():
        m = re.match(r"^O\d+\s+O\s+(\S+)\s+(\S+)\s+(\S+)", line.strip())
        if m:
            sites.append(tuple(parse_frac(g) for g in m.groups()))
    uniq: list[tuple[float, float, float]] = []
    for sx, sy, sz in sites:
        for op in OPS:
            p = tuple(wrap(u) for u in apply_op(op, sx, sy, sz))
            if not any(near(p, q) for q in uniq):
                uniq.append(p)
    if len(uniq) != 152:
        raise SystemExit(f"expected 152 oxygens, got {len(uniq)}")
    a = 20.1966
    c = 7.8912
    cart = sorted((fx * a, fy * a, fz * c) for fx, fy, fz in uniq)
    lines = [
        "ITEM: TIMESTEP",
        "0",
        "ITEM: NUMBER OF ATOMS",
        "152",
        "ITEM: BOX BOUNDS pp pp pp",
        f"0.0 {a:.10f}",
        f"0.0 {a:.10f}",
        f"0.0 {c:.10f}",
        "ITEM: ATOMS id type x y z",
    ]
    for i, (x, y, z) in enumerate(cart, 1):
        lines.append(f"{i} 1 {x:.8f} {y:.8f} {z:.8f}")
    dump_path.write_text("\n".join(lines) + "\n")
    print(dump_path, "n=152")


if __name__ == "__main__":
    main()
