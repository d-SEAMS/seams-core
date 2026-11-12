#!/usr/bin/env python3
"""Ideal ice Ih and Ic, plus the liquid frames already in the tree.

The oxygen bond is 2.75 A, the same length as the local-symmetry cell.
Ice Ic is the 8-atom cubic cell. Ice Ih is the 8-atom orthorhombic
cell built from two hexagonal primitives.
128, 4096 and 110592 molecules are rectangular supercells of those
cells. Coordinates stay in memory. Pass --out only when a caller wants
one frame written, and keep that file out of git.
"""

from __future__ import annotations

import argparse
import math
from pathlib import Path

import numpy as np

BOND_A = 2.75
ROOT = Path(__file__).resolve().parents[1]

# (nx, ny, nz). Ice Ic uses the 8-atom cubic cell. Ice Ih uses the
# 8-atom orthorhombic cell made from two hexagonal primitives, so the
# box has no tilt and a three-length box is the same cell.
IC_REPS = {
    128: (4, 2, 2),
    4096: (8, 8, 8),
    110592: (24, 24, 24),
}
IH_REPS = {
    128: (2, 2, 4),
    4096: (8, 8, 8),
    110592: (24, 24, 24),
}
SIZES = (128, 4096, 110592)


def _ic_basis():
    fcc = np.array(
        [[0.0, 0.0, 0.0], [0.5, 0.5, 0.0], [0.5, 0.0, 0.5], [0.0, 0.5, 0.5]],
        dtype=np.float64,
    )
    return np.vstack([fcc, fcc + 0.25])


def _ih_basis():
    z = 1.0 / 16.0
    return np.array(
        [
            [1.0 / 3.0, 2.0 / 3.0, z],
            [2.0 / 3.0, 1.0 / 3.0, 0.5 + z],
            [2.0 / 3.0, 1.0 / 3.0, 1.0 - z],
            [1.0 / 3.0, 2.0 / 3.0, 0.5 - z],
        ],
        dtype=np.float64,
    )


def ice_ic(molecules):
    """Orthorhombic supercell of cubic ice. Lattice rows, LAMMPS box9."""
    nx, ny, nz = IC_REPS[molecules]
    a = 4.0 * BOND_A / math.sqrt(3.0)
    lengths = np.array([nx, ny, nz], dtype=np.float64) * a
    basis = _ic_basis()
    axes = [np.arange(n, dtype=np.float64) for n in (nx, ny, nz)]
    ijk = np.stack(np.meshgrid(*axes, indexing="ij"), axis=-1).reshape(-1, 3)
    frac = ((ijk[:, None, :] + basis[None, :, :]) / np.array([nx, ny, nz])).reshape(-1, 3)
    pos = frac * lengths
    lattice = np.diag(lengths)
    box9 = np.array(
        [0.0, 0.0, 0.0, lengths[0], lengths[1], lengths[2], 0.0, 0.0, 0.0],
        dtype=np.float64,
    )
    return _frame(f"Ic-{molecules}", "Ic", pos, lattice, box9)


def _ih_ortho_basis():
    """Eight oxygens in the orthorhombic cell, wrapped into [0, L)."""
    c = BOND_A / 0.375
    a = c / math.sqrt(8.0 / 3.0)
    hex_vectors = np.array(
        [[a, 0.0, 0.0], [-0.5 * a, a * math.sqrt(3.0) / 2.0, 0.0], [0.0, 0.0, c]],
        dtype=np.float64,
    )
    primitive = _ih_basis() @ hex_vectors
    both = np.vstack([primitive, primitive + hex_vectors[1]])
    lengths = np.array([a, a * math.sqrt(3.0), c], dtype=np.float64)
    # A coordinate that is only a roundoff below 0 must stay at 0, not jump
    # to the far side of the cell.
    wrapped = both - np.floor(both / lengths + 1e-8) * lengths
    return wrapped, lengths


def ice_ih(molecules):
    """Orthorhombic supercell of ice Ih. Eight molecules per cell."""
    nx, ny, nz = IH_REPS[molecules]
    basis, lengths = _ih_ortho_basis()
    if len(basis) != 8:
        raise RuntimeError("ice Ih orthorhombic cell is not 8 molecules")
    axes = [np.arange(n, dtype=np.float64) for n in (nx, ny, nz)]
    ijk = np.stack(np.meshgrid(*axes, indexing="ij"), axis=-1).reshape(-1, 3)
    supercell = lengths * np.array([nx, ny, nz], dtype=np.float64)
    frac = (ijk[:, None, :] * lengths + basis[None, :, :]) / supercell
    pos = frac.reshape(-1, 3) * supercell
    lattice = np.diag(supercell)
    box9 = np.array(
        [0.0, 0.0, 0.0, supercell[0], supercell[1], supercell[2], 0.0, 0.0, 0.0],
        dtype=np.float64,
    )
    return _frame(f"Ih-{molecules}", "Ih", pos, lattice, box9)


def _frame(name, phase, pos, lattice, box9):
    pos = np.ascontiguousarray(pos, dtype=np.float64)
    return {
        "name": name,
        "phase": phase,
        "positions": pos,
        "lattice": np.ascontiguousarray(lattice, dtype=np.float64),
        "box9": np.ascontiguousarray(box9, dtype=np.float64),
        "molecules": int(pos.shape[0]),
    }


def read_lammps_first(path, type_keep, name, phase):
    """First frame of a LAMMPS dump. Orthorhombic bounds only."""
    path = Path(path)
    atoms = []
    bounds = []
    with path.open() as handle:
        line = handle.readline()
        while line and not line.startswith("ITEM: NUMBER OF ATOMS"):
            line = handle.readline()
        n = int(handle.readline())
        line = handle.readline()
        if not line.startswith("ITEM: BOX BOUNDS"):
            raise ValueError(f"{path} has no box bounds")
        if "xy" in line:
            raise ValueError(f"{path} is triclinic; this reader keeps ortho frames")
        for _ in range(3):
            lo, hi = handle.readline().split()[:2]
            bounds.append((float(lo), float(hi)))
        header = handle.readline()
        if not header.startswith("ITEM: ATOMS"):
            raise ValueError(f"{path} has no atom header")
        cols = header.split()[2:]
        ix = cols.index("x")
        iy = cols.index("y")
        iz = cols.index("z")
        it = cols.index("type")
        for _ in range(n):
            parts = handle.readline().split()
            if int(parts[it]) != type_keep:
                continue
            atoms.append((float(parts[ix]), float(parts[iy]), float(parts[iz])))
    if not atoms:
        raise ValueError(f"{path} has no atoms of type {type_keep}")
    pos = np.asarray(atoms, dtype=np.float64)
    lengths = np.array([hi - lo for lo, hi in bounds], dtype=np.float64)
    origin = np.array([lo for lo, _hi in bounds], dtype=np.float64)
    lattice = np.diag(lengths)
    box9 = np.array(
        [origin[0], origin[1], origin[2], lengths[0], lengths[1], lengths[2], 0.0, 0.0, 0.0],
        dtype=np.float64,
    )
    return _frame(name, phase, pos, lattice, box9)


def liquid_frames(root=None):
    root = Path(root) if root is not None else ROOT
    traj = root / "input" / "traj"
    frames = []
    mw = traj / "dump-l.lammpstrj"
    tip = traj / "exampleTraj.lammpstrj"
    if mw.is_file():
        frames.append(read_lammps_first(mw, 1, "liquid-mw", "liquid-mw"))
    if tip.is_file():
        frames.append(read_lammps_first(tip, 2, "tip4p-o", "tip4p-o"))
    return frames


def ideal_frames(sizes=SIZES):
    frames = []
    for n in sizes:
        if n not in IC_REPS:
            raise ValueError(f"no ice cell for {n} molecules")
        frames.append(ice_ic(n))
        frames.append(ice_ih(n))
    return frames


def campaign_frames(sizes=SIZES, root=None):
    return ideal_frames(sizes) + liquid_frames(root)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--sizes", default="128,4096,110592")
    parser.add_argument("--out", help="write one named frame as an npy archive")
    parser.add_argument("--name", default="Ic-128")
    args = parser.parse_args()
    sizes = tuple(int(part) for part in args.sizes.split(",") if part)
    frames = campaign_frames(sizes)
    for frame in frames:
        box = frame["box9"]
        print(
            f"{frame['name']} molecules={frame['molecules']} "
            f"lx={box[3]:.6f} ly={box[4]:.6f} lz={box[5]:.6f} xy={box[6]:.6f}"
        )
    if args.out:
        chosen = next(frame for frame in frames if frame["name"] == args.name)
        np.savez(args.out, positions=chosen["positions"], lattice=chosen["lattice"], box9=chosen["box9"])


if __name__ == "__main__":
    main()
