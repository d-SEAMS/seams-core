"""Tilted ice Ih and replicated sH cells.

Ice Ih is the orthorhombic generator with a LAMMPS tilt (xy, xz, yz).
sH is the genice triclinic cell, replicated in fractional space. 136 does
not divide 4096, 32768 or 110592, so the supercell is the nearest
nx*ny*nz multiple.
"""

from __future__ import annotations

from pathlib import Path

import numpy as np

from mixed_phase_frames import ice_ih

ROOT = Path(__file__).resolve().parents[1]
SII_PATH = ROOT / "input" / "traj" / "genice_sH.lammpstrj"


def lammps_matrix(lx, ly, lz, xy, xz, yz, origin=(0.0, 0.0, 0.0)):
    """Columns are the cell vectors a, b, c. Origin is separate."""
    matrix = np.array(
        [[lx, xy, xz], [0.0, ly, yz], [0.0, 0.0, lz]],
        dtype=np.float64,
    )
    return matrix, np.asarray(origin, dtype=np.float64)


def box9_of(lx, ly, lz, xy, xz, yz, origin=(0.0, 0.0, 0.0)):
    ox, oy, oz = origin
    return (ox, oy, oz, lx, ly, lz, xy, xz, yz)


def tilted_ice(nx: int, ny: int, nz: int, xy_frac=0.15, xz_frac=0.05, yz_frac=0.08):
    pos, (lx, ly, lz) = ice_ih(nx, ny, nz)
    sx = pos[:, 0] / lx
    sy = pos[:, 1] / ly
    sz = pos[:, 2] / lz
    xy = xy_frac * ly
    xz = xz_frac * lz
    yz = yz_frac * lz
    cart = np.column_stack((lx * sx + xy * sy + xz * sz, ly * sy + yz * sz, lz * sz))
    return np.ascontiguousarray(cart, dtype=np.float64), box9_of(lx, ly, lz, xy, xz, yz)


def _parse_sH():
    text = SII_PATH.read_text().splitlines()
    n = int(text[3])
    bounds = []
    tilts = []
    for line in text[5:8]:
        parts = line.split()
        lo, hi = float(parts[0]), float(parts[1])
        tilt = float(parts[2]) if len(parts) > 2 else 0.0
        bounds.append((lo, hi))
        tilts.append(tilt)
    xy, xz, yz = tilts
    xlo_b, xhi_b = bounds[0]
    ylo_b, yhi_b = bounds[1]
    zlo_b, zhi_b = bounds[2]
    xmin = min(0.0, xy, xz, xy + xz)
    xmax = max(0.0, xy, xz, xy + xz)
    ymin = min(0.0, yz)
    ymax = max(0.0, yz)
    xlo = xlo_b - xmin
    ylo = ylo_b - ymin
    zlo = zlo_b
    lx = (xhi_b - xmax) - xlo
    ly = (yhi_b - ymax) - ylo
    lz = zhi_b - zlo
    coords = []
    for line in text[9 : 9 + n]:
        parts = line.split()
        coords.append((float(parts[-3]), float(parts[-2]), float(parts[-1])))
    pos = np.asarray(coords, dtype=np.float64)
    return pos, (xlo, ylo, zlo), (lx, ly, lz), (xy, xz, yz)


def replicate_sH(nx: int, ny: int, nz: int):
    pos, origin, (lx, ly, lz), (xy, xz, yz) = _parse_sH()
    matrix, origin = lammps_matrix(lx, ly, lz, xy, xz, yz, origin)
    inv = np.linalg.inv(matrix)
    frac = (pos - origin) @ inv.T
    frac -= np.floor(frac)
    blocks = []
    for ix in range(nx):
        for iy in range(ny):
            for iz in range(nz):
                shifted = (frac + np.array([ix, iy, iz])) / np.array([nx, ny, nz])
                blocks.append(shifted)
    all_frac = np.vstack(blocks)
    lx2, ly2, lz2 = nx * lx, ny * ly, nz * lz
    xy2, xz2, yz2 = ny * xy, nz * xz, nz * yz
    big, origin2 = lammps_matrix(lx2, ly2, lz2, xy2, xz2, yz2, (0.0, 0.0, 0.0))
    cart = all_frac @ big.T
    return np.ascontiguousarray(cart, dtype=np.float64), box9_of(lx2, ly2, lz2, xy2, xz2, yz2)


def nearest_sH(target: int):
    """Supercell nearest the target count, preferring a compact cell.

    136 does not divide 4096, 32768 or 110592. Among replications with the
    same molecule count, the one with the least elongated edges is kept.
    A 1x1xN rod makes the cell list do extra image work.
    """
    unit = 136
    _pos, _origin, (lx, ly, lz), _tilt = _parse_sH()
    best = None
    goal = target / unit
    limit = max(4, int(round(goal ** (1.0 / 3.0))) + 8)
    for nx in range(1, limit + 1):
        for ny in range(1, limit + 1):
            for nz in range(1, limit + 1):
                n = unit * nx * ny * nz
                err = abs(n - target)
                # Stay near the requested count. A rod with a smaller
                # count error is a worse cell than a compact one a
                # fraction of a percent away.
                if err > max(136, int(0.02 * target)):
                    continue
                edges = (nx * lx, ny * ly, nz * lz)
                aspect = max(edges) / min(edges)
                rank = (aspect, err)
                if best is None or rank < best[0]:
                    best = (rank, nx, ny, nz, n)
    _rank, nx, ny, nz, n = best
    pos, box = replicate_sH(nx, ny, nz)
    return pos, box, (nx, ny, nz, n)
