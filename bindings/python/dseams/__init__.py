# Copyright (c) d-SEAMS developers
"""d-SEAMS from the stable C API.

Positions are Angstrom, C-contiguous, and are not copied on the way in.
``box`` is three edge lengths or the nine LAMMPS numbers
``ox, oy, oz, lx, ly, lz, xy, xz, yz``.
"""

from __future__ import annotations

import numpy as np

from dseams import _core

PUBLISHED = _core.PUBLISHED
OVITO = _core.OVITO
CUTOFF = _core.CUTOFF
MUTUAL = _core.MUTUAL

version = _core.version
set_num_threads = _core.set_num_threads
num_threads = _core.num_threads
Frame = _core.Frame


def _as_box(box):
    raw = np.ascontiguousarray(box, dtype=np.float64).reshape(-1)
    if raw.size not in (3, 9):
        raise ValueError("box is 3 lengths or 9 LAMMPS cell numbers")
    return raw


def chill_plus(positions, box, rule=PUBLISHED, graph=CUTOFF, cutoff=3.5):
    """CHILL+ labels and the nine population counts."""
    pos = np.ascontiguousarray(positions, dtype=np.float64)
    frame = Frame()
    frame.set_frame(pos, _as_box(box))
    return frame.chill_plus(cutoff, graph, rule)


def q4q6(positions, box, cutoff=3.5):
    pos = np.ascontiguousarray(positions, dtype=np.float64)
    frame = Frame()
    frame.set_frame(pos, _as_box(box))
    return frame.q4q6(cutoff)


def neighbours(positions, box, cutoff=3.5, graph=CUTOFF):
    """CSR neighbour list, self entries omitted."""
    pos = np.ascontiguousarray(positions, dtype=np.float64)
    frame = Frame()
    frame.set_frame(pos, _as_box(box))
    return frame.neighbours(cutoff, graph)


def rings(positions, box, cutoff=3.5, max_depth=6):
    """Ring sizes, symmetry labels, size histogram, and the six-ring atom fraction."""
    pos = np.ascontiguousarray(positions, dtype=np.float64)
    frame = Frame()
    frame.set_frame(pos, _as_box(box))
    return frame.rings(cutoff, max_depth)


def _reduced_cell(vectors):
    """LAMMPS lower-triangular cell and the rotation that puts coordinates in it."""
    cell = np.asarray(vectors, dtype=np.float64).reshape(3, 3)
    a, b, c = cell
    lx = float(np.linalg.norm(a))
    if lx <= 0.0:
        raise ValueError("cell vector a has zero length")
    ax = a / lx
    xy = float(np.dot(b, ax))
    by = b - xy * ax
    ly = float(np.linalg.norm(by))
    if ly <= 0.0:
        raise ValueError("cell vector b has zero height")
    bx = by / ly
    cx = np.cross(ax, bx)
    if np.dot(cx, c) < 0.0:
        bx = -bx
        xy = -xy
        cx = np.cross(ax, bx)
    xz = float(np.dot(c, ax))
    yz = float(np.dot(c, bx))
    lz = float(np.dot(c, cx))
    rot = np.vstack((ax, bx, cx))
    box = np.array([0.0, 0.0, 0.0, lx, ly, lz, xy, xz, yz], dtype=np.float64)
    return box, rot


def from_ase(atoms):
    """Positions and a LAMMPS cell from an ASE ``Atoms`` object."""
    box, rot = _reduced_cell(atoms.cell.array)
    pos = np.asarray(atoms.get_positions(), dtype=np.float64) @ rot.T
    return np.ascontiguousarray(pos), box


def from_mdanalysis(universe, select="all"):
    """Positions and a LAMMPS cell from an MDAnalysis universe or atom group."""
    group = universe.select_atoms(select) if hasattr(universe, "select_atoms") else universe
    dims = np.asarray(universe.dimensions if hasattr(universe, "dimensions") else group.dimensions,
                      dtype=np.float64)
    if dims.shape[0] < 6:
        raise ValueError("universe.dimensions needs lengths and angles")
    a, b, c, alpha, beta, gamma = (float(v) for v in dims[:6])
    rad = np.pi / 180.0
    alpha, beta, gamma = alpha * rad, beta * rad, gamma * rad
    sin_g = np.sin(gamma)
    if a <= 0.0 or b <= 0.0 or c <= 0.0 or abs(sin_g) < 1e-18:
        raise ValueError("box lengths and gamma must be positive")
    xy = b * np.cos(gamma)
    ly = b * sin_g
    xz = c * np.cos(beta)
    yz = c * (np.cos(alpha) - np.cos(beta) * np.cos(gamma)) / sin_g
    tol = 1e-10 * max(a, b, c)
    if abs(xy) < tol:
        xy = 0.0
    if abs(xz) < tol:
        xz = 0.0
    if abs(yz) < tol:
        yz = 0.0
    lz2 = c * c - xz * xz - yz * yz
    if ly <= 0.0 or lz2 <= 0.0:
        raise ValueError("cell volume is not positive")
    box = np.array([0.0, 0.0, 0.0, a, ly, float(np.sqrt(lz2)), xy, xz, yz], dtype=np.float64)
    pos = np.ascontiguousarray(group.positions, dtype=np.float64)
    return pos, box
