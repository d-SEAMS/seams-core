"""LAMMPS dump bounds to the cell edge.

A triclinic dump stores the bounding box, not the cell edges.
``xlo_bound = xlo + min(0, xy, xz, xy+xz)``. The tilt factors sit in
the third column, in the order xy, xz, yz.
"""

from __future__ import annotations


def cell_from_bounds(rows):
    """Three bound rows ``(lo, hi[, tilt])`` to origin, edges, and tilts.

    Origin is ``(xlo, ylo, zlo)``. Edges are ``(lx, ly, lz)`` with
    ``lx = xhi - xlo``. Tilts are ``(xy, xz, yz)``.
    """
    xy = float(rows[0][2]) if len(rows[0]) > 2 else 0.0
    xz = float(rows[1][2]) if len(rows[1]) > 2 else 0.0
    yz = float(rows[2][2]) if len(rows[2]) > 2 else 0.0
    xmin = min(0.0, xy, xz, xy + xz)
    xmax = max(0.0, xy, xz, xy + xz)
    ymin = min(0.0, yz)
    ymax = max(0.0, yz)
    xlo_b, xhi_b = float(rows[0][0]), float(rows[0][1])
    ylo_b, yhi_b = float(rows[1][0]), float(rows[1][1])
    zlo_b, zhi_b = float(rows[2][0]), float(rows[2][1])
    origin = (xlo_b - xmin, ylo_b - ymin, zlo_b)
    edges = ((xhi_b - xlo_b) - xmax + xmin, (yhi_b - ylo_b) - ymax + ymin, zhi_b - zlo_b)
    return origin, edges, (xy, xz, yz)


def read_dump_cell(path, frame=1):
    """Origin, edges, and tilts of a 1-based LAMMPS dump frame."""
    seen = 0
    with open(path) as handle:
        while True:
            line = handle.readline()
            if not line:
                break
            if not line.startswith("ITEM: TIMESTEP"):
                continue
            handle.readline()
            seen += 1
            if handle.readline().strip() != "ITEM: NUMBER OF ATOMS":
                raise ValueError("dump is missing the atom count")
            count = int(handle.readline())
            bounds = handle.readline()
            if not bounds.startswith("ITEM: BOX BOUNDS"):
                raise ValueError("dump is missing box bounds")
            rows = [handle.readline().split() for _ in range(3)]
            for _ in range(count + 1):
                handle.readline()
            if seen == frame:
                return cell_from_bounds(rows)
    raise ValueError(f"{path} has {seen} frames, not {frame}")
