"""Tetrahedral q, LSI, and ring winding against a pure-Python reference.

The reference below shares no code with libseams. It reads the LAMMPS dump
itself, finds neighbours with its own cell list and orthorhombic minimum
image, and evaluates

- Errington--Debenedetti q = 1 - 9/(2N(N-1)) sum_{j<k} (cos psi_jk + 1/3)^2,
- Shiratani--Sasai LSI: sorted distances r_1 < ... < r_N <= cutoff < r_{N+1},
  population variance of the N gaps r_{i+1} - r_i,
- winding: the minimum-image bond steps summed around a ring, in cell
  lengths.

``nucleation_f1.lammpstrj`` is 4096 supercooled liquid water sites;
``ice1h.lammpstrj`` is a 128-molecule ice Ih cell short enough for
primitive rings to wrap.
"""

from __future__ import annotations

import ctypes
import itertools
import math
import os
from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]
WATER = ROOT / "input" / "traj" / "nucleation_f1.lammpstrj"
ICE = ROOT / "input" / "traj" / "ice1h.lammpstrj"

lib = ctypes.CDLL(os.environ["SEAMS_LANG_LIB"])
c_int_p = ctypes.POINTER(ctypes.c_int)
c_double_p = ctypes.POINTER(ctypes.c_double)

lib.seams_frame_from_arrays.argtypes = [
    c_double_p,
    ctypes.c_int,
    c_double_p,
    c_int_p,
    c_int_p,
]
lib.seams_frame_from_arrays.restype = ctypes.c_void_p
lib.seams_frame_free.argtypes = [ctypes.c_void_p]
lib.seams_frame_tetrahedral.argtypes = [
    ctypes.c_void_p,
    c_int_p,
    ctypes.c_int,
    c_int_p,
    ctypes.c_int,
    ctypes.c_int,
    ctypes.c_double,
    c_double_p,
]
lib.seams_frame_tetrahedral.restype = ctypes.c_int
lib.seams_frame_lsi.argtypes = [
    ctypes.c_void_p,
    c_int_p,
    ctypes.c_int,
    c_int_p,
    ctypes.c_int,
    ctypes.c_double,
    c_double_p,
]
lib.seams_frame_lsi.restype = ctypes.c_int
lib.seams_graph_rings.argtypes = [
    c_int_p,
    ctypes.c_int,
    c_int_p,
    ctypes.c_int,
    ctypes.c_int,
    c_int_p,
    ctypes.c_int,
    c_int_p,
    ctypes.c_int,
    c_int_p,
    ctypes.c_int,
    c_int_p,
    c_int_p,
]
lib.seams_graph_rings.restype = ctypes.c_int
lib.seams_frame_ring_winding.argtypes = [
    ctypes.c_void_p,
    c_int_p,
    ctypes.c_int,
    c_int_p,
    ctypes.c_int,
    ctypes.c_int,
    c_int_p,
    c_int_p,
]
lib.seams_frame_ring_winding.restype = ctypes.c_int

TOLERANCE = 1e-12
Q_CUTOFF = 3.5
LSI_CUTOFF = 3.7
SHELL_RADIUS = 6.0


def read_dump(path: Path) -> tuple[list[tuple[float, float, float]], list[float]]:
    lines = path.read_text().splitlines()
    count = int(lines[lines.index("ITEM: NUMBER OF ATOMS") + 1])
    bounds_at = next(i for i, row in enumerate(lines) if row.startswith("ITEM: BOX BOUNDS"))
    bounds = [tuple(map(float, lines[bounds_at + 1 + k].split()[:2])) for k in range(3)]
    atoms_at = next(i for i, row in enumerate(lines) if row.startswith("ITEM: ATOMS"))
    columns = lines[atoms_at].split()[2:]
    ix, iy, iz = (columns.index(name) for name in ("x", "y", "z"))
    points = []
    for row in lines[atoms_at + 1 : atoms_at + 1 + count]:
        fields = row.split()
        points.append((float(fields[ix]), float(fields[iy]), float(fields[iz])))
    origin = [lo for lo, _ in bounds]
    lengths = [hi - lo for lo, hi in bounds]
    return points, origin + lengths


def minimum_image(a, b, lengths):
    step = []
    for k in range(3):
        d = b[k] - a[k]
        d -= lengths[k] * round(d / lengths[k])
        step.append(d)
    return step


def neighbour_shells(points, box, radius, centres):
    """Sorted (distance, unit vector, index) within radius for each centre."""
    origin, lengths = box[:3], box[3:]
    cells = [int(length // radius) for length in lengths]
    assert min(cells) >= 3

    def cell_of(p):
        return tuple(
            int(((p[k] - origin[k]) % lengths[k]) / lengths[k] * cells[k]) % cells[k]
            for k in range(3)
        )

    bins: dict[tuple[int, int, int], list[int]] = {}
    for j, p in enumerate(points):
        bins.setdefault(cell_of(p), []).append(j)
    shells = []
    for i in centres:
        home = cell_of(points[i])
        shell = []
        for offset in itertools.product((-1, 0, 1), repeat=3):
            key = tuple((home[k] + offset[k]) % cells[k] for k in range(3))
            for j in bins.get(key, ()):
                if j == i:
                    continue
                d = minimum_image(points[i], points[j], lengths)
                r = math.sqrt(d[0] * d[0] + d[1] * d[1] + d[2] * d[2])
                if r <= radius:
                    shell.append((r, [d[0] / r, d[1] / r, d[2] / r], j))
        shell.sort(key=lambda entry: entry[0])
        shells.append(shell)
    return shells


def reference_q(shell) -> float:
    n = len(shell)
    if n < 4:
        return math.nan
    total = 0.0
    for j in range(n):
        for k in range(j + 1, n):
            u, v = shell[j][1], shell[k][1]
            total += (u[0] * v[0] + u[1] * v[1] + u[2] * v[2] + 1.0 / 3.0) ** 2
    return 1.0 - 9.0 / (2.0 * n * (n - 1)) * total


def reference_lsi(shell, cutoff) -> float:
    inside = sum(1 for entry in shell if entry[0] <= cutoff)
    if inside == 0 or inside == len(shell):
        return math.nan
    gaps = [shell[i + 1][0] - shell[i][0] for i in range(inside)]
    mean = sum(gaps) / inside
    return sum((gap - mean) ** 2 for gap in gaps) / inside


def make_frame(points, box):
    xyz = (ctypes.c_double * (3 * len(points)))(*[c for p in points for c in p])
    cell = (ctypes.c_double * 9)(*box, 0.0, 0.0, 0.0)
    frame = lib.seams_frame_from_arrays(xyz, len(points), cell, None, None)
    assert frame
    return frame


def same(a: float, b: float) -> bool:
    if math.isnan(a) or math.isnan(b):
        return math.isnan(a) and math.isnan(b)
    return abs(a - b) <= TOLERANCE


def test_liquid_water_tetrahedral_order_and_lsi_match_reference() -> None:
    points, box = read_dump(WATER)
    centres = list(range(0, len(points), 8))
    shells = neighbour_shells(points, box, SHELL_RADIUS, centres)
    # Every shell reaches past the LSI cutoff, and no tie at the fourth
    # neighbour makes the four-nearest shell ambiguous.
    assert all(shell[-1][0] > LSI_CUTOFF for shell in shells)
    assert min(shell[4][0] - shell[3][0] for shell in shells) > 1e-6

    selected = (ctypes.c_int * len(centres))(*centres)
    nearest = (ctypes.c_double * len(centres))()
    within = (ctypes.c_double * len(centres))()
    lsi = (ctypes.c_double * len(centres))()
    frame = make_frame(points, box)
    try:
        assert lib.seams_frame_tetrahedral(frame, selected, len(centres), None, 0, 0, 0.0, nearest) == 0
        assert lib.seams_frame_tetrahedral(frame, selected, len(centres), None, 0, 1, Q_CUTOFF, within) == 0
        assert lib.seams_frame_lsi(frame, selected, len(centres), None, 0, LSI_CUTOFF, lsi) == 0
    finally:
        lib.seams_frame_free(frame)

    ref_nearest = [reference_q(shell[:4]) for shell in shells]
    ref_within = [reference_q([e for e in shell if e[0] <= Q_CUTOFF]) for shell in shells]
    ref_lsi = [reference_lsi(shell, LSI_CUTOFF) for shell in shells]

    assert all(same(a, b) for a, b in zip(nearest, ref_nearest))
    assert all(same(a, b) for a, b in zip(within, ref_within))
    assert all(same(a, b) for a, b in zip(lsi, ref_lsi))
    assert all(math.isfinite(value) for value in nearest)
    assert all(math.isfinite(value) for value in lsi)
    assert sum(math.isfinite(value) for value in within) > len(centres) // 2
    # A tetrahedral liquid lies between the ideal gas (0) and ice (1).
    assert 0.3 < sum(nearest) / len(centres) < 0.95


def bond_graph(points, box, cutoff):
    shells = neighbour_shells(points, box, cutoff, range(len(points)))
    offsets = [0]
    columns = []
    for shell in shells:
        columns.extend(sorted(entry[2] for entry in shell))
        offsets.append(len(columns))
    return offsets, columns


def library_rings(offsets, columns, max_depth):
    off = (ctypes.c_int * len(offsets))(*offsets)
    col = (ctypes.c_int * len(columns))(*columns)
    histogram = (ctypes.c_int * (max_depth + 1))()
    n_rings = ctypes.c_int()
    n_vertices = ctypes.c_int()
    args = (off, len(offsets), col, len(columns), max_depth, histogram, len(histogram))
    assert lib.seams_graph_rings(*args, None, 0, None, 0, ctypes.byref(n_rings), ctypes.byref(n_vertices)) == 0
    ring_offsets = (ctypes.c_int * (n_rings.value + 1))()
    vertices = (ctypes.c_int * n_vertices.value)()
    assert (
        lib.seams_graph_rings(
            *args,
            ring_offsets,
            len(ring_offsets),
            vertices,
            len(vertices),
            ctypes.byref(n_rings),
            ctypes.byref(n_vertices),
        )
        == 0
    )
    return ring_offsets, vertices, n_rings.value


def reference_winding(points, lengths, ring):
    total = [0.0, 0.0, 0.0]
    for i, a in enumerate(ring):
        b = ring[(i + 1) % len(ring)]
        step = minimum_image(points[a], points[b], lengths)
        for k in range(3):
            total[k] += step[k]
    winding = [round(total[k] / lengths[k]) for k in range(3)]
    for k in range(3):
        assert abs(total[k] - winding[k] * lengths[k]) < 1e-9
    return winding


def library_windings(frame, ring_offsets, vertices, n_rings):
    per_axis = []
    contractible = (ctypes.c_int * n_rings)()
    for axis in range(3):
        winding = (ctypes.c_int * n_rings)()
        rc = lib.seams_frame_ring_winding(
            frame, ring_offsets, n_rings, vertices, len(vertices), axis, winding, contractible
        )
        assert rc == 0
        per_axis.append(list(winding))
    return [[per_axis[k][r] for k in range(3)] for r in range(n_rings)], list(contractible)


def check_windings(path: Path, max_depth: int) -> list[list[int]]:
    points, box = read_dump(path)
    offsets, columns = bond_graph(points, box, Q_CUTOFF)
    ring_offsets, vertices, n_rings = library_rings(offsets, columns, max_depth)
    assert n_rings > 0
    frame = make_frame(points, box)
    try:
        windings, contractible = library_windings(frame, ring_offsets, vertices, n_rings)
    finally:
        lib.seams_frame_free(frame)
    rings = [list(vertices[ring_offsets[r] : ring_offsets[r + 1]]) for r in range(n_rings)]
    expected = [reference_winding(points, box[3:], ring) for ring in rings]
    assert windings == expected
    assert contractible == [int(w == [0, 0, 0]) for w in expected]
    return expected


def test_liquid_water_ring_winding_matches_reference() -> None:
    check_windings(WATER, 7)


def test_wrapping_ice_ring_winding_matches_reference() -> None:
    expected = check_windings(ICE, 10)
    wrapping = [w for w in expected if w != [0, 0, 0]]
    assert wrapping
    assert any(component < 0 for w in wrapping for component in w)
    assert any(component > 0 for w in wrapping for component in w)
