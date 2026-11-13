#!/usr/bin/env python3
"""C API against freud and OVITO.

The six frames are ice Ih, ice Ic, sI, the TIP4P/ice liquid, the brine
slab, and the nucleation frame. Neighbour rows, local q4 and q6, the
Lechner-Dellago q6-bar, CHILL+ labels, diamond labels, the RDF, and the
coordination are read from libseams through ctypes.

A missing freud or ovito import is a failed run. The numeric match is
1e-6 where both sides are float64. A looser residual is recorded only
where the other library's own storage or cutoff is the cause, and the
analytic block names which convention that is.
"""

from __future__ import annotations

import argparse
import ctypes
import math
import os
import sys
from pathlib import Path

os.environ.setdefault("OMP_NUM_THREADS", "1")
os.environ.setdefault("RAYON_NUM_THREADS", "1")
os.environ.setdefault("OVITO_THREAD_COUNT", "1")

try:
    import freud
    import numpy as np
    from ovito.data import DataCollection, NearestNeighborFinder, Particles, SimulationCell
    from ovito.modifiers import (
        ChillPlusModifier,
        CoordinationAnalysisModifier,
        CreateBondsModifier,
        PolyhedralTemplateMatchingModifier,
        SteinhardtParametersModifier,
    )
    from ovito.pipeline import Pipeline, StaticSource
    from scipy.special import sph_harm_y
except ImportError as exc:
    sys.stderr.write(f"freud or ovito is not installed ({exc})\n")
    raise SystemExit(1)

CUTOFF = 3.5
K_NEIGH = 4
ABS_TOL = 1e-6
# freud 3.6 stores q_l in float32. On a disordered shell the residual
# against the float64 sum sits near 2e-5.
FREUD_Q_TOL = 3e-5
RDF_BINS = 40
CHILL_FROM_OVITO = {2: 0, 1: 1, 0: 2, 3: 3, 4: 4, 5: 5}
SYM_TO_PTM = {0: 7, 1: 6}
PTM_TO_SYM = {7: 0, 6: 1}

FAILURES: list[str] = []


def fail(message: str) -> None:
    FAILURES.append(message)
    print(f"fail: {message}", file=sys.stderr)


def expect(cond: bool, message: str) -> None:
    if not cond:
        fail(message)


def load_library() -> ctypes.CDLL:
    path = os.environ.get("SEAMS_LIB", "")
    if not path:
        fail("SEAMS_LIB is not set")
        raise SystemExit(1)
    lib = ctypes.CDLL(path)
    lib.seams_system_create.restype = ctypes.c_void_p
    lib.seams_system_destroy.argtypes = [ctypes.c_void_p]
    lib.seams_system_set_frame.restype = ctypes.c_int
    lib.seams_system_set_frame.argtypes = [
        ctypes.c_void_p,
        ctypes.POINTER(ctypes.c_double),
        ctypes.c_int,
        ctypes.POINTER(ctypes.c_double),
        ctypes.POINTER(ctypes.c_int),
        ctypes.POINTER(ctypes.c_int),
    ]
    lib.seams_last_error.restype = ctypes.c_char_p
    lib.seams_set_num_threads.restype = ctypes.c_int
    lib.seams_set_num_threads.argtypes = [ctypes.c_int]
    lib.seams_ql.restype = ctypes.c_int
    lib.seams_ql.argtypes = [
        ctypes.c_void_p,
        ctypes.c_double,
        ctypes.c_int,
        ctypes.c_int,
        ctypes.c_int,
        ctypes.POINTER(ctypes.c_double),
        ctypes.c_size_t,
        ctypes.POINTER(ctypes.c_double),
        ctypes.c_size_t,
    ]
    lib.seams_q6bar.restype = ctypes.c_int
    lib.seams_q6bar.argtypes = [
        ctypes.c_void_p,
        ctypes.c_double,
        ctypes.c_int,
        ctypes.c_double,
        ctypes.POINTER(ctypes.c_double),
        ctypes.c_size_t,
        ctypes.POINTER(ctypes.c_double),
        ctypes.POINTER(ctypes.c_double),
        ctypes.POINTER(ctypes.c_int),
    ]
    lib.seams_neighbours_build.restype = ctypes.c_int
    lib.seams_neighbours_build.argtypes = [
        ctypes.c_void_p,
        ctypes.c_double,
        ctypes.c_int,
        ctypes.c_int,
        ctypes.POINTER(ctypes.c_void_p),
    ]
    lib.seams_neighbours_knn.restype = ctypes.c_int
    lib.seams_neighbours_knn.argtypes = [
        ctypes.c_void_p,
        ctypes.c_int,
        ctypes.c_double,
        ctypes.c_int,
        ctypes.c_int,
        ctypes.POINTER(ctypes.c_void_p),
    ]
    lib.seams_neighbours_copy.restype = ctypes.c_int
    lib.seams_neighbours_copy.argtypes = [
        ctypes.c_void_p,
        ctypes.POINTER(ctypes.c_int),
        ctypes.c_size_t,
        ctypes.POINTER(ctypes.c_size_t),
        ctypes.POINTER(ctypes.c_int),
        ctypes.c_size_t,
        ctypes.POINTER(ctypes.c_size_t),
    ]
    lib.seams_neighbours_destroy.argtypes = [ctypes.c_void_p]
    lib.seams_chill_plus.restype = ctypes.c_int
    lib.seams_chill_plus.argtypes = [
        ctypes.c_void_p,
        ctypes.c_double,
        ctypes.c_int,
        ctypes.c_int,
        ctypes.c_int,
        ctypes.POINTER(ctypes.c_int),
        ctypes.POINTER(ctypes.c_int),
    ]
    lib.seams_system_symmetry.restype = ctypes.c_int
    lib.seams_system_symmetry.argtypes = [
        ctypes.c_void_p,
        ctypes.c_double,
        ctypes.POINTER(ctypes.c_int),
        ctypes.c_size_t,
    ]
    lib.seams_rdf.restype = ctypes.c_int
    lib.seams_rdf.argtypes = [
        ctypes.c_void_p,
        ctypes.c_double,
        ctypes.c_int,
        ctypes.c_int,
        ctypes.c_int,
        ctypes.POINTER(ctypes.c_double),
        ctypes.POINTER(ctypes.c_double),
        ctypes.POINTER(ctypes.c_int),
        ctypes.c_size_t,
        ctypes.POINTER(ctypes.c_size_t),
        ctypes.POINTER(ctypes.c_double),
    ]
    lib.seams_coordination.restype = ctypes.c_int
    lib.seams_coordination.argtypes = [
        ctypes.c_void_p,
        ctypes.c_double,
        ctypes.c_int,
        ctypes.c_int,
        ctypes.c_int,
        ctypes.POINTER(ctypes.c_double),
    ]
    if lib.seams_set_num_threads(1) != 0:
        fail("seams_set_num_threads")
        raise SystemExit(1)
    freud.parallel.set_num_threads(1)
    return lib


def status_ok(lib: ctypes.CDLL, code: int, what: str) -> bool:
    if code == 0:
        return True
    detail = lib.seams_last_error()
    text = detail.decode() if detail else ""
    fail(f"{what} status {code} {text}")
    return False


def open_system(lib: ctypes.CDLL, pts: np.ndarray, box: np.ndarray):
    n = len(pts)
    flat = np.ascontiguousarray(pts.reshape(-1), dtype=np.float64)
    cell = (ctypes.c_double * 9)(*[float(v) for v in box])
    handle = lib.seams_system_create()
    code = lib.seams_system_set_frame(
        handle,
        flat.ctypes.data_as(ctypes.POINTER(ctypes.c_double)),
        n,
        cell,
        None,
        None,
    )
    if not status_ok(lib, code, "set_frame"):
        lib.seams_system_destroy(handle)
        return None
    return handle


def copy_csr(lib: ctypes.CDLL, handle) -> tuple[np.ndarray, np.ndarray]:
    n_off = ctypes.c_size_t()
    nnz = ctypes.c_size_t()
    lib.seams_neighbours_copy(handle, None, 0, ctypes.byref(n_off), None, 0, ctypes.byref(nnz))
    off = (ctypes.c_int * n_off.value)()
    cols = (ctypes.c_int * max(int(nnz.value), 1))()
    lib.seams_neighbours_copy(
        handle,
        off,
        n_off.value,
        ctypes.byref(n_off),
        cols,
        max(int(nnz.value), 1),
        ctypes.byref(nnz),
    )
    return np.ctypeslib.as_array(off).copy(), np.ctypeslib.as_array(cols)[: int(nnz.value)].copy()


def rows_of(off: np.ndarray, cols: np.ndarray) -> list[tuple[int, ...]]:
    return [tuple(sorted(int(v) for v in cols[off[i] : off[i + 1]])) for i in range(len(off) - 1)]


def freud_box(box: np.ndarray) -> freud.box.Box:
    return freud.box.Box(
        Lx=float(box[3]),
        Ly=float(box[4]),
        Lz=float(box[5]),
        xy=float(box[6]),
        xz=float(box[7]),
        yz=float(box[8]),
    )


def ovito_data(pts: np.ndarray, box: np.ndarray) -> DataCollection:
    data = DataCollection()
    cell = SimulationCell(pbc=(True, True, True))
    lx, ly, lz = [float(v) for v in box[3:6]]
    xy, xz, yz = [float(v) for v in box[6:9]]
    cell[:, :3] = np.array([[lx, xy, xz], [0.0, ly, yz], [0.0, 0.0, lz]], dtype=float)
    cell[:, 3] = np.asarray(box[:3], dtype=float)
    data.objects.append(cell)
    part = Particles()
    part.create_property("Position", data=np.ascontiguousarray(pts, dtype=np.float64))
    data.objects.append(part)
    return data


def freud_directed(pts: np.ndarray, box: np.ndarray, query: dict) -> list[list[int]]:
    cloud = freud.locality.AABBQuery(freud_box(box), pts)
    found = cloud.query(pts, query).toNeighborList()
    rows: list[list[int]] = [[] for _ in range(len(pts))]
    for i, j in zip(found.query_point_indices, found.point_indices):
        rows[int(i)].append(int(j))
    return rows


def close_directed(directed: list[list[int]], mutual: bool) -> list[tuple[int, ...]]:
    sets = [set(row) for row in directed]
    closed = [set() for _ in sets]
    for i, row in enumerate(sets):
        for j in row:
            if mutual and i not in sets[j]:
                continue
            closed[i].add(j)
            closed[j].add(i)
    return [tuple(sorted(row)) for row in closed]


def ovito_bonds(pts: np.ndarray, box: np.ndarray, cutoff: float) -> list[tuple[int, ...]]:
    pipe = Pipeline(source=StaticSource(data=ovito_data(pts, box)))
    pipe.modifiers.append(CreateBondsModifier(mode=CreateBondsModifier.Mode.Uniform, cutoff=cutoff))
    out = pipe.compute()
    bonds = out.particles.bonds
    rows: list[set[int]] = [set() for _ in range(len(pts))]
    if bonds is not None and bonds.count:
        for a, b in np.asarray(bonds.topology):
            rows[int(a)].add(int(b))
            rows[int(b)].add(int(a))
    return [tuple(sorted(row)) for row in rows]


def ovito_knn(pts: np.ndarray, box: np.ndarray, k: int) -> list[list[int]]:
    finder = NearestNeighborFinder(k, ovito_data(pts, box))
    return [[int(neigh.index) for neigh in finder.find(i)] for i in range(len(pts))]


def steinhardt_column(pts: np.ndarray, box: np.ndarray, order: int, averaged: bool) -> np.ndarray:
    pipe = Pipeline(source=StaticSource(data=ovito_data(pts, box)))
    mod = SteinhardtParametersModifier(cutoff=CUTOFF, max_l=order, compute_averaged=averaged)
    mod.neighbor_mode = SteinhardtParametersModifier.NeighborMode.Cutoff
    pipe.modifiers.append(mod)
    out = pipe.compute()
    key = "Smooth Steinhardt" if averaged else "Steinhardt"
    column = np.asarray(out.particles[key], dtype=np.float64)
    return np.ascontiguousarray(column[:, -1])


def chill_labels(pts: np.ndarray, box: np.ndarray) -> np.ndarray:
    pipe = Pipeline(source=StaticSource(data=ovito_data(pts, box)))
    pipe.modifiers.append(ChillPlusModifier(cutoff=CUTOFF))
    kinds = np.asarray(pipe.compute().particles["Structure Type"])
    return np.array([CHILL_FROM_OVITO[int(v)] for v in kinds], dtype=np.int32)


def ptm_labels(pts: np.ndarray, box: np.ndarray) -> tuple[np.ndarray, np.ndarray]:
    mod = PolyhedralTemplateMatchingModifier(rmsd_cutoff=0.1, output_rmsd=True)
    for structure in mod.structures:
        if int(structure.id) in (6, 7):
            structure.enabled = True
    pipe = Pipeline(source=StaticSource(data=ovito_data(pts, box)))
    pipe.modifiers.append(mod)
    out = pipe.compute()
    return (
        np.asarray(out.particles["Structure Type"]).copy(),
        np.asarray(out.particles["RMSD"], dtype=np.float64).copy(),
    )


def max_abs(left: np.ndarray, right: np.ndarray) -> float:
    other = np.asarray(right, dtype=np.float64)
    delta = np.abs(np.asarray(left, dtype=np.float64) - other)
    live = np.isfinite(other)
    if not np.any(live):
        return 0.0
    return float(np.max(delta[live]))


def mic_rows(pts: np.ndarray, box: np.ndarray, index: int, neigh: list[int]) -> np.ndarray:
    lengths = np.asarray(box[3:6], dtype=np.float64)
    delta = pts[neigh] - pts[index]
    delta -= lengths * np.round(delta / lengths)
    return delta


def qlm_mean(vecs: np.ndarray, order: int) -> np.ndarray:
    radius = np.linalg.norm(vecs, axis=1)
    theta = np.arccos(np.clip(vecs[:, 2] / radius, -1.0, 1.0))
    phi = np.arctan2(vecs[:, 1], vecs[:, 0])
    means = []
    for moment in range(-order, order + 1):
        means.append(np.mean(sph_harm_y(order, moment, theta, phi)))
    return np.array(means)


def q_from_qlm(qlm: np.ndarray, order: int) -> float:
    pref = 4.0 * math.pi / (2.0 * order + 1.0)
    return math.sqrt(pref * float(np.sum(np.abs(qlm) ** 2)))


def shell_volume(bin_index: int, width: float) -> float:
    r0 = width * bin_index
    r1 = width * (bin_index + 1)
    return 4.0 * math.pi * (r1**3 - r0**3) / 3.0


def integrate_cn(g_values: np.ndarray, width: float, rho: float, r_cut: float) -> float:
    total = 0.0
    for index, value in enumerate(g_values):
        r0 = width * index
        if r0 >= r_cut:
            break
        r1 = min(width * (index + 1), r_cut)
        total += 4.0 * math.pi * rho * float(value) * (r1**3 - r0**3) / 3.0
    return total


def load_dump(path: Path) -> tuple[np.ndarray, np.ndarray]:
    lines = path.read_text(encoding="ascii").splitlines()
    index = 0
    count = None
    atoms: list[tuple[float, float, float]] = []
    lo = hi = None
    while index < len(lines):
        line = lines[index]
        if line.startswith("ITEM: NUMBER"):
            count = int(lines[index + 1])
            index += 2
            continue
        if line.startswith("ITEM: BOX"):
            bounds = []
            for step in range(3):
                bounds.append([float(part) for part in lines[index + 1 + step].split()])
            lo = [bounds[axis][0] for axis in range(3)]
            hi = [bounds[axis][1] for axis in range(3)]
            index += 4
            continue
        if line.startswith("ITEM: ATOMS"):
            index += 1
            while count is not None and len(atoms) < count and index < len(lines):
                parts = lines[index].split()
                index += 1
                if len(parts) < 5:
                    continue
                atoms.append((float(parts[2]), float(parts[3]), float(parts[4])))
            break
        index += 1
    if count is None or lo is None or hi is None or len(atoms) != count:
        raise RuntimeError(f"could not read {path}")
    pts = np.array(atoms, dtype=np.float64)
    box = np.array(
        [lo[0], lo[1], lo[2], hi[0] - lo[0], hi[1] - lo[1], hi[2] - lo[2], 0.0, 0.0, 0.0],
        dtype=np.float64,
    )
    return pts, box


def load_oxygen_text(path: Path) -> tuple[np.ndarray, np.ndarray]:
    text = path.read_text(encoding="ascii").split()
    count = int(text[0])
    lengths = [float(text[1]), float(text[2]), float(text[3])]
    body = np.array([float(value) for value in text[4 : 4 + 3 * count]], dtype=np.float64)
    pts = body.reshape(count, 3)
    box = np.array([0.0, 0.0, 0.0, *lengths, 0.0, 0.0, 0.0], dtype=np.float64)
    return pts, box


def ice_ic() -> tuple[np.ndarray, np.ndarray]:
    repeat = 2
    bond = 2.75
    lattice = 4.0 * bond / math.sqrt(3.0)
    basis = np.array(
        [
            [0.0, 0.0, 0.0],
            [0.5, 0.5, 0.0],
            [0.5, 0.0, 0.5],
            [0.0, 0.5, 0.5],
            [0.25, 0.25, 0.25],
            [0.75, 0.75, 0.25],
            [0.75, 0.25, 0.75],
            [0.25, 0.75, 0.75],
        ]
    )
    grid = np.stack(np.meshgrid(*[np.arange(repeat)] * 3, indexing="ij"), -1).reshape(-1, 3)
    pts = (grid[:, None, :] + basis).reshape(-1, 3) * lattice
    side = repeat * lattice
    box = np.array([0.0, 0.0, 0.0, side, side, side, 0.0, 0.0, 0.0], dtype=np.float64)
    return np.ascontiguousarray(pts, dtype=np.float64), box


def frames_of(root: Path) -> list[tuple[str, np.ndarray, np.ndarray]]:
    return [
        ("ice_ih", *load_dump(root / "input/traj/ice1h.lammpstrj")),
        ("ice_ic", *ice_ic()),
        ("si", *load_dump(root / "input/traj/genice_sI.lammpstrj")),
        ("liquid", *load_oxygen_text(root / "benchmarks/paper-cells/tip4pice-liquid-4096.txt")),
        ("brine", *load_dump(root / "input/traj/nacl_front.lammpstrj")),
        ("nucleation", *load_dump(root / "input/traj/nucleation_f1.lammpstrj")),
    ]


def seams_rows(lib: ctypes.CDLL, handle, kind: str, k: int = K_NEIGH) -> list[tuple[int, ...]]:
    owned = ctypes.c_void_p()
    if kind == "cutoff":
        code = lib.seams_neighbours_build(handle, CUTOFF, 0, 0, ctypes.byref(owned))
    elif kind == "knn":
        code = lib.seams_neighbours_knn(handle, k, CUTOFF, 0, 0, ctypes.byref(owned))
    else:
        code = lib.seams_neighbours_knn(handle, k, CUTOFF, 0, 1, ctypes.byref(owned))
    if not status_ok(lib, code, kind):
        return []
    off, cols = copy_csr(lib, owned)
    lib.seams_neighbours_destroy(owned)
    return rows_of(off, cols)


def check_ptm(name: str, sym: np.ndarray, kinds: np.ndarray, rmsd: np.ndarray) -> int:
    explained = 0
    for index in range(len(sym)):
        kernel = int(sym[index])
        other = int(kinds[index])
        distance = float(rmsd[index])
        kernel_ptm = SYM_TO_PTM.get(kernel)
        if kernel_ptm == other:
            continue
        if kernel_ptm is None and other not in PTM_TO_SYM:
            continue
        if kernel_ptm is not None and other not in PTM_TO_SYM and distance > 0.1:
            explained += 1
            continue
        if other in PTM_TO_SYM and kernel != PTM_TO_SYM[other] and distance < 0.1:
            explained += 1
            continue
        fail(f"{name} ptm atom {index} kernel {kernel} ovito {other} rmsd {distance:.6g}")
        return explained
    if name in ("ice_ih", "ice_ic"):
        expect(explained == 0, f"{name} ideal crystal disagrees with PTM ({explained})")
    # The generated ice Ic lattice sits on the diamond template. The ice Ih
    # file is the same label with a finite template RMSD, still under 0.1.
    if name == "ice_ic":
        expect(float(np.max(rmsd)) < ABS_TOL, f"{name} PTM rmsd {float(np.max(rmsd)):.3e}")
    if name == "ice_ih":
        expect(float(np.max(rmsd)) < 0.1, f"{name} PTM rmsd {float(np.max(rmsd)):.3e}")
    return explained


def compare_frame(lib: ctypes.CDLL, name: str, pts: np.ndarray, box: np.ndarray) -> dict:
    n = len(pts)
    handle = open_system(lib, pts, box)
    record = {"frame": name, "n": n}
    if handle is None:
        return record
    q4 = (ctypes.c_double * n)()
    q6 = (ctypes.c_double * n)()
    q6bar = (ctypes.c_double * n)()
    q6bar_api = (ctypes.c_double * n)()
    ok = status_ok(lib, lib.seams_ql(handle, CUTOFF, 0, 4, 0, q4, n, None, 0), "q4")
    ok = status_ok(lib, lib.seams_ql(handle, CUTOFF, 0, 6, 0, q6, n, q6bar, n), "q6") and ok
    mean = ctypes.c_double()
    median = ctypes.c_double()
    n_ice = ctypes.c_int()
    ok = (
        status_ok(
            lib,
            lib.seams_q6bar(
                handle,
                CUTOFF,
                0,
                0.36,
                q6bar_api,
                n,
                ctypes.byref(mean),
                ctypes.byref(median),
                ctypes.byref(n_ice),
            ),
            "q6bar",
        )
        and ok
    )
    if not ok:
        lib.seams_system_destroy(handle)
        return record
    q4a = np.ctypeslib.as_array(q4).copy()
    q6a = np.ctypeslib.as_array(q6).copy()
    q6ba = np.ctypeslib.as_array(q6bar).copy()
    q6b2 = np.ctypeslib.as_array(q6bar_api).copy()
    expect(max_abs(q6ba, q6b2) <= ABS_TOL, f"{name} q6bar api {max_abs(q6ba, q6b2):.3e}")

    cutoff_rows = seams_rows(lib, handle, "cutoff")
    knn_rows = seams_rows(lib, handle, "knn")
    mutual_rows = seams_rows(lib, handle, "mutual")
    freud_ball = freud_directed(pts, box, {"mode": "ball", "r_max": CUTOFF, "exclude_ii": True})
    freud_near = freud_directed(
        pts, box, {"mode": "nearest", "num_neighbors": K_NEIGH, "exclude_ii": True}
    )
    ovito_near = ovito_knn(pts, box, K_NEIGH)
    ball_rows = [tuple(sorted(row)) for row in freud_ball]
    record["cutoff"] = sum(left != right for left, right in zip(cutoff_rows, ball_rows))
    record["cutoff_ovito"] = sum(
        left != right for left, right in zip(cutoff_rows, ovito_bonds(pts, box, CUTOFF))
    )
    record["knn"] = sum(
        left != right for left, right in zip(knn_rows, close_directed(freud_near, False))
    )
    record["knn_ovito"] = sum(
        left != right for left, right in zip(knn_rows, close_directed(ovito_near, False))
    )
    record["mutual"] = sum(
        left != right for left, right in zip(mutual_rows, close_directed(freud_near, True))
    )
    record["mutual_ovito"] = sum(
        left != right for left, right in zip(mutual_rows, close_directed(ovito_near, True))
    )
    for key in ("cutoff", "cutoff_ovito", "knn", "knn_ovito", "mutual", "mutual_ovito"):
        expect(record[key] == 0, f"{name} {key} mismatches {record[key]}")

    order4 = freud.order.Steinhardt(l=4)
    order6 = freud.order.Steinhardt(l=6)
    order6b = freud.order.Steinhardt(l=6, average=True)
    neigh = {"r_max": CUTOFF, "exclude_ii": True}
    order4.compute((freud_box(box), pts), neighbors=neigh)
    order6.compute((freud_box(box), pts), neighbors=neigh)
    order6b.compute((freud_box(box), pts), neighbors=neigh)
    record["q4_ovito"] = max_abs(q4a, steinhardt_column(pts, box, 4, False))
    record["q6_ovito"] = max_abs(q6a, steinhardt_column(pts, box, 6, False))
    record["q6bar_ovito"] = max_abs(q6ba, steinhardt_column(pts, box, 6, True))
    record["q4_freud"] = max_abs(q4a, order4.ql)
    record["q6_freud"] = max_abs(q6a, order6.ql)
    record["q6bar_freud"] = max_abs(q6ba, order6b.ql)
    for key in ("q4_ovito", "q6_ovito", "q6bar_ovito"):
        expect(record[key] <= ABS_TOL, f"{name} {key} {record[key]:.3e}")
    for key in ("q4_freud", "q6_freud", "q6bar_freud"):
        expect(record[key] <= FREUD_Q_TOL, f"{name} {key} {record[key]:.3e}")
    empty = ~np.isfinite(np.asarray(order6.ql, dtype=np.float64))
    if np.any(empty):
        expect(float(np.max(np.abs(q6a[empty]))) == 0.0, f"{name} empty shell is not 0")

    labels = (ctypes.c_int * n)()
    counts = (ctypes.c_int * 9)()
    if status_ok(
        lib,
        lib.seams_chill_plus(handle, CUTOFF, 0, 0, 1, labels, counts),
        "chill+",
    ):
        got = np.ctypeslib.as_array(labels).copy()
        record["chill"] = int(np.sum(got != chill_labels(pts, box)))
        expect(record["chill"] == 0, f"{name} chill+ mismatches {record['chill']}")

    sym = (ctypes.c_int * n)()
    if status_ok(lib, lib.seams_system_symmetry(handle, 0.15, sym, n), "symmetry"):
        kernel = np.ctypeslib.as_array(sym).copy()
        kinds, rmsd = ptm_labels(pts, box)
        record["ptm_split"] = check_ptm(name, kernel, kinds, rmsd)
        if name == "ice_ih":
            expect(int(np.all(kernel == 0)) == 1, f"{name} kernel is not hexagonal")
        if name == "ice_ic":
            expect(int(np.all(kernel == 1)) == 1, f"{name} kernel is not cubic")

    r_max = 0.45 * float(np.min(box[3:6]))
    radius = (ctypes.c_double * RDF_BINS)()
    gram = (ctypes.c_double * RDF_BINS)()
    hist = (ctypes.c_int * RDF_BINS)()
    n_out = ctypes.c_size_t()
    volume = ctypes.c_double()
    if status_ok(
        lib,
        lib.seams_rdf(
            handle,
            r_max,
            RDF_BINS,
            0,
            0,
            radius,
            gram,
            hist,
            RDF_BINS,
            ctypes.byref(n_out),
            ctypes.byref(volume),
        ),
        "rdf",
    ):
        g_seams = np.ctypeslib.as_array(gram).copy()
        r_seams = np.ctypeslib.as_array(radius).copy()
        counts = np.ctypeslib.as_array(hist).copy()
        scale = (n - 1) / n
        pipe = Pipeline(source=StaticSource(data=ovito_data(pts, box)))
        pipe.modifiers.append(
            CoordinationAnalysisModifier(cutoff=r_max, number_of_bins=RDF_BINS)
        )
        table = pipe.compute().tables["coordination-rdf"]
        g_ovito = np.asarray(table["g(r)"], dtype=np.float64)
        r_ovito = np.asarray(table.xy(), dtype=np.float64)[:, 0]
        record["rdf_ovito"] = float(np.max(np.abs(g_seams * scale - g_ovito)))
        record["rdf_r"] = float(np.max(np.abs(r_seams - r_ovito)))
        expect(record["rdf_ovito"] <= ABS_TOL, f"{name} rdf {record['rdf_ovito']:.3e}")
        expect(record["rdf_r"] <= ABS_TOL, f"{name} rdf r {record['rdf_r']:.3e}")
        rdf = freud.density.RDF(bins=RDF_BINS, r_max=r_max)
        rdf.compute((freud_box(box), pts))
        g_freud = np.asarray(rdf.rdf, dtype=np.float64)
        record["rdf_freud"] = float(np.max(np.abs(g_ovito - g_freud)))
        width = r_max / RDF_BINS
        rho = n / float(volume.value)
        ideal = np.array(
            [rho * n * shell_volume(i, width) for i in range(RDF_BINS)], dtype=np.float64
        )
        ovito_counts = np.rint(g_ovito * ideal).astype(np.int64)
        # seams stores each unordered pair once. OVITO stores it twice.
        # freud's g(r) is that same N^2 histogram in float32, so a pair
        # whose float32 distance crosses a bin edge moves by one bin.
        directed = 2 * counts.astype(np.int64)
        record["rdf_pairs"] = int(np.sum(np.abs(ovito_counts - directed)))
        expect(record["rdf_pairs"] == 0, f"{name} rdf pair counts {record['rdf_pairs']}")
        freud_limit = 1e-5 if name in ("ice_ih", "ice_ic", "si") else 1e-3
        expect(
            record["rdf_freud"] <= freud_limit,
            f"{name} freud rdf {record['rdf_freud']:.3e}",
        )

    degree = float(np.mean([len(row) for row in cutoff_rows])) if cutoff_rows else 0.0
    record["degree"] = degree
    cn = ctypes.c_double()
    if status_ok(lib, lib.seams_coordination(handle, CUTOFF, 100, 0, 0, ctypes.byref(cn)), "cn"):
        scaled = float(cn.value) * (n - 1) / n
        record["cn"] = abs(scaled - degree)
        expect(record["cn"] <= ABS_TOL, f"{name} coordination {scaled:.6g} degree {degree:.6g}")

    lib.seams_system_destroy(handle)
    return record


def analytic_ql(lib: ctypes.CDLL, pts: np.ndarray, box: np.ndarray) -> None:
    handle = open_system(lib, pts, box)
    if handle is None:
        return
    rows = seams_rows(lib, handle, "cutoff")
    lib.seams_system_destroy(handle)
    populated = next(i for i, row in enumerate(rows) if row)
    vecs = mic_rows(pts, box, populated, list(rows[populated]))
    qlm = qlm_mean(vecs, 6)
    analytic = q_from_qlm(qlm, 6)
    qlm32 = qlm.astype(np.complex64).astype(np.complex128)
    rounded = q_from_qlm(qlm32, 6)
    n = len(pts)
    handle = open_system(lib, pts, box)
    q6 = (ctypes.c_double * n)()
    status_ok(lib, lib.seams_ql(handle, CUTOFF, 0, 6, 0, q6, n, None, 0), "analytic q6")
    lib.seams_system_destroy(handle)
    seams = float(np.ctypeslib.as_array(q6)[populated])
    ovito = float(steinhardt_column(pts, box, 6, False)[populated])
    op = freud.order.Steinhardt(l=6)
    op.compute((freud_box(box), pts), neighbors={"r_max": CUTOFF, "exclude_ii": True})
    fr = float(np.asarray(op.ql, dtype=np.float64)[populated])
    expect(abs(seams - analytic) <= ABS_TOL, f"analytic q6 seams {abs(seams - analytic):.3e}")
    expect(abs(ovito - analytic) <= ABS_TOL, f"analytic q6 ovito {abs(ovito - analytic):.3e}")
    expect(abs(fr - rounded) <= abs(fr - analytic) + 1e-8, "freud q6 is not the float32 sum")
    print(
        f"analytic q6 atom {populated} value {analytic:.8f} "
        f"seams {abs(seams - analytic):.3e} ovito {abs(ovito - analytic):.3e} "
        f"freud {abs(fr - analytic):.3e} float32 {abs(fr - rounded):.3e}"
    )


def analytic_average(lib: ctypes.CDLL) -> None:
    """Three oxygens. q6-bar is the norm of the mean q_6m, self included."""
    side = 30.0
    pts = np.array([[0.0, 0.0, 0.0], [2.0, 0.0, 0.0], [0.0, 2.0, 0.0]], dtype=np.float64)
    box = np.array([0.0, 0.0, 0.0, side, side, side, 0.0, 0.0, 0.0])
    cutoff = 2.5
    handle = open_system(lib, pts, box)
    if handle is None:
        return
    owned = ctypes.c_void_p()
    status_ok(lib, lib.seams_neighbours_build(handle, cutoff, 0, 0, ctypes.byref(owned)), "tri")
    off, cols = copy_csr(lib, owned)
    lib.seams_neighbours_destroy(owned)
    rows = rows_of(off, cols)
    q6 = (ctypes.c_double * 3)()
    bar = (ctypes.c_double * 3)()
    status_ok(lib, lib.seams_ql(handle, cutoff, 0, 6, 0, q6, 3, bar, 3), "tri q6")
    lib.seams_system_destroy(handle)
    local = []
    comps = []
    for index, row in enumerate(rows):
        comp = qlm_mean(mic_rows(pts, box, index, list(row)), 6)
        comps.append(comp)
        local.append(q_from_qlm(comp, 6))
    seams_local = np.ctypeslib.as_array(q6).copy()
    seams_bar = np.ctypeslib.as_array(bar).copy()
    expect(max_abs(seams_local, np.array(local)) <= ABS_TOL, "triangle local q6")
    # Lechner-Dellago: mean of q_lm over the atom and the atoms it lists.
    for index, row in enumerate(rows):
        group = [comps[index]] + [comps[j] for j in row]
        averaged = q_from_qlm(np.mean(np.stack(group), axis=0), 6)
        scalar = float(np.mean([local[index]] + [local[j] for j in row]))
        expect(abs(float(seams_bar[index]) - averaged) <= ABS_TOL, f"q6bar atom {index}")
        expect(abs(averaged - scalar) > 1e-4, "scalar mean of q6 collided with q_lm mean")
    one = np.array([[0.0, 0.0, 0.0], [2.0, 0.0, 0.0]], dtype=np.float64)
    handle = open_system(lib, one, box)
    q6 = (ctypes.c_double * 2)()
    status_ok(lib, lib.seams_ql(handle, cutoff, 0, 6, 0, q6, 2, None, 0), "pair q6")
    lib.seams_system_destroy(handle)
    expect(max_abs(np.ctypeslib.as_array(q6), np.array([1.0, 1.0])) <= ABS_TOL, "one-neighbour q6")
    print("analytic q6bar Lechner-Dellago one-neighbour q6=1")


def analytic_cutoff(lib: ctypes.CDLL) -> None:
    """OVITO keeps a pair at the cutoff. d-SEAMS and freud keep a pair strictly inside it."""
    side = 20.0
    box = np.array([0.0, 0.0, 0.0, side, side, side, 0.0, 0.0, 0.0])
    exact = np.array([[0.0, 0.0, 0.0], [CUTOFF, 0.0, 0.0]], dtype=np.float64)
    under = np.array([[0.0, 0.0, 0.0], [CUTOFF - 1e-4, 0.0, 0.0]], dtype=np.float64)

    def degree(pts: np.ndarray) -> tuple[int, int, int]:
        handle = open_system(lib, pts, box)
        rows = seams_rows(lib, handle, "cutoff") if handle else []
        if handle:
            lib.seams_system_destroy(handle)
        freud_rows = freud_directed(pts, box, {"mode": "ball", "r_max": CUTOFF, "exclude_ii": True})
        ovito_rows = ovito_bonds(pts, box, CUTOFF)
        return (
            sum(len(row) for row in rows),
            sum(len(row) for row in freud_rows),
            sum(len(row) for row in ovito_rows),
        )

    seams_n, freud_n, ovito_n = degree(exact)
    expect(seams_n == 0 and freud_n == 0 and ovito_n == 2, f"exact cutoff {seams_n} {freud_n} {ovito_n}")
    seams_n, freud_n, ovito_n = degree(under)
    expect(seams_n == 2 and freud_n == 2 and ovito_n == 2, f"inner cutoff {seams_n} {freud_n} {ovito_n}")
    print(f"analytic cutoff exact ovito-only under all-three")


def analytic_knn(lib: ctypes.CDLL) -> None:
    """k = 1 on a line. The C rows are the undirected closure of the directed list."""
    pts = np.array(
        [[0.0, 0.0, 0.0], [1.0, 0.0, 0.0], [3.0, 0.0, 0.0], [6.0, 0.0, 0.0]],
        dtype=np.float64,
    )
    box = np.array([0.0, 0.0, 0.0, 40.0, 40.0, 40.0, 0.0, 0.0, 0.0])
    handle = open_system(lib, pts, box)
    if handle is None:
        return
    one_way = seams_rows(lib, handle, "knn", k=1)
    mutual = seams_rows(lib, handle, "mutual", k=1)
    lib.seams_system_destroy(handle)
    directed = freud_directed(pts, box, {"mode": "nearest", "num_neighbors": 1, "exclude_ii": True})
    ovito = ovito_knn(pts, box, 1)
    expect(one_way == close_directed(directed, False), f"knn one-way {one_way}")
    expect(mutual == close_directed(directed, True), f"knn mutual {mutual}")
    expect(one_way == close_directed(ovito, False), "knn ovito one-way")
    expect(mutual == [(1,), (0,), (), ()], f"knn mutual rows {mutual}")
    expect(one_way[1] == (0, 2), f"knn asymmetric row {one_way[1]}")
    print(f"analytic knn one-way {one_way} mutual {mutual}")


def analytic_origin(lib: ctypes.CDLL) -> None:
    """A pair across the boundary of a cell whose origin is not zero."""
    origin = 2.0
    side = 10.0
    box = np.array([origin, origin, origin, side, side, side, 0.0, 0.0, 0.0])
    pts = np.array(
        [[origin + 0.1, 5.0, 5.0], [origin + side - 0.1, 5.0, 5.0]],
        dtype=np.float64,
    )
    handle = open_system(lib, pts, box)
    rows = seams_rows(lib, handle, "cutoff") if handle else []
    if handle:
        lib.seams_system_destroy(handle)
    freud_rows = [
        tuple(sorted(row))
        for row in freud_directed(pts, box, {"mode": "ball", "r_max": CUTOFF, "exclude_ii": True})
    ]
    ovito_rows = ovito_bonds(pts, box, CUTOFF)
    expect(rows == [(1,), (0,)] and rows == freud_rows and rows == ovito_rows, f"wrapped pair {rows}")
    print("analytic periodic image distance 0.2")


def analytic_rdf(lib: ctypes.CDLL, pts: np.ndarray, box: np.ndarray) -> None:
    """N(N-1) makes the shell integral N/(N-1) times the neighbour count."""
    n = len(pts)
    handle = open_system(lib, pts, box)
    if handle is None:
        return
    rows = seams_rows(lib, handle, "cutoff")
    degree = float(np.mean([len(row) for row in rows]))
    cn = ctypes.c_double()
    status_ok(lib, lib.seams_coordination(handle, CUTOFF, 100, 0, 0, ctypes.byref(cn)), "ice cn")
    lib.seams_system_destroy(handle)
    scaled = float(cn.value) * (n - 1) / n
    expect(abs(scaled - degree) <= ABS_TOL, f"scaled cn {scaled} degree {degree}")
    expect(abs(float(cn.value) - degree) > 1e-3, "rdf normalisation already matched N^2")
    small = np.array(
        [[0.0, 0.0, 0.0], [1.0, 0.0, 0.0], [0.0, 1.0, 0.0], [0.0, 0.0, 1.0]],
        dtype=np.float64,
    )
    small_box = np.array([0.0, 0.0, 0.0, 10.0, 10.0, 10.0, 0.0, 0.0, 0.0])
    handle = open_system(lib, small, small_box)
    radius = (ctypes.c_double * 10)()
    gram = (ctypes.c_double * 10)()
    hist = (ctypes.c_int * 10)()
    n_out = ctypes.c_size_t()
    volume = ctypes.c_double()
    code = lib.seams_rdf(
        handle, 6.0, 10, 0, 0, radius, gram, hist, 10, ctypes.byref(n_out), ctypes.byref(volume)
    )
    lib.seams_system_destroy(handle)
    expect(code == 0, "rdf past half an edge")
    raised = False
    try:
        freud.density.RDF(bins=10, r_max=6.0).compute((freud.box.Box.cube(10.0), small))
    except RuntimeError:
        raised = True
    expect(raised, "freud accepted r_max past half the shortest edge")
    print(
        f"analytic rdf cn {float(cn.value):.6f} scaled {scaled:.6f} degree {degree:.6f} "
        f"freud refuses r_max > L/2"
    )


def print_record(record: dict) -> None:
    keys = (
        "n",
        "cutoff",
        "cutoff_ovito",
        "knn",
        "knn_ovito",
        "mutual",
        "q4_ovito",
        "q6_ovito",
        "q6bar_ovito",
        "q4_freud",
        "q6_freud",
        "q6bar_freud",
        "chill",
        "ptm_split",
        "rdf_ovito",
        "rdf_freud",
        "cn",
        "degree",
    )
    parts = [record["frame"]]
    for key in keys:
        if key not in record:
            continue
        value = record[key]
        if isinstance(value, float):
            parts.append(f"{key} {value:.3e}")
        else:
            parts.append(f"{key} {value}")
    print(" ".join(parts))


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--root", required=True)
    args = parser.parse_args()
    root = Path(args.root)
    lib = load_library()
    loaded = frames_of(root)
    ic_pts, ic_box = loaded[1][1], loaded[1][2]
    analytic_cutoff(lib)
    analytic_knn(lib)
    analytic_origin(lib)
    analytic_average(lib)
    analytic_ql(lib, loaded[3][1], loaded[3][2])
    analytic_rdf(lib, ic_pts, ic_box)
    for name, pts, box in loaded:
        print_record(compare_frame(lib, name, pts, box))
    if FAILURES:
        print(f"{len(FAILURES)} crosscheck comparisons failed", file=sys.stderr)
        raise SystemExit(1)
    print("capi crosscheck ok")


if __name__ == "__main__":
    main()
