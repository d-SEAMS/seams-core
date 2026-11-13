#!/usr/bin/env python3
"""Tilted ice Ih and sH against freud, OVITO and mdapy.

Pinned ABAB, 20 samples. Cutoff pairs, k nearest neighbours, q6, both
CHILL+ modes, and primitive rings. A later scale run is 442368 molecules.
"""

from __future__ import annotations

import ctypes
import os
import sys
import time
from pathlib import Path

# Pin before libraries start worker threads. taskset on the process is the
# same mask; this keeps a later import from being the first word on it.
if len(sys.argv) >= 3 and sys.argv[1] in ("time", "scale"):
    try:
        os.sched_setaffinity(0, set(range(int(sys.argv[2]))))
    except OSError:
        pass

import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parent))

from h2h_cells import nearest_sH, tilted_ice
from h2h_freud_ovito import (
    CUTOFF,
    KNN,
    SEAMS_CHILL_OVITO,
    SEAMS_GRAPH_CUTOFF,
    copy_csr,
    core_count,
    library_path,
    load_seams,
    pairs_from_csr,
    require,
    same_edges,
    set_ovito_threads,
    undirected,
)

SEAMS_CHILL_PAPER = 0
N_ABAB = 20
N_SCALE = 8


def quiet():
    from PySide6.QtCore import QThreadPool

    QThreadPool.globalInstance().waitForDone(1000)


def time_once(fn) -> float:
    quiet()
    t0 = time.perf_counter_ns()
    fn()
    return (time.perf_counter_ns() - t0) / 1e6


def stats(samples):
    arr = np.asarray(samples, dtype=np.float64)
    q1, med, q3 = np.percentile(arr, [25, 50, 75])
    return float(med), float(q3 - q1), float(arr.min())


def open_frame(lib, pos, box9):
    held = np.ascontiguousarray(pos, dtype=np.float64)
    n = int(held.shape[0])
    box = (ctypes.c_double * 9)(*[float(v) for v in box9])
    types = (ctypes.c_int * n)(*([1] * n))
    ids = (ctypes.c_int * n)(*range(1, n + 1))
    ptr = held.ctypes.data_as(ctypes.POINTER(ctypes.c_double))
    sys = lib.seams_system_create()
    require(lib, lib.seams_system_set_frame(sys, ptr, n, box, types, ids), "set_frame")
    return sys, held


def cell_matrix(box9):
    ox, oy, oz, lx, ly, lz, xy, xz, yz = [float(v) for v in box9]
    # Columns are a, b, c.
    matrix = np.array([[lx, xy, xz], [0.0, ly, yz], [0.0, 0.0, lz]], dtype=np.float64)
    origin = np.array([ox, oy, oz], dtype=np.float64)
    return matrix, origin


def mdapy_vectors(box9):
    matrix, _origin = cell_matrix(box9)
    # mdapy rows are a, b, c.
    return np.ascontiguousarray(matrix.T)


def freud_box(box9):
    import freud

    matrix, _origin = cell_matrix(box9)
    return freud.box.Box.from_matrix(matrix)


def ovito_data(pos, box9):
    from ovito.data import DataCollection

    matrix, origin = cell_matrix(box9)
    cell = np.zeros((3, 4), dtype=np.float64)
    cell[:, :3] = matrix
    cell[:, 3] = origin
    data = DataCollection()
    data.create_cell(cell, pbc=(True, True, True))
    parts = data.create_particles(count=len(pos))
    parts.create_property("Position", data=np.ascontiguousarray(pos))
    return data


def dseams_pairs(lib, sys):
    handle = ctypes.c_void_p()
    require(
        lib,
        lib.seams_neighbours_build(sys, CUTOFF, 1, SEAMS_GRAPH_CUTOFF, ctypes.byref(handle)),
        "cutoff",
    )
    offsets, cols = copy_csr(lib, handle)
    return pairs_from_csr(offsets, cols)


def dseams_knn(lib, sys, n):
    handle = ctypes.c_void_p()
    require(lib, lib.seams_neighbours_knn(sys, KNN, CUTOFF, 1, 0, ctypes.byref(handle)), "knn")
    offsets, cols = copy_csr(lib, handle)
    rows = []
    for i in range(n):
        rows.append(tuple(int(v) for v in cols[offsets[i] : offsets[i + 1]]))
    return rows


def dseams_q6(lib, sys, n):
    q = (ctypes.c_double * n)()
    bar = (ctypes.c_double * n)()
    require(lib, lib.seams_steinhardt(sys, CUTOFF, 1, 6, q, n, bar, n), "q6")
    return np.ctypeslib.as_array(q).copy(), np.ctypeslib.as_array(bar).copy()


def dseams_chill(lib, sys, n, rule):
    labels = (ctypes.c_int * n)()
    counts = (ctypes.c_int * 9)()
    require(
        lib,
        lib.seams_chill_plus(sys, CUTOFF, 1, SEAMS_GRAPH_CUTOFF, rule, labels, counts),
        "chill+",
    )
    return np.ctypeslib.as_array(labels).copy()


def bind_rings(lib):
    lib.seams_rings_find.argtypes = [
        ctypes.c_void_p, ctypes.c_double, ctypes.c_int, ctypes.c_int, ctypes.c_int,
        ctypes.POINTER(ctypes.c_void_p),
    ]
    lib.seams_rings_find.restype = ctypes.c_int
    lib.seams_rings_destroy.argtypes = [ctypes.c_void_p]
    lib.seams_rings_count.argtypes = [ctypes.c_void_p]
    lib.seams_rings_count.restype = ctypes.c_int
    lib.seams_rings_vertices.argtypes = [
        ctypes.c_void_p, ctypes.c_int, ctypes.POINTER(ctypes.c_int),
        ctypes.c_size_t, ctypes.POINTER(ctypes.c_size_t),
    ]
    lib.seams_rings_vertices.restype = ctypes.c_int
    lib.seams_rings_histogram.argtypes = [
        ctypes.c_void_p, ctypes.POINTER(ctypes.c_int), ctypes.c_int,
    ]
    lib.seams_rings_histogram.restype = ctypes.c_int


def dseams_ring_sets(lib, sys):
    rings = ctypes.c_void_p()
    require(lib, lib.seams_rings_find(sys, CUTOFF, 1, 0, 6, ctypes.byref(rings)), "rings")
    n_rings = lib.seams_rings_count(rings)
    found = set()
    hist = [0] * 7
    buf = (ctypes.c_int * 12)()
    n_out = ctypes.c_size_t()
    for index in range(n_rings):
        require(lib, lib.seams_rings_vertices(rings, index, buf, 12, ctypes.byref(n_out)), "vertices")
        members = tuple(sorted(int(buf[k]) for k in range(n_out.value)))
        found.add(members)
        if len(members) < len(hist):
            hist[len(members)] += 1
    lib.seams_rings_destroy(rings)
    return found, hist


def freud_pairs(pos, box9):
    import freud

    box = freud_box(box9)
    found = freud.locality.AABBQuery(box, pos).query(
        pos, {"r_max": CUTOFF, "exclude_ii": True}
    ).toNeighborList()
    return undirected(np.asarray(found[:, 0]), np.asarray(found[:, 1]))


def ovito_pairs(data):
    from ovito.modifiers import CreateBondsModifier
    from ovito.pipeline import Pipeline, StaticSource

    pipe = Pipeline(source=StaticSource(data=data))
    pipe.modifiers.append(CreateBondsModifier(cutoff=CUTOFF))
    out = pipe.compute()
    bonds = out.particles.bonds
    topo = np.asarray(bonds.topology)
    return undirected(topo[:, 0], topo[:, 1])


def mdapy_pairs(system):
    import mdapy

    neigh = mdapy.NeighborList.build(system, cutoff=CUTOFF)
    return pairs_from_csr(np.asarray(neigh.offsets), np.asarray(neigh.j))


def mic_distance(pos, box9, i, j):
    matrix, origin = cell_matrix(box9)
    inv = np.linalg.inv(matrix)
    delta = pos[j] - pos[i]
    frac = delta @ inv.T
    frac -= np.round(frac)
    cart = frac @ matrix.T
    return float(np.linalg.norm(cart))


def explain_pairs(name, pos, box9, reference, other, other_name):
    if same_edges(reference, other):
        print(f"check {name} cutoff {other_name} pairs {reference.shape[0]} match exact", flush=True)
        return True
    ref = {tuple(row) for row in reference}
    got = {tuple(row) for row in other}
    only_ref = sorted(ref - got)
    only_got = sorted(got - ref)
    print(
        f"check {name} cutoff {other_name} dseams {len(ref)} other {len(got)} "
        f"only_dseams {len(only_ref)} only_other {len(only_got)}",
        flush=True,
    )
    shown = 0
    for i, j in (only_ref + only_got)[:6]:
        dist = mic_distance(pos, box9, int(i), int(j))
        side = "dseams" if (i, j) in ref and (i, j) not in got else other_name
        print(f"  pair {i} {j} mic_A {dist:.8f} only {side}", flush=True)
        shown += 1
    return False


def frames_for_check():
    out = []
    for nx, ny, nz, label in ((4, 4, 4, "ice512"), (8, 8, 8, "ice4096")):
        pos, box = tilted_ice(nx, ny, nz)
        out.append((label, pos, box))
    pos, box, rep = nearest_sH(4096)
    out.append((f"sH{rep[3]}", pos, box))
    return out


def check_frame(lib, name, pos, box9):
    import freud
    import mdapy
    import rings_statistics_matter
    from ase import Atoms
    from ovito.data import NearestNeighborFinder
    from ovito.modifiers import ChillPlusModifier, SteinhardtParametersModifier
    from ovito.pipeline import Pipeline, StaticSource

    n = len(pos)
    sys, _held = open_frame(lib, pos, box9)
    require(lib, lib.seams_set_num_threads(1), "threads")
    d_pairs = dseams_pairs(lib, sys)
    f_pairs = freud_pairs(pos, box9)
    data = ovito_data(pos, box9)
    o_pairs = ovito_pairs(data)
    system = mdapy.System(pos, box=mdapy_vectors(box9))
    m_pairs = mdapy_pairs(system)
    ok = True
    ok &= explain_pairs(name, pos, box9, d_pairs, f_pairs, "freud")
    ok &= explain_pairs(name, pos, box9, d_pairs, o_pairs, "ovito")
    ok &= explain_pairs(name, pos, box9, d_pairs, m_pairs, "mdapy")

    d_knn = dseams_knn(lib, sys, n)
    fbox = freud_box(box9)
    f_nn = freud.locality.AABBQuery(fbox, pos).query(
        pos, {"mode": "nearest", "num_neighbors": KNN, "r_max": CUTOFF + 1.5, "exclude_ii": True}
    ).toNeighborList()
    f_knn = [[] for _ in range(n)]
    for i, j in zip(f_nn[:, 0], f_nn[:, 1]):
        f_knn[int(i)].append(int(j))
    f_knn = [tuple(row) for row in f_knn]
    o_idx, _vec = NearestNeighborFinder(KNN, data).find_all()
    o_knn = [tuple(int(v) for v in row) for row in np.asarray(o_idx)]
    m_nn = mdapy.NeighborList.build(system, num_neighbors=KNN)
    m_off = np.asarray(m_nn.offsets)
    m_col = np.asarray(m_nn.j)
    m_knn = [tuple(int(v) for v in m_col[m_off[i] : m_off[i + 1]]) for i in range(n)]
    for other_name, other in (("freud", f_knn), ("ovito", o_knn), ("mdapy", m_knn)):
        miss = sum(tuple(sorted(a)) != tuple(sorted(b)) for a, b in zip(d_knn, other))
        print(f"check {name} knn {other_name} atom_mismatches {miss}", flush=True)
        if miss:
            ok = False
            for i, (a, b) in enumerate(zip(d_knn, other)):
                if tuple(sorted(a)) != tuple(sorted(b)):
                    da = [mic_distance(pos, box9, i, j) for j in a]
                    db = [mic_distance(pos, box9, i, j) for j in b]
                    print(f"  atom {i} dseams {a} dist {da} {other_name} {b} dist {db}", flush=True)
                    break

    d_q, d_bar = dseams_q6(lib, sys, n)
    local_q = freud.order.Steinhardt(l=6, average=False)
    local_q.compute((fbox, pos), neighbors={"r_max": CUTOFF, "exclude_ii": True})
    avg_q = freud.order.Steinhardt(l=6, average=True)
    avg_q.compute((fbox, pos), neighbors={"r_max": CUTOFF, "exclude_ii": True})
    print(
        f"check {name} q6 freud max_abs {np.max(np.abs(d_q - np.asarray(local_q.ql))):.3e} "
        f"qbar_max_abs {np.max(np.abs(d_bar - np.asarray(avg_q.ql))):.3e}",
        flush=True,
    )
    stein = SteinhardtParametersModifier(
        neighbor_mode=SteinhardtParametersModifier.NeighborMode.Cutoff,
        cutoff=CUTOFF,
        max_l=6,
        compute_averaged=True,
    )
    stein_pipe = Pipeline(source=StaticSource(data=data))
    stein_pipe.modifiers.append(stein)
    sout = stein_pipe.compute()
    st = np.asarray(sout.particles["Steinhardt"])
    sm = np.asarray(sout.particles["Smooth Steinhardt"])
    q6_ov = float(np.max(np.abs(d_q - st[:, 4])))
    b6_ov = float(np.max(np.abs(d_bar - sm[:, 4])))
    print(f"check {name} q6 ovito max_abs {q6_ov:.3e} qbar_max_abs {b6_ov:.3e}", flush=True)
    if q6_ov > 1e-8 or b6_ov > 1e-8:
        ok = False
    local = mdapy.analysis.SteinhardtParameters((6,), cutoff=CUTOFF, average=False).compute(system)
    averaged = mdapy.analysis.SteinhardtParameters((6,), cutoff=CUTOFF, average=True).compute(system)
    mq = np.asarray(local.q, dtype=np.float64).reshape(n, -1)[:, 0]
    mb = np.asarray(averaged.q, dtype=np.float64).reshape(n, -1)[:, 0]
    q6_m = float(np.max(np.abs(d_q - mq)))
    b6_m = float(np.max(np.abs(d_bar - mb)))
    print(f"check {name} q6 mdapy max_abs {q6_m:.3e} qbar_max_abs {b6_m:.3e}", flush=True)
    if q6_m > 1e-8 or b6_m > 1e-8:
        ok = False

    paper = dseams_chill(lib, sys, n, SEAMS_CHILL_PAPER)
    ovito_mode = dseams_chill(lib, sys, n, SEAMS_CHILL_OVITO)
    chill = ChillPlusModifier(cutoff=CUTOFF)
    chill_pipe = Pipeline(source=StaticSource(data=data))
    chill_pipe.modifiers.append(chill)
    cl = np.asarray(chill_pipe.compute().particles["Structure Type"], dtype=np.int32)
    ice = mdapy.identify.ChillPlus(cutoff=CUTOFF).compute(system)
    ml = np.asarray(ice.structure, dtype=np.int32)
    # mdapy and OVITO number the same labels. Compare ovito mode directly.
    ov_miss = int(np.count_nonzero(ovito_mode != cl))
    md_miss = int(np.count_nonzero(ovito_mode != ml))
    paper_miss = int(np.count_nonzero(paper != ovito_mode))
    print(
        f"check {name} chill ovito_mode_vs_ovito {ov_miss} "
        f"ovito_mode_vs_mdapy {md_miss} paper_vs_ovito_mode {paper_miss}",
        flush=True,
    )
    if ov_miss or md_miss:
        ok = False

    d_rings, d_hist = dseams_ring_sets(lib, sys)
    atoms = Atoms(symbols=["O"] * n, positions=pos, cell=mdapy_vectors(box9), pbc=True)
    _counts, blocks = rings_statistics_matter.ring_statistics(
        atoms, cutoff=CUTOFF, maxlvl=8, supercell=False, deduplicate=True,
        backend="rust", num_threads=1,
    )
    r_rings = set()
    for block in blocks:
        arr = np.asarray(block)
        if arr.size == 0 or arr.ndim != 2:
            continue
        if arr.shape[1] > 6:
            continue
        for row in arr:
            r_rings.add(tuple(sorted(int(v) for v in row)))
    only_d = len(d_rings - r_rings)
    only_r = len(r_rings - d_rings)
    print(
        f"check {name} rings dseams {len(d_rings)} ringsstat {len(r_rings)} "
        f"only_dseams {only_d} only_ringsstat {only_r} hist {d_hist}",
        flush=True,
    )
    if only_d or only_r:
        ok = False
    lib.seams_system_destroy(sys)
    return ok


def report_pair(op, threads, molecules, rival, a, b):
    am, ai, amin = stats(a)
    bm, bi, bmin = stats(b)
    win = am < bm and amin < bmin
    print(
        f"robust op {op} threads {threads} molecules {molecules} "
        f"dseams_median {am:.4f} dseams_iqr {ai:.4f} dseams_min {amin:.4f} "
        f"{rival}_median {bm:.4f} {rival}_iqr {bi:.4f} {rival}_min {bmin:.4f} "
        f"{'WIN' if win else 'LOSE'}",
        flush=True,
    )
    return win


def time_frame(lib, name, pos, box9, width, samples, only=None):
    import freud
    import mdapy
    from ovito.modifiers import (
        ChillPlusModifier,
        CreateBondsModifier,
        SteinhardtParametersModifier,
    )
    from ovito.pipeline import Pipeline, StaticSource
    from ovito.data import NearestNeighborFinder

    try:
        os.sched_setaffinity(0, set(range(width if width > 1 else 1)))
    except OSError:
        pass
    n = len(pos)
    require(lib, lib.seams_set_num_threads(width), "threads")
    mdapy.set_num_threads(width)
    freud.parallel.set_num_threads(width)
    set_ovito_threads(width)
    sys, _held = open_frame(lib, pos, box9)
    fbox = freud_box(box9)
    system = mdapy.System(pos, box=mdapy_vectors(box9))
    data = ovito_data(pos, box9)
    bonds = CreateBondsModifier(cutoff=CUTOFF)
    bond_pipe = Pipeline(source=StaticSource(data=data))
    bond_pipe.modifiers.append(bonds)
    stein = SteinhardtParametersModifier(
        neighbor_mode=SteinhardtParametersModifier.NeighborMode.Cutoff,
        cutoff=CUTOFF, max_l=6, compute_averaged=True,
    )
    stein_pipe = Pipeline(source=StaticSource(data=data))
    stein_pipe.modifiers.append(stein)
    chill = ChillPlusModifier(cutoff=CUTOFF)
    chill_pipe = Pipeline(source=StaticSource(data=data))
    chill_pipe.modifiers.append(chill)
    q = (ctypes.c_double * n)()
    bar = (ctypes.c_double * n)()
    labels = (ctypes.c_int * n)()
    counts = (ctypes.c_int * 9)()

    def d_cutoff():
        handle = ctypes.c_void_p()
        require(
            lib,
            lib.seams_neighbours_build(sys, CUTOFF, 1, SEAMS_GRAPH_CUTOFF, ctypes.byref(handle)),
            "cutoff",
        )
        lib.seams_neighbours_destroy(handle)

    def d_knn():
        handle = ctypes.c_void_p()
        require(lib, lib.seams_neighbours_knn(sys, KNN, CUTOFF, 1, 0, ctypes.byref(handle)), "knn")
        lib.seams_neighbours_destroy(handle)

    def d_q6():
        require(lib, lib.seams_steinhardt(sys, CUTOFF, 1, 6, q, n, bar, n), "q6")

    def d_paper():
        require(
            lib,
            lib.seams_chill_plus(sys, CUTOFF, 1, SEAMS_GRAPH_CUTOFF, SEAMS_CHILL_PAPER, labels, counts),
            "paper",
        )

    def d_chill():
        require(
            lib,
            lib.seams_chill_plus(sys, CUTOFF, 1, SEAMS_GRAPH_CUTOFF, SEAMS_CHILL_OVITO, labels, counts),
            "ovito-chill",
        )

    def d_rings():
        rings = ctypes.c_void_p()
        require(lib, lib.seams_rings_find(sys, CUTOFF, 1, 0, 6, ctypes.byref(rings)), "rings")
        lib.seams_rings_destroy(rings)

    def f_cutoff():
        freud.locality.AABBQuery(fbox, pos).query(
            pos, {"r_max": CUTOFF, "exclude_ii": True}
        ).toNeighborList()

    def f_knn():
        freud.locality.AABBQuery(fbox, pos).query(
            pos, {"mode": "nearest", "num_neighbors": KNN, "r_max": CUTOFF + 1.5, "exclude_ii": True}
        ).toNeighborList()

    def f_q6():
        op = freud.order.Steinhardt(l=6, average=True)
        op.compute((fbox, pos), neighbors={"r_max": CUTOFF, "exclude_ii": True})

    def o_cutoff():
        bonds.notify_object_changed()
        bond_pipe.compute()

    def o_knn():
        NearestNeighborFinder(KNN, data).find_all()

    def o_q6():
        stein.notify_object_changed()
        stein_pipe.compute()

    def o_chill():
        chill.notify_object_changed()
        chill_pipe.compute()

    def o_rings():
        # Built inside the call so the cache cannot skip the search.
        pipe = Pipeline(source=StaticSource(data=data))
        pipe.modifiers.append(CreateBondsModifier(cutoff=CUTOFF))
        from ovito.modifiers import FindRingsModifier
        pipe.modifiers.append(FindRingsModifier(maximum_ring_size=6))
        pipe.compute()

    def m_cutoff():
        mdapy.NeighborList.build(system, cutoff=CUTOFF)

    def m_knn():
        mdapy.NeighborList.build(system, num_neighbors=KNN)

    def m_q6():
        mdapy.analysis.SteinhardtParameters((6,), cutoff=CUTOFF, average=False).compute(system)
        mdapy.analysis.SteinhardtParameters((6,), cutoff=CUTOFF, average=True).compute(system)

    def m_chill():
        mdapy.identify.ChillPlus(cutoff=CUTOFF).compute(system)

    def m_rings():
        from ase import Atoms
        import rings_statistics_matter
        atoms = Atoms(symbols=["O"] * n, positions=pos, cell=mdapy_vectors(box9), pbc=True)
        rings_statistics_matter.ring_statistics(
            atoms, cutoff=CUTOFF, maxlvl=8, supercell=False, deduplicate=True,
            backend="rust", num_threads=width,
        )

    jobs = {
        "cutoff": (d_cutoff, {"freud": f_cutoff, "ovito": o_cutoff, "mdapy": m_cutoff}),
        "knn": (d_knn, {"freud": f_knn, "ovito": o_knn, "mdapy": m_knn}),
        "q6": (d_q6, {"freud": f_q6, "ovito": o_q6, "mdapy": m_q6}),
        "chill_paper": (d_paper, {}),
        "chill_plus": (d_chill, {"ovito": o_chill, "mdapy": m_chill}),
        "rings": (d_rings, {"ovito": o_rings, "ringsstat": m_rings}),
    }
    wins = True
    n_samples = samples
    for op, (ours, rivals) in jobs.items():
        if only is not None and op not in only:
            continue
        ours()
        for fn in rivals.values():
            fn()
    for op, (ours, rivals) in jobs.items():
        if only is not None and op not in only:
            continue
        if not rivals:
            got = [time_once(ours) for _ in range(n_samples)]
            med, iqr, amin = stats(got)
            print(
                f"robust op {op} threads {width} molecules {n} frame {name} "
                f"dseams_median {med:.4f} dseams_iqr {iqr:.4f} dseams_min {amin:.4f}",
                flush=True,
            )
            continue
        for rival, fn in rivals.items():
            a, b = [], []
            for _ in range(n_samples):
                a.append(time_once(ours))
                b.append(time_once(fn))
            if not report_pair(f"{name} {op}", width, n, rival, a, b):
                wins = False
    lib.seams_system_destroy(sys)
    return wins


def main():
    mode = sys.argv[1] if len(sys.argv) > 1 else "all"
    os.environ.setdefault("QT_QPA_PLATFORM", "offscreen")
    os.environ.setdefault("OMP_PROC_BIND", "close")
    os.environ.setdefault("OMP_PLACES", "cores")
    os.environ.setdefault("OMP_WAIT_POLICY", "PASSIVE")
    import freud
    import mdapy
    import ovito

    lib = load_seams(library_path())
    bind_rings(lib)
    if len(sys.argv) >= 3 and sys.argv[1] in ("time", "scale"):
        width = int(sys.argv[2])
        try:
            os.sched_setaffinity(0, set(range(width)))
        except OSError:
            pass
        require(lib, lib.seams_set_num_threads(width), "threads")
    print(f"cores {os.cpu_count()} affinity {core_count()}", flush=True)
    print(f"freud {freud.__version__} ovito {'.'.join(str(v) for v in ovito.version)} mdapy {mdapy.__version__}", flush=True)
    print("cell triclinic LAMMPS tilt xy xz yz", flush=True)

    if mode in ("check", "all"):
        ok = True
        # Small cells first, then the three requested sizes.
        small = []
        pos, box = tilted_ice(2, 2, 2)
        small.append((f"ice{len(pos)}", pos, box))
        pos, box = replicate_one()
        small.append((f"sH{len(pos)}", pos, box))
        for name, pos, box in small:
            print(f"frame {name} molecules {len(pos)}", flush=True)
            ok &= check_frame(lib, name, pos, box)
        if mode == "check":
            raise SystemExit(0 if ok else 1)
        if not ok:
            raise SystemExit("triclinic check failed")

    if mode in ("time", "all"):
        width = int(sys.argv[2]) if len(sys.argv) > 2 else 1
        catalog = []
        for ncell, label in ((8, "ice4096"), (16, "ice32768"), (24, "ice110592")):
            pos, box = tilted_ice(ncell, ncell, ncell)
            catalog.append((label, pos, box))
        for target in (4096, 32768, 110592):
            pos, box, rep = nearest_sH(target)
            catalog.append((f"sH{rep[3]}", pos, box))
            print(f"sH target {target} replication {rep[0]}x{rep[1]}x{rep[2]} molecules {rep[3]}", flush=True)
        print(f"pinned_threads {width}", flush=True)
        losses = []
        if width == 1:
            for name, pos, box in catalog:
                print(f"agree {name} molecules {len(pos)}", flush=True)
                if not check_frame(lib, name, pos, box):
                    raise SystemExit(f"disagreement on {name}")
        for name, pos, box in catalog:
            print(f"frame {name} molecules {len(pos)}", flush=True)
            if not time_frame(lib, name, pos, box, width, N_ABAB):
                losses.append(name)
        if losses:
            print(f"losses {losses}", flush=True)
            raise SystemExit(1)

    if mode in ("scale", "all"):
        width = int(sys.argv[2]) if len(sys.argv) > 2 else 1
        pos, box = tilted_ice(24, 24, 96)
        assert len(pos) == 442368
        print(f"scale molecules {len(pos)} pinned_threads {width}", flush=True)
        # Cutoff, both CHILL+ modes, and rings. Drop knn and q6.
        keep = {"cutoff", "chill_paper", "chill_plus", "rings"}
        if not time_frame(lib, "ice442368", pos, box, width, N_SCALE, only=keep):
            raise SystemExit(f"scale loss threads {width}")


def replicate_one():
    from h2h_cells import replicate_sH
    return replicate_sH(1, 1, 1)


if __name__ == "__main__":
    main()
