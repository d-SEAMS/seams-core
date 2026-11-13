#!/usr/bin/env python3
"""One process: d-SEAMS, freud, and OVITO on ice and the mixed frames.

Eight samples at one thread and at eight. The record is the core count.
The host name is not written. CHILL+ is the ovito definition. Freud has
no CHILL+ and no diamond identifier. OVITO's polyhedral template
matching is run with the cubic and hexagonal diamond templates enabled,
which the modifier leaves off by default.
"""

from __future__ import annotations

import ctypes
import importlib.util
import math
import os
import statistics
import time
from pathlib import Path

os.environ.setdefault("OPENBLAS_NUM_THREADS", "1")
os.environ.setdefault("MKL_NUM_THREADS", "1")
os.environ.setdefault("NUMEXPR_NUM_THREADS", "1")

import numpy as np

ROOT = Path(__file__).resolve().parents[1]
SAMPLES = 8
SIZES = (8, 16, 24)
CUTOFF = 3.5
KNN = 4
SEAMS_GRAPH_CUTOFF = 0
SEAMS_CHILL_OVITO = 1

SYM = {0: "hex", 1: "cubic", 2: "stack", -1: "other"}
PTM = {0: "other", 6: "cubic", 7: "hex"}
DIA = {
    0: "other",
    1: "cubic",
    2: "cubic_1",
    3: "cubic_2",
    4: "hex",
    5: "hex_1",
    6: "hex_2",
}


def core_count() -> int:
    try:
        return len(os.sched_getaffinity(0))
    except (AttributeError, OSError):
        return os.cpu_count() or 1


def library_path() -> Path:
    env = os.environ.get("SEAMS_LIB", "")
    if env:
        return Path(env)
    for candidate in (
        ROOT / "bbdir-release" / "src" / "libseams.so",
        ROOT / "build" / "src" / "libseams.so",
    ):
        if candidate.is_file():
            return candidate
    raise SystemExit("libseams.so was not found; set SEAMS_LIB")


def load_seams(path: Path):
    yoda = path.parent / "libyodaLib.so"
    if yoda.is_file():
        ctypes.CDLL(str(yoda), mode=ctypes.RTLD_GLOBAL)
    lib = ctypes.CDLL(str(path))
    lib.seams_set_num_threads.argtypes = [ctypes.c_int]
    lib.seams_set_num_threads.restype = ctypes.c_int
    lib.seams_system_create.restype = ctypes.c_void_p
    lib.seams_system_destroy.argtypes = [ctypes.c_void_p]
    lib.seams_system_set_frame.argtypes = [
        ctypes.c_void_p,
        ctypes.POINTER(ctypes.c_double),
        ctypes.c_int,
        ctypes.POINTER(ctypes.c_double),
        ctypes.POINTER(ctypes.c_int),
        ctypes.POINTER(ctypes.c_int),
    ]
    lib.seams_system_set_frame.restype = ctypes.c_int
    lib.seams_neighbours_build.argtypes = [
        ctypes.c_void_p, ctypes.c_double, ctypes.c_int, ctypes.c_int,
        ctypes.POINTER(ctypes.c_void_p),
    ]
    lib.seams_neighbours_build.restype = ctypes.c_int
    lib.seams_neighbours_knn.argtypes = [
        ctypes.c_void_p, ctypes.c_int, ctypes.c_double, ctypes.c_int, ctypes.c_int,
        ctypes.POINTER(ctypes.c_void_p),
    ]
    lib.seams_neighbours_knn.restype = ctypes.c_int
    lib.seams_neighbours_destroy.argtypes = [ctypes.c_void_p]
    lib.seams_neighbours_copy.argtypes = [
        ctypes.c_void_p,
        ctypes.POINTER(ctypes.c_int), ctypes.c_size_t, ctypes.POINTER(ctypes.c_size_t),
        ctypes.POINTER(ctypes.c_int), ctypes.c_size_t, ctypes.POINTER(ctypes.c_size_t),
    ]
    lib.seams_neighbours_copy.restype = ctypes.c_int
    lib.seams_steinhardt.argtypes = [
        ctypes.c_void_p, ctypes.c_double, ctypes.c_int, ctypes.c_int,
        ctypes.POINTER(ctypes.c_double), ctypes.c_size_t,
        ctypes.POINTER(ctypes.c_double), ctypes.c_size_t,
    ]
    lib.seams_steinhardt.restype = ctypes.c_int
    lib.seams_chill_plus.argtypes = [
        ctypes.c_void_p, ctypes.c_double, ctypes.c_int, ctypes.c_int, ctypes.c_int,
        ctypes.POINTER(ctypes.c_int), ctypes.POINTER(ctypes.c_int),
    ]
    lib.seams_chill_plus.restype = ctypes.c_int
    lib.seams_system_symmetry.argtypes = [
        ctypes.c_void_p, ctypes.c_double, ctypes.POINTER(ctypes.c_int), ctypes.c_size_t,
    ]
    lib.seams_system_symmetry.restype = ctypes.c_int
    lib.seams_last_error.restype = ctypes.c_char_p
    return lib


def ice_positions(ncell: int):
    bond = 2.76
    a = bond * math.sqrt(8.0 / 3.0)
    c_len = a * math.sqrt(8.0 / 3.0)
    b = a * math.sqrt(3.0)
    basis = (
        (1.0 / 3.0, 2.0 / 3.0, 0.0),
        (2.0 / 3.0, 1.0 / 3.0, 0.5),
        (1.0 / 3.0, 2.0 / 3.0, 3.0 / 8.0),
        (2.0 / 3.0, 1.0 / 3.0, 7.0 / 8.0),
    )
    lx, ly, lz = ncell * a, ncell * b, ncell * c_len
    coords = []
    for ix in range(ncell):
        for iy in range(ncell):
            for iz in range(ncell):
                for shift in (0.0, 0.5):
                    for site in basis:
                        coords.append((
                            (shift * a + site[0] * a - site[1] * a / 2.0 + ix * a) % lx,
                            (shift * b + site[1] * b / 2.0 + iy * b) % ly,
                            (site[2] + iz) * c_len,
                        ))
    return np.asarray(coords, dtype=np.float64), (lx, ly, lz)


def mixed_frames():
    spec = importlib.util.spec_from_file_location(
        "mixed_phase_frames", ROOT / "scripts" / "mixed_phase_frames.py"
    )
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return [(name, builder()) for name, builder in module.FRAMES.items()]


def require(lib, status: int, what: str) -> None:
    if status == 0:
        return
    err = lib.seams_last_error()
    detail = err.decode() if err else "unknown"
    raise SystemExit(f"d-SEAMS {what} failed: {detail}")


def emit(lines, echo, text: str) -> None:
    lines.append(text)
    if echo:
        print(text, flush=True)


def samples_of(fn, n=SAMPLES) -> list[float]:
    fn()
    out = []
    for _ in range(n):
        t0 = time.perf_counter_ns()
        fn()
        out.append((time.perf_counter_ns() - t0) / 1e6)
    return out


def emit_time(lines, echo, molecules, code, op, threads, got) -> None:
    raw = " ".join(f"{v:.4f}" for v in got)
    emit(
        lines, echo,
        f"arm {code} op {op} threads {threads} molecules {molecules} samples_ms {raw}",
    )
    emit(
        lines, echo,
        f"arm {code} op {op} threads {threads} molecules {molecules} "
        f"median_ms {statistics.median(got):.4f}",
    )


def undirected(left, right) -> np.ndarray:
    a = np.minimum(left, right).astype(np.int64, copy=False)
    b = np.maximum(left, right).astype(np.int64, copy=False)
    pairs = np.unique(np.stack((a, b), axis=1), axis=0)
    return pairs


def pairs_from_csr(offsets, cols) -> np.ndarray:
    counts = np.diff(offsets)
    sources = np.repeat(np.arange(counts.shape[0], dtype=np.int64), counts)
    return undirected(sources, np.asarray(cols, dtype=np.int64))


def same_edges(left, right) -> bool:
    if left.shape != right.shape:
        return False
    order_l = np.lexsort((left[:, 1], left[:, 0]))
    order_r = np.lexsort((right[:, 1], right[:, 0]))
    return bool(np.array_equal(left[order_l], right[order_r]))


def copy_csr(lib, nb):
    n_offsets = ctypes.c_size_t()
    nnz = ctypes.c_size_t()
    require(
        lib,
        lib.seams_neighbours_copy(nb, None, 0, ctypes.byref(n_offsets), None, 0, ctypes.byref(nnz)),
        "neighbour length",
    )
    offsets = (ctypes.c_int * n_offsets.value)()
    cols = (ctypes.c_int * max(nnz.value, 1))()
    require(
        lib,
        lib.seams_neighbours_copy(
            nb, offsets, n_offsets.value, ctypes.byref(n_offsets),
            cols, max(nnz.value, 1), ctypes.byref(nnz),
        ),
        "neighbour copy",
    )
    off = np.ctypeslib.as_array(offsets).copy()
    col = np.ctypeslib.as_array(cols)[: nnz.value].copy()
    lib.seams_neighbours_destroy(nb)
    return off, col


def frame_system(lib, pos, lengths):
    held = np.ascontiguousarray(pos, dtype=np.float64)
    n = int(held.shape[0])
    lx, ly, lz = lengths
    box9 = (ctypes.c_double * 9)(0, 0, 0, lx, ly, lz, 0, 0, 0)
    types = (ctypes.c_int * n)(*([1] * n))
    ids = (ctypes.c_int * n)(*range(1, n + 1))
    ptr = held.ctypes.data_as(ctypes.POINTER(ctypes.c_double))
    sys = lib.seams_system_create()
    require(lib, lib.seams_system_set_frame(sys, ptr, n, box9, types, ids), "set_frame")
    return sys


def confusion(left, right) -> str:
    counts = {}
    for a, b in zip(left, right):
        key = (a, b)
        counts[key] = counts.get(key, 0) + 1
    return " ".join(f"{a}/{b}={counts[a, b]}" for a, b in sorted(counts))


def set_ovito_threads(width: int) -> None:
    from PySide6.QtCore import QThreadPool

    QThreadPool.globalInstance().setMaxThreadCount(width)


def run_session(echo: bool = False) -> list[str]:
    import freud
    from ovito.data import DataCollection, NearestNeighborFinder, Particles, SimulationCell
    from ovito.modifiers import (
        ChillPlusModifier,
        CreateBondsModifier,
        IdentifyDiamondModifier,
        PolyhedralTemplateMatchingModifier,
        SteinhardtParametersModifier,
    )
    from ovito.pipeline import Pipeline, StaticSource

    lib = load_seams(library_path())
    lines: list[str] = []
    emit(lines, echo, f"cores {core_count()}")
    emit(lines, echo, f"samples {SAMPLES}")
    emit(lines, echo, f"freud {freud.__version__}")
    import ovito
    emit(lines, echo, f"ovito {'.'.join(str(v) for v in ovito.version)}")
    emit(lines, echo, f"cutoff_A {CUTOFF}")
    emit(lines, echo, f"knn {KNN}")
    emit(lines, echo, "note freud ql is float32 and freud has no CHILL+ or diamond identifier")
    emit(lines, echo, "note OVITO diamond templates are enabled for polyhedral template matching")

    frames = []
    for ncell in SIZES:
        pos, lengths = ice_positions(ncell)
        frames.append((f"ice{pos.shape[0]}", pos, lengths))
    frames.extend((name, pos, lengths) for name, (pos, lengths) in mixed_frames())

    for name, pos, lengths in frames:
        pos = np.ascontiguousarray(pos, dtype=np.float64)
        n = int(pos.shape[0])
        lx, ly, lz = lengths
        emit(lines, echo, f"frame {name} molecules {n}")
        sys = frame_system(lib, pos, lengths)
        require(lib, lib.seams_set_num_threads(1), "threads")
        nb = ctypes.c_void_p()
        require(lib, lib.seams_neighbours_build(sys, CUTOFF, 1, SEAMS_GRAPH_CUTOFF, ctypes.byref(nb)), "cutoff")
        off, col = copy_csr(lib, nb)
        d_pairs = pairs_from_csr(off, col)
        require(lib, lib.seams_neighbours_knn(sys, KNN, CUTOFF, 1, 0, ctypes.byref(nb)), "knn")
        koff, kcol = copy_csr(lib, nb)
        d_knn = [tuple(sorted(int(v) for v in kcol[koff[i]:koff[i + 1]])) for i in range(n)]

        q6 = (ctypes.c_double * n)()
        b6 = (ctypes.c_double * n)()
        q3 = (ctypes.c_double * n)()
        b3 = (ctypes.c_double * n)()
        require(lib, lib.seams_steinhardt(sys, CUTOFF, 1, 6, q6, n, b6, n), "q6")
        require(lib, lib.seams_steinhardt(sys, CUTOFF, 1, 3, q3, n, b3, n), "q3")
        d_q6 = np.ctypeslib.as_array(q6).copy()
        d_b6 = np.ctypeslib.as_array(b6).copy()
        d_q3 = np.ctypeslib.as_array(q3).copy()
        d_b3 = np.ctypeslib.as_array(b3).copy()
        labels = (ctypes.c_int * n)()
        counts = (ctypes.c_int * 9)()
        require(
            lib,
            lib.seams_chill_plus(sys, CUTOFF, 1, SEAMS_GRAPH_CUTOFF, SEAMS_CHILL_OVITO, labels, counts),
            "chill+",
        )
        d_chill = np.ctypeslib.as_array(labels).copy()
        sym = (ctypes.c_int * n)()
        require(lib, lib.seams_system_symmetry(sys, 0.0, sym, n), "symmetry")
        d_sym = np.ctypeslib.as_array(sym).copy()

        fbox = freud.box.Box(lx, ly, lz)
        freud.parallel.set_num_threads(1)
        query = freud.locality.AABBQuery(fbox, pos)
        flist = query.query(pos, {"r_max": CUTOFF, "exclude_ii": True}).toNeighborList()
        f_pairs = undirected(
            np.asarray(flist.query_point_indices), np.asarray(flist.point_indices)
        )
        if same_edges(d_pairs, f_pairs):
            emit(lines, echo, f"check cutoff freud {name} pairs {d_pairs.shape[0]} match exact")
        else:
            emit(
                lines, echo,
                f"check cutoff freud {name} dseams {d_pairs.shape[0]} freud {f_pairs.shape[0]} differ",
            )
        nlist = query.query(
            pos, {"mode": "nearest", "num_neighbors": KNN, "r_max": CUTOFF + 1.5, "exclude_ii": True}
        ).toNeighborList()
        f_knn = [[] for _ in range(n)]
        for i, j in zip(nlist.query_point_indices, nlist.point_indices):
            f_knn[int(i)].append(int(j))
        f_knn = [tuple(sorted(row)) for row in f_knn]
        knn_miss = sum(a != b for a, b in zip(d_knn, f_knn))
        emit(lines, echo, f"check knn freud {name} atom_mismatches {knn_miss}")

        def freud_q(order):
            local = freud.order.Steinhardt(l=order, average=False)
            local.compute((fbox, pos), neighbors={"r_max": CUTOFF, "exclude_ii": True})
            averaged = freud.order.Steinhardt(l=order, average=True)
            averaged.compute((fbox, pos), neighbors={"r_max": CUTOFF, "exclude_ii": True})
            return np.asarray(local.ql, dtype=np.float64), np.asarray(averaged.ql, dtype=np.float64)

        for order, local, bar in ((6, d_q6, d_b6), (3, d_q3, d_b3)):
            fq, fb = freud_q(order)
            emit(
                lines, echo,
                f"check q{order} freud {name} max_abs {np.max(np.abs(local - fq)):.3e} "
                f"qbar_max_abs {np.max(np.abs(bar - fb)):.3e} freud_dtype float32",
            )

        data = DataCollection()
        cell = SimulationCell(pbc=(True, True, True))
        matrix = np.zeros((3, 4))
        matrix[0, 0], matrix[1, 1], matrix[2, 2] = lx, ly, lz
        cell[...] = matrix
        data.objects.append(cell)
        particles = Particles()
        particles.create_property("Position", data=pos)
        data.objects.append(particles)
        bonds_mod = CreateBondsModifier(cutoff=CUTOFF)
        bond_pipe = Pipeline(source=StaticSource(data=data))
        bond_pipe.modifiers.append(bonds_mod)
        bonded = bond_pipe.compute()
        topo = np.asarray(bonded.particles.bonds.topology)
        o_pairs = undirected(topo[:, 0], topo[:, 1])
        if same_edges(d_pairs, o_pairs):
            emit(lines, echo, f"check cutoff ovito {name} pairs {d_pairs.shape[0]} match exact")
        else:
            emit(
                lines, echo,
                f"check cutoff ovito {name} dseams {d_pairs.shape[0]} ovito {o_pairs.shape[0]} differ",
            )
        idx, _vecs = NearestNeighborFinder(KNN, data).find_all()
        o_knn = [tuple(sorted(int(v) for v in row)) for row in np.asarray(idx)]
        ov_knn_miss = sum(a != b for a, b in zip(d_knn, o_knn))
        emit(lines, echo, f"check knn ovito {name} atom_mismatches {ov_knn_miss}")

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
        # Columns are q2, q3, q4, q5, q6 when max_l is 6.
        emit(
            lines, echo,
            f"check q3 ovito {name} max_abs {np.max(np.abs(d_q3 - st[:, 1])):.3e} "
            f"qbar_max_abs {np.max(np.abs(d_b3 - sm[:, 1])):.3e}",
        )
        emit(
            lines, echo,
            f"check q6 ovito {name} max_abs {np.max(np.abs(d_q6 - st[:, 4])):.3e} "
            f"qbar_max_abs {np.max(np.abs(d_b6 - sm[:, 4])):.3e}",
        )

        chill = ChillPlusModifier(cutoff=CUTOFF)
        chill_pipe = Pipeline(source=StaticSource(data=data))
        chill_pipe.modifiers.append(chill)
        ov_chill = np.asarray(chill_pipe.compute().particles["Structure Type"])
        # d-SEAMS: 0 cubic, 1 hexagonal, 2 water, 3 interfacial, 4 clathrate, 5 interclathrate
        # OVITO: 0 other, 1 hexagonal, 2 cubic, 3 interfacial, 4 hydrate, 5 interfacial hydrate
        mapped = np.array([2, 1, 0, 3, 4, 5, 0, 0, 0], dtype=np.int64)[d_chill]
        chill_miss = int(np.count_nonzero(mapped != ov_chill))
        emit(lines, echo, f"check chill_plus ovito {name} atom_mismatches {chill_miss}")
        if chill_miss:
            raise SystemExit(f"{name}: CHILL+ ovito mode disagrees with the modifier")

        ptm = PolyhedralTemplateMatchingModifier(rmsd_cutoff=0.1)
        ptm.structures[PolyhedralTemplateMatchingModifier.Type.CUBIC_DIAMOND].enabled = True
        ptm.structures[PolyhedralTemplateMatchingModifier.Type.HEX_DIAMOND].enabled = True
        ptm_pipe = Pipeline(source=StaticSource(data=data))
        ptm_pipe.modifiers.append(ptm)
        ov_ptm = np.asarray(ptm_pipe.compute().particles["Structure Type"])
        dia = IdentifyDiamondModifier()
        dia_pipe = Pipeline(source=StaticSource(data=data))
        dia_pipe.modifiers.append(dia)
        ov_dia = np.asarray(dia_pipe.compute().particles["Structure Type"])
        sym_name = [SYM.get(int(v), "other") for v in d_sym]
        ptm_name = [PTM.get(int(v), "other") for v in ov_ptm]
        dia_name = [DIA.get(int(v), "other") for v in ov_dia]
        emit(lines, echo, f"check ptm {name} " + confusion(sym_name, ptm_name))
        emit(lines, echo, f"check diamond {name} " + confusion(sym_name, dia_name))

        q6b = (ctypes.c_double * n)()
        b6b = (ctypes.c_double * n)()
        q3b = (ctypes.c_double * n)()
        b3b = (ctypes.c_double * n)()
        lab_b = (ctypes.c_int * n)()
        cnt_b = (ctypes.c_int * 9)()
        sym_b = (ctypes.c_int * n)()

        def time_dseams(width: int) -> None:
            require(lib, lib.seams_set_num_threads(width), "threads")

            def cutoff():
                handle = ctypes.c_void_p()
                require(
                    lib,
                    lib.seams_neighbours_build(sys, CUTOFF, 1, SEAMS_GRAPH_CUTOFF, ctypes.byref(handle)),
                    "cutoff",
                )
                lib.seams_neighbours_destroy(handle)

            def knn():
                handle = ctypes.c_void_p()
                require(lib, lib.seams_neighbours_knn(sys, KNN, CUTOFF, 1, 0, ctypes.byref(handle)), "knn")
                lib.seams_neighbours_destroy(handle)

            def q6_call():
                require(lib, lib.seams_steinhardt(sys, CUTOFF, 1, 6, q6b, n, b6b, n), "q6")

            def q3_call():
                require(lib, lib.seams_steinhardt(sys, CUTOFF, 1, 3, q3b, n, b3b, n), "q3")

            def chill_call():
                require(
                    lib,
                    lib.seams_chill_plus(
                        sys, CUTOFF, 1, SEAMS_GRAPH_CUTOFF, SEAMS_CHILL_OVITO, lab_b, cnt_b
                    ),
                    "chill+",
                )

            def sym_call():
                require(lib, lib.seams_system_symmetry(sys, 0.0, sym_b, n), "symmetry")

            def q3q6_call():
                require(lib, lib.seams_steinhardt(sys, CUTOFF, 1, 3, q3b, n, b3b, n), "q3")
                require(lib, lib.seams_steinhardt(sys, CUTOFF, 1, 6, q6b, n, b6b, n), "q6")

            emit_time(lines, echo, n, "dseams", "cutoff", width, samples_of(cutoff))
            emit_time(lines, echo, n, "dseams", "knn", width, samples_of(knn))
            emit_time(lines, echo, n, "dseams", "q6", width, samples_of(q6_call))
            emit_time(lines, echo, n, "dseams", "q3", width, samples_of(q3_call))
            emit_time(lines, echo, n, "dseams", "q3q6", width, samples_of(q3q6_call))
            emit_time(lines, echo, n, "dseams", "chill_plus_ovito", width, samples_of(chill_call))
            emit_time(lines, echo, n, "dseams", "symmetry", width, samples_of(sym_call))

        def time_freud(width: int) -> None:
            freud.parallel.set_num_threads(width)

            def cutoff():
                found = freud.locality.AABBQuery(fbox, pos).query(
                    pos, {"r_max": CUTOFF, "exclude_ii": True}
                ).toNeighborList()
                if int(found.neighbor_counts.shape[0]) != n:
                    raise SystemExit("freud cutoff shape")

            def knn():
                found = freud.locality.AABBQuery(fbox, pos).query(
                    pos,
                    {"mode": "nearest", "num_neighbors": KNN, "r_max": CUTOFF + 1.5, "exclude_ii": True},
                ).toNeighborList()
                if int(found.neighbor_counts.shape[0]) != n:
                    raise SystemExit("freud knn shape")

            def q6_call():
                op = freud.order.Steinhardt(l=6, average=True)
                op.compute((fbox, pos), neighbors={"r_max": CUTOFF, "exclude_ii": True})
                if op.ql.shape[0] != n:
                    raise SystemExit("freud q6 shape")

            def q3_call():
                op = freud.order.Steinhardt(l=3, average=True)
                op.compute((fbox, pos), neighbors={"r_max": CUTOFF, "exclude_ii": True})
                if op.ql.shape[0] != n:
                    raise SystemExit("freud q3 shape")

            emit_time(lines, echo, n, "freud", "cutoff", width, samples_of(cutoff))
            emit_time(lines, echo, n, "freud", "knn", width, samples_of(knn))
            emit_time(lines, echo, n, "freud", "q6", width, samples_of(q6_call))
            emit_time(lines, echo, n, "freud", "q3", width, samples_of(q3_call))

        def time_ovito(width: int) -> None:
            set_ovito_threads(width)

            def cutoff():
                bonds_mod.notify_object_changed()
                out = bond_pipe.compute()
                if out.particles.bonds.count < 1 and n > 1:
                    raise SystemExit("ovito cutoff empty")

            def knn():
                found, _vecs = NearestNeighborFinder(KNN, data).find_all()
                if found.shape != (n, KNN):
                    raise SystemExit("ovito knn shape")

            def q_call():
                stein.notify_object_changed()
                out = stein_pipe.compute()
                if out.particles["Steinhardt"].shape[0] != n:
                    raise SystemExit("ovito q shape")

            def chill_call():
                chill.notify_object_changed()
                out = chill_pipe.compute()
                if out.particles["Structure Type"].shape[0] != n:
                    raise SystemExit("ovito chill shape")

            def ptm_call():
                ptm.notify_object_changed()
                out = ptm_pipe.compute()
                if out.particles["Structure Type"].shape[0] != n:
                    raise SystemExit("ovito ptm shape")

            def dia_call():
                dia.notify_object_changed()
                out = dia_pipe.compute()
                if out.particles["Structure Type"].shape[0] != n:
                    raise SystemExit("ovito diamond shape")

            emit_time(lines, echo, n, "ovito", "cutoff", width, samples_of(cutoff))
            emit_time(lines, echo, n, "ovito", "knn", width, samples_of(knn))
            emit_time(lines, echo, n, "ovito", "steinhardt_l6", width, samples_of(q_call))
            emit_time(lines, echo, n, "ovito", "chill_plus", width, samples_of(chill_call))
            emit_time(lines, echo, n, "ovito", "ptm", width, samples_of(ptm_call))
            emit_time(lines, echo, n, "ovito", "diamond", width, samples_of(dia_call))

        for width in (1, 8):
            time_dseams(width)
            time_freud(width)
            time_ovito(width)
        lib.seams_system_destroy(sys)
    return lines


def main() -> None:
    for line in run_session(echo=True):
        pass


if __name__ == "__main__":
    main()
