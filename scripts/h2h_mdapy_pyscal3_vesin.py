#!/usr/bin/env python3
"""One process: d-SEAMS, mdapy, pyscal3, and vesin on the same ice Ih cells.

Eight samples at 4096, 32768, and 110592 oxygen molecules. d-SEAMS, mdapy,
pyscal3, and vesin run at one thread and at eight where the tool has a
thread count. The record is the core count. The host name is not written.

Neighbour pairs are unordered and exclude the atom itself. q6 is the
cutoff Steinhardt parameter and q6bar is the Lechner-Dellago average of
q_lm, including the atom itself. CHILL+ has two definitions. paper, the
default, averages q_lm over the four closest neighbours. ovito averages
q_lm over every neighbour inside 3.5 A. mdapy and the OVITO modifier are
compared in ovito mode, atom for atom. On this ice the cutoff shell is
already those four neighbours.
"""

from __future__ import annotations

import ctypes
import math
import os
import statistics
import time
from pathlib import Path

# Keep BLAS off the cores the neighbour kernels are measured on.
os.environ.setdefault("OPENBLAS_NUM_THREADS", "1")
os.environ.setdefault("MKL_NUM_THREADS", "1")
os.environ.setdefault("NUMEXPR_NUM_THREADS", "1")

import numpy as np
from ase import Atoms

import mdapy
import pyscal3
from pyscal3._bridge import neighbor_arrays
import vesin

ROOT = Path(__file__).resolve().parents[1]
SAMPLES = 8
SIZES = (8, 16, 24)  # orthorhombic cells: 4096, 32768, 110592 molecules
CUTOFF = 3.5
KNN = 4
Q_ATOL = 1e-8

SEAMS_GRAPH_CUTOFF = 0
SEAMS_GRAPH_MUTUAL = 1
SEAMS_CHILL_PUBLISHED = 0
SEAMS_CHILL_OVITO = 1

# d-SEAMS label index -> mdapy / OVITO class name.
DSEAMS_AS_OVITO = {
    0: "cubic_ice",
    1: "hexagonal_ice",
    2: "other",
    3: "interfacial_ice",
    4: "hydrate",
    5: "interfacial_hydrate",
    6: "other",
    7: "other",
    8: "other",
}
OVITO_CODE = {
    0: "other",
    1: "hexagonal_ice",
    2: "cubic_ice",
    3: "interfacial_ice",
    4: "hydrate",
    5: "interfacial_hydrate",
}

DSEAMS_ICE = (
    "cubic",
    "hexagonal",
    "water",
    "interfacial",
    "clathrate",
    "interclathrate",
    "unclassified",
    "recubic",
    "rehex",
)

_PAUSE = None


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
        ctypes.c_void_p,
        ctypes.c_double,
        ctypes.c_int,
        ctypes.c_int,
        ctypes.POINTER(ctypes.c_void_p),
    ]
    lib.seams_neighbours_build.restype = ctypes.c_int
    lib.seams_neighbours_knn.argtypes = [
        ctypes.c_void_p,
        ctypes.c_int,
        ctypes.c_double,
        ctypes.c_int,
        ctypes.c_int,
        ctypes.POINTER(ctypes.c_void_p),
    ]
    lib.seams_neighbours_knn.restype = ctypes.c_int
    lib.seams_neighbours_destroy.argtypes = [ctypes.c_void_p]
    lib.seams_neighbours_copy.argtypes = [
        ctypes.c_void_p,
        ctypes.POINTER(ctypes.c_int),
        ctypes.c_size_t,
        ctypes.POINTER(ctypes.c_size_t),
        ctypes.POINTER(ctypes.c_int),
        ctypes.c_size_t,
        ctypes.POINTER(ctypes.c_size_t),
    ]
    lib.seams_neighbours_copy.restype = ctypes.c_int
    lib.seams_steinhardt.argtypes = [
        ctypes.c_void_p,
        ctypes.c_double,
        ctypes.c_int,
        ctypes.c_int,
        ctypes.POINTER(ctypes.c_double),
        ctypes.c_size_t,
        ctypes.POINTER(ctypes.c_double),
        ctypes.c_size_t,
    ]
    lib.seams_steinhardt.restype = ctypes.c_int
    lib.seams_chill.argtypes = [
        ctypes.c_void_p,
        ctypes.c_double,
        ctypes.c_int,
        ctypes.c_int,
        ctypes.c_int,
        ctypes.POINTER(ctypes.c_int),
        ctypes.c_size_t,
        ctypes.POINTER(ctypes.c_int),
    ]
    lib.seams_chill.restype = ctypes.c_int
    lib.seams_chill_plus.argtypes = [
        ctypes.c_void_p,
        ctypes.c_double,
        ctypes.c_int,
        ctypes.c_int,
        ctypes.c_int,
        ctypes.POINTER(ctypes.c_int),
        ctypes.POINTER(ctypes.c_int),
    ]
    lib.seams_chill_plus.restype = ctypes.c_int
    lib.seams_last_error.restype = ctypes.c_char_p
    return lib


def park_openmp() -> None:
    """Stop a spinning OpenMP team so the next library gets the cores."""
    global _PAUSE
    if _PAUSE is False:
        return
    if _PAUSE is None:
        try:
            lib = ctypes.CDLL("libgomp.so.1")
            fn = lib.omp_pause_resource_all
            fn.argtypes = [ctypes.c_int]
            fn.restype = None
            _PAUSE = fn
        except (OSError, AttributeError):
            _PAUSE = False
            return
    _PAUSE(1)


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
    pos = np.asarray(coords, dtype=np.float64)
    return pos, (lx, ly, lz)


def undirected(left: np.ndarray, right: np.ndarray) -> np.ndarray:
    a = np.minimum(left, right)
    b = np.maximum(left, right)
    mask = a != b
    edges = np.stack((a[mask], b[mask]), axis=1)
    if edges.shape[0] == 0:
        return edges.astype(np.int64, copy=False)
    order = np.lexsort((edges[:, 1], edges[:, 0]))
    edges = edges[order]
    keep = np.empty(edges.shape[0], dtype=bool)
    keep[0] = True
    keep[1:] = np.any(edges[1:] != edges[:-1], axis=1)
    return edges[keep]


def pairs_from_csr(offsets: np.ndarray, cols: np.ndarray) -> np.ndarray:
    counts = np.diff(offsets)
    sources = np.repeat(np.arange(counts.shape[0], dtype=np.int64), counts)
    return undirected(sources, np.asarray(cols, dtype=np.int64))


def knn_matrix(offsets: np.ndarray, cols: np.ndarray, k: int):
    degree = np.diff(offsets.astype(np.int64, copy=False))
    if degree.size == 0 or not np.all(degree == k):
        return None
    matrix = np.asarray(cols, dtype=np.int64).reshape(degree.shape[0], k)
    matrix.sort(axis=1)
    return matrix


def samples_of(fn, n=SAMPLES) -> list[float]:
    fn()
    out = []
    for _ in range(n):
        t0 = time.perf_counter_ns()
        fn()
        out.append((time.perf_counter_ns() - t0) / 1e6)
    return out


def emit(lines: list[str], echo, text: str) -> None:
    lines.append(text)
    if echo:
        print(text, flush=True)


def emit_time(lines, echo, molecules, code, op, threads, got) -> None:
    raw = " ".join(f"{v:.4f}" for v in got)
    emit(
        lines, echo,
        f"arm {code} op {op} threads {threads} molecules {molecules} "
        f"samples_ms {raw}",
    )
    emit(
        lines, echo,
        f"arm {code} op {op} threads {threads} molecules {molecules} "
        f"median_ms {statistics.median(got):.4f}",
    )


def require(lib, status: int, what: str) -> None:
    if status == 0:
        return
    err = lib.seams_last_error()
    detail = err.decode() if err else "unknown"
    raise SystemExit(f"d-SEAMS {what} failed: {detail}")


def copy_csr(lib, nb) -> tuple[np.ndarray, np.ndarray]:
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


def dseams_cutoff(lib, sys):
    nb = ctypes.c_void_p()
    require(lib, lib.seams_neighbours_build(sys, CUTOFF, 1, SEAMS_GRAPH_CUTOFF, ctypes.byref(nb)), "cutoff")
    return copy_csr(lib, nb)


def dseams_knn(lib, sys):
    nb = ctypes.c_void_p()
    require(
        lib,
        lib.seams_neighbours_knn(sys, KNN, CUTOFF, 1, 0, ctypes.byref(nb)),
        "knn",
    )
    return copy_csr(lib, nb)


def dseams_q6(lib, sys, q, bar):
    n = len(q)
    require(lib, lib.seams_steinhardt(sys, CUTOFF, 1, 6, q, n, bar, n), "q6")


def dseams_chill(lib, sys, labels, counts):
    n = len(labels)
    require(
        lib,
        lib.seams_chill(
            sys, CUTOFF, 1, SEAMS_GRAPH_CUTOFF, SEAMS_CHILL_PUBLISHED,
            labels, n, counts,
        ),
        "chill",
    )


def dseams_chill_plus(lib, sys, labels, counts, rule=SEAMS_CHILL_PUBLISHED):
    require(
        lib,
        lib.seams_chill_plus(
            sys, CUTOFF, 1, SEAMS_GRAPH_CUTOFF, rule,
            labels, counts,
        ),
        "chill+",
    )


def ovito_names(pos, lengths):
    from ovito.data import DataCollection, Particles, SimulationCell
    from ovito.modifiers import ChillPlusModifier
    from ovito.pipeline import Pipeline, StaticSource

    lx, ly, lz = lengths
    data = DataCollection()
    cell = SimulationCell(pbc=(True, True, True))
    matrix = np.zeros((3, 4))
    matrix[0, 0], matrix[1, 1], matrix[2, 2] = lx, ly, lz
    cell[...] = matrix
    data.objects.append(cell)
    particles = Particles()
    particles.create_property("Position", data=np.ascontiguousarray(pos, dtype=np.float64))
    data.objects.append(particles)
    pipe = Pipeline(source=StaticSource(data=data))
    pipe.modifiers.append(ChillPlusModifier(cutoff=CUTOFF))
    out = pipe.compute()
    codes = np.asarray(out.particles["Structure Type"], dtype=np.int64)
    return np.array([OVITO_CODE[int(v)] for v in codes])


def require_atom_match(left, right, what: str) -> None:
    if left.shape != right.shape or not np.array_equal(left, right):
        n = int(np.count_nonzero(left != right)) if left.shape == right.shape else -1
        raise SystemExit(f"{what}: {n} atoms differ")


def histogram(labels: np.ndarray, names: tuple[str, ...]) -> str:
    counts = np.bincount(labels.astype(np.int64), minlength=len(names))
    parts = [f"{name}={int(counts[i])}" for i, name in enumerate(names) if counts[i]]
    return " ".join(parts) if parts else "empty"


def max_abs(a: np.ndarray, b: np.ndarray) -> float:
    return float(np.max(np.abs(a - b))) if a.size else 0.0


def vesin_pairs(pos, box, threads: int) -> np.ndarray:
    calc = vesin.NeighborList(cutoff=CUTOFF, full_list=True, n_threads=threads)
    left, right = calc.compute(points=pos, box=box, periodic=True, quantities="ij", copy=False)
    return undirected(np.asarray(left, dtype=np.int64), np.asarray(right, dtype=np.int64))


def mdapy_cutoff(system, threads: int):
    mdapy.set_num_threads(threads)
    neigh = mdapy.NeighborList.build(system, cutoff=CUTOFF)
    return pairs_from_csr(np.asarray(neigh.offsets), np.asarray(neigh.j))


def mdapy_knn(system, threads: int):
    mdapy.set_num_threads(threads)
    neigh = mdapy.NeighborList.build(system, num_neighbors=KNN)
    return knn_matrix(np.asarray(neigh.offsets), np.asarray(neigh.j), KNN)


def pyscal_neighbors(atoms, threads: int, method: str):
    pyscal3.set_num_threads(threads)
    if method == "cutoff":
        pyscal3.find_neighbors(atoms, method="cutoff", cutoff=CUTOFF, store_rows=False)
    else:
        pyscal3.find_neighbors(atoms, method="number", nmax=KNN, store_rows=False)
    offsets, values = neighbor_arrays(atoms, "neighbors")
    return np.asarray(offsets), np.asarray(values["neighbors"])


def same_edges(a: np.ndarray, b: np.ndarray) -> bool:
    return a.shape == b.shape and np.array_equal(a, b)


def mic_distances(pos: np.ndarray, lengths, edges: np.ndarray) -> np.ndarray:
    if edges.size == 0:
        return np.empty(0)
    delta = pos[edges[:, 1]] - pos[edges[:, 0]]
    box = np.asarray(lengths, dtype=np.float64)
    delta -= box * np.round(delta / box)
    return np.linalg.norm(delta, axis=1)


def explain_pairs(lines, echo, pos, lengths, reference, other, other_name: str) -> None:
    if same_edges(reference, other):
        emit(lines, echo, f"check cutoff {other_name} pairs {reference.shape[0]} match exact")
        return
    ref = {tuple(row) for row in reference}
    got = {tuple(row) for row in other}
    only_ref = np.asarray(sorted(ref - got), dtype=np.int64).reshape(-1, 2) if ref - got else np.empty((0, 2), np.int64)
    only_got = np.asarray(sorted(got - ref), dtype=np.int64).reshape(-1, 2) if got - ref else np.empty((0, 2), np.int64)
    d_ref = mic_distances(pos, lengths, only_ref) if only_ref.size else np.empty(0)
    d_got = mic_distances(pos, lengths, only_got) if only_got.size else np.empty(0)

    def span(dist: np.ndarray) -> str:
        if dist.size == 0:
            return "none"
        return f"n={dist.size} min={dist.min():.6f} max={dist.max():.6f}"

    emit(
        lines, echo,
        f"check cutoff {other_name} pairs dseams {reference.shape[0]} "
        f"other {other.shape[0]} only_dseams {span(d_ref)} only_other {span(d_got)}",
    )
    near = 0
    if d_got.size:
        near += int(np.count_nonzero(np.abs(d_got - CUTOFF) <= 1e-8))
    if d_ref.size:
        near += int(np.count_nonzero(np.abs(d_ref - CUTOFF) <= 1e-8))
    if near:
        emit(
            lines, echo,
            f"note cutoff {other_name} disagreement is {near} pairs sitting on {CUTOFF} A; "
            "d-SEAMS keeps a pair only when the minimum image is strictly inside the cutoff",
        )
    else:
        emit(
            lines, echo,
            f"note cutoff {other_name} disagreement is not a boundary pair at {CUTOFF} A",
        )


def explain_q(lines, echo, kind: str, tool: str, ref: np.ndarray, got: np.ndarray, same_shell: bool) -> None:
    diff = max_abs(ref, got)
    emit(lines, echo, f"check {kind} {tool} max_abs {diff:.3e}")
    if diff <= Q_ATOL:
        return
    if not same_shell:
        emit(
            lines, echo,
            f"note {kind} {tool} differs because the neighbour shell is not the same set of atoms",
        )
        return
    scale = float(np.median(got / ref)) if np.all(ref > 0) else float("nan")
    emit(
        lines, echo,
        f"note {kind} {tool} uses the same neighbours and still differs; "
        f"median ratio other/dseams {scale:.6f}. "
        "d-SEAMS q_l is sqrt(4*pi/(2l+1) * sum |q_lm|^2) and q6bar averages q_lm "
        "over the atom and its neighbours before that norm",
    )


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


def check_mixed(lib, lines, echo) -> None:
    """Ovito mode matches mdapy and the OVITO modifier on the mixed frames."""
    import importlib.util

    spec = importlib.util.spec_from_file_location(
        "mixed_phase_frames", ROOT / "scripts" / "mixed_phase_frames.py"
    )
    mixed = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(mixed)
    mdapy.set_num_threads(1)
    for name, builder in mixed.FRAMES.items():
        pos, lengths = builder()
        n = int(pos.shape[0])
        sys = frame_system(lib, pos, lengths)
        labels = (ctypes.c_int * n)()
        counts = (ctypes.c_int * 9)()
        require(lib, lib.seams_set_num_threads(1), "threads")
        dseams_chill_plus(lib, sys, labels, counts, SEAMS_CHILL_OVITO)
        one = np.ctypeslib.as_array(labels).copy()
        require(lib, lib.seams_set_num_threads(8), "threads")
        dseams_chill_plus(lib, sys, labels, counts, SEAMS_CHILL_OVITO)
        if not np.array_equal(one, np.ctypeslib.as_array(labels)):
            raise SystemExit(f"mixed {name}: ovito mode changed with the thread count")
        dseams_chill_plus(lib, sys, labels, counts, SEAMS_CHILL_PUBLISHED)
        paper8 = np.ctypeslib.as_array(labels).copy()
        require(lib, lib.seams_set_num_threads(1), "threads")
        dseams_chill_plus(lib, sys, labels, counts, SEAMS_CHILL_PUBLISHED)
        if not np.array_equal(paper8, np.ctypeslib.as_array(labels)):
            raise SystemExit(f"mixed {name}: paper mode changed with the thread count")
        ours = np.array([DSEAMS_AS_OVITO[int(v)] for v in one])
        system = mdapy.System(pos, box=lengths)
        ice = mdapy.identify.ChillPlus(cutoff=CUTOFF).compute(system)
        md = np.array([tuple(ice.names)[int(v)] for v in np.asarray(ice.structure)])
        real = ovito_names(pos, lengths)
        require_atom_match(ours, md, f"mixed {name} d-SEAMS ovito vs mdapy")
        require_atom_match(ours, real, f"mixed {name} d-SEAMS ovito vs OVITO")
        lib.seams_system_destroy(sys)
        emit(lines, echo, f"check mixed {name} ovito mode matches mdapy and OVITO on {n} oxygens")


def run_session(echo: bool = False) -> list[str]:
    lib = load_seams(library_path())
    lines: list[str] = []
    emit(lines, echo, f"cores {core_count()}")
    emit(lines, echo, f"samples {SAMPLES}")
    emit(lines, echo, "cell orthorhombic ice Ih oxygen molecules")
    emit(lines, echo, f"cutoff_A {CUTOFF}")
    emit(lines, echo, f"knn {KNN}")
    emit(
        lines, echo,
        "note vesin builds a cutoff list and has no k-nearest operation",
    )
    emit(
        lines, echo,
        "note pyscal3 has no CHILL; its ice identifier here is extended CNA diamond structure",
    )
    emit(
        lines, echo,
        "note mdapy and OVITO CHILL+ are compared to d-SEAMS ovito mode, which sums q_lm over every neighbour inside 3.5 A",
    )
    emit(
        lines, echo,
        "note d-SEAMS paper mode sums q_lm over the four closest neighbours and is timed separately",
    )
    emit(
        lines, echo,
        "note an OpenMP team is paused before mdapy, pyscal3, and vesin so it does not stay in active wait",
    )

    check_mixed(lib, lines, echo)

    for ncell in SIZES:
        pos, lengths = ice_positions(ncell)
        n = int(pos.shape[0])
        lx, ly, lz = lengths
        emit(lines, echo, f"molecules {n}")
        box = np.diag(np.asarray(lengths, dtype=np.float64))
        system = mdapy.System(pos, box=lengths)
        atoms = Atoms(positions=pos, cell=lengths, pbc=True)
        box9 = (ctypes.c_double * 9)(0, 0, 0, lx, ly, lz, 0, 0, 0)
        types = (ctypes.c_int * n)(*([1] * n))
        ids = (ctypes.c_int * n)(*range(1, n + 1))
        ptr = pos.ctypes.data_as(ctypes.POINTER(ctypes.c_double))
        sys = lib.seams_system_create()
        require(lib, lib.seams_system_set_frame(sys, ptr, n, box9, types, ids), "set_frame")
        require(lib, lib.seams_set_num_threads(1), "threads")

        d_off, d_col = dseams_cutoff(lib, sys)
        d_pairs = pairs_from_csr(d_off, d_col)
        emit(lines, echo, f"check cutoff dseams pairs {d_pairs.shape[0]}")
        park_openmp()
        explain_pairs(lines, echo, pos, lengths, d_pairs, vesin_pairs(pos, box, 1), "vesin")
        explain_pairs(lines, echo, pos, lengths, d_pairs, mdapy_cutoff(system, 1), "mdapy")
        p_off, p_col = pyscal_neighbors(atoms, 1, "cutoff")
        explain_pairs(lines, echo, pos, lengths, d_pairs, pairs_from_csr(p_off, p_col), "pyscal3")

        require(lib, lib.seams_set_num_threads(1), "threads")
        k_off, k_col = dseams_knn(lib, sys)
        d_knn = knn_matrix(k_off, k_col, KNN)
        if d_knn is None:
            degree = np.diff(k_off)
            emit(
                lines, echo,
                f"check knn dseams degree min {int(degree.min())} max {int(degree.max())} "
                "not a constant k",
            )
        else:
            emit(lines, echo, f"check knn dseams rows {d_knn.shape[0]} k {KNN}")
        park_openmp()
        m_knn = mdapy_knn(system, 1)
        k_off_p, k_col_p = pyscal_neighbors(atoms, 1, "number")
        p_knn = knn_matrix(k_off_p, k_col_p, KNN)
        if d_knn is not None and m_knn is not None and same_edges(d_knn, m_knn):
            emit(lines, echo, "check knn mdapy match exact")
        elif d_knn is None or m_knn is None:
            emit(lines, echo, "check knn mdapy degree differs from k")
        else:
            mismatch = int(np.count_nonzero(np.any(d_knn != m_knn, axis=1)))
            emit(lines, echo, f"check knn mdapy atoms_differ {mismatch}")
            emit(
                lines, echo,
                "note knn mdapy is the k nearest in the periodic box with ties broken by atom index; "
                "d-SEAMS takes k nearest inside cutoff + 1.5 A",
            )
        if d_knn is not None and p_knn is not None and same_edges(d_knn, p_knn):
            emit(lines, echo, "check knn pyscal3 match exact")
        elif d_knn is None or p_knn is None:
            emit(lines, echo, "check knn pyscal3 degree differs from k")
        else:
            mismatch = int(np.count_nonzero(np.any(d_knn != p_knn, axis=1)))
            emit(lines, echo, f"check knn pyscal3 atoms_differ {mismatch}")
            emit(
                lines, echo,
                "note knn pyscal3 method number keeps nmax nearest neighbours; "
                "d-SEAMS takes k nearest inside cutoff + 1.5 A",
            )
        shell_pairs = pairs_from_csr(k_off, k_col)
        if same_edges(d_pairs, shell_pairs):
            emit(
                lines, echo,
                "note on this ice the four nearest neighbours are exactly the 3.5 A cutoff shell",
            )
        else:
            emit(
                lines, echo,
                f"note knn undirected pairs {shell_pairs.shape[0]} differ from cutoff pairs {d_pairs.shape[0]}",
            )

        require(lib, lib.seams_set_num_threads(1), "threads")
        q = (ctypes.c_double * n)()
        bar = (ctypes.c_double * n)()
        dseams_q6(lib, sys, q, bar)
        q_np = np.ctypeslib.as_array(q).copy()
        bar_np = np.ctypeslib.as_array(bar).copy()
        park_openmp()
        mdapy.set_num_threads(1)
        shell = mdapy.NeighborList.build(system, cutoff=CUTOFF)
        local = mdapy.analysis.SteinhardtParameters((6,), cutoff=CUTOFF, average=False).compute(
            system, neighbors=shell
        )
        averaged = mdapy.analysis.SteinhardtParameters((6,), cutoff=CUTOFF, average=True).compute(
            system, neighbors=shell
        )
        mq = np.asarray(local.q, dtype=np.float64).reshape(n, -1)[:, 0]
        mbar = np.asarray(averaged.q, dtype=np.float64).reshape(n, -1)[:, 0]
        same_mdapy_shell = same_edges(d_pairs, pairs_from_csr(np.asarray(shell.offsets), np.asarray(shell.j)))
        explain_q(lines, echo, "q6", "mdapy", q_np, mq, same_mdapy_shell)
        explain_q(lines, echo, "q6bar", "mdapy", bar_np, mbar, same_mdapy_shell)

        pyscal3.set_num_threads(1)
        pyscal3.find_neighbors(atoms, method="cutoff", cutoff=CUTOFF, store_rows=False)
        pyscal3.steinhardt_parameter(atoms, 6, averaged=True)
        pq = np.asarray(atoms.arrays["pyscal_q6"], dtype=np.float64)
        pbar = np.asarray(atoms.arrays["pyscal_avg_q6"], dtype=np.float64)
        p_off, p_col = neighbor_arrays(atoms, "neighbors")
        same_pyscal_shell = same_edges(d_pairs, pairs_from_csr(np.asarray(p_off), np.asarray(p_col["neighbors"])))
        explain_q(lines, echo, "q6", "pyscal3", q_np, pq, same_pyscal_shell)
        explain_q(lines, echo, "q6bar", "pyscal3", bar_np, pbar, same_pyscal_shell)

        require(lib, lib.seams_set_num_threads(1), "threads")
        labels = (ctypes.c_int * n)()
        counts = (ctypes.c_int * 9)()
        dseams_chill(lib, sys, labels, counts)
        chill_labels = np.ctypeslib.as_array(labels).copy()
        emit(lines, echo, "check chill dseams " + histogram(chill_labels, DSEAMS_ICE))
        dseams_chill_plus(lib, sys, labels, counts, SEAMS_CHILL_PUBLISHED)
        plus_labels = np.ctypeslib.as_array(labels).copy()
        emit(lines, echo, "check chill_plus paper " + histogram(plus_labels, DSEAMS_ICE))
        dseams_chill_plus(lib, sys, labels, counts, SEAMS_CHILL_OVITO)
        ovito_mode = np.ctypeslib.as_array(labels).copy()
        emit(lines, echo, "check chill_plus ovito " + histogram(ovito_mode, DSEAMS_ICE))
        park_openmp()
        mdapy.set_num_threads(1)
        ice = mdapy.identify.ChillPlus(cutoff=CUTOFF).compute(system)
        ice_labels = np.asarray(ice.structure, dtype=np.int64)
        emit(lines, echo, "check ice mdapy " + histogram(ice_labels, tuple(ice.names)))
        ours = np.array([DSEAMS_AS_OVITO[int(v)] for v in ovito_mode])
        md_names = np.array([tuple(ice.names)[int(v)] for v in ice_labels])
        real = ovito_names(pos, lengths)
        require_atom_match(ours, md_names, f"ice {n} d-SEAMS ovito vs mdapy")
        require_atom_match(ours, real, f"ice {n} d-SEAMS ovito vs OVITO")
        emit(lines, echo, f"check ice {n} ovito mode matches mdapy and OVITO on every oxygen")
        d_hex = int(np.count_nonzero(chill_labels == 1))
        p_hex = int(np.count_nonzero(plus_labels == 1))
        m_hex = int(np.count_nonzero(ice_labels == 1))
        if d_hex == n and p_hex == n and m_hex == n:
            emit(
                lines, echo,
                "note published CHILL, published CHILL+, and mdapy CHILL+ each label every oxygen hexagonal on this ice",
            )
        else:
            emit(
                lines, echo,
                f"note ice labels differ: dseams CHILL hexagonal {d_hex}, "
                f"dseams CHILL+ hexagonal {p_hex}, mdapy CHILL+ hexagonal {m_hex} of {n}. "
                "Thresholds were not changed. CHILL is Moore's four-neighbour correlation; "
                "CHILL+ adds the interfacial and hydrate classes",
            )
        pyscal3.set_num_threads(1)
        diamond = pyscal3.diamond_structure(atoms)
        parts = " ".join(f"{key.replace(' ', '_')}={value}" for key, value in diamond.items() if value)
        emit(lines, echo, "check ice pyscal3 " + parts)
        hex_d = int(diamond.get("hex diamond", 0))
        if hex_d != n:
            emit(
                lines, echo,
                "note pyscal3 diamond_structure is extended CNA with a fixed second argument of 2.0 "
                "and is not CHILL; a perfect ice Ih cell is not required to match its hex-diamond count",
            )

        require(lib, lib.seams_set_num_threads(8), "threads")
        off8, col8 = dseams_cutoff(lib, sys)
        if not same_edges(d_pairs, pairs_from_csr(off8, col8)):
            raise SystemExit("d-SEAMS cutoff pairs changed with the thread count")
        k8_off, k8_col = dseams_knn(lib, sys)
        k8 = knn_matrix(k8_off, k8_col, KNN)
        if d_knn is None or k8 is None or not same_edges(d_knn, k8):
            raise SystemExit("d-SEAMS knn changed with the thread count")
        q8 = (ctypes.c_double * n)()
        bar8 = (ctypes.c_double * n)()
        dseams_q6(lib, sys, q8, bar8)
        if max_abs(q_np, np.ctypeslib.as_array(q8)) > 0.0 or max_abs(bar_np, np.ctypeslib.as_array(bar8)) > 0.0:
            raise SystemExit("d-SEAMS q6 changed with the thread count")
        labels8 = (ctypes.c_int * n)()
        dseams_chill(lib, sys, labels8, counts)
        if not np.array_equal(chill_labels, np.ctypeslib.as_array(labels8)):
            raise SystemExit("d-SEAMS CHILL changed with the thread count")
        dseams_chill_plus(lib, sys, labels8, counts, SEAMS_CHILL_PUBLISHED)
        if not np.array_equal(plus_labels, np.ctypeslib.as_array(labels8)):
            raise SystemExit("d-SEAMS CHILL+ paper mode changed with the thread count")
        dseams_chill_plus(lib, sys, labels8, counts, SEAMS_CHILL_OVITO)
        if not np.array_equal(ovito_mode, np.ctypeslib.as_array(labels8)):
            raise SystemExit("d-SEAMS CHILL+ ovito mode changed with the thread count")
        emit(lines, echo, "check dseams threads 1 and 8 agree for paper and ovito")

        q_buf = (ctypes.c_double * n)()
        bar_buf = (ctypes.c_double * n)()
        lab_buf = (ctypes.c_int * n)()
        cnt_buf = (ctypes.c_int * 9)()

        def time_dseams(width: int) -> None:
            require(lib, lib.seams_set_num_threads(width), "threads")

            def cutoff():
                nb = ctypes.c_void_p()
                require(
                    lib,
                    lib.seams_neighbours_build(sys, CUTOFF, 1, SEAMS_GRAPH_CUTOFF, ctypes.byref(nb)),
                    "cutoff",
                )
                if lib.seams_neighbours_copy(nb, None, 0, ctypes.byref(ctypes.c_size_t()), None, 0, ctypes.byref(ctypes.c_size_t())) != 0:
                    raise SystemExit("d-SEAMS cutoff length")
                lib.seams_neighbours_destroy(nb)

            def knn():
                nb = ctypes.c_void_p()
                require(lib, lib.seams_neighbours_knn(sys, KNN, CUTOFF, 1, 0, ctypes.byref(nb)), "knn")
                lib.seams_neighbours_destroy(nb)

            def q6():
                dseams_q6(lib, sys, q_buf, bar_buf)

            def chill():
                dseams_chill(lib, sys, lab_buf, cnt_buf)

            def chill_plus():
                dseams_chill_plus(lib, sys, lab_buf, cnt_buf, SEAMS_CHILL_PUBLISHED)

            def chill_plus_ovito():
                dseams_chill_plus(lib, sys, lab_buf, cnt_buf, SEAMS_CHILL_OVITO)

            emit_time(lines, echo, n, "dseams", "cutoff", width, samples_of(cutoff))
            emit_time(lines, echo, n, "dseams", "knn", width, samples_of(knn))
            emit_time(lines, echo, n, "dseams", "q6", width, samples_of(q6))
            emit_time(lines, echo, n, "dseams", "chill", width, samples_of(chill))
            emit_time(lines, echo, n, "dseams", "chill_plus_paper", width, samples_of(chill_plus))
            emit_time(lines, echo, n, "dseams", "chill_plus_ovito", width, samples_of(chill_plus_ovito))

        time_dseams(1)
        time_dseams(8)

        def time_others(width: int) -> None:
            park_openmp()
            calc = vesin.NeighborList(cutoff=CUTOFF, full_list=True, n_threads=width)

            def vesin_cutoff():
                left, right = calc.compute(points=pos, box=box, periodic=True, quantities="ij", copy=False)
                if left.shape[0] < n:
                    raise SystemExit("vesin cutoff short")

            emit_time(lines, echo, n, "vesin", "cutoff", width, samples_of(vesin_cutoff))

            mdapy.set_num_threads(width)

            def mdapy_cut():
                neigh = mdapy.NeighborList.build(system, cutoff=CUTOFF)
                if int(neigh.offsets.shape[0]) != n + 1:
                    raise SystemExit("mdapy cutoff shape")

            def mdapy_k():
                neigh = mdapy.NeighborList.build(system, num_neighbors=KNN)
                if int(neigh.j.shape[0]) != n * KNN:
                    raise SystemExit("mdapy knn shape")

            def mdapy_q():
                neigh = mdapy.NeighborList.build(system, cutoff=CUTOFF)
                local_q = mdapy.analysis.SteinhardtParameters((6,), cutoff=CUTOFF, average=False).compute(
                    system, neighbors=neigh
                )
                avg_q = mdapy.analysis.SteinhardtParameters((6,), cutoff=CUTOFF, average=True).compute(
                    system, neighbors=neigh
                )
                if local_q.q.shape[0] != n or avg_q.q.shape[0] != n:
                    raise SystemExit("mdapy q6 shape")

            def mdapy_ice():
                result = mdapy.identify.ChillPlus(cutoff=CUTOFF).compute(system)
                if int(result.structure.shape[0]) != n:
                    raise SystemExit("mdapy ice shape")

            emit_time(lines, echo, n, "mdapy", "cutoff", width, samples_of(mdapy_cut))
            emit_time(lines, echo, n, "mdapy", "knn", width, samples_of(mdapy_k))
            emit_time(lines, echo, n, "mdapy", "q6", width, samples_of(mdapy_q))
            emit_time(lines, echo, n, "mdapy", "chill_plus", width, samples_of(mdapy_ice))

            pyscal3.set_num_threads(width)

            def pyscal_cut():
                pyscal3.find_neighbors(atoms, method="cutoff", cutoff=CUTOFF, store_rows=False)
                if "pyscal_bond_offsets" not in atoms.info and "pyscal_neighbors" not in atoms.arrays:
                    raise SystemExit("pyscal3 cutoff empty")

            def pyscal_k():
                pyscal3.find_neighbors(atoms, method="number", nmax=KNN, store_rows=False)
                offsets, values = neighbor_arrays(atoms, "neighbors")
                if int(np.asarray(values["neighbors"]).shape[0]) != n * KNN:
                    raise SystemExit("pyscal3 knn shape")

            def pyscal_q():
                pyscal3.find_neighbors(atoms, method="cutoff", cutoff=CUTOFF, store_rows=False)
                got = pyscal3.steinhardt_parameter(atoms, 6, averaged=True)
                if got[0].shape[0] != n:
                    raise SystemExit("pyscal3 q6 shape")

            def pyscal_ice():
                got = pyscal3.diamond_structure(atoms)
                if sum(got.values()) != n:
                    raise SystemExit("pyscal3 diamond count")

            emit_time(lines, echo, n, "pyscal3", "cutoff", width, samples_of(pyscal_cut))
            emit_time(lines, echo, n, "pyscal3", "knn", width, samples_of(pyscal_k))
            emit_time(lines, echo, n, "pyscal3", "q6", width, samples_of(pyscal_q))
            emit_time(lines, echo, n, "pyscal3", "diamond", width, samples_of(pyscal_ice))

        time_others(1)
        time_others(8)
        lib.seams_system_destroy(sys)
    return lines


def main() -> None:
    run_session(echo=True)


if __name__ == "__main__":
    main()
