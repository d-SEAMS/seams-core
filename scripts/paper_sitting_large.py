#!/usr/bin/env python3
"""One sitting at 884736 molecules, eight times the 110592 ice Ic cell.

The operations are the ice rows of sitting.csv plus the symmetry
classifier and both PTM template sets: rings, local q4, local q6, the
q4+q6 library call, the cutoff neighbour list, CHILL+, the symmetry
classifier, PTM with its default templates, and PTM with only the
cubic and hexagonal diamond templates. Rivals are
rings-statistics-matter 0.3.0, mdapy, vesin, freud and OVITO.

One warm-up, then eight calls. A rival call over 30 s is kept at six
calls and marked. A rival call over 120 s, or a process killed for
memory, is recorded as that cell and is not replaced with an estimate.
Each width is its own process, pinned with taskset. The median of an
even count is the mean of the two central samples. The interquartile
range is type-7. The host name is not read.
"""

from __future__ import annotations

import argparse
import csv
import json
import os
import re
import subprocess
import sys
import threading
import time
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "scripts"))
import paper_sitting as sit
import paper_stats as stats

MOLECULES = 884736
CELLS = 48
SIX = 1769472
CUTOFF = 3.5
WARMUP = 1
REPEATS = 8
SLOW_REPEATS = 6
SLOW_MS = 30_000.0
LIMIT_S = 120.0
WIDTHS = (1, 2, 4, 8)

FIELDS = [
    "source_commit", "cpu_model", "cpu_cores", "affinity",
    "phase", "operation", "code", "version",
    "molecules", "threads", "warmup", "n",
    "median_ms", "iqr_ms", "q1_ms", "q3_ms", "samples_ms", "detail",
]

# code key, published name, rival (the 30 s / 120 s rules apply)
CODES = (
    ("dseams", "d-SEAMS", False),
    ("rsm", "rings-statistics-matter", True),
    ("mdapy", "mdapy", True),
    ("vesin", "vesin", True),
    ("freud", "freud", True),
    ("ovito", "OVITO", True),
)

LSCPU_FIELDS = (
    "Architecture",
    "CPU(s)",
    "On-line CPU(s) list",
    "Vendor ID",
    "Model name",
    "CPU family",
    "Model",
    "Thread(s) per core",
    "Core(s) per socket",
    "Socket(s)",
    "BogoMIPS",
    "L1d cache",
    "L1i cache",
    "L2 cache",
    "L3 cache",
)


def emit(payload: dict) -> None:
    sys.stdout.write(json.dumps(payload) + "\n")
    sys.stdout.flush()


def skip(operation: str) -> bool:
    """True when a restarted worker has not yet reached --from-op."""
    allow = os.environ.get("SITTING_ONLY", "")
    if allow and operation not in set(allow.split(",")):
        return True
    start = os.environ.get("SITTING_FROM", "")
    if not start or getattr(skip, "open", False):
        return False
    if operation == start:
        skip.open = True
        return False
    return True


def project_version() -> str:
    text = (ROOT / "meson.build").read_text(encoding="utf-8")
    match = re.search(r"version:\s*'([0-9]+\.[0-9]+\.[0-9]+)'", text)
    if not match:
        return "unrecorded"
    return f"seams {match.group(1)}"


def lscpu_details() -> list[str]:
    text = subprocess.check_output(["lscpu"], text=True)
    wanted = {label: "" for label in LSCPU_FIELDS}
    for line in text.splitlines():
        if ":" not in line:
            continue
        label, value = line.split(":", 1)
        label = label.strip()
        if label in wanted:
            wanted[label] = value.strip()
    return [f"{label}: {wanted[label]}" for label in LSCPU_FIELDS]


def ice():
    import h2h_protocol as proto
    pos, side = proto.ice_ic(CELLS)
    if len(pos) != MOLECULES:
        raise RuntimeError(f"ice Ic has {len(pos)} molecules, need {MOLECULES}")
    import numpy as np
    xyz = np.ascontiguousarray(pos, dtype=np.float64)
    box = np.array([0, 0, 0, side, side, side, 0, 0, 0], dtype=np.float64)
    return xyz, box, float(side)


def hold_affinity(threads: int) -> int:
    """Put the process back on the taskset cores.

    OMP_PROC_BIND narrows the main thread to one core after a runtime
    loads. Rayon, TBB and mdapy then inherit that mask and stay serial.
    """
    allowed = set(range(threads))
    try:
        os.sched_setaffinity(0, allowed)
    except (AttributeError, OSError):
        return threads
    return len(os.sched_getaffinity(0))


def run_samples(call, rival: bool, threads: int):
    """Return (samples, last, note). A rival over 30 s stops at six calls."""
    repeats = REPEATS
    note = ""
    samples = []
    last = None
    index = 0
    while len(samples) < repeats:
        if rival:
            pinned = hold_affinity(threads)
            if pinned != threads:
                raise RuntimeError(f"affinity has {pinned} cores, width is {threads}")
        emit({"event": "sample-begin"})
        t0 = time.perf_counter()
        last = call()
        elapsed = (time.perf_counter() - t0) * 1000.0
        emit({"event": "sample-end", "elapsed_ms": elapsed})
        if rival and elapsed > LIMIT_S * 1000.0:
            raise TimeoutError("exceeded 120 s")
        if index >= WARMUP:
            samples.append(elapsed)
        if rival and elapsed > SLOW_MS:
            repeats = SLOW_REPEATS
            note = "n=6 because a call exceeded 30 s"
            if len(samples) > SLOW_REPEATS:
                del samples[SLOW_REPEATS:]
        index += 1
        if len(samples) >= repeats:
            break
        if index > WARMUP + REPEATS + 2:
            raise RuntimeError("sample loop did not finish")
    return samples, last, note


def save_agree(threads: int, code: str, operation: str, values) -> None:
    root = os.environ.get("SIT_AGREE")
    if not root:
        return
    import numpy as np
    path = Path(root) / f"t{threads}-{code}-{operation}.npy"
    path.parent.mkdir(parents=True, exist_ok=True)
    np.save(path, np.asarray(values))


def time_dseams(lib, core, xyz, box, version: str, threads: int) -> None:
    import ctypes
    import numpy as np
    lib.seams_frame_neighbour_pairs.restype = ctypes.c_int
    lib.seams_frame_neighbour_pairs.argtypes = [ctypes.c_void_p]
    handle, n = sit.frame_from(lib, xyz, box)
    try:
        counts = (ctypes.c_int * 7)()

        def rings():
            sit.fail_lang(lib, lib.seams_frame_rings(
                handle, CUTOFF, 1, b"cutoff", 6, counts, 7), "rings")
            return int(counts[6])

        if not skip("rings"):
            emit({"event": "begin", "operation": "rings"})
            samples, six, note = run_samples(rings, False, 0)
            if six != SIX:
                raise RuntimeError(f"d-SEAMS six-rings {six} != {SIX}")
            _done("rings", "d-SEAMS", version, n, samples, f"six={six}", note)

        q4 = (ctypes.c_double * n)()
        q6 = (ctypes.c_double * n)()

        def local(order, dest):
            def run():
                sit.fail_lang(lib, lib.seams_frame_steinhardt(
                    handle, CUTOFF, 1, order, dest, None), f"q{order}")
                return float(dest[0])
            return run

        for order, dest, name in ((4, q4, "q4"), (6, q6, "q6")):
            if skip(name):
                continue
            emit({"event": "begin", "operation": name})
            samples, _probe, note = run_samples(local(order, dest), False, 0)
            if name == "q6":
                save_agree(threads, "d-SEAMS", "q6", np.ctypeslib.as_array(dest))
            _done(name, "d-SEAMS", version, n, samples,
                  "seams_frame_steinhardt; local q only", note)

        def both():
            sit.fail_lang(lib, lib.seams_frame_q4_q6(
                handle, CUTOFF, 1, q4, q6, None, None), "q4+q6")
            return float(q4[0]) + float(q6[0])

        if not skip("q4+q6"):
            emit({"event": "begin", "operation": "q4+q6"})
            samples, _probe, note = run_samples(both, False, 0)
            _done("q4+q6", "d-SEAMS", version, n, samples,
                  "seams_frame_q4_q6; both averages null", note)

        def neighbours():
            sit.fail_lang(lib, lib.seams_frame_neighbours(
                handle, CUTOFF, 1, b"cutoff"), "neighbours")
            return int(lib.seams_frame_neighbour_pairs(handle))

        if not skip("neighbours"):
            emit({"event": "begin", "operation": "neighbours"})
            samples, pairs, note = run_samples(neighbours, False, 0)
            save_agree(threads, "d-SEAMS", "neighbours", [pairs])
            _done("neighbours", "d-SEAMS", version, n, samples,
                  f"seams_frame_neighbours cutoff pairs={pairs}", note)

        labels = (ctypes.c_int * n)()
        counts9 = (ctypes.c_int * 9)()

        def chill():
            sit.fail_lang(lib, lib.seams_frame_chill_plus(
                handle, CUTOFF, 1, b"cutoff", labels, counts9), "chill+")
            return int(counts9[0])

        if not skip("chill+"):
            emit({"event": "begin", "operation": "chill+"})
            samples, cubic, note = run_samples(chill, False, 0)
            if cubic != n:
                raise RuntimeError(f"CHILL+ cubic {cubic} != {n}")
            save_agree(threads, "d-SEAMS", "chill+", np.ctypeslib.as_array(labels))
            _done("chill+", "d-SEAMS", version, n, samples, f"cubic={cubic}", note)

        core.seams_symmetry_classify.restype = ctypes.c_int
        core.seams_symmetry_classify.argtypes = [
            ctypes.POINTER(ctypes.c_double), ctypes.c_int,
            ctypes.POINTER(ctypes.c_double), ctypes.POINTER(ctypes.c_int)]
        sym = np.empty(n, dtype=np.int32)
        lat = np.zeros(9, dtype=np.float64)
        lat[0] = lat[4] = lat[8] = box[3]
        xyz_ptr = xyz.ctypes.data_as(ctypes.POINTER(ctypes.c_double))
        lat_ptr = lat.ctypes.data_as(ctypes.POINTER(ctypes.c_double))
        lab_ptr = sym.ctypes.data_as(ctypes.POINTER(ctypes.c_int))

        def symmetry():
            code = core.seams_symmetry_classify(xyz_ptr, n, lat_ptr, lab_ptr)
            if code != 0:
                raise RuntimeError(f"seams_symmetry_classify returned {code}")
            return int(np.count_nonzero(sym == 1))

        if not skip("classifier"):
            emit({"event": "begin", "operation": "classifier"})
            samples, ic, note = run_samples(symmetry, False, 0)
            if ic != n:
                raise RuntimeError(f"classifier ice Ic {ic} != {n}")
            save_agree(threads, "d-SEAMS", "classifier", sym)
            _done("classifier", "d-SEAMS", version, n, samples,
                  f"ice-Ic={ic} seams_symmetry_classify", note)
    finally:
        lib.seams_frame_free(handle)


def _done(operation, code, version, molecules, samples, detail, note) -> None:
    if note:
        detail = f"{detail}; {note}" if detail else note
    emit({
        "event": "done",
        "operation": operation,
        "code": code,
        "version": version,
        "molecules": molecules,
        "samples_ms": samples,
        "detail": detail,
        "slow": bool(note),
    })


def _fail(operation, code, reason: str) -> None:
    emit({
        "event": "fail",
        "operation": operation,
        "code": code,
        "reason": reason,
    })


def time_rsm(xyz, side: float, threads: int) -> None:
    version = sit.require_version("rings-statistics-matter", "0.3.0")
    import numpy as np
    from ase import Atoms
    from rings_statistics_matter import ring_statistics
    atoms = Atoms(
        symbols=["O"] * len(xyz),
        positions=np.ascontiguousarray(xyz),
        cell=[side, side, side],
        pbc=True,
    )
    os.environ["RAYON_NUM_THREADS"] = str(threads)

    def run():
        hold_affinity(threads)
        counts, _rings = ring_statistics(
            atoms, cutoff=CUTOFF, maxlvl=8, supercell=False, deduplicate=True,
            backend="rust", num_threads=threads)
        values = np.asarray(counts, dtype=np.float64)
        six = int(round(float(values[5])))
        if six != SIX:
            raise RuntimeError(f"rings-statistics-matter six-rings {six} != {SIX}")
        return six

    if skip("rings"):
        return
    emit({"event": "begin", "operation": "rings"})
    try:
        samples, six, note = run_samples(run, True, threads)
    except TimeoutError:
        _fail("rings", "rings-statistics-matter", "exceeded 120 s")
        return
    except MemoryError:
        _fail("rings", "rings-statistics-matter", "out of memory")
        return
    _done("rings", "rings-statistics-matter", version, MOLECULES, samples,
          f"six={six} maxlvl=8 num_threads={threads} RAYON_NUM_THREADS={threads}", note)


def rival_env(threads: int) -> dict[str, str]:
    """Thread controls for rivals. OpenMP binding is left unset.

    OMP_PROC_BIND=close pins the main thread to one core. Rayon, TBB and
    mdapy inherit that mask, so a requested width of 8 still runs on one
    core. taskset is the cap. Each library then gets its own count.
    """
    env = dict(os.environ)
    width = str(threads)
    env.pop("OMP_PROC_BIND", None)
    env.pop("OMP_PLACES", None)
    env["OMP_NUM_THREADS"] = width
    env["OMP_DYNAMIC"] = "FALSE"
    env["RAYON_NUM_THREADS"] = width
    env["MDAPY_NUM_THREADS"] = width
    env["TBB_NUM_THREADS"] = width
    env["OVITO_THREAD_COUNT"] = width
    env["OPENBLAS_NUM_THREADS"] = "1"
    env["MKL_NUM_THREADS"] = "1"
    env["NUMEXPR_NUM_THREADS"] = "1"
    return env


def arm_mdapy(threads: int) -> str:
    """Thread count for mdapy's kernels, and for Taichi if that runtime is present."""
    os.environ["MDAPY_NUM_THREADS"] = str(threads)
    import mdapy
    from mdapy import get_num_threads, set_num_threads
    set_num_threads(threads)
    note = f"mdapy_threads={get_num_threads()} MDAPY_NUM_THREADS={threads}"
    try:
        import taichi as ti
    except ModuleNotFoundError:
        return note + " taichi=absent"
    ti.init(arch=ti.cpu, cpu_max_num_threads=threads)
    return note + f" taichi_cpu_max_num_threads={threads}"


def time_mdapy(xyz, side: float, threads: int) -> None:
    os.environ["MDAPY_NUM_THREADS"] = str(threads)
    import numpy as np
    import mdapy
    from mdapy import System, get_num_threads, set_num_threads
    from mdapy.analysis import RingStatistics, SteinhardtParameters
    from mdapy.identify import ChillPlus
    from mdapy.neighbor import NeighborList
    version = mdapy.__version__
    thread_note = arm_mdapy(threads)
    set_num_threads(threads)
    system = System(xyz, box=side)
    jobs = []

    def neighbours():
        hold_affinity(threads)
        set_num_threads(threads)
        built = NeighborList.build(system, cutoff=CUTOFF)
        return int(sit.np_len(built))

    jobs.append(("neighbours", neighbours, "pairs"))

    for order, name in ((4, "q4"), (6, "q6")):
        tool = SteinhardtParameters(degrees=(order,), cutoff=CUTOFF)

        def run(tool=tool, order=order):
            hold_affinity(threads)
            set_num_threads(threads)
            result = tool.compute(system)
            if order == 6:
                save_agree(threads, "mdapy", "q6", np.asarray(result.q).reshape(-1))
            return int(result.q.shape[0])

        jobs.append((name, run, "rows"))

    def rings_real():
        import numpy as np
        hold_affinity(threads)
        set_num_threads(threads)
        result = RingStatistics(max_size=6, min_size=3, cutoff=CUTOFF).compute(system)
        values = np.asarray(result.counts)
        return int(values.reshape(-1)[-1])

    jobs.append(("rings", rings_real, "last_count"))

    def chill():
        hold_affinity(threads)
        set_num_threads(threads)
        result = ChillPlus(cutoff=CUTOFF).compute(system)
        return int(getattr(result, "structure", result).shape[0])

    jobs.append(("chill+", chill, "rows"))
    _rival_jobs("mdapy", version, jobs, threads, thread_note + f" get_num_threads={get_num_threads()}")


def time_vesin(xyz, side: float, threads: int) -> None:
    import vesin
    import numpy as np
    version = vesin.__version__
    box = np.diag([side, side, side])
    finder = vesin.NeighborList(
        cutoff=CUTOFF, full_list=True, n_threads=threads, algorithm="cell_list")

    def run():
        hold_affinity(threads)
        left, _right = finder.compute(
            points=xyz, box=box, periodic=True, quantities="ij", copy=True)
        pairs = int(np.asarray(left).shape[0])
        save_agree(threads, "vesin", "neighbours", [pairs])
        return pairs

    _rival_jobs(
        "vesin", version, [("neighbours", run, "pairs")], threads,
        f"n_threads={threads} algorithm=cell_list")


def time_freud(xyz, side: float, threads: int) -> None:
    import freud
    import numpy as np
    version = freud.__version__
    os.environ["TBB_NUM_THREADS"] = str(threads)
    freud.parallel.set_num_threads(threads)
    box = freud.box.Box.cube(side)
    pts = np.ascontiguousarray(xyz - side / 2.0)

    def neighbours():
        hold_affinity(threads)
        freud.parallel.set_num_threads(threads)
        query = freud.locality.AABBQuery(box, pts)
        found = query.query(pts, {"r_max": CUTOFF, "exclude_ii": True}).toNeighborList()
        pairs = int(len(found))
        save_agree(threads, "freud", "neighbours", [pairs])
        return pairs

    jobs = [("neighbours", neighbours, "pairs")]
    for order, name in ((4, "q4"), (6, "q6")):
        def run(order=order):
            hold_affinity(threads)
            freud.parallel.set_num_threads(threads)
            op = freud.order.Steinhardt(l=order, average=False)
            op.compute((box, pts), neighbors={"r_max": CUTOFF, "exclude_ii": True})
            values = np.asarray(op.ql).reshape(-1)
            if order == 6:
                save_agree(threads, "freud", "q6", values)
            return int(values.shape[0])

        jobs.append((name, run, "rows"))
    _rival_jobs(
        "freud", version, jobs, threads,
        f"freud_threads={threads} TBB_NUM_THREADS={threads}")


def time_ovito(xyz, side: float, threads: int) -> None:
    os.environ.setdefault("QT_QPA_PLATFORM", "offscreen")
    import ovito
    import numpy as np
    from ovito.data import DataCollection
    from ovito.modifiers import (
        ChillPlusModifier,
        PolyhedralTemplateMatchingModifier,
        SteinhardtParametersModifier,
    )
    from ovito.pipeline import Pipeline, StaticSource
    version = ".".join(str(part) for part in tuple(ovito.version)[:3])
    data = DataCollection()
    data.create_cell(np.diag([side, side, side]), pbc=(True, True, True))
    parts = data.create_particles(count=len(xyz))
    parts.create_property("Position", data=xyz)

    def chill():
        hold_affinity(threads)
        pipe = Pipeline(source=StaticSource(data=data))
        pipe.modifiers.append(ChillPlusModifier(cutoff=CUTOFF))
        result = pipe.compute()
        save_agree(threads, "OVITO", "chill+", np.asarray(result.particles["Structure Type"]))
        return int(result.particles.count)

    jobs = [("chill+", chill, "rows")]
    for order, name in ((4, "q4"), (6, "q6")):
        def run(order=order):
            hold_affinity(threads)
            pipe = Pipeline(source=StaticSource(data=data))
            mod = SteinhardtParametersModifier(
                cutoff=CUTOFF, max_l=order, compute_averaged=False)
            mod.neighbor_mode = SteinhardtParametersModifier.NeighborMode.Cutoff
            pipe.modifiers.append(mod)
            result = pipe.compute()
            return int(np.asarray(result.particles["Steinhardt"]).shape[0])

        jobs.append((name, run, "rows"))

    def ptm(kind: str):
        def run():
            hold_affinity(threads)
            pipe = Pipeline(source=StaticSource(data=data))
            mod = PolyhedralTemplateMatchingModifier(rmsd_cutoff=0.1)
            if kind == "ice":
                for structure in mod.structures:
                    structure.enabled = structure.id in (6, 7)
            pipe.modifiers.append(mod)
            result = pipe.compute()
            values = np.asarray(result.particles["Structure Type"])
            if kind == "ice":
                save_agree(threads, "OVITO", "ptm-ice", values)
            return int(values.shape[0])
        return run

    jobs.append(("ptm-default", ptm("default"), "rows"))
    jobs.append(("ptm-ice", ptm("ice"), "rows"))
    _rival_jobs("OVITO", version, jobs, threads, f"OVITO_THREAD_COUNT={threads}")


def _rival_jobs(code: str, version: str, jobs, threads: int, thread_note: str) -> None:
    import gc
    for operation, call, kind in jobs:
        if skip(operation):
            continue
        emit({"event": "begin", "operation": operation})
        try:
            samples, value, note = run_samples(call, True, threads)
        except TimeoutError:
            _fail(operation, code, "exceeded 120 s")
            gc.collect()
            continue
        except MemoryError:
            _fail(operation, code, "out of memory")
            gc.collect()
            continue
        detail = f"{kind}={value} {thread_note}"
        if operation.startswith("ptm"):
            detail = f"{detail} rmsd=0.1"
            if operation == "ptm-ice":
                detail += " templates=cubic+hexagonal-diamond"
            else:
                detail += " templates=default"
        _done(operation, code, version, MOLECULES, samples, detail, note)
        gc.collect()


def time_pyscal(xyz, side: float, threads: int) -> None:
    import numpy as np
    import pyscal3
    from ase import Atoms
    version = pyscal3.__version__
    os.environ["PYSCAL_NUM_THREADS"] = str(threads)
    pyscal3.set_num_threads(threads)
    atoms = Atoms(
        symbols=["O"] * len(xyz),
        positions=np.ascontiguousarray(xyz),
        cell=[side, side, side],
        pbc=True,
    )

    def neighbours():
        hold_affinity(threads)
        pyscal3.set_num_threads(threads)
        pyscal3.find_neighbors(atoms, method="cutoff", cutoff=CUTOFF)
        rows = atoms.arrays["pyscal_neighbors"]
        pairs = int(np.asarray(rows).size)
        save_agree(threads, "pyscal3", "neighbours", [pairs])
        return pairs

    def q6():
        hold_affinity(threads)
        pyscal3.set_num_threads(threads)
        pyscal3.find_neighbors(atoms, method="cutoff", cutoff=CUTOFF)
        values = np.asarray(pyscal3.steinhardt_parameter(atoms, 6)[0], dtype=np.float64).reshape(-1)
        save_agree(threads, "pyscal3", "q6", values)
        return int(values.shape[0])

    _rival_jobs(
        "pyscal3", version,
        [("neighbours", neighbours, "pairs"), ("q6", q6, "rows")],
        threads, f"pyscal_threads={threads} PYSCAL_NUM_THREADS={threads}",
    )


def time_plumed(xyz, side: float, threads: int) -> None:
    import math
    import ctypes
    import numpy as np
    os.environ["PLUMED_NUM_THREADS"] = str(threads)
    os.environ["OMP_NUM_THREADS"] = str(threads)
    library = os.environ.get("PLUMED_LIB", "")
    if not library:
        raise RuntimeError("set PLUMED_LIB to libplumed.so")
    kernel = library.replace("libplumed.so", "libplumedKernel.so")
    ctypes.CDLL(kernel, mode=ctypes.RTLD_GLOBAL)
    lib = ctypes.CDLL(library, mode=ctypes.RTLD_GLOBAL)
    lib.plumed_create.restype = ctypes.c_void_p
    lib.plumed_cmd.argtypes = [ctypes.c_void_p, ctypes.c_char_p, ctypes.c_void_p]
    handle = lib.plumed_create()
    keep = []

    def cmd(key, val=None):
        if isinstance(val, str):
            buf = ctypes.create_string_buffer(val.encode())
            keep.append(buf)
            lib.plumed_cmd(handle, key.encode(), ctypes.cast(buf, ctypes.c_void_p))
        elif val is None:
            lib.plumed_cmd(handle, key.encode(), None)
        else:
            lib.plumed_cmd(handle, key.encode(), val)

    n = len(xyz)
    positions = np.ascontiguousarray(xyz / 10.0)
    box = np.zeros((3, 3))
    box[0, 0] = box[1, 1] = box[2, 2] = side / 10.0
    mass = np.ones(n)
    charge = np.zeros(n)
    forces = np.zeros_like(positions)
    virial = np.zeros((3, 3))
    cmd("setMDEngine", "python")
    cmd("setNatoms", ctypes.byref(ctypes.c_int(n)))
    cmd("setTimestep", ctypes.byref(ctypes.c_double(0.001)))
    cmd("setKbT", ctypes.byref(ctypes.c_double(1.0)))
    cmd("setLogFile", "/dev/null")
    cmd("init")
    cmd("readInputLine",
        f"obs: Q6 SPECIES=1-{n} SWITCH={{RATIONAL R_0=0.35 D_MAX=0.35 NN=100 MM=200}}")
    rank = np.zeros(1, dtype=np.intc)
    cmd("getDataRank obs", rank.ctypes.data_as(ctypes.c_void_p))
    shape = np.zeros(int(rank[0]), dtype=np.intc)
    cmd("getDataShape obs", shape.ctypes.data_as(ctypes.c_void_p))
    data = np.zeros(tuple(int(v) for v in shape) or (n,), dtype=np.float64)
    cmd("setMemoryForData obs", data.ctypes.data_as(ctypes.c_void_p))
    step = ctypes.c_int(1)
    scale = math.sqrt(4.0 * math.pi / (2 * 6 + 1))

    def run():
        hold_affinity(threads)
        step.value += 1
        cmd("setStep", ctypes.byref(step))
        cmd("setBox", box.ctypes.data_as(ctypes.c_void_p))
        cmd("setPositions", positions.ctypes.data_as(ctypes.c_void_p))
        cmd("setMasses", mass.ctypes.data_as(ctypes.c_void_p))
        cmd("setCharges", charge.ctypes.data_as(ctypes.c_void_p))
        cmd("setForces", forces.ctypes.data_as(ctypes.c_void_p))
        cmd("setVirial", virial.ctypes.data_as(ctypes.c_void_p))
        cmd("calc")
        values = np.array(data, copy=True).reshape(-1) * scale
        save_agree(threads, "PLUMED", "q6", values)
        return int(values.shape[0])

    _rival_jobs(
        "PLUMED", "2.10.0", [("q6", run, "rows")], threads,
        f"PLUMED_NUM_THREADS={threads}",
    )


def time_ovito_coord(xyz, side: float, threads: int) -> None:
    import numpy as np
    import ovito
    from ovito.data import DataCollection
    from ovito.modifiers import CoordinationAnalysisModifier
    version = ".".join(str(part) for part in tuple(ovito.version)[:3])

    def run():
        hold_affinity(threads)
        fresh = DataCollection()
        fresh.create_cell(np.diag([side, side, side]), pbc=(True, True, True))
        atoms = fresh.create_particles(count=len(xyz))
        atoms.create_property("Position", data=xyz)
        fresh.apply(CoordinationAnalysisModifier(cutoff=CUTOFF))
        coord = np.asarray(fresh.particles["Coordination"], dtype=np.int32)
        save_agree(threads, "OVITO", "coordination", coord)
        return int(np.sum(coord))

    _rival_jobs(
        "OVITO", version, [("coordination", run, "sum")], threads,
        f"OVITO_THREAD_COUNT={threads} CoordinationAnalysisModifier cutoff={CUTOFF}",
    )


def worker(args) -> int:
    try:
        os.sched_setaffinity(0, set(range(args.threads)))
    except (AttributeError, OSError):
        pass
    os.environ["OVITO_THREAD_COUNT"] = str(args.threads)
    os.environ["QT_QPA_PLATFORM"] = "offscreen"
    cores, affinity = stats.cores_now()
    emit({"event": "hello", "cores": cores, "affinity": affinity})
    if cores != args.threads:
        raise RuntimeError(f"affinity has {cores} cores, width is {args.threads}")
    if args.from_op:
        os.environ["SITTING_FROM"] = args.from_op
    xyz, box, side = ice()
    key = args.code
    if key == "dseams":
        import ctypes
        lib = sit.open_lang(args.lib)
        core = ctypes.CDLL(args.core)
        time_dseams(lib, core, xyz, box, args.version, args.threads)
        return 0
    if key == "rsm":
        time_rsm(xyz, side, args.threads)
        return 0
    if key == "mdapy":
        time_mdapy(xyz, side, args.threads)
        return 0
    if key == "vesin":
        time_vesin(xyz, side, args.threads)
        return 0
    if key == "freud":
        time_freud(xyz, side, args.threads)
        return 0
    if key == "ovito":
        time_ovito(xyz, side, args.threads)
        return 0
    if key == "pyscal":
        time_pyscal(xyz, side, args.threads)
        return 0
    if key == "plumed":
        time_plumed(xyz, side, args.threads)
        return 0
    if key == "ovito-coord":
        time_ovito_coord(xyz, side, args.threads)
        return 0
    raise RuntimeError(f"unknown code {key}")


def drive(code_key: str, threads: int, rival: bool, args, env) -> list[dict]:
    """Run one code at one width. Restart after a killed call."""
    rows = []
    start = ""
    attempts = 0
    while True:
        attempts += 1
        if attempts > 12:
            raise RuntimeError(f"{code_key} at {threads} threads restarted too many times")
        cmd = [
            "taskset", "-c", stats.taskset_cores(threads),
            sys.executable, str(Path(__file__).resolve()),
            "--worker", "--code", code_key, "--threads", str(threads),
            "--lib", args.lib, "--core", args.core, "--version", args.version,
        ]
        if start:
            cmd.extend(["--from-op", start])
        proc = subprocess.Popen(
            cmd, cwd=ROOT, env=env, text=True, bufsize=1,
            stdout=subprocess.PIPE, stderr=subprocess.PIPE)
        state = {"sample": False, "op": "", "killed": ""}
        stderr_parts: list[str] = []

        def read_err():
            assert proc.stderr is not None
            stderr_parts.append(proc.stderr.read())

        err_thread = threading.Thread(target=read_err, daemon=True)
        err_thread.start()

        def watch():
            while proc.poll() is None:
                time.sleep(0.25)
                if rival and state["sample"] and state["killed"] == "":
                    # sample-begin stores the deadline on the state
                    deadline = state.get("deadline", 0.0)
                    if deadline and time.monotonic() > deadline:
                        state["killed"] = "timeout"
                        proc.kill()
                        return

        threading.Thread(target=watch, daemon=True).start()
        assert proc.stdout is not None
        finished = False
        failed = ""
        hello = False
        try:
            for line in proc.stdout:
                line = line.strip()
                if not line.startswith("{"):
                    continue
                msg = json.loads(line)
                event = msg.get("event")
                if event == "hello":
                    hello = True
                    if int(msg["cores"]) != threads:
                        proc.kill()
                        raise RuntimeError(
                            f"{code_key} pinned {msg['cores']} cores for width {threads}")
                    continue
                if event == "begin":
                    state["op"] = msg["operation"]
                    continue
                if event == "sample-begin":
                    state["sample"] = True
                    state["deadline"] = time.monotonic() + LIMIT_S
                    continue
                if event == "sample-end":
                    state["sample"] = False
                    state["deadline"] = 0.0
                    continue
                if event == "done":
                    rows.append(msg)
                    state["op"] = ""
                    continue
                if event == "fail":
                    rows.append(msg)
                    failed = msg["operation"]
                    state["op"] = ""
                    continue
        finally:
            proc.wait(timeout=30)
            err_thread.join(timeout=5)
        if proc.returncode == 0 and hello:
            return rows
        op = state["op"] or failed
        if state["killed"] == "timeout":
            rows.append({
                "event": "fail", "operation": op, "code": _name(code_key),
                "reason": "exceeded 120 s",
            })
            start = _next_op(code_key, op)
            if not start:
                return rows
            continue
        err = "".join(stderr_parts)
        if proc.returncode and proc.returncode < 0 and not hello:
            raise RuntimeError(f"{code_key} failed before measuring\n{err[-2000:]}")
        out_of_memory = (proc.returncode is not None and proc.returncode < 0) or "MemoryError" in err
        if out_of_memory:
            rows.append({
                "event": "fail", "operation": op, "code": _name(code_key),
                "reason": "out of memory",
            })
            start = _next_op(code_key, op)
            if not start:
                return rows
            continue
        raise RuntimeError(f"{code_key} threads {threads} exit {proc.returncode}\n{err[-4000:]}")


def _name(code_key: str) -> str:
    names = {key: name for key, name, _rival in CODES}
    names.update({"pyscal": "pyscal3", "plumed": "PLUMED", "ovito-coord": "OVITO"})
    return names.get(code_key, code_key)


OPS = {
    "dseams": ["rings", "q4", "q6", "q4+q6", "neighbours", "chill+", "classifier"],
    "rsm": ["rings"],
    "mdapy": ["neighbours", "q4", "q6", "rings", "chill+"],
    "vesin": ["neighbours"],
    "freud": ["neighbours", "q4", "q6"],
    "ovito": ["chill+", "q4", "q6", "ptm-default", "ptm-ice"],
    "pyscal": ["neighbours", "q6"],
    "plumed": ["q6"],
    "ovito-coord": ["coordination"],
}


def _next_op(code_key: str, operation: str) -> str:
    ops = OPS[code_key]
    if operation not in ops:
        return ""
    index = ops.index(operation) + 1
    if index >= len(ops):
        return ""
    return ops[index]


def row_from(msg: dict, threads: int, affinity: str, info: dict) -> dict:
    if msg.get("event") == "fail" or "samples_ms" not in msg:
        return {
            "source_commit": info["source_commit"],
            "cpu_model": info["cpu_model"],
            "cpu_cores": info["cpu_cores"],
            "affinity": affinity,
            "phase": "ice",
            "operation": msg.get("operation", ""),
            "code": msg.get("code", ""),
            "version": msg.get("version", ""),
            "molecules": MOLECULES,
            "threads": threads,
            "warmup": WARMUP,
            "n": 0,
            "median_ms": msg.get("reason", "out of memory"),
            "iqr_ms": "",
            "q1_ms": "",
            "q3_ms": "",
            "samples_ms": "",
            "detail": msg.get("reason", ""),
        }
    summary = stats.summarize(msg["samples_ms"])
    return {
        "source_commit": info["source_commit"],
        "cpu_model": info["cpu_model"],
        "cpu_cores": info["cpu_cores"],
        "affinity": affinity,
        "phase": "ice",
        "operation": msg["operation"],
        "code": msg["code"],
        "version": msg["version"],
        "molecules": msg["molecules"],
        "threads": threads,
        "warmup": WARMUP,
        "n": summary["n"],
        "median_ms": f"{summary['median_ms']:.6f}",
        "iqr_ms": f"{summary['iqr_ms']:.6f}",
        "q1_ms": f"{summary['q1_ms']:.6f}",
        "q3_ms": f"{summary['q3_ms']:.6f}",
        "samples_ms": stats.fmt_samples(msg["samples_ms"]),
        "detail": msg.get("detail", ""),
    }


def write_outputs(out: Path, rows: list[dict], info: dict) -> None:
    out.mkdir(parents=True, exist_ok=True)
    dest = out / "sitting-large.csv"
    with dest.open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=FIELDS)
        writer.writeheader()
        writer.writerows(rows)
    log = [
        *stats.machine_lines(ROOT),
        f"cores {info['cores']}",
        "lscpu",
        *lscpu_details(),
        "warmup 1",
        "n 8",
        "n_if_rival_call_exceeds_30s 6",
        "rival_call_limit_s 120",
        f"molecules {MOLECULES}",
        "cell ice Ic 48^3 conventional cells, eight times 24^3",
        "cutoff_A 3.5",
        "widths 1 2 4 8",
        "one sitting; each width and each code is one process",
        "host name is not recorded",
    ]
    for row in rows:
        log.append(
            f"{row['operation']} {row['code']} threads {row['threads']} "
            f"n {row['n']} median_ms {row['median_ms']} iqr_ms {row['iqr_ms']} "
            f"affinity {row['affinity']} {row['detail']}"
        )
    (out / "raw-sitting-large.txt").write_text("\n".join(log) + "\n", encoding="utf-8")
    print(f"wrote {dest} ({len(rows)} rows)", flush=True)


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--out", type=Path, default=ROOT / "benchmarks/paper-cells")
    parser.add_argument("--lib", default=str(ROOT / "bbdir-release/bindings/libseams_lang.so"))
    parser.add_argument("--core", default=str(ROOT / "bbdir-release/src/libyodaLib.so"))
    parser.add_argument("--seams", default=str(ROOT / "bbdir-release/seams"))
    parser.add_argument("--worker", action="store_true")
    parser.add_argument("--code", default="dseams")
    parser.add_argument("--threads", type=int, default=8)
    parser.add_argument("--version", default="")
    parser.add_argument("--from-op", default="")
    parser.add_argument("--only", default="", help="comma-separated code keys")
    parser.add_argument("--widths", default="", help="comma-separated thread counts")
    parser.add_argument("--extra", action="store_true",
                        help="time pyscal3, PLUMED q6 and OVITO coordination and merge them in")
    args = parser.parse_args()
    if args.worker:
        return worker(args)
    if not Path(args.lib).is_file() or not Path(args.core).is_file():
        raise SystemExit("missing the release libraries")
    info = stats.machine(ROOT)
    args.version = args.version or sit.seams_version(Path(args.seams))
    if args.version == "unrecorded":
        args.version = project_version()
    only = set(args.only.split(",")) if args.only else set()
    widths = tuple(int(part) for part in args.widths.split(",") if part) or WIDTHS
    catalog = list(CODES)
    if args.extra:
        catalog = [
            ("dseams", "d-SEAMS", False),
            ("pyscal", "pyscal3", True),
            ("plumed", "PLUMED", True),
            ("ovito-coord", "OVITO", True),
        ]
    if args.extra:
        agree = args.out / ".sitting-agree"
        if agree.is_dir():
            for path in agree.glob("*.npy"):
                path.unlink()
    collected = []
    for threads in widths:
        affinity = stats.taskset_cores(threads)
        for key, _name, rival in catalog:
            if only and key not in only:
                continue
            print(f"threads {threads} code {key}", flush=True)
            env = rival_env(threads) if rival else stats.pin_env(threads)
            if args.extra and key == "dseams":
                env["SITTING_ONLY"] = "rings,q6,neighbours"
            prior = env.get("LD_LIBRARY_PATH", "")
            libs = sit.library_path(args.lib)
            system = "/usr/lib/x86_64-linux-gnu"
            parts = [libs]
            if prior:
                parts.append(prior)
            parts.append(system)
            env["LD_LIBRARY_PATH"] = os.pathsep.join(parts)
            env["QT_QPA_PLATFORM"] = "offscreen"
            env["PYTHONUNBUFFERED"] = "1"
            env["SIT_AGREE"] = str(args.out / ".sitting-agree")
            messages = drive(key, threads, rival, args, env)
            for msg in messages:
                if not msg.get("code"):
                    msg["code"] = _name
                collected.append(row_from(msg, threads, affinity, info))
            print(
                f"threads {threads} code {key} rows "
                f"{sum(1 for row in collected if row['threads'] == threads and row['code'] == _name)}",
                flush=True)
            if not args.extra:
                write_outputs(args.out, collected, info)
    if not collected:
        raise SystemExit("the sitting produced no rows")
    problems = record_agreement(collected, args.out / ".sitting-agree")
    if args.extra:
        collected = merge_existing(args.out / "sitting-large.csv", collected)
    write_outputs(args.out, collected, info)
    require_rival_scaling([row for row in collected if row["code"] in {
        "pyscal3", "PLUMED", "OVITO", "rings-statistics-matter", "mdapy", "vesin", "freud",
    }])
    if problems:
        raise SystemExit("results do not match the rivals\n" + "\n".join(problems))
    return 0


def merge_existing(path: Path, new_rows: list[dict]) -> list[dict]:
    if not path.is_file():
        return new_rows
    old = list(csv.DictReader(path.open()))
    index = {
        (row["code"], row["operation"], str(row["threads"])): i
        for i, row in enumerate(old)
    }
    for row in new_rows:
        key = (row["code"], row["operation"], str(row["threads"]))
        if key in index:
            old[index[key]] = row
        else:
            old.append(row)
    return old


def _note(row: dict, text: str) -> None:
    detail = row.get("detail", "")
    row["detail"] = f"{detail} {text}".strip()


def record_agreement(rows: list[dict], root: Path) -> None:
    """Stamp each paired row with how it compares to the other codes.

    q6 values, neighbour-pair counts, CHILL+ cubic labels and the
    classifier versus ice-only PTM are compared on the same frame.
    """
    import numpy as np
    by = {(str(r["threads"]), r["code"], r["operation"]): r for r in rows}
    problems = []

    def load(threads: str, code: str, operation: str):
        path = root / f"t{threads}-{code}-{operation}.npy"
        if not path.is_file():
            return None
        return np.load(path)

    for threads in ("1", "2", "4", "8"):
        ours = load(threads, "d-SEAMS", "q6")
        freud = load(threads, "freud", "q6")
        mdapy = load(threads, "mdapy", "q6")
        gaps = []
        if ours is not None and freud is not None:
            gaps.append(("freud", float(np.max(np.abs(ours - freud)))))
        if ours is not None and mdapy is not None:
            gaps.append(("mdapy", float(np.max(np.abs(ours - mdapy)))))
        pyscal = load(threads, "pyscal3", "q6")
        plumed = load(threads, "PLUMED", "q6")
        if ours is not None and pyscal is not None:
            gaps.append(("pyscal3", float(np.max(np.abs(ours - pyscal)))))
        if ours is not None and plumed is not None:
            gaps.append(("PLUMED", float(np.max(np.abs(ours - plumed)))))
        if gaps:
            text = " ".join(f"agree_{name}_max_abs={gap:.3e}" for name, gap in gaps)
            codes = ["d-SEAMS"] + [name for name, _gap in gaps]
            for code in codes:
                row = by.get((threads, code, "q6"))
                if row:
                    _note(row, text)
            for name, gap in gaps:
                if gap > 1e-4:
                    problems.append(f"q6 threads {threads} {name} {gap:.3e}")
        counts = []
        names = ("d-SEAMS", "vesin", "freud", "pyscal3")
        present = []
        for code in names:
            values = load(threads, code, "neighbours")
            count = None if values is None else int(np.asarray(values).reshape(-1)[0])
            if count is not None:
                present.append(code)
                counts.append(count)
        if "d-SEAMS" in present and len(present) >= 2:
            text = "pairs=" + "/".join(
                f"{code}:{count}" for code, count in zip(present, counts)
            )
            matched = len(set(counts)) == 1
            text += " agree=" + ",".join(present) if matched else " disagree=" + ",".join(present)
            if not matched:
                problems.append(f"neighbours threads {threads} {list(zip(present, counts))}")
            for code in present:
                row = by.get((threads, code, "neighbours"))
                if row:
                    _note(row, text)
        dchill = load(threads, "d-SEAMS", "chill+")
        ochill = load(threads, "OVITO", "chill+")
        if dchill is not None and ochill is not None:
            d_cubic = int(np.count_nonzero(dchill == 0))
            o_cubic = int(np.count_nonzero(ochill == 2))
            text = f"cubic_dseams={d_cubic} cubic_ovito={o_cubic}"
            for code in ("d-SEAMS", "OVITO"):
                row = by.get((threads, code, "chill+"))
                if row:
                    _note(row, text)
            if d_cubic != len(dchill) or o_cubic != len(ochill):
                problems.append(f"chill+ threads {threads} d-SEAMS {d_cubic} OVITO {o_cubic}")
        labels = load(threads, "d-SEAMS", "classifier")
        ptm = load(threads, "OVITO", "ptm-ice")
        if labels is not None and ptm is not None:
            matched = int(np.count_nonzero((labels == 1) & (ptm == 6)))
            text = f"agree_ice_ic_cubic_diamond={matched}/{len(labels)}"
            for code, operation in (("d-SEAMS", "classifier"), ("OVITO", "ptm-ice")):
                row = by.get((threads, code, operation))
                if row:
                    _note(row, text)
            if matched != len(labels):
                problems.append(f"classifier threads {threads} {matched}/{len(labels)}")
        coord = load(threads, "OVITO", "coordination")
        ours_pairs = load(threads, "d-SEAMS", "neighbours")
        if coord is not None and ours_pairs is not None:
            total = int(np.sum(coord))
            pairs = int(np.asarray(ours_pairs).reshape(-1)[0])
            text = f"coordination_sum={total} dseams_pairs={pairs}"
            for code, operation in (("OVITO", "coordination"), ("d-SEAMS", "neighbours")):
                row = by.get((threads, code, operation))
                if row:
                    _note(row, text)
            if total != pairs:
                problems.append(f"coordination threads {threads} {total} vs {pairs}")
        for code in ("d-SEAMS", "rings-statistics-matter"):
            row = by.get((threads, code, "rings"))
            if row and "six=1769472" in row.get("detail", ""):
                _note(row, "agree_six=1769472")
    return problems


def require_rival_scaling(rows: list[dict]) -> None:
    """A rival row is trusted only when its time falls from 1 thread to 8."""
    rivals = {name for _key, name, rival in CODES if rival}
    grouped: dict[tuple[str, str], dict[int, float]] = {}
    for row in rows:
        if row["code"] not in rivals:
            continue
        try:
            value = float(row["median_ms"])
        except ValueError:
            raise SystemExit(
                f"{row['code']} {row['operation']} at {row['threads']} threads "
                f"is {row['median_ms']}, not a time")
        grouped.setdefault((row["code"], row["operation"]), {})[int(row["threads"])] = value
    problems = []
    for key in sorted(grouped):
        times = grouped[key]
        if 1 not in times or 8 not in times:
            problems.append(f"{key[0]} {key[1]} missing a width: {sorted(times)}")
            continue
        if not times[8] < times[1]:
            text = " ".join(f"{width}:{times[width]:.1f}" for width in sorted(times))
            problems.append(f"{key[0]} {key[1]} does not fall with threads: {text}")
    if problems:
        raise SystemExit("rival thread counts did not take effect\n" + "\n".join(problems))


if __name__ == "__main__":
    sys.exit(main())
