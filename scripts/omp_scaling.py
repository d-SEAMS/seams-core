#!/usr/bin/env python3
# Copyright (c) d-SEAMS developers
"""OpenMP scaling of the d-SEAMS kernels, idle and beside a busy neighbour.

Ideal ice Ic at 110592 and 884736 molecules, cutoff 3.5 A. Each d-SEAMS
call does its whole job, neighbour list included, at 1, 2, 4 and 8
threads:

  neighbours, q4, q6, rings   flat C API in libyodaLib (seams_cutoff_csr,
                              seams_ql_cutoff, seams_primitive_rings)
  chill, chill_ovito,         libseams frame handle (seams_chill,
  chill_plus, chill_plus_ovito  seams_chill_plus, published and OVITO)
  symmetry                    seams_system_symmetry, the shell match
  knn                         seams_neighbours_knn, four nearest, one-way
  rings_find                  seams_rings_find, rings with symmetry codes
  cages                       seams_cages_find
  prisms                      seams_prisms
  monolayer                   seams_monolayer, the ring census
  rdf2d                       seams_frame_rdf2d in libseams_lang, 8 A

mdapy and pyscal3 q4 and q6, OVITO q6, and the four-nearest graph from
mdapy, pyscal3, freud and OVITO, run at the rival widths (1, 2, 4 and 8
by default) on both frames in the same sitting. Every rival call builds
its own neighbour list, as the d-SEAMS calls do.

Two placements. A bound row runs under taskset on cpus 0 .. width-1 with
OMP_PROC_BIND=close and OMP_PLACES=cores; a rival's bound row keeps the
taskset but not the OpenMP binding, which would pin its main thread and
every thread its own pool starts later. An unbound row is the default
environment: the whole process mask, OMP_NUM_THREADS (and the rival's
own thread count) and nothing else, with OMP_PROC_BIND, OMP_PLACES,
OMP_WAIT_POLICY and GOMP_SPINCOUNT removed. d-SEAMS runs the unbound
placement for --unbound-ops.

The contended rows repeat the 8-thread rows while scripts/omp_spinner.c
holds an 8-thread OpenMP team in its default wait (OMP_WAIT_POLICY unset)
on the same cpus. The build labelled "<label>-fixed" is the after build
with SEAMS_FAIR_WIDTH=0, which keeps the requested width under contention.

One process per code, frame, operation, width, placement and load. One
warm-up call, then n timed calls (n >= 10). The median of an even sample
is the mean of the two central values and the IQR is type 7. Each row
records the core count and the affinity of the worker. The header
records the tree of the measuring checkout and no commit hash. No host
name is read or written.

--before times a second libyodaLib.so, for example a build of the parent
commit, under the label given by --before-label. libseams.so is read from
the same directory and libseams_lang.so from ../bindings. The library
paths are not written to the outputs.
"""

from __future__ import annotations

import argparse
import csv
import json
import os
import shutil
import subprocess
import sys
import tempfile
import time
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "scripts"))

import paper_stats as stats  # noqa: E402

CUTOFF = 3.5
RDF_CUTOFF = 8.0
RDF_BIN = 0.05
REPETITIONS = (24, 48)
MAX_RING = 6
WARMUP = 1
MIN_REPEATS = 10
WIDTHS = (1, 2, 4, 8)
FLAT_OPS = ("neighbours", "q4", "q6", "rings")
SYSTEM_OPS = ("chill", "chill_ovito", "chill_plus", "chill_plus_ovito", "symmetry", "knn",
              "rings_find", "cages", "prisms", "monolayer")
LANG_OPS = ("rdf2d",)
DSEAMS_OPS = FLAT_OPS + SYSTEM_OPS + LANG_OPS
KNN = 4
RIVAL_OPS = {"mdapy": ("q4", "q6", "knn"), "pyscal3": ("q4", "q6", "knn"), "freud": ("knn",),
             "ovito": ("q6", "knn")}
PLACEMENTS = ("bound", "unbound")
UNBOUND_OPS = ("neighbours", "q6", "knn")
UNBOUND_DROP = ("OMP_PROC_BIND", "OMP_PLACES", "OMP_WAIT_POLICY", "GOMP_SPINCOUNT",
                "GOMP_CPU_AFFINITY", "KMP_AFFINITY", "OMP_DISPLAY_ENV", "OMP_SCHEDULE",
                "SEAMS_PIN")
SPINNER = ROOT / "scripts" / "omp_spinner.c"
TRACKED = ("src", "include", "bindings", "scripts")


def frame(repetitions: int):
    import h2h_protocol as proto
    return proto.ice_ic(repetitions)


START_MASK = set(os.sched_getaffinity(0)) if hasattr(os, "sched_getaffinity") else set()


def hold_affinity(threads: int) -> int:
    """Put the process back on the mask it started with after a runtime narrows it."""
    try:
        os.sched_setaffinity(0, START_MASK)
        return len(os.sched_getaffinity(0))
    except (AttributeError, OSError):
        return threads


def sibling_libs(lib_path: str) -> tuple[Path, Path]:
    lib = Path(lib_path).resolve()
    return lib.parent / "libseams.so", lib.parent.parent / "bindings" / "libseams_lang.so"


def flat_call(lib_path: str, op: str, xyz, side: float):
    import ctypes
    import numpy as np
    n = len(xyz)
    lib = ctypes.CDLL(lib_path)
    dp = ctypes.POINTER(ctypes.c_double)
    ip = ctypes.POINTER(ctypes.c_int)
    box = np.array([0.0, 0.0, 0.0, side, side, side, 0.0, 0.0, 0.0])
    xp = xyz.ctypes.data_as(dp)
    bp = box.ctypes.data_as(dp)
    if op == "neighbours":
        lib.seams_cutoff_csr.argtypes = [dp, ctypes.c_int, dp, ctypes.c_double, ip, ip,
                                         ctypes.c_int, ip]
        offsets = np.empty(n + 1, dtype=np.int32)
        cols = np.empty(16 * n, dtype=np.int32)
        count = ctypes.c_int(0)

        def run():
            rc = lib.seams_cutoff_csr(xp, n, bp, CUTOFF, offsets.ctypes.data_as(ip),
                                      cols.ctypes.data_as(ip), len(cols), ctypes.byref(count))
            if rc != 0:
                raise RuntimeError(f"seams_cutoff_csr returned {rc}")
            return f"pairs={count.value}"
        return run
    if op in ("q4", "q6"):
        lib.seams_ql_cutoff.argtypes = [dp, ctypes.c_int, dp, ctypes.c_double, ctypes.c_int, dp]
        ql = np.empty(n, dtype=np.float64)
        order = int(op[1])

        def run():
            rc = lib.seams_ql_cutoff(xp, n, bp, CUTOFF, order, ql.ctypes.data_as(dp))
            if rc != 0:
                raise RuntimeError(f"seams_ql_cutoff returned {rc}")
            return f"mean={float(ql.mean()):.12f}"
        return run
    if op == "rings":
        lib.seams_primitive_rings.argtypes = [dp, ctypes.c_int, dp, ctypes.c_double,
                                              ctypes.c_int, ip, ctypes.c_int]
        counts = np.zeros(MAX_RING + 1, dtype=np.int32)

        def run():
            rc = lib.seams_primitive_rings(xp, n, bp, CUTOFF, MAX_RING,
                                           counts.ctypes.data_as(ip), len(counts))
            if rc != 0:
                raise RuntimeError(f"seams_primitive_rings returned {rc}")
            return f"six={int(counts[6])}"
        return run
    raise ValueError(op)


def system_call(lib_path: str, op: str, xyz, side: float):
    import ctypes
    import numpy as np
    n = len(xyz)
    api_path, _ = sibling_libs(lib_path)
    lib = ctypes.CDLL(str(api_path))
    dp = ctypes.POINTER(ctypes.c_double)
    ip = ctypes.POINTER(ctypes.c_int)
    vp = ctypes.c_void_p
    size = ctypes.c_size_t
    num = ctypes.c_double
    lib.seams_system_create.restype = vp
    lib.seams_system_set_ortho.argtypes = [vp, dp, ctypes.c_int] + [num] * 6 + [ip, ip]
    handle = lib.seams_system_create()
    rc = lib.seams_system_set_ortho(handle, xyz.ctypes.data_as(dp), n, 0.0, 0.0, 0.0,
                                    side, side, side, None, None)
    if rc != 0:
        raise RuntimeError(f"seams_system_set_ortho returned {rc}")
    labels = np.zeros(n, dtype=np.int32)
    counts9 = np.zeros(9, dtype=np.int32)
    lp = labels.ctypes.data_as(ip)
    cp = counts9.ctypes.data_as(ip)

    def ok(rc, name):
        if rc != 0:
            raise RuntimeError(f"{name} returned {rc}")

    if op in ("chill", "chill_ovito"):
        lib.seams_chill.argtypes = [vp, num, ctypes.c_int, ctypes.c_int, ctypes.c_int, ip,
                                    size, ip]
        rule = 1 if op.endswith("ovito") else 0

        def run():
            ok(lib.seams_chill(handle, CUTOFF, 0, 0, rule, lp, n, cp), "seams_chill")
            return "counts=" + ":".join(str(int(c)) for c in counts9)
        return run
    if op in ("chill_plus", "chill_plus_ovito"):
        lib.seams_chill_plus.argtypes = [vp, num, ctypes.c_int, ctypes.c_int, ctypes.c_int,
                                         ip, ip]
        rule = 1 if op.endswith("ovito") else 0

        def run():
            ok(lib.seams_chill_plus(handle, CUTOFF, 0, 0, rule, lp, cp), "seams_chill_plus")
            return "counts=" + ":".join(str(int(c)) for c in counts9)
        return run
    if op == "symmetry":
        lib.seams_system_symmetry.argtypes = [vp, num, ip, size]

        def run():
            ok(lib.seams_system_symmetry(handle, 0.0, lp, n), "seams_system_symmetry")
            tally = np.bincount(labels + 1, minlength=4)
            return "labels=" + ":".join(str(int(c)) for c in tally)
        return run
    if op == "knn":
        lib.seams_neighbours_knn.argtypes = [vp, ctypes.c_int, num, ctypes.c_int, ctypes.c_int,
                                             ctypes.POINTER(vp)]
        lib.seams_neighbours_copy.argtypes = [vp, ip, size, ctypes.POINTER(size), ip, size,
                                              ctypes.POINTER(size)]
        lib.seams_neighbours_destroy.argtypes = [vp]

        def run():
            graph = vp()
            ok(lib.seams_neighbours_knn(handle, KNN, CUTOFF, 0, 0, ctypes.byref(graph)),
               "seams_neighbours_knn")
            rows = size()
            pairs = size()
            lib.seams_neighbours_copy(graph, None, 0, ctypes.byref(rows), None, 0,
                                      ctypes.byref(pairs))
            lib.seams_neighbours_destroy(graph)
            return f"pairs={pairs.value}"
        return run
    if op == "rings_find":
        lib.seams_rings_find.argtypes = [vp, num, ctypes.c_int, ctypes.c_int, ctypes.c_int,
                                         ctypes.POINTER(vp)]
        lib.seams_rings_histogram.argtypes = [vp, ip, ctypes.c_int]
        lib.seams_rings_destroy.argtypes = [vp]
        histogram = np.zeros(MAX_RING + 1, dtype=np.int32)

        def run():
            rings = vp()
            ok(lib.seams_rings_find(handle, CUTOFF, 0, 0, MAX_RING, ctypes.byref(rings)),
               "seams_rings_find")
            lib.seams_rings_histogram(rings, histogram.ctypes.data_as(ip), len(histogram))
            lib.seams_rings_destroy(rings)
            return f"six={int(histogram[6])}"
        return run
    if op == "cages":
        lib.seams_cages_find.argtypes = [vp, num, ctypes.c_int, ctypes.c_int, ctypes.c_int,
                                         ctypes.c_int, ctypes.POINTER(vp)]
        lib.seams_cages_counts.argtypes = [vp, ip, ip, ip, ip]
        lib.seams_cages_destroy.argtypes = [vp]

        def run():
            cages = vp()
            ok(lib.seams_cages_find(handle, CUTOFF, 0, 0, 4, 0, ctypes.byref(cages)),
               "seams_cages_find")
            found = [ctypes.c_int() for _ in range(4)]
            lib.seams_cages_counts(cages, *(ctypes.byref(v) for v in found))
            lib.seams_cages_destroy(cages)
            return "cages=" + ":".join(str(v.value) for v in found)
        return run
    if op == "prisms":
        lib.seams_prisms.argtypes = [vp, num, ctypes.c_int, ip, ip, size]
        perfect = np.zeros(7, dtype=np.int32)
        imperfect = np.zeros(7, dtype=np.int32)

        def run():
            ok(lib.seams_prisms(handle, CUTOFF, 0, perfect.ctypes.data_as(ip),
                                imperfect.ctypes.data_as(ip), 7), "seams_prisms")
            return f"prisms={int(perfect.sum())}:{int(imperfect.sum())}"
        return run
    if op == "monolayer":
        lib.seams_monolayer.argtypes = [vp, num, ctypes.c_int, ctypes.c_int, num, ip, size]
        census = np.zeros(MAX_RING + 1, dtype=np.int32)

        def run():
            ok(lib.seams_monolayer(handle, CUTOFF, 0, MAX_RING, side * side,
                                   census.ctypes.data_as(ip), len(census)), "seams_monolayer")
            return f"six={int(census[6])}"
        return run
    raise ValueError(op)


def lang_call(lib_path: str, op: str, xyz, side: float):
    import ctypes
    import numpy as np
    if op != "rdf2d":
        raise ValueError(op)
    _, lang_path = sibling_libs(lib_path)
    lang = ctypes.CDLL(str(lang_path))
    dp = ctypes.POINTER(ctypes.c_double)
    ip = ctypes.POINTER(ctypes.c_int)
    vp = ctypes.c_void_p
    lang.seams_frame_from_arrays.restype = vp
    lang.seams_frame_from_arrays.argtypes = [dp, ctypes.c_int, dp, ip, ip]
    lang.seams_frame_rdf2d.argtypes = [vp, ctypes.c_double, ctypes.c_double, ctypes.c_char_p,
                                       dp, dp]
    box = np.array([0.0, 0.0, 0.0, side, side, side, 0.0, 0.0, 0.0])
    handle = lang.seams_frame_from_arrays(xyz.ctypes.data_as(dp), len(xyz),
                                          box.ctypes.data_as(dp), None, None)
    if not handle:
        raise RuntimeError("seams_frame_from_arrays failed")
    nbin = int(RDF_CUTOFF / RDF_BIN)
    r = np.zeros(nbin)
    g = np.zeros(nbin)

    def run():
        rc = lang.seams_frame_rdf2d(handle, RDF_CUTOFF, RDF_BIN, None, r.ctypes.data_as(dp),
                                    g.ctypes.data_as(dp))
        if rc != 0:
            raise RuntimeError(f"seams_frame_rdf2d returned {rc}")
        return f"gsum={float(g.sum()):.9f}"
    return run


def dseams_call(lib_path: str, op: str, repetitions: int = REPETITIONS[0]):
    xyz, side = frame(repetitions)
    if op in FLAT_OPS:
        return flat_call(lib_path, op, xyz, side)
    if op in SYSTEM_OPS:
        return system_call(lib_path, op, xyz, side)
    if op in LANG_OPS:
        return lang_call(lib_path, op, xyz, side)
    raise ValueError(op)


def mdapy_call(op: str, threads: int, repetitions: int):
    os.environ["MDAPY_NUM_THREADS"] = str(threads)
    import mdapy
    from mdapy import System, set_num_threads
    from mdapy.analysis import SteinhardtParameters
    xyz, side = frame(repetitions)
    set_num_threads(threads)
    system = System(xyz, box=side)
    if op == "knn":
        def run():
            hold_affinity(threads)
            set_num_threads(threads)
            pairs = mdapy.NeighborList.build(system, num_neighbors=KNN)
            return f"pairs={len(pairs.j)}"
        return run, f"mdapy {mdapy.__version__}"
    tool = SteinhardtParameters(degrees=(int(op[1]),), cutoff=CUTOFF)

    def run():
        hold_affinity(threads)
        set_num_threads(threads)
        result = tool.compute(system)
        return f"rows={int(result.q.shape[0])}"
    return run, f"mdapy {mdapy.__version__}"


def pyscal_call(op: str, threads: int, repetitions: int):
    import numpy as np
    import pyscal3
    from ase import Atoms
    os.environ["PYSCAL_NUM_THREADS"] = str(threads)
    pyscal3.set_num_threads(threads)
    xyz, side = frame(repetitions)
    atoms = Atoms(symbols=["O"] * len(xyz), positions=np.ascontiguousarray(xyz),
                  cell=[side, side, side], pbc=True)
    if op == "knn":
        def run():
            hold_affinity(threads)
            pyscal3.set_num_threads(threads)
            pyscal3.find_neighbors(atoms, method="number", nmax=KNN, store_rows=False)
            return f"pairs={len(atoms.info['pyscal_bond_neighbors'])}"
        return run, f"pyscal3 {pyscal3.__version__}"
    order = int(op[1])

    def run():
        hold_affinity(threads)
        pyscal3.set_num_threads(threads)
        pyscal3.find_neighbors(atoms, method="cutoff", cutoff=CUTOFF)
        values = np.asarray(pyscal3.steinhardt_parameter(atoms, order)[0]).reshape(-1)
        return f"rows={int(values.shape[0])}"
    return run, f"pyscal3 {pyscal3.__version__}"


def freud_call(op: str, threads: int, repetitions: int):
    import freud
    if op != "knn":
        raise ValueError(op)
    freud.parallel.set_num_threads(threads)
    xyz, side = frame(repetitions)
    box = freud.box.Box.cube(side)
    points = box.wrap(xyz - side / 2.0)
    query = freud.locality.AABBQuery(box, points)

    def run():
        hold_affinity(threads)
        found = query.query(points, {"num_neighbors": KNN, "exclude_ii": True})
        return f"pairs={len(found.toNeighborList())}"
    return run, f"freud {freud.__version__}"


def ovito_call(op: str, threads: int, repetitions: int):
    import numpy as np
    import ovito
    from ovito.data import DataCollection, NearestNeighborFinder, Particles, SimulationCell
    from ovito.modifiers import SteinhardtParametersModifier
    xyz, side = frame(repetitions)
    data = DataCollection()
    cell = SimulationCell(pbc=(True, True, True))
    cell[...] = [[side, 0.0, 0.0, 0.0], [0.0, side, 0.0, 0.0], [0.0, 0.0, side, 0.0]]
    data.objects.append(cell)
    particles = Particles()
    particles.create_property("Position", data=np.ascontiguousarray(xyz))
    data.objects.append(particles)
    version = f"ovito {ovito.version_string}"
    if op == "knn":
        def run():
            hold_affinity(threads)
            found, _ = NearestNeighborFinder(KNN, data).find_all()
            return f"pairs={int(found.size)}"
        return run, version
    if op != "q6":
        raise ValueError(op)
    tool = SteinhardtParametersModifier(
        neighbor_mode=SteinhardtParametersModifier.NeighborMode.Cutoff, cutoff=CUTOFF, max_l=6)

    def run():
        hold_affinity(threads)
        out = data.clone()
        out.apply(tool)
        q6 = np.asarray(out.particles["Steinhardt"])[:, -1]
        return f"mean={float(q6.mean()):.12f}"
    return run, version


def worker(args) -> int:
    # Read before a runtime loads: OMP_PROC_BIND narrows the main thread.
    cores, affinity = stats.cores_now()
    threads = args.threads
    rival = args.code != "dseams"
    if args.code == "dseams":
        call = dseams_call(args.lib, args.op, args.repetitions)
        version = args.label
    elif args.code == "mdapy":
        call, version = mdapy_call(args.op, threads, args.repetitions)
    elif args.code == "pyscal3":
        call, version = pyscal_call(args.op, threads, args.repetitions)
    elif args.code == "freud":
        call, version = freud_call(args.op, threads, args.repetitions)
    elif args.code == "ovito":
        call, version = ovito_call(args.op, threads, args.repetitions)
    else:
        raise ValueError(args.code)
    samples = []
    check = ""
    for index in range(WARMUP + args.repeats):
        if rival:
            hold_affinity(threads)
        t0 = time.perf_counter()
        check = call()
        elapsed = (time.perf_counter() - t0) * 1000.0
        if index >= WARMUP:
            samples.append(elapsed)
    print(json.dumps({"samples": samples, "check": check, "version": version,
                      "cores": cores, "affinity": affinity}), flush=True)
    return 0


def worker_env(code: str, threads: int, lib: str, extra: dict[str, str],
               placement: str) -> dict[str, str]:
    width = str(threads)
    if code == "dseams" and placement == "bound":
        env = stats.pin_env(threads)
    elif code == "dseams":
        env = dict(os.environ)
        for key in UNBOUND_DROP + ("OMP_DYNAMIC", "RAYON_NUM_THREADS"):
            env.pop(key, None)
        env.update({"OMP_NUM_THREADS": width, "OPENBLAS_NUM_THREADS": "1",
                    "MKL_NUM_THREADS": "1", "NUMEXPR_NUM_THREADS": "1"})
    else:
        # OMP_PROC_BIND would pin the rival's main thread to one cpu, and
        # its own thread pool would inherit that mask.
        env = dict(os.environ)
        for key in UNBOUND_DROP:
            env.pop(key, None)
        env.update({"OMP_NUM_THREADS": width, "OMP_DYNAMIC": "FALSE",
                    "RAYON_NUM_THREADS": width, "MDAPY_NUM_THREADS": width,
                    "PYSCAL_NUM_THREADS": width, "OVITO_THREAD_COUNT": width,
                    "OPENBLAS_NUM_THREADS": "1", "MKL_NUM_THREADS": "1",
                    "NUMEXPR_NUM_THREADS": "1"})
    env.pop("SEAMS_FAIR_WIDTH", None)
    env.update(extra)
    extra_paths = [str(Path(lib).resolve().parent),
                   str(ROOT / "bbdir-release/subprojects/readcon-core"),
                   str(ROOT / "bbdir-release/subprojects/readcon-db")]
    if env.get("LD_LIBRARY_PATH"):
        extra_paths.append(env["LD_LIBRARY_PATH"])
    env["LD_LIBRARY_PATH"] = os.pathsep.join(extra_paths)
    return env


def drive(code: str, op: str, threads: int, lib: str, label: str, repeats: int,
          repetitions: int, extra: dict[str, str], placement: str) -> dict:
    cmd = [sys.executable, __file__, "--worker", "--code", code, "--op", op,
           "--threads", str(threads), "--repeats", str(repeats), "--lib", lib,
           "--label", label, "--repetitions", str(repetitions)]
    if placement == "bound":
        cmd = ["taskset", "-c", stats.taskset_cores(threads)] + cmd
    done = subprocess.run(cmd, env=worker_env(code, threads, lib, extra, placement),
                          capture_output=True, text=True, check=False)
    if done.returncode != 0:
        raise RuntimeError(f"{code} {op} {threads}: {done.stderr.strip()[-2000:]}")
    lines = [line for line in done.stdout.splitlines() if line.startswith("{")]
    if not lines:
        raise RuntimeError(f"{code} {op} {threads}: no result line")
    return json.loads(lines[-1])


class Spinner:
    """scripts/omp_spinner.c on cpus 0 .. width-1, built in a scratch dir."""

    def __init__(self, width: int):
        self.width = width
        self.scratch = tempfile.TemporaryDirectory()
        self.proc = None

    def __enter__(self):
        compiler = os.environ.get("CC") or shutil.which("cc") or shutil.which("gcc")
        if compiler is None:
            raise RuntimeError("no C compiler for the spinner")
        exe = Path(self.scratch.name) / "omp_spinner"
        subprocess.run([compiler, "-O2", "-fopenmp", str(SPINNER), "-o", str(exe)],
                       check=True)
        env = dict(os.environ)
        for key in ("OMP_WAIT_POLICY", "OMP_PROC_BIND", "OMP_PLACES", "GOMP_SPINCOUNT"):
            env.pop(key, None)
        env["OMP_NUM_THREADS"] = str(self.width)
        self.proc = subprocess.Popen(["taskset", "-c", stats.taskset_cores(self.width),
                                      str(exe)], env=env, stdout=subprocess.PIPE, text=True)
        line = self.proc.stdout.readline()
        if not line.startswith("spinner"):
            raise RuntimeError("the spinner did not start")
        time.sleep(1.0)
        return self

    def __exit__(self, *exc):
        if self.proc is not None:
            self.proc.kill()
            self.proc.wait()
        self.scratch.cleanup()
        return False


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__,
                                     formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("--out", type=Path, default=ROOT / "bench/omp-scaling")
    parser.add_argument("--lib", default=str(ROOT / "bbdir-release/src/libyodaLib.so"))
    parser.add_argument("--label", default="after")
    parser.add_argument("--before", default="", help="a second libyodaLib.so to time")
    parser.add_argument("--before-label", default="before")
    parser.add_argument("--before-tree", default="", help="tree the --before build is from")
    parser.add_argument("--widths", default=",".join(str(w) for w in WIDTHS))
    parser.add_argument("--repetitions", default=",".join(str(r) for r in REPETITIONS),
                        help="ice Ic cells per edge; 8 * r**3 molecules")
    parser.add_argument("--ops", default=",".join(DSEAMS_OPS))
    parser.add_argument("--rival-widths", default=",".join(str(w) for w in WIDTHS))
    parser.add_argument("--rivals", default=",".join(RIVAL_OPS))
    parser.add_argument("--placements", default=",".join(PLACEMENTS))
    parser.add_argument("--unbound-ops", default=",".join(UNBOUND_OPS),
                        help="d-SEAMS operations that also run unbound")
    parser.add_argument("--contend-width", type=int, default=8,
                        help="width of the contended rows and the spinner; 0 skips them")
    parser.add_argument("--repeats", type=int, default=MIN_REPEATS)
    parser.add_argument("--worker", action="store_true")
    parser.add_argument("--code", default="dseams")
    parser.add_argument("--op", default="q4")
    parser.add_argument("--threads", type=int, default=1)
    args = parser.parse_args()
    if args.worker:
        args.repetitions = int(args.repetitions)
        return worker(args)
    if args.repeats < MIN_REPEATS:
        parser.error(f"--repeats must be at least {MIN_REPEATS}")

    widths = [int(w) for w in args.widths.split(",") if w]
    repetitions = [int(r) for r in args.repetitions.split(",") if r]
    ops = [op for op in args.ops.split(",") if op]
    unknown = sorted(set(ops) - set(DSEAMS_OPS))
    if unknown:
        parser.error(f"unknown operations {unknown}")
    rival_widths = [int(w) for w in args.rival_widths.split(",") if w]
    rivals = [r for r in args.rivals.split(",") if r]
    placements = [p for p in args.placements.split(",") if p]
    if set(placements) - set(PLACEMENTS):
        parser.error(f"unknown placements {sorted(set(placements) - set(PLACEMENTS))}")
    unbound_ops = [op for op in args.unbound_ops.split(",") if op in ops]
    contend = args.contend_width
    available = len(os.sched_getaffinity(0))
    if max(widths + rival_widths + [contend]) > available:
        parser.error(f"widths need {max(widths + rival_widths + [contend])} cpus, "
                     f"{available} allowed")

    builds = []
    if args.before:
        builds.append((args.before_label, args.before, {}))
    builds.append((args.label, args.lib, {}))
    fixed = (f"{args.label}-fixed", args.lib, {"SEAMS_FAIR_WIDTH": "0"})
    idle = []
    busy = []
    for place in placements:
        for reps in repetitions:
            for op in ops:
                if place == "unbound" and op not in unbound_ops:
                    continue
                for threads in widths:
                    for label, lib, extra in builds:
                        idle.append(("dseams", reps, op, threads, label, lib, extra, place))
                    if threads == contend and place == "bound":
                        idle.append(("dseams", reps, op, threads, *fixed, place))
                if contend:
                    for label, lib, extra in builds + ([fixed] if place == "bound" else []):
                        busy.append(("dseams", reps, op, contend, label, lib, extra, place))
        for code in rivals:
            for reps in repetitions:
                for op in RIVAL_OPS[code]:
                    for threads in rival_widths:
                        idle.append((code, reps, op, threads, code, args.lib, {}, place))
                    if contend:
                        busy.append((code, reps, op, contend, code, args.lib, {}, place))

    args.out.mkdir(parents=True, exist_ok=True)
    molecules = ",".join(str(8 * r ** 3) for r in repetitions)
    header = [
        "# omp-scaling sitting",
        f"# source_tree {stats.source_tree(ROOT)}",
        f"# worktree {stats.worktree_state(ROOT, TRACKED)}",
        f"# frame ice Ic molecules={molecules} cutoff={CUTOFF} max_ring={MAX_RING} "
        f"rdf2d_cutoff={RDF_CUTOFF} rdf2d_bin={RDF_BIN}",
        f"# warmup={WARMUP} repeats={args.repeats} one process per row",
        "# bound: taskset 0..width-1; d-SEAMS also OMP_PROC_BIND=close OMP_PLACES=cores",
        "# unbound: whole process mask; OMP_NUM_THREADS and the code's own thread count only; "
        + " ".join(UNBOUND_DROP) + " unset",
        f"# contended rows: scripts/omp_spinner.c on cpus 0-{contend - 1}, "
        f"OMP_NUM_THREADS={contend}, OMP_WAIT_POLICY unset"
        if contend else "# contended rows: none",
        f"# build {args.label}-fixed: build {args.label} with SEAMS_FAIR_WIDTH=0",
    ] + [f"# {line}" for line in stats.machine_lines(ROOT, commit=False)]
    if args.before:
        header.append(f"# build {args.before_label} tree {args.before_tree or 'unrecorded'}")
    header.append(f"# build {args.label} tree {stats.source_tree(ROOT)}")
    raw_lines = list(header)
    rows = []

    def run(jobs, load):
        for code, reps, op, threads, label, lib, extra, place in jobs:
            got = drive(code, op, threads, lib, label, args.repeats, reps, extra, place)
            summary = stats.summarize(got["samples"])
            row = {
                "code": code,
                "build": label,
                "version": got["version"],
                "molecules": 8 * reps ** 3,
                "operation": op,
                "threads": threads,
                "placement": place,
                "load": load,
                "cores": got["cores"],
                "affinity": got["affinity"],
                "n": summary["n"],
                "median_ms": f"{summary['median_ms']:.3f}",
                "q1_ms": f"{summary['q1_ms']:.3f}",
                "q3_ms": f"{summary['q3_ms']:.3f}",
                "iqr_ms": f"{summary['iqr_ms']:.3f}",
                "check": got["check"],
            }
            rows.append(row)
            raw_lines.append(
                f"{code} {label} molecules={8 * reps ** 3} {op} threads={threads} "
                f"placement={place} load={load} "
                f"cores={got['cores']} affinity={got['affinity']} {got['check']} samples_ms "
                f"{stats.fmt_samples(got['samples'])}")
            print(f"{code:8s} {label:12s} {8 * reps ** 3:7d} {op:16s} {threads} {place:7s} "
                  f"{load:9s} "
                  f"{row['median_ms']:>10s} ms", flush=True)

    run(idle, "idle")
    if busy:
        with Spinner(contend):
            run(busy, "contended")

    with open(args.out / "scaling.csv", "w", newline="") as handle:
        handle.write(f"# source_tree {stats.source_tree(ROOT)}\n")
        writer = csv.DictWriter(handle, fieldnames=list(rows[0].keys()), lineterminator="\n")
        writer.writeheader()
        writer.writerows(rows)
    (args.out / "raw-scaling.txt").write_text("\n".join(raw_lines) + "\n")
    return 0


if __name__ == "__main__":
    sys.exit(main())
