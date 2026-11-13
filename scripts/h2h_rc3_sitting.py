#!/usr/bin/env python3
# Copyright (c) d-SEAMS developers
"""One sitting of the 2.14.0 rc3 head-to-head.

Ideal ice Ic at 110592 and 884736 molecules, one thread and eight.
One warm-up, then six timed calls. The median of an even sample is the
mean of the two central values. The interquartile range is type-7.
Cores are the allowed list. The host name is not recorded.

Cages are the GenIce sI and sII frames, where cages exist. King's
shortest-path count is timed at every size whose first call returns
within 10 minutes.
"""

from __future__ import annotations

import argparse
import csv
import ctypes
import json
import math
import os
import shutil
import subprocess
import sys
import time
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))
import h2h_protocol as proto

ROOT = Path(__file__).resolve().parents[1]
LIB = ROOT / "bbdir-release" / "bindings" / "libseams_lang.so"
def plumed_library():
    """libplumed.so for the symfunc build.

    PLUMED_LIB names the library. Otherwise it is taken from the
    plumed executable on PATH, so a conda kernel earlier on the
    loader path is not used by accident.
    """
    named = os.environ.get("PLUMED_LIB")
    if named:
        return Path(named)
    exe = shutil.which("plumed")
    if exe:
        candidate = Path(exe).resolve().parent.parent / "lib" / "libplumed.so"
        if candidate.is_file():
            return candidate
    raise RuntimeError("set PLUMED_LIB to the symfunc libplumed.so")


PLUGIN = ROOT / "bbdir-release" / "engines" / "libSeamsPlumed.so"
WARMUP = 1
REPEATS = 6
CUTOFF = proto.CUTOFF_A


def cores_of(threads):
    return "0" if threads <= 1 else "0-" + str(threads - 1)


def ice(repetitions):
    pos, side = proto.ice_ic(repetitions)
    return np_pos(pos), float(side)


def np_pos(pos):
    import numpy as np
    return np.ascontiguousarray(pos, dtype=np.float64)


def read_dump(path):
    positions = []
    side = None
    with open(path) as handle:
        line = handle.readline()
        while line:
            if line.startswith("ITEM: NUMBER"):
                n = int(handle.readline())
            elif line.startswith("ITEM: BOX"):
                lo, hi = handle.readline().split()[:2]
                side = float(hi) - float(lo)
                handle.readline()
                handle.readline()
            elif line.startswith("ITEM: ATOMS"):
                for _ in range(n):
                    parts = handle.readline().split()
                    positions.append(tuple(float(v) for v in parts[-3:]))
                break
            line = handle.readline()
    return np_pos(positions), float(side)


def summarize(samples):
    return proto.summarize(samples)


def time_calls(fn, warmup=WARMUP, repeats=REPEATS, max_call_s=None):
    for _ in range(warmup):
        t0 = time.perf_counter()
        fn()
        if max_call_s is not None and time.perf_counter() - t0 > max_call_s:
            raise TimeoutError("warm-up exceeded the 10 min cap")
    samples = []
    last = None
    for _ in range(repeats):
        t0 = time.perf_counter()
        last = fn()
        elapsed = time.perf_counter() - t0
        if max_call_s is not None and elapsed > max_call_s:
            raise TimeoutError("timed call exceeded the 10 min cap")
        samples.append(elapsed * 1e3)
    return samples, last


class Seams:
    def __init__(self, threads):
        os.environ["OMP_NUM_THREADS"] = str(threads)
        os.environ["OMP_DYNAMIC"] = "FALSE"
        self.lib = ctypes.CDLL(str(LIB))
        lib = self.lib
        lib.seams_lang_version.restype = ctypes.c_char_p
        lib.seams_last_error.restype = ctypes.c_char_p
        lib.seams_frame_from_arrays.restype = ctypes.c_void_p
        lib.seams_frame_from_arrays.argtypes = [
            ctypes.POINTER(ctypes.c_double), ctypes.c_int,
            ctypes.POINTER(ctypes.c_double),
            ctypes.POINTER(ctypes.c_int), ctypes.POINTER(ctypes.c_int)]
        lib.seams_frame_free.argtypes = [ctypes.c_void_p]
        lib.seams_frame_neighbours.argtypes = [
            ctypes.c_void_p, ctypes.c_double, ctypes.c_int, ctypes.c_char_p]
        lib.seams_frame_neighbour_pairs.restype = ctypes.c_int
        lib.seams_frame_neighbour_pairs.argtypes = [ctypes.c_void_p]
        lib.seams_frame_steinhardt.argtypes = [
            ctypes.c_void_p, ctypes.c_double, ctypes.c_int, ctypes.c_int,
            ctypes.POINTER(ctypes.c_double), ctypes.POINTER(ctypes.c_double)]
        lib.seams_frame_chill_plus.argtypes = [
            ctypes.c_void_p, ctypes.c_double, ctypes.c_int, ctypes.c_char_p,
            ctypes.POINTER(ctypes.c_int), ctypes.POINTER(ctypes.c_int)]
        lib.seams_frame_rings.argtypes = [
            ctypes.c_void_p, ctypes.c_double, ctypes.c_int, ctypes.c_char_p,
            ctypes.c_int, ctypes.POINTER(ctypes.c_int), ctypes.c_int]
        lib.seams_frame_symmetry.argtypes = [
            ctypes.c_void_p, ctypes.c_double, ctypes.POINTER(ctypes.c_int)]
        lib.seams_frame_signature.argtypes = [
            ctypes.c_void_p, ctypes.c_double, ctypes.c_int, ctypes.c_char_p,
            ctypes.c_char_p, ctypes.POINTER(ctypes.c_int), ctypes.POINTER(ctypes.c_int)]
        self.version = lib.seams_lang_version().decode()

    def frame(self, pos, side):
        import numpy as np
        n = len(pos)
        box = np.array([0, 0, 0, side, side, side, 0, 0, 0], dtype=np.float64)
        types = np.ones(n, dtype=np.int32)
        ids = np.arange(1, n + 1, dtype=np.int32)
        handle = self.lib.seams_frame_from_arrays(
            pos.ctypes.data_as(ctypes.POINTER(ctypes.c_double)), n,
            box.ctypes.data_as(ctypes.POINTER(ctypes.c_double)),
            types.ctypes.data_as(ctypes.POINTER(ctypes.c_int)),
            ids.ctypes.data_as(ctypes.POINTER(ctypes.c_int)))
        if not handle:
            raise RuntimeError(self.lib.seams_last_error().decode())
        return handle, n

    def neighbours(self, pos, side):
        handle, _n = self.frame(pos, side)

        def run():
            if self.lib.seams_frame_neighbours(handle, CUTOFF, 1, b"cutoff") != 0:
                raise RuntimeError(self.lib.seams_last_error().decode())
            return int(self.lib.seams_frame_neighbour_pairs(handle))

        samples, last = time_calls(run)
        self.lib.seams_frame_free(handle)
        return samples, last, self.version

    def ql(self, pos, side, order):
        import numpy as np
        handle, n = self.frame(pos, side)
        ql = np.empty(n)
        bar = np.empty(n)

        def run():
            if self.lib.seams_frame_steinhardt(
                    handle, CUTOFF, 1, order,
                    ql.ctypes.data_as(ctypes.POINTER(ctypes.c_double)),
                    bar.ctypes.data_as(ctypes.POINTER(ctypes.c_double))) != 0:
                raise RuntimeError(self.lib.seams_last_error().decode())
            return ql.copy()

        samples, last = time_calls(run)
        self.lib.seams_frame_free(handle)
        return samples, last, self.version

    def chill(self, pos, side):
        import numpy as np
        handle, n = self.frame(pos, side)
        labels = np.empty(n, dtype=np.int32)
        counts = (ctypes.c_int * 9)()

        def run():
            if self.lib.seams_frame_chill_plus(
                    handle, CUTOFF, 1, b"cutoff",
                    labels.ctypes.data_as(ctypes.POINTER(ctypes.c_int)), counts) != 0:
                raise RuntimeError(self.lib.seams_last_error().decode())
            return int(np.sum(labels == 0))

        samples, last = time_calls(run)
        self.lib.seams_frame_free(handle)
        return samples, last, self.version

    def rings(self, pos, side):
        handle, _n = self.frame(pos, side)
        counts = (ctypes.c_int * 7)()

        def run():
            if self.lib.seams_frame_rings(
                    handle, CUTOFF, 1, b"cutoff", 6, counts, 7) != 0:
                raise RuntimeError(self.lib.seams_last_error().decode())
            return [int(counts[i]) for i in range(3, 7)]

        samples, last = time_calls(run)
        self.lib.seams_frame_free(handle)
        return samples, last, self.version

    def symmetry(self, pos, side):
        import numpy as np
        handle, n = self.frame(pos, side)
        labels = np.empty(n, dtype=np.int32)

        def run():
            if self.lib.seams_frame_symmetry(
                    handle, 0.12, labels.ctypes.data_as(ctypes.POINTER(ctypes.c_int))) != 0:
                raise RuntimeError(self.lib.seams_last_error().decode())
            return int(np.sum(labels > 0))

        samples, last = time_calls(run)
        self.lib.seams_frame_free(handle)
        return samples, last, self.version

    def cages(self, pos, side):
        handle, _n = self.frame(pos, side)
        specs = (b"512", b"51262", b"51264")

        def one(spec):
            nc = ctypes.c_int()
            na = ctypes.c_int()
            if self.lib.seams_frame_signature(
                    handle, CUTOFF, 1, b"cutoff", spec, ctypes.byref(nc), ctypes.byref(na)) != 0:
                raise RuntimeError(self.lib.seams_last_error().decode())
            return int(nc.value)

        def run():
            return tuple(one(spec) for spec in specs)

        samples, last = time_calls(run)
        self.lib.seams_frame_free(handle)
        return samples, last, self.version


def worker(kind, molecules, threads, frame_path, calls=0):
    import numpy as np
    os.environ["FASTATOMSTRUCT_NO_TELEMETRY"] = "1"
    os.environ["FASTATOMSTRUCT_SUPPRESS_TELEMETRY_NOTICE"] = "1"
    os.environ["OMP_NUM_THREADS"] = str(threads)
    data = np.load(frame_path)
    pos = data["positions"]
    side = float(data["side"])
    seams = Seams(threads)
    detail = ""
    agree = ""
    if kind == "dseams-neighbours":
        samples, last, version = seams.neighbours(pos, side)
        detail = f"pairs={last}"
    elif kind == "vesin-neighbours":
        from vesin import NeighborList
        import vesin
        calc = NeighborList(cutoff=CUTOFF, full_list=True, n_threads=threads, algorithm="cell_list")
        cell = np.eye(3) * side

        def run():
            left, right = calc.compute(points=pos, box=cell, periodic=True, quantities="ij", copy=True)
            return int(len(left))

        samples, last = time_calls(run)
        version = vesin.__version__
        detail = f"pairs={last} full_list=1 algorithm=cell_list n_threads={threads} copy=1"
    elif kind in ("dseams-q4", "dseams-q6"):
        order = 4 if kind.endswith("q4") else 6
        samples, last, version = seams.ql(pos, side, order)
        detail = f"mean={float(np.mean(last)):.8f}"
        agree = last
    elif kind in ("mdapy-q4", "mdapy-q6"):
        import mdapy
        from mdapy import System
        from mdapy.analysis import SteinhardtParameters
        mdapy.set_num_threads(threads)
        order = 4 if kind.endswith("q4") else 6
        tool = SteinhardtParameters(degrees=(order,), cutoff=CUTOFF, average=False)
        system = System(pos, box=side)

        def run():
            return np.asarray(tool.compute(system).q).reshape(-1)

        samples, last = time_calls(run)
        version = mdapy.__version__
        detail = f"mean={float(np.mean(last)):.8f} threads={mdapy.get_num_threads()}"
        agree = last
    elif kind in ("freud-q4", "freud-q6"):
        import freud
        freud.parallel.set_num_threads(threads)
        order = 4 if kind.endswith("q4") else 6
        box = freud.box.Box.cube(side)
        ql = freud.order.Steinhardt(order)

        def run():
            ql.compute((box, pos), neighbors={"r_max": CUTOFF, "exclude_ii": True})
            return np.asarray(ql.particle_order, dtype=np.float64)

        samples, last = time_calls(run)
        version = freud.__version__
        detail = f"mean={float(np.mean(last)):.8f}"
        agree = last
    elif kind in ("plumed-q4", "plumed-q6"):
        samples, last, version = plumed_ql(pos, side, 4 if kind.endswith("q4") else 6, threads)
        detail = f"mean={float(np.mean(last)):.8f} norm=sqrt(4pi/(2l+1))"
        agree = last
    elif kind == "dseams-chill":
        samples, last, version = seams.chill(pos, side)
        detail = f"cubic={last}"
    elif kind == "mdapy-chill":
        import mdapy
        from mdapy import System
        from mdapy.identify import ChillPlus
        mdapy.set_num_threads(threads)
        tool = ChillPlus(cutoff=CUTOFF)
        system = System(pos, box=side)

        def run():
            result = tool.compute(system)
            return int(np.sum(np.asarray(result.structure) == 2))

        samples, last = time_calls(run)
        version = mdapy.__version__
        detail = f"cubic={last}"
    elif kind == "ovito-chill":
        samples, last, version = ovito_mod(pos, side, threads, "chill")
        detail = f"cubic={last}"
    elif kind == "dseams-rings":
        samples, last, version = seams.rings(pos, side)
        detail = f"census={last}"
    elif kind == "mdapy-rings":
        import mdapy
        from mdapy import System
        from mdapy.analysis import RingStatistics
        mdapy.set_num_threads(threads)
        system = System(pos, box=side)

        def run():
            result = RingStatistics(max_size=6, min_size=3, cutoff=CUTOFF).compute(system)
            return [int(result.counts[k]) for k in range(3, 7)]

        samples, last = time_calls(run)
        version = mdapy.__version__
        detail = f"census={last}"
    elif kind == "ovito-rings":
        samples, last, version = ovito_mod(pos, side, threads, "rings")
        detail = f"sink={last}"
    elif kind == "ringsstat-rings":
        from ase import Atoms
        import rings_statistics_matter
        atoms = Atoms(symbols=["O"] * len(pos), positions=pos, cell=[side, side, side], pbc=True)

        def run():
            # maxlvl 8: a six-ring midpoint is at depth 3, and the kernel
            # only tries midpoints closer than maxlvl/2.
            counts, _rings = rings_statistics_matter.ring_statistics(
                atoms, cutoff=CUTOFF, maxlvl=8, supercell=False,
                deduplicate=True, backend="rust", num_threads=threads)
            values = [int(round(float(counts[k]))) for k in range(2, 6)]
            return values

        samples, last = time_calls(run)
        version = "0.3.0"
        detail = f"census={last}"
    elif kind == "dseams-ptm":
        samples, last, version = seams.symmetry(pos, side)
        detail = f"labelled={last}"
    elif kind == "ovito-ptm":
        samples, last, version = ovito_mod(pos, side, threads, "ptm")
        detail = f"labelled={last}"
    elif kind == "dseams-cages":
        samples, last, version = seams.cages(pos, side)
        detail = f"512/51262/51264={last}"
    elif kind == "clavis-cages":
        import clavis_engine as engine
        # Hydrogens are an ice-rule orientation of the four nearest
        # oxygens, built before the timer. TRACE reads that graph.
        text, box_nm = ice_rule_waters(pos, side)

        def run():
            raw = engine.analyze_waters(text, "[]", box_nm, "trace", 0.35, 35.0)
            parsed = json.loads(raw)
            return (len(parsed["cages_512"]), len(parsed["cages_51262"]), len(parsed["cages_51264"]))

        samples, last = time_calls(run)
        version = str(engine.version())
        detail = f"512/51262/51264={last}"
    elif kind == "king":
        import fastatomstruct as fs
        fs.disable_telemetry()
        fs.set_num_threads(threads)
        elements = [fs.Element.O] * len(pos)
        atoms = fs.Atoms(pos, np.eye(3) * side, elements)
        cutoff = fs.analysis.distance.CutoffMode.Fixed(CUTOFF)
        method = fs.analysis.rings.RingMethod.King

        def run():
            sizes, counts = fs.analysis.rings.ring_size_distribution(atoms, cutoff, 6, method)
            return [int(round(float(v))) for v in counts]

        if calls:
            samples, last = time_calls(run, warmup=0, repeats=calls, max_call_s=600)
        else:
            samples, last = time_calls(run, max_call_s=600)
        version = fs.__version__
        detail = f"census={last}"
    else:
        raise SystemExit(kind)
    payload = {
        "kind": kind,
        "version": version,
        "molecules": int(len(pos)),
        "threads": threads,
        "cpus": proto.cpus_allowed(),
        "samples_ms": samples,
        "detail": detail,
    }
    if isinstance(agree, np.ndarray):
        payload["values"] = agree.tolist()
    json.dump(payload, sys.stdout)
    sys.stdout.write("\n")


def plumed_ql(pos, side, order, threads):
    import numpy as np
    os.environ["OMP_NUM_THREADS"] = str(threads)
    # The conda prefix's RPATH beats LD_LIBRARY_PATH and would load an
    # older kernel, which has no Q6. Open the matching kernel by path first.
    library = plumed_library()
    kernel = library.with_name("libplumedKernel.so")
    ctypes.CDLL(str(kernel), mode=ctypes.RTLD_GLOBAL)
    lib = ctypes.CDLL(str(library), mode=ctypes.RTLD_GLOBAL)
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

    n = len(pos)
    positions = np.ascontiguousarray(pos / 10.0)
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
    cmd("readInputLine", f"obs: Q{order} SPECIES=1-{n} SWITCH={{RATIONAL R_0=0.35 D_MAX=0.35 NN=100 MM=200}}")
    rank = np.zeros(1, dtype=np.intc)
    cmd("getDataRank obs", rank.ctypes.data_as(ctypes.c_void_p))
    shape = np.zeros(int(rank[0]), dtype=np.intc)
    cmd("getDataShape obs", shape.ctypes.data_as(ctypes.c_void_p))
    data = np.zeros(tuple(int(v) for v in shape) or (n,), dtype=np.float64)
    cmd("setMemoryForData obs", data.ctypes.data_as(ctypes.c_void_p))
    step = ctypes.c_int(1)
    scale = math.sqrt(4.0 * math.pi / (2 * order + 1))

    def run():
        step.value += 1
        cmd("setStep", ctypes.byref(step))
        cmd("setBox", box.ctypes.data_as(ctypes.c_void_p))
        cmd("setPositions", positions.ctypes.data_as(ctypes.c_void_p))
        cmd("setMasses", mass.ctypes.data_as(ctypes.c_void_p))
        cmd("setCharges", charge.ctypes.data_as(ctypes.c_void_p))
        cmd("setForces", forces.ctypes.data_as(ctypes.c_void_p))
        cmd("setVirial", virial.ctypes.data_as(ctypes.c_void_p))
        cmd("calc")
        return np.array(data, copy=True).reshape(-1) * scale

    samples, last = time_calls(run)
    return samples, last, "2.10.2+symfunc"


def ovito_mod(pos, side, threads, which):
    import numpy as np
    os.environ["OVITO_THREAD_COUNT"] = str(threads)
    os.environ["QT_QPA_PLATFORM"] = "offscreen"
    from ovito.data import DataCollection
    from ovito.modifiers import (
        ChillPlusModifier, CreateBondsModifier, FindRingsModifier,
        PolyhedralTemplateMatchingModifier)
    import ovito

    def run():
        data = DataCollection()
        data.create_cell(
            matrix=[[side, 0, 0, 0], [0, side, 0, 0], [0, 0, side, 0]],
            pbc=(True, True, True))
        parts = data.create_particles(count=len(pos))
        parts.create_property("Position", data=np.ascontiguousarray(pos))
        if which == "chill":
            data.apply(ChillPlusModifier(cutoff=CUTOFF))
            values = np.asarray(data.particles["Structure Type"], dtype=np.int32)
            return int(np.sum(values == int(ChillPlusModifier.Type.CUBIC_ICE)))
        if which == "rings":
            data.apply(CreateBondsModifier(cutoff=CUTOFF))
            data.apply(FindRingsModifier(minimum_ring_size=3, maximum_ring_size=6))
            hist = data.tables["ring-size-histogram"]
            found = dict(zip(
                [int(v) for v in np.asarray(hist["Ring Size"])],
                [int(v) for v in np.asarray(hist["Count"])]))
            return [found.get(size, 0) for size in range(3, 7)]
        mod = PolyhedralTemplateMatchingModifier()
        for structure in mod.structures:
            structure.enabled = structure.id in {6, 7}
        data.apply(mod)
        values = np.asarray(data.particles["Structure Type"], dtype=np.int32)
        return int(np.sum(values == 6))

    samples, last = time_calls(run)
    version = getattr(ovito, "version_string", None)
    if version is None:
        version = ".".join(str(v) for v in tuple(ovito.version)[:3])
    return samples, last, version


def launch(kind, molecules, threads, frame, timeout, calls=0):
    cmd = [
        "taskset", "-c", cores_of(threads), sys.executable, __file__,
        "--worker", kind, "--molecules", str(molecules), "--threads", str(threads),
        "--frame", str(frame),
    ]
    if calls:
        cmd.extend(["--calls", str(calls)])
    env = os.environ.copy()
    env["OMP_NUM_THREADS"] = str(threads)
    env["OMP_DYNAMIC"] = "FALSE"
    env["FASTATOMSTRUCT_NO_TELEMETRY"] = "1"
    done = subprocess.run(cmd, cwd=ROOT, env=env, text=True, capture_output=True, timeout=timeout)
    if done.returncode != 0:
        raise RuntimeError(done.stderr[-2000:] or done.stdout[-2000:])
    return json.loads(done.stdout.strip().splitlines()[-1])


def row_of(payload, operation, code):
    summary = summarize(payload["samples_ms"])
    return {
        "code": code,
        "version": payload["version"],
        "molecules": payload["molecules"],
        "operation": operation,
        "threads": payload["threads"],
        "cores": payload["cpus"],
        "n": summary["n"],
        "median_ms": f"{summary['median_ms']:.3f}",
        "iqr_ms": f"{summary['iqr_ms']:.3f}",
        "detail": payload["detail"],
        "samples_ms": " ".join(f"{v:.3f}" for v in payload["samples_ms"]),
    }


def max_abs(left, right):
    import numpy as np
    a = np.asarray(left, dtype=np.float64).reshape(-1)
    b = np.asarray(right, dtype=np.float64).reshape(-1)
    n = min(len(a), len(b))
    if n == 0:
        return float("nan")
    return float(np.max(np.abs(a[:n] - b[:n])))


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--worker", default="")
    parser.add_argument("--molecules", type=int, default=0)
    parser.add_argument("--threads", type=int, default=1)
    parser.add_argument("--frame", default="")
    parser.add_argument("--out", default=str(ROOT / "benchmarks" / "rc3-sitting.csv"))
    parser.add_argument("--log", default=str(ROOT / "benchmarks" / "rc3-sitting.txt"))
    parser.add_argument("--vesin-gap", action="store_true")
    parser.add_argument("--calls", type=int, default=0)
    args = parser.parse_args()
    if args.vesin_gap:
        vesin_gap(ROOT / "benchmarks" / "rc3-vesin-gap.txt")
        return
    if args.worker:
        worker(args.worker, args.molecules, args.threads, args.frame, args.calls)
        return
    import numpy as np
    scratch = ROOT / "bbdir" / "rc3-sitting"
    scratch.mkdir(parents=True, exist_ok=True)
    frames = {}
    for rep, count in ((24, 110592), (48, 884736)):
        pos, side = proto.ice_ic(rep)
        path = scratch / f"ic-{count}.npz"
        np.savez(path, positions=np.ascontiguousarray(pos), side=np.float64(side))
        frames[count] = path
    for name in ("genice_sI", "genice_sII"):
        pos, side = read_dump(ROOT / "input" / "traj" / f"{name}.lammpstrj")
        path = scratch / f"{name}.npz"
        np.savez(path, positions=pos, side=np.float64(side))
        frames[name] = path
    # TRACE drops a ring that crosses the periodic boundary. The sI unit
    # cell is smaller than its cages, so the 2x2x2 replicate is the frame
    # on which both codes count the same polyhedra.
    base, base_side = read_dump(ROOT / "input" / "traj" / "genice_sI.lammpstrj")
    shifts = np.stack(np.meshgrid(*[np.arange(2)] * 3, indexing="ij"), axis=-1).reshape(-1, 3)
    tiled = np.ascontiguousarray((base[None, :, :] + shifts[:, None, :] * base_side).reshape(-1, 3))
    path = scratch / "genice_sI_222.npz"
    np.savez(path, positions=tiled, side=np.float64(base_side * 2))
    frames["genice_sI_222"] = path
    fields = ["code", "version", "molecules", "operation", "threads", "cores",
              "n", "median_ms", "iqr_ms", "detail", "samples_ms"]
    out = Path(args.out)
    log = Path(args.log)
    out.parent.mkdir(parents=True, exist_ok=True)
    pairs = []
    for molecules in (110592, 884736):
        for threads in (1, 8):
            pairs.append((molecules, threads, "neighbours",
                          [("d-SEAMS", "dseams-neighbours"), ("vesin", "vesin-neighbours")]))
            for op, order in (("q4", "q4"), ("q6", "q6")):
                pairs.append((molecules, threads, op, [
                    ("d-SEAMS", f"dseams-{order}"),
                    ("mdapy", f"mdapy-{order}"),
                    ("freud", f"freud-{order}"),
                    ("PLUMED", f"plumed-{order}"),
                ]))
            pairs.append((molecules, threads, "CHILL+", [
                ("d-SEAMS", "dseams-chill"), ("mdapy", "mdapy-chill"), ("OVITO", "ovito-chill")]))
            pairs.append((molecules, threads, "rings", [
                ("d-SEAMS", "dseams-rings"), ("OVITO", "ovito-rings"),
                ("mdapy", "mdapy-rings"), ("rings-statistics-matter", "ringsstat-rings")]))
            pairs.append((molecules, threads, "PTM", [
                ("d-SEAMS", "dseams-ptm"), ("OVITO", "ovito-ptm")]))
    gap = subprocess.run(
        ["taskset", "-c", "0-7", sys.executable, __file__, "--vesin-gap"],
        cwd=ROOT, env=os.environ.copy(), check=False)
    if gap.returncode != 0:
        raise SystemExit("vesin gap probe failed")
    with out.open("w", newline="") as handle, log.open("w") as raw:
        writer = csv.DictWriter(handle, fieldnames=fields)
        writer.writeheader()
        raw.write(f"cores {os.cpu_count()} warmup {WARMUP} n {REPEATS}\n")
        for molecules, threads, operation, codes in pairs:
            got = {}
            for code, kind in codes:
                raw.write(f"begin {kind} molecules={molecules} threads={threads}\n")
                raw.flush()
                timeout = 2400 if kind in (
                    "ringsstat-rings", "plumed-q4", "plumed-q6", "mdapy-rings", "ovito-ptm"
                ) else 900
                try:
                    payload = launch(kind, molecules, threads, frames[molecules], timeout=timeout)
                except Exception as exc:
                    raw.write(f"FAIL {kind} {exc}\n")
                    raw.flush()
                    continue
                record = row_of(payload, operation, code)
                writer.writerow(record)
                handle.flush()
                got[code] = payload
                raw.write(
                    f"row {code} {operation} mol={record['molecules']} t={record['threads']} "
                    f"cores={record['cores']} n={record['n']} median={record['median_ms']} "
                    f"iqr={record['iqr_ms']} {record['detail']}\n")
                raw.flush()
            if operation in ("q4", "q6") and "d-SEAMS" in got and "mdapy" in got:
                gap = max_abs(got["d-SEAMS"].get("values", []), got["mdapy"].get("values", []))
                raw.write(f"agree {operation} mdapy max_abs={gap:.3e}\n")
        for name in ("genice_sI", "genice_sI_222", "genice_sII"):
            for threads in (1, 8):
                for code, kind in (("d-SEAMS", "dseams-cages"), ("ClavisCage", "clavis-cages")):
                    raw.write(f"begin {kind} frame={name} threads={threads}\n")
                    raw.flush()
                    try:
                        payload = launch(kind, 0, threads, frames[name], timeout=900)
                    except Exception as exc:
                        raw.write(f"FAIL {kind} {name} {exc}\n")
                        continue
                    record = row_of(payload, "cages", code)
                    record["molecules"] = payload["molecules"]
                    writer.writerow(record)
                    handle.flush()
                    raw.write(
                        f"row {code} cages mol={record['molecules']} t={threads} "
                        f"cores={payload['cpus']} median={record['median_ms']} {record['detail']}\n")
        # King, growing sizes, 10 minute cap on the whole timed set.
        for threads in (1, 8):
            for rep in (5, 8, 11, 14, 16, 20, 24):
                count = proto.molecules_for(rep)
                path = scratch / f"king-{count}.npz"
                if not path.is_file():
                    pos, side = proto.ice_ic(rep)
                    np.savez(path, positions=np.ascontiguousarray(pos), side=np.float64(side))
                raw.write(f"begin king molecules={count} threads={threads}\n")
                raw.flush()
                try:
                    ref = launch("dseams-rings", count, threads, path, timeout=600)
                    # One process per call so a size past the cap is killed
                    # at 10 minutes and is not written down.
                    launch("king", count, threads, path, timeout=600, calls=1)
                    parts = []
                    last = None
                    for _ in range(REPEATS):
                        last = launch("king", count, threads, path, timeout=600, calls=1)
                        parts.append(last["samples_ms"][0])
                    last["samples_ms"] = parts
                except Exception as exc:
                    raw.write(f"king stopped molecules={count} threads={threads} within the 10 min cap: {exc}\n")
                    raw.flush()
                    break
                for code, payload, operation in (
                        ("d-SEAMS", ref, "rings"), ("fastatomstruct", last, "king")):
                    record = row_of(payload, operation, code)
                    writer.writerow(record)
                    handle.flush()
                    raw.write(
                        f"row {code} {operation} mol={count} t={threads} "
                        f"cores={payload['cpus']} median={record['median_ms']} {record['detail']}\n")
    print(f"wrote {out}")


def ice_rule_waters(pos, side):
    """JSON waters for TRACE.

    The four nearest oxygens are a symmetric graph on these crystals.
    An Eulerian circuit orients each water toward two of them, which is
    an ice rule. Each hydrogen sits 0.96 A along that bond. Coordinates
    in the JSON are nanometres.
    """
    import numpy as np
    from vesin import NeighborList
    cell = np.eye(3) * side
    left, right, dist = NeighborList(cutoff=CUTOFF, full_list=True).compute(
        points=pos, box=cell, periodic=True, quantities="ijd", copy=True)
    buckets = [[] for _ in pos]
    for i, j, distance in zip(left, right, dist):
        buckets[int(i)].append((float(distance), int(j)))
    adj = []
    for rows in buckets:
        rows.sort()
        chosen = [j for _distance, j in rows[:4]]
        if len(chosen) != 4:
            raise RuntimeError("hydrate oxygen does not have four neighbours inside 3.5 A")
        adj.append(chosen)
    n = len(pos)
    unused = [set() for _ in range(n)]
    edges = {}
    eid = 0
    for i, others in enumerate(adj):
        for j in others:
            if i < j:
                edges[eid] = (i, j)
                unused[i].add(eid)
                unused[j].add(eid)
                eid += 1
    donated = [[] for _ in range(n)]
    for start in range(n):
        if not unused[start]:
            continue
        stack = [start]
        circuit = []
        while stack:
            vertex = stack[-1]
            if unused[vertex]:
                edge = unused[vertex].pop()
                a, b = edges[edge]
                other = b if a == vertex else a
                unused[other].discard(edge)
                stack.append(other)
            else:
                circuit.append(stack.pop())
        for a, b in zip(circuit, circuit[1:]):
            donated[a].append(b)
    waters = []
    for i, point in enumerate(pos):
        if len(donated[i]) != 2:
            raise RuntimeError(f"ice-rule orientation failed at water {i}")
        hydrogens = []
        for j in donated[i]:
            delta = pos[j] - point
            delta -= side * np.round(delta / side)
            hydrogens.append(point + delta / np.linalg.norm(delta) * 0.96)
        waters.append({
            "oxygen": (point / 10.0).tolist(),
            "hydrogen1": (hydrogens[0] / 10.0).tolist(),
            "hydrogen2": (hydrogens[1] / 10.0).tolist(),
        })
    return json.dumps(waters), [side / 10.0, side / 10.0, side / 10.0]


def vesin_gap(log_path):
    """Why one sitting saw ~20 ms and another ~50 ms at 110592 molecules."""
    import numpy as np
    from vesin import NeighborList
    import vesin
    pos, side = proto.ice_ic(24)
    pos = np.ascontiguousarray(pos)
    cell = np.eye(3) * side
    fortran = np.asfortranarray(pos)
    lines = [f"vesin {vesin.__version__} molecules {len(pos)} cutoff {CUTOFF}"]
    configs = [
        ("full cell copy ij t8", dict(full_list=True, n_threads=8, algorithm="cell_list"), True, "ij", pos),
        ("full cell nocopy ij t8", dict(full_list=True, n_threads=8, algorithm="cell_list"), False, "ij", pos),
        ("half cell copy ij t8", dict(full_list=False, n_threads=8, algorithm="cell_list"), True, "ij", pos),
        ("full auto copy ij t8", dict(full_list=True, n_threads=8, algorithm="auto"), True, "ij", pos),
        ("full cell copy ijS t8", dict(full_list=True, n_threads=8, algorithm="cell_list"), True, "ijS", pos),
        ("full cell copy ij t1", dict(full_list=True, n_threads=1, algorithm="cell_list"), True, "ij", pos),
        ("full cell copy ij t0", dict(full_list=True, n_threads=0, algorithm="cell_list"), True, "ij", pos),
        ("full cell fortran t8", dict(full_list=True, n_threads=8, algorithm="cell_list"), True, "ij", fortran),
    ]
    for name, opts, copy, quantities, points in configs:
        try:
            calc = NeighborList(cutoff=CUTOFF, **opts)

            def run(calc=calc, points=points, quantities=quantities, copy=copy):
                left, _right = calc.compute(
                    points=points, box=cell, periodic=True, quantities=quantities, copy=copy)[:2]
                return len(left)

            run()
            samples = []
            for _ in range(6):
                t0 = time.perf_counter()
                n = run()
                samples.append((time.perf_counter() - t0) * 1e3)
            summary = proto.summarize(samples)
            lines.append(
                f"{name} pairs={n} median_ms={summary['median_ms']:.3f} "
                f"iqr={summary['iqr_ms']:.3f} cpus={proto.cpus_allowed()} "
                f"samples={' '.join(f'{v:.3f}' for v in samples)}")
        except Exception as exc:
            lines.append(f"{name} FAILED {exc}")
    Path(log_path).write_text("\n".join(lines) + "\n")
    print("\n".join(lines))


if __name__ == "__main__":
    main()
