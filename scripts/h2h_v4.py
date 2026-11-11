#!/usr/bin/env python3
"""Release head-to-head on the fast tree, one core and eight cores.

CHILL+ against mdapy and OVITO. Rings against OVITO Find rings, matscipy
and cycless. Ideal ice Ic, cutoff 3.5 A, rings of size 3-6. One warm-up,
six timed calls, median and type-7 IQR. Each code uses its own thread
count: d-SEAMS sets OMP_NUM_THREADS and RAYON_NUM_THREADS, mdapy calls
set_num_threads, OVITO reads OVITO_THREAD_COUNT before import. matscipy
and cycless have no parallel ring API; they run with OpenMP set to the
same width.
"""

from __future__ import annotations

import argparse
import csv
import json
import os
import subprocess
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
HERE = Path(__file__).resolve()
FLEXI = (
    "/cvmfs/software.eessi.io/versions/2026.06/software/linux/x86_64/intel/"
    "sapphirerapids/software/FlexiBLAS/3.5.0-GCC-15.2.0/lib"
)
HOST_LIB = "/usr/lib/x86_64-linux-gnu"

FIELDS = [
    "threads", "code", "version", "operation", "molecules", "repetitions",
    "warmup", "n", "median_ms", "q1_ms", "q3_ms", "iqr_ms", "agreement",
    "detail", "samples_ms", "cpus",
]


def lib_env(lib: str) -> str:
    build = Path(lib).resolve().parents[1]
    parts = [
        str(build / "bindings"),
        str(build / "src"),
        str(build / "subprojects/readcon-core"),
        str(build / "subprojects/readcon-db"),
        str(build / "subprojects/linkcell"),
        str(build / "subprojects/vesin"),
        str(build / "subprojects/minimage"),
        FLEXI,
    ]
    return os.pathsep.join(parts)


def parse_args(argv):
    parser = argparse.ArgumentParser()
    parser.add_argument("--lib")
    parser.add_argument("--out", default="benchmarks/h2h-watch-v4.csv")
    parser.add_argument("--python")
    parser.add_argument("--threads", type=int, default=1)
    parser.add_argument("--worker", choices=(
        "dseams", "mdapy", "ovito", "matscipy", "cycless"))
    parser.add_argument("--operation", choices=("chill+", "rings"))
    parser.add_argument("--frame")
    parser.add_argument("--skip-cycless", action="store_true")
    parser.add_argument("--only", choices=("mdapy", "ovito", "matscipy", "cycless"))
    parser.add_argument("--sizes", default="5,11,24")
    return parser.parse_args(argv)


def child_env(kind, threads, lib):
    env = os.environ.copy()
    width = str(threads)
    env["OMP_NUM_THREADS"] = width
    env["OMP_DYNAMIC"] = "FALSE"
    env["OMP_PROC_BIND"] = "close"
    env["MKL_NUM_THREADS"] = "1"
    env["OPENBLAS_NUM_THREADS"] = "1"
    env["NUMEXPR_NUM_THREADS"] = "1"
    env["RAYON_NUM_THREADS"] = width
    if kind == "dseams":
        env["LD_LIBRARY_PATH"] = lib_env(lib)
        env.pop("OVITO_THREAD_COUNT", None)
    elif kind == "ovito":
        env["LD_LIBRARY_PATH"] = HOST_LIB
        env["OVITO_THREAD_COUNT"] = width
        env["QT_QPA_PLATFORM"] = "offscreen"
    else:
        env.pop("OVITO_THREAD_COUNT", None)
        env.pop("LD_LIBRARY_PATH", None)
    return env


def pinned(python, threads, kind, extra, lib):
    cores = "0" if threads == 1 else "0-7"
    cmd = ["taskset", "-c", cores, python, str(HERE),
           "--worker", kind, "--threads", str(threads), *extra]
    done = subprocess.run(
        cmd, cwd=ROOT, env=child_env(kind, threads, lib),
        text=True, capture_output=True)
    if done.returncode != 0:
        tail = (done.stderr or done.stdout or "").strip().splitlines()
        return {"error": tail[-1] if tail else f"exit {done.returncode}",
                "returncode": done.returncode, "code": kind}
    line = [row for row in done.stdout.splitlines() if row.startswith("{")][-1]
    payload = json.loads(line)
    if payload.get("cpus") != cores:
        return {"error": f"not pinned: {payload.get('cpus')}", "code": kind}
    return payload


def emit(payload):
    json.dump(payload, sys.stdout)
    sys.stdout.write("\n")


def cpus_now():
    try:
        with open("/proc/self/status", encoding="utf-8") as handle:
            for line in handle:
                if line.startswith("Cpus_allowed_list:"):
                    return line.split(":", 1)[1].strip()
    except OSError:
        return "unknown"
    return "unknown"


def run_worker(args):
    # OpenMP proc-bind narrows the main thread after the first parallel
    # region. Read the mask taskset set, before numpy or the library.
    pinned_cpus = cpus_now()
    sys.path.insert(0, str(HERE.parent))
    os.environ["OMP_NUM_THREADS"] = str(args.threads)
    os.environ["RAYON_NUM_THREADS"] = str(args.threads)
    os.environ["OMP_DYNAMIC"] = "FALSE"
    import numpy as np
    import h2h_protocol as proto

    data = np.load(args.frame)
    pos = data["positions"]
    side = float(data["side"])
    n = len(pos)

    if args.worker == "dseams":
        import ctypes
        lib = ctypes.CDLL(args.lib)
        lib.seams_lang_version.restype = ctypes.c_char_p
        lib.seams_last_error.restype = ctypes.c_char_p
        lib.seams_frame_from_arrays.restype = ctypes.c_void_p
        lib.seams_frame_from_arrays.argtypes = [
            ctypes.POINTER(ctypes.c_double), ctypes.c_int,
            ctypes.POINTER(ctypes.c_double),
            ctypes.POINTER(ctypes.c_int), ctypes.POINTER(ctypes.c_int)]
        lib.seams_frame_free.argtypes = [ctypes.c_void_p]
        lib.seams_frame_chill_plus.argtypes = [
            ctypes.c_void_p, ctypes.c_double, ctypes.c_int, ctypes.c_char_p,
            ctypes.POINTER(ctypes.c_int), ctypes.POINTER(ctypes.c_int)]
        lib.seams_frame_rings.argtypes = [
            ctypes.c_void_p, ctypes.c_double, ctypes.c_int, ctypes.c_char_p,
            ctypes.c_int, ctypes.POINTER(ctypes.c_int), ctypes.c_int]
        version = lib.seams_lang_version().decode()
        xyz = np.ascontiguousarray(pos, dtype=np.float64)
        box = np.array([0, 0, 0, side, side, side, 0, 0, 0], dtype=np.float64)
        types = np.ones(n, dtype=np.int32)
        ids = np.arange(1, n + 1, dtype=np.int32)
        handle = lib.seams_frame_from_arrays(
            xyz.ctypes.data_as(ctypes.POINTER(ctypes.c_double)), n,
            box.ctypes.data_as(ctypes.POINTER(ctypes.c_double)),
            types.ctypes.data_as(ctypes.POINTER(ctypes.c_int)),
            ids.ctypes.data_as(ctypes.POINTER(ctypes.c_int)))
        if not handle:
            raise RuntimeError(lib.seams_last_error().decode())

        def prep():
            return handle

        if args.operation == "chill+":
            def call(_frame):
                labels = np.empty(n, dtype=np.int32)
                counts = (ctypes.c_int * 9)()
                status = lib.seams_frame_chill_plus(
                    handle, proto.CUTOFF_A, 1, b"cutoff",
                    labels.ctypes.data_as(ctypes.POINTER(ctypes.c_int)), counts)
                if status != 0:
                    raise RuntimeError(lib.seams_last_error().decode())
                return [proto.DSEAMS_CLASS[int(v)] for v in labels]
        else:
            def call(_frame):
                counts = (ctypes.c_int * (proto.MAX_RING + 1))()
                status = lib.seams_frame_rings(
                    handle, proto.CUTOFF_A, 1, b"cutoff", proto.MAX_RING, counts,
                    proto.MAX_RING + 1)
                if status != 0:
                    raise RuntimeError(lib.seams_last_error().decode())
                return [int(counts[i]) for i in range(3, proto.MAX_RING + 1)]
        samples, result = proto.time_calls(prep, call)
        payload = {
            "code": "d-SEAMS", "version": version, "operation": args.operation,
            "molecules": int(n), "samples_ms": samples,
            "cpus": pinned_cpus, "threads": args.threads,
        }
        if args.operation == "chill+":
            payload["labels"] = result
        else:
            payload["census"] = result
        lib.seams_frame_free(handle)
        emit(payload)
        return

    if args.worker == "mdapy":
        import mdapy
        from mdapy import System
        from mdapy.identify import ChillPlus
        mdapy.set_num_threads(args.threads)

        def prep():
            return System(pos, box=side)

        if args.operation == "chill+":
            def call(system):
                labels = np.asarray(
                    ChillPlus(cutoff=proto.CUTOFF_A).compute(system).structure,
                    dtype=np.int32)
                return [proto.MDAPY_CLASS[int(v)] for v in labels]
        else:
            from mdapy.analysis import RingStatistics

            def call(system):
                counts = np.asarray(RingStatistics(
                    max_size=proto.MAX_RING, cutoff=proto.CUTOFF_A).compute(system).counts)
                return [int(counts[size]) if size < len(counts) else 0
                        for size in range(3, proto.MAX_RING + 1)]

        samples, result = proto.time_calls(prep, call)
        payload = {
            "code": "mdapy", "version": mdapy.__version__,
            "operation": args.operation, "molecules": int(n),
            "samples_ms": samples, "cpus": pinned_cpus,
            "threads": args.threads,
        }
        if args.operation == "chill+":
            payload["labels"] = result
        else:
            payload["census"] = result
        emit(payload)
        return

    if args.worker == "ovito":
        os.environ["OVITO_THREAD_COUNT"] = str(args.threads)
        os.environ["QT_QPA_PLATFORM"] = "offscreen"
        from ovito.data import DataCollection
        from ovito.modifiers import (
            ChillPlusModifier, CreateBondsModifier, FindRingsModifier)
        import ovito
        ovito_class = {
            0: "other", 1: "hexagonal", 2: "cubic", 3: "interfacial",
            4: "hydrate", 5: "interfacial_hydrate",
        }

        def prep():
            data = DataCollection()
            data.create_cell(
                matrix=[[side, 0, 0, 0], [0, side, 0, 0], [0, 0, side, 0]],
                pbc=(True, True, True))
            particles = data.create_particles(count=n)
            particles.create_property("Position", data=np.ascontiguousarray(pos))
            return data

        if args.operation == "chill+":
            def call(data):
                data.apply(ChillPlusModifier(cutoff=proto.CUTOFF_A))
                labels = np.asarray(data.particles["Structure Type"], dtype=np.int32)
                return [ovito_class[int(v)] for v in labels]
        else:
            def call(data):
                data.apply(CreateBondsModifier(cutoff=proto.CUTOFF_A))
                data.apply(FindRingsModifier(
                    minimum_ring_size=3, maximum_ring_size=proto.MAX_RING))
                hist = data.tables["ring-size-histogram"]
                found = dict(zip(
                    [int(v) for v in np.asarray(hist["Ring Size"])],
                    [int(v) for v in np.asarray(hist["Count"])]))
                return [found.get(size, 0) for size in range(3, proto.MAX_RING + 1)]
        samples, result = proto.time_calls(prep, call)
        payload = {
            "code": "OVITO", "version": ovito.version_string,
            "operation": args.operation, "molecules": int(n),
            "samples_ms": samples, "cpus": pinned_cpus,
            "threads": args.threads,
        }
        if args.operation == "chill+":
            payload["labels"] = result
        else:
            payload["census"] = result
        emit(payload)
        return

    if args.worker == "matscipy":
        import matscipy
        from ase import Atoms
        from matscipy.rings import ring_statistics
        atoms = Atoms(
            symbols=["O"] * n, positions=pos, cell=[side, side, side], pbc=True)

        def call(_atoms):
            counts = np.asarray(ring_statistics(
                atoms, proto.CUTOFF_A, maxlength=proto.MAX_RING + 1))
            return [int(counts[size]) if size < len(counts) else 0
                    for size in range(3, proto.MAX_RING + 1)]

        samples, census = proto.time_calls(lambda: atoms, call)
        emit({
            "code": "matscipy", "version": matscipy.__version__,
            "operation": "rings", "molecules": int(n), "samples_ms": samples,
            "census": census, "cpus": pinned_cpus, "threads": args.threads,
        })
        return

    import importlib.metadata
    import networkx as nx
    import pairlist
    from cycless.cycles import cycles_iter
    version = importlib.metadata.version("cycless")

    def call(points):
        fractional = points / side
        graph = nx.Graph()
        graph.add_nodes_from(range(n))
        graph.add_edges_from(pairlist.pairs_iter(
            fractional, proto.CUTOFF_A, np.diag([side, side, side]),
            distance=False))
        census = [0, 0, 0, 0]
        for ring in cycles_iter(graph, proto.MAX_RING, pos=fractional):
            size = len(ring)
            if 3 <= size <= proto.MAX_RING:
                census[size - 3] += 1
        return census

    samples, census = proto.time_calls(lambda: pos, call)
    emit({
        "code": "cycless", "version": version, "operation": "rings",
        "molecules": int(n), "samples_ms": samples, "census": census,
        "cpus": pinned_cpus, "threads": args.threads,
    })


def agreement(reference, other, operation):
    if "error" in other or "error" in reference:
        return "", other.get("error") or reference.get("error")
    if operation == "chill+":
        left, right = reference["labels"], other["labels"]
        matched = sum(a == b for a, b in zip(left, right))
        return f"{matched / len(left):.6f}", f"matched={matched}/{len(left)}"
    same = reference["census"] == other["census"]
    return f"{float(same):.6f}", (
        f"dseams={reference['census']};other={other['census']}")


def write_row(writer, payload, operation, repetitions, fraction, detail):
    if "error" in payload:
        writer.writerow({
            "threads": payload.get("threads", ""),
            "code": payload.get("code", "unknown"),
            "version": payload.get("version", ""),
            "operation": operation,
            "molecules": payload.get("molecules", ""),
            "repetitions": repetitions,
            "warmup": "",
            "n": 0,
            "median_ms": "",
            "q1_ms": "",
            "q3_ms": "",
            "iqr_ms": "",
            "agreement": "",
            "detail": payload["error"],
            "samples_ms": "",
            "cpus": "",
        })
        return
    sys.path.insert(0, str(HERE.parent))
    import h2h_protocol as proto
    summary = proto.summarize(payload["samples_ms"])
    writer.writerow({
        "threads": payload["threads"],
        "code": payload["code"],
        "version": payload["version"],
        "operation": operation,
        "molecules": payload["molecules"],
        "repetitions": repetitions,
        "warmup": proto.WARMUP,
        "n": summary["n"],
        "median_ms": f"{summary['median_ms']:.3f}",
        "q1_ms": f"{summary['q1_ms']:.3f}",
        "q3_ms": f"{summary['q3_ms']:.3f}",
        "iqr_ms": f"{summary['iqr_ms']:.3f}",
        "agreement": fraction if payload["code"] != "d-SEAMS" else "1.000000",
        "detail": detail,
        "samples_ms": " ".join(f"{v:.3f}" for v in payload["samples_ms"]),
        "cpus": payload["cpus"],
    })


def campaign(args):
    sys.path.insert(0, str(HERE.parent))
    os.environ.pop("LD_LIBRARY_PATH", None)
    import numpy as np
    import h2h_protocol as proto
    import release_build

    lib = args.lib
    if not lib:
        lib, _build = release_build.ensure_release()
    python = args.python or os.environ.get("H2H_PYTHON") or sys.executable
    out = ROOT / args.out
    out.parent.mkdir(parents=True, exist_ok=True)
    sizes = tuple(int(v) for v in args.sizes.split(",") if v)
    jobs = [
        ("chill+", "mdapy"),
        ("chill+", "ovito"),
        ("rings", "mdapy"),
        ("rings", "ovito"),
        ("rings", "matscipy"),
    ]
    if not args.skip_cycless:
        jobs.append(("rings", "cycless"))
    if args.only:
        jobs = [job for job in jobs if job[1] == args.only]
    with out.open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=FIELDS)
        writer.writeheader()
        for threads in (1, 8):
            for reps in sizes:
                pos, side = proto.ice_ic(reps)
                frame = ROOT / "benchmarks" / f".ice-ic-{len(pos)}.npz"
                np.savez(frame, positions=pos, side=side)
                by_op = {}
                for operation, other in jobs:
                    if operation not in by_op:
                        by_op[operation] = pinned(
                            python, threads, "dseams",
                            ["--operation", operation, "--frame", str(frame),
                             "--lib", lib], lib)
                    reference = by_op[operation]
                    rival = pinned(
                        python, threads, other,
                        ["--operation", operation, "--frame", str(frame),
                         "--lib", lib], lib)
                    fraction, detail = agreement(reference, rival, operation)
                    reference["molecules"] = len(pos)
                    rival["molecules"] = len(pos)
                    reference["threads"] = threads
                    rival["threads"] = threads
                    if other == jobs[0][1] or operation != jobs[0][0]:
                        # d-SEAMS is written once per operation. The first
                        # chill+ rival and every later operation's first rival
                        # are the points where the reference is new.
                        first = [name for op, name in jobs if op == operation][0]
                        if other == first:
                            write_row(writer, reference, operation, reps,
                                      "1.000000", detail)
                    write_row(writer, rival, operation, reps, fraction, detail)
                    handle.flush()
                    print(
                        f"threads={threads} N={len(pos)} {operation} "
                        f"d-SEAMS={reference.get('samples_ms', ['err'])} "
                        f"{other}={rival.get('error') or rival.get('samples_ms')}",
                        flush=True)
                frame.unlink(missing_ok=True)


def main(argv=None):
    args = parse_args(sys.argv[1:] if argv is None else argv)
    if args.worker:
        run_worker(args)
        return
    campaign(args)


if __name__ == "__main__":
    main()
