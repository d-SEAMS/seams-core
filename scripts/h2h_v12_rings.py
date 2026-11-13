#!/usr/bin/env python3
"""Primitive rings against rings-statistics-matter 0.2.0 and 0.3.0 on ideal ice Ic.

One warm-up, six timed calls, median and type-7 interquartile range.
One core is taskset -c 0. Eight cores are taskset -c 0-7. The ice is
built in the worker and not written out.

The rival is the public ring_statistics call. Its Rust kernel tries
midpoints closer than maxlvl/2, so a six-ring (midpoint at distance 3)
needs maxlvl 8. The package docstring's "maxlvl - 1" would pick 7, and
that depth returns no rings on this crystal. supercell is off: the
periodic cell already has two six-rings per molecule. The default
supercell repeats only the 1000-molecule cell, and that call is a
separate operation.
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
MAXLVL = 8

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
    parser.add_argument("--out", default="benchmarks/h2h-watch-v12.csv")
    parser.add_argument("--python")
    parser.add_argument("--threads", type=int, default=1)
    parser.add_argument("--worker", choices=("dseams", "rsm"))
    parser.add_argument("--repetitions", type=int)
    parser.add_argument("--operation", choices=("rings", "rings-default"), default="rings")
    parser.add_argument("--sizes", default="5,11,24")
    return parser.parse_args(argv)


def child_env(kind, threads, lib):
    env = os.environ.copy()
    width = str(threads)
    env["OMP_NUM_THREADS"] = width
    env["OMP_DYNAMIC"] = "FALSE"
    env["MKL_NUM_THREADS"] = "1"
    env["OPENBLAS_NUM_THREADS"] = "1"
    env["NUMEXPR_NUM_THREADS"] = "1"
    env["RAYON_NUM_THREADS"] = width
    env.pop("OMP_PROC_BIND", None)
    if kind == "dseams":
        env["LD_LIBRARY_PATH"] = lib_env(lib)
    else:
        env.pop("LD_LIBRARY_PATH", None)
    return env


def cpus_now():
    try:
        with open("/proc/self/status", encoding="utf-8") as handle:
            for line in handle:
                if line.startswith("Cpus_allowed_list:"):
                    return line.split(":", 1)[1].strip()
    except OSError:
        return "unknown"
    return "unknown"


def emit(payload):
    json.dump(payload, sys.stdout)
    sys.stdout.write("\n")


def census_from_counts(counts):
    """Ring counts of size 3 through 6. counts[n - 1] is the n-ring count."""
    import numpy as np
    values = np.asarray(counts, dtype=np.float64)
    if values.shape[0] < 6:
        raise RuntimeError(f"ring count vector is shorter than 6: {values.shape}")
    chosen = values[2:6]
    if not np.allclose(chosen, np.round(chosen)):
        raise RuntimeError(f"non-integer ring counts: {chosen.tolist()}")
    return [int(round(float(v))) for v in chosen]


def run_worker(args):
    pinned_cpus = cpus_now()
    sys.path.insert(0, str(HERE.parent))
    os.environ["OMP_NUM_THREADS"] = str(args.threads)
    os.environ["RAYON_NUM_THREADS"] = str(args.threads)
    os.environ["OMP_DYNAMIC"] = "FALSE"
    import numpy as np
    import h2h_protocol as proto

    pos, side = proto.ice_ic(args.repetitions)
    n = len(pos)
    supercell = args.operation == "rings-default"

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

        def call(_frame):
            counts = (ctypes.c_int * (proto.MAX_RING + 1))()
            status = lib.seams_frame_rings(
                handle, proto.CUTOFF_A, 1, b"cutoff", proto.MAX_RING, counts,
                proto.MAX_RING + 1)
            if status != 0:
                raise RuntimeError(lib.seams_last_error().decode())
            return [int(counts[i]) for i in range(3, proto.MAX_RING + 1)]

        samples, result = proto.time_calls(lambda: handle, call)
        lib.seams_frame_free(handle)
        emit({
            "code": "d-SEAMS", "version": version, "operation": args.operation,
            "molecules": int(n), "samples_ms": samples, "census": result,
            "cpus": pinned_cpus, "threads": args.threads,
        })
        return

    from ase import Atoms
    from importlib.metadata import version as pkg_version
    from rings_statistics_matter import ring_statistics

    atoms = Atoms(
        symbols=["O"] * n, positions=np.ascontiguousarray(pos, dtype=np.float64),
        cell=[side, side, side], pbc=True)

    def call(frame):
        counts, _rings = ring_statistics(
            frame, cutoff=proto.CUTOFF_A, maxlvl=MAXLVL, supercell=supercell,
            deduplicate=True, backend="rust", num_threads=args.threads)
        return census_from_counts(counts)

    samples, result = proto.time_calls(lambda: atoms, call)
    emit({
        "code": "rings-statistics-matter",
        "version": pkg_version("rings-statistics-matter"),
        "operation": args.operation, "molecules": int(n),
        "samples_ms": samples, "census": result,
        "cpus": pinned_cpus, "threads": args.threads,
        "maxlvl": MAXLVL, "supercell": supercell,
    })


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
                "code": kind, "threads": threads}
    line = [row for row in done.stdout.splitlines() if row.startswith("{")][-1]
    payload = json.loads(line)
    if payload.get("cpus") != cores:
        return {"error": f"not pinned: {payload.get('cpus')}",
                "code": kind, "threads": threads}
    return payload


def agreement(reference, other):
    if "error" in other or "error" in reference:
        return "", other.get("error") or reference.get("error")
    left, right = reference["census"], other["census"]
    same = left == right
    detail = f"dseams={left};rings-statistics-matter={right};maxlvl={MAXLVL}"
    if other.get("supercell"):
        detail += ";supercell=true"
    return ("1.000000" if same else "0.000000"), detail


def write_row(writer, payload, repetitions, fraction, detail):
    if "error" in payload:
        writer.writerow({
            "threads": payload.get("threads", ""),
            "code": payload.get("code", "unknown"),
            "version": "",
            "operation": payload.get("operation", "rings"),
            "molecules": "",
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
        "operation": payload["operation"],
        "molecules": payload["molecules"],
        "repetitions": repetitions,
        "warmup": proto.WARMUP,
        "n": summary["n"],
        "median_ms": f"{summary['median_ms']:.3f}",
        "q1_ms": f"{summary['q1_ms']:.3f}",
        "q3_ms": f"{summary['q3_ms']:.3f}",
        "iqr_ms": f"{summary['iqr_ms']:.3f}",
        "agreement": fraction,
        "detail": detail,
        "samples_ms": " ".join(f"{v:.3f}" for v in payload["samples_ms"]),
        "cpus": payload["cpus"],
    })


def campaign(args):
    sys.path.insert(0, str(HERE.parent))
    os.environ.pop("LD_LIBRARY_PATH", None)
    import h2h_protocol as proto
    import release_build

    lib = args.lib
    if not lib:
        lib, _build = release_build.ensure_release()
    else:
        release_build.ensure_release(lib)
    python = args.python or os.environ.get("H2H_PYTHON") or sys.executable
    out = ROOT / args.out
    out.parent.mkdir(parents=True, exist_ok=True)
    sizes = tuple(int(v) for v in args.sizes.split(",") if v)
    with out.open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=FIELDS)
        writer.writeheader()
        for threads in (1, 8):
            for reps in sizes:
                operations = ["rings"]
                if reps == 5:
                    operations.append("rings-default")
                for operation in operations:
                    extra = [
                        "--repetitions", str(reps), "--lib", lib,
                        "--operation", operation,
                    ]
                    reference = pinned(python, threads, "dseams", extra, lib)
                    rival = pinned(python, threads, "rsm", extra, lib)
                    fraction, detail = agreement(reference, rival)
                    molecules = proto.molecules_for(reps)
                    for payload in (reference, rival):
                        payload["molecules"] = molecules
                        payload["threads"] = threads
                        payload["operation"] = operation
                    write_row(writer, reference, reps, fraction, detail)
                    write_row(writer, rival, reps, fraction, detail)
                    handle.flush()
                    print(
                        f"threads={threads} N={molecules} op={operation} "
                        f"d-SEAMS={reference.get('samples_ms', ['err'])} "
                        f"rsm={rival.get('error') or rival.get('samples_ms')} "
                        f"agreement={fraction} {detail}",
                        flush=True)


def main(argv=None):
    args = parse_args(sys.argv[1:] if argv is None else argv)
    if args.worker:
        run_worker(args)
        return
    campaign(args)


if __name__ == "__main__":
    main()
