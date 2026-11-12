#!/usr/bin/env python3
"""Cutoff CHILL+ against ClavisCage 0.3.18 on ideal ice Ic.

One warm-up, six timed calls, median and type-7 interquartile range.
One core is taskset -c 0. Eight cores are taskset -c 0-7. The ice is
built in the worker and not written out. Hydrogens sit on the oxygen,
which is the placement ClavisCage's own tests use; its CHILL+ reads
only the oxygen coordinates. ClavisCage's public call parses a JSON
water list, and that parse is inside the timed call.
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

FIELDS = [
    "threads", "code", "version", "operation", "molecules", "repetitions",
    "warmup", "n", "median_ms", "q1_ms", "q3_ms", "iqr_ms", "agreement",
    "detail", "samples_ms", "cpus",
]

CLAVIS_CLASS = {
    0: "other",
    1: "hexagonal",
    2: "cubic",
    3: "interfacial",
    4: "hydrate",
    5: "interfacial_hydrate",
}


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
    parser.add_argument("--out", default="benchmarks/h2h-watch-v11.csv")
    parser.add_argument("--python")
    parser.add_argument("--threads", type=int, default=1)
    parser.add_argument("--worker", choices=("dseams", "clavis"))
    parser.add_argument("--repetitions", type=int)
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


def waters_json(points_nm):
    """One JSON array. Hydrogens coincide with the oxygen."""
    waters = []
    for oxygen in points_nm:
        xyz = [float(oxygen[0]), float(oxygen[1]), float(oxygen[2])]
        waters.append({
            "oxygen": xyz,
            "hydrogen1": xyz,
            "hydrogen2": xyz,
        })
    return json.dumps(waters)


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
            labels = np.empty(n, dtype=np.int32)
            counts = (ctypes.c_int * 9)()
            status = lib.seams_frame_chill_plus(
                handle, proto.CUTOFF_A, 1, b"cutoff",
                labels.ctypes.data_as(ctypes.POINTER(ctypes.c_int)), counts)
            if status != 0:
                raise RuntimeError(lib.seams_last_error().decode())
            return [proto.DSEAMS_CLASS[int(v)] for v in labels]

        samples, result = proto.time_calls(lambda: handle, call)
        lib.seams_frame_free(handle)
        emit({
            "code": "d-SEAMS", "version": version, "operation": "chill+",
            "molecules": int(n), "samples_ms": samples, "labels": result,
            "cpus": pinned_cpus, "threads": args.threads,
        })
        return

    import clavis_engine as engine
    points_nm = np.ascontiguousarray(pos, dtype=np.float64) / 10.0
    payload = waters_json(points_nm)
    side_nm = side / 10.0
    box_nm = [side_nm, side_nm, side_nm]
    empty = "[]"

    def call(text):
        raw = engine.analyze_waters(text, empty, box_nm, "chill", 0.35)
        parsed = json.loads(raw)
        return [CLAVIS_CLASS[int(v)] for v in parsed["per_atom_class"]]

    samples, result = proto.time_calls(lambda: payload, call)
    emit({
        "code": "ClavisCage", "version": engine.version(),
        "operation": "chill+", "molecules": int(n), "samples_ms": samples,
        "labels": result, "cpus": pinned_cpus, "threads": args.threads,
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


def class_counts(labels):
    counts = {}
    for name in labels:
        counts[name] = counts.get(name, 0) + 1
    return counts


def agreement(reference, other):
    if "error" in other or "error" in reference:
        return "", other.get("error") or reference.get("error")
    left, right = reference["labels"], other["labels"]
    matched = sum(a == b for a, b in zip(left, right))
    detail = (
        f"matched={matched}/{len(left)};"
        f"dseams={class_counts(left)};clavis={class_counts(right)}"
    )
    return f"{matched / len(left):.6f}", detail


def write_row(writer, payload, repetitions, fraction, detail):
    if "error" in payload:
        writer.writerow({
            "threads": payload.get("threads", ""),
            "code": payload.get("code", "unknown"),
            "version": "",
            "operation": "chill+",
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
        "operation": "chill+",
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
                extra = ["--repetitions", str(reps), "--lib", lib]
                reference = pinned(python, threads, "dseams", extra, lib)
                rival = pinned(python, threads, "clavis", extra, lib)
                fraction, detail = agreement(reference, rival)
                reference["molecules"] = proto.molecules_for(reps)
                reference["threads"] = threads
                rival["threads"] = threads
                if "molecules" not in rival:
                    rival["molecules"] = proto.molecules_for(reps)
                write_row(writer, reference, reps, "1.000000", detail)
                write_row(writer, rival, reps, fraction, detail)
                handle.flush()
                print(
                    f"threads={threads} N={proto.molecules_for(reps)} "
                    f"d-SEAMS={reference.get('samples_ms', ['err'])} "
                    f"clavis={rival.get('error') or rival.get('samples_ms')} "
                    f"agreement={fraction}",
                    flush=True)


def main(argv=None):
    args = parse_args(sys.argv[1:] if argv is None else argv)
    if args.worker:
        run_worker(args)
        return
    campaign(args)


if __name__ == "__main__":
    main()
