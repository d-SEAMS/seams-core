#!/usr/bin/env python3
"""Neighbours and q6 with every OpenMP placement variable unset.

One fresh process per surface and operation. Three warm-up calls, then
one timed call. The only OpenMP variable set in the child is
OMP_NUM_THREADS. OMP_PROC_BIND, OMP_PLACES, OMP_WAIT_POLICY and
GOMP_SPINCOUNT are removed, as is SEAMS_FAIR_WIDTH. n is at least 8.
The header records the tree of this checkout and the core count. It
does not record a commit hash or a host name.
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
sys.path.insert(0, str(ROOT / "scripts"))

import paper_stats as stats  # noqa: E402

CUTOFF = 3.5
WARMUP = 3
MIN_N = 8
CELLS = 24
DROP = (
    "OMP_PROC_BIND",
    "OMP_PLACES",
    "OMP_WAIT_POLICY",
    "GOMP_SPINCOUNT",
    "OMP_DISPLAY_ENV",
    "OMP_SCHEDULE",
    "SEAMS_FAIR_WIDTH",
    "SEAMS_NEIGH_PROFILE",
    "SEAMS_QL_TIME",
)


def ice():
    import numpy as np
    import h2h_protocol as proto
    xyz, side = proto.ice_ic(CELLS)
    xyz = np.ascontiguousarray(xyz, dtype=np.float64)
    box = np.array([0.0, 0.0, 0.0, side, side, side, 0.0, 0.0, 0.0], dtype=np.float64)
    return xyz, box


def c_call(lib_path: str, op: str):
    import ctypes
    import numpy as np
    lib = ctypes.CDLL(lib_path)
    xyz, box = ice()
    n = len(xyz)
    dp = ctypes.POINTER(ctypes.c_double)
    ip = ctypes.POINTER(ctypes.c_int)
    xp = xyz.ctypes.data_as(dp)
    bp = box.ctypes.data_as(dp)
    if op == "neighbours":
        lib.seams_cutoff_csr.argtypes = [
            dp, ctypes.c_int, dp, ctypes.c_double, ip, ip, ctypes.c_int, ip]
        offsets = np.empty(n + 1, dtype=np.int32)
        cols = np.empty(16 * n, dtype=np.int32)
        count = ctypes.c_int(0)

        def run():
            rc = lib.seams_cutoff_csr(
                xp, n, bp, CUTOFF, offsets.ctypes.data_as(ip),
                cols.ctypes.data_as(ip), len(cols), ctypes.byref(count))
            if rc != 0:
                raise RuntimeError(f"seams_cutoff_csr returned {rc}")
            return f"pairs={count.value}"
        return run
    lib.seams_ql_cutoff.argtypes = [
        dp, ctypes.c_int, dp, ctypes.c_double, ctypes.c_int, dp]
    ql = np.empty(n, dtype=np.float64)

    def run():
        rc = lib.seams_ql_cutoff(xp, n, bp, CUTOFF, 6, ql.ctypes.data_as(dp))
        if rc != 0:
            raise RuntimeError(f"seams_ql_cutoff returned {rc}")
        return f"mean={float(ql.mean()):.12f}"
    return run


def frame_call(lang_path: str, op: str):
    import ctypes
    import numpy as np
    lib = ctypes.CDLL(lang_path)
    xyz, box = ice()
    n = len(xyz)
    dp = ctypes.POINTER(ctypes.c_double)
    ip = ctypes.POINTER(ctypes.c_int)
    lib.seams_frame_from_arrays.restype = ctypes.c_void_p
    lib.seams_frame_from_arrays.argtypes = [dp, ctypes.c_int, dp, ip, ip]
    types = np.ones(n, dtype=np.int32)
    ids = np.arange(1, n + 1, dtype=np.int32)
    frame = lib.seams_frame_from_arrays(
        xyz.ctypes.data_as(dp), n, box.ctypes.data_as(dp),
        types.ctypes.data_as(ip), ids.ctypes.data_as(ip))
    if not frame:
        raise RuntimeError("could not build the frame")
    if op == "neighbours":
        lib.seams_frame_neighbours.restype = ctypes.c_int
        lib.seams_frame_neighbours.argtypes = [
            ctypes.c_void_p, ctypes.c_double, ctypes.c_int, ctypes.c_char_p]

        def run():
            rc = lib.seams_frame_neighbours(frame, CUTOFF, 1, b"cutoff")
            if rc != 0:
                raise RuntimeError(f"seams_frame_neighbours returned {rc}")
            return "frame=neighbours"
        return run
    lib.seams_frame_steinhardt.restype = ctypes.c_int
    lib.seams_frame_steinhardt.argtypes = [
        ctypes.c_void_p, ctypes.c_double, ctypes.c_int, ctypes.c_int, dp, dp]
    ql = np.empty(n, dtype=np.float64)

    def run():
        rc = lib.seams_frame_steinhardt(
            frame, CUTOFF, 1, 6, ql.ctypes.data_as(dp), None)
        if rc != 0:
            raise RuntimeError(f"seams_frame_steinhardt returned {rc}")
        return f"mean={float(ql.mean()):.12f}"
    return run


def worker(args) -> int:
    import time
    cores, affinity = stats.cores_now()
    if args.surface == "c":
        call = c_call(args.lib, args.op)
    else:
        call = frame_call(args.lang, args.op)
    timed = None
    check = ""
    for index in range(WARMUP + 1):
        t0 = time.perf_counter()
        check = call()
        elapsed = (time.perf_counter() - t0) * 1000.0
        if index == WARMUP:
            timed = elapsed
    print(json.dumps({
        "ms": timed, "check": check, "cores": cores, "affinity": affinity,
    }), flush=True)
    return 0


def child_env(lib: str, lang: str) -> dict[str, str]:
    env = dict(os.environ)
    for key in DROP:
        env.pop(key, None)
    env["OMP_NUM_THREADS"] = "8"
    env["OPENBLAS_NUM_THREADS"] = "1"
    env["MKL_NUM_THREADS"] = "1"
    env["NUMEXPR_NUM_THREADS"] = "1"
    paths = [str(Path(lib).resolve().parent), str(Path(lang).resolve().parent)]
    if env.get("LD_LIBRARY_PATH"):
        paths.append(env["LD_LIBRARY_PATH"])
    env["LD_LIBRARY_PATH"] = os.pathsep.join(paths)
    return env


def drive(surface: str, op: str, lib: str, lang: str) -> dict:
    cmd = [sys.executable, __file__, "--worker", "--surface", surface, "--op", op,
           "--lib", lib, "--lang", lang]
    done = subprocess.run(cmd, env=child_env(lib, lang), capture_output=True, text=True,
                          check=False)
    if done.returncode != 0:
        raise RuntimeError(f"{surface} {op}: {done.stderr.strip()[-2000:]}")
    lines = [line for line in done.stdout.splitlines() if line.startswith("{")]
    if not lines:
        raise RuntimeError(f"{surface} {op}: no result")
    return json.loads(lines[-1])


def display_env(lib: str, lang: str) -> str:
    env = child_env(lib, lang)
    env["OMP_DISPLAY_ENV"] = "VERBOSE"
    code = "import ctypes,sys; ctypes.CDLL(sys.argv[1])"
    done = subprocess.run([sys.executable, "-c", code, lib], env=env,
                          capture_output=True, text=True, check=False)
    text = done.stderr.strip()
    if "OPENMP DISPLAY ENVIRONMENT" not in text:
        text = (done.stderr + done.stdout).strip()
    return text


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--out", type=Path, default=ROOT / "bench/omp-unbound")
    parser.add_argument("--lib", default=str(ROOT / "bbdir-release/src/libyodaLib.so"))
    parser.add_argument("--lang", default=str(ROOT / "bbdir-release/bindings/libseams_lang.so"))
    parser.add_argument("--n", type=int, default=MIN_N)
    parser.add_argument("--worker", action="store_true")
    parser.add_argument("--surface", choices=("c", "frame"), default="c")
    parser.add_argument("--op", choices=("neighbours", "q6"), default="neighbours")
    args = parser.parse_args()
    if args.worker:
        return worker(args)
    if args.n < MIN_N:
        parser.error(f"--n must be at least {MIN_N}")
    if not Path(args.lib).is_file() or not Path(args.lang).is_file():
        parser.error("release libraries are not built")

    jobs = [("c", "neighbours"), ("c", "q6"), ("frame", "neighbours"), ("frame", "q6")]
    rows = []
    raw = [
        "# unbound omp sitting",
        f"# source_tree {stats.source_tree(ROOT)}",
        "# molecules 110592 cutoff 3.5 ice Ic",
        f"# warmup {WARMUP} timed 1 processes {args.n} one process per sample",
        "# OMP_NUM_THREADS=8",
        "# OMP_PROC_BIND OMP_PLACES OMP_WAIT_POLICY GOMP_SPINCOUNT unset",
        f"# cores {len(os.sched_getaffinity(0))}",
    ]
    shown = display_env(args.lib, args.lang)
    if shown:
        raw.extend(f"# {line}" for line in shown.splitlines())
    for surface, op in jobs:
        samples = []
        check = ""
        cores = ""
        affinity = ""
        for _ in range(args.n):
            got = drive(surface, op, args.lib, args.lang)
            samples.append(got["ms"])
            check = got["check"]
            cores = str(got["cores"])
            affinity = got["affinity"]
        summary = stats.summarize(samples)
        row = {
            "surface": surface,
            "operation": op,
            "threads": 8,
            "cores": cores,
            "affinity": affinity,
            "n": summary["n"],
            "median_ms": f"{summary['median_ms']:.3f}",
            "q1_ms": f"{summary['q1_ms']:.3f}",
            "q3_ms": f"{summary['q3_ms']:.3f}",
            "check": check,
        }
        rows.append(row)
        raw.append(
            f"{surface} {op} threads=8 cores={cores} affinity={affinity} {check} "
            f"samples_ms {stats.fmt_samples(samples)}")
        print(f"{surface:6} {op:12} {row['median_ms']:>8} ms", flush=True)

    args.out.mkdir(parents=True, exist_ok=True)
    with open(args.out / "unbound.csv", "w", newline="") as handle:
        handle.write(f"# source_tree {stats.source_tree(ROOT)}\n")
        writer = csv.DictWriter(handle, fieldnames=list(rows[0].keys()), lineterminator="\n")
        writer.writeheader()
        writer.writerows(rows)
    (args.out / "raw-unbound.txt").write_text("\n".join(raw) + "\n")
    return 0


if __name__ == "__main__":
    sys.exit(main())
