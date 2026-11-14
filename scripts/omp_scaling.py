#!/usr/bin/env python3
"""OpenMP scaling of the cutoff neighbour list, q4, q6 and the ring census.

Ideal ice Ic, 110592 molecules, cutoff 3.5 A, rings up to six. Every
timed call builds its own neighbour list. d-SEAMS is called through the
C API (seams_cutoff_csr, seams_ql_cutoff, seams_primitive_rings) at 1,
2, 4 and 8 threads; mdapy and pyscal3 q4 and q6 run at the rival widths
(8 by default) in the same sitting.

One process per code, operation and width, pinned with taskset to cpus
0 .. width-1. One warm-up call, then n timed calls (n >= 8). The median
of an even sample is the mean of the two central values and the IQR is
type 7. Each row records the affinity of the worker. No host name is
read or written.

--before times a second libyodaLib.so, for example a build of the parent
commit, under the label given by --before-label. The library path is
not written to the outputs.
"""

from __future__ import annotations

import argparse
import csv
import json
import os
import subprocess
import sys
import time
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "scripts"))

import paper_stats as stats  # noqa: E402

CUTOFF = 3.5
REPETITIONS = 24
MAX_RING = 6
WARMUP = 1
WIDTHS = (1, 2, 4, 8)
DSEAMS_OPS = ("neighbours", "q4", "q6", "rings")
RIVAL_OPS = ("q4", "q6")


def frame():
    import h2h_protocol as proto
    return proto.ice_ic(REPETITIONS)


def hold_affinity(threads: int) -> int:
    """Put the process back on its taskset cpus after a runtime narrows it."""
    try:
        os.sched_setaffinity(0, set(range(threads)))
        return len(os.sched_getaffinity(0))
    except (AttributeError, OSError):
        return threads


def dseams_call(lib_path: str, op: str):
    import ctypes
    import numpy as np
    xyz, side = frame()
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


def mdapy_call(op: str, threads: int):
    os.environ["MDAPY_NUM_THREADS"] = str(threads)
    import mdapy
    from mdapy import System, set_num_threads
    from mdapy.analysis import SteinhardtParameters
    xyz, side = frame()
    set_num_threads(threads)
    system = System(xyz, box=side)
    tool = SteinhardtParameters(degrees=(int(op[1]),), cutoff=CUTOFF)

    def run():
        hold_affinity(threads)
        set_num_threads(threads)
        result = tool.compute(system)
        return f"rows={int(result.q.shape[0])}"
    return run, f"mdapy {mdapy.__version__}"


def pyscal_call(op: str, threads: int):
    import numpy as np
    import pyscal3
    from ase import Atoms
    os.environ["PYSCAL_NUM_THREADS"] = str(threads)
    pyscal3.set_num_threads(threads)
    xyz, side = frame()
    atoms = Atoms(symbols=["O"] * len(xyz), positions=np.ascontiguousarray(xyz),
                  cell=[side, side, side], pbc=True)
    order = int(op[1])

    def run():
        hold_affinity(threads)
        pyscal3.set_num_threads(threads)
        pyscal3.find_neighbors(atoms, method="cutoff", cutoff=CUTOFF)
        values = np.asarray(pyscal3.steinhardt_parameter(atoms, order)[0]).reshape(-1)
        return f"rows={int(values.shape[0])}"
    return run, f"pyscal3 {pyscal3.__version__}"


def worker(args) -> int:
    # Read before a runtime loads: OMP_PROC_BIND narrows the main thread.
    cores, affinity = stats.cores_now()
    threads = args.threads
    rival = args.code != "dseams"
    if args.code == "dseams":
        call = dseams_call(args.lib, args.op)
        version = args.label
    elif args.code == "mdapy":
        call, version = mdapy_call(args.op, threads)
    elif args.code == "pyscal3":
        call, version = pyscal_call(args.op, threads)
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


def worker_env(code: str, threads: int, lib: str) -> dict[str, str]:
    if code == "dseams":
        env = stats.pin_env(threads)
    else:
        # OMP_PROC_BIND would pin the rival's main thread to one cpu, and
        # its own thread pool would inherit that mask.
        env = dict(os.environ)
        env.pop("OMP_PROC_BIND", None)
        env.pop("OMP_PLACES", None)
        width = str(threads)
        env.update({"OMP_NUM_THREADS": width, "OMP_DYNAMIC": "FALSE",
                    "RAYON_NUM_THREADS": width, "MDAPY_NUM_THREADS": width,
                    "PYSCAL_NUM_THREADS": width, "OPENBLAS_NUM_THREADS": "1",
                    "MKL_NUM_THREADS": "1", "NUMEXPR_NUM_THREADS": "1"})
    extra = [str(Path(lib).resolve().parent),
             str(ROOT / "bbdir-release/subprojects/readcon-core")]
    if env.get("LD_LIBRARY_PATH"):
        extra.append(env["LD_LIBRARY_PATH"])
    env["LD_LIBRARY_PATH"] = os.pathsep.join(extra)
    return env


def drive(code: str, op: str, threads: int, lib: str, label: str, repeats: int) -> dict:
    cmd = ["taskset", "-c", stats.taskset_cores(threads), sys.executable, __file__,
           "--worker", "--code", code, "--op", op, "--threads", str(threads),
           "--repeats", str(repeats), "--lib", lib, "--label", label]
    done = subprocess.run(cmd, env=worker_env(code, threads, lib), capture_output=True,
                          text=True, check=False)
    if done.returncode != 0:
        raise RuntimeError(f"{code} {op} {threads}: {done.stderr.strip()[-2000:]}")
    lines = [line for line in done.stdout.splitlines() if line.startswith("{")]
    if not lines:
        raise RuntimeError(f"{code} {op} {threads}: no result line")
    return json.loads(lines[-1])


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__,
                                     formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("--out", type=Path, default=ROOT / "bench/omp-scaling")
    parser.add_argument("--lib", default=str(ROOT / "bbdir-release/src/libyodaLib.so"))
    parser.add_argument("--label", default="after")
    parser.add_argument("--before", default="", help="a second libyodaLib.so to time")
    parser.add_argument("--before-label", default="before")
    parser.add_argument("--before-commit", default="", help="commit the --before build is from")
    parser.add_argument("--widths", default=",".join(str(w) for w in WIDTHS))
    parser.add_argument("--rival-widths", default="8")
    parser.add_argument("--rivals", default="mdapy,pyscal3")
    parser.add_argument("--repeats", type=int, default=8)
    parser.add_argument("--worker", action="store_true")
    parser.add_argument("--code", default="dseams")
    parser.add_argument("--op", default="q4")
    parser.add_argument("--threads", type=int, default=1)
    args = parser.parse_args()
    if args.worker:
        return worker(args)
    if args.repeats < 8:
        parser.error("--repeats must be at least 8")

    widths = [int(w) for w in args.widths.split(",") if w]
    rival_widths = [int(w) for w in args.rival_widths.split(",") if w]
    rivals = [r for r in args.rivals.split(",") if r]
    available = len(os.sched_getaffinity(0))
    if max(widths + rival_widths) > available:
        parser.error(f"widths need {max(widths + rival_widths)} cpus, {available} allowed")

    builds = []
    if args.before:
        builds.append((args.before_label, args.before))
    builds.append((args.label, args.lib))
    jobs = []
    for op in DSEAMS_OPS:
        for threads in widths:
            for label, lib in builds:
                jobs.append(("dseams", op, threads, label, lib))
    for code in rivals:
        for op in RIVAL_OPS:
            for threads in rival_widths:
                jobs.append((code, op, threads, code, args.lib))

    args.out.mkdir(parents=True, exist_ok=True)
    header = [
        "# omp-scaling sitting",
        f"# frame ice Ic molecules={8 * REPETITIONS ** 3} cutoff={CUTOFF} max_ring={MAX_RING}",
        f"# warmup={WARMUP} repeats={args.repeats} taskset=0..width-1 one process per row",
    ] + [f"# {line}" for line in stats.machine_lines(ROOT)]
    if args.before:
        header.append(f"# build {args.before_label} commit {args.before_commit or 'unrecorded'}")
    header.append(f"# build {args.label} commit {stats.source_commit(ROOT)}")
    raw_lines = list(header)
    rows = []
    for code, op, threads, label, lib in jobs:
        got = drive(code, op, threads, lib, label, args.repeats)
        summary = stats.summarize(got["samples"])
        row = {
            "code": code,
            "build": label,
            "version": got["version"],
            "operation": op,
            "threads": threads,
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
            f"{code} {label} {op} threads={threads} cores={got['cores']} "
            f"affinity={got['affinity']} {got['check']} samples_ms "
            f"{stats.fmt_samples(got['samples'])}")
        print(f"{code:8s} {label:8s} {op:10s} {threads} {row['median_ms']:>9s} ms",
              flush=True)

    with open(args.out / "scaling.csv", "w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(rows[0].keys()), lineterminator="\n")
        writer.writeheader()
        writer.writerows(rows)
    (args.out / "raw-scaling.txt").write_text("\n".join(raw_lines) + "\n")
    return 0


if __name__ == "__main__":
    sys.exit(main())
