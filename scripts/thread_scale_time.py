#!/usr/bin/env python3
"""Time CHILL+, q6, rings and neighbours at 1, 2, 4 and 8 threads.

Release -O3 only. One warm-up, six timed calls, median and type-7 IQR.
d-SEAMS takes its width from OMP_NUM_THREADS and gives that width to
Rayon. mdapy uses its own pool. OVITO reads OVITO_THREAD_COUNT before
import. freud uses freud.parallel.set_num_threads.
"""

from __future__ import annotations

import csv
import os
import platform
import statistics
import subprocess
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "scripts"))

import release_build

THREADS = (1, 2, 4, 8)
SIZES = (11, 24)  # 10648 and 110592
WARMUP = 1
REPEATS = 6
FLEXI = (
    "/cvmfs/software.eessi.io/versions/2026.06/software/linux/x86_64/intel/"
    "sapphirerapids/software/FlexiBLAS/3.5.0-GCC-15.2.0/lib"
)


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
    return ":".join(parts)


def summarize(samples):
    ordered = sorted(samples)
    med = statistics.median(ordered)
    q1, _, q3 = statistics.quantiles(ordered, n=4, method="inclusive")
    return med, q1, q3, q3 - q1


def hardware():
    nproc = subprocess.check_output(["nproc"], text=True).strip()
    nproc_all = subprocess.check_output(["nproc", "--all"], text=True).strip()
    cpu = platform.processor() or platform.machine()
    try:
        for line in Path("/proc/cpuinfo").read_text().splitlines():
            if line.startswith("model name"):
                cpu = line.split(":", 1)[1].strip()
                break
    except OSError:
        pass
    allowed = "unknown"
    try:
        for line in Path("/proc/self/status").read_text().splitlines():
            if line.startswith("Cpus_allowed_list:"):
                allowed = line.split(":", 1)[1].strip()
    except OSError:
        pass
    return f"{cpu}; nproc={nproc}; nproc_all={nproc_all}; affinity={allowed}"


def seams_env(threads, lib):
    env = os.environ.copy()
    env.pop("LD_PRELOAD", None)
    env["LD_LIBRARY_PATH"] = lib_env(lib)
    env["OMP_NUM_THREADS"] = str(threads)
    env["OMP_DYNAMIC"] = "FALSE"
    env["OMP_PROC_BIND"] = "close"
    env["RAYON_NUM_THREADS"] = str(threads)
    env["OPENBLAS_NUM_THREADS"] = "1"
    env["MKL_NUM_THREADS"] = "1"
    env["SEAMS_LIB"] = lib
    return env


def run_seams(threads, molecules_m, mode, lib):
    # h2h_dseams warms up once, then prints one line per timed call.
    cmd = ["taskset", "-c", "0-7", "/tmp/h2h_dseams", str(molecules_m), str(REPEATS), mode]
    done = subprocess.run(cmd, env=seams_env(threads, lib), text=True, capture_output=True, check=True)
    return [float(line) for line in done.stdout.split() if line.strip()]


def py_env():
    env = os.environ.copy()
    env.pop("LD_PRELOAD", None)
    env["LD_LIBRARY_PATH"] = "/usr/lib/x86_64-linux-gnu"
    env["QT_QPA_PLATFORM"] = "offscreen"
    env["FASTATOMSTRUCT_NO_TELEMETRY"] = "1"
    env["OPENBLAS_NUM_THREADS"] = "1"
    env["MKL_NUM_THREADS"] = "1"
    return env


OTHER = r"""
import json, math, os, sys, time
sys.path.insert(0, "/agent/repos/seams-core/scripts")
import h2h_protocol as proto
import numpy as np
kind, m, threads, repeats = sys.argv[1], int(sys.argv[2]), int(sys.argv[3]), int(sys.argv[4])
pos, side = proto.ice_ic(m)

def samples(prep, call):
    call(prep())
    out = []
    for _ in range(repeats):
        frame = prep()
        t0 = time.perf_counter()
        call(frame)
        out.append((time.perf_counter() - t0) * 1000.0)
    return out

if kind == "mdapy":
    import mdapy
    from mdapy import System
    from mdapy.identify import ChillPlus
    mdapy.set_num_threads(threads)
    if mdapy.get_num_threads() != threads:
        raise SystemExit(f"mdapy threads {mdapy.get_num_threads()} != {threads}")
    def prep():
        return System(pos, box=side)
    def call(system):
        return ChillPlus(cutoff=3.5).compute(system)
    ver = mdapy.__version__
elif kind == "ovito":
    os.environ["OVITO_THREAD_COUNT"] = str(threads)
    import ovito
    from ovito.data import DataCollection
    from ovito.modifiers import ChillPlusModifier
    from ovito.pipeline import Pipeline, StaticSource
    if tuple(ovito.version)[:3] != (3, 17, 2):
        raise SystemExit(f"ovito {ovito.version}")
    def prep():
        data = DataCollection()
        data.create_cell(np.diag([side, side, side]), pbc=(True, True, True))
        parts = data.create_particles(count=len(pos))
        parts.create_property("Position", data=pos)
        return data
    def call(data):
        pipe = Pipeline(source=StaticSource(data=data))
        pipe.modifiers.append(ChillPlusModifier(cutoff=3.5))
        return pipe.compute()
    ver = ".".join(str(v) for v in tuple(ovito.version)[:3])
elif kind == "freud":
    import freud
    freud.parallel.set_num_threads(threads)
    box = freud.box.Box.cube(side)
    pts = pos - side / 2.0
    query = freud.locality.AABBQuery(box, pts)
    def prep():
        return query
    def call(q):
        return q.query(pts, {"r_max": 3.5, "exclude_ii": True}).toNeighborList()
    ver = freud.__version__
else:
    raise SystemExit(kind)
print(ver)
for v in samples(prep, call):
    print(f"{v:.6f}")
"""


def run_other(kind, threads, molecules_m):
    # OVITO reads its thread count at import, so each width is a new process.
    env = py_env()
    if kind == "ovito":
        env["OVITO_THREAD_COUNT"] = str(threads)
    done = subprocess.run(
        ["taskset", "-c", "0-7", "/tmp/mdapy-venv/bin/python", "-c", OTHER,
         kind, str(molecules_m), str(threads), str(REPEATS)],
        env=env, text=True, capture_output=True, check=False,
    )
    if done.returncode != 0:
        raise RuntimeError(done.stderr)
    lines = [ln for ln in done.stdout.splitlines() if ln.strip()]
    version = lines[0]
    samples = [float(ln) for ln in lines[1:]]
    return version, samples


def main():
    lib, build = release_build.ensure_release(
        "/tmp/builds/chill-rel/bindings/libseams_lang.so")
    hw = hardware()
    out = ROOT / "benchmarks" / "thread-scale.csv"
    fields = [
        "hardware", "build", "code", "version", "threads", "operation",
        "molecules", "warmup", "n", "median_ms", "q1_ms", "q3_ms", "iqr_ms",
        "samples_ms",
    ]
    rows = []
    jobs = []
    for m in SIZES:
        molecules = 8 * m ** 3
        for threads in THREADS:
            for mode, op in (
                ("chill", "published"),
                ("ovito", "ovito-chill+"),
                ("q6", "q6"),
                ("rings", "rings"),
                ("neigh", "neighbours"),
            ):
                jobs.append(("dseams", threads, m, molecules, mode, op))
            jobs.append(("mdapy", threads, m, molecules, "mdapy", "chill+"))
            jobs.append(("ovito", threads, m, molecules, "ovito", "chill+"))
            jobs.append(("freud", threads, m, molecules, "freud", "neighbours"))
    for threads in (1, 8):
        molecules = 8 * 48 ** 3
        for mode, op in (
            ("chill", "published"),
            ("ovito", "ovito-chill+"),
            ("q6", "q6"),
            ("rings", "rings"),
            ("neigh", "neighbours"),
        ):
            jobs.append(("dseams", threads, 48, molecules, mode, op))
        jobs.append(("mdapy", threads, 48, molecules, "mdapy", "chill+"))
        jobs.append(("ovito", threads, 48, molecules, "ovito", "chill+"))
        jobs.append(("freud", threads, 48, molecules, "freud", "neighbours"))
    with out.open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields)
        writer.writeheader()
        for kind, threads, m, molecules, mode, op in jobs:
            print(f"{kind} threads={threads} m={m} {op}", flush=True)
            if kind == "dseams":
                samples = run_seams(threads, m, mode, lib)
                version = "2.14.0"
                code = "d-SEAMS"
            else:
                version, samples = run_other(mode, threads, m)
                code = {"mdapy": "mdapy", "ovito": "OVITO", "freud": "freud"}[kind]
            if len(samples) < REPEATS:
                raise SystemExit(f"short sample {kind} {threads} {m}: {samples}")
            med, q1, q3, iqr = summarize(samples[:REPEATS])
            writer.writerow({
                "hardware": hw,
                "build": build,
                "code": code,
                "version": version,
                "threads": threads,
                "operation": op,
                "molecules": molecules,
                "warmup": WARMUP,
                "n": REPEATS,
                "median_ms": f"{med:.6f}",
                "q1_ms": f"{q1:.6f}",
                "q3_ms": f"{q3:.6f}",
                "iqr_ms": f"{iqr:.6f}",
                "samples_ms": ",".join(f"{v:.6f}" for v in samples[:REPEATS]),
            })
            handle.flush()
            print(f"  {med:.3f} ms  iqr {iqr:.3f}", flush=True)
    print(f"wrote {out}", flush=True)


if __name__ == "__main__":
    main()
