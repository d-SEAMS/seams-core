#!/usr/bin/env python3
"""Strong scaling of the spatial-domain split.

Ranks and threads are each 1, 2, 4 or 8, and their product is at most 8.
One warm-up, six timed calls. Speedup and efficiency are against one rank
and one thread. Each sample line from the benchmark is
wall, ghost-exchange, compute, halo fraction.
"""

from __future__ import annotations

import csv
import os
import statistics
import subprocess
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
BENCH = "/tmp/rank_bench"
LIBS = ":".join(
    [
        "/tmp/builds/chill-mpi/src",
        "/tmp/builds/chill-mpi/subprojects/readcon-core",
        "/tmp/builds/chill-mpi/subprojects/readcon-db",
        "/tmp/builds/chill-mpi/subprojects/linkcell",
        "/tmp/builds/chill-mpi/subprojects/vesin",
        "/tmp/builds/chill-mpi/subprojects/minimage",
        "/usr/lib/x86_64-linux-gnu/openmpi/lib",
        "/cvmfs/software.eessi.io/versions/2026.06/software/linux/x86_64/intel/"
        "sapphirerapids/software/FlexiBLAS/3.5.0-GCC-15.2.0/lib",
    ]
)
OPS = ("neigh", "chill", "q6", "rings")
SIZES = (24, 48)
REPS = 6


def configs():
    for ranks in (1, 2, 4, 8):
        for threads in (1, 2, 4, 8):
            if ranks * threads <= 8:
                yield ranks, threads


def run(cells, ranks, threads, op):
    env = os.environ.copy()
    env.pop("LD_PRELOAD", None)
    env["LD_LIBRARY_PATH"] = LIBS
    env["OMP_NUM_THREADS"] = str(threads)
    env["RAYON_NUM_THREADS"] = str(threads)
    env["OMP_WAIT_POLICY"] = "PASSIVE"
    env["OMP_DYNAMIC"] = "FALSE"
    env["OMPI_MCA_btl"] = "^openib"
    env["OMPI_MCA_pml"] = "ob1"
    mapping = ["--bind-to", "core", "--map-by", f"node:PE={threads}"]
    cmd = [
        "mpirun",
        *mapping,
        "-n",
        str(ranks),
        BENCH,
        str(cells),
        str(REPS),
        op,
    ]
    done = subprocess.run(cmd, env=env, text=True, capture_output=True, check=False)
    if done.returncode != 0:
        raise RuntimeError(done.stderr or done.stdout)
    rows = []
    for line in done.stdout.splitlines():
        parts = line.split()
        if len(parts) == 4:
            rows.append(tuple(float(v) for v in parts))
    if len(rows) < REPS:
        raise RuntimeError(f"short {op} {cells} {ranks}x{threads}: {done.stdout} {done.stderr}")
    return rows[:REPS]


def main():
    out = ROOT / "benchmarks" / "rank-scale.csv"
    fields = [
        "molecules", "operation", "ranks", "threads_per_rank", "cores",
        "median_ms", "q1_ms", "q3_ms", "iqr_ms", "speedup", "efficiency",
        "halo_fraction", "comm_ms", "compute_ms", "samples_ms",
    ]
    measured = []
    for cells in SIZES:
        for op in OPS:
            for ranks, threads in configs():
                print(f"{op} cells={cells} ranks={ranks} threads={threads}", flush=True)
                samples = run(cells, ranks, threads, op)
                measured.append((cells, op, ranks, threads, samples))
                print(f"  {statistics.median(s[0] for s in samples):.1f} ms", flush=True)
    base = {}
    for cells, op, ranks, threads, samples in measured:
        if ranks == 1 and threads == 1:
            base[(cells, op)] = statistics.median(s[0] for s in samples)
    with out.open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields)
        writer.writeheader()
        for cells, op, ranks, threads, samples in measured:
            walls = [s[0] for s in samples]
            ordered = sorted(walls)
            q1, _, q3 = statistics.quantiles(ordered, n=4, method="inclusive")
            med = statistics.median(ordered)
            cores = ranks * threads
            speed = base[(cells, op)] / med
            writer.writerow({
                "molecules": 8 * cells ** 3,
                "operation": op,
                "ranks": ranks,
                "threads_per_rank": threads,
                "cores": cores,
                "median_ms": f"{med:.6f}",
                "q1_ms": f"{q1:.6f}",
                "q3_ms": f"{q3:.6f}",
                "iqr_ms": f"{q3 - q1:.6f}",
                "speedup": f"{speed:.3f}",
                "efficiency": f"{speed / cores:.3f}",
                "halo_fraction": f"{statistics.median(s[3] for s in samples):.6f}",
                "comm_ms": f"{statistics.median(s[1] for s in samples):.6f}",
                "compute_ms": f"{statistics.median(s[2] for s in samples):.6f}",
                "samples_ms": ",".join(f"{v:.6f}" for v in walls),
            })
    print(f"wrote {out}", flush=True)


if __name__ == "__main__":
    main()
