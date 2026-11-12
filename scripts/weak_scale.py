#!/usr/bin/env python3
"""Weak scaling at 110592 atoms per rank.

One thread per rank. The ice is 24x24x24, 24x24x48, 24x48x48 and
48x48x48 unit cells so the atom count doubles with the rank count.
"""

from __future__ import annotations

import csv
import os
import statistics
import subprocess
import tempfile
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
_TMP = Path(tempfile.gettempdir())
BENCH = str(_TMP / "rank_bench")
_CHILL = _TMP / "builds" / "chill-mpi"
LIBS = ":".join(
    [
        str(_CHILL / "src"),
        str(_CHILL / "subprojects" / "readcon-core"),
        str(_CHILL / "subprojects" / "readcon-db"),
        str(_CHILL / "subprojects" / "linkcell"),
        str(_CHILL / "subprojects" / "vesin"),
        str(_CHILL / "subprojects" / "minimage"),
        "/usr/lib/x86_64-linux-gnu/openmpi/lib",
        "/cvmfs/software.eessi.io/versions/2026.06/software/linux/x86_64/intel/"
        "sapphirerapids/software/FlexiBLAS/3.5.0-GCC-15.2.0/lib",
    ]
)
SHAPES = ((1, 24, 24, 24), (2, 24, 24, 48), (4, 24, 48, 48), (8, 48, 48, 48))
OPS = ("neigh", "chill", "q6", "rings")
REPS = 6


def run(ranks, mx, my, mz, op):
    env = os.environ.copy()
    env.pop("LD_PRELOAD", None)
    env["LD_LIBRARY_PATH"] = LIBS
    env["OMP_NUM_THREADS"] = "1"
    env["RAYON_NUM_THREADS"] = "1"
    env["OMP_WAIT_POLICY"] = "PASSIVE"
    env["OMP_DYNAMIC"] = "FALSE"
    env["OMPI_MCA_btl"] = "^openib"
    env["OMPI_MCA_pml"] = "ob1"
    cmd = [
        "mpirun", "--bind-to", "core", "--map-by", "node:PE=1",
        "-n", str(ranks), BENCH, str(mx), str(my), str(mz), str(REPS), op,
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
        raise RuntimeError(done.stdout + done.stderr)
    return rows[:REPS]


def main():
    out = ROOT / "benchmarks" / "weak-scale.csv"
    fields = [
        "ranks", "molecules", "atoms_per_rank", "operation", "median_ms",
        "iqr_ms", "efficiency", "halo_fraction", "comm_ms", "compute_ms",
        "samples_ms",
    ]
    measured = []
    for ranks, mx, my, mz in SHAPES:
        for op in OPS:
            print(f"weak {op} ranks={ranks} {mx}x{my}x{mz}", flush=True)
            samples = run(ranks, mx, my, mz, op)
            measured.append((ranks, mx * my * mz * 8, op, samples))
            print(f"  {statistics.median(s[0] for s in samples):.1f} ms", flush=True)
    base = {}
    for ranks, molecules, op, samples in measured:
        if ranks == 1:
            base[op] = statistics.median(s[0] for s in samples)
    with out.open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields)
        writer.writeheader()
        for ranks, molecules, op, samples in measured:
            walls = [s[0] for s in samples]
            ordered = sorted(walls)
            q1, _, q3 = statistics.quantiles(ordered, n=4, method="inclusive")
            med = statistics.median(ordered)
            writer.writerow({
                "ranks": ranks,
                "molecules": molecules,
                "atoms_per_rank": molecules // ranks,
                "operation": op,
                "median_ms": f"{med:.6f}",
                "iqr_ms": f"{q3 - q1:.6f}",
                "efficiency": f"{base[op] / med:.3f}",
                "halo_fraction": f"{statistics.median(s[3] for s in samples):.6f}",
                "comm_ms": f"{statistics.median(s[1] for s in samples):.6f}",
                "compute_ms": f"{statistics.median(s[2] for s in samples):.6f}",
                "samples_ms": ",".join(f"{v:.6f}" for v in walls),
            })
    print(f"wrote {out}", flush=True)


if __name__ == "__main__":
    main()
