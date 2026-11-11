#!/usr/bin/env python3
"""Best threaded time of mdapy, OVITO and freud at 1 and 8 threads.

The minimum median is the code's best parallel configuration on this
machine. d-SEAMS is not timed here; its distributed times are in
benchmarks/rank-scale.csv.
"""

from __future__ import annotations

import csv
import statistics
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "scripts"))
import thread_scale_time as bench

SIZES = (24, 48)
THREADS = (1, 8)
KINDS = ("mdapy", "ovito", "freud")


def main():
    out = ROOT / "benchmarks" / "h2h-parallel.csv"
    fields = ["molecules", "code", "threads", "median_ms", "iqr_ms", "best", "samples_ms"]
    rows = []
    for m in SIZES:
        molecules = 8 * m ** 3
        for kind in KINDS:
            for threads in THREADS:
                print(f"{kind} m={m} threads={threads}", flush=True)
                version, samples = bench.run_other(kind, threads, m)
                ordered = sorted(samples)
                q1, _, q3 = statistics.quantiles(ordered, n=4, method="inclusive")
                med = statistics.median(ordered)
                rows.append({
                    "molecules": molecules,
                    "code": kind,
                    "threads": threads,
                    "median_ms": f"{med:.6f}",
                    "iqr_ms": f"{q3 - q1:.6f}",
                    "best": "0",
                    "samples_ms": ",".join(f"{v:.6f}" for v in samples),
                    "version": version,
                })
                print(f"  {med:.1f} ms ({version})", flush=True)
    best = {}
    for row in rows:
        key = (row["molecules"], row["code"])
        if key not in best or float(row["median_ms"]) < float(best[key]["median_ms"]):
            best[key] = row
    for row in rows:
        key = (row["molecules"], row["code"])
        row["best"] = "1" if row is best[key] else "0"
    with out.open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields + ["version"])
        writer.writeheader()
        writer.writerows(rows)
    print(f"wrote {out}", flush=True)


if __name__ == "__main__":
    main()
