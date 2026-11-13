#!/usr/bin/env python3
"""One sitting of the hydrate signature search at 1, 2, 4 and 8 threads.

The cells are GenIce sI and sII replicated 2x2x2, 4x4x4, 6x6x6 and 8x8x8.
Each sample is one linkcell neighbour list, one depth-6 ring network, and
one combined search of 5^12, 5^12 6^2 and 5^12 6^4. The replicate is built
in the benchmark process and is not written out.
One warm-up and at least eight repeats. The median of an even sample is
the mean of the two central values, and the interquartile range is the
type-7 upper quartile minus the lower one. The log records the core
count, the lscpu model name, and the source commit.
"""

from __future__ import annotations

import argparse
import csv
import subprocess
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "scripts"))
import paper_stats as stats

CELLS = (
    ("sI", "input/traj/genice_sI.lammpstrj", (2, 4, 6, 8)),
    ("sII", "input/traj/genice_sII.lammpstrj", (2, 4, 6, 8)),
)
THREADS = (1, 2, 4, 8)


def model_name() -> str:
    done = subprocess.run(["lscpu"], text=True, capture_output=True, check=True)
    for line in done.stdout.splitlines():
        if line.startswith("Model name:"):
            return line.split(":", 1)[1].strip()
    raise SystemExit("lscpu did not report a model name")


def source_commit() -> str:
    done = subprocess.run(
        ["git", "rev-parse", "HEAD"],
        cwd=ROOT, text=True, capture_output=True, check=True,
    )
    return done.stdout.strip()


def run_cell(bench: Path, traj: Path, tile: int, threads: int, warmup: int,
             repeats: int) -> tuple[int, str, list[float]]:
    env = stats.pin_env(threads)
    cmd = [
        "taskset", "-c", stats.taskset_cores(threads),
        str(bench), str(traj), str(tile), str(warmup), str(repeats),
    ]
    done = subprocess.run(cmd, cwd=ROOT, env=env, text=True, capture_output=True)
    if done.returncode != 0:
        tail = (done.stderr or done.stdout or "").strip().splitlines()
        raise RuntimeError(tail[-1] if tail else f"exit {done.returncode}")
    atoms = None
    census = ""
    samples: list[float] = []
    for line in done.stdout.splitlines():
        if line.startswith("atoms "):
            atoms = int(line.split()[1])
        elif line.startswith("census "):
            census = line.split(" ", 1)[1]
        elif line.startswith("ms "):
            samples.append(float(line.split()[1]))
    if atoms is None or census == "" or len(samples) != repeats:
        raise RuntimeError(f"short bench output for {traj.name} x{tile} threads {threads}")
    return atoms, census, samples


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--bench", type=Path,
                        default=ROOT / "bbdir-release/tests/bench_trace_threads")
    parser.add_argument("--out", type=Path, default=ROOT / "benchmarks")
    parser.add_argument("--warmup", type=int, default=1)
    parser.add_argument("--repeats", type=int, default=8)
    args = parser.parse_args()
    if args.repeats < 8:
        raise SystemExit("n must be at least 8")
    if not args.bench.is_file():
        raise SystemExit(f"missing {args.bench}")
    args.out.mkdir(parents=True, exist_ok=True)
    cores, affinity = stats.cores_now()
    model = model_name()
    commit = source_commit()
    log_lines = [
        f"cores {cores}",
        f"model {model}",
        f"commit {commit}",
        f"affinity {affinity}",
        "threads 1 2 4 8",
        f"warmup {args.warmup}",
        f"n {args.repeats}",
        "operation one linkcell list, one ring network, signatures 512 51262 51264",
        "driver bench_trace_threads",
        "replicates are built in memory from the GenIce unit cells",
        "one sitting, serial: each cell at 1, 2, 4 and 8 threads before the next cell",
    ]
    rows = []
    for name, rel, tiles in CELLS:
        traj = ROOT / rel
        for tile in tiles:
            for threads in THREADS:
                atoms, census, samples = run_cell(
                    args.bench, traj, tile, threads, args.warmup, args.repeats)
                summary = stats.summarize(samples)
                label = f"{name} {tile}x{tile}x{tile}"
                log_lines.append(
                    f"{label} threads {threads} atoms {atoms} census {census} "
                    f"median_ms {summary['median_ms']:.6f} iqr_ms {summary['iqr_ms']:.6f} "
                    f"samples {stats.fmt_samples(samples)}"
                )
                rows.append({
                    "cell": label,
                    "crystal": name,
                    "tile": tile,
                    "molecules": atoms,
                    "operation": "signature-cages",
                    "threads": threads,
                    "cores": cores,
                    "model": model,
                    "commit": commit,
                    "affinity": affinity,
                    "warmup": args.warmup,
                    "n": args.repeats,
                    "median_ms": f"{summary['median_ms']:.6f}",
                    "iqr_ms": f"{summary['iqr_ms']:.6f}",
                    "q1_ms": f"{summary['q1_ms']:.6f}",
                    "q3_ms": f"{summary['q3_ms']:.6f}",
                    "census": census,
                    "samples_ms": stats.fmt_samples(samples),
                })
                print(
                    f"{label} threads {threads} {summary['median_ms']:.3f} ms "
                    f"iqr {summary['iqr_ms']:.3f} census {census}",
                    flush=True,
                )
    raw = args.out / "trace-threads.txt"
    raw.write_text("\n".join(log_lines) + "\n", encoding="utf-8")
    dest = args.out / "trace-threads.csv"
    with dest.open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(rows[0]))
        writer.writeheader()
        writer.writerows(rows)
    print(f"wrote {dest}", flush=True)
    return 0


if __name__ == "__main__":
    sys.exit(main())
