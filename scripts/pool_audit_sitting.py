#!/usr/bin/env python3
"""One sitting of the pool-sized kernels at 1, 2, 4 and 8 threads.

Ice Ic at 4096, 32768, 110592 and 884736 molecules is built in the
benchmark and is not written out. The kernels are primitive rings,
local and cutoff q4 and q6, CHILL+, the shell label, and one carried
hop-ball update. One warm-up and at least eight repeats. The median of
an even sample is the mean of the two central values, and the
interquartile range is the type-7 upper quartile minus the lower one.
The log records the core count, the lscpu model name, and the source
commit.
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

KERNELS = (
    "rings",
    "q4-local",
    "q6-local",
    "q4-cutoff",
    "q6-cutoff",
    "chill",
    "symmetry",
    "hop",
)
MOLECULES = (4096, 32768, 110592, 884736)
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


def run_kernel(bench: Path, kernel: str, molecules: int, threads: int,
               warmup: int, repeats: int) -> tuple[int, list[float]]:
    env = stats.pin_env(threads)
    cmd = [
        "taskset", "-c", stats.taskset_cores(threads),
        str(bench), kernel, str(molecules), str(warmup), str(repeats),
    ]
    done = subprocess.run(cmd, cwd=ROOT, env=env, text=True, capture_output=True)
    if done.returncode != 0:
        tail = (done.stderr or done.stdout or "").strip().splitlines()
        raise RuntimeError(
            f"{kernel} {molecules} threads {threads}: "
            + (tail[-1] if tail else f"exit {done.returncode}"))
    fingerprint = 0
    samples: list[float] = []
    for line in done.stdout.splitlines():
        if line.startswith("fingerprint "):
            fingerprint = int(line.split()[1])
        elif line.startswith("ms "):
            samples.append(float(line.split()[1]))
    if fingerprint == 0 or len(samples) != repeats:
        raise RuntimeError(f"short bench output for {kernel} {molecules} threads {threads}")
    return fingerprint, samples


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--bench", type=Path,
                        default=ROOT / "bbdir-release/tests/bench_pool_audit")
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
        "ice Ic built in memory, bond 2.75 A, cutoff 3.5 A",
        "kernels rings q4-local q6-local q4-cutoff q6-cutoff chill symmetry hop",
        "hop times one carried update after a prime on the unperturbed crystal",
        "one sitting, serial: each kernel and size at 1, 2, 4 and 8 threads",
    ]
    rows = []
    for kernel in KERNELS:
        for molecules in MOLECULES:
            for threads in THREADS:
                fingerprint, samples = run_kernel(
                    args.bench, kernel, molecules, threads, args.warmup, args.repeats)
                summary = stats.summarize(samples)
                log_lines.append(
                    f"{kernel} molecules {molecules} threads {threads} "
                    f"median_ms {summary['median_ms']:.6f} iqr_ms {summary['iqr_ms']:.6f} "
                    f"fingerprint {fingerprint} samples {stats.fmt_samples(samples)}"
                )
                rows.append({
                    "kernel": kernel,
                    "crystal": "Ic",
                    "molecules": molecules,
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
                    "fingerprint": fingerprint,
                    "samples_ms": stats.fmt_samples(samples),
                })
                print(
                    f"{kernel} {molecules} threads {threads} {summary['median_ms']:.3f} ms "
                    f"iqr {summary['iqr_ms']:.3f}",
                    flush=True,
                )
    raw = args.out / "pool-audit.txt"
    raw.write_text("\n".join(log_lines) + "\n", encoding="utf-8")
    dest = args.out / "pool-audit.csv"
    with dest.open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(rows[0]))
        writer.writeheader()
        writer.writerows(rows)
    print(f"wrote {dest}", flush=True)
    return 0


if __name__ == "__main__":
    sys.exit(main())
