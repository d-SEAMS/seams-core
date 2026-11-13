#!/usr/bin/env python3
"""Carried-cage times on the 0.1, 1 and 10 ps mW strides.

Each dump is the growth window continued from published frame 104:
2000 steps of Nose-Hoover NPT, then 200 transitions. One warm-up and
eight repeats, eight OpenMP threads. The driver is
bench_trajectory_incremental. The carried time is the median of the
per-pass mean affilInc_ms, and the ratio is the median of
sum(affilBatch_ms) / sum(affilInc_ms) on the steady frames of each pass.
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

STRIDES = {10: 0.1, 100: 1.0, 1000: 10.0}


def parse_pass(text: str) -> dict:
    batch = 0.0
    inc = 0.0
    steady = 0
    equal = True
    steady_line = ""
    for line in text.splitlines():
        if line.startswith("# steady"):
            steady_line = line
            continue
        if not line or line.startswith("#"):
            continue
        parts = line.split()
        if len(parts) < 18:
            raise RuntimeError(f"short frame line: {line}")
        frame = int(parts[0])
        if parts[17] != "yes" or parts[16] != "yes":
            equal = False
        if frame <= 1:
            continue
        batch += float(parts[9])
        inc += float(parts[10])
        steady += 1
    if steady == 0 or steady_line == "":
        raise RuntimeError("bench_trajectory_incremental wrote no steady frames")
    if not steady_line.endswith("yes"):
        equal = False
    return {
        "steady_frames": steady,
        "sum_batch_ms": batch,
        "sum_inc_ms": inc,
        "mean_inc_ms": inc / steady,
        "ratio": batch / inc if inc else float("inf"),
        "equal": equal,
        "steady": steady_line,
    }


def run_pass(bench: Path, dump: Path, frames: int, threads: int) -> str:
    env = stats.pin_env(threads)
    cmd = [
        "taskset", "-c", stats.taskset_cores(threads),
        str(bench), str(dump), str(frames), "1", "6", "2.0", "cutoff",
    ]
    done = subprocess.run(cmd, cwd=ROOT, env=env, text=True, capture_output=True)
    if done.returncode not in (0, 1):
        tail = (done.stderr or done.stdout or "").strip().splitlines()
        raise RuntimeError(tail[-1] if tail else f"exit {done.returncode}")
    return done.stdout


def count_frames(path: Path) -> int:
    count = 0
    with path.open() as handle:
        for line in handle:
            if line.startswith("ITEM: TIMESTEP"):
                count += 1
    return count


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--bench", type=Path,
                        default=ROOT / "bbdir-release/tests/bench_trajectory_incremental")
    parser.add_argument("--out", type=Path, default=ROOT / "benchmarks/paper-cells")
    parser.add_argument("--warmup", type=int, default=1)
    parser.add_argument("--repeats", type=int, default=8)
    parser.add_argument("--threads", type=int, default=8)
    args = parser.parse_args()
    if args.repeats < 8:
        raise SystemExit("n must be at least 8")
    if not args.bench.is_file():
        raise SystemExit(f"missing {args.bench}")
    args.out.mkdir(parents=True, exist_ok=True)
    cores, affinity = stats.cores_now()
    log_lines = [
        *stats.machine_lines(ROOT),
        f"cores {cores}",
        f"affinity {affinity}",
        f"threads {args.threads}",
        f"warmup {args.warmup}",
        f"n {args.repeats}",
        "graph cutoff",
        "driver bench_trajectory_incremental",
        "carried_ms median of per-pass mean affilInc_ms on frames after the first",
        "ratio median of sum(affilBatch_ms)/sum(affilInc_ms) on those frames",
    ]
    rows = []
    for stride, picoseconds in STRIDES.items():
        dump = ROOT / "data" / "hop" / f"stride-{stride}.lammpstrj"
        frames = count_frames(dump)
        if frames < 201:
            raise SystemExit(f"{dump.name} has {frames} frames, need 201")
        passes = []
        for index in range(args.warmup + args.repeats):
            text = run_pass(args.bench, dump, frames, args.threads)
            parsed = parse_pass(text)
            role = "warmup" if index < args.warmup else "sample"
            log_lines.append(
                f"stride_ps {picoseconds} {role} "
                f"mean_inc_ms {parsed['mean_inc_ms']:.6f} "
                f"sum_batch_ms {parsed['sum_batch_ms']:.6f} "
                f"sum_inc_ms {parsed['sum_inc_ms']:.6f} "
                f"ratio {parsed['ratio']:.6f} "
                f"steady_frames {parsed['steady_frames']} "
                f"equal {'yes' if parsed['equal'] else 'NO'}"
            )
            if not parsed["equal"]:
                raise SystemExit(f"stride {stride} pass {index} changed the cage set")
            if role == "sample":
                passes.append(parsed)
        carried = [item["mean_inc_ms"] for item in passes]
        ratios = [item["ratio"] for item in passes]
        carried_stats = stats.summarize(carried)
        ratio_stats = stats.summarize(ratios)
        if carried_stats["n"] != args.repeats:
            raise SystemExit("short sample")
        rows.append({
            "cell": f"carried cages, {picoseconds:g} ps",
            "stride_steps": stride,
            "stride_ps": f"{picoseconds:g}",
            "graph": "cutoff",
            "molecules": 4096,
            "frames": frames,
            "transitions": frames - 1,
            "cores": cores,
            "affinity": affinity,
            "threads": args.threads,
            "warmup": args.warmup,
            "n": args.repeats,
            "median_carried_ms": f"{carried_stats['median_ms']:.6f}",
            "iqr_carried_ms": f"{carried_stats['iqr_ms']:.6f}",
            "q1_carried_ms": f"{carried_stats['q1_ms']:.6f}",
            "q3_carried_ms": f"{carried_stats['q3_ms']:.6f}",
            "median_ratio": f"{ratio_stats['median_ms']:.6f}",
            "iqr_ratio": f"{ratio_stats['iqr_ms']:.6f}",
            "samples_carried_ms": stats.fmt_samples(carried),
            "samples_ratio": stats.fmt_samples(ratios),
            "discussion": "carried-cage time is median_carried_ms, not the ratio",
        })
        print(
            f"stride {picoseconds:g} ps carried {carried_stats['median_ms']:.4f} ms "
            f"ratio {ratio_stats['median_ms']:.4f}",
            flush=True,
        )
    raw = args.out / "raw-hop-stride.txt"
    raw.write_text("\n".join(log_lines) + "\n", encoding="utf-8")
    dest = args.out / "tab-hop-stride.csv"
    with dest.open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(rows[0]))
        writer.writeheader()
        writer.writerows(rows)
    print(f"wrote {dest}", flush=True)
    return 0


if __name__ == "__main__":
    sys.exit(main())
