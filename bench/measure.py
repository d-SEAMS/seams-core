#!/usr/bin/env python3
"""Interleaved, order-flipped comparison of two seams_bench binaries.

Each round runs one binary then the other. Even rounds start with the
base binary; odd rounds start with the other, so a slow drift does not
land on one side. A same-binary A/A run uses the same schedule. A case
is reported faster or slower only when its bootstrap 95% CI of the
median ratio lies entirely outside the A/A CI.
"""

from __future__ import annotations

import argparse
import os
import random
import statistics
import subprocess
import sys
from pathlib import Path


CASES = ("neigh", "q3", "q6", "rings", "ice")


def parse_args() -> argparse.Namespace:
    here = Path(__file__).resolve().parent
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("--base", required=True, type=Path)
    p.add_argument("--patch", required=True, type=Path)
    p.add_argument("--rounds", type=int, default=20)
    p.add_argument("--boot", type=int, default=10000)
    p.add_argument("--seed", type=int, default=1)
    p.add_argument("--cells", type=int, default=8)
    p.add_argument("--case", action="append", dest="cases")
    p.add_argument("--out", type=Path, default=here / "local" / "compare.bench.csv")
    p.add_argument("--taskset-cpu", default="1")
    p.add_argument(
        "--lib-path",
        default="",
        help="Extra LD_LIBRARY_PATH entries, colon-separated",
    )
    return p.parse_args()


def pin_prefix(cpu: str) -> list[str]:
    if not cpu:
        return []
    from shutil import which

    if which("taskset") is None:
        return []
    return ["taskset", "-c", cpu]


def run_case(binary: Path, case: str, cells: int, cpu: str, lib_path: str) -> tuple[float, str]:
    env = os.environ.copy()
    env["OMP_NUM_THREADS"] = "1"
    env["OMP_DYNAMIC"] = "FALSE"
    env["OMP_PROC_BIND"] = "TRUE"
    env["SEAMS_BENCH_THREADS"] = "1"
    extra = lib_path or os.environ.get("SEAMS_BENCH_LIB_PATH", "")
    if extra:
        prev = env.get("LD_LIBRARY_PATH", "")
        env["LD_LIBRARY_PATH"] = extra if not prev else extra + ":" + prev
    cmd = pin_prefix(cpu) + [
        str(binary),
        "--cells",
        str(cells),
        "--case",
        case,
        "--warmup",
        "1",
        "--reps",
        "1",
    ]
    proc = subprocess.run(cmd, env=env, check=True, capture_output=True, text=True)
    seconds = None
    checksum = None
    for line in proc.stdout.splitlines():
        if not line.startswith("case "):
            continue
        parts = line.split()
        # case NAME rep I seconds S checksum H
        seconds = float(parts[5])
        checksum = parts[7]
    if seconds is None or checksum is None:
        sys.stderr.write(proc.stdout)
        sys.stderr.write(proc.stderr)
        raise RuntimeError(f"no timing line from {binary} case {case}")
    return seconds, checksum


def bootstrap_median_ci(ratios: list[float], nboot: int, seed: int) -> tuple[float, float, float]:
    rng = random.Random(seed)
    n = len(ratios)
    medians = []
    for _ in range(nboot):
        sample = [ratios[rng.randrange(n)] for _ in range(n)]
        medians.append(statistics.median(sample))
    medians.sort()
    lo_i = int(0.025 * nboot)
    hi_i = int(0.975 * nboot) - 1
    return statistics.median(ratios), medians[lo_i], medians[hi_i]


def paired_ratios(
    left: Path,
    right: Path,
    case: str,
    rounds: int,
    cells: int,
    cpu: str,
    check_equal: bool,
    lib_path: str,
) -> list[float]:
    ratios = []
    for i in range(rounds):
        if i % 2 == 0:
            t_left, c_left = run_case(left, case, cells, cpu, lib_path)
            t_right, c_right = run_case(right, case, cells, cpu, lib_path)
        else:
            t_right, c_right = run_case(right, case, cells, cpu, lib_path)
            t_left, c_left = run_case(left, case, cells, cpu, lib_path)
        if check_equal and c_left != c_right:
            raise RuntimeError(
                f"checksum mismatch on {case} round {i}: {c_left} vs {c_right}"
            )
        ratios.append(t_right / t_left)
        print(
            f"{case} round {i} left {t_left:.6f} right {t_right:.6f} "
            f"ratio {ratios[-1]:.6f}",
            flush=True,
        )
    return ratios


def outside(lo: float, hi: float, band_lo: float, band_hi: float) -> str:
    if hi < band_lo:
        return "faster"
    if lo > band_hi:
        return "slower"
    return "within"


def main() -> int:
    args = parse_args()
    if args.rounds < 20:
        print("need at least 20 rounds", file=sys.stderr)
        return 2
    cases = args.cases if args.cases else list(CASES)
    args.out.parent.mkdir(parents=True, exist_ok=True)
    lines = ["case,median_ratio,ci_lo,ci_hi,aa_lo,aa_hi,call"]
    for idx, case in enumerate(cases):
        print(f"== A/B {case}", flush=True)
        ab = paired_ratios(
            args.base,
            args.patch,
            case,
            args.rounds,
            args.cells,
            args.taskset_cpu,
            True,
            args.lib_path,
        )
        print(f"== A/A {case}", flush=True)
        aa = paired_ratios(
            args.base,
            args.base,
            case,
            args.rounds,
            args.cells,
            args.taskset_cpu,
            True,
            args.lib_path,
        )
        med, lo, hi = bootstrap_median_ci(ab, args.boot, args.seed + idx)
        _, aa_lo, aa_hi = bootstrap_median_ci(aa, args.boot, args.seed + 1000 + idx)
        call = outside(lo, hi, aa_lo, aa_hi)
        lines.append(
            f"{case},{med:.6f},{lo:.6f},{hi:.6f},{aa_lo:.6f},{aa_hi:.6f},{call}"
        )
        print(
            f"SUMMARY {case} median {med:.4f} CI [{lo:.4f}, {hi:.4f}] "
            f"A/A [{aa_lo:.4f}, {aa_hi:.4f}] {call}",
            flush=True,
        )
    args.out.write_text("\n".join(lines) + "\n")
    print(f"wrote {args.out}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
