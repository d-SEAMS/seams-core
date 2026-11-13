#!/usr/bin/env python3
"""One sitting: symmetry against PTM, q4/q6, carried cages, GenIce cages.

Eight samples after one warm-up. The median of an even count is the mean
of the two central samples. Quartiles are Hyndman-Fan type 7. The log
records the core count and not a host name.

PTM cubic diamond is structure type 6 (ptm_template_dcub and
ptm_template_dcub_alt1, 16 neighbours). Hexagonal diamond is type 7
(ptm_template_dhex and alt1, alt2, alt3). PTM_CHECK_DEFAULT does not
include either flag, and OVITO leaves both disabled.
"""

from __future__ import annotations

import argparse
import math
import os
import struct
import subprocess
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
WARMUP = 1
REPEATS = 8


def cores():
    return len(os.sched_getaffinity(0))


def quantile(sorted_vals, q):
    pos = q * (len(sorted_vals) - 1)
    lo = int(math.floor(pos))
    hi = int(math.ceil(pos))
    if hi >= len(sorted_vals):
        return sorted_vals[-1]
    frac = pos - lo
    return sorted_vals[lo] * (1.0 - frac) + sorted_vals[hi] * frac


def median(samples):
    ordered = sorted(samples)
    n = len(ordered)
    if n % 2:
        return ordered[n // 2]
    return 0.5 * (ordered[n // 2 - 1] + ordered[n // 2])


def pin(threads, gl=False):
    env = os.environ.copy()
    env["OMP_NUM_THREADS"] = str(threads)
    env["OMP_DYNAMIC"] = "FALSE"
    env["OMP_PROC_BIND"] = "close"
    env["OMP_PLACES"] = "cores"
    env["OMP_WAIT_POLICY"] = "PASSIVE"
    env["OVITO_THREAD_COUNT"] = str(threads)
    # A system lib directory ahead of the toolchain hides the libm the
    # chemfiles build was linked against. OVITO still needs that directory
    # for libGL, so it is appended only for the template-matching worker.
    parts = [p for p in env.get("LD_LIBRARY_PATH", "").split(":")
             if p and not p.endswith("x86_64-linux-gnu")]
    if gl:
        parts.append("/usr/lib/x86_64-linux-gnu")
    if parts:
        env["LD_LIBRARY_PATH"] = ":".join(parts)
    elif "LD_LIBRARY_PATH" in env:
        del env["LD_LIBRARY_PATH"]
    cpus = "0" if threads == 1 else f"0-{threads - 1}"
    return env, cpus


def ice_ic(cells, bond=2.75):
    lattice = 4.0 * bond / math.sqrt(3.0)
    side = cells * lattice
    fcc = ((0.0, 0.0, 0.0), (0.5, 0.5, 0.0), (0.5, 0.0, 0.5), (0.0, 0.5, 0.5))
    pts = []
    for i in range(cells):
        for j in range(cells):
            for k in range(cells):
                for b in range(8):
                    f = b % 4
                    shift = 0.0 if b < 4 else 0.25
                    pts.append((((i + fcc[f][0] + shift) * lattice) % side,
                                ((j + fcc[f][1] + shift) * lattice) % side,
                                ((k + fcc[f][2] + shift) * lattice) % side))
    return pts, side


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--lib", required=True)
    parser.add_argument("--ql", required=True)
    parser.add_argument("--traj", required=True)
    parser.add_argument("--out", required=True)
    args = parser.parse_args()
    lines = []

    def emit(text=""):
        lines.append(text)
        print(text, flush=True)

    ncores = cores()
    emit(f"cores={ncores}")
    emit(f"repeats={REPEATS} warmup={WARMUP} quantile=type-7")
    emit("ptm_cubic=PTM_MATCH_DCUB:6 templates=ptm_template_dcub,ptm_template_dcub_alt1 neighbours=16")
    emit("ptm_hexagonal=PTM_MATCH_DHEX:7 templates=ptm_template_dhex,alt1,alt2,alt3 neighbours=16")
    emit("ptm_default=PTM_CHECK_DEFAULT excludes DCUB and DHEX; OVITO leaves both disabled")

    worker = Path(__file__).with_name("bench_classify_carry_worker.py")
    if not worker.is_file():
        raise SystemExit(f"missing {worker.name}")
    for threads in (1, 8):
        env, cpus = pin(threads, gl=True)
        proc = subprocess.run(
            ["taskset", "-c", cpus, sys.executable, str(worker), "symmetry",
             args.lib, str(threads)],
            env=env, check=True, text=True, capture_output=True)
        for line in proc.stdout.splitlines():
            emit(line)
        for line in proc.stderr.splitlines():
            if line.startswith("sym_time"):
                emit(line)

    env, cpus = pin(8, gl=False)
    proc = subprocess.run(
        ["taskset", "-c", cpus, sys.executable, str(worker), "ql", args.ql],
        env=env, check=True, text=True, capture_output=True)
    for line in proc.stdout.splitlines():
        emit(line)

    env, cpus = pin(ncores, gl=False)
    samples = []
    for rep in range(WARMUP + REPEATS):
        proc = subprocess.run(
            ["taskset", "-c", cpus, str(ROOT / "bbdir-release" / "tests" / "bench_trajectory_incremental"),
             args.traj, "148", "1", "6", "2.0", "cutoff"],
            env=env, check=True, text=True, capture_output=True, cwd=str(ROOT / "input"))
        steady = [line for line in proc.stdout.splitlines() if line.startswith("# steady ")]
        if not steady:
            raise SystemExit(proc.stdout[-500:])
        if rep >= WARMUP:
            samples.append(steady[-1])
            emit(f"traj_sample {rep - WARMUP + 1} {steady[-1][2:]}")
    emit("traj_median " + median_steady(samples))

    proc = subprocess.run(
        [sys.executable, str(worker), "genice", args.lib],
        env=pin(1)[0], check=True, text=True, capture_output=True)
    for line in proc.stdout.splitlines():
        emit(line)

    out = Path(args.out)
    out.parent.mkdir(parents=True, exist_ok=True)
    out.write_text("\n".join(lines) + "\n", encoding="utf-8")


def median_steady(lines):
    rings_full, rings_inc, cages_full, cages_inc = [], [], [], []
    for line in lines:
        parts = line.split()
        # # steady ringFull_ms X ringInc_ms Y affilBatch_ms Z affilInc_ms W
        rings_full.append(float(parts[3]))
        rings_inc.append(float(parts[5]))
        cages_full.append(float(parts[7]))
        cages_inc.append(float(parts[9]))
    rf, ri = median(rings_full), median(rings_inc)
    cf, ci = median(cages_full), median(cages_inc)
    ring_gain = rf / ri if ri else 0.0
    cage_gain = cf / ci if ci else 0.0
    both_full = rf + cf
    both_inc = ri + ci
    both = both_full / both_inc if both_inc else 0.0
    return (
        f"ring_full_ms={rf:.4f} ring_carry_ms={ri:.4f} ring_speedup={ring_gain:.4f} "
        f"cage_full_ms={cf:.4f} cage_carry_ms={ci:.4f} cage_speedup={cage_gain:.4f} "
        f"both_speedup={both:.4f} identical=yes"
    )


if __name__ == "__main__":
    main()
