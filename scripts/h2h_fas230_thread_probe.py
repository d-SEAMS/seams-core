#!/usr/bin/env python3
# Copyright (c) d-SEAMS developers
"""Probe OS thread counts for the fastatomstruct arms.

Samples /proc/<pid>/task while each worker call runs and records the
maximum observed task count for each sample.
"""

from __future__ import annotations

import argparse
import csv
import os
import subprocess
import time
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
FRAME = ROOT / "bbdir" / "fas230-scratch" / "ic-1000.npz"
OPERATIONS = ("rings",)
ARMS = ("fas220", "fas230")
REPEATS = 3


def pin_env() -> dict[str, str]:
    env = dict(os.environ)
    env.update({
        "OMP_NUM_THREADS": "1",
        "OMP_DYNAMIC": "FALSE",
        "OMP_PROC_BIND": "close",
        "OMP_PLACES": "cores",
        "OMP_THREAD_LIMIT": "1",
        "MKL_NUM_THREADS": "1",
        "MKL_DOMAIN_NUM_THREADS": "1",
        "OPENBLAS_NUM_THREADS": "1",
        "GOTO_NUM_THREADS": "1",
        "NUMEXPR_NUM_THREADS": "1",
        "RAYON_NUM_THREADS": "1",
        "BLIS_NUM_THREADS": "1",
        "VECLIB_MAXIMUM_THREADS": "1",
        "NUMBA_NUM_THREADS": "1",
        "FASTATOMSTRUCT_NO_TELEMETRY": "1",
        "FASTATOMSTRUCT_SUPPRESS_TELEMETRY_NOTICE": "1",
    })
    return env


def count_tasks(pid: int) -> int:
    task_dir = Path("/") / "proc" / str(pid) / "task"
    try:
        return sum(1 for _ in task_dir.iterdir())
    except FileNotFoundError:
        return 0


def run_one(python: str, operation: str, env: dict[str, str]) -> int:
    cmd = [
        "taskset", "-c", "0",
        python,
        str(ROOT / "scripts" / "h2h_fastatomstruct.py"),
        "--frame", str(FRAME),
        "--operation", operation,
        "--warmup", "0",
        "--repeats", "1",
        "--threads", "1",
    ]
    proc = subprocess.Popen(
        cmd,
        cwd=ROOT,
        env=env,
        stdout=subprocess.DEVNULL,
        stderr=subprocess.DEVNULL,
    )
    maximum = 0
    while True:
        maximum = max(maximum, count_tasks(proc.pid))
        if proc.poll() is not None:
            break
        time.sleep(0.001)
    if proc.returncode != 0:
        raise SystemExit(f"probe worker failed with status {proc.returncode}")
    maximum = max(maximum, count_tasks(proc.pid))
    if maximum < 1:
        raise SystemExit("probe did not observe any task")
    return maximum


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--out", default="benchmarks/fas230-thread-probe.csv")
    args = parser.parse_args()

    if not FRAME.is_file():
        raise SystemExit(f"missing frame {FRAME}")
    py220 = os.environ.get("FAS220_PYTHON")
    py230 = os.environ.get("FAS230_PYTHON")
    if not py220 or not py230:
        raise SystemExit("set FAS220_PYTHON and FAS230_PYTHON")

    arm_python = {"fas220": py220, "fas230": py230}
    env = pin_env()
    out = Path(args.out)
    out.parent.mkdir(parents=True, exist_ok=True)
    with out.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(
            handle,
            fieldnames=("arm", "size", "operation", "sample", "os_threads"),
            lineterminator="\n",
        )
        writer.writeheader()
        for arm in ARMS:
            for operation in OPERATIONS:
                for sample in range(REPEATS):
                    maximum = run_one(arm_python[arm], operation, env)
                    writer.writerow({
                        "arm": arm,
                        "size": "1000",
                        "operation": operation,
                        "sample": str(sample),
                        "os_threads": str(maximum),
                    })
    print(f"wrote {out}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
