#!/usr/bin/env python3
# Copyright (c) d-SEAMS developers
"""Run test_serial_affinity under taskset and check the library keeps the mask.

Each case starts the binary with a clean thread environment: OMP_NUM_THREADS,
OMP_PROC_BIND, OMP_PLACES and RAYON_NUM_THREADS are removed and only the
case's own settings are put back. The binary asserts that every call leaves
the process on the cpus taskset gave it, that no thread can run or did run
elsewhere, and that a one-thread width has a single thread. The first case is
taskset -c <first cpu> with OMP_NUM_THREADS=1 and RAYON_NUM_THREADS unset.

    python3 tests/test_serial_affinity.py bbdir/tests/test_serial_affinity
"""

from __future__ import annotations

import os
import shutil
import subprocess
import sys

CLEAR = ("OMP_NUM_THREADS", "OMP_PROC_BIND", "OMP_PLACES", "OMP_DYNAMIC",
         "RAYON_NUM_THREADS", "SEAMS_FAIR_WIDTH", "GOMP_CPU_AFFINITY", "KMP_AFFINITY")


def cases(cpus: list[int]) -> list[dict]:
    one = str(cpus[0])
    rows = [
        {"name": "taskset one cpu, OMP_NUM_THREADS=1", "taskset": one, "cpus": one,
         "width": 1, "env": {"OMP_NUM_THREADS": "1"}},
        {"name": "taskset one cpu, OMP_NUM_THREADS=1, OMP_PROC_BIND=close",
         "taskset": one, "cpus": one, "width": 1,
         "env": {"OMP_NUM_THREADS": "1", "OMP_PROC_BIND": "close", "OMP_PLACES": "cores"}},
    ]
    if len(cpus) >= 2:
        two = f"{cpus[0]},{cpus[1]}"
        rows += [
            {"name": "taskset two cpus, seams_set_num_threads(1)", "taskset": two,
             "cpus": two, "width": 1, "env": {}, "args": ["--set-threads", "1"]},
            {"name": "taskset two cpus, OMP_NUM_THREADS=1, RAYON_NUM_THREADS=8",
             "taskset": two, "cpus": two, "width": 1,
             "env": {"OMP_NUM_THREADS": "1", "RAYON_NUM_THREADS": "8"}},
            {"name": "taskset two cpus, two threads, OMP_PROC_BIND=close", "taskset": two,
             "cpus": two, "width": 2,
             "env": {"OMP_NUM_THREADS": "2", "OMP_PROC_BIND": "close", "OMP_PLACES": "cores"}},
        ]
    return rows


def main() -> int:
    if len(sys.argv) != 2:
        print(__doc__)
        return 2
    binary = sys.argv[1]
    taskset = shutil.which("taskset")
    if taskset is None:
        print("taskset is not on PATH")
        return 1
    cpus = sorted(os.sched_getaffinity(0))
    failed = []
    for case in cases(cpus):
        env = {key: value for key, value in os.environ.items() if key not in CLEAR}
        env.update(case["env"])
        cmd = [taskset, "-c", case["taskset"], binary, "--cpus", case["cpus"],
               "--width", str(case["width"]), *case.get("args", [])]
        print(f"== {case['name']}: {' '.join(cmd)}", flush=True)
        run = subprocess.run(cmd, env=env, capture_output=True, text=True, timeout=600)
        sys.stdout.write(run.stdout)
        sys.stdout.write(run.stderr)
        if run.returncode != 0:
            failed.append(case["name"])
    if len(cpus) < 2:
        print(f"only cpu {cpus[0]} is allowed; the two-cpu cases need a second cpu")
    for name in failed:
        print(f"FAILED {name}")
    return 1 if failed else 0


if __name__ == "__main__":
    sys.exit(main())
