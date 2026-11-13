#!/usr/bin/env python3
"""Median and interquartile range for one paper-cell sitting.

An even sample count uses the mean of the two central samples. The
interquartile range is the type-7 upper quartile minus the lower one.
The core count is the size of this process's CPU affinity. Nothing
here reads a host name.
"""

from __future__ import annotations

import os
import statistics


def cores_now() -> tuple[int, str]:
    try:
        allowed = sorted(os.sched_getaffinity(0))
    except AttributeError:
        allowed = list(range(os.cpu_count() or 1))
    if not allowed:
        return 0, ""
    if allowed == list(range(allowed[0], allowed[-1] + 1)):
        label = f"{allowed[0]}-{allowed[-1]}" if len(allowed) > 1 else str(allowed[0])
    else:
        label = ",".join(str(cpu) for cpu in allowed)
    return len(allowed), label


def summarize(samples: list[float]) -> dict[str, float]:
    ordered = sorted(float(v) for v in samples)
    count = len(ordered)
    if count < 2:
        raise ValueError(f"need at least two samples, got {count}")
    if count % 2 == 0:
        median = 0.5 * (ordered[count // 2 - 1] + ordered[count // 2])
    else:
        median = ordered[count // 2]
    q1, _q2, q3 = statistics.quantiles(ordered, n=4, method="inclusive")
    return {
        "n": count,
        "median_ms": median,
        "q1_ms": q1,
        "q3_ms": q3,
        "iqr_ms": q3 - q1,
    }


def fmt_samples(samples: list[float]) -> str:
    return " ".join(f"{float(v):.6f}" for v in samples)


def pin_env(threads: int, base: dict[str, str] | None = None) -> dict[str, str]:
    env = dict(os.environ if base is None else base)
    width = str(threads)
    env["OMP_NUM_THREADS"] = width
    env["OMP_DYNAMIC"] = "FALSE"
    env["OMP_PROC_BIND"] = "close"
    env["OMP_PLACES"] = "cores"
    env["RAYON_NUM_THREADS"] = width
    env["OPENBLAS_NUM_THREADS"] = "1"
    env["MKL_NUM_THREADS"] = "1"
    env["NUMEXPR_NUM_THREADS"] = "1"
    return env


def taskset_cores(threads: int) -> str:
    if threads < 1:
        raise ValueError("threads")
    return "0" if threads == 1 else f"0-{threads - 1}"


def _lscpu_field(label: str) -> str:
    import subprocess
    text = subprocess.check_output(["lscpu"], text=True)
    prefix = label + ":"
    for line in text.splitlines():
        if line.startswith(prefix):
            return line.split(":", 1)[1].strip()
    raise RuntimeError(f"lscpu has no {label}")


def source_commit(root) -> str:
    import subprocess
    done = subprocess.run(
        ["git", "-C", str(root), "rev-parse", "HEAD"],
        check=True, text=True, capture_output=True)
    return done.stdout.strip()


def machine(root) -> dict[str, str]:
    """Commit, CPU model and core count. The host name is not read."""
    cores, affinity = cores_now()
    return {
        "source_commit": source_commit(root),
        "cpu_model": _lscpu_field("Model name"),
        "cpu_cores": _lscpu_field("CPU(s)"),
        "cores": str(cores),
        "affinity": affinity,
    }


def machine_lines(root) -> list[str]:
    info = machine(root)
    return [
        f"source_commit {info['source_commit']}",
        f"cpu_model {info['cpu_model']}",
        f"cpu_cores {info['cpu_cores']}",
        f"affinity {info['affinity']}",
    ]
