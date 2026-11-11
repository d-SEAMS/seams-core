"""Shared head-to-head protocol for the related-work watch.

Serial, one core (taskset -c 0), one thread. One warm-up, then at least
six timed calls. The median of an even sample is the mean of the two
central values. Quartiles are Hyndman-Fan type 7, and the IQR is Q3-Q1.
Frames are ideal ice Ic built in memory. Nothing here writes a trajectory.
"""

from __future__ import annotations

import math
import os
import time

import numpy as np

CUTOFF_A = 3.5
BOND_A = 2.75
MAX_RING = 6
WARMUP = 1
REPEATS = 6
CORE = "0"
# Conventional cells of cubic diamond. 5^3, 11^3 and 24^3 cells are
# 1000, 10648 and 110592 molecules.
REPETITIONS = (5, 11, 24)

# d-SEAMS molSys::atom_state_type, merged the way the paper merges
# reclassified cubic and hexagonal into cubic and hexagonal, and water
# into other.
DSEAMS_CLASS = {
    0: "cubic",
    1: "hexagonal",
    2: "other",
    3: "interfacial",
    4: "hydrate",
    5: "interfacial_hydrate",
    6: "other",
    7: "cubic",
    8: "hexagonal",
}

# mdapy ChillPlus codes. Hexagonal is checked before cubic; the two
# bond tests do not overlap.
MDAPY_CLASS = {
    0: "other",
    1: "hexagonal",
    2: "cubic",
    3: "interfacial",
    4: "hydrate",
    5: "interfacial_hydrate",
}


def pin_env(base=None):
    env = dict(os.environ if base is None else base)
    env.update({
        "OMP_NUM_THREADS": "1",
        "OMP_DYNAMIC": "FALSE",
        "OMP_PROC_BIND": "close",
        "MKL_NUM_THREADS": "1",
        "OPENBLAS_NUM_THREADS": "1",
        "NUMEXPR_NUM_THREADS": "1",
        "RAYON_NUM_THREADS": "1",
        "FASTATOMSTRUCT_NO_TELEMETRY": "1",
        "FASTATOMSTRUCT_SUPPRESS_TELEMETRY_NOTICE": "1",
    })
    return env


def apply_pin_env():
    for key, value in pin_env().items():
        if key not in os.environ or key.endswith("THREADS") or key.startswith("FASTATOMSTRUCT"):
            os.environ[key] = value
    os.environ["OMP_NUM_THREADS"] = "1"
    os.environ["RAYON_NUM_THREADS"] = "1"


def ice_ic(repetitions, bond=BOND_A):
    """Oxygen sites of ideal ice Ic. Second neighbours sit beyond 3.5 A."""
    lattice = 4.0 * bond / math.sqrt(3.0)
    fcc = np.array(
        [[0.0, 0.0, 0.0], [0.5, 0.5, 0.0], [0.5, 0.0, 0.5], [0.0, 0.5, 0.5]],
        dtype=np.float64,
    )
    basis = np.vstack([fcc, fcc + 0.25])
    axes = [np.arange(repetitions, dtype=np.float64) for _ in range(3)]
    ijk = np.stack(np.meshgrid(*axes, indexing="ij"), axis=-1).reshape(-1, 3)
    points = (ijk[:, None, :] + basis[None, :, :]).reshape(-1, 3) * lattice
    side = repetitions * lattice
    np.mod(points, side, out=points)
    return np.ascontiguousarray(points), float(side)


def molecules_for(repetitions):
    return 8 * repetitions ** 3


def quantile_linear(sorted_vals, q):
    if len(sorted_vals) == 1:
        return sorted_vals[0]
    pos = q * (len(sorted_vals) - 1)
    lo = int(pos)
    frac = pos - lo
    if lo + 1 >= len(sorted_vals):
        return sorted_vals[-1]
    return sorted_vals[lo] * (1.0 - frac) + sorted_vals[lo + 1] * frac


def summarize(samples_ms):
    if len(samples_ms) < REPEATS:
        raise ValueError(f"need n >= {REPEATS}, got {len(samples_ms)}")
    ordered = sorted(float(v) for v in samples_ms)
    mid = len(ordered) // 2
    if len(ordered) % 2 == 0:
        median = 0.5 * (ordered[mid - 1] + ordered[mid])
    else:
        median = ordered[mid]
    q1 = quantile_linear(ordered, 0.25)
    q3 = quantile_linear(ordered, 0.75)
    return {
        "n": len(ordered),
        "median_ms": median,
        "q1_ms": q1,
        "q3_ms": q3,
        "iqr_ms": q3 - q1,
        "min_ms": ordered[0],
        "max_ms": ordered[-1],
    }


def time_calls(prep, run, warmup=WARMUP, repeats=REPEATS):
    """prep() is untimed. run(prep()) is timed. The warm-up calls are dropped."""
    samples = []
    last = None
    for i in range(warmup + repeats):
        ready = prep()
        t0 = time.perf_counter()
        last = run(ready)
        dt = (time.perf_counter() - t0) * 1e3
        if i >= warmup:
            samples.append(dt)
    return samples, last


def cpus_allowed():
    try:
        with open("/proc/self/status", encoding="utf-8") as handle:
            for line in handle:
                if line.startswith("Cpus_allowed_list:"):
                    return line.split(":", 1)[1].strip()
    except OSError:
        return "unknown"
    return "unknown"
