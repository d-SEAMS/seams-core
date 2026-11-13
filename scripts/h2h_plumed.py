#!/usr/bin/env python3
"""PLUMED 2.10 COORDINATIONNUMBER and Q6 against d-SEAMS.

The switch is a step at 3.5 A, the same strict cutoff as the neighbour
list. PLUMED Q6 is scaled by sqrt(4*pi/(2*6+1)) so it is the local
Steinhardt parameter. Pinned ABAB, 20 samples.
"""

from __future__ import annotations

import ctypes
import math
import os
import sys
import time
from pathlib import Path

import numpy as np

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "scripts"))

from h2h_freud_ovito import (  # noqa: E402
    CUTOFF,
    SEAMS_GRAPH_CUTOFF,
    copy_csr,
    core_count,
    frame_system,
    ice_positions,
    library_path,
    load_seams,
    require,
)

N = 20
SCALE = math.sqrt(4.0 * math.pi / 13.0)
# D_MAX is the neighbour-list cutoff. Without it PLUMED builds a dense
# adjacency matrix and faults past a few thousand atoms. On this ice the
# first oxygen shell is at 2.76 A and the next is at 4.51 A, so 4.2 A
# contains the same neighbours as the 3.5 A step.
SWITCH = "SWITCH={CUSTOM FUNC=step(3.5-x) R_0=1 D_MAX=4.2}"
KERNEL = (
    "/cvmfs/software.eessi.io/versions/2026.06/software/linux/x86_64/intel/"
    "sapphirerapids/software/PLUMED/2.10.0-foss-2026.1/lib/libplumedKernel.so"
)


def pin(width: int) -> None:
    try:
        os.sched_setaffinity(0, set(range(width)))
    except OSError:
        pass


def stats(samples):
    arr = np.asarray(samples, dtype=np.float64)
    q1, med, q3 = np.percentile(arr, [25, 50, 75])
    return float(med), float(q3 - q1), float(arr.min())


def time_once(fn) -> float:
    t0 = time.perf_counter_ns()
    fn()
    return (time.perf_counter_ns() - t0) / 1e6


def degree(lib, sys) -> np.ndarray:
    handle = ctypes.c_void_p()
    require(
        lib,
        lib.seams_neighbours_build(sys, CUTOFF, 1, SEAMS_GRAPH_CUTOFF, ctypes.byref(handle)),
        "cutoff",
    )
    offsets, _cols = copy_csr(lib, handle)
    return np.diff(offsets).astype(np.float64)


def local_q6(lib, sys, n) -> np.ndarray:
    q = (ctypes.c_double * n)()
    bar = (ctypes.c_double * n)()
    require(lib, lib.seams_steinhardt(sys, CUTOFF, 1, 6, q, n, bar, n), "q6")
    return np.ctypeslib.as_array(q).copy()


class PlumedCV:
    def __init__(self, n, lengths, kind):
        import plumed

        if kind == "coordination":
            line = f"obs: COORDINATIONNUMBER SPECIES=1-{n} {SWITCH}"
        elif kind == "q6":
            line = f"obs: Q6 SPECIES=1-{n} {SWITCH}"
        else:
            raise ValueError(kind)
        self.kind = kind
        self.calc = plumed.Plumed()
        self.calc.cmd("setMDEngine", "python")
        self.calc.cmd("setNatoms", n)
        self.calc.cmd("setTimestep", 1.0)
        self.calc.cmd("setKbT", 1.0)
        self.calc.cmd("setLogFile", "/dev/null")
        self.calc.cmd("init")
        self.calc.cmd("readInputLine", line)
        rank = np.zeros(1, dtype=np.int64)
        self.calc.cmd("getDataRank obs", rank)
        shape = np.zeros(int(rank[0]), dtype=np.int64)
        self.calc.cmd("getDataShape obs", shape)
        self.data = np.zeros(tuple(int(v) for v in shape), dtype=np.float64)
        self.calc.cmd("setMemoryForData obs", self.data)
        self.box = np.diag(np.asarray(lengths, dtype=np.float64))
        self.mass = np.ones(n, dtype=np.float64)
        self.charge = np.zeros(n, dtype=np.float64)
        self.forces = np.zeros((n, 3), dtype=np.float64)
        self.virial = np.zeros((3, 3), dtype=np.float64)
        self.step = 0

    def values(self, pos) -> np.ndarray:
        self.step += 1
        self.calc.cmd("setStep", self.step)
        self.calc.cmd("setBox", self.box)
        self.calc.cmd("setPositions", pos)
        self.calc.cmd("setMasses", self.mass)
        self.calc.cmd("setCharges", self.charge)
        self.calc.cmd("setForces", self.forces)
        self.calc.cmd("setVirial", self.virial)
        self.calc.cmd("calc")
        out = np.array(self.data, copy=True).reshape(-1)
        if self.kind == "q6":
            out = out * SCALE
        return out


def report(op, width, n, a, b):
    am, ai, amin = stats(a)
    bm, bi, bmin = stats(b)
    win = am < bm and amin < bmin
    print(
        f"robust op {op} threads {width} molecules {n} "
        f"dseams_median {am:.4f} dseams_iqr {ai:.4f} dseams_min {amin:.4f} "
        f"plumed_median {bm:.4f} plumed_iqr {bi:.4f} plumed_min {bmin:.4f} "
        f"{'WIN' if win else 'LOSE'}",
        flush=True,
    )
    return win


def main():
    width = int(sys.argv[1]) if len(sys.argv) > 1 else 1
    os.environ.setdefault("PLUMED_KERNEL", KERNEL)
    pin(width)
    import plumed  # noqa: F401

    pin(width)
    lib = load_seams(library_path())
    require(lib, lib.seams_set_num_threads(width), "threads")
    pin(width)
    print(f"cores {os.cpu_count()} affinity {core_count()}", flush=True)
    print(f"pinned_threads {width}", flush=True)
    print("plumed 2.10.0", flush=True)
    print(f"samples {N}", flush=True)
    print(f"switch {SWITCH}", flush=True)
    wins = True
    for ncell in (8, 16, 24):
        pos, lengths = ice_positions(ncell)
        pos = np.ascontiguousarray(pos)
        n = len(pos)
        cloud = frame_system(lib, pos, lengths)
        d_cn = degree(lib, cloud)
        d_q = local_q6(lib, cloud, n)
        p_cn = PlumedCV(n, lengths, "coordination")
        p_q = PlumedCV(n, lengths, "q6")
        cn = p_cn.values(pos)
        qq = p_q.values(pos)
        print(
            f"check ice{n} coordination max_abs {np.max(np.abs(d_cn - cn)):.3e}",
            flush=True,
        )
        print(
            f"check ice{n} q6 max_abs {np.max(np.abs(d_q - qq)):.3e}",
            flush=True,
        )
        if np.max(np.abs(d_cn - cn)) > 1e-8 or np.max(np.abs(d_q - qq)) > 1e-8:
            raise SystemExit(f"ice{n} values differ")

        def dcn():
            pin(width)
            degree(lib, cloud)

        def pcn():
            pin(width)
            p_cn.values(pos)

        def dq():
            pin(width)
            local_q6(lib, cloud, n)

        def pq():
            pin(width)
            p_q.values(pos)

        dcn(); pcn(); dq(); pq()
        for op, ours, other in (("coordination", dcn, pcn), ("q6", dq, pq)):
            a, b = [], []
            for _ in range(N):
                a.append(time_once(ours))
                b.append(time_once(other))
            if not report(f"ice{n} {op}", width, n, a, b):
                wins = False
        lib.seams_system_destroy(cloud)
    if not wins:
        raise SystemExit(1)


if __name__ == "__main__":
    main()
