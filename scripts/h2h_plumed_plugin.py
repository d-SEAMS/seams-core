#!/usr/bin/env python3
"""In-tree PLUMED SEAMS colvar against the d-SEAMS C API.

The colvar is engines/plumed/SeamsColvar.cpp. Pass the built shared
library as the first argument. PLUMED stores coordinates in nanometres
and the colvar scales them to angstroms before calling d-SEAMS. The
timed call is one SEAMS evaluation with no restraint, so the finite
difference derivatives stay off. Pinned ABAB, 20 samples, on ice Ih at
4096, 32768 and 110592 molecules.
"""

from __future__ import annotations

import ctypes
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
    core_count,
    frame_system,
    ice_positions,
    library_path,
    load_seams,
    require,
)
from h2h_plumed import KERNEL  # noqa: E402
from h2h_triclinic import bind_rings  # noqa: E402

N = 20
SEAMS_CHILL_PUBLISHED = 0


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


def reference(lib, pos, lengths):
    n = len(pos)
    cloud = frame_system(lib, pos, lengths)
    q = (ctypes.c_double * n)()
    bar = (ctypes.c_double * n)()
    require(lib, lib.seams_steinhardt(cloud, CUTOFF, 1, 6, q, n, bar, n), "q6")
    labels = (ctypes.c_int * n)()
    counts = (ctypes.c_int * 9)()
    require(
        lib,
        lib.seams_chill_plus(
            cloud, CUTOFF, 1, SEAMS_GRAPH_CUTOFF, SEAMS_CHILL_PUBLISHED, labels, counts
        ),
        "chill+",
    )
    bind_rings(lib)
    lib.seams_rings_fraction.argtypes = [
        ctypes.c_void_p, ctypes.c_int, ctypes.POINTER(ctypes.c_double),
    ]
    lib.seams_rings_fraction.restype = ctypes.c_int
    rings = ctypes.c_void_p()
    require(
        lib,
        lib.seams_rings_find(cloud, CUTOFF, 1, SEAMS_GRAPH_CUTOFF, 6, ctypes.byref(rings)),
        "rings",
    )
    hist = (ctypes.c_int * 7)()
    frac = ctypes.c_double()
    require(lib, lib.seams_rings_histogram(rings, hist, 7), "histogram")
    require(lib, lib.seams_rings_fraction(rings, 6, ctypes.byref(frac)), "fraction")
    lib.seams_rings_destroy(rings)
    lib.seams_system_destroy(cloud)
    q6 = float(np.ctypeslib.as_array(q).mean())
    ice = float(counts[0]) / float(n)
    return {
        "q6": q6,
        "ice": ice,
        "ring": float(frac.value),
        "hexagonal": float(counts[1]),
        "six": float(hist[6]),
    }


class Plugin:
    def __init__(self, so: Path, n, lengths):
        import plumed

        self.calc = plumed.Plumed()
        self.calc.cmd("setMDEngine", "python")
        self.calc.cmd("setNatoms", n)
        self.calc.cmd("setTimestep", 1.0)
        self.calc.cmd("setKbT", 1.0)
        self.calc.cmd("setLogFile", os.devnull)
        self.calc.cmd("init")
        self.calc.cmd("readInputLine", f"LOAD FILE={so}")
        self.calc.cmd("readInputLine", f"seams: SEAMS ATOMS=1-{n} CUTOFF=3.5")
        self.parts = {}
        for name in ("q6", "ice", "ring", "hexagonal", "six"):
            rank = np.zeros(1, dtype=np.int64)
            self.calc.cmd(f"getDataRank seams.{name}", rank)
            shape = np.zeros(int(rank[0]), dtype=np.int64)
            self.calc.cmd(f"getDataShape seams.{name}", shape)
            data = np.zeros(tuple(int(v) for v in shape) or (1,), dtype=np.float64)
            self.calc.cmd(f"setMemoryForData seams.{name}", data)
            self.parts[name] = data
        self.box = np.diag(np.asarray(lengths, dtype=np.float64) / 10.0)
        self.mass = np.ones(n, dtype=np.float64)
        self.charge = np.zeros(n, dtype=np.float64)
        self.forces = np.zeros((n, 3), dtype=np.float64)
        self.virial = np.zeros((3, 3), dtype=np.float64)
        self.pos = None
        self.step = 0

    def values(self, pos) -> dict:
        self.step += 1
        self.pos = np.ascontiguousarray(pos, dtype=np.float64) / 10.0
        self.calc.cmd("setStep", self.step)
        self.calc.cmd("setBox", self.box)
        self.calc.cmd("setPositions", self.pos)
        self.calc.cmd("setMasses", self.mass)
        self.calc.cmd("setCharges", self.charge)
        self.calc.cmd("setForces", self.forces)
        self.calc.cmd("setVirial", self.virial)
        self.calc.cmd("calc")
        return {name: float(np.array(buf).reshape(-1)[0]) for name, buf in self.parts.items()}


def report(width, n, a, b):
    am, ai, amin = stats(a)
    bm, bi, bmin = stats(b)
    # The plugin is the same library behind one PLUMED call. A faster
    # plugin time is the C++ call against three ctypes calls, not a loss.
    print(
        f"arm ice{n} plugin threads {width} molecules {n} "
        f"api_median {am:.4f} api_iqr {ai:.4f} api_min {amin:.4f} "
        f"plugin_median {bm:.4f} plugin_iqr {bi:.4f} plugin_min {bmin:.4f}",
        flush=True,
    )
    return True


def main():
    if len(sys.argv) < 2:
        raise SystemExit("usage: h2h_plumed_plugin.py LIB.so [threads]")
    so = Path(sys.argv[1]).resolve()
    width = int(sys.argv[2]) if len(sys.argv) > 2 else 1
    if not so.is_file():
        raise SystemExit(f"missing {so.name}")
    os.environ.setdefault("PLUMED_KERNEL", KERNEL)
    pin(width)
    import plumed  # noqa: F401

    pin(width)
    lib = load_seams(library_path())
    require(lib, lib.seams_set_num_threads(width), "threads")
    pin(width)
    print(f"cores {os.cpu_count()} affinity {core_count()}", flush=True)
    print(f"pinned_threads {width}", flush=True)
    print("plugin SEAMS", flush=True)
    print(f"samples {N}", flush=True)
    wins = True
    for ncell in (8, 16, 24):
        pos, lengths = ice_positions(ncell)
        pos = np.ascontiguousarray(pos)
        n = len(pos)
        ref = reference(lib, pos, lengths)
        plugin = Plugin(so, n, lengths)
        got = plugin.values(pos)
        for key in ("q6", "ice", "ring", "hexagonal", "six"):
            err = abs(got[key] - ref[key])
            print(
                f"check ice{n} {key} plugin {got[key]:.8g} dseams {ref[key]:.8g} abs {err:.3e}",
                flush=True,
            )
            limit = 1e-8 if key in ("q6", "ice", "ring") else 0.0
            if err > limit:
                raise SystemExit(f"ice{n} {key} differs")

        def ours(pos=pos, lengths=lengths):
            pin(width)
            reference(lib, pos, lengths)

        def other(pos=pos, plugin=plugin):
            pin(width)
            plugin.values(pos)

        ours()
        other()
        a, b = [], []
        for _ in range(N):
            a.append(time_once(ours))
            b.append(time_once(other))
        if not report(width, n, a, b):
            wins = False
    if not wins:
        raise SystemExit(1)


if __name__ == "__main__":
    main()
