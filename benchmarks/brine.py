"""Airspeed Velocity timings for brine oxygen CHILL+ and neighbours.

The fair table is scripts/h2h_brine.py (warmup 1, n >= 6, median and IQR).
This file times one serial call so asv can track the same arm. Frames are
produced by input/brine/regenerate.py and are not in the tree.
"""

import os
import sys
from pathlib import Path

os.environ.setdefault("OMP_NUM_THREADS", "1")

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "scripts"))


def _frame():
    work = ROOT / "input" / "brine" / "work"
    preferred = work / "nacl-n1000-m1.0" / "dump.lammpstrj"
    if preferred.is_file():
        return preferred
    found = sorted(work.glob("*/dump.lammpstrj"))
    return found[0] if found else None


class BrineOxygen:
    def setup(self):
        path = _frame()
        if path is None:
            raise NotImplementedError("brine frames are not generated")
        import h2h_brine as bench

        self.bench = bench
        self.points, self.box = bench.read_dump(path)
        self.lib = bench.load_dseams()
        import numpy as np

        n = self.points.shape[0]
        self.xyz = np.ascontiguousarray(self.points.reshape(-1))
        self.box9 = np.array(
            [0, 0, 0, self.box[0], self.box[1], self.box[2], 0, 0, 0],
            dtype=np.float64,
        )
        self.types = np.ones(n, dtype=np.int32)
        self.labels = np.empty(n, dtype=np.int32)
        self.degree = np.empty(n, dtype=np.int32)
        self.n = n

    def time_chill_plus(self):
        import ctypes

        rc = self.lib.seams_chill_plus_species(
            self.xyz.ctypes.data_as(ctypes.POINTER(ctypes.c_double)),
            self.n,
            self.box9.ctypes.data_as(ctypes.POINTER(ctypes.c_double)),
            self.types.ctypes.data_as(ctypes.POINTER(ctypes.c_int)),
            None,
            1,
            0,
            None,
            0,
            3.5,
            b"cutoff",
            self.labels.ctypes.data_as(ctypes.POINTER(ctypes.c_int)),
            None,
            None,
            None,
        )
        if rc != 0:
            raise RuntimeError("chill+")

    def time_neighbours(self):
        import ctypes

        rc = self.lib.seams_oxygen_cutoff_degree(
            self.xyz.ctypes.data_as(ctypes.POINTER(ctypes.c_double)),
            self.n,
            self.box9.ctypes.data_as(ctypes.POINTER(ctypes.c_double)),
            self.types.ctypes.data_as(ctypes.POINTER(ctypes.c_int)),
            None,
            1,
            0,
            None,
            0,
            3.5,
            self.degree.ctypes.data_as(ctypes.POINTER(ctypes.c_int)),
        )
        if rc != 0:
            raise RuntimeError("neighbours")
