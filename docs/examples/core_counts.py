"""CHILL+ and ring counts for Ice Ih and sII, through libseams.

``SEAMS_LIBDIR`` is the directory that holds ``libseams.so``. The numbers
are ``seams_chill_plus`` and ``seams_rings_histogram`` at 3.5 angstroms
on ``input/traj/ice1h.lammpstrj`` and ``input/traj/genice_sII.lammpstrj``.
"""

import ctypes
import os
from pathlib import Path

ROOT = Path(__file__).resolve().parents[2]


def _load_dump(path):
    xyz = []
    types = []
    box = None
    reading = False
    with open(path) as handle:
        for line in handle:
            if line.startswith("ITEM: BOX BOUNDS"):
                bounds = []
                for _ in range(3):
                    bounds.append([float(item) for item in handle.readline().split()])
                box = [
                    bounds[0][0],
                    bounds[1][0],
                    bounds[2][0],
                    bounds[0][1] - bounds[0][0],
                    bounds[1][1] - bounds[1][0],
                    bounds[2][1] - bounds[2][0],
                    0.0,
                    0.0,
                    0.0,
                ]
            elif line.startswith("ITEM: ATOMS"):
                reading = True
            elif reading:
                parts = line.split()
                if len(parts) < 5:
                    break
                types.append(int(parts[1]))
                xyz.extend(float(value) for value in parts[2:5])
    return xyz, types, box


def _classify(lib, path):
    xyz, types, box = _load_dump(path)
    n = len(xyz) // 3
    sys = lib.seams_system_create()
    points = (ctypes.c_double * len(xyz))(*xyz)
    box9 = (ctypes.c_double * 9)(*box)
    kinds = (ctypes.c_int * n)(*types)
    if lib.seams_system_set_frame(sys, points, n, box9, kinds, None) != 0:
        raise SystemExit("set_frame {}".format(path))
    labels = (ctypes.c_int * n)()
    counts = (ctypes.c_int * 9)()
    if lib.seams_chill_plus(sys, 3.5, 0, 0, 0, labels, counts) != 0:
        raise SystemExit("chill+ {}".format(path))
    rings = ctypes.c_void_p()
    if lib.seams_rings_find(sys, 3.5, 0, 0, 6, ctypes.byref(rings)) != 0:
        raise SystemExit("rings {}".format(path))
    hist = (ctypes.c_int * 13)()
    if lib.seams_rings_histogram(rings, hist, 13) != 0:
        raise SystemExit("histogram {}".format(path))
    lib.seams_rings_destroy(rings)
    lib.seams_system_destroy(sys)
    return n, [counts[i] for i in range(9)], [hist[i] for i in range(13)]


def check(libdir=None):
    """Print the Ice Ih and sII census and refuse a mismatch.

    Ice Ih is 128 hexagonal waters and 256 six-rings. sII is 136
    clathrate waters, 144 five-rings and 16 six-rings.

    >>> check()
    hexagonal 128
    rings6 256
    clathrate 136
    rings5 144
    """
    libdir = libdir or os.environ.get("SEAMS_LIBDIR", "")
    if not libdir:
        raise SystemExit("set SEAMS_LIBDIR to the directory that holds libseams.so")
    lib = ctypes.CDLL(os.path.join(libdir, "libseams.so"))
    lib.seams_system_create.restype = ctypes.c_void_p
    lib.seams_system_set_frame.argtypes = [
        ctypes.c_void_p,
        ctypes.POINTER(ctypes.c_double),
        ctypes.c_int,
        ctypes.POINTER(ctypes.c_double),
        ctypes.POINTER(ctypes.c_int),
        ctypes.c_void_p,
    ]
    lib.seams_system_set_frame.restype = ctypes.c_int
    lib.seams_chill_plus.argtypes = [
        ctypes.c_void_p,
        ctypes.c_double,
        ctypes.c_int,
        ctypes.c_int,
        ctypes.c_int,
        ctypes.POINTER(ctypes.c_int),
        ctypes.POINTER(ctypes.c_int),
    ]
    lib.seams_chill_plus.restype = ctypes.c_int
    lib.seams_rings_find.argtypes = [
        ctypes.c_void_p,
        ctypes.c_double,
        ctypes.c_int,
        ctypes.c_int,
        ctypes.c_int,
        ctypes.POINTER(ctypes.c_void_p),
    ]
    lib.seams_rings_find.restype = ctypes.c_int
    lib.seams_rings_histogram.argtypes = [
        ctypes.c_void_p,
        ctypes.POINTER(ctypes.c_int),
        ctypes.c_int,
    ]
    lib.seams_rings_histogram.restype = ctypes.c_int
    lib.seams_rings_destroy.argtypes = [ctypes.c_void_p]
    lib.seams_system_destroy.argtypes = [ctypes.c_void_p]
    lib.seams_set_num_threads.argtypes = [ctypes.c_int]
    lib.seams_set_num_threads.restype = ctypes.c_int
    if lib.seams_set_num_threads(1) != 0:
        raise SystemExit("threads")
    n, counts, rings = _classify(lib, str(ROOT / "input/traj/ice1h.lammpstrj"))
    if n != 128 or counts[1] != 128 or rings[6] != 256:
        raise SystemExit("ice Ih {} {} {}".format(n, counts, rings))
    n, counts, rings = _classify(lib, str(ROOT / "input/traj/genice_sII.lammpstrj"))
    if n != 136 or counts[4] != 136 or rings[5] != 144 or rings[6] != 16:
        raise SystemExit("sII {} {} {}".format(n, counts, rings))
    print("hexagonal 128")
    print("rings6 256")
    print("clathrate 136")
    print("rings5 144")


if __name__ == "__main__":
    check()
