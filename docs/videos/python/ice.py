#!/usr/bin/env python3
# Copyright (c) d-SEAMS developers
"""Ice Ih census through libseams. Usage: ice.py DUMP LIBDIR"""

import ctypes
import os
import sys


def read_dump(path):
    xyz = []
    edges = []
    reading = False
    with open(path) as handle:
        while True:
            line = handle.readline()
            if line == "":
                break
            if line.startswith("ITEM: BOX BOUNDS"):
                for _ in range(3):
                    lo, hi = handle.readline().split()[:2]
                    edges.append(float(hi) - float(lo))
            elif line.startswith("ITEM: ATOMS"):
                reading = True
            elif reading:
                parts = line.split()
                if len(parts) < 5:
                    break
                xyz.extend(float(value) for value in parts[2:5])
    return xyz, edges


def main():
    path, libdir = sys.argv[1], sys.argv[2]
    xyz, edges = read_dump(path)
    n = len(xyz) // 3
    if n != 128 or len(edges) != 3:
        raise SystemExit("atom count {}".format(n))
    lib = ctypes.CDLL(os.path.join(libdir, "libseams.so"))
    lib.seams_system_create.restype = ctypes.c_void_p
    lib.seams_system_destroy.argtypes = [ctypes.c_void_p]
    lib.seams_system_set_ortho.restype = ctypes.c_int
    lib.seams_system_set_ortho.argtypes = [
        ctypes.c_void_p,
        ctypes.POINTER(ctypes.c_double),
        ctypes.c_int,
        ctypes.c_double,
        ctypes.c_double,
        ctypes.c_double,
        ctypes.c_double,
        ctypes.c_double,
        ctypes.c_double,
        ctypes.c_void_p,
        ctypes.c_void_p,
    ]
    lib.seams_neighbours_build.restype = ctypes.c_int
    lib.seams_neighbours_build.argtypes = [
        ctypes.c_void_p,
        ctypes.c_double,
        ctypes.c_int,
        ctypes.c_int,
        ctypes.POINTER(ctypes.c_void_p),
    ]
    lib.seams_neighbours_count.restype = ctypes.c_int
    lib.seams_neighbours_count.argtypes = [ctypes.c_void_p]
    lib.seams_neighbours_destroy.argtypes = [ctypes.c_void_p]
    lib.seams_q4q6.restype = ctypes.c_int
    lib.seams_q4q6.argtypes = [
        ctypes.c_void_p,
        ctypes.c_double,
        ctypes.c_int,
        ctypes.POINTER(ctypes.c_double),
        ctypes.POINTER(ctypes.c_double),
    ]
    lib.seams_chill_plus.restype = ctypes.c_int
    lib.seams_chill_plus.argtypes = [
        ctypes.c_void_p,
        ctypes.c_double,
        ctypes.c_int,
        ctypes.c_int,
        ctypes.c_int,
        ctypes.POINTER(ctypes.c_int),
        ctypes.POINTER(ctypes.c_int),
    ]
    lib.seams_rings_find.restype = ctypes.c_int
    lib.seams_rings_find.argtypes = [
        ctypes.c_void_p,
        ctypes.c_double,
        ctypes.c_int,
        ctypes.c_int,
        ctypes.c_int,
        ctypes.POINTER(ctypes.c_void_p),
    ]
    lib.seams_rings_histogram.restype = ctypes.c_int
    lib.seams_rings_histogram.argtypes = [ctypes.c_void_p, ctypes.POINTER(ctypes.c_int), ctypes.c_int]
    lib.seams_rings_destroy.argtypes = [ctypes.c_void_p]
    coords = (ctypes.c_double * len(xyz))(*xyz)
    sys_handle = lib.seams_system_create()
    if lib.seams_system_set_ortho(sys_handle, coords, n, 0, 0, 0, edges[0], edges[1], edges[2], None, None) != 0:
        raise SystemExit("frame")
    nb = ctypes.c_void_p()
    if lib.seams_neighbours_build(sys_handle, 3.5, 1, 0, ctypes.byref(nb)) != 0:
        raise SystemExit("neighbours")
    nI = lib.seams_neighbours_count(nb)
    lib.seams_neighbours_destroy(nb)
    q4 = (ctypes.c_double * n)()
    q6 = (ctypes.c_double * n)()
    if lib.seams_q4q6(sys_handle, 3.5, 1, q4, q6) != 0:
        raise SystemExit("q6")
    mean = sum(q6) / n
    labels = (ctypes.c_int * n)()
    counts = (ctypes.c_int * 9)()
    if lib.seams_chill_plus(sys_handle, 3.5, 1, 0, 0, labels, counts) != 0:
        raise SystemExit("chill")
    rings = ctypes.c_void_p()
    if lib.seams_rings_find(sys_handle, 3.5, 1, 0, 6, ctypes.byref(rings)) != 0:
        raise SystemExit("rings")
    hist = (ctypes.c_int * 13)()
    if lib.seams_rings_histogram(rings, hist, 13) != 0:
        raise SystemExit("histogram")
    print("nop {} nI {} q6 {:.6f} hexagonal {} rings {}".format(n, nI, mean, counts[1], hist[6]))
    lib.seams_rings_destroy(rings)
    lib.seams_system_destroy(sys_handle)


if __name__ == "__main__":
    main()
