#!/usr/bin/env python3
"""Eight-atom CHILL+ and q4/q6 through dseams, or through libseams."""

import ctypes
import math
import os
import sys


def diamond():
    a = 2.75 * 4.0 / math.sqrt(3.0)
    fcc = ((0, 0, 0), (0.5, 0.5, 0), (0.5, 0, 0.5), (0, 0.5, 0.5))
    basis = ((0, 0, 0), (0.25, 0.25, 0.25))
    xyz = []
    for shift in fcc:
        for site in basis:
            xyz.extend((shift[i] + site[i]) * a for i in range(3))
    return xyz, a


def via_dseams(xyz, edge):
    import numpy as np

    import dseams

    pos = np.asarray(xyz, dtype=np.float64).reshape(8, 3)
    box = np.asarray([edge, edge, edge], dtype=np.float64)
    dseams.set_num_threads(1)
    if dseams.num_threads() != 1:
        raise SystemExit("dseams thread width")
    labels, counts = dseams.chill_plus(pos, box)
    q4, q6 = dseams.q4q6(pos, box)
    if int(counts[0]) != 8 or len(labels) != 8 or len(q6) != 8:
        raise SystemExit("dseams counts")
    if float(np.mean(q6)) <= 0.3:
        raise SystemExit("dseams q6")
    print("python api ok dseams", dseams.version())


def via_ctypes(libdir, xyz, edge):
    lib = ctypes.CDLL(os.path.join(libdir, "libseams.so"))
    lib.seams_version_string.restype = ctypes.c_char_p
    lib.seams_set_num_threads.argtypes = [ctypes.c_int]
    lib.seams_set_num_threads.restype = ctypes.c_int
    lib.seams_get_num_threads.restype = ctypes.c_int
    lib.seams_system_create.restype = ctypes.c_void_p
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
    lib.seams_system_set_ortho.restype = ctypes.c_int
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
    lib.seams_q4q6.argtypes = [
        ctypes.c_void_p,
        ctypes.c_double,
        ctypes.c_int,
        ctypes.POINTER(ctypes.c_double),
        ctypes.POINTER(ctypes.c_double),
    ]
    lib.seams_q4q6.restype = ctypes.c_int
    lib.seams_system_destroy.argtypes = [ctypes.c_void_p]
    if lib.seams_set_num_threads(1) != 0 or lib.seams_get_num_threads() != 1:
        raise SystemExit("ctypes thread width")
    points = (ctypes.c_double * len(xyz))(*xyz)
    labels = (ctypes.c_int * 8)()
    counts = (ctypes.c_int * 9)()
    q4 = (ctypes.c_double * 8)()
    q6 = (ctypes.c_double * 8)()
    sys_handle = lib.seams_system_create()
    if lib.seams_system_set_ortho(sys_handle, points, 8, 0, 0, 0, edge, edge, edge, None, None) != 0:
        raise SystemExit("ctypes frame")
    if lib.seams_chill_plus(sys_handle, 3.5, 0, 0, 0, labels, counts) != 0:
        raise SystemExit("ctypes chill+")
    if lib.seams_q4q6(sys_handle, 3.5, 0, q4, q6) != 0:
        raise SystemExit("ctypes q4q6")
    mean = sum(q6) / 8.0
    lib.seams_system_destroy(sys_handle)
    if counts[0] != 8 or mean <= 0.3:
        raise SystemExit("ctypes result")
    print("python api ok ctypes", lib.seams_version_string().decode())


def main():
    xyz, edge = diamond()
    try:
        via_dseams(xyz, edge)
    except ImportError:
        if len(sys.argv) < 2:
            raise SystemExit("pass the directory that holds libseams.so")
        via_ctypes(sys.argv[1], xyz, edge)


if __name__ == "__main__":
    main()
