#!/usr/bin/env python3
# Copyright (c) d-SEAMS developers
"""Eight-atom cubic cell through the ctypes package."""

import math
import os
import sys

import numpy as np

from seams import PUBLISHED, System, set_num_threads, version


def diamond():
    a = 2.75 * 4.0 / math.sqrt(3.0)
    fcc = ((0, 0, 0), (0.5, 0.5, 0), (0.5, 0, 0.5), (0, 0.5, 0.5))
    basis = ((0, 0, 0), (0.25, 0.25, 0.25))
    xyz = [
        tuple((shift[i] + site[i]) * a for i in range(3))
        for shift in fcc
        for site in basis
    ]
    return np.asarray(xyz, dtype=np.float64), a


def main():
    if len(sys.argv) > 1:
        os.environ["SEAMS_LIB"] = sys.argv[1]
    pos, edge = diamond()
    set_num_threads(1)
    with System() as sys_handle:
        sys_handle.set_ortho(pos, (edge, edge, edge))
        ice = sys_handle.chill_plus(rule=PUBLISHED)
        order = sys_handle.q4q6()
        with sys_handle.neighbours() as nb:
            degree = int(np.diff(nb.copy().offsets).mean())
        with sys_handle.rings() as rings:
            six = int(rings.histogram(7)[6])
            fraction = rings.fraction(6)
        cn = sys_handle.coordination()
    print("version", version())
    print("n", pos.shape[0])
    print("cubic", int(ice.counts[0]))
    print("q6", round(float(order.q6.mean()), 6))
    print("degree", degree)
    print("six", six)
    print("fraction", round(fraction, 6))
    print("cn", round(cn, 6))


if __name__ == "__main__":
    main()
