#!/usr/bin/env python3
"""Neighbours, rings, and CHILL+ from Python.

The frame is a context manager. Labels, the ring histogram, and the
CHILL+ counts are NumPy arrays that the library fills in place.
"""

import ctypes
import os
import sys

import numpy as np


class Frame:
    def __init__(self, path):
        lib = os.environ.get("SEAMS_LANG_LIB", "libseams_lang.so")
        self._lib = ctypes.CDLL(lib)
        self._lib.seams_frame_read.restype = ctypes.c_void_p
        self._lib.seams_frame_read.argtypes = [
            ctypes.c_char_p,
            ctypes.c_int,
            ctypes.c_int,
        ]
        self._lib.seams_frame_free.argtypes = [ctypes.c_void_p]
        self._lib.seams_frame_nop.restype = ctypes.c_int
        self._lib.seams_frame_nop.argtypes = [ctypes.c_void_p]
        self._lib.seams_frame_neighbours.restype = ctypes.c_int
        self._lib.seams_frame_neighbours.argtypes = [
            ctypes.c_void_p,
            ctypes.c_double,
            ctypes.c_int,
            ctypes.c_char_p,
        ]
        self._lib.seams_frame_neighbour_pairs.restype = ctypes.c_int
        self._lib.seams_frame_neighbour_pairs.argtypes = [ctypes.c_void_p]
        self._lib.seams_frame_rings.restype = ctypes.c_int
        self._lib.seams_frame_rings.argtypes = [
            ctypes.c_void_p,
            ctypes.c_double,
            ctypes.c_int,
            ctypes.c_char_p,
            ctypes.c_int,
            ctypes.POINTER(ctypes.c_int),
            ctypes.c_int,
        ]
        self._lib.seams_frame_chill_plus.restype = ctypes.c_int
        self._lib.seams_frame_chill_plus.argtypes = [
            ctypes.c_void_p,
            ctypes.c_double,
            ctypes.c_int,
            ctypes.c_char_p,
            ctypes.POINTER(ctypes.c_int),
            ctypes.POINTER(ctypes.c_int),
        ]
        self._lib.seams_last_error.restype = ctypes.c_char_p
        self.ptr = self._lib.seams_frame_read(os.fsencode(path), 1, -1)
        if not self.ptr:
            raise RuntimeError(self._error())

    def _error(self):
        text = self._lib.seams_last_error()
        return text.decode() if text else "seams call failed"

    def __enter__(self):
        return self

    def __exit__(self, exc_type, exc, tb):
        if self.ptr:
            self._lib.seams_frame_free(self.ptr)
            self.ptr = None
        return False

    def pairs(self):
        if self._lib.seams_frame_neighbours(self.ptr, 3.5, -1, b"cutoff") != 0:
            raise RuntimeError(self._error())
        return int(self._lib.seams_frame_neighbour_pairs(self.ptr))

    def rings(self):
        counts = np.zeros(7, dtype=np.int32)
        if (
            self._lib.seams_frame_rings(
                self.ptr,
                3.5,
                -1,
                b"cutoff",
                6,
                counts.ctypes.data_as(ctypes.POINTER(ctypes.c_int)),
                7,
            )
            != 0
        ):
            raise RuntimeError(self._error())
        return counts

    def chill(self):
        n = int(self._lib.seams_frame_nop(self.ptr))
        labels = np.zeros(n, dtype=np.int32)
        counts = np.zeros(9, dtype=np.int32)
        if (
            self._lib.seams_frame_chill_plus(
                self.ptr,
                3.5,
                -1,
                b"cutoff",
                labels.ctypes.data_as(ctypes.POINTER(ctypes.c_int)),
                counts.ctypes.data_as(ctypes.POINTER(ctypes.c_int)),
            )
            != 0
        ):
            raise RuntimeError(self._error())
        return labels, counts


def main():
    root = os.environ["DSEAMS_ROOT"]
    path = os.path.join(root, "input", "traj", "ice1h.lammpstrj")
    with Frame(path) as frame:
        pairs = frame.pairs()
        rings = frame.rings()
        _labels, counts = frame.chill()
    ring_text = " ".join(str(int(value)) for value in rings)
    chill_text = " ".join(str(int(value)) for value in counts)
    print(f"pairs {pairs} rings {ring_text} chill {chill_text}")


if __name__ == "__main__":
    try:
        main()
    except Exception as err:
        print(err, file=sys.stderr)
        sys.exit(1)
