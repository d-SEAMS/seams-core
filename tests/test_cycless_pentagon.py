#!/usr/bin/env python3
"""The interface pentagon cycless misses is in the d-SEAMS census."""

from __future__ import annotations

import ctypes
import sys
from pathlib import Path

import numpy as np

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "scripts"))

from mixed_phase_frames import interface

WITNESS = frozenset((123, 298, 260, 337, 383))


def library_path() -> Path:
    if len(sys.argv) > 1:
        return Path(sys.argv[1])
    for candidate in (
        ROOT / "bbdir-release" / "src" / "libseams.so",
        ROOT / "build" / "src" / "libseams.so",
    ):
        if candidate.is_file():
            return candidate
    raise SystemExit("libseams.so was not found")


def main() -> None:
    path = library_path()
    yoda = path.with_name("libyodaLib.so")
    if yoda.is_file():
        ctypes.CDLL(str(yoda), mode=ctypes.RTLD_GLOBAL)
    lib = ctypes.CDLL(str(path))
    lib.seams_system_create.restype = ctypes.c_void_p
    lib.seams_system_destroy.argtypes = [ctypes.c_void_p]
    lib.seams_system_set_frame.argtypes = [
        ctypes.c_void_p, ctypes.POINTER(ctypes.c_double), ctypes.c_int,
        ctypes.POINTER(ctypes.c_double), ctypes.POINTER(ctypes.c_int),
        ctypes.POINTER(ctypes.c_int),
    ]
    lib.seams_system_set_frame.restype = ctypes.c_int
    lib.seams_rings_find.argtypes = [
        ctypes.c_void_p, ctypes.c_double, ctypes.c_int, ctypes.c_int, ctypes.c_int,
        ctypes.POINTER(ctypes.c_void_p),
    ]
    lib.seams_rings_find.restype = ctypes.c_int
    lib.seams_rings_destroy.argtypes = [ctypes.c_void_p]
    lib.seams_rings_count.argtypes = [ctypes.c_void_p]
    lib.seams_rings_count.restype = ctypes.c_int
    lib.seams_rings_vertices.argtypes = [
        ctypes.c_void_p, ctypes.c_int, ctypes.POINTER(ctypes.c_int),
        ctypes.c_size_t, ctypes.POINTER(ctypes.c_size_t),
    ]
    lib.seams_rings_vertices.restype = ctypes.c_int

    pos, (lx, ly, lz) = interface()
    held = np.ascontiguousarray(pos, dtype=np.float64)
    n = int(held.shape[0])
    box = (ctypes.c_double * 9)(0, 0, 0, lx, ly, lz, 0, 0, 0)
    types = (ctypes.c_int * n)(*([1] * n))
    ids = (ctypes.c_int * n)(*range(1, n + 1))
    sys = lib.seams_system_create()
    if lib.seams_system_set_frame(
        sys, held.ctypes.data_as(ctypes.POINTER(ctypes.c_double)), n, box, types, ids
    ) != 0:
        raise SystemExit("set_frame failed")
    rings = ctypes.c_void_p()
    if lib.seams_rings_find(sys, 3.5, 1, 0, 6, ctypes.byref(rings)) != 0:
        raise SystemExit("rings failed")
    found = False
    buf = (ctypes.c_int * 12)()
    n_out = ctypes.c_size_t()
    for index in range(lib.seams_rings_count(rings)):
        if lib.seams_rings_vertices(rings, index, buf, 12, ctypes.byref(n_out)) != 0:
            raise SystemExit("vertices failed")
        if frozenset(int(buf[k]) for k in range(n_out.value)) == WITNESS:
            found = True
            break
    lib.seams_rings_destroy(rings)
    lib.seams_system_destroy(sys)
    if not found:
        raise SystemExit(f"pentagon {sorted(WITNESS)} is missing from the census")
    print(f"pentagon {sorted(WITNESS)} is in the census")


if __name__ == "__main__":
    main()
