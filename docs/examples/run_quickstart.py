#!/usr/bin/env python3
"""Run docs/examples/quickstart.py against dseams or libseams."""

import ctypes
import os
import runpy
import sys
import types
from pathlib import Path

ROOT = Path(__file__).resolve().parents[2]


def facade(libdir):
    os.environ["LD_LIBRARY_PATH"] = str(libdir) + os.pathsep + os.environ.get("LD_LIBRARY_PATH", "")
    lib = ctypes.CDLL(str(Path(libdir) / "libseams.so"))
    lib.seams_system_create.restype = ctypes.c_void_p
    lib.seams_system_set_ortho.restype = ctypes.c_int
    lib.seams_system_set_ortho.argtypes = [
        ctypes.c_void_p, ctypes.POINTER(ctypes.c_double), ctypes.c_int,
        ctypes.c_double, ctypes.c_double, ctypes.c_double,
        ctypes.c_double, ctypes.c_double, ctypes.c_double,
        ctypes.c_void_p, ctypes.c_void_p,
    ]
    lib.seams_chill_plus.restype = ctypes.c_int
    lib.seams_chill_plus.argtypes = [
        ctypes.c_void_p, ctypes.c_double, ctypes.c_int, ctypes.c_int, ctypes.c_int,
        ctypes.POINTER(ctypes.c_int), ctypes.POINTER(ctypes.c_int),
    ]
    lib.seams_system_destroy.argtypes = [ctypes.c_void_p]

    def chill_plus(xyz, box):
        flat = [float(value) for row in xyz for value in row]
        count = len(flat) // 3
        points = (ctypes.c_double * len(flat))(*flat)
        edge = float(box[0])
        handle = lib.seams_system_create()
        if lib.seams_system_set_ortho(handle, points, count, 0, 0, 0, edge, edge, edge, None, None) != 0:
            raise SystemExit("frame")
        labels = (ctypes.c_int * count)()
        counts = (ctypes.c_int * 9)()
        if lib.seams_chill_plus(handle, 3.5, 0, 0, 0, labels, counts) != 0:
            raise SystemExit("chill+")
        lib.seams_system_destroy(handle)
        return list(labels), list(counts)

    module = types.ModuleType("dseams")
    module.chill_plus = chill_plus
    sys.modules["dseams"] = module


def main():
    try:
        import dseams  # noqa: F401
    except ImportError:
        if len(sys.argv) < 2:
            raise SystemExit("pass the directory that holds libseams.so")
        facade(sys.argv[1])
    runpy.run_path(str(ROOT / "docs" / "examples" / "quickstart.py"), run_name="__main__")
    print("quickstart python ok")


if __name__ == "__main__":
    main()
