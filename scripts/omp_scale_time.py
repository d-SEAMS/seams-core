#!/usr/bin/env python3
"""Time the C API at one thread and at eight, in one process.

Eight samples follow one warm-up. The record is the core count of the
machine. The host name is not part of it.
"""

from __future__ import annotations

import ctypes
import math
import os
import statistics
import tempfile
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
SAMPLES = 8


def library_path() -> Path:
    env = os.environ.get("SEAMS_LIB", "")
    if env:
        return Path(env)
    for candidate in (
        ROOT / "bbdir-release" / "src" / "libseams.so",
        ROOT / "build" / "src" / "libseams.so",
    ):
        if candidate.is_file():
            return candidate
    raise SystemExit("libseams.so was not found; set SEAMS_LIB")


def load_library(path: Path):
    lib = ctypes.CDLL(str(path))
    lib.seams_set_num_threads.argtypes = [ctypes.c_int]
    lib.seams_set_num_threads.restype = ctypes.c_int
    lib.seams_system_create.restype = ctypes.c_void_p
    lib.seams_system_destroy.argtypes = [ctypes.c_void_p]
    lib.seams_system_set_frame.argtypes = [
        ctypes.c_void_p,
        ctypes.POINTER(ctypes.c_double),
        ctypes.c_int,
        ctypes.POINTER(ctypes.c_double),
        ctypes.POINTER(ctypes.c_int),
        ctypes.POINTER(ctypes.c_int),
    ]
    lib.seams_system_set_frame.restype = ctypes.c_int
    lib.seams_system_read.argtypes = [ctypes.c_void_p, ctypes.c_char_p, ctypes.c_int, ctypes.c_int]
    lib.seams_system_read.restype = ctypes.c_int
    lib.seams_neighbours_knn.argtypes = [
        ctypes.c_void_p,
        ctypes.c_int,
        ctypes.c_double,
        ctypes.c_int,
        ctypes.c_int,
        ctypes.POINTER(ctypes.c_void_p),
    ]
    lib.seams_neighbours_knn.restype = ctypes.c_int
    lib.seams_neighbours_destroy.argtypes = [ctypes.c_void_p]
    lib.seams_steinhardt.argtypes = [
        ctypes.c_void_p,
        ctypes.c_double,
        ctypes.c_int,
        ctypes.c_int,
        ctypes.POINTER(ctypes.c_double),
        ctypes.c_size_t,
        ctypes.POINTER(ctypes.c_double),
        ctypes.c_size_t,
    ]
    lib.seams_steinhardt.restype = ctypes.c_int
    lib.seams_chill.argtypes = [
        ctypes.c_void_p,
        ctypes.c_double,
        ctypes.c_int,
        ctypes.c_int,
        ctypes.c_int,
        ctypes.POINTER(ctypes.c_int),
        ctypes.c_size_t,
        ctypes.POINTER(ctypes.c_int),
    ]
    lib.seams_chill.restype = ctypes.c_int
    lib.seams_last_error.restype = ctypes.c_char_p
    return lib


def ice_ih(ncell: int):
    bond = 2.76
    a = bond * math.sqrt(8.0 / 3.0)
    c = a * math.sqrt(8.0 / 3.0)
    b = a * math.sqrt(3.0)
    basis = (
        (1.0 / 3.0, 2.0 / 3.0, 0.0),
        (2.0 / 3.0, 1.0 / 3.0, 0.5),
        (1.0 / 3.0, 2.0 / 3.0, 3.0 / 8.0),
        (2.0 / 3.0, 1.0 / 3.0, 7.0 / 8.0),
    )
    lx, ly, lz = ncell * a, ncell * b, ncell * c
    xyz = []
    for ix in range(ncell):
        for iy in range(ncell):
            for iz in range(ncell):
                for shift in (0.0, 0.5):
                    for site in basis:
                        x = (shift * a + site[0] * a - site[1] * a / 2.0 + ix * a) % lx
                        y = (shift * b + site[1] * b / 2.0 + iy * b) % ly
                        z = (site[2] + iz) * c
                        xyz.extend((x, y, z))
    n = len(xyz) // 3
    arr = (ctypes.c_double * len(xyz))(*xyz)
    box = (ctypes.c_double * 9)(0, 0, 0, lx, ly, lz, 0, 0, 0)
    types = (ctypes.c_int * n)(*([1] * n))
    ids = (ctypes.c_int * n)(*range(1, n + 1))
    return n, arr, box, types, ids


def write_dump(path: Path, n: int, xyz, box) -> None:
    lines = [
        "ITEM: TIMESTEP",
        "0",
        "ITEM: NUMBER OF ATOMS",
        str(n),
        "ITEM: BOX BOUNDS pp pp pp",
        f"0 {box[3]}",
        f"0 {box[4]}",
        f"0 {box[5]}",
        "ITEM: ATOMS id type x y z",
    ]
    for i in range(n):
        lines.append(
            f"{i + 1} 1 {xyz[3 * i]:.8f} {xyz[3 * i + 1]:.8f} {xyz[3 * i + 2]:.8f}"
        )
    path.write_text("\n".join(lines) + "\n", encoding="ascii")


def median_of(samples: list[float]) -> float:
    return statistics.median(samples)


def session(lib, threads: int, samples: int = SAMPLES) -> dict[str, float]:
    if lib.seams_set_num_threads(threads) != 0:
        raise SystemExit(f"seams_set_num_threads({threads}) failed")
    n, xyz, box, types, ids = ice_ih(24)
    sys = lib.seams_system_create()
    if lib.seams_system_set_frame(sys, xyz, n, box, types, ids) != 0:
        raise SystemExit("set_frame failed")
    q = (ctypes.c_double * n)()
    bar = (ctypes.c_double * n)()
    labels = (ctypes.c_int * n)()
    counts = (ctypes.c_int * 9)()
    nb = ctypes.c_void_p()

    def knn():
        if lib.seams_neighbours_knn(sys, 4, 3.5, 1, 1, ctypes.byref(nb)) != 0:
            raise SystemExit("knn failed")
        lib.seams_neighbours_destroy(nb)

    def q6():
        if lib.seams_steinhardt(sys, 3.5, 1, 6, q, n, bar, n) != 0:
            raise SystemExit("q6 failed")

    def chill():
        if lib.seams_chill(sys, 3.5, 1, 0, 0, labels, n, counts) != 0:
            raise SystemExit("chill failed")

    with tempfile.TemporaryDirectory(prefix="seams-omp-") as folder:
        dump = Path(folder) / "ice.lammpstrj"
        write_dump(dump, n, xyz, box)
        encoded = os.fsencode(dump)

        def read():
            fresh = lib.seams_system_create()
            if lib.seams_system_read(fresh, encoded, 1, 1) != 0:
                raise SystemExit("read failed")
            lib.seams_system_destroy(fresh)

        out = {}
        for name, fn in (("knn", knn), ("q6", q6), ("chill", chill), ("read", read)):
            fn()
            ticks = []
            for _ in range(samples):
                t0 = time_ns()
                fn()
                ticks.append((time_ns() - t0) / 1e6)
            out[name] = median_of(ticks)
    lib.seams_system_destroy(sys)
    return out


def time_ns() -> int:
    return __import__("time").perf_counter_ns()


def main() -> None:
    lib = load_library(library_path())
    cores = os.cpu_count() or 1
    print(f"cores {cores}")
    print(f"samples {SAMPLES}")
    for threads in (1, 8):
        got = session(lib, threads)
        for name, ms in got.items():
            print(f"{name} threads {threads} median_ms {ms:.3f}")


if __name__ == "__main__":
    main()
