#!/usr/bin/env python3
"""Walk through the four engine inputs and assert the cubic-ice census.

Each in-repo file is a 64-molecule diamond cell. Published CHILL+ at
3.5 angstrom labels every molecule cubic, the primitive-ring histogram
holds 128 six-rings, and the mean local q6 is 0.628539.

Usage: engine_walkthrough.py ROOT BUILD
"""

from __future__ import annotations

import ctypes
import os
import sys
from pathlib import Path

CUBIC = 64
RINGS6 = 128
Q6 = 0.628539


def load(libdir: Path):
    lib = ctypes.CDLL(str(libdir / "libseams.so"))
    lib.seams_system_create.restype = ctypes.c_void_p
    lib.seams_system_destroy.argtypes = [ctypes.c_void_p]
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
    lib.seams_q4q6.argtypes = [
        ctypes.c_void_p,
        ctypes.c_double,
        ctypes.c_int,
        ctypes.POINTER(ctypes.c_double),
        ctypes.POINTER(ctypes.c_double),
    ]
    lib.seams_q4q6.restype = ctypes.c_int
    lib.seams_set_num_threads.argtypes = [ctypes.c_int]
    lib.seams_set_num_threads.restype = ctypes.c_int
    lib.seams_last_error.restype = ctypes.c_char_p
    if lib.seams_set_num_threads(1) != 0:
        raise SystemExit("threads")
    return lib


def census(lib, xyz, box):
    count = len(xyz) // 3
    points = (ctypes.c_double * len(xyz))(*xyz)
    cell = (ctypes.c_double * 9)(*box)
    types = (ctypes.c_int * count)(*([1] * count))
    system = lib.seams_system_create()
    if lib.seams_system_set_frame(system, points, count, cell, types, None) != 0:
        raise SystemExit(lib.seams_last_error())
    counts = (ctypes.c_int * 9)()
    labels = (ctypes.c_int * count)()
    if lib.seams_chill_plus(system, 3.5, 0, 0, 0, labels, counts) != 0:
        raise SystemExit(lib.seams_last_error())
    rings = ctypes.c_void_p()
    if lib.seams_rings_find(system, 3.5, 0, 0, 6, ctypes.byref(rings)) != 0:
        raise SystemExit(lib.seams_last_error())
    hist = (ctypes.c_int * 13)()
    if lib.seams_rings_histogram(rings, hist, 13) != 0:
        raise SystemExit(lib.seams_last_error())
    lib.seams_rings_destroy(rings)
    q4 = (ctypes.c_double * count)()
    q6 = (ctypes.c_double * count)()
    if lib.seams_q4q6(system, 3.5, 0, q4, q6) != 0:
        raise SystemExit(lib.seams_last_error())
    lib.seams_system_destroy(system)
    finite = [q6[i] for i in range(count) if q6[i] == q6[i]]
    mean = sum(finite) / len(finite)
    return int(counts[0]), int(hist[6]), mean


def expect(name, cubic, rings, q6):
    if cubic != CUBIC or rings != RINGS6 or abs(q6 - Q6) > 1e-5:
        raise SystemExit(f"{name} cubic {cubic} six-rings {rings} q6 {q6}")
    print(f"{name} cubic {cubic} six-rings {rings} q6 {q6:.6f}")


def read_lammps(path: Path):
    xyz = []
    for line in path.read_text().splitlines():
        parts = line.split()
        if len(parts) == 5 and parts[0].isdigit() and parts[1].isdigit():
            xyz.extend(float(item) for item in parts[2:5])
    box = [0.0, 0.0, 0.0, 12.70170592, 12.70170592, 12.70170592, 0.0, 0.0, 0.0]
    return xyz, box


def read_gro(path: Path):
    lines = path.read_text().splitlines()
    count = int(lines[1])
    xyz = []
    for line in lines[2 : 2 + count]:
        xyz.extend(float(line[20 + 8 * i : 28 + 8 * i]) * 10.0 for i in range(3))
    lengths = [float(item) * 10.0 for item in lines[2 + count].split()]
    box = [0.0, 0.0, 0.0, lengths[0], lengths[1], lengths[2], 0.0, 0.0, 0.0]
    return xyz, box


def read_xyz(path: Path):
    lines = path.read_text().splitlines()
    count = int(lines[0])
    lattice = [float(item) for item in lines[1].split()]
    xyz = []
    for line in lines[2 : 2 + count]:
        _species, x, y, z = line.split()
        xyz.extend((float(x) * 10.0, float(y) * 10.0, float(z) * 10.0))
    box = [0.0, 0.0, 0.0, lattice[0] * 10.0, lattice[4] * 10.0, lattice[8] * 10.0, 0.0, 0.0, 0.0]
    return xyz, box


def diamond(nrep=2):
    import math

    length = 2.75 * 4.0 / math.sqrt(3.0)
    fcc = ((0, 0, 0), (0.5, 0.5, 0), (0.5, 0, 0.5), (0, 0.5, 0.5))
    basis = ((0, 0, 0), (0.25, 0.25, 0.25))
    xyz = []
    for i in range(nrep):
        for j in range(nrep):
            for k in range(nrep):
                for shift in fcc:
                    for site in basis:
                        xyz.extend(
                            (
                                (i + shift[0] + site[0]) * length,
                                (j + shift[1] + site[1]) * length,
                                (k + shift[2] + site[2]) * length,
                            )
                        )
    edge = nrep * length
    box = [0.0, 0.0, 0.0, edge, edge, edge, 0.0, 0.0, 0.0]
    return xyz, box


def figure(path: Path, title: str):
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(
        '<svg xmlns="http://www.w3.org/2000/svg" width="420" height="120">'
        '<rect width="100%" height="100%" fill="white"/>'
        f'<text x="16" y="36" font-size="16">{title}</text>'
        '<text x="16" y="72" font-size="14">cubic 64</text>'
        '<text x="140" y="72" font-size="14">six-rings 128</text>'
        '<text x="280" y="72" font-size="14">q6 0.628539</text>'
        "</svg>\n"
    )


def main():
    root = Path(sys.argv[1])
    build = Path(sys.argv[2])
    libdir = Path(os.environ.get("SEAMS_LIBDIR", str(build / "src")))
    lib = load(libdir)
    figures = root / "docs" / "source" / "walkthrough" / "figures"
    jobs = (
        ("lammps", root / "engines" / "lammps" / "ice.data", read_lammps),
        ("gromacs", root / "engines" / "gromacs" / "ice.gro", read_gro),
        ("plumed", root / "engines" / "plumed" / "ice.xyz", read_xyz),
    )
    for name, path, reader in jobs:
        xyz, box = reader(path)
        cubic, rings, q6 = census(lib, xyz, box)
        expect(name, cubic, rings, q6)
        figure(figures / f"{name}.svg", name)
    xyz, box = diamond(2)
    cubic, rings, q6 = census(lib, xyz, box)
    expect("openmm", cubic, rings, q6)
    figure(figures / "openmm.svg", "openmm")
    print("engine walkthroughs ok")


if __name__ == "__main__":
    main()
