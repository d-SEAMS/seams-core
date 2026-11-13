#!/usr/bin/env python3
# Copyright (c) d-SEAMS developers
"""Classify the bundled eight-atom cell the way the engine plugins do.

LAMMPS ``run_hook`` reads positions in angstroms and calls
``seams_chill_plus``. The OpenMM plugin stores nanometres and multiplies
by 10 before the same call. This script does both conversions on
``engines/tutorial/ice.data`` and requires every label to be cubic.
``meson test engines_tutorial`` runs it against the build's libseams.
"""

import ctypes
import pathlib
import sys

ROOT = pathlib.Path(__file__).resolve().parents[1]
DATA = pathlib.Path(__file__).resolve().parent / "tutorial" / "ice.data"
CUTOFF = 3.5
CUBIC = 0


def read_data(path):
    """LAMMPS atomic data: id type x y z, and the orthorhombic edge."""
    text = path.read_text().splitlines()
    natoms = None
    length = None
    atoms = []
    mode = None
    for line in text:
        parts = line.split()
        if not parts:
            continue
        if len(parts) >= 2 and parts[1] == "atoms":
            natoms = int(parts[0])
            continue
        if len(parts) >= 4 and parts[2:] == ["xlo", "xhi"]:
            length = float(parts[1]) - float(parts[0])
            continue
        if parts[0] == "Atoms":
            mode = "atoms"
            continue
        if mode == "atoms" and len(parts) >= 5 and parts[0].isdigit():
            atoms.append((int(parts[0]), float(parts[2]), float(parts[3]), float(parts[4])))
    atoms.sort()
    if natoms is None or length is None or len(atoms) != natoms:
        raise SystemExit("ice.data did not parse")
    xyz = []
    for _ident, x, y, z in atoms:
        xyz.extend((x, y, z))
    return natoms, length, xyz


def load_library(build):
    lib = pathlib.Path(build) / "src" / "libseams.so"
    if not lib.is_file():
        raise SystemExit("libseams.so is not at {}".format(lib))
    handle = ctypes.CDLL(str(lib))
    handle.seams_system_create.restype = ctypes.c_void_p
    handle.seams_system_set_ortho.argtypes = [
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
    handle.seams_chill_plus.argtypes = [
        ctypes.c_void_p,
        ctypes.c_double,
        ctypes.c_int,
        ctypes.c_int,
        ctypes.c_int,
        ctypes.POINTER(ctypes.c_int),
        ctypes.POINTER(ctypes.c_int),
    ]
    handle.seams_system_destroy.argtypes = [ctypes.c_void_p]
    return handle


def chill_plus(lib, xyz, length):
    n = len(xyz) // 3
    coords = (ctypes.c_double * len(xyz))(*xyz)
    labels = (ctypes.c_int * n)()
    counts = (ctypes.c_int * 9)()
    system = lib.seams_system_create()
    if not system:
        raise SystemExit("seams_system_create failed")
    code = lib.seams_system_set_ortho(
        system, coords, n, 0.0, 0.0, 0.0, length, length, length, None, None
    )
    if code != 0:
        lib.seams_system_destroy(system)
        raise SystemExit("seams_system_set_ortho returned {}".format(code))
    code = lib.seams_chill_plus(system, CUTOFF, 0, 0, 0, labels, counts)
    lib.seams_system_destroy(system)
    if code != 0:
        raise SystemExit("seams_chill_plus returned {}".format(code))
    return [labels[i] for i in range(n)], counts[CUBIC]


def openmm_angstroms(xyz):
    """The plugin reads OpenMM nanometres and multiplies by 10."""
    nanometres = [value / 10.0 for value in xyz]
    return [value * 10.0 for value in nanometres]


def main():
    build = sys.argv[1] if len(sys.argv) > 1 else ROOT / "build"
    natoms, length, xyz = read_data(DATA)
    if natoms != 8:
        raise SystemExit("the bundled cell has {} atoms".format(natoms))
    lib = load_library(build)
    labels, cubic = chill_plus(lib, xyz, length)
    if cubic != natoms or any(label != CUBIC for label in labels):
        raise SystemExit("LAMMPS path labels {}".format(labels))
    scaled = openmm_angstroms(xyz)
    again, cubic_nm = chill_plus(lib, scaled, length)
    if again != labels or cubic_nm != natoms:
        raise SystemExit("OpenMM scaling changed the labels")
    print("engines tutorial ok atoms {} cubic {}".format(natoms, cubic))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
