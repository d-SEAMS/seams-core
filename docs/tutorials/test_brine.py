#!/usr/bin/env python3
"""Ions per layer on a growing ice front.

Water is the interface frame. Nine ions were added: one chloride in a
fully crystalline slab and eight ions in the liquid.
"""

import ctypes
import sys
from pathlib import Path

WIDTH = 3.7
ICE = {0, 1, 7, 8}
# bin -> ion count. Bin 3 is ice. The others are liquid.
IONS = {3: 1, 9: 2, 10: 2, 12: 2, 13: 2}


def load(path):
    lines = Path(path).read_text().splitlines()
    box = None
    atoms = []
    index = 0
    while index < len(lines):
        if lines[index].startswith("ITEM: BOX"):
            rows = []
            for _ in range(3):
                index += 1
                rows.append([float(value) for value in lines[index].split()])
            box = (
                rows[0][0],
                rows[1][0],
                rows[2][0],
                rows[0][1] - rows[0][0],
                rows[1][1] - rows[1][0],
                rows[2][1] - rows[2][0],
            )
        elif lines[index].startswith("ITEM: ATOMS"):
            cols = lines[index].split()[2:]
            index += 1
            while index < len(lines) and not lines[index].startswith("ITEM:"):
                atoms.append(dict(zip(cols, lines[index].split())))
                index += 1
            break
        index += 1
    return box, atoms


def water_labels(libdir, box, water):
    lib = ctypes.CDLL(str(Path(libdir) / "libseams.so"))
    lib.seams_system_create.restype = ctypes.c_void_p
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
    lib.seams_system_destroy.argtypes = [ctypes.c_void_p]
    count = len(water)
    xyz = (ctypes.c_double * (3 * count))()
    for index, atom in enumerate(water):
        xyz[3 * index] = float(atom["x"])
        xyz[3 * index + 1] = float(atom["y"])
        xyz[3 * index + 2] = float(atom["z"])
    ox, oy, oz, lx, ly, lz = box
    handle = lib.seams_system_create()
    if lib.seams_system_set_ortho(handle, xyz, count, ox, oy, oz, lx, ly, lz, None, None) != 0:
        raise SystemExit("frame")
    labels = (ctypes.c_int * count)()
    counts = (ctypes.c_int * 9)()
    if lib.seams_chill_plus(handle, 3.5, 0, 0, 0, labels, counts) != 0:
        raise SystemExit("chill+")
    lib.seams_system_destroy(handle)
    return [labels[index] for index in range(count)], counts[1], oz


def main():
    root = Path(sys.argv[1])
    libdir = sys.argv[2]
    box, atoms = load(root / "input" / "traj" / "nacl_front.lammpstrj")
    water = [atom for atom in atoms if atom["type"] == "1"]
    ions = [atom for atom in atoms if atom["type"] in {"3", "4"}]
    if len(water) != 864 or len(ions) != 9:
        raise SystemExit(f"atoms water {len(water)} ions {len(ions)}")
    labels, hexagonal, oz = water_labels(libdir, box, water)
    if hexagonal != 352:
        raise SystemExit(f"hexagonal {hexagonal}")
    ice_fraction = {}
    for index, atom in enumerate(water):
        slot = int((float(atom["z"]) - oz) / WIDTH)
        ice, total = ice_fraction.get(slot, (0, 0))
        ice_fraction[slot] = (ice + (1 if labels[index] in ICE else 0), total + 1)
    got = {}
    ice_ions = 0
    liquid_ions = 0
    for atom in ions:
        slot = int((float(atom["z"]) - oz) / WIDTH)
        got[slot] = got.get(slot, 0) + 1
        ice, total = ice_fraction[slot]
        if ice / total >= 0.5:
            ice_ions += 1
        else:
            liquid_ions += 1
    if got != IONS or ice_ions != 1 or liquid_ions != 8:
        raise SystemExit(f"ions {got} ice {ice_ions} liquid {liquid_ions}")
    print("brine ok ice_ions 1 liquid_ions 8")


if __name__ == "__main__":
    main()
