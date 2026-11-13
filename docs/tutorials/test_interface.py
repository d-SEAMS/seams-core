#!/usr/bin/env python3
"""CHILL+ layer profile and the ice front on one interface frame."""

import ctypes
import os
import re
import subprocess
import sys
from pathlib import Path

WIDTH = 3.7
# Crystal labels: cubic, hexagonal, reCubic, reHex.
ICE = {0, 1, 7, 8}
# ice molecules in each 3.7 A slab along z, then the atoms in that slab.
LAYERS = (
    (30, 69),
    (65, 69),
    (39, 39),
    (72, 72),
    (69, 69),
    (39, 41),
    (34, 71),
    (4, 62),
    (0, 57),
    (0, 66),
    (0, 61),
    (0, 58),
    (1, 65),
    (0, 62),
    (0, 3),
)


def load_dump(path):
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


def labels_of(libdir, box, atoms):
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
    count = len(atoms)
    xyz = (ctypes.c_double * (3 * count))()
    for index, atom in enumerate(atoms):
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
    return [labels[index] for index in range(count)], [counts[index] for index in range(9)], oz


def front_z(layers):
    """z where the ice fraction crosses one half, between two slabs."""
    for index in range(1, len(layers)):
        ice0, total0 = layers[index - 1]
        ice1, total1 = layers[index]
        if total0 == 0 or total1 == 0:
            continue
        left = ice0 / total0
        right = ice1 / total1
        if left >= 0.5 > right:
            span = left - right
            frac = (left - 0.5) / span
            return (index - 1 + frac) * WIDTH
    raise SystemExit("no front")


def main():
    root = Path(sys.argv[1])
    seams = sys.argv[2]
    libdir = sys.argv[3]
    frame = root / "input" / "traj" / "interface_f0.lammpstrj"
    env = os.environ.copy()
    env["LD_LIBRARY_PATH"] = str(libdir) + os.pathsep + env.get("LD_LIBRARY_PATH", "")
    text = subprocess.check_output(
        [seams, "chill-plus", str(frame), "--type", "1", "--cutoff", "3.5", "--layers", "--axis", "2"],
        env=env,
        text=True,
    )
    if "hexagonal 352" not in text or "water 434" not in text:
        raise SystemExit(text)
    stack = re.search(r"stack (\S+)", text)
    if stack is None or stack.group(1) != "HHHHHHHH....C.":
        raise SystemExit(text)
    box, atoms = load_dump(frame)
    labels, counts, oz = labels_of(libdir, box, atoms)
    if counts[0] != 1 or counts[1] != 352 or counts[2] != 434:
        raise SystemExit(f"counts {counts}")
    layers = []
    for ice, total in LAYERS:
        layers.append([0, 0])
    for index, atom in enumerate(atoms):
        slot = int((float(atom["z"]) - oz) / WIDTH)
        if slot < 0 or slot >= len(layers):
            raise SystemExit(f"z bin {slot}")
        layers[slot][1] += 1
        if labels[index] in ICE:
            layers[slot][0] += 1
    if [tuple(row) for row in layers] != list(LAYERS):
        raise SystemExit(f"layers {layers}")
    position = front_z(layers)
    if not 22.0 < position < 22.1:
        raise SystemExit(f"front {position}")
    print(f"interface ok front {position:.3f} stack {stack.group(1)}")


if __name__ == "__main__":
    main()
