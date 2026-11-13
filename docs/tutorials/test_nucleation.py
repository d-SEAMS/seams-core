#!/usr/bin/env python3
"""Largest ice cluster and the cubic and hexagonal fractions.

The file is the first three frames of the 148-frame mW nucleation
trajectory. Crystal atoms are cubic, hexagonal, reCubic and reHex.
"""

import ctypes
import sys
from pathlib import Path

# cubic, hexagonal, largest crystal cluster
EXPECT = ((47, 62, 59), (50, 60, 64), (44, 69, 70))
CRYSTAL = {0, 1, 7, 8}


def frames_of(path):
    text = Path(path).read_text()
    parts = [part for part in text.split("ITEM: TIMESTEP\n") if part.strip()]
    return ["ITEM: TIMESTEP\n" + part for part in parts]


def one_frame(body):
    lines = body.splitlines()
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


def cluster(libdir, box, atoms):
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
    lib.seams_neighbours_build.restype = ctypes.c_int
    lib.seams_neighbours_build.argtypes = [
        ctypes.c_void_p,
        ctypes.c_double,
        ctypes.c_int,
        ctypes.c_int,
        ctypes.POINTER(ctypes.c_void_p),
    ]
    lib.seams_neighbours_csr.restype = ctypes.c_int
    lib.seams_neighbours_csr.argtypes = [
        ctypes.c_void_p,
        ctypes.POINTER(ctypes.POINTER(ctypes.c_int)),
        ctypes.POINTER(ctypes.c_int),
        ctypes.POINTER(ctypes.POINTER(ctypes.c_int)),
        ctypes.POINTER(ctypes.c_int),
    ]
    lib.seams_neighbours_destroy.argtypes = [ctypes.c_void_p]
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
    neighbours = ctypes.c_void_p()
    if lib.seams_neighbours_build(handle, 3.5, 0, 0, ctypes.byref(neighbours)) != 0:
        raise SystemExit("neighbours")
    offsets = ctypes.POINTER(ctypes.c_int)()
    cols = ctypes.POINTER(ctypes.c_int)()
    n_off = ctypes.c_int()
    nnz = ctypes.c_int()
    if lib.seams_neighbours_csr(neighbours, ctypes.byref(offsets), ctypes.byref(n_off), ctypes.byref(cols), ctypes.byref(nnz)) != 0:
        raise SystemExit("csr")
    parent = list(range(count))

    def find(atom):
        while parent[atom] != atom:
            parent[atom] = parent[parent[atom]]
            atom = parent[atom]
        return atom

    for index in range(count):
        if labels[index] not in CRYSTAL:
            continue
        for edge in range(offsets[index], offsets[index + 1]):
            other = cols[edge]
            if labels[other] in CRYSTAL:
                left, right = find(index), find(other)
                if left != right:
                    parent[right] = left
    sizes = {}
    for index in range(count):
        if labels[index] in CRYSTAL:
            root = find(index)
            sizes[root] = sizes.get(root, 0) + 1
    lib.seams_neighbours_destroy(neighbours)
    lib.seams_system_destroy(handle)
    largest = max(sizes.values()) if sizes else 0
    return counts[0], counts[1], largest, count


def main():
    root = Path(sys.argv[1])
    libdir = sys.argv[2]
    path = root / "input" / "traj" / "nucleation_three.lammpstrj"
    bodies = frames_of(path)
    if len(bodies) != 3:
        raise SystemExit(f"frames {len(bodies)}")
    for number, body in enumerate(bodies):
        cubic, hexagonal, largest, count = cluster(libdir, *one_frame(body))
        want_c, want_h, want_largest = EXPECT[number]
        if (cubic, hexagonal, largest) != (want_c, want_h, want_largest):
            raise SystemExit(f"frame {number + 1} {(cubic, hexagonal, largest)}")
        if count != 4096:
            raise SystemExit(f"nop {count}")
        if cubic / count != want_c / 4096 or hexagonal / count != want_h / 4096:
            raise SystemExit("fraction")
    print("nucleation ok", EXPECT)


if __name__ == "__main__":
    main()
