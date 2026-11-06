"""Worked example of the manuscript: closed-form counts against the library.

Builds ideal ice Ic and Ih (the 1,000-site frames of the head-to-head
campaign) and periodic square, pentagonal and hexagonal ice nanotubes, then
calls the C interface of seams-core (libseams_lang) for the ring census, the
per-atom cage labels and the prism blocks. Writes a LAMMPS dump of every
frame and one tab-separated table.

Usage: python3 worked_example.py LIBSEAMS_LANG OUTDIR
"""

import ctypes
import math
import pathlib
import sys

BOND_A = 2.75
PRISM_HEIGHT_A = 2.845
CUTOFF_A = 3.5
TUBE_LAYERS = 12
TUBE_BOX_A = 25.0
# The ring census counts the rings of MIN_RING to MAX_RING members, and the
# prism blocks are those of the square, pentagonal and hexagonal tubes.
MIN_RING = 3
MAX_RING = 6
PRISM_SIZES = (4, 5, 6)
# Atom type of every site of the built frames, the type the C interface keeps.
ATOM_TYPE = 1
# Nearest neighbours of the k-nearest graphs of seams_frame_cages, its default.
CAGE_NEIGHBOURS = 4
# The cage labels run without the ring completion of the seeded graph.
RING_COMPLETION = 0


def cubic_diamond(repetitions):
    lattice = 4.0 * BOND_A / math.sqrt(3.0)
    fcc = ((0.0, 0.0, 0.0), (0.5, 0.5, 0.0), (0.5, 0.0, 0.5), (0.0, 0.5, 0.5))
    basis = fcc + tuple((x + 0.25, y + 0.25, z + 0.25) for x, y, z in fcc)
    side = repetitions * lattice
    points = [
        (((i + x) * lattice) % side, ((j + y) * lattice) % side,
         ((k + z) * lattice) % side)
        for i in range(repetitions)
        for j in range(repetitions)
        for k in range(repetitions)
        for x, y, z in basis
    ]
    return points, (side, side, side)


def lonsdaleite(repetitions):
    lattice = BOND_A / ((3.0 / 8.0) * math.sqrt(8.0 / 3.0))
    height = lattice * math.sqrt(8.0 / 3.0)
    hexagonal = (
        (1.0 / 3.0, 2.0 / 3.0, 0.0),
        (2.0 / 3.0, 1.0 / 3.0, 0.5),
        (1.0 / 3.0, 2.0 / 3.0, 3.0 / 8.0),
        (2.0 / 3.0, 1.0 / 3.0, 7.0 / 8.0),
    )
    x_length = lattice
    y_length = lattice * math.sqrt(3.0)
    basis = []
    for hx, hy, hz in hexagonal:
        x = (hx - 0.5 * hy) * lattice
        y = hy * math.sqrt(3.0) * lattice / 2.0
        for shift_x, shift_y in ((0.0, 0.0), (x_length / 2.0, y_length / 2.0)):
            basis.append(((x + shift_x) % x_length, (y + shift_y) % y_length,
                          hz * height))
    box = (repetitions * x_length, repetitions * y_length, repetitions * height)
    points = [
        ((i * x_length + x) % box[0], (j * y_length + y) % box[1],
         (k * height + z) % box[2])
        for i in range(repetitions)
        for j in range(repetitions)
        for k in range(repetitions)
        for x, y, z in basis
    ]
    return points, box


def nanotube(sides):
    radius = BOND_A / (2.0 * math.sin(math.pi / sides))
    centre = TUBE_BOX_A / 2.0
    points = []
    for layer in range(TUBE_LAYERS):
        for vertex in range(sides):
            angle = 2.0 * math.pi * vertex / sides
            points.append((centre + radius * math.cos(angle),
                           centre + radius * math.sin(angle),
                           (layer + 0.5) * PRISM_HEIGHT_A))
    return points, (TUBE_BOX_A, TUBE_BOX_A, TUBE_LAYERS * PRISM_HEIGHT_A)


def write_lammps(path, points, box):
    with open(path, "w") as handle:
        handle.write("ITEM: TIMESTEP\n0\nITEM: NUMBER OF ATOMS\n")
        handle.write(f"{len(points)}\nITEM: BOX BOUNDS pp pp pp\n")
        for length in box:
            handle.write(f"0.0 {length:.12f}\n")
        handle.write("ITEM: ATOMS id type x y z\n")
        for index, (x, y, z) in enumerate(points, 1):
            handle.write(f"{index} 1 {x:.12f} {y:.12f} {z:.12f}\n")


class Library:
    def __init__(self, path):
        lib = ctypes.CDLL(path)
        lib.seams_last_error.restype = ctypes.c_char_p
        lib.seams_lang_version.restype = ctypes.c_char_p
        lib.seams_frame_from_arrays.restype = ctypes.c_void_p
        lib.seams_frame_from_arrays.argtypes = [
            ctypes.POINTER(ctypes.c_double), ctypes.c_int,
            ctypes.POINTER(ctypes.c_double), ctypes.POINTER(ctypes.c_int),
            ctypes.POINTER(ctypes.c_int)]
        lib.seams_frame_free.argtypes = [ctypes.c_void_p]
        lib.seams_frame_rings.argtypes = [
            ctypes.c_void_p, ctypes.c_double, ctypes.c_int, ctypes.c_char_p,
            ctypes.c_int, ctypes.POINTER(ctypes.c_int), ctypes.c_int]
        lib.seams_frame_cages.argtypes = [
            ctypes.c_void_p, ctypes.c_double, ctypes.c_int, ctypes.c_int,
            ctypes.c_char_p, ctypes.c_int, ctypes.POINTER(ctypes.c_int),
            ctypes.POINTER(ctypes.c_int), ctypes.POINTER(ctypes.c_int),
            ctypes.POINTER(ctypes.c_int), ctypes.POINTER(ctypes.c_int)]
        lib.seams_frame_prisms.argtypes = [
            ctypes.c_void_p, ctypes.c_double, ctypes.c_int, ctypes.c_int,
            ctypes.c_char_p, ctypes.c_void_p, ctypes.POINTER(ctypes.c_int)]
        self.lib = lib

    def check(self, status, what):
        if status != 0:
            raise RuntimeError(f"{what}: {self.lib.seams_last_error().decode()}")

    def frame(self, points, box):
        n = len(points)
        xyz = (ctypes.c_double * (3 * n))(*[c for p in points for c in p])
        box9 = (ctypes.c_double * 9)(0.0, 0.0, 0.0, box[0], box[1], box[2],
                                     0.0, 0.0, 0.0)
        types = (ctypes.c_int * n)(*([ATOM_TYPE] * n))
        ids = (ctypes.c_int * n)(*range(1, n + 1))
        handle = self.lib.seams_frame_from_arrays(xyz, n, box9, types, ids)
        if not handle:
            raise RuntimeError(self.lib.seams_last_error().decode())
        return handle

    def rings(self, frame, max_depth):
        counts = (ctypes.c_int * (max_depth + 1))()
        self.check(self.lib.seams_frame_rings(frame, CUTOFF_A, ATOM_TYPE,
                                              b"cutoff", max_depth, counts,
                                              max_depth + 1),
                   "rings")
        return list(counts)

    def cages(self, frame, n):
        per_atom = (ctypes.c_int * n)()
        out = [ctypes.c_int() for _ in range(4)]
        self.check(self.lib.seams_frame_cages(
            frame, CUTOFF_A, ATOM_TYPE, CAGE_NEIGHBOURS, b"cutoff",
            RING_COMPLETION, per_atom,
            *[ctypes.byref(o) for o in out]), "cages")
        return [o.value for o in out]

    def prisms(self, frame, max_depth):
        rows = (ctypes.c_int * (3 * (max_depth - 2)))()
        n_rows = ctypes.c_int(max_depth - 2)
        self.check(self.lib.seams_frame_prisms(
            frame, CUTOFF_A, ATOM_TYPE, max_depth, b"",
            ctypes.cast(rows, ctypes.c_void_p), ctypes.byref(n_rows)),
            "prisms")
        return {rows[3 * i]: (rows[3 * i + 1], rows[3 * i + 2])
                for i in range(n_rows.value)}


def main():
    lib = Library(sys.argv[1])
    out = pathlib.Path(sys.argv[2])
    out.mkdir(parents=True, exist_ok=True)
    frames = {
        "ice-Ic": cubic_diamond(5),
        "ice-Ih": lonsdaleite(5),
        "tube-4": nanotube(4),
        "tube-5": nanotube(5),
        "tube-6": nanotube(6),
    }
    rows = ["frame\tsites\tbox_max_A\trings_3_to_6\tcage_hex_cubic_water_both"
            "\tprisms_4_5_6\theight_percent_4_5_6"]
    for name, (points, box) in frames.items():
        write_lammps(out / f"{name}.lammpstrj", points, box)
        frame = lib.frame(points, box)
        try:
            census = lib.rings(frame, MAX_RING)[MIN_RING:MAX_RING + 1]
            cages = lib.cages(frame, len(points))
            prisms = lib.prisms(frame, MAX_RING)
        finally:
            lib.lib.seams_frame_free(frame)
        n_max = max(box) / PRISM_HEIGHT_A
        blocks = [sum(prisms.get(size, (0, 0))) for size in PRISM_SIZES]
        heights = [f"{100.0 * b / n_max:.3f}" for b in blocks]
        rows.append("\t".join([
            name, str(len(points)), f"{max(box):.6f}",
            ",".join(map(str, census)), ",".join(map(str, cages)),
            ",".join(map(str, blocks)), ",".join(heights)]))
    table = "\n".join(rows) + "\n"
    (out / "worked_example.tsv").write_text(table)
    sys.stdout.write(f"libseams_lang {lib.lib.seams_lang_version().decode()}\n")
    sys.stdout.write(table)


if __name__ == "__main__":
    main()
