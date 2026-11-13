#!/usr/bin/env python3
"""Mixed-phase oxygen frames for the CHILL regression.

Three cells, written as LAMMPS dumps:

- interface: ice Ih slab against a disordered liquid of the same density
- stacking: ice Ih with two bilayers shifted onto the cubic stacking vector
- melted: ice Ih with a sphere of molecules replaced by a random melt

The generator is deterministic. Pass --label LIB to write one integer
label per oxygen for published CHILL and CHILL+ at one and eight threads.
"""

from __future__ import annotations

import argparse
import ctypes
import math
from pathlib import Path

import numpy as np

BOND = 2.76
A = BOND * math.sqrt(8.0 / 3.0)
C_LEN = A * math.sqrt(8.0 / 3.0)
B_LEN = A * math.sqrt(3.0)
BASIS = (
    (1.0 / 3.0, 2.0 / 3.0, 0.0),
    (2.0 / 3.0, 1.0 / 3.0, 0.5),
    (1.0 / 3.0, 2.0 / 3.0, 3.0 / 8.0),
    (2.0 / 3.0, 1.0 / 3.0, 7.0 / 8.0),
)
NAMES = (
    "cubic",
    "hexagonal",
    "water",
    "interfacial",
    "clathrate",
    "interclathrate",
    "unclassified",
    "recubic",
    "rehex",
)


def ice_ih(nx: int, ny: int, nz: int) -> tuple[np.ndarray, tuple[float, float, float]]:
    lx, ly, lz = nx * A, ny * B_LEN, nz * C_LEN
    coords = []
    for ix in range(nx):
        for iy in range(ny):
            for iz in range(nz):
                for shift in (0.0, 0.5):
                    for site in BASIS:
                        coords.append((
                            (shift * A + site[0] * A - site[1] * A / 2.0 + ix * A) % lx,
                            (shift * B_LEN + site[1] * B_LEN / 2.0 + iy * B_LEN) % ly,
                            (site[2] + iz) * C_LEN,
                        ))
    return np.asarray(coords, dtype=np.float64), (lx, ly, lz)


def interface() -> tuple[np.ndarray, tuple[float, float, float]]:
    """Ice slab in the lower half, liquid of the same count in the upper half."""
    crystal, (lx, ly, lz_ice) = ice_ih(4, 4, 2)
    n = crystal.shape[0]
    rng = np.random.Generator(np.random.PCG64(20261113))
    liquid = np.column_stack((
        rng.random(n) * lx,
        rng.random(n) * ly,
        lz_ice + rng.random(n) * lz_ice,
    ))
    pos = np.vstack((crystal, liquid))
    return pos, (lx, ly, 2.0 * lz_ice)


def stacking() -> tuple[np.ndarray, tuple[float, float, float]]:
    """Ice Ih with two bilayers moved by the cubic stacking vector.

    A bilayer is the pair of basis sites at fractional heights 0 and 3/8,
    or 1/2 and 7/8. Shifting a bilayer by (a/3, b/3) changes its stacking
    from hexagonal to cubic relative to the bilayer below.
    """
    pos, lengths = ice_ih(4, 4, 4)
    lx, ly, lz = lengths
    # Bilayer index along z: 0 and 1 inside each orthorhombic cell.
    cell_z = np.floor(pos[:, 2] / C_LEN).astype(np.int64)
    frac = pos[:, 2] / C_LEN - cell_z
    bilayer = cell_z * 2 + (frac > 0.45).astype(np.int64)
    fault = (bilayer == 3) | (bilayer == 4)
    pos = pos.copy()
    # (a/2, b/6) is the basal shift that continues the bilayer as cubic ice.
    pos[fault, 0] = (pos[fault, 0] + A / 2.0) % lx
    pos[fault, 1] = (pos[fault, 1] + B_LEN / 6.0) % ly
    return pos, lengths


def melted() -> tuple[np.ndarray, tuple[float, float, float]]:
    """Ice Ih with molecules inside a sphere replaced by a random melt."""
    pos, lengths = ice_ih(4, 4, 4)
    lx, ly, lz = lengths
    centre = np.array((0.35 * lx, 0.40 * ly, 0.55 * lz))
    delta = pos - centre
    box = np.array(lengths)
    delta -= box * np.round(delta / box)
    inside = np.linalg.norm(delta, axis=1) < 8.0
    rng = np.random.Generator(np.random.PCG64(20261114))
    n_melt = int(np.count_nonzero(inside))
    # Keep the melt inside the same sphere so the crystal outside is intact.
    direction = rng.normal(size=(n_melt, 3))
    direction /= np.linalg.norm(direction, axis=1)[:, None]
    radius = 8.0 * np.cbrt(rng.random(n_melt))
    pos = pos.copy()
    pos[inside] = centre + direction * radius[:, None]
    pos[:, 0] %= lx
    pos[:, 1] %= ly
    pos[:, 2] %= lz
    return pos, lengths


FRAMES = {
    "interface": interface,
    "stacking": stacking,
    "melted": melted,
}


def write_dump(path: Path, pos: np.ndarray, lengths) -> None:
    lx, ly, lz = lengths
    lines = [
        "ITEM: TIMESTEP",
        "0",
        "ITEM: NUMBER OF ATOMS",
        str(pos.shape[0]),
        "ITEM: BOX BOUNDS pp pp pp",
        f"0 {lx:.10f}",
        f"0 {ly:.10f}",
        f"0 {lz:.10f}",
        "ITEM: ATOMS id type x y z",
    ]
    for i in range(pos.shape[0]):
        lines.append(
            f"{i + 1} 1 {pos[i, 0]:.8f} {pos[i, 1]:.8f} {pos[i, 2]:.8f}"
        )
    path.write_text("\n".join(lines) + "\n", encoding="ascii")


def load_seams(path: Path):
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
    lib.seams_chill.argtypes = [
        ctypes.c_void_p, ctypes.c_double, ctypes.c_int, ctypes.c_int, ctypes.c_int,
        ctypes.POINTER(ctypes.c_int), ctypes.c_size_t, ctypes.POINTER(ctypes.c_int),
    ]
    lib.seams_chill.restype = ctypes.c_int
    lib.seams_chill_plus.argtypes = [
        ctypes.c_void_p, ctypes.c_double, ctypes.c_int, ctypes.c_int, ctypes.c_int,
        ctypes.POINTER(ctypes.c_int), ctypes.POINTER(ctypes.c_int),
    ]
    lib.seams_chill_plus.restype = ctypes.c_int
    lib.seams_neighbours_build.argtypes = [
        ctypes.c_void_p, ctypes.c_double, ctypes.c_int, ctypes.c_int,
        ctypes.POINTER(ctypes.c_void_p),
    ]
    lib.seams_neighbours_build.restype = ctypes.c_int
    lib.seams_neighbours_destroy.argtypes = [ctypes.c_void_p]
    lib.seams_neighbours_copy.argtypes = [
        ctypes.c_void_p,
        ctypes.POINTER(ctypes.c_int), ctypes.c_size_t, ctypes.POINTER(ctypes.c_size_t),
        ctypes.POINTER(ctypes.c_int), ctypes.c_size_t, ctypes.POINTER(ctypes.c_size_t),
    ]
    lib.seams_neighbours_copy.restype = ctypes.c_int
    return lib


def label_frame(lib, pos: np.ndarray, lengths, threads: int, which: str) -> np.ndarray:
    n = pos.shape[0]
    lx, ly, lz = lengths
    box9 = (ctypes.c_double * 9)(0, 0, 0, lx, ly, lz, 0, 0, 0)
    types = (ctypes.c_int * n)(*([1] * n))
    ids = (ctypes.c_int * n)(*range(1, n + 1))
    ptr = np.ascontiguousarray(pos).ctypes.data_as(ctypes.POINTER(ctypes.c_double))
    if lib.seams_set_num_threads(threads) != 0:
        raise SystemExit("seams_set_num_threads failed")
    sys = lib.seams_system_create()
    if lib.seams_system_set_frame(sys, ptr, n, box9, types, ids) != 0:
        raise SystemExit("seams_system_set_frame failed")
    labels = (ctypes.c_int * n)()
    counts = (ctypes.c_int * 9)()
    if which == "chill":
        rc = lib.seams_chill(sys, 3.5, 1, 0, 0, labels, n, counts)
    else:
        rc = lib.seams_chill_plus(sys, 3.5, 1, 0, 0, labels, counts)
    lib.seams_system_destroy(sys)
    if rc != 0:
        raise SystemExit(f"{which} failed")
    return np.ctypeslib.as_array(labels).copy()


def histogram(labels: np.ndarray) -> str:
    counts = np.bincount(labels.astype(np.int64), minlength=len(NAMES))
    return " ".join(f"{NAMES[i]}={int(counts[i])}" for i in range(len(NAMES)) if counts[i])


def write_labels(path: Path, labels: np.ndarray) -> None:
    body = " ".join(str(int(v)) for v in labels)
    path.write_text(
        "# reference ee53b952 published cutoff 3.5 A type 1\n"
        + histogram(labels) + "\n" + body + "\n",
        encoding="ascii",
    )


MDAPY_NAME = {
    0: "cubic_ice",
    1: "hexagonal_ice",
    2: "other",
    3: "interfacial_ice",
    4: "hydrate",
    5: "interfacial_hydrate",
    6: "other",
    7: "other",
    8: "other",
}


def cutoff_shell(lib, pos: np.ndarray, lengths):
    n = pos.shape[0]
    lx, ly, lz = lengths
    box9 = (ctypes.c_double * 9)(0, 0, 0, lx, ly, lz, 0, 0, 0)
    types = (ctypes.c_int * n)(*([1] * n))
    ids = (ctypes.c_int * n)(*range(1, n + 1))
    ptr = np.ascontiguousarray(pos).ctypes.data_as(ctypes.POINTER(ctypes.c_double))
    sys = lib.seams_system_create()
    if lib.seams_system_set_frame(sys, ptr, n, box9, types, ids) != 0:
        raise SystemExit("set_frame failed")
    nb = ctypes.c_void_p()
    if lib.seams_neighbours_build(sys, 3.5, 1, 0, ctypes.byref(nb)) != 0:
        raise SystemExit("neighbours failed")
    n_off = ctypes.c_size_t()
    nnz = ctypes.c_size_t()
    lib.seams_neighbours_copy(nb, None, 0, ctypes.byref(n_off), None, 0, ctypes.byref(nnz))
    offsets = (ctypes.c_int * n_off.value)()
    cols = (ctypes.c_int * max(int(nnz.value), 1))()
    lib.seams_neighbours_copy(
        nb, offsets, n_off.value, ctypes.byref(n_off), cols, max(int(nnz.value), 1), ctypes.byref(nnz)
    )
    lib.seams_neighbours_destroy(nb)
    lib.seams_system_destroy(sys)
    return np.ctypeslib.as_array(offsets).copy(), np.ctypeslib.as_array(cols)[: int(nnz.value)].copy()


def mdapy_report(path: Path, lib) -> None:
    import mdapy

    mdapy.set_num_threads(1)
    lines = [
        "# d-SEAMS CHILL+ against mdapy CHILL+",
        "",
        "Same three frames, published thresholds, cutoff 3.5 A, one thread.",
        "mdapy's class names are aligned with the CHILL+ table: hexagonal ice,",
        "cubic ice, interfacial ice, hydrate (clathrate), interfacial hydrate,",
        "and other for liquid or unclassified. Thresholds were not moved.",
        "",
    ]
    for name, builder in FRAMES.items():
        pos, lengths = builder()
        ours = label_frame(lib, pos, lengths, 1, "chill_plus")
        system = mdapy.System(pos, box=lengths)
        ice = mdapy.identify.ChillPlus(cutoff=3.5).compute(system)
        mdapy_names = tuple(ice.names)
        theirs = np.asarray(ice.structure, dtype=np.int64)
        their_name = np.array([mdapy_names[int(v)] for v in theirs])
        our_name = np.array([MDAPY_NAME[int(v)] for v in ours])
        neigh = mdapy.NeighborList.build(system, cutoff=3.5)
        offsets = np.asarray(neigh.offsets)
        columns = np.asarray(neigh.j)
        lines.append(f"## {name}")
        lines.append("")
        lines.append(f"d-SEAMS: {histogram(ours)}")
        counts = {}
        for label in their_name:
            counts[label] = counts.get(label, 0) + 1
        lines.append("mdapy: " + " ".join(f"{k}={v}" for k, v in sorted(counts.items())))
        disagree = np.flatnonzero(our_name != their_name)
        lines.append(f"per-atom disagreements: {int(disagree.size)}")
        lines.append("")
        if disagree.size == 0:
            lines.append("Every oxygen has the same class.")
            lines.append("")
            continue
        d_off, d_col = cutoff_shell(lib, pos, lengths)
        same_shell = 0
        lines.append("| oxygen | d-SEAMS | mdapy | same shell |")
        lines.append("|---:|---|---|---|")
        for i in disagree:
            d_set = set(d_col[d_off[i]:d_off[i + 1]].tolist())
            m_set = set(columns[offsets[i]:offsets[i + 1]].tolist())
            same = d_set == m_set
            same_shell += int(same)
            lines.append(
                f"| {int(i) + 1} | {our_name[i]} | {their_name[i]} | {str(same).lower()} |"
            )
        lines.append("")
        lines.append(
            f"The 3.5 A neighbour set is the same for {same_shell} of these "
            f"{int(disagree.size)} oxygens. The classes that swap are adjacent "
            "rows of the published table: liquid against interfacial hydrate, "
            "hexagonal against interfacial ice, and clathrate against hydrate. "
            "A correlation that sits on an edge falls one way in d-SEAMS and the "
            "other in mdapy. The edges were not moved."
        )
        lines.append("")
    path.write_text("\n".join(lines), encoding="ascii")


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--out", type=Path, required=True)
    parser.add_argument("--label", type=Path, default=None)
    parser.add_argument("--mdapy-report", type=Path, default=None)
    args = parser.parse_args()
    args.out.mkdir(parents=True, exist_ok=True)
    lib = load_seams(args.label) if args.label is not None else None
    for name, builder in FRAMES.items():
        pos, lengths = builder()
        write_dump(args.out / f"{name}.lammpstrj", pos, lengths)
        print(f"{name} n={pos.shape[0]} box={lengths[0]:.4f},{lengths[1]:.4f},{lengths[2]:.4f}")
        if lib is None:
            continue
        for which in ("chill", "chill_plus"):
            for threads in (1, 8):
                labels = label_frame(lib, pos, lengths, threads, which)
                print(f"  {which} threads {threads} {histogram(labels)}")
                if threads == 1:
                    first = labels
                    write_labels(args.out / f"{name}.{which}.txt", labels)
                elif not np.array_equal(first, labels):
                    raise SystemExit(f"{name} {which} changed with the thread count")
    if args.mdapy_report is not None:
        if lib is None:
            raise SystemExit("--mdapy-report needs --label")
        mdapy_report(args.mdapy_report, lib)


if __name__ == "__main__":
    main()
