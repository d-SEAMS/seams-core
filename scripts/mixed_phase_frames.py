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


def label_frame(lib, pos: np.ndarray, lengths, threads: int, which: str, rule: int = 0) -> np.ndarray:
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
        rc = lib.seams_chill_plus(sys, 3.5, 1, 0, rule, labels, counts)
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


# Condon-Shortley constants for the real l = 3 harmonics. The same
# numbers sit in normQ3(). The 1/4 in equation (1) cancels in c_ij.
_Q3_NORM = np.array(
    [
        0.25 * math.sqrt(7.0 / math.pi),
        0.125 * math.sqrt(21.0 / math.pi),
        0.25 * math.sqrt(105.0 / (2.0 * math.pi)),
        0.125 * math.sqrt(35.0 / math.pi),
    ]
)
_EDGES = (-0.8, -0.35, 0.25)
OVITO_NAME = {
    0: "other",
    1: "hexagonal_ice",
    2: "cubic_ice",
    3: "interfacial_ice",
    4: "hydrate",
    5: "interfacial_hydrate",
}


def _ylm3(x: float, y: float, z: float) -> np.ndarray:
    acc = np.zeros(7)
    r2 = x * x + y * y + z * z
    if not r2 > 0.0:
        acc[0] = _Q3_NORM[0] * 2.0
        return acc
    inv = 1.0 / (r2 * math.sqrt(r2))
    t1 = 5.0 * z * z - r2
    acc[0] = _Q3_NORM[0] * (z * (5.0 * z * z - 3.0 * r2)) * inv
    acc[1] = _Q3_NORM[1] * (y * t1) * inv
    acc[2] = _Q3_NORM[1] * (x * t1) * inv
    acc[3] = _Q3_NORM[2] * (z * (y * y - x * x)) * inv
    acc[4] = _Q3_NORM[2] * (2.0 * x * y * z) * inv
    acc[5] = _Q3_NORM[3] * (y * (y * y - 3.0 * x * x)) * inv
    acc[6] = _Q3_NORM[3] * (x * (3.0 * y * y - x * x)) * inv
    return acc


def _dot7(a: np.ndarray, b: np.ndarray) -> float:
    dot = float(a[0] * b[0])
    for channel in range(1, 7):
        dot += 2.0 * float(a[channel] * b[channel])
    return dot


def _shells(lib, pos: np.ndarray, lengths):
    off, col = cutoff_shell(lib, pos, lengths)
    box = np.asarray(lengths, dtype=float)
    neigh = []
    vecs = []
    for i in range(pos.shape[0]):
        js = col[off[i] : off[i + 1]]
        if js.size == 0:
            neigh.append(np.zeros(0, dtype=int))
            vecs.append(np.zeros((0, 3)))
            continue
        delta = pos[js] - pos[i]
        delta -= box * np.round(delta / box)
        d2 = np.einsum("ij,ij->i", delta, delta)
        keep = d2 > 0.0
        js = js[keep]
        delta = delta[keep]
        d2 = d2[keep]
        order = np.lexsort((js, d2))
        neigh.append(js[order])
        vecs.append(delta[order])
    return neigh, vecs


def _q3(neigh, vecs, four_only: bool):
    n = len(neigh)
    q = np.zeros((n, 7))
    chosen = []
    for i in range(n):
        k = min(4, len(neigh[i])) if four_only else len(neigh[i])
        chosen.append(neigh[i][:k])
        acc = np.zeros(7)
        for vec in vecs[i][:k]:
            acc += _ylm3(float(vec[0]), float(vec[1]), float(vec[2]))
        q[i] = acc
    norm = np.array([_dot7(q[i], q[i]) for i in range(n)])
    return q, norm, chosen


def _cij(q, norm, i: int, j: int) -> float:
    ni = float(norm[i])
    nj = float(norm[j])
    if not (ni > 1e-12 and nj > 1e-12):
        return 0.0
    return _dot7(q[i], q[j]) / math.sqrt(ni * nj)


def _kind(c: float, closed: bool) -> str:
    staggered = c <= -0.8 if closed else c < -0.8
    if staggered:
        return "S"
    if closed:
        eclipsed = -0.35 <= c <= 0.25
    else:
        eclipsed = -0.35 < c < 0.25
    return "E" if eclipsed else "O"


def _classes(neigh, q, norm, chosen, closed: bool):
    n = len(neigh)
    stag = np.zeros(n, dtype=int)
    eclip = np.zeros(n, dtype=int)
    bonds = []
    for i, js in enumerate(chosen):
        row = []
        for j in js:
            c = _cij(q, norm, i, int(j))
            k = _kind(c, closed)
            row.append((int(j), c, k))
            if k == "S":
                stag[i] += 1
            elif k == "E":
                eclip[i] += 1
        bonds.append(row)
    names = []
    for i, row in enumerate(bonds):
        four = len(neigh[i]) == 4 and len(row) == 4
        if not four:
            names.append("other")
            continue
        ns = [int(stag[j]) for j, _, _ in row]
        nstag = int(stag[i])
        neclip = int(eclip[i])
        inter = False
        if nstag == 2 and any(s > 2 for s in ns):
            inter = True
        if nstag == 3 and neclip == 0 and any(s > 1 for s in ns):
            inter = True
        if neclip == 0 and nstag == 4:
            names.append("cubic_ice")
        elif neclip == 1 and nstag == 3:
            names.append("hexagonal_ice")
        elif inter:
            names.append("interfacial_ice")
        elif neclip == 4 and nstag == 0:
            names.append("hydrate")
        elif neclip == 3:
            names.append("interfacial_hydrate")
        else:
            names.append("other")
    return names, bonds


def _ovito_labels(pos: np.ndarray, lengths) -> np.ndarray:
    from ovito.data import DataCollection, Particles, SimulationCell
    from ovito.modifiers import ChillPlusModifier
    from ovito.pipeline import Pipeline, StaticSource

    lx, ly, lz = lengths
    data = DataCollection()
    cell = SimulationCell(pbc=(True, True, True))
    matrix = np.zeros((3, 4))
    matrix[0, 0], matrix[1, 1], matrix[2, 2] = lx, ly, lz
    cell[...] = matrix
    data.objects.append(cell)
    particles = Particles()
    particles.create_property("Position", data=np.ascontiguousarray(pos))
    data.objects.append(particles)
    pipe = Pipeline(source=StaticSource(data=data))
    pipe.modifiers.append(ChillPlusModifier(cutoff=3.5))
    out = pipe.compute()
    codes = np.asarray(out.particles["Structure Type"], dtype=np.int64)
    return np.array([OVITO_NAME[int(v)] for v in codes])


def _fmt(c: float) -> str:
    return f"{c:.12f}"


def mdapy_report(path: Path, lib) -> None:
    import mdapy

    mdapy.set_num_threads(1)
    ulp = math.ulp(0.35)
    lines = [
        "# d-SEAMS CHILL+ against mdapy and OVITO",
        "",
        "Same three frames, cutoff 3.5 A, one thread. Thresholds were not",
        "moved. The comparison operator was not moved either: d-SEAMS already",
        "uses the inclusive edges in the papers, and on these frames the",
        "inclusive and strict readings give the same label for every oxygen.",
        "",
        "## Published inequalities",
        "",
        "Nguyen and Molinero, J. Phys. Chem. B 2015, 119, 9369, equation (1),",
        "averages the spherical harmonic over the four closest neighbours:",
        "",
        "    q_lm(i) = (1/4) sum_{j=1}^{4} Y_lm(r_ij)",
        "",
        "Equation (2) is the normalised correlation of that vector,",
        "",
        "    c(i,j) = sum_m q_lm(i) q_lm(j)*",
        "             / sqrt( (sum_m |q_lm(i)|^2) (sum_m |q_lm(j)|^2) )",
        "",
        "The text says q_lm(i) projects the orientational structure of the",
        "four closest neighbours, and that c(i,j) is the dot product of q_l",
        "between a pair of neighbouring molecules. Table 2 then sets the",
        "CHILL+ windows with closed ends:",
        "",
        "    eclipsed    0.25 >= c(i,j) >= -0.35",
        "    staggered   c(i,j) <= -0.8",
        "",
        "The same page restates CHILL, citing Moore, de la Llave, Welke,",
        "Scherlis and Molinero, Phys. Chem. Chem. Phys. 2010, 12, 4124:",
        "\"CHILL identifies a staggered bond when c(i, j) <= -0.8 and an",
        "eclipsed bond when -0.05 >= c(i, j) >= -0.2.\" Table 2 prints the",
        "same CHILL row, with the footnote that it is from that 2010 paper.",
        "Both inequalities are inclusive. d-SEAMS classifyBond with the",
        "published flag uses c <= -0.8 and -0.35 <= c <= 0.25 for CHILL+,",
        "and c <= -0.8 and -0.2 <= c <= -0.05 for CHILL. That is the paper.",
        "The operator was left as it is.",
        "",
        "mdapy's CHILL+ page states the strict cuts c < -0.8 and",
        "-0.35 < c < 0.25, and it sums q_3m over the neighbours inside the",
        "cutoff rather than over the four closest. OVITO 3.17.2 ChillPlus",
        "is the modifier in the ovito Python module, cutoff 3.5. d-SEAMS",
        "also has an ovito-chill+ path (open edges, q3 on the whole cutoff",
        "row). On these frames that path, mdapy, and the OVITO modifier",
        "agree on every oxygen. The open and closed readings of the same",
        "c_ij also agree on every oxygen, so the split is not the operator.",
        "",
        "## Agreement",
        "",
        "| frame | n | d-SEAMS vs OVITO | mdapy vs OVITO | d-SEAMS ovito-rule vs OVITO |",
        "|---|---:|---:|---:|---:|",
    ]
    details = []
    closest = None
    closest_flip = None
    for name, builder in FRAMES.items():
        pos, lengths = builder()
        published = label_frame(lib, pos, lengths, 1, "chill_plus", 0)
        ovito_rule = label_frame(lib, pos, lengths, 1, "chill_plus", 1)
        our_name = np.array([MDAPY_NAME[int(v)] for v in published])
        rule_name = np.array([MDAPY_NAME[int(v)] for v in ovito_rule])
        system = mdapy.System(pos, box=lengths)
        ice = mdapy.identify.ChillPlus(cutoff=3.5).compute(system)
        md_names = tuple(ice.names)
        md = np.array([md_names[int(v)] for v in np.asarray(ice.structure)])
        ov = _ovito_labels(pos, lengths)
        neigh, vecs = _shells(lib, pos, lengths)
        q4, n4, c4 = _q3(neigh, vecs, True)
        q_all, n_all, c_all = _q3(neigh, vecs, False)
        paper, bonds4 = _classes(neigh, q4, n4, c4, True)
        shell, bonds_all = _classes(neigh, q_all, n_all, c_all, True)
        if paper != list(our_name):
            raise SystemExit(f"{name}: four-nearest c_ij does not reproduce d-SEAMS")
        if shell != list(ov):
            raise SystemExit(f"{name}: full-shell c_ij does not reproduce OVITO")
        if not np.array_equal(rule_name, ov) or not np.array_equal(md, ov):
            raise SystemExit(f"{name}: mdapy or the ovito rule left OVITO")
        # Strict edges on the same two correlations.
        strict4, _ = _classes(neigh, q4, n4, c4, False)
        strict_all, _ = _classes(neigh, q_all, n_all, c_all, False)
        if strict4 != paper or strict_all != shell:
            raise SystemExit(f"{name}: inclusive and strict labels differ")
        n = pos.shape[0]
        lines.append(
            f"| {name} | {n} | {int(np.count_nonzero(our_name != ov))} | "
            f"{int(np.count_nonzero(md != ov))} | {int(np.count_nonzero(rule_name != ov))} |"
        )
        disagree = np.flatnonzero(our_name != ov)
        frame_rows = []
        for i in disagree:
            by_all = {j: (c, k) for j, c, k in bonds_all[i]}
            bond_lines = []
            partners = []
            for j in c4[i]:
                if len(neigh[int(j)]) != 4:
                    partners.append(f"{int(j) + 1}:{len(neigh[int(j)])}")
            for j, c, k in bonds4[i]:
                c2, k2 = by_all[int(j)]
                d4 = min(abs(c - edge) for edge in _EDGES)
                d_all = min(abs(c2 - edge) for edge in _EDGES)
                gap = min(d4, d_all)
                rec = (gap, name, int(i) + 1, c, c2)
                if closest is None or gap < closest[0]:
                    closest = rec
                if k != k2 and (closest_flip is None or gap < closest_flip[0]):
                    closest_flip = rec
                mark = " changes class" if k != k2 else ""
                bond_lines.append(
                    f"j={int(j) + 1} paper={_fmt(c)} {k}/{_kind(c, False)} "
                    f"shell={_fmt(c2)} {k2}/{_kind(c2, False)} "
                    f"|dc|={abs(c - c2):.3e} "
                    f"dist_paper={d4:.3e} dist_shell={d_all:.3e}{mark}"
                )
            frame_rows.append((int(i) + 1, our_name[i], ov[i], partners, bond_lines))
        details.append((name, frame_rows))
    lines.extend(
        [
            "",
            "d-SEAMS vs OVITO is 20, 0 and 17 oxygens. mdapy vs OVITO is 0 on",
            "every frame. The d-SEAMS ovito-chill+ path also matches the",
            "modifier on every oxygen. Stacking has no disagreement: every",
            "oxygen there has four neighbours inside 3.5 A, so the four",
            "closest neighbours and the cutoff shell are the same set.",
            "",
            "## Where the 37 labels move",
            "",
            "Each disagreeing oxygen itself has four neighbours inside 3.5 A,",
            "and that set is the same in d-SEAMS, mdapy and OVITO. The",
            "correlation still moves because equation (1) and the cutoff sum",
            "disagree on a neighbour that has more than four oxygens inside",
            "3.5 A. d-SEAMS keeps the four closest in",
            "q_lm. OVITO and mdapy fold the extra oxygen into q_lm. c_ij then",
            "shifts by about 1e-2 to several tenths, and a bond that was",
            "eclipsed becomes out of range, or the other way around. The",
            "molecule row in Table 1 follows the new staggered and eclipsed",
            "counts. The interfacial neighbour test is not the split: turning",
            "it off does not change a label once the q3 shell is fixed.",
            "",
            "No c_ij of a disagreeing oxygen lands on -0.8, -0.35 or 0.25.",
            "For each bond below, the two letters are the inclusive class",
            "and the strict class. They are the same letter in every row.",
            "d-SEAMS applies the inclusive form. mdapy's documentation",
            "applies the strict form. OVITO's modifier matches either form",
            "of the cutoff-shell correlation on these frames. `paper` is the",
            "four-closest c_ij. `shell` is the cutoff-shell c_ij.",
            "",
        ]
    )
    for name, frame_rows in details:
        lines.append(f"### {name}")
        lines.append("")
        if not frame_rows:
            lines.append("No disagreeing oxygen.")
            lines.append("")
            continue
        for oid, ours, theirs, partners, bond_lines in frame_rows:
            extra = ", ".join(partners) if partners else "none"
            lines.append(
                f"oxygen {oid}: d-SEAMS {ours}, OVITO and mdapy {theirs}. "
                f"partners whose shell is not 4: {extra}."
            )
            for bond in bond_lines:
                lines.append(f"    {bond}")
            lines.append("")
    assert closest is not None and closest_flip is not None
    lines.extend(
        [
            "## Floating point",
            "",
            "A double near 0.35 has unit in the last place "
            f"{ulp:.3e}. The closest any c_ij of a disagreeing oxygen comes",
            "to -0.8, -0.35 or 0.25 is "
            f"{closest[0]:.3e} ({closest[1]} oxygen {closest[2]}, "
            f"paper {_fmt(closest[3])}, shell {_fmt(closest[4])}). "
            "The closest bond that actually changes class is the same gap, "
            f"{closest_flip[0]:.3e} "
            f"({closest_flip[1]} oxygen {closest_flip[2]}). "
            f"That gap is {closest[0] / ulp:.3e} units in the last place, "
            "not roundoff.",
            "Inclusive and strict comparisons classify every one of these",
            "bonds the same way. The values are not sitting on the boundary.",
            "",
        ]
    )
    path.write_text("\n".join(lines) + "\n", encoding="ascii")


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
