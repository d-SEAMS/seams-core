"""Per-atom and averaged q4 and q6 against freud and pyscal3.

OVITO 3.15.5 has no SteinhardtParametersModifier, so this comparison
does not call OVITO. The three codes use the same rule: a neighbour
is a pair whose minimum-image distance is strictly below 3.5 Å.
Ideal ice Ih is input/traj/ice1h.lammpstrj. Ideal ice Ic is the
diamond cell at a 2.75 Å bond. The liquid is that cell with Gaussian
noise of 0.9 Å, seed 91001.

pyscal3 is float64. Local q_l and the Lechner-Dellago average, and
the mean of each over atoms that have a shell, must agree to 1e-10.

freud 3.6 stores q_l as float32. On the two ideal crystals the gap is
at most eight units in the last place of that float32 value. On the liquid
the spherical-harmonic sum is accumulated in float32, so the residual
is larger. The test prints that residual and requires it to stay
under 2e-5. An empty shell is 0 here and NaN in freud and pyscal3.
"""

from __future__ import annotations

import ctypes
import math
import os
import subprocess
import sys
from pathlib import Path

try:
    import numpy as np
except ImportError:
    sys.stderr.write("numpy is not installed\n")
    raise SystemExit(77)

try:
    import freud
    import pyscal3
    from ase import Atoms
except ImportError:
    sys.stderr.write("freud, pyscal3, or ase is not installed\n")
    raise SystemExit(77)

LIB = os.environ.get("SEAMS_LANG_LIB")
if not LIB:
    sys.stderr.write("SEAMS_LANG_LIB is not set\n")
    raise SystemExit(77)

CUTOFF = 3.5
PYSCAL_TOL = 1e-10
# Units in the last place of freud's float32 q_l on a crystal.
FREUD_ULP = 8
# float32 accumulation of a disordered shell. Not a 1e-10 comparison.
FREUD_LIQUID = 2e-5
ROOT = Path(__file__).resolve().parents[1]
IH_DUMP = ROOT / "input" / "traj" / "ice1h.lammpstrj"


def ice_ic(cells, bond=2.75):
    lattice = 4.0 * bond / math.sqrt(3.0)
    side = cells * lattice
    fcc = np.array([[0.0, 0.0, 0.0], [0.5, 0.5, 0.0], [0.5, 0.0, 0.5], [0.0, 0.5, 0.5]])
    basis = np.vstack([fcc, fcc + 0.25])
    ijk = np.stack(np.meshgrid(*[np.arange(cells)] * 3, indexing="ij"), -1).reshape(-1, 3)
    pts = (ijk[:, None, :] + basis).reshape(-1, 3) * lattice
    return np.ascontiguousarray(pts % side), (float(side), float(side), float(side))


def liquid(cells):
    pts, edges = ice_ic(cells)
    rng = np.random.default_rng(91001)
    noisy = np.mod(pts + rng.normal(0.0, 0.9, pts.shape), edges[0])
    return np.ascontiguousarray(noisy), edges


def ice_ih():
    lines = IH_DUMP.read_text().splitlines()
    index = 0
    while "BOX BOUNDS" not in lines[index]:
        index += 1
    index += 1
    edges = []
    for _ in range(3):
        lo, hi = (float(part) for part in lines[index].split()[:2])
        edges.append(hi - lo)
        index += 1
    while not lines[index].startswith("ITEM: ATOMS"):
        index += 1
    columns = {name: col for col, name in enumerate(lines[index].split()[2:])}
    index += 1
    pts = []
    while index < len(lines) and lines[index] and not lines[index].startswith("ITEM"):
        parts = lines[index].split()
        pts.append([float(parts[columns[axis]]) for axis in "xyz"])
        index += 1
    return np.ascontiguousarray(pts, dtype=np.float64), tuple(edges)


def seams_ql(lib, pts, edges, order):
    n = len(pts)
    xyz = np.ascontiguousarray(pts, dtype=np.float64).reshape(-1)
    ox, oy, oz = 0.0, 0.0, 0.0
    lx, ly, lz = edges
    box = (ctypes.c_double * 9)(ox, oy, oz, lx, ly, lz, 0.0, 0.0, 0.0)
    types = (ctypes.c_int * n)(*([1] * n))
    ids = (ctypes.c_int * n)(*range(1, n + 1))
    frame = lib.seams_frame_from_arrays(
        xyz.ctypes.data_as(ctypes.POINTER(ctypes.c_double)), n, box, types, ids)
    if not frame:
        raise RuntimeError("frame")
    ql = (ctypes.c_double * n)()
    bar = (ctypes.c_double * n)()
    try:
        if lib.seams_frame_steinhardt(frame, CUTOFF, 1, order, ql, bar) != 0:
            raise RuntimeError("steinhardt")
        return np.ctypeslib.as_array(ql).copy(), np.ctypeslib.as_array(bar).copy()
    finally:
        lib.seams_frame_free(frame)


def freud_ql(pts, edges, order, averaged):
    freud.parallel.set_num_threads(int(os.environ.get("OMP_NUM_THREADS", "1")))
    box = freud.box.Box(*edges)
    op = freud.order.Steinhardt(l=order, average=averaged)
    op.compute((box, pts), neighbors={"r_max": CUTOFF})
    return np.asarray(op.ql, dtype=np.float64).copy()


def pyscal_ql(pts, edges, order, averaged):
    if hasattr(pyscal3, "set_num_threads"):
        pyscal3.set_num_threads(int(os.environ.get("OMP_NUM_THREADS", "1")))
    atoms = Atoms(
        symbols=["H"] * len(pts),
        positions=np.ascontiguousarray(pts, dtype=np.float64),
        cell=list(edges),
        pbc=True,
    )
    pyscal3.find_neighbors(atoms, method="cutoff", cutoff=CUTOFF)
    values = pyscal3.steinhardt_parameter(atoms, order, averaged=averaged)
    return np.asarray(values[0], dtype=np.float64).copy()


def live_gap(ours, theirs):
    mask = np.isfinite(ours) & np.isfinite(theirs)
    if not np.any(mask):
        raise SystemExit("no shell is finite in both codes")
    return float(np.max(np.abs(ours[mask] - theirs[mask]))), mask


def freud_gap(name, kind, ours, theirs):
    """Return a note. Crystal gaps are float32 ulps. Liquid gaps are the float32 sum."""
    gap, mask = live_gap(ours, theirs)
    if name == "liquid":
        if gap > FREUD_LIQUID:
            raise SystemExit(
                f"liquid {kind} freud residual {gap:.3e} is past the float32 accumulation bound"
            )
        return f"freud {kind} {gap:.3e} (float32 sum)"
    rounded = ours[mask].astype(np.float32)
    ulp = np.spacing(rounded).astype(np.float64)
    steps = float(np.max(np.abs(theirs[mask] - rounded.astype(np.float64)) / ulp))
    if steps > FREUD_ULP:
        raise SystemExit(f"{name} {kind} freud {steps:.2f} ulp, gap {gap:.3e}")
    return f"freud {kind} {steps:.2f} ulp"


def check_order(lib, name, pts, edges, order):
    local, bar = seams_ql(lib, pts, edges, order)
    py_local = pyscal_ql(pts, edges, order, averaged=False)
    py_bar = pyscal_ql(pts, edges, order, averaged=True)
    fr_local = freud_ql(pts, edges, order, averaged=False)
    fr_bar = freud_ql(pts, edges, order, averaged=True)

    empty = ~np.isfinite(fr_local)
    if not np.array_equal(empty, ~np.isfinite(py_local)):
        raise SystemExit(f"{name} l={order} empty shells differ between freud and pyscal3")
    if np.any(empty):
        if float(np.max(np.abs(local[empty]))) != 0.0 or float(np.max(np.abs(bar[empty]))) != 0.0:
            raise SystemExit(f"{name} l={order} empty shell is not 0")
        if np.any(np.isfinite(py_bar[empty])) or np.any(np.isfinite(fr_bar[empty])):
            raise SystemExit(f"{name} l={order} averaged empty shell is finite")
        print(
            f"{name} l={order} empty {int(empty.sum())}: "
            "d-SEAMS writes 0, freud and pyscal3 write NaN"
        )

    d_local, live = live_gap(local, py_local)
    d_bar, live_bar = live_gap(bar, py_bar)
    if d_local > PYSCAL_TOL or d_bar > PYSCAL_TOL:
        raise SystemExit(f"{name} l={order} pyscal3 local {d_local:.3e} bar {d_bar:.3e}")
    mean_local = abs(float(np.mean(local[live])) - float(np.mean(py_local[live])))
    mean_bar = abs(float(np.mean(bar[live_bar])) - float(np.mean(py_bar[live_bar])))
    if mean_local > PYSCAL_TOL or mean_bar > PYSCAL_TOL:
        raise SystemExit(f"{name} l={order} mean local {mean_local:.3e} bar {mean_bar:.3e}")

    note_local = freud_gap(name, "q", local, fr_local)
    note_bar = freud_gap(name, "qbar", bar, fr_bar)
    print(
        f"{name} n={len(pts)} l={order} threads={os.environ.get('OMP_NUM_THREADS')} "
        f"pyscal3 q {d_local:.3e} qbar {d_bar:.3e} "
        f"mean {mean_local:.3e} {mean_bar:.3e} {note_local} {note_bar}"
    )


def check(lib):
    frames = (
        ("Ih",) + ice_ih(),
        ("Ic",) + ice_ic(4),
        ("liquid",) + liquid(4),
    )
    for name, pts, edges in frames:
        for order in (4, 6):
            check_order(lib, name, pts, edges, order)


def load():
    lib = ctypes.CDLL(LIB)
    lib.seams_frame_from_arrays.restype = ctypes.c_void_p
    lib.seams_frame_from_arrays.argtypes = [
        ctypes.POINTER(ctypes.c_double), ctypes.c_int,
        ctypes.POINTER(ctypes.c_double),
        ctypes.POINTER(ctypes.c_int), ctypes.POINTER(ctypes.c_int)]
    lib.seams_frame_free.argtypes = [ctypes.c_void_p]
    lib.seams_frame_steinhardt.restype = ctypes.c_int
    lib.seams_frame_steinhardt.argtypes = [
        ctypes.c_void_p, ctypes.c_double, ctypes.c_int, ctypes.c_int,
        ctypes.POINTER(ctypes.c_double), ctypes.POINTER(ctypes.c_double)]
    return lib


def main():
    if os.environ.get("QL_EQUAL_CHILD") != "1":
        for threads in ("1", "8"):
            env = os.environ.copy()
            env["OMP_NUM_THREADS"] = threads
            env["QL_EQUAL_CHILD"] = "1"
            done = subprocess.run([sys.executable, __file__], env=env)
            if done.returncode != 0:
                raise SystemExit(done.returncode)
        return
    check(load())


if __name__ == "__main__":
    main()
