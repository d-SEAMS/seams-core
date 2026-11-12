"""q-bar and normalised w on ice Ih, ice Ic and a liquid, against pyscal3."""

import ctypes
import math
import os
import sys

try:
    import numpy as np
except ImportError:
    sys.stderr.write("numpy is not installed\n")
    raise SystemExit(77)

try:
    import pyscal3
    from ase import Atoms
except ImportError:
    sys.stderr.write("pyscal3 is not installed\n")
    raise SystemExit(77)

LIB = os.environ.get("SEAMS_LANG_LIB")
if not LIB:
    sys.stderr.write("SEAMS_LANG_LIB is not set\n")
    raise SystemExit(77)

CUTOFF = 3.5
TOL = 1e-6


def load():
    lib = ctypes.CDLL(LIB)
    lib.seams_frame_from_arrays.restype = ctypes.c_void_p
    lib.seams_frame_from_arrays.argtypes = [
        ctypes.POINTER(ctypes.c_double), ctypes.c_int,
        ctypes.POINTER(ctypes.c_double), ctypes.POINTER(ctypes.c_int),
        ctypes.POINTER(ctypes.c_int)]
    lib.seams_frame_free.argtypes = [ctypes.c_void_p]
    lib.seams_frame_steinhardt.restype = ctypes.c_int
    lib.seams_frame_steinhardt.argtypes = [
        ctypes.c_void_p, ctypes.c_double, ctypes.c_int, ctypes.c_int,
        ctypes.POINTER(ctypes.c_double), ctypes.POINTER(ctypes.c_double)]
    lib.seams_frame_wl.restype = ctypes.c_int
    lib.seams_frame_wl.argtypes = [
        ctypes.c_void_p, ctypes.c_double, ctypes.c_int, ctypes.c_int,
        ctypes.POINTER(ctypes.c_double)]
    lib.seams_frame_q4_q6_chill.restype = ctypes.c_int
    lib.seams_frame_q4_q6_chill.argtypes = [
        ctypes.c_void_p, ctypes.c_double, ctypes.c_int,
        ctypes.POINTER(ctypes.c_double), ctypes.POINTER(ctypes.c_double),
        ctypes.POINTER(ctypes.c_double), ctypes.POINTER(ctypes.c_double),
        ctypes.POINTER(ctypes.c_int), ctypes.POINTER(ctypes.c_int)]
    return lib


def ice_ic(cells=3):
    bond = 2.75
    lattice = 4.0 * bond / math.sqrt(3.0)
    side = cells * lattice
    fcc = np.array([[0.0, 0, 0], [0.5, 0.5, 0], [0.5, 0, 0.5], [0, 0.5, 0.5]])
    basis = np.vstack([fcc, fcc + 0.25])
    ijk = np.stack(np.meshgrid(*[np.arange(cells)] * 3, indexing="ij"), -1).reshape(-1, 3)
    pts = ((ijk[:, None, :] + basis).reshape(-1, 3) * lattice) % side
    return np.ascontiguousarray(pts), np.array([side, side, side])


def ice_ih(nx=4, ny=2, nz=2):
    # Ideal P6_3/mmc oxygen framework. a is scaled so the four
    # nearest neighbours sit at 2.75 A. The orthorhombic cell holds
    # two hexagonal cells.
    a = 4.5 * 2.75 / 2.755675956773843
    c = a * math.sqrt(8.0 / 3.0)
    frac = np.array([
        [1.0 / 3.0, 2.0 / 3.0, 1.0 / 16.0],
        [2.0 / 3.0, 1.0 / 3.0, 9.0 / 16.0],
        [2.0 / 3.0, 1.0 / 3.0, 15.0 / 16.0],
        [1.0 / 3.0, 2.0 / 3.0, 7.0 / 16.0],
    ])
    a1 = np.array([a, 0.0, 0.0])
    a2 = np.array([-a / 2.0, a * math.sqrt(3.0) / 2.0, 0.0])
    a3 = np.array([0.0, 0.0, c])
    cell_box = np.array([a, a * math.sqrt(3.0), c])
    sites = []
    for f in frac:
        p = f[0] * a1 + f[1] * a2 + f[2] * a3
        sites.extend((p, p + a2))
    sites = np.array(sites)
    sites = sites - cell_box * np.floor(sites / cell_box)
    key = np.round(sites, 6)
    _, index = np.unique(key, axis=0, return_index=True)
    sites = sites[np.sort(index)]
    blocks = []
    for ix in range(nx):
        for iy in range(ny):
            for iz in range(nz):
                blocks.append(sites + np.array([ix, iy, iz]) * cell_box)
    return np.ascontiguousarray(np.vstack(blocks)), cell_box * np.array([nx, ny, nz], float)


def liquid():
    rng = np.random.default_rng(91001)
    n = 200
    side = (n / 0.033) ** (1.0 / 3.0)
    return np.ascontiguousarray(rng.random((n, 3)) * side), np.array([side, side, side])


def frame(lib, pts, box):
    n = len(pts)
    flat = np.ascontiguousarray(pts.reshape(-1))
    box9 = np.array([0.0, 0.0, 0.0, box[0], box[1], box[2], 0.0, 0.0, 0.0])
    types = np.ones(n, dtype=np.int32)
    ids = np.arange(1, n + 1, dtype=np.int32)
    return lib.seams_frame_from_arrays(
        flat.ctypes.data_as(ctypes.POINTER(ctypes.c_double)), n,
        box9.ctypes.data_as(ctypes.POINTER(ctypes.c_double)),
        types.ctypes.data_as(ctypes.POINTER(ctypes.c_int)),
        ids.ctypes.data_as(ctypes.POINTER(ctypes.c_int))), n


def seams_params(lib, pts, box):
    fr, n = frame(lib, pts, box)
    out = {}
    for order in (4, 6):
        ql = (ctypes.c_double * n)()
        qb = (ctypes.c_double * n)()
        wl = (ctypes.c_double * n)()
        if lib.seams_frame_steinhardt(fr, CUTOFF, 1, order, ql, qb) != 0:
            raise RuntimeError("steinhardt failed")
        if lib.seams_frame_wl(fr, CUTOFF, 1, order, wl) != 0:
            raise RuntimeError("wl failed")
        out[order] = (
            np.ctypeslib.as_array(ql).copy(),
            np.ctypeslib.as_array(qb).copy(),
            np.ctypeslib.as_array(wl).copy(),
        )
    lib.seams_frame_free(fr)
    return out


def pyscal_params(pts, box):
    atoms = Atoms(positions=pts, cell=box, pbc=True)
    pyscal3.find_neighbors(atoms, method="cutoff", cutoff=CUTOFF)
    out = {}
    for order in (4, 6):
        local = np.asarray(pyscal3.steinhardt_parameter(atoms, l=order, averaged=False)[0], float)
        averaged = np.asarray(pyscal3.steinhardt_parameter(atoms, l=order, averaged=True)[0], float)
        invariant = np.asarray(
            pyscal3.wigner_w_parameter(atoms, l=order, averaged=False, normalized=True)[0],
            float,
        )
        out[order] = (local, averaged, invariant)
    return out


def max_abs(left, right):
    if left.shape != right.shape:
        raise AssertionError(f"shape {left.shape} vs {right.shape}")
    if np.any(~np.isfinite(left)) or np.any(~np.isfinite(right)):
        raise AssertionError("non-finite order parameter")
    return float(np.max(np.abs(left - right)))


def labels_of(lib, pts, box):
    fr, n = frame(lib, pts, box)
    q4 = (ctypes.c_double * n)()
    q6 = (ctypes.c_double * n)()
    b4 = (ctypes.c_double * n)()
    b6 = (ctypes.c_double * n)()
    labels = (ctypes.c_int * n)()
    counts = (ctypes.c_int * 9)()
    if lib.seams_frame_q4_q6_chill(fr, CUTOFF, 1, q4, q6, b4, b6, labels, counts) != 0:
        raise RuntimeError("fused call failed")
    lib.seams_frame_free(fr)
    return np.ctypeslib.as_array(labels).copy()


def main():
    lib = load()
    report = []
    for name, builder, label in (
        ("Ic", ice_ic, 0),
        ("Ih", ice_ih, 1),
    ):
        pts, box = builder()
        got = seams_params(lib, pts, box)
        ref = pyscal_params(pts, box)
        kind = labels_of(lib, pts, box)
        if not np.all(kind == label):
            raise SystemExit(f"{name} CHILL+ labels are not all {label}")
        for order in (4, 6):
            for piece, title in enumerate(("q", "qbar", "w")):
                err = max_abs(got[order][piece], ref[order][piece])
                report.append(f"{name} {title}{order} max {err:.3e} value {got[order][piece][0]:.8f}")
                if err > TOL:
                    raise SystemExit(report[-1])
    pts, box = liquid()
    got = seams_params(lib, pts, box)
    ref = pyscal_params(pts, box)
    for order in (4, 6):
        for piece, title in enumerate(("q", "qbar", "w")):
            err = max_abs(got[order][piece], ref[order][piece])
            report.append(f"liquid {title}{order} max {err:.3e}")
            if err > TOL:
                raise SystemExit(report[-1])
    for line in report:
        print(line)


if __name__ == "__main__":
    main()
