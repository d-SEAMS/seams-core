"""Per-atom q4 and q6 against OVITO and freud.

Ideal ice Ic and one noisy liquid, at one thread and at eight. The
library call leaves the averaged buffer null, so this is the local
sum. OVITO 3.17 is float64 and must agree to 1e-10. freud 3.6 stores
q_l as float32: on ice that is a few ulps, and on a disordered shell
the float32 sum sits near 1e-5. An empty shell is NaN in freud and 0
here and in OVITO.
"""

from __future__ import annotations

import ctypes
import math
import os
import subprocess
import sys

try:
    import numpy as np
except ImportError:
    sys.stderr.write("numpy is not installed\n")
    raise SystemExit(77)

try:
    from ovito.data import DataCollection, Particles, SimulationCell
    from ovito.modifiers import SteinhardtParametersModifier
    from ovito.pipeline import Pipeline, StaticSource
    import freud
except ImportError:
    sys.stderr.write("ovito or freud is not installed\n")
    raise SystemExit(77)

LIB = os.environ.get("SEAMS_LANG_LIB")
if not LIB:
    sys.stderr.write("SEAMS_LANG_LIB is not set\n")
    raise SystemExit(77)

CUTOFF = 3.5
OVITO_TOL = 1e-10
# A few units in the last place of freud's float32 q_l.
FREUD_ULP = 3


def ice_ic(cells, bond=2.75):
    lattice = 4.0 * bond / math.sqrt(3.0)
    side = cells * lattice
    fcc = np.array([[0.0, 0.0, 0.0], [0.5, 0.5, 0.0], [0.5, 0.0, 0.5], [0.0, 0.5, 0.5]])
    basis = np.vstack([fcc, fcc + 0.25])
    ijk = np.stack(np.meshgrid(*[np.arange(cells)] * 3, indexing="ij"), -1).reshape(-1, 3)
    pts = (ijk[:, None, :] + basis).reshape(-1, 3) * lattice
    return np.ascontiguousarray(pts % side), float(side)


def liquid(cells):
    pts, side = ice_ic(cells)
    rng = np.random.default_rng(91001)
    return np.ascontiguousarray(np.mod(pts + rng.normal(0.0, 0.9, pts.shape), side)), side


def seams_ql(lib, pts, side, order, averaged):
    n = len(pts)
    xyz = np.ascontiguousarray(pts, dtype=np.float64).reshape(-1)
    box = (ctypes.c_double * 9)(0, 0, 0, side, side, side, 0, 0, 0)
    types = (ctypes.c_int * n)(*([1] * n))
    ids = (ctypes.c_int * n)(*range(1, n + 1))
    frame = lib.seams_frame_from_arrays(
        xyz.ctypes.data_as(ctypes.POINTER(ctypes.c_double)), n, box, types, ids)
    if not frame:
        raise RuntimeError("frame")
    ql = (ctypes.c_double * n)()
    bar = None if not averaged else (ctypes.c_double * n)()
    try:
        if lib.seams_frame_steinhardt(frame, CUTOFF, 1, order, ql, bar) != 0:
            raise RuntimeError("steinhardt")
        return np.ctypeslib.as_array(ql).copy()
    finally:
        lib.seams_frame_free(frame)


def ovito_ql(pts, side, order):
    data = DataCollection()
    cell = SimulationCell(pbc=(True, True, True))
    cell[:, :3] = np.diag([side, side, side])
    data.objects.append(cell)
    part = Particles()
    part.create_property("Position", data=pts)
    data.objects.append(part)
    pipe = Pipeline(source=StaticSource(data=data))
    mod = SteinhardtParametersModifier(cutoff=CUTOFF, max_l=order, compute_averaged=False)
    mod.neighbor_mode = SteinhardtParametersModifier.NeighborMode.Cutoff
    pipe.modifiers.append(mod)
    out = pipe.compute()
    arr = np.asarray(out.particles["Steinhardt"], dtype=np.float64)
    return np.ascontiguousarray(arr[:, -1])


def freud_ql(pts, side, order):
    freud.parallel.set_num_threads(int(os.environ.get("OMP_NUM_THREADS", "1")))
    box = freud.box.Box.cube(side)
    op = freud.order.Steinhardt(l=order, average=False)
    op.compute((box, pts), neighbors={"r_max": CUTOFF})
    return np.asarray(op.ql, dtype=np.float64).copy()


def check(lib):
    frames = (("ice",) + ice_ic(7), ("liquid",) + liquid(7))
    for name, pts, side in frames:
        for order in (4, 6):
            local = seams_ql(lib, pts, side, order, averaged=False)
            full = seams_ql(lib, pts, side, order, averaged=True)
            ov = ovito_ql(pts, side, order)
            fr = freud_ql(pts, side, order)
            d_ov = float(np.max(np.abs(local - ov)))
            d_full = float(np.max(np.abs(local - full)))
            if d_ov > OVITO_TOL or d_full > OVITO_TOL:
                raise SystemExit(
                    f"{name} l={order} ovito {d_ov:.3e} local-vs-bar {d_full:.3e}")
            # freud's q_l is float32, and an empty shell is NaN there.
            # OVITO and this library write 0 for that shell.
            empty = ~np.isfinite(fr)
            if np.any(empty) and float(np.max(np.abs(local[empty]))) != 0.0:
                raise SystemExit(f"{name} l={order} empty shell is not 0")
            if np.any(empty) and float(np.max(np.abs(ov[empty]))) > OVITO_TOL:
                raise SystemExit(f"{name} l={order} ovito empty shell")
            live = ~empty
            if name == "ice":
                rounded = local[live].astype(np.float32)
                ulp = np.spacing(rounded).astype(np.float64)
                steps = float(np.max(np.abs(fr[live] - rounded.astype(np.float64)) / ulp))
                if steps > FREUD_ULP:
                    raise SystemExit(f"{name} l={order} freud {steps:.2f} ulp")
                freud_note = f"freud_ulp {steps:.2f}"
            else:
                # freud accumulates q_l in float32. On a disordered shell that
                # residual is about 1e-5, which is not a 1e-10 comparison.
                d_fr = float(np.max(np.abs(local[live] - fr[live])))
                if d_fr > 2e-5:
                    raise SystemExit(f"{name} l={order} freud {d_fr:.3e}")
                freud_note = f"freud {d_fr:.3e}"
            print(f"{name} n={len(pts)} l={order} threads={os.environ.get('OMP_NUM_THREADS')} "
                  f"ovito {d_ov:.3e} {freud_note} empty {int(empty.sum())}")


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
            env["OVITO_THREAD_COUNT"] = threads
            env["QL_EQUAL_CHILD"] = "1"
            done = subprocess.run([sys.executable, __file__], env=env)
            if done.returncode == 77:
                raise SystemExit(77)
            if done.returncode != 0:
                raise SystemExit(done.returncode)
        return
    check(load())


if __name__ == "__main__":
    main()
