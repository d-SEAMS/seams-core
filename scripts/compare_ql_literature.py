#!/usr/bin/env python3
"""Print freud / pyscal3 / DScribe numbers on the same FCC and BCC cells
as tests/compare_structure_desc.cpp. Missing imports print MISSING.
"""

from __future__ import annotations

import importlib
import importlib.util
import sys

import numpy as np


def lattice(basis, reps: int, a: float):
    pos = []
    for i in range(reps):
        for j in range(reps):
            for k in range(reps):
                for b in basis:
                    pos.append(((i + b[0]) * a, (j + b[1]) * a, (k + b[2]) * a))
    xyz = np.asarray(pos, dtype=np.float64)
    L = reps * a
    return xyz, L


def fcc():
    return lattice(
        ((0.0, 0.0, 0.0), (0.5, 0.5, 0.0), (0.5, 0.0, 0.5), (0.0, 0.5, 0.5)),
        3,
        4.0,
    )


def bcc():
    return lattice(((0.0, 0.0, 0.0), (0.5, 0.5, 0.5)), 3, 4.0)


def have(name: str):
    return importlib.util.find_spec(name) is not None


def run_freud(name: str, xyz, L: float, r_max: float):
    import freud

    box = freud.box.Box.cube(L)
    print(f"tool=freud lattice={name} n={len(xyz)} L={L} cut_ql={r_max}")
    for ell in (4, 6, 8):
        op = freud.order.Steinhardt(l=ell)
        op.compute((box, xyz), neighbors={"r_max": r_max, "exclude_ii": True})
        ql = np.asarray(op.particle_order, dtype=np.float64)
        print(f"  ql{ell}={float(np.mean(ql)):.6f}")


def run_pyscal(name: str, xyz, L: float, r_max: float):
    from pyscal3 import System

    # Same cell as compare_structure_desc.cpp. atoms+box ctor fills ghosts.
    sys_ = System(
        atoms={"positions": xyz},
        box=[[L, 0.0, 0.0], [0.0, L, 0.0], [0.0, 0.0, L]],
    )
    sys_.find.neighbors(method="cutoff", cutoff=r_max)
    print(f"tool=pyscal3 lattice={name} n={len(xyz)} L={L} cut_ql={r_max}")
    qs = sys_.calculate.steinhardt_parameter([4, 6, 8])
    arr = np.asarray(qs)
    if arr.ndim == 2 and arr.shape[0] == 3:
        for ell, row in zip((4, 6, 8), arr):
            print(f"  ql{ell}={float(np.mean(row)):.6f}")
    elif arr.ndim == 2 and arr.shape[1] == 3:
        for ell, col in zip((4, 6, 8), arr.T):
            print(f"  ql{ell}={float(np.mean(col)):.6f}")
    else:
        print(f"  raw_shape={arr.shape} raw={arr!r}")


def run_dscribe(name: str, xyz, L: float, r_max: float):
    from dscribe.descriptors import SOAP

    soap = SOAP(
        species=["X"],
        r_cut=r_max,
        n_max=3,
        l_max=6,
        periodic=True,
        sparse=False,
    )
    # DScribe wants an ASE Atoms
    from ase import Atoms

    atoms = Atoms(symbols=["X"] * len(xyz), positions=xyz, cell=[L, L, L], pbc=True)
    spec = soap.create(atoms)
    print(
        f"tool=dscribe lattice={name} n={len(xyz)} L={L} soap_shape={tuple(spec.shape)}"
    )
    print(f"  soap_norm0={float(np.linalg.norm(spec[0])):.6f}")


def run_ptm(name: str, xyz, L: float):
    print(f"tool=ptm lattice={name} n={len(xyz)} L={L} status=MISSING")


def main():
    cases = (
        ("fcc", fcc(), 3.2),
        ("bcc", bcc(), 4.0),
    )
    print("python", sys.version.replace("\n", " "))
    for pkg in ("freud", "pyscal3", "dscribe", "ovito", "featomic"):
        spec = importlib.util.find_spec(pkg)
        ver = "?"
        if spec is not None:
            try:
                mod = importlib.import_module(pkg)
                ver = getattr(mod, "__version__", "?")
            except Exception as exc:  # noqa: BLE001
                ver = f"import_failed:{exc}"
        print(f"pkg {pkg} present={spec is not None} version={ver}")

    for name, (xyz, L), cut in cases:
        if have("freud"):
            try:
                run_freud(name, xyz, L, cut)
            except Exception as exc:  # noqa: BLE001
                print(f"tool=freud lattice={name} error={exc}")
        else:
            print(f"tool=freud lattice={name} status=MISSING")
        if have("pyscal3"):
            try:
                run_pyscal(name, xyz, L, cut)
            except Exception as exc:  # noqa: BLE001
                print(f"tool=pyscal3 lattice={name} error={exc}")
        else:
            print(f"tool=pyscal3 lattice={name} status=MISSING")
        if have("dscribe"):
            try:
                run_dscribe(name, xyz, L, cut)
            except Exception as exc:  # noqa: BLE001
                print(f"tool=dscribe lattice={name} error={exc}")
        else:
            print(f"tool=dscribe lattice={name} status=MISSING")
        run_ptm(name, xyz, L)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
