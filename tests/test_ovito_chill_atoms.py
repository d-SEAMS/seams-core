#!/usr/bin/env python3
"""OVITO-rule CHILL+ against OVITO 3.17 ChillPlusModifier, every atom.

Exit 77 when OVITO cannot be imported. The published four-neighbour
rule is not this modifier; the ovito rule is.
"""

from __future__ import annotations

import ctypes
import os
import sys
from pathlib import Path

import numpy as np

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "scripts"))

CUTOFF = 3.5
# d-SEAMS label -> OVITO ChillPlus structure type
TO_OVITO = np.array([2, 1, 0, 3, 4, 5, 0, 0, 0], dtype=np.int64)


def prepare_ovito() -> bool:
    site = os.environ.get("SEAMS_OVITO_SITE", "")
    if site and os.path.isdir(site) and site not in sys.path:
        sys.path.insert(0, site)
    extra = []
    gl = "/usr/lib/x86_64-linux-gnu"
    qt = os.path.join(site, "PySide6", "Qt", "lib")
    for path in (gl, qt):
        if os.path.isdir(path):
            extra.append(path)
    if extra and os.environ.get("SEAMS_OVITO_LIBS") != "1":
        current = os.environ.get("LD_LIBRARY_PATH", "")
        os.environ["LD_LIBRARY_PATH"] = ":".join(extra + ([current] if current else []))
        os.environ["SEAMS_OVITO_LIBS"] = "1"
        os.environ.setdefault("QT_QPA_PLATFORM", "offscreen")
        os.execv(sys.executable, [sys.executable, *sys.argv])
    os.environ.setdefault("QT_QPA_PLATFORM", "offscreen")
    try:
        import ovito  # noqa: F401

        return True
    except ImportError:
        return False


def load_seams(path: Path):
    yoda = path.parent / "libyodaLib.so"
    if yoda.is_file():
        ctypes.CDLL(str(yoda), mode=ctypes.RTLD_GLOBAL)
    lib = ctypes.CDLL(str(path))
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
    lib.seams_set_num_threads.argtypes = [ctypes.c_int]
    lib.seams_set_num_threads.restype = ctypes.c_int
    lib.seams_chill_plus.argtypes = [
        ctypes.c_void_p, ctypes.c_double, ctypes.c_int, ctypes.c_int, ctypes.c_int,
        ctypes.POINTER(ctypes.c_int), ctypes.POINTER(ctypes.c_int),
    ]
    lib.seams_chill_plus.restype = ctypes.c_int
    lib.seams_last_error.restype = ctypes.c_char_p
    return lib


def frame(lib, pos, lengths):
    n = len(pos)
    xyz = np.ascontiguousarray(pos, dtype=np.float64)
    box = (ctypes.c_double * 9)(0, 0, 0, lengths[0], lengths[1], lengths[2], 0, 0, 0)
    types = (ctypes.c_int * n)(*([1] * n))
    ids = (ctypes.c_int * n)(*range(1, n + 1))
    sys = lib.seams_system_create()
    rc = lib.seams_system_set_frame(
        sys, xyz.ctypes.data_as(ctypes.POINTER(ctypes.c_double)), n, box, types, ids
    )
    if rc != 0:
        err = lib.seams_last_error()
        raise SystemExit(err.decode() if err else "set_frame")
    return sys


def dseams_ovito(lib, pos, lengths) -> np.ndarray:
    sys = frame(lib, pos, lengths)
    n = len(pos)
    labels = (ctypes.c_int * n)()
    counts = (ctypes.c_int * 9)()
    rc = lib.seams_chill_plus(sys, CUTOFF, 1, 0, 1, labels, counts)
    if rc != 0:
        err = lib.seams_last_error()
        raise SystemExit(err.decode() if err else "chill+")
    out = np.ctypeslib.as_array(labels).copy()
    lib.seams_system_destroy(sys)
    return out


def ovito_chill(pos, lengths) -> np.ndarray:
    from ovito.data import DataCollection
    from ovito.modifiers import ChillPlusModifier
    from ovito.pipeline import Pipeline, StaticSource

    data = DataCollection()
    data.create_cell(np.diag(lengths), pbc=(True, True, True))
    parts = data.create_particles(count=len(pos))
    parts.create_property("Position", data=np.ascontiguousarray(pos, dtype=np.float64))
    pipe = Pipeline(source=StaticSource(data=data))
    pipe.modifiers.append(ChillPlusModifier(cutoff=CUTOFF))
    return np.asarray(pipe.compute().particles["Structure Type"])


def check(name, lib, pos, lengths) -> None:
    ours = TO_OVITO[dseams_ovito(lib, pos, lengths)]
    theirs = ovito_chill(pos, lengths)
    miss = int(np.count_nonzero(ours != theirs))
    print(f"{name} atoms {len(pos)} mismatches {miss}")
    if miss:
        raise SystemExit(f"{name} disagrees with OVITO on {miss} atoms")


def main() -> None:
    if not prepare_ovito():
        print("ovito is not importable", file=sys.stderr)
        raise SystemExit(77)
    if len(sys.argv) < 2:
        print("usage: test_ovito_chill_atoms.py libseams.so", file=sys.stderr)
        raise SystemExit(2)
    lib = load_seams(Path(sys.argv[1]))
    if lib.seams_set_num_threads(1) != 0:
        raise SystemExit("threads")
    from gen_ice_ic import ice_ic
    from h2h_freud_ovito import ice_positions, mixed_frames

    pos, lengths = ice_positions(4)
    check("Ih512", lib, pos, lengths)
    pts, side = ice_ic(4)
    check("Ic512", lib, np.ascontiguousarray(pts), (side, side, side))
    for name, (pos, lengths) in mixed_frames():
        check(name, lib, pos, lengths)
    print("ovito-rule CHILL+ matches OVITO 3.17 on every atom")


if __name__ == "__main__":
    main()
