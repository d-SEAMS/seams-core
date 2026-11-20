# Copyright (c) d-SEAMS developers
"""ctypes bridge onto libseams for trajectory and time-correlation tools."""

from __future__ import annotations

import ctypes
import ctypes.util
import os
from pathlib import Path

import numpy as np

_REPO_ROOT = Path(__file__).resolve().parents[3]
_LIB_CACHE: ctypes.CDLL | None = None

SEAMS_OK = 0
SEAMS_ERR_NULL = 1
SEAMS_ERR_RANGE = 2
SEAMS_ERR_BOX = 3
SEAMS_ERR_ARG = 4
SEAMS_ERR_STATE = 5
SEAMS_ERR_NOMEM = 6
SEAMS_ERR_INTERNAL = 7


def _find_libseams() -> str:
    env_lib = os.environ.get("SEAMS_LIB")
    if env_lib and os.path.isfile(env_lib):
        return env_lib

    candidates: list[Path] = []
    for sub in ("bbdir-head", "bbdir", "build", "bbdir-develop"):
        for name in ("libseams.so", "libseams.dylib", "seams.dll"):
            candidates.append(_REPO_ROOT / sub / "src" / name)
    if _REPO_ROOT.is_dir():
        for folder in sorted(_REPO_ROOT.glob("bbdir*")):
            for name in ("libseams.so", "libseams.dylib", "seams.dll"):
                candidates.append(folder / "src" / name)
    conda = os.environ.get("CONDA_PREFIX")
    if conda:
        for name in ("libseams.so", "libseams.dylib", "seams.dll"):
            candidates.append(Path(conda) / "lib" / name)
    for cand in candidates:
        if cand.is_file():
            return str(cand)

    found = ctypes.util.find_library("seams")
    if found:
        return found
    return env_lib or "libseams.so"


def get_lib() -> ctypes.CDLL:
    """Return the configured libseams CDLL singleton."""
    global _LIB_CACHE
    if _LIB_CACHE is not None:
        return _LIB_CACHE

    lib = ctypes.CDLL(_find_libseams())
    c_int = ctypes.c_int
    c_size_t = ctypes.c_size_t
    c_double = ctypes.c_double
    c_void_p = ctypes.c_void_p
    c_char_p = ctypes.c_char_p
    p_int = ctypes.POINTER(c_int)
    p_size_t = ctypes.POINTER(c_size_t)
    p_double = ctypes.POINTER(c_double)
    p_void_p = ctypes.POINTER(c_void_p)

    lib.seams_last_error_message.argtypes = []
    lib.seams_last_error_message.restype = c_char_p
    lib.seams_handle_error.argtypes = [c_void_p]
    lib.seams_handle_error.restype = c_char_p
    lib.seams_status_string.argtypes = [c_int]
    lib.seams_status_string.restype = c_char_p

    lib.seams_path_nframes.argtypes = [c_char_p, p_size_t]
    lib.seams_path_nframes.restype = c_int

    lib.seams_system_create.argtypes = []
    lib.seams_system_create.restype = c_void_p
    lib.seams_system_destroy.argtypes = [c_void_p]
    lib.seams_system_destroy.restype = None
    lib.seams_system_read.argtypes = [c_void_p, c_char_p, c_int, c_int]
    lib.seams_system_read.restype = c_int
    lib.seams_system_set_atoms.argtypes = [
        c_void_p,
        p_double,
        c_size_t,
        p_double,
        c_size_t,
        p_int,
        p_int,
    ]
    lib.seams_system_set_atoms.restype = c_int
    lib.seams_system_natoms.argtypes = [c_void_p, p_size_t]
    lib.seams_system_natoms.restype = c_int
    lib.seams_system_box.argtypes = [c_void_p, p_double, c_size_t]
    lib.seams_system_box.restype = c_int
    lib.seams_system_positions.argtypes = [c_void_p, p_double, c_size_t]
    lib.seams_system_positions.restype = c_int
    lib.seams_system_types.argtypes = [c_void_p, p_int, c_size_t]
    lib.seams_system_types.restype = c_int
    lib.seams_system_ids.argtypes = [c_void_p, p_int, c_size_t]
    lib.seams_system_ids.restype = c_int
    lib.seams_chill_plus_n.argtypes = [
        c_void_p,
        c_double,
        c_int,
        c_int,
        c_int,
        p_int,
        c_size_t,
        p_int,
    ]
    lib.seams_chill_plus_n.restype = c_int

    lib.seams_mols.argtypes = [c_void_p, p_int, c_size_t]
    lib.seams_mols.restype = c_int
    lib.seams_set_mols.argtypes = [c_void_p, p_int, c_size_t]
    lib.seams_set_mols.restype = c_int

    lib.seams_molecule_centres.argtypes = [
        c_void_p,
        p_int,
        c_size_t,
        p_double,
        c_size_t,
        p_int,
        p_double,
        c_size_t,
        p_size_t,
    ]
    lib.seams_molecule_centres.restype = c_int

    lib.seams_traj_unwrap_positions.argtypes = [
        p_double,
        c_size_t,
        c_size_t,
        p_double,
        c_size_t,
        p_double,
        c_size_t,
    ]
    lib.seams_traj_unwrap_positions.restype = c_int

    lib.seams_traj_unwrap_centres.argtypes = [
        p_double,
        c_size_t,
        c_size_t,
        p_double,
        c_size_t,
        p_int,
        c_size_t,
        p_double,
        c_size_t,
        p_int,
        c_size_t,
        p_double,
        c_size_t,
        p_size_t,
    ]
    lib.seams_traj_unwrap_centres.restype = c_int

    lib.seams_water_vectors.argtypes = [
        c_void_p,
        c_int,
        c_int,
        p_int,
        p_double,
        p_double,
        p_double,
        p_double,
        c_size_t,
        p_size_t,
    ]
    lib.seams_water_vectors.restype = c_int

    lib.seams_molecule_dipoles.argtypes = [
        c_void_p,
        p_int,
        c_size_t,
        p_double,
        c_size_t,
        p_int,
        p_double,
        p_double,
        c_size_t,
        p_size_t,
    ]
    lib.seams_molecule_dipoles.restype = c_int

    lib.seams_orientation_order.argtypes = [
        p_double,
        c_size_t,
        p_double,
        p_double,
        p_double,
        c_size_t,
        p_double,
        p_double,
    ]
    lib.seams_orientation_order.restype = c_int

    lib.seams_msd_create.argtypes = [
        c_size_t,
        c_size_t,
        c_size_t,
        c_double,
        p_double,
        p_void_p,
    ]
    lib.seams_msd_create.restype = c_int
    lib.seams_msd_destroy.argtypes = [c_void_p]
    lib.seams_msd_destroy.restype = None
    lib.seams_msd_push.argtypes = [c_void_p, p_double, c_size_t, p_double, c_size_t]
    lib.seams_msd_push.restype = c_int
    lib.seams_msd_curve.argtypes = [
        c_void_p,
        p_double,
        p_double,
        p_double,
        p_double,
        p_size_t,
        c_size_t,
        p_size_t,
    ]
    lib.seams_msd_curve.restype = c_int
    lib.seams_msd_fit.argtypes = [
        c_void_p,
        c_double,
        c_double,
        p_double,
        p_double,
        p_double,
        p_double,
    ]
    lib.seams_msd_fit.restype = c_int

    lib.seams_hbond_tcf_create.argtypes = [
        c_size_t,
        c_size_t,
        c_double,
        c_int,
        p_void_p,
    ]
    lib.seams_hbond_tcf_create.restype = c_int
    lib.seams_hbond_tcf_destroy.argtypes = [c_void_p]
    lib.seams_hbond_tcf_destroy.restype = None
    lib.seams_hbond_tcf_push.argtypes = [c_void_p, p_int, p_int, c_size_t]
    lib.seams_hbond_tcf_push.restype = c_int
    lib.seams_hbond_tcf_curve.argtypes = [
        c_void_p,
        p_double,
        p_double,
        p_double,
        p_size_t,
        c_size_t,
        p_size_t,
    ]
    lib.seams_hbond_tcf_curve.restype = c_int
    lib.seams_hbond_tcf_fit.argtypes = [
        c_void_p,
        c_double,
        c_double,
        p_double,
        p_double,
        p_double,
        p_double,
    ]
    lib.seams_hbond_tcf_fit.restype = c_int

    lib.seams_reorient_create.argtypes = [
        c_size_t,
        c_size_t,
        c_size_t,
        c_double,
        p_void_p,
    ]
    lib.seams_reorient_create.restype = c_int
    lib.seams_reorient_destroy.argtypes = [c_void_p]
    lib.seams_reorient_destroy.restype = None
    lib.seams_reorient_push.argtypes = [c_void_p, p_double, c_size_t]
    lib.seams_reorient_push.restype = c_int
    lib.seams_reorient_curve.argtypes = [
        c_void_p,
        p_double,
        p_double,
        p_double,
        p_size_t,
        c_size_t,
        p_size_t,
    ]
    lib.seams_reorient_curve.restype = c_int
    lib.seams_reorient_fit.argtypes = [
        c_void_p,
        c_double,
        c_double,
        p_double,
        p_double,
        p_double,
        p_double,
    ]
    lib.seams_reorient_fit.restype = c_int

    _LIB_CACHE = lib
    return lib


def check_status(status: int, handle: int | ctypes.c_void_p | None = None) -> None:
    """Raise an informative Python exception when status != SEAMS_OK."""
    if status == SEAMS_OK:
        return
    lib = get_lib()
    raw_msg = (
        lib.seams_handle_error(handle)
        if handle is not None
        else lib.seams_last_error_message()
    )
    if not raw_msg:
        raw_msg = lib.seams_status_string(status)
    text = raw_msg.decode("utf-8", errors="replace") if raw_msg else f"seams error {status}"
    if status == SEAMS_ERR_NOMEM:
        raise MemoryError(text)
    if status in (SEAMS_ERR_RANGE, SEAMS_ERR_BOX, SEAMS_ERR_ARG, SEAMS_ERR_NULL):
        raise ValueError(text)
    raise RuntimeError(text)


def double_ptr(arr: np.ndarray):
    return arr.ctypes.data_as(ctypes.POINTER(ctypes.c_double))


def int_ptr(arr: np.ndarray):
    return arr.ctypes.data_as(ctypes.POINTER(ctypes.c_int))


def sizet_ptr(arr: np.ndarray):
    return arr.ctypes.data_as(ctypes.POINTER(ctypes.c_size_t))


def normalize_box9(cell, box_low=None) -> np.ndarray:
    """Convert a (3,), (6,), (9,), or (3, 3) periodic box specification to box9.

    Returns ``[ox, oy, oz, lx, ly, lz, xy, xz, yz]`` as float64.
    A 6-element vector follows the d-SEAMS ``PointCloud.box`` / ``Frame.from_arrays``
    convention ``[xspan, yspan, zspan, xy, xz, yz]``.
    A (3, 3) array is treated as row lattice vectors ``[a; b; c]`` in the
    restricted-triclinic orientation ``a = (lx, 0, 0)``, ``b = (xy, ly, 0)``,
    ``c = (xz, yz, lz)``.
    """
    arr = np.asarray(cell, dtype=np.float64)
    ox, oy, oz = (0.0, 0.0, 0.0) if box_low is None else map(float, box_low)
    if arr.shape == (9,):
        out = np.array(arr, dtype=np.float64, copy=True)
        if box_low is not None:
            out[:3] = (ox, oy, oz)
        return out
    if arr.shape == (3,):
        return np.array(
            [ox, oy, oz, arr[0], arr[1], arr[2], 0.0, 0.0, 0.0], dtype=np.float64
        )
    if arr.shape == (6,):
        xspan, yspan, zspan, xy, xz, yz = map(float, arr)
        xmin = min(0.0, xy, xz, xy + xz)
        xmax = max(0.0, xy, xz, xy + xz)
        ymin = min(0.0, yz)
        ymax = max(0.0, yz)
        lx = xspan - xmax + xmin
        ly = yspan - ymax + ymin
        lz = zspan
        return np.array(
            [ox - xmin, oy - ymin, oz, lx, ly, lz, xy, xz, yz], dtype=np.float64
        )
    if arr.shape == (3, 3):
        if (
            abs(arr[0, 1]) > 1e-10
            or abs(arr[0, 2]) > 1e-10
            or abs(arr[1, 2]) > 1e-10
        ):
            raise ValueError(
                "3x3 cell matrix must be restricted triclinic lower-triangular "
                "[[lx, 0, 0], [xy, ly, 0], [xz, yz, lz]]"
            )
        return np.array(
            [
                ox,
                oy,
                oz,
                arr[0, 0],
                arr[1, 1],
                arr[2, 2],
                arr[1, 0],
                arr[2, 0],
                arr[2, 1],
            ],
            dtype=np.float64,
        )
    raise ValueError(
        f"cell must have shape (3,), (6,), (9,), or (3, 3); got {arr.shape}"
    )


def box9_to_dump_bounds(box9: np.ndarray) -> tuple[list[float], list[float]]:
    """Convert box9 to PointCloud (box, boxLow) vectors."""
    ox, oy, oz, lx, ly, lz, xy, xz, yz = map(float, box9)
    if abs(xy) == 0.0 and abs(xz) == 0.0 and abs(yz) == 0.0:
        return [lx, ly, lz], [ox, oy, oz]
    xmin = min(0.0, xy, xz, xy + xz)
    xmax = max(0.0, xy, xz, xy + xz)
    ymin = min(0.0, yz)
    ymax = max(0.0, yz)
    return [lx + xmax - xmin, ly + ymax - ymin, lz, xy, xz, yz], [
        ox + xmin,
        oy + ymin,
        oz,
    ]
