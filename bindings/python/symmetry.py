# Copyright (c) d-SEAMS developers
"""Local symmetry from libyodaLib, loaded with ctypes.

Set ``SEAMS_CORE_LIB`` to the built ``libyodaLib``. Template ids are
``ice-Ih``, ``ice-Ic``, ``ice-Isd``, ``5^12``, ``5^12 6^2``,
``5^12 6^4``, ``ring-5`` and ``ring-6``.
"""

from __future__ import annotations

import ctypes
import os

import numpy as np

_NAMES = (
    "ice-Ih",
    "ice-Ic",
    "ice-Isd",
    "5^12",
    "5^12 6^2",
    "5^12 6^4",
    "ring-5",
    "ring-6",
)


def _library():
    path = os.environ.get("SEAMS_CORE_LIB", "libyodaLib")
    lib = ctypes.CDLL(path)
    lib.seams_symmetry_match.argtypes = [
        ctypes.POINTER(ctypes.c_double),
        ctypes.c_int,
        ctypes.POINTER(ctypes.c_int),
        ctypes.POINTER(ctypes.c_double),
        ctypes.POINTER(ctypes.c_double),
        ctypes.POINTER(ctypes.c_double),
        ctypes.POINTER(ctypes.c_double),
        ctypes.POINTER(ctypes.c_int),
    ]
    lib.seams_symmetry_match.restype = ctypes.c_int
    lib.seams_symmetry_point_group.argtypes = [
        ctypes.POINTER(ctypes.c_double),
        ctypes.c_int,
        ctypes.c_double,
        ctypes.c_char_p,
        ctypes.c_int,
        ctypes.POINTER(ctypes.c_int),
    ]
    lib.seams_symmetry_point_group.restype = ctypes.c_int
    lib.seams_symmetry_space_group.argtypes = [
        ctypes.POINTER(ctypes.c_double),
        ctypes.POINTER(ctypes.c_double),
        ctypes.c_int,
        ctypes.c_double,
        ctypes.POINTER(ctypes.c_int),
        ctypes.c_char_p,
        ctypes.c_int,
    ]
    lib.seams_symmetry_space_group.restype = ctypes.c_int
    lib.seams_symmetry_template.argtypes = [
        ctypes.c_int,
        ctypes.POINTER(ctypes.c_double),
        ctypes.c_int,
        ctypes.POINTER(ctypes.c_int),
    ]
    lib.seams_symmetry_template.restype = ctypes.c_int
    lib.seams_symmetry_crystal.argtypes = [
        ctypes.c_int,
        ctypes.POINTER(ctypes.c_double),
        ctypes.POINTER(ctypes.c_double),
        ctypes.c_int,
        ctypes.POINTER(ctypes.c_int),
    ]
    lib.seams_symmetry_crystal.restype = ctypes.c_int
    lib.seams_symmetry_backend.argtypes = [ctypes.c_char_p, ctypes.c_int]
    lib.seams_symmetry_backend.restype = ctypes.c_int
    lib.seams_symmetry_classify_tol.argtypes = [
        ctypes.POINTER(ctypes.c_double),
        ctypes.c_int,
        ctypes.POINTER(ctypes.c_double),
        ctypes.c_double,
        ctypes.POINTER(ctypes.c_int),
    ]
    lib.seams_symmetry_classify_tol.restype = ctypes.c_int
    return lib


def backend():
    lib = _library()
    buf = ctypes.create_string_buffer(32)
    lib.seams_symmetry_backend(buf, 32)
    return buf.value.decode()


def match(points):
    """Best template for an ``n x 3`` environment.

    The rotation and translation take the permuted environment onto the
    template.
    """
    env = np.ascontiguousarray(points, dtype=np.float64)
    if env.ndim != 2 or env.shape[1] != 3:
        raise ValueError("points must be n x 3")
    n = env.shape[0]
    ident = ctypes.c_int(-1)
    rmsd = ctypes.c_double(0)
    haus = ctypes.c_double(0)
    rotation = np.zeros(9, dtype=np.float64)
    translation = np.zeros(3, dtype=np.float64)
    perm = np.empty(n, dtype=np.int32)
    lib = _library()
    code = lib.seams_symmetry_match(
        env.ctypes.data_as(ctypes.POINTER(ctypes.c_double)),
        ctypes.c_int(n),
        ctypes.byref(ident),
        ctypes.byref(rmsd),
        ctypes.byref(haus),
        rotation.ctypes.data_as(ctypes.POINTER(ctypes.c_double)),
        translation.ctypes.data_as(ctypes.POINTER(ctypes.c_double)),
        perm.ctypes.data_as(ctypes.POINTER(ctypes.c_int)),
    )
    if code != 0:
        raise RuntimeError("symmetry match failed")
    return {
        "name": _NAMES[ident.value],
        "rmsd": rmsd.value,
        "hausdorff": haus.value,
        "rotation": rotation.reshape(3, 3),
        "translation": translation,
        "permutation": perm,
    }


def point_group(points, sym_thr=0.05):
    env = np.ascontiguousarray(points, dtype=np.float64)
    n = env.shape[0]
    symbol = ctypes.create_string_buffer(16)
    nops = ctypes.c_int(0)
    lib = _library()
    code = lib.seams_symmetry_point_group(
        env.ctypes.data_as(ctypes.POINTER(ctypes.c_double)),
        ctypes.c_int(n),
        ctypes.c_double(sym_thr),
        symbol,
        16,
        ctypes.byref(nops),
    )
    if code == 2:
        return None
    if code != 0:
        raise RuntimeError("point group failed")
    return symbol.value.decode(), nops.value


def space_group(lattice, cartesian, symprec=1e-3):
    """``lattice`` is three row vectors. ``cartesian`` is ``n x 3``."""
    vec = np.ascontiguousarray(lattice, dtype=np.float64).reshape(9)
    xyz = np.ascontiguousarray(cartesian, dtype=np.float64)
    number = ctypes.c_int(0)
    symbol = ctypes.create_string_buffer(16)
    lib = _library()
    code = lib.seams_symmetry_space_group(
        vec.ctypes.data_as(ctypes.POINTER(ctypes.c_double)),
        xyz.ctypes.data_as(ctypes.POINTER(ctypes.c_double)),
        ctypes.c_int(xyz.shape[0]),
        ctypes.c_double(symprec),
        ctypes.byref(number),
        symbol,
        16,
    )
    if code == 2:
        return None
    if code != 0:
        raise RuntimeError("space group failed")
    return number.value, symbol.value.decode()


def template(ident):
    if isinstance(ident, str):
        ident = _NAMES.index(ident)
    lib = _library()
    count = ctypes.c_int(0)
    if lib.seams_symmetry_template(ctypes.c_int(ident), None, 0, ctypes.byref(count)) != 0:
        raise RuntimeError("unknown template")
    buf = np.empty((count.value, 3), dtype=np.float64)
    code = lib.seams_symmetry_template(
        ctypes.c_int(ident),
        buf.ctypes.data_as(ctypes.POINTER(ctypes.c_double)),
        count.value,
        ctypes.byref(count),
    )
    if code != 0:
        raise RuntimeError("template copy failed")
    return buf


def classify(cartesian, lattice, shell_rmsd=0.12):
    """Per-atom ice label. 0 ice Ih, 1 ice Ic, 2 ice Isd, -1 otherwise.

    ``shell_rmsd`` is the scale-free ceiling on the second-shell residual
    after the optimal rotation and scale. Zero is a perfect shape.
    """
    xyz = np.ascontiguousarray(cartesian, dtype=np.float64)
    vec = np.ascontiguousarray(lattice, dtype=np.float64).reshape(9)
    labels = np.empty(xyz.shape[0], dtype=np.int32)
    code = _library().seams_symmetry_classify_tol(
        xyz.ctypes.data_as(ctypes.POINTER(ctypes.c_double)),
        ctypes.c_int(xyz.shape[0]),
        vec.ctypes.data_as(ctypes.POINTER(ctypes.c_double)),
        ctypes.c_double(shell_rmsd),
        labels.ctypes.data_as(ctypes.POINTER(ctypes.c_int)),
    )
    if code != 0:
        raise RuntimeError("classify failed")
    return labels


def crystal(which):
    """``which`` is ``ih``, ``ic`` or ``fault``. Returns lattice rows and coordinates."""
    index = {"ih": 0, "ic": 1, "fault": 2}[which]
    lib = _library()
    count = ctypes.c_int(0)
    lattice = np.empty(9, dtype=np.float64)
    if lib.seams_symmetry_crystal(index, lattice.ctypes.data_as(ctypes.POINTER(ctypes.c_double)), None, 0, ctypes.byref(count)) != 0:
        raise RuntimeError("crystal query failed")
    xyz = np.empty((count.value, 3), dtype=np.float64)
    code = lib.seams_symmetry_crystal(
        index,
        lattice.ctypes.data_as(ctypes.POINTER(ctypes.c_double)),
        xyz.ctypes.data_as(ctypes.POINTER(ctypes.c_double)),
        count.value,
        ctypes.byref(count),
    )
    if code != 0:
        raise RuntimeError("crystal copy failed")
    return lattice.reshape(3, 3), xyz
