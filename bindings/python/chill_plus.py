# Copyright (c) d-SEAMS developers
"""CHILL+ from libseams_lang, loaded with ctypes.

``rule`` is ``chill+``, ``ovito-chill+``, or the text of a document
that contains a ``rule`` key. A filesystem path is not opened here.
``seams chill-plus --rule FILE`` is the command that reads a file.
Set ``SEAMS_LANG_LIB`` to the built ``libseams_lang``.
"""

from __future__ import annotations

import ctypes
import os

import numpy as np

_NAMES = (
    "cubic",
    "hexagonal",
    "water",
    "interfacial",
    "clathrate",
    "interClathrate",
    "unclassified",
    "reCubic",
    "reHex",
)


def _library():
    path = os.environ.get("SEAMS_LANG_LIB", "libseams_lang")
    lib = ctypes.CDLL(path)
    lib.seams_frame_from_arrays.restype = ctypes.c_void_p
    lib.seams_frame_from_arrays.argtypes = [
        ctypes.POINTER(ctypes.c_double),
        ctypes.c_int,
        ctypes.POINTER(ctypes.c_double),
        ctypes.POINTER(ctypes.c_int),
        ctypes.POINTER(ctypes.c_int),
    ]
    lib.seams_frame_chill_plus_rule.restype = ctypes.c_int
    lib.seams_frame_q4_q6.restype = ctypes.c_int
    lib.seams_frame_q4_q6.argtypes = [
        ctypes.c_void_p,
        ctypes.c_double,
        ctypes.c_int,
        ctypes.POINTER(ctypes.c_double),
        ctypes.POINTER(ctypes.c_double),
        ctypes.POINTER(ctypes.c_double),
        ctypes.POINTER(ctypes.c_double),
    ]
    lib.seams_frame_q4_q6_chill.restype = ctypes.c_int
    lib.seams_frame_q4_q6_chill.argtypes = [
        ctypes.c_void_p,
        ctypes.c_double,
        ctypes.c_int,
        ctypes.POINTER(ctypes.c_double),
        ctypes.POINTER(ctypes.c_double),
        ctypes.POINTER(ctypes.c_double),
        ctypes.POINTER(ctypes.c_double),
        ctypes.POINTER(ctypes.c_int),
        ctypes.POINTER(ctypes.c_int),
    ]
    lib.seams_frame_free.argtypes = [ctypes.c_void_p]
    lib.seams_last_error.restype = ctypes.c_char_p
    return lib


def chill_plus(positions, box, rule="chill+", cutoff=3.5, graph="cutoff"):
    """Label one orthorhombic frame.

    ``box`` is ``(ox, oy, oz, lx, ly, lz, xy, xz, yz)`` or three lengths.
    """
    pos = np.ascontiguousarray(positions, dtype=np.float64)
    n = pos.shape[0]
    raw = np.zeros(9, dtype=np.float64)
    lengths = np.asarray(box, dtype=np.float64).reshape(-1)
    if lengths.size == 3:
        raw[3:6] = lengths
    elif lengths.size == 9:
        raw[:] = lengths
    else:
        raise ValueError("box must be 3 lengths or 9 LAMMPS cell numbers")
    types = np.ones(n, dtype=np.int32)
    ids = np.arange(1, n + 1, dtype=np.int32)
    labels = np.empty(n, dtype=np.int32)
    counts = np.zeros(9, dtype=np.int32)
    lib = _library()
    frame = lib.seams_frame_from_arrays(
        pos.ctypes.data_as(ctypes.POINTER(ctypes.c_double)),
        ctypes.c_int(n),
        raw.ctypes.data_as(ctypes.POINTER(ctypes.c_double)),
        types.ctypes.data_as(ctypes.POINTER(ctypes.c_int)),
        ids.ctypes.data_as(ctypes.POINTER(ctypes.c_int)),
    )
    if not frame:
        err = lib.seams_last_error()
        raise RuntimeError(err.decode() if err else "null frame")
    try:
        code = lib.seams_frame_chill_plus_rule(
            ctypes.c_void_p(frame),
            ctypes.c_double(cutoff),
            ctypes.c_int(0),
            ctypes.c_char_p(graph.encode()),
            ctypes.c_char_p(str(rule).encode()),
            labels.ctypes.data_as(ctypes.POINTER(ctypes.c_int)),
            counts.ctypes.data_as(ctypes.POINTER(ctypes.c_int)),
        )
    finally:
        lib.seams_frame_free(ctypes.c_void_p(frame))
    if code != 0:
        err = lib.seams_last_error()
        raise RuntimeError(err.decode() if err else "CHILL+ failed")
    return labels, counts


def q4_q6(positions, box, cutoff=3.5, average=True):
    """Local q4 and q6 on the cutoff shell.

    Degrees 4 and 6 use the solid-harmonic kernel. ``average`` fills the
    Lechner-Dellago means. ``box`` is three lengths or the nine-number cell.
    """
    pos = np.ascontiguousarray(positions, dtype=np.float64)
    n = pos.shape[0]
    raw = np.zeros(9, dtype=np.float64)
    lengths = np.asarray(box, dtype=np.float64).reshape(-1)
    if lengths.size == 3:
        raw[3:6] = lengths
    elif lengths.size == 9:
        raw[:] = lengths
    else:
        raise ValueError("box must be 3 lengths or 9 LAMMPS cell numbers")
    types = np.ones(n, dtype=np.int32)
    ids = np.arange(1, n + 1, dtype=np.int32)
    q4 = np.empty(n, dtype=np.float64)
    q6 = np.empty(n, dtype=np.float64)
    q4bar = np.empty(n, dtype=np.float64) if average else None
    q6bar = np.empty(n, dtype=np.float64) if average else None
    lib = _library()
    frame = lib.seams_frame_from_arrays(
        pos.ctypes.data_as(ctypes.POINTER(ctypes.c_double)),
        ctypes.c_int(n),
        raw.ctypes.data_as(ctypes.POINTER(ctypes.c_double)),
        types.ctypes.data_as(ctypes.POINTER(ctypes.c_int)),
        ids.ctypes.data_as(ctypes.POINTER(ctypes.c_int)),
    )
    if not frame:
        err = lib.seams_last_error()
        raise RuntimeError(err.decode() if err else "null frame")
    try:
        code = lib.seams_frame_q4_q6(
            ctypes.c_void_p(frame),
            ctypes.c_double(cutoff),
            ctypes.c_int(0),
            q4.ctypes.data_as(ctypes.POINTER(ctypes.c_double)),
            q6.ctypes.data_as(ctypes.POINTER(ctypes.c_double)),
            None
            if q4bar is None
            else q4bar.ctypes.data_as(ctypes.POINTER(ctypes.c_double)),
            None
            if q6bar is None
            else q6bar.ctypes.data_as(ctypes.POINTER(ctypes.c_double)),
        )
    finally:
        lib.seams_frame_free(ctypes.c_void_p(frame))
    if code != 0:
        err = lib.seams_last_error()
        raise RuntimeError(err.decode() if err else "q4 q6 failed")
    out = {"q4": q4, "q6": q6}
    if average:
        out["q4bar"] = q4bar
        out["q6bar"] = q6bar
    return out


def q4_q6_chill(positions, box, cutoff=3.5):
    """Local q4, q6, Lechner-Dellago q-bar, and published CHILL+ labels.

    One cutoff list. ``box`` is three lengths or the nine-number cell.
    """
    pos = np.ascontiguousarray(positions, dtype=np.float64)
    n = pos.shape[0]
    raw = np.zeros(9, dtype=np.float64)
    lengths = np.asarray(box, dtype=np.float64).reshape(-1)
    if lengths.size == 3:
        raw[3:6] = lengths
    elif lengths.size == 9:
        raw[:] = lengths
    else:
        raise ValueError("box must be 3 lengths or 9 LAMMPS cell numbers")
    types = np.ones(n, dtype=np.int32)
    ids = np.arange(1, n + 1, dtype=np.int32)
    q4 = np.empty(n, dtype=np.float64)
    q6 = np.empty(n, dtype=np.float64)
    q4bar = np.empty(n, dtype=np.float64)
    q6bar = np.empty(n, dtype=np.float64)
    labels = np.empty(n, dtype=np.int32)
    counts = np.zeros(9, dtype=np.int32)
    lib = _library()
    frame = lib.seams_frame_from_arrays(
        pos.ctypes.data_as(ctypes.POINTER(ctypes.c_double)),
        ctypes.c_int(n),
        raw.ctypes.data_as(ctypes.POINTER(ctypes.c_double)),
        types.ctypes.data_as(ctypes.POINTER(ctypes.c_int)),
        ids.ctypes.data_as(ctypes.POINTER(ctypes.c_int)),
    )
    if not frame:
        err = lib.seams_last_error()
        raise RuntimeError(err.decode() if err else "null frame")
    try:
        code = lib.seams_frame_q4_q6_chill(
            ctypes.c_void_p(frame),
            ctypes.c_double(cutoff),
            ctypes.c_int(0),
            q4.ctypes.data_as(ctypes.POINTER(ctypes.c_double)),
            q6.ctypes.data_as(ctypes.POINTER(ctypes.c_double)),
            q4bar.ctypes.data_as(ctypes.POINTER(ctypes.c_double)),
            q6bar.ctypes.data_as(ctypes.POINTER(ctypes.c_double)),
            labels.ctypes.data_as(ctypes.POINTER(ctypes.c_int)),
            counts.ctypes.data_as(ctypes.POINTER(ctypes.c_int)),
        )
    finally:
        lib.seams_frame_free(ctypes.c_void_p(frame))
    if code != 0:
        err = lib.seams_last_error()
        raise RuntimeError(err.decode() if err else "fused q4 q6 CHILL+ failed")
    return {
        "q4": q4,
        "q6": q6,
        "q4bar": q4bar,
        "q6bar": q6bar,
        "labels": labels,
        "counts": counts,
    }
