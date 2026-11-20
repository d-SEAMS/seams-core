# Copyright (c) d-SEAMS developers
"""Directed hydrogen bonds from a heavy cloud and a hydrogen cloud.

The geometry matches ``seams_directed_hbonds``. The angle is acceptor-first:
0 degrees is a linear donor-hydrogen-acceptor geometry. The acceptor-hydrogen
distance is strictly below ``dist_cutoff``. The angle is kept at or below
``angle_cutoff``. The donor-acceptor pair is strictly inside ``pair_cutoff``.
"""

from __future__ import annotations

import ctypes
import os
from dataclasses import dataclass


class DirectedHbond(ctypes.Structure):
    _fields_ = [
        ("size", ctypes.c_uint32),
        ("donor", ctypes.c_uint32),
        ("hydrogen", ctypes.c_uint32),
        ("acceptor", ctypes.c_uint32),
        ("d_AH", ctypes.c_double),
        ("d_DA", ctypes.c_double),
        ("angle", ctypes.c_double),
    ]


@dataclass(frozen=True)
class DirectedBond:
    donor: int
    hydrogen: int
    acceptor: int
    d_AH: float
    d_DA: float
    angle: float


class SeamsError(RuntimeError):
    def __init__(self, status, message):
        super().__init__(message or "seams status {}".format(status))
        self.status = status


def load_library(path=None):
    """Open libseams. ``path`` defaults to the ``SEAMS_LIB`` environment variable."""
    lib_path = path or os.environ.get("SEAMS_LIB")
    if not lib_path:
        raise SeamsError(-1, "SEAMS_LIB is not set")
    lib = ctypes.CDLL(lib_path)
    lib.seams_system_create.restype = ctypes.c_void_p
    lib.seams_system_destroy.argtypes = [ctypes.c_void_p]
    lib.seams_system_count.argtypes = [ctypes.c_void_p]
    lib.seams_system_count.restype = ctypes.c_int
    lib.seams_last_error.restype = ctypes.c_char_p
    lib.seams_system_set_hydrogen_frame.argtypes = [
        ctypes.c_void_p,
        ctypes.POINTER(ctypes.c_double),
        ctypes.c_size_t,
        ctypes.POINTER(ctypes.c_int),
        ctypes.POINTER(ctypes.c_int),
        ctypes.POINTER(ctypes.c_int),
        ctypes.POINTER(ctypes.c_double),
        ctypes.c_size_t,
        ctypes.POINTER(ctypes.c_int),
        ctypes.POINTER(ctypes.c_int),
        ctypes.POINTER(ctypes.c_int),
        ctypes.POINTER(ctypes.c_double),
        ctypes.c_size_t,
    ]
    lib.seams_system_set_hydrogen_frame.restype = ctypes.c_int
    lib.seams_directed_hbonds.argtypes = [
        ctypes.c_void_p,
        ctypes.c_double,
        ctypes.c_double,
        ctypes.c_double,
        ctypes.POINTER(ctypes.c_int),
        ctypes.c_size_t,
        ctypes.POINTER(ctypes.c_int),
        ctypes.c_size_t,
        ctypes.c_int,
        ctypes.c_int,
        ctypes.POINTER(DirectedHbond),
        ctypes.c_size_t,
        ctypes.POINTER(ctypes.c_size_t),
    ]
    lib.seams_directed_hbonds.restype = ctypes.c_int
    return lib


class Frame:
    """One frame built without a LAMMPS path."""

    def __init__(self, handle, library):
        self._handle = handle
        self._lib = library

    @classmethod
    def from_clouds(
        cls,
        heavy_positions,
        hydrogen_positions,
        box,
        heavy_types,
        heavy_mol_ids,
        hydrogen_types,
        hydrogen_mol_ids,
        heavy_ids=None,
        hydrogen_ids=None,
        library=None,
    ):
        """Store heavies, then hydrogens, with one type and one molecule id each."""
        lib = library or load_library()
        heavy = _doubles(heavy_positions)
        hydrogen = _doubles(hydrogen_positions)
        if len(heavy) % 3 != 0 or len(hydrogen) % 3 != 0:
            raise SeamsError(4, "coordinates are not a multiple of 3")
        n_heavy = len(heavy) // 3
        n_hydrogen = len(hydrogen) // 3
        box9 = _doubles(box)
        if len(box9) != 9:
            raise SeamsError(2, "box length is 9")
        heavy_types_array = _ints(heavy_types)
        heavy_mol_ids_array = _ints(heavy_mol_ids)
        heavy_ids_array = _ints(heavy_ids) if heavy_ids is not None else None
        hydrogen_types_array = _ints(hydrogen_types) if hydrogen_types is not None else None
        hydrogen_mol_ids_array = _ints(hydrogen_mol_ids) if hydrogen_mol_ids is not None else None
        hydrogen_ids_array = _ints(hydrogen_ids) if hydrogen_ids is not None else None
        lengths = (
            ("heavy_types", heavy_types_array, n_heavy),
            ("heavy_mol_ids", heavy_mol_ids_array, n_heavy),
            ("heavy_ids", heavy_ids_array, n_heavy),
            ("hydrogen_types", hydrogen_types_array, n_hydrogen),
            ("hydrogen_mol_ids", hydrogen_mol_ids_array, n_hydrogen),
            ("hydrogen_ids", hydrogen_ids_array, n_hydrogen),
        )
        for name, values, expected in lengths:
            if values is not None and len(values) != expected:
                raise SeamsError(4, "{} length must equal {}".format(name, expected))
            if values is None and expected and name in ("hydrogen_types", "hydrogen_mol_ids"):
                raise SeamsError(4, "{} length must equal {}".format(name, expected))
        handle = lib.seams_system_create()
        if not handle:
            raise SeamsError(1, "frame handle was not created")
        status = lib.seams_system_set_hydrogen_frame(
            handle,
            heavy,
            n_heavy,
            heavy_types_array,
            heavy_ids_array,
            heavy_mol_ids_array,
            hydrogen if n_hydrogen else None,
            n_hydrogen,
            hydrogen_types_array if n_hydrogen else None,
            hydrogen_ids_array,
            hydrogen_mol_ids_array if n_hydrogen else None,
            box9,
            9,
        )
        if status != 0:
            message = lib.seams_last_error()
            lib.seams_system_destroy(handle)
            raise SeamsError(status, message.decode() if message else "")
        return cls(handle, lib)

    def count(self):
        return int(self._lib.seams_system_count(self._handle))

    def directed_hbonds(
        self,
        donor_types,
        acceptor_types,
        hydrogen_type=2,
        pair_cutoff=3.5,
        dist_cutoff=2.42,
        angle_cutoff=30.0,
        exclude_same_mol=True,
    ):
        """Return directed edges ``(donor, hydrogen, acceptor)`` with geometry."""
        donors = _ints(donor_types)
        acceptors = _ints(acceptor_types)
        n_out = ctypes.c_size_t(0)
        status = self._lib.seams_directed_hbonds(
            self._handle,
            pair_cutoff,
            dist_cutoff,
            angle_cutoff,
            donors,
            len(donor_types),
            acceptors,
            len(acceptor_types),
            int(hydrogen_type),
            1 if exclude_same_mol else 0,
            None,
            0,
            ctypes.byref(n_out),
        )
        if status != 0:
            message = self._lib.seams_last_error()
            raise SeamsError(status, message.decode() if message else "")
        count = n_out.value
        if count == 0:
            return []
        rows = (DirectedHbond * count)()
        for row in rows:
            row.size = ctypes.sizeof(DirectedHbond)
        status = self._lib.seams_directed_hbonds(
            self._handle,
            pair_cutoff,
            dist_cutoff,
            angle_cutoff,
            donors,
            len(donor_types),
            acceptors,
            len(acceptor_types),
            int(hydrogen_type),
            1 if exclude_same_mol else 0,
            rows,
            count,
            ctypes.byref(n_out),
        )
        if status != 0:
            message = self._lib.seams_last_error()
            raise SeamsError(status, message.decode() if message else "")
        return [
            DirectedBond(int(row.donor), int(row.hydrogen), int(row.acceptor), row.d_AH, row.d_DA, row.angle)
            for row in rows
        ]

    def close(self):
        if self._handle:
            self._lib.seams_system_destroy(self._handle)
            self._handle = None

    def __del__(self):
        self.close()

    def __enter__(self):
        return self

    def __exit__(self, *_exc):
        self.close()


def _doubles(values):
    return (ctypes.c_double * len(values))(*values)


def _ints(values):
    return (ctypes.c_int * len(values))(*[int(value) for value in values])
