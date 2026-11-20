#!/usr/bin/env python3
import collections
import ctypes
import pathlib
import sys


class LoadOptions(ctypes.Structure):
    _fields_ = [
        ("resnames", ctypes.POINTER(ctypes.c_char_p)),
        ("resname_count", ctypes.c_int),
        ("excluded_resnames", ctypes.POINTER(ctypes.c_char_p)),
        ("excluded_resname_count", ctypes.c_int),
        ("atomnames", ctypes.POINTER(ctypes.c_char_p)),
        ("atomname_count", ctypes.c_int),
        ("excluded_atomnames", ctypes.POINTER(ctypes.c_char_p)),
        ("excluded_atomname_count", ctypes.c_int),
        ("type_map", ctypes.POINTER(ctypes.c_char_p)),
        ("type_map_count", ctypes.c_int),
    ]


def strings(values):
    return (ctypes.c_char_p * len(values))(*values)


def main() -> int:
    if len(sys.argv) != 3:
        raise SystemExit("usage: test_chemfiles_options.py LIB FIXTURE.gro")
    lib = ctypes.CDLL(sys.argv[1])
    lib.seams_frame_read_options.restype = ctypes.c_void_p
    lib.seams_frame_read_options.argtypes = [
        ctypes.c_char_p,
        ctypes.c_int,
        ctypes.c_int,
        ctypes.POINTER(LoadOptions),
    ]
    lib.seams_frame_nop.restype = ctypes.c_int
    lib.seams_frame_nop.argtypes = [ctypes.c_void_p]
    lib.seams_frame_types.restype = ctypes.POINTER(ctypes.c_int)
    lib.seams_frame_types.argtypes = [ctypes.c_void_p]
    lib.seams_frame_free.argtypes = [ctypes.c_void_p]
    lib.seams_last_error.restype = ctypes.c_char_p

    excluded_resnames = strings([b"SOLX"])
    excluded_atomnames = strings([b"MW"])
    type_map = strings(
        [b"O_2=21", b"N_2=22", b"C_3=23", b"C_2=24", b"C_R=25", b"H_=26"]
    )
    options = LoadOptions(
        None,
        0,
        excluded_resnames,
        len(excluded_resnames),
        None,
        0,
        excluded_atomnames,
        len(excluded_atomnames),
        type_map,
        len(type_map),
    )
    path = str(pathlib.Path(sys.argv[2])).encode()
    frame = lib.seams_frame_read_options(path, 1, -1, ctypes.byref(options))
    if not frame:
        sys.stderr.write(lib.seams_last_error().decode() + "\n")
        return 1
    try:
        count = lib.seams_frame_nop(frame)
        if count != 1320:
            sys.stderr.write(f"expected 1320 atoms, got {count}\n")
            return 1
        types = lib.seams_frame_types(frame)
        census = collections.Counter(types[index] for index in range(count))
        expected = {
            1: 80,
            8: 40,
            21: 60,
            22: 60,
            23: 120,
            24: 120,
            25: 360,
            26: 480,
        }
        if dict(census) != expected:
            sys.stderr.write(f"unexpected type census: {dict(census)}\n")
            return 1
    finally:
        lib.seams_frame_free(frame)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
