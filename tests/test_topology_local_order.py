"""Python-facing checks for explicit rings, winding, tetrahedral q, and LSI."""

from __future__ import annotations

import ctypes
import math
import os


lib = ctypes.CDLL(os.environ["SEAMS_LANG_LIB"])
c_int_p = ctypes.POINTER(ctypes.c_int)
c_double_p = ctypes.POINTER(ctypes.c_double)

lib.seams_frame_from_arrays.argtypes = [
    c_double_p,
    ctypes.c_int,
    c_double_p,
    c_int_p,
    c_int_p,
]
lib.seams_frame_from_arrays.restype = ctypes.c_void_p
lib.seams_frame_free.argtypes = [ctypes.c_void_p]
lib.seams_graph_rings.argtypes = [
    c_int_p,
    ctypes.c_int,
    c_int_p,
    ctypes.c_int,
    ctypes.c_int,
    c_int_p,
    ctypes.c_int,
    c_int_p,
    ctypes.c_int,
    c_int_p,
    ctypes.c_int,
    c_int_p,
    c_int_p,
]
lib.seams_graph_rings.restype = ctypes.c_int
lib.seams_frame_ring_winding.argtypes = [
    ctypes.c_void_p,
    c_int_p,
    ctypes.c_int,
    c_int_p,
    ctypes.c_int,
    ctypes.c_int,
    c_int_p,
    c_int_p,
]
lib.seams_frame_ring_winding.restype = ctypes.c_int
lib.seams_frame_tetrahedral.argtypes = [
    ctypes.c_void_p,
    c_int_p,
    ctypes.c_int,
    c_int_p,
    ctypes.c_int,
    ctypes.c_int,
    ctypes.c_double,
    c_double_p,
]
lib.seams_frame_tetrahedral.restype = ctypes.c_int
lib.seams_frame_lsi.argtypes = [
    ctypes.c_void_p,
    c_int_p,
    ctypes.c_int,
    c_int_p,
    ctypes.c_int,
    ctypes.c_double,
    c_double_p,
]
lib.seams_frame_lsi.restype = ctypes.c_int


def ints(values: list[int]) -> ctypes.Array[ctypes.c_int]:
    return (ctypes.c_int * len(values))(*values)


def doubles(values: list[float]) -> ctypes.Array[ctypes.c_double]:
    return (ctypes.c_double * len(values))(*values)


def rings_from_csr(
    offsets: list[int], columns: list[int], max_depth: int
) -> tuple[list[list[int]], list[int]]:
    off = ints(offsets)
    col = ints(columns)
    histogram = (ctypes.c_int * (max_depth + 1))()
    n_rings = ctypes.c_int()
    n_vertices = ctypes.c_int()
    rc = lib.seams_graph_rings(
        off,
        len(offsets),
        col,
        len(columns),
        max_depth,
        histogram,
        len(histogram),
        None,
        0,
        None,
        0,
        ctypes.byref(n_rings),
        ctypes.byref(n_vertices),
    )
    assert rc == 0
    ring_offsets = (ctypes.c_int * (n_rings.value + 1))()
    vertices = (ctypes.c_int * n_vertices.value)()
    rc = lib.seams_graph_rings(
        off,
        len(offsets),
        col,
        len(columns),
        max_depth,
        histogram,
        len(histogram),
        ring_offsets,
        len(ring_offsets),
        vertices,
        len(vertices),
        ctypes.byref(n_rings),
        ctypes.byref(n_vertices),
    )
    assert rc == 0
    rings = [
        list(vertices[ring_offsets[i] : ring_offsets[i + 1]])
        for i in range(n_rings.value)
    ]
    return rings, list(histogram)


def make_frame(xyz: list[float], box: list[float]) -> ctypes.c_void_p:
    frame = lib.seams_frame_from_arrays(
        doubles(xyz), len(xyz) // 3, doubles(box), None, None
    )
    assert frame
    return frame


def test_explicit_graph_ring_histogram() -> None:
    rings, histogram = rings_from_csr(
        [0, 1, 2, 3, 4, 5, 6], [1, 2, 3, 4, 5, 0], 6
    )
    assert rings == [[0, 1, 2, 3, 4, 5]]
    assert histogram[6] == 1


def test_periodic_axis_winding() -> None:
    rings, _ = rings_from_csr([0, 1, 2, 3], [1, 2, 0], 3)
    frame = make_frame(
        [1, 4, 4, 4, 4, 4, 7, 4, 4],
        [0, 0, 0, 9, 9, 9, 0, 0, 0],
    )
    try:
        packed = ints(rings[0])
        offsets = ints([0, len(rings[0])])
        winding = (ctypes.c_int * 1)()
        contractible = (ctypes.c_int * 1)()
        rc = lib.seams_frame_ring_winding(
            frame,
            offsets,
            1,
            packed,
            len(packed),
            0,
            winding,
            contractible,
        )
        assert rc == 0
        assert winding[0] != 0
        assert contractible[0] == 0
    finally:
        lib.seams_frame_free(frame)


def test_tetrahedral_order_and_lsi_oracles() -> None:
    box = [0, 0, 0, 100, 100, 100, 0, 0, 0]
    tetra = make_frame(
        [
            50,
            50,
            50,
            51,
            51,
            51,
            51,
            49,
            49,
            49,
            51,
            49,
            49,
            49,
            51,
        ],
        box,
    )
    centre = ints([0])
    partners = ints([1, 2, 3, 4])
    q = (ctypes.c_double * 1)()
    try:
        assert (
            lib.seams_frame_tetrahedral(
                tetra, centre, 1, partners, 4, 1, 2.0, q
            )
            == 0
        )
        assert math.isclose(q[0], 1.0, abs_tol=1e-12)
        assert (
            lib.seams_frame_tetrahedral(
                tetra, centre, 1, partners, 3, 0, 0.0, q
            )
            == 0
        )
        assert math.isnan(q[0])
    finally:
        lib.seams_frame_free(tetra)

    radial = make_frame(
        [50, 50, 50, 51, 50, 50, 50, 51.5, 50, 50, 50, 52.5, 47, 50, 50],
        box,
    )
    lsi = (ctypes.c_double * 1)()
    try:
        assert lib.seams_frame_lsi(radial, centre, 1, partners, 4, 3.1, lsi) == 0
        assert math.isclose(lsi[0], 1.0 / 18.0, abs_tol=1e-12)
    finally:
        lib.seams_frame_free(radial)
