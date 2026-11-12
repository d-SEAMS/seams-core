# Copyright (c) d-SEAMS developers
"""In-situ CHILL+ and primitive rings through libseams_lang."""

from __future__ import annotations

import ctypes
import os
import time

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


class _PhaseTimes(ctypes.Structure):
    _fields_ = [
        ("neighbours", ctypes.c_double),
        ("chill", ctypes.c_double),
        ("rings", ctypes.c_double),
        ("reused", ctypes.c_int),
    ]


class InSitu:
    """One library handle. Each observe() is one frame.

    A skin above zero keeps the neighbour list while every oxygen has
    moved by at most half that skin. The bonds are still the pairs
    inside the cutoff.
    """

    def __init__(self, library=None, cutoff=3.5, graph="cutoff", max_ring=6, skin=0.0):
        path = library or os.environ.get("SEAMS_LANG_LIB", "libseams_lang")
        self.lib = ctypes.CDLL(path)
        self.cutoff = float(cutoff)
        self.graph = graph
        self.max_ring = int(max_ring)
        self.skin = float(skin)
        self.frames = []
        self.keep_frames = True
        self._frame = None
        self._frame_n = 0
        self._oxy = None
        self._box = np.zeros(9, dtype=np.float64)
        self._labels = None
        self._counts = np.zeros(9, dtype=np.int32)
        self._rings = np.zeros(self.max_ring + 1, dtype=np.int32)
        self._phase = _PhaseTimes()
        self.profile = {
            "copy": 0.0,
            "neighbours": 0.0,
            "chill": 0.0,
            "rings": 0.0,
            "reports": 0,
            "reused": 0,
        }
        lib = self.lib
        lib.seams_last_error.restype = ctypes.c_char_p
        lib.seams_frame_from_arrays.restype = ctypes.c_void_p
        lib.seams_frame_from_arrays.argtypes = [
            ctypes.POINTER(ctypes.c_double),
            ctypes.c_int,
            ctypes.POINTER(ctypes.c_double),
            ctypes.POINTER(ctypes.c_int),
            ctypes.POINTER(ctypes.c_int),
        ]
        lib.seams_frame_free.argtypes = [ctypes.c_void_p]
        lib.seams_frame_chill_plus.restype = ctypes.c_int
        lib.seams_frame_chill_plus.argtypes = [
            ctypes.c_void_p,
            ctypes.c_double,
            ctypes.c_int,
            ctypes.c_char_p,
            ctypes.POINTER(ctypes.c_int),
            ctypes.POINTER(ctypes.c_int),
        ]
        lib.seams_frame_rings.restype = ctypes.c_int
        lib.seams_frame_rings.argtypes = [
            ctypes.c_void_p,
            ctypes.c_double,
            ctypes.c_int,
            ctypes.c_char_p,
            ctypes.c_int,
            ctypes.POINTER(ctypes.c_int),
            ctypes.c_int,
        ]
        lib.seams_frame_read.restype = ctypes.c_void_p
        lib.seams_frame_read.argtypes = [ctypes.c_char_p, ctypes.c_int, ctypes.c_int]
        lib.seams_frame_nop.restype = ctypes.c_int
        lib.seams_frame_nop.argtypes = [ctypes.c_void_p]
        lib.seams_frame_set_positions.restype = ctypes.c_int
        lib.seams_frame_set_positions.argtypes = [
            ctypes.c_void_p,
            ctypes.POINTER(ctypes.c_double),
            ctypes.c_int,
            ctypes.POINTER(ctypes.c_double),
        ]
        lib.seams_frame_chill_rings.restype = ctypes.c_int
        lib.seams_frame_chill_rings.argtypes = [
            ctypes.c_void_p,
            ctypes.c_double,
            ctypes.c_int,
            ctypes.c_char_p,
            ctypes.c_int,
            ctypes.c_double,
            ctypes.POINTER(ctypes.c_int),
            ctypes.POINTER(ctypes.c_int),
            ctypes.POINTER(ctypes.c_int),
            ctypes.c_int,
            ctypes.POINTER(_PhaseTimes),
        ]

    def _fail(self, code, what):
        if code == 0:
            return
        err = self.lib.seams_last_error()
        text = err.decode() if err else what
        raise RuntimeError(text)

    def box_from_rows(self, rows):
        """OpenMM box vectors, in angstrom, as nine LAMMPS cell numbers."""
        vecs = np.asarray(rows, dtype=np.float64).reshape(3, 3)
        a, b, c = vecs
        raw = self._box
        raw[:] = 0.0
        raw[3] = a[0]
        raw[4] = b[1]
        raw[5] = c[2]
        raw[6] = b[0]
        raw[7] = c[0]
        raw[8] = c[1]
        return raw

    def reset_profile(self):
        for key in self.profile:
            self.profile[key] = 0

    def _box9(self, box):
        """Three edge lengths or the nine LAMMPS numbers.

        A length-3 box is orthorhombic: slots 3-5 are lx, ly, lz and the
        origin and tilts stay zero. Any other length is rejected.
        """
        lengths = np.asarray(box, dtype=np.float64).reshape(-1).copy()
        raw = self._box
        raw[:] = 0.0
        if lengths.size == 3:
            raw[3:6] = lengths
        elif lengths.size == 9:
            raw[:] = lengths
        else:
            raise ValueError("box must be 3 lengths or 9 LAMMPS cell numbers")
        return raw

    def _hold(self, positions, box):
        """One frame object. Later reports write into its coordinate buffer."""
        pos = np.ascontiguousarray(positions, dtype=np.float64)
        n = int(pos.shape[0])
        raw = self._box9(box)
        if self._frame is not None and n == self._frame_n:
            self._fail(
                self.lib.seams_frame_set_positions(
                    ctypes.c_void_p(self._frame),
                    pos.ctypes.data_as(ctypes.POINTER(ctypes.c_double)),
                    ctypes.c_int(n),
                    raw.ctypes.data_as(ctypes.POINTER(ctypes.c_double)),
                ),
                "positions",
            )
            return
        if self._frame is not None:
            self.lib.seams_frame_free(ctypes.c_void_p(self._frame))
            self._frame = None
        types = np.ones(n, dtype=np.int32)
        ids = np.arange(1, n + 1, dtype=np.int32)
        frame = self.lib.seams_frame_from_arrays(
            pos.ctypes.data_as(ctypes.POINTER(ctypes.c_double)),
            ctypes.c_int(n),
            raw.ctypes.data_as(ctypes.POINTER(ctypes.c_double)),
            types.ctypes.data_as(ctypes.POINTER(ctypes.c_int)),
            ids.ctypes.data_as(ctypes.POINTER(ctypes.c_int)),
        )
        if not frame:
            self._fail(1, "null frame")
        self._frame = frame
        self._frame_n = n
        self._labels = np.empty(n, dtype=np.int32)

    def _classify(self):
        self._fail(
            self.lib.seams_frame_chill_rings(
                ctypes.c_void_p(self._frame),
                ctypes.c_double(self.cutoff),
                ctypes.c_int(1),
                ctypes.c_char_p(self.graph.encode()),
                ctypes.c_int(self.max_ring),
                ctypes.c_double(self.skin),
                self._labels.ctypes.data_as(ctypes.POINTER(ctypes.c_int)),
                self._counts.ctypes.data_as(ctypes.POINTER(ctypes.c_int)),
                self._rings.ctypes.data_as(ctypes.POINTER(ctypes.c_int)),
                ctypes.c_int(self.max_ring + 1),
                ctypes.byref(self._phase),
            ),
            "CHILL+ and rings",
        )
        prof = self.profile
        prof["neighbours"] += self._phase.neighbours
        prof["chill"] += self._phase.chill
        prof["rings"] += self._phase.rings
        prof["reports"] += 1
        prof["reused"] += int(self._phase.reused)

    def analyse(self, positions, box):
        """CHILL+ labels and primitive-ring counts for one orthorhombic frame."""
        self._hold(positions, box)
        self._classify()
        return self._labels.copy(), self._counts.copy(), self._rings.copy()

    def analyse_file(self, path, frame_index):
        """Same classifiers, reading a dump from disk."""
        lib = self.lib
        held = lib.seams_frame_read(str(path).encode(), ctypes.c_int(frame_index), ctypes.c_int(1))
        if not held:
            self._fail(1, "unreadable frame")
        n = lib.seams_frame_nop(ctypes.c_void_p(held))
        labels = np.empty(n, dtype=np.int32)
        counts = np.zeros(9, dtype=np.int32)
        rings = np.zeros(self.max_ring + 1, dtype=np.int32)
        try:
            self._fail(
                lib.seams_frame_chill_plus(
                    ctypes.c_void_p(held),
                    ctypes.c_double(self.cutoff),
                    ctypes.c_int(1),
                    ctypes.c_char_p(self.graph.encode()),
                    labels.ctypes.data_as(ctypes.POINTER(ctypes.c_int)),
                    counts.ctypes.data_as(ctypes.POINTER(ctypes.c_int)),
                ),
                "CHILL+",
            )
            self._fail(
                lib.seams_frame_rings(
                    ctypes.c_void_p(held),
                    ctypes.c_double(self.cutoff),
                    ctypes.c_int(1),
                    ctypes.c_char_p(self.graph.encode()),
                    ctypes.c_int(self.max_ring),
                    rings.ctypes.data_as(ctypes.POINTER(ctypes.c_int)),
                    ctypes.c_int(self.max_ring + 1),
                ),
                "rings",
            )
        finally:
            lib.seams_frame_free(ctypes.c_void_p(held))
        return labels.copy(), counts.copy(), rings.copy()

    def observe(self, step, positions, box, ids, types=None, copy_s=0.0):
        """Classify one OpenMM frame. copy_s is the host-copy time already spent."""
        del ids
        t0 = time.perf_counter()
        pos = np.asarray(positions, dtype=np.float64)
        if pos.ndim != 2:
            pos = pos.reshape(-1, 3)
        if types is None:
            kept = np.ascontiguousarray(pos, dtype=np.float64)
        else:
            kind = np.asarray(types)
            if getattr(self, "_index_n", -1) != kind.shape[0]:
                self._index = np.flatnonzero(kind == 1)
                self._index_n = int(kind.shape[0])
            n = int(self._index.shape[0])
            if self._oxy is None or self._oxy.shape != (n, 3):
                self._oxy = np.empty((n, 3), dtype=np.float64)
            np.take(pos, self._index, axis=0, out=self._oxy)
            kept = self._oxy
        self.profile["copy"] += float(copy_s) + (time.perf_counter() - t0)
        self._hold(kept, box)
        self._classify()
        if self.keep_frames:
            self.frames.append(
                {
                    "step": int(step),
                    "positions": np.array(kept, copy=True),
                    "box": np.array(self._box, copy=True),
                    "labels": self._labels.copy(),
                    "counts": self._counts.copy(),
                    "rings": self._rings.copy(),
                }
            )
        return self._labels, self._counts, self._rings
