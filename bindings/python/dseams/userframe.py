# Copyright (c) d-SEAMS developers
"""Frame handle over libseams_lang.

The classifiers, rings, cages, Steinhardt parameters, prisms, monolayer
rings, and the in-plane RDF are the C API. This module only owns the
frame pointer and the buffers those calls fill.
"""

from __future__ import annotations

import ctypes
import ctypes.util
import os
import sys
from pathlib import Path

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


class IceCounts(dict):
    """CHILL+ histogram. Missing names read as 0."""

    def __missing__(self, key):
        return 0

    def __getattr__(self, name):
        if name.startswith("_"):
            raise AttributeError(name)
        return self[name]

    def __repr__(self):
        parts = [f"{key}={value}" for key, value in sorted(self.items()) if value]
        return "IceCounts(" + ", ".join(parts) + ")"


class CageScore:
    """Per-atom hexagonal-cage and double-diamond flags.

    ``n_ih`` counts atoms in a hexagonal cage, ``n_ic`` atoms in a
    double-diamond cage, and ``n_water`` atoms in neither. An atom in
    both contributes to ``n_ih`` and ``n_ic``.
    """

    def __init__(self, hc, ddc):
        self.hc = [bool(flag) for flag in hc]
        self.ddc = [bool(flag) for flag in ddc]

    @property
    def n_ih(self):
        return sum(1 for flag in self.hc if flag)

    @property
    def n_ic(self):
        return sum(1 for flag in self.ddc if flag)

    @property
    def n_water(self):
        return sum(1 for hc, ddc in zip(self.hc, self.ddc) if not hc and not ddc)

    @property
    def n_both(self):
        return sum(1 for hc, ddc in zip(self.hc, self.ddc) if hc and ddc)

    def __repr__(self):
        return f"CageScore(n_ih={self.n_ih}, n_ic={self.n_ic}, n_water={self.n_water})"


class _PrismRow(ctypes.Structure):
    _fields_ = (
        ("size", ctypes.c_int),
        ("perfect", ctypes.c_int),
        ("imperfect", ctypes.c_int),
    )


def _lang_library():
    """``SEAMS_LANG_LIB``, else the copy installed with this package, else the loader's."""
    path = os.environ.get("SEAMS_LANG_LIB")
    if path:
        return path
    name = {"darwin": "libseams_lang.dylib", "win32": "seams_lang.dll"}.get(sys.platform, "libseams_lang.so")
    here = Path(__file__).resolve().parent
    for folder in (here, here.parent / ".pydseamslib.mesonpy.libs", Path(sys.prefix) / "lib"):
        if (folder / name).is_file():
            return str(folder / name)
    return ctypes.util.find_library("seams_lang") or name


_LIB = None


def _library():
    global _LIB
    if _LIB is not None:
        return _LIB
    lib = ctypes.CDLL(_lang_library())
    lib.seams_last_error.restype = ctypes.c_char_p
    lib.seams_frame_read.restype = ctypes.c_void_p
    lib.seams_frame_read.argtypes = [ctypes.c_char_p, ctypes.c_int, ctypes.c_int]
    lib.seams_frame_from_arrays.restype = ctypes.c_void_p
    lib.seams_frame_from_arrays.argtypes = [
        ctypes.POINTER(ctypes.c_double),
        ctypes.c_int,
        ctypes.POINTER(ctypes.c_double),
        ctypes.POINTER(ctypes.c_int),
        ctypes.POINTER(ctypes.c_int),
    ]
    lib.seams_frame_free.argtypes = [ctypes.c_void_p]
    lib.seams_frame_nop.restype = ctypes.c_int
    lib.seams_frame_nop.argtypes = [ctypes.c_void_p]
    lib.seams_frame_positions.restype = ctypes.POINTER(ctypes.c_double)
    lib.seams_frame_positions.argtypes = [ctypes.c_void_p]
    lib.seams_frame_box.argtypes = [ctypes.c_void_p, ctypes.POINTER(ctypes.c_double)]
    lib.seams_frame_types.restype = ctypes.POINTER(ctypes.c_int)
    lib.seams_frame_types.argtypes = [ctypes.c_void_p]
    lib.seams_frame_chill_plus.restype = ctypes.c_int
    lib.seams_frame_chill_plus.argtypes = [
        ctypes.c_void_p,
        ctypes.c_double,
        ctypes.c_int,
        ctypes.c_char_p,
        ctypes.POINTER(ctypes.c_int),
        ctypes.POINTER(ctypes.c_int),
    ]
    lib.seams_frame_cages.restype = ctypes.c_int
    lib.seams_frame_cages.argtypes = [
        ctypes.c_void_p,
        ctypes.c_double,
        ctypes.c_int,
        ctypes.c_int,
        ctypes.c_char_p,
        ctypes.c_int,
        ctypes.POINTER(ctypes.c_int),
        ctypes.POINTER(ctypes.c_int),
        ctypes.POINTER(ctypes.c_int),
        ctypes.POINTER(ctypes.c_int),
        ctypes.POINTER(ctypes.c_int),
    ]
    lib.seams_frame_signature.restype = ctypes.c_int
    lib.seams_frame_signature.argtypes = [
        ctypes.c_void_p,
        ctypes.c_double,
        ctypes.c_int,
        ctypes.c_char_p,
        ctypes.c_char_p,
        ctypes.POINTER(ctypes.c_int),
        ctypes.POINTER(ctypes.c_int),
    ]
    lib.seams_frame_steinhardt.restype = ctypes.c_int
    lib.seams_frame_steinhardt.argtypes = [
        ctypes.c_void_p,
        ctypes.c_double,
        ctypes.c_int,
        ctypes.c_int,
        ctypes.POINTER(ctypes.c_double),
        ctypes.POINTER(ctypes.c_double),
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
    lib.seams_frame_rdf2d_nbin.restype = ctypes.c_int
    lib.seams_frame_rdf2d_nbin.argtypes = [ctypes.c_double, ctypes.c_double]
    lib.seams_frame_rdf2d.restype = ctypes.c_int
    lib.seams_frame_rdf2d.argtypes = [
        ctypes.c_void_p,
        ctypes.c_double,
        ctypes.c_double,
        ctypes.c_char_p,
        ctypes.POINTER(ctypes.c_double),
        ctypes.POINTER(ctypes.c_double),
    ]
    lib.seams_frame_prisms.restype = ctypes.c_int
    lib.seams_frame_prisms.argtypes = [
        ctypes.c_void_p,
        ctypes.c_double,
        ctypes.c_int,
        ctypes.c_int,
        ctypes.c_char_p,
        ctypes.POINTER(_PrismRow),
        ctypes.POINTER(ctypes.c_int),
    ]
    lib.seams_frame_monolayer.restype = ctypes.c_int
    lib.seams_frame_monolayer.argtypes = [
        ctypes.c_void_p,
        ctypes.c_double,
        ctypes.c_int,
        ctypes.c_int,
        ctypes.c_double,
        ctypes.c_char_p,
        ctypes.POINTER(ctypes.c_int),
        ctypes.c_int,
    ]
    _LIB = lib
    return lib


def _check(lib, code):
    if code != 0:
        err = lib.seams_last_error()
        raise RuntimeError(err.decode() if err else "d-SEAMS call failed")


class Frame:
    """One configuration. Classification calls libseams_lang."""

    def __init__(self, ptr, symbols=None, cutoff=3.5):
        if not ptr:
            raise RuntimeError("null frame")
        self._ptr = ptr
        self._symbols = list(symbols) if symbols is not None else None
        self.cutoff = float(cutoff)
        self._labels = None
        self._cages = None
        self._owned = True

    def close(self):
        if self._owned and self._ptr:
            _library().seams_frame_free(self._ptr)
            self._ptr = None
            self._owned = False

    def __del__(self):
        self.close()

    @classmethod
    def _adopt(cls, ptr, symbols=None, cutoff=3.5):
        if not ptr:
            err = _library().seams_last_error()
            raise RuntimeError(err.decode() if err else "could not read the frame")
        return cls(ptr, symbols=symbols, cutoff=cutoff)

    @property
    def n_atoms(self):
        return int(_library().seams_frame_nop(self._ptr))

    @property
    def positions(self):
        n = self.n_atoms
        raw = _library().seams_frame_positions(self._ptr)
        return np.ctypeslib.as_array(raw, shape=(n, 3)).copy()

    @property
    def box(self):
        raw = np.zeros(9, dtype=np.float64)
        _library().seams_frame_box(self._ptr, raw.ctypes.data_as(ctypes.POINTER(ctypes.c_double)))
        if np.any(raw[6:] != 0.0):
            return raw.copy()
        return raw[3:6].copy()

    def chill_plus(self, graph="cutoff"):
        """Published CHILL+ histogram. Index names match the C label order."""
        n = self.n_atoms
        labels = np.zeros(n, dtype=np.int32)
        counts = np.zeros(9, dtype=np.int32)
        code = _library().seams_frame_chill_plus(
            self._ptr,
            ctypes.c_double(self.cutoff),
            ctypes.c_int(0),
            graph.encode(),
            labels.ctypes.data_as(ctypes.POINTER(ctypes.c_int)),
            counts.ctypes.data_as(ctypes.POINTER(ctypes.c_int)),
        )
        _check(_library(), code)
        self._labels = [ _NAMES[int(label)] if 0 <= int(label) < len(_NAMES) else "unclassified" for label in labels ]
        named = IceCounts({name: int(counts[i]) for i, name in enumerate(_NAMES)})
        return named

    def cages(self, seeded=True, k=4, ring_adjacent=False):
        """Hexagonal and double-diamond membership from seams_frame_cages."""
        n = self.n_atoms
        per = np.zeros(n, dtype=np.int32)
        n_hex = ctypes.c_int()
        n_cubic = ctypes.c_int()
        n_water = ctypes.c_int()
        n_both = ctypes.c_int()
        graph = b"seeded" if seeded else b"cutoff"
        code = _library().seams_frame_cages(
            self._ptr,
            ctypes.c_double(self.cutoff),
            ctypes.c_int(0),
            ctypes.c_int(int(k)),
            graph,
            ctypes.c_int(1 if ring_adjacent else 0),
            per.ctypes.data_as(ctypes.POINTER(ctypes.c_int)),
            ctypes.byref(n_hex),
            ctypes.byref(n_cubic),
            ctypes.byref(n_water),
            ctypes.byref(n_both),
        )
        _check(_library(), code)
        hc = [int(code_) in (1, 3) for code_ in per]
        ddc = [int(code_) in (2, 3) for code_ in per]
        self._cages = CageScore(hc, ddc)
        return self._cages

    def cages_by_signature(self, signature, max_ring_size=None):
        """Closed cages of one ring census.

        The C API returns the cage count and the number of distinct
        vertices, not one vertex list per cage. The list length is the
        cage count so ``len`` matches the documented print. ``n_atoms``
        is the vertex count.
        """
        del max_ring_size  # the C search picks the depth from the signature
        n_cages = ctypes.c_int()
        n_atoms = ctypes.c_int()
        code = _library().seams_frame_signature(
            self._ptr,
            ctypes.c_double(self.cutoff),
            ctypes.c_int(0),
            b"cutoff",
            str(signature).encode(),
            ctypes.byref(n_cages),
            ctypes.byref(n_atoms),
        )
        _check(_library(), code)
        found = [
            {"signature": str(signature), "vertices": [], "n_atoms": int(n_atoms.value)}
            for _ in range(int(n_cages.value))
        ]
        found_n = int(n_atoms.value)

        class _Found(list):
            n_atoms = found_n

        out = _Found(found)
        out.n_atoms = found_n
        return out

    def steinhardt(self, order_l=6):
        """Local q_l and the Lechner-Dellago average from seams_frame_steinhardt."""
        n = self.n_atoms
        ql = np.zeros(n, dtype=np.float64)
        qlbar = np.zeros(n, dtype=np.float64)
        code = _library().seams_frame_steinhardt(
            self._ptr,
            ctypes.c_double(self.cutoff),
            ctypes.c_int(0),
            ctypes.c_int(int(order_l)),
            ql.ctypes.data_as(ctypes.POINTER(ctypes.c_double)),
            qlbar.ctypes.data_as(ctypes.POINTER(ctypes.c_double)),
        )
        _check(_library(), code)
        return {"ql": ql.tolist(), "ql_bar": qlbar.tolist()}

    def ring_census(self, max_size=8):
        """Primitive rings by size. ``counts[i]`` from seams_frame_rings."""
        depth = int(max_size)
        counts = np.zeros(depth + 1, dtype=np.int32)
        code = _library().seams_frame_rings(
            self._ptr,
            ctypes.c_double(self.cutoff),
            ctypes.c_int(0),
            b"cutoff",
            ctypes.c_int(depth),
            counts.ctypes.data_as(ctypes.POINTER(ctypes.c_int)),
            ctypes.c_int(counts.size),
        )
        _check(_library(), code)
        return {size: int(counts[size]) for size in range(3, depth + 1) if int(counts[size])}

    def rdf_2d(self, output_dir, cutoff=12.0, binwidth=0.05):
        """In-plane g(r) from seams_frame_rdf2d."""
        nbin = int(_library().seams_frame_rdf2d_nbin(ctypes.c_double(cutoff), ctypes.c_double(binwidth)))
        radius = np.zeros(nbin, dtype=np.float64)
        values = np.zeros(nbin, dtype=np.float64)
        code = _library().seams_frame_rdf2d(
            self._ptr,
            ctypes.c_double(cutoff),
            ctypes.c_double(binwidth),
            str(output_dir).encode(),
            radius.ctypes.data_as(ctypes.POINTER(ctypes.c_double)),
            values.ctypes.data_as(ctypes.POINTER(ctypes.c_double)),
        )
        _check(_library(), code)
        return radius.tolist(), values.tolist()

    def monolayer_rings(self, output_dir, sheet_area, max_depth=4):
        """Quasi-2D ring counts from seams_frame_monolayer.

        Coverage is the number the C writer stores in
        ``topoMonolayer/coverageAreaXY.dat``.
        """
        depth = int(max_depth)
        counts = np.zeros(depth + 1, dtype=np.int32)
        code = _library().seams_frame_monolayer(
            self._ptr,
            ctypes.c_double(self.cutoff),
            ctypes.c_int(0),
            ctypes.c_int(depth),
            ctypes.c_double(float(sheet_area)),
            str(output_dir).encode(),
            counts.ctypes.data_as(ctypes.POINTER(ctypes.c_int)),
            ctypes.c_int(counts.size),
        )
        _check(_library(), code)
        report = Path(output_dir) / "topoMonolayer" / "coverageAreaXY.dat"
        if not report.is_file():
            return {
                size: {"count": int(counts[size]), "coverage_xy": 0.0}
                for size in range(3, depth + 1)
            }
        fields = report.read_text(encoding="utf-8").splitlines()[-1].split()[1:]
        found = {}
        for size, number, area in zip(fields[::3], fields[1::3], fields[2::3]):
            found[int(size)] = {"count": int(number), "coverage_xy": float(area)}
        return found

    def find_prisms(self, output_dir="output/", max_depth=6, shape_matching=False):
        """Prism census from seams_frame_prisms. The return is the C rows."""
        del shape_matching
        depth = int(max_depth)
        n_rows = ctypes.c_int(depth - 2)
        rows = (_PrismRow * (depth - 2))()
        code = _library().seams_frame_prisms(
            self._ptr,
            ctypes.c_double(self.cutoff),
            ctypes.c_int(0),
            ctypes.c_int(depth),
            str(output_dir).encode(),
            rows,
            ctypes.byref(n_rows),
        )
        _check(_library(), code)
        return [
            {"size": int(row.size), "perfect": int(row.perfect), "imperfect": int(row.imperfect)}
            for row in rows
        ]

    def to_ase(self):
        """ASE Atoms for this frame, with ice labels and cage flags when set."""
        import ase
        from ase import Atoms

        del ase
        n = self.n_atoms
        symbols = self._symbols if self._symbols is not None else ["O"] * n
        cell = np.asarray(self.box, dtype=float).reshape(-1)
        if cell.size == 9:
            lx, ly, lz, xy, xz, yz = cell[3:9]
            matrix = [[lx, 0.0, 0.0], [xy, ly, 0.0], [xz, yz, lz]]
        else:
            matrix = cell
        atoms = Atoms(symbols=symbols, positions=self.positions, cell=matrix, pbc=True)
        if self._labels is not None:
            atoms.arrays["ice_type"] = np.array(self._labels, dtype=object)
        if self._cages is not None:
            atoms.arrays["hc"] = np.array(self._cages.hc)
            atoms.arrays["ddc"] = np.array(self._cages.ddc)
        atoms.info["dseams_n_atoms"] = n
        return atoms

    def to_solvis(self, expand_box=True):
        """solvis System around :meth:`to_ase`. solvis-tools must be installed."""
        try:
            from solvis.system import System
        except ImportError as exc:
            raise ImportError("solvis interop needs solvis-tools") from exc
        return System(self.to_ase(), expand_box=expand_box)


def read(filename, frame=1, atom_type=None, all_atoms=False, **kwargs):
    """Read one frame with seams_frame_read.

    ``atom_type`` keeps that LAMMPS species. ``all_atoms`` keeps every
    atom. Otherwise the library guesses oxygen, which is what the
    tutorials do when they omit the type.
    """
    del kwargs
    if all_atoms:
        filt = -1
    elif atom_type is None:
        filt = 0
    else:
        filt = int(atom_type)
    ptr = _library().seams_frame_read(os.fsencode(filename), int(frame), filt)
    return Frame._adopt(ptr)


def _select_mask(atoms, select):
    if select is None:
        return np.ones(len(atoms), dtype=bool)
    symbols = np.array(atoms.get_chemical_symbols())
    if isinstance(select, str):
        return symbols == select
    if isinstance(select, (list, tuple, set, frozenset)):
        wanted = {item if isinstance(item, str) else int(item) for item in select}
        return np.array([sym in wanted or int(z) in wanted for sym, z in zip(symbols, atoms.numbers)])
    return np.array([int(z) == int(select) for z in atoms.numbers])


def from_ase(atoms, select="O", cutoff=None, bonded="auto"):
    """Frame for one ASE species. The coordinates go through seams_frame_from_arrays."""
    del bonded
    if not hasattr(atoms, "get_positions"):
        raise TypeError("from_ase expects an ASE Atoms object")
    keep = _select_mask(atoms, select)
    if not np.any(keep):
        raise ValueError(f"no atoms matched select={select!r}")
    import dseams

    box, rot = dseams._reduced_cell(atoms.cell.array)
    positions = np.ascontiguousarray(
        np.asarray(atoms.get_positions(), dtype=np.float64) @ rot.T, dtype=np.float64
    )[keep]
    symbols = [sym for sym, yes in zip(atoms.get_chemical_symbols(), keep) if yes]
    n = positions.shape[0]
    types = np.ones(n, dtype=np.int32)
    ids = np.arange(1, n + 1, dtype=np.int32)
    ptr = _library().seams_frame_from_arrays(
        positions.ctypes.data_as(ctypes.POINTER(ctypes.c_double)),
        ctypes.c_int(n),
        box.ctypes.data_as(ctypes.POINTER(ctypes.c_double)),
        types.ctypes.data_as(ctypes.POINTER(ctypes.c_int)),
        ids.ctypes.data_as(ctypes.POINTER(ctypes.c_int)),
    )
    frame = Frame._adopt(ptr, symbols=symbols, cutoff=3.5 if cutoff is None else cutoff)
    return frame


def to_solvis(frame, expand_box=True):
    return frame.to_solvis(expand_box=expand_box)
