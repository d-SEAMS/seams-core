# Copyright (c) d-SEAMS developers
"""ctypes signatures for include/seams.h. Importing this module does not load the library."""

from __future__ import annotations

import ctypes
import os
from typing import Any

import numpy as np

from seams.errors import SeamsError

_C = ctypes.c_int
_D = ctypes.c_double
_S = ctypes.c_size_t
_P = ctypes.c_void_p
_PD = ctypes.POINTER(_D)
_PI = ctypes.POINTER(_C)
_PS = ctypes.POINTER(_S)
_PP = ctypes.POINTER(_P)


def _text(raw: bytes | None) -> str:
    if not raw:
        return ""
    return raw.decode()


class Api:
    """One loaded ``libseams``.

    Parameters
    ----------
    path : str
        Filesystem path or a loader name. ``SEAMS_LIB`` is the usual path.
    """

    def __init__(self, path: str) -> None:
        self.path = path
        try:
            lib = ctypes.CDLL(path)
        except OSError as exc:
            raise OSError(f"cannot load libseams from {path}") from exc
        self._lib = lib
        self._bind()

    def _fn(self, name: str) -> Any:
        return getattr(self._lib, name)

    def _bind(self) -> None:
        void = None
        specs: dict[str, tuple[Any, list[Any]]] = {
            "seams_abi": (ctypes.c_int, []),
            "seams_version": (ctypes.c_int, []),
            "seams_version_major": (ctypes.c_int, []),
            "seams_version_minor": (ctypes.c_int, []),
            "seams_version_patch": (ctypes.c_int, []),
            "seams_version_string": (ctypes.c_char_p, []),
            "seams_status_string": (ctypes.c_char_p, [_C]),
            "seams_last_error": (ctypes.c_char_p, []),
            "seams_handle_error": (ctypes.c_char_p, [_P]),
            "seams_label_name": (ctypes.c_char_p, [_C]),
            "seams_ring_symmetry_name": (ctypes.c_char_p, [_C]),
            "seams_set_num_threads": (_C, [_C]),
            "seams_get_num_threads": (_C, []),
            "seams_system_create": (_P, []),
            "seams_system_destroy": (void, [_P]),
            "seams_system_set_atoms": (_C, [_P, _PD, _S, _PD, _S, _PI, _PI]),
            "seams_system_set_ortho": (_C, [_P, _PD, _C, _D, _D, _D, _D, _D, _D, _PI, _PI]),
            "seams_system_read": (_C, [_P, ctypes.c_char_p, _C, _C]),
            "seams_path_nframes": (_C, [ctypes.c_char_p, _PS]),
            "seams_system_count": (_C, [_P]),
            "seams_neighbours_build": (_C, [_P, _D, _C, _C, _PP]),
            "seams_neighbours_destroy": (void, [_P]),
            "seams_neighbours_count": (_C, [_P]),
            "seams_neighbours_copy": (_C, [_P, _PI, _S, _PS, _PI, _S, _PS]),
            "seams_chill": (_C, [_P, _D, _C, _C, _C, _PI, _S, _PI]),
            "seams_chill_plus": (_C, [_P, _D, _C, _C, _C, _PI, _PI]),
            "seams_ql": (_C, [_P, _D, _C, _C, _C, _PD, _S, _PD, _S]),
            "seams_q4q6": (_C, [_P, _D, _C, _PD, _PD]),
            "seams_q6bar": (_C, [_P, _D, _C, _D, _PD, _S, _PD, _PD, _PI]),
            "seams_rings_find": (_C, [_P, _D, _C, _C, _C, _PP]),
            "seams_rings_destroy": (void, [_P]),
            "seams_rings_count": (_C, [_P]),
            "seams_rings_histogram": (_C, [_P, _PI, _C]),
            "seams_rings_labels": (_C, [_P, _PI, _PI]),
            "seams_rings_vertices": (_C, [_P, _C, _PI, _S, _PS]),
            "seams_rings_fraction": (_C, [_P, _C, _PD]),
            "seams_cages": (_C, [_P, _D, _C, _C, _C, _C, _PP]),
            "seams_cages_destroy": (void, [_P]),
            "seams_cages_counts": (_C, [_P, _PI, _PI, _PI, _PI]),
            "seams_cages_labels": (_C, [_P, _PI, _S]),
            "seams_signature_cages": (_C, [_P, _D, _C, _C, ctypes.c_char_p, _PS, _PS]),
            "seams_system_symmetry": (_C, [_P, _D, _PI, _S]),
            "seams_brine_analyse": (_C, [_P, _PI, _S, _C, _C, _C, _C, _PP]),
            "seams_brine_destroy": (void, [_P]),
            "seams_brine_fronts": (_C, [_P, _PD, _PD, _PD, _PI]),
            "seams_brine_rejection": (_C, [_P, _PD, _PD, _PD, _PI]),
            "seams_brine_profile": (_C, [_P, _PD, _PD, _PD, _S, _PS]),
            "seams_brine_pockets": (_C, [_P, _PI, _PI]),
            "seams_steinhardt": (_C, [_P, _D, _C, _C, _PD, _S, _PD, _S]),
            "seams_neighbours_knn": (_C, [_P, _C, _D, _C, _C, _PP]),
            "seams_rdf": (_C, [_P, _D, _C, _C, _C, _PD, _PD, _PI, _S, _PS, _PD]),
            "seams_running_cn": (_C, [_P, _D, _C, _C, _C, _PD, _PD, _S, _PS]),
            "seams_coordination": (_C, [_P, _D, _C, _C, _C, _PD]),
            "seams_hbonds": (_C, [_P, _D, _C, _C, _D, _D, _PI]),
            "seams_density_axis": (_C, [_P, _C, _C, _C, _PD, _PD, _S, _PS]),
            "seams_fingerprint": (_C, [_P, _D, _C, _C, _C, _C, _C, _PP]),
            "seams_fingerprint_destroy": (void, [_P]),
            "seams_fingerprint_key": (_C, [_P, ctypes.c_char_p, _S, _PS]),
            "seams_fingerprint_atom": (_C, [_P, _C, ctypes.c_char_p, _S, _PS]),
            "seams_fingerprint_rings": (_C, [_P, _PI, _S]),
            "seams_contact_pairs": (_C, [_P, ctypes.c_char_p, _PS, _PI, _PI]),
            "seams_domains": (_C, [_P, _D, ctypes.c_char_p, _C, _PI, _PI, _PD]),
            "seams_prisms": (_C, [_P, _D, _C, _PI, _PI, _S]),
            "seams_monolayer": (_C, [_P, _D, _C, _C, _D, _PI, _S]),
            "seams_site_table": (_C, [_P, ctypes.c_char_p, _PI, _S]),
        }
        for name, (restype, argtypes) in specs.items():
            fn = self._fn(name)
            fn.restype = restype
            fn.argtypes = argtypes

    def check(self, code: int, handle: int | None) -> None:
        """Raise :class:`seams.errors.SeamsError` when ``code`` is not zero.

        The message is ``seams_handle_error(handle)``. A null handle reads
        the thread-local string.
        """
        if code == 0:
            return
        raw = self._fn("seams_handle_error")(None if handle is None else handle)
        text = _text(raw)
        if not text:
            text = _text(self._fn("seams_status_string")(code)) or "seams call failed"
        raise SeamsError(int(code), text)

    def version(self) -> str:
        return _text(self._fn("seams_version_string")())

    def abi(self) -> int:
        return int(self._fn("seams_abi")())

    def version_parts(self) -> tuple[int, int, int]:
        return (
            int(self._fn("seams_version_major")()),
            int(self._fn("seams_version_minor")()),
            int(self._fn("seams_version_patch")()),
        )

    def set_threads(self, n: int) -> None:
        self.check(int(self._fn("seams_set_num_threads")(n)), None)

    def threads(self) -> int:
        return int(self._fn("seams_get_num_threads")())

    def status_string(self, code: int) -> str:
        return _text(self._fn("seams_status_string")(code))

    def label_name(self, code: int) -> str:
        return _text(self._fn("seams_label_name")(code))

    def ring_name(self, code: int) -> str:
        return _text(self._fn("seams_ring_symmetry_name")(code))

    def create(self) -> int:
        ptr = self._fn("seams_system_create")()
        if not ptr:
            self.check(1, None)
            raise SeamsError(6, "null system")
        return int(ptr)

    def destroy_system(self, ptr: int) -> None:
        self._fn("seams_system_destroy")(ptr)

    def destroy_neighbours(self, ptr: int) -> None:
        self._fn("seams_neighbours_destroy")(ptr)

    def destroy_rings(self, ptr: int) -> None:
        self._fn("seams_rings_destroy")(ptr)

    def destroy_cages(self, ptr: int) -> None:
        self._fn("seams_cages_destroy")(ptr)

    def destroy_brine(self, ptr: int) -> None:
        self._fn("seams_brine_destroy")(ptr)

    def destroy_fingerprint(self, ptr: int) -> None:
        self._fn("seams_fingerprint_destroy")(ptr)

    def count(self, ptr: int) -> int:
        return int(self._fn("seams_system_count")(ptr))

    def nframes(self, path: bytes) -> int:
        out = _S()
        self.check(int(self._fn("seams_path_nframes")(path, ctypes.byref(out))), None)
        return int(out.value)

    def read(self, ptr: int, path: bytes, frame: int, type_filter: int) -> None:
        self.check(int(self._fn("seams_system_read")(ptr, path, frame, type_filter)), ptr)

    def set_atoms(
        self,
        ptr: int,
        xyz: np.ndarray[Any, np.dtype[np.float64]],
        n: int,
        box9: np.ndarray[Any, np.dtype[np.float64]],
        types: np.ndarray[Any, np.dtype[np.int32]] | None,
        ids: np.ndarray[Any, np.dtype[np.int32]] | None,
    ) -> None:
        self.check(
            int(
                self._fn("seams_system_set_atoms")(
                    ptr,
                    _f64(xyz),
                    n,
                    _f64(box9),
                    9,
                    None if types is None else _i32(types),
                    None if ids is None else _i32(ids),
                )
            ),
            ptr,
        )

    def set_ortho(
        self,
        ptr: int,
        xyz: np.ndarray[Any, np.dtype[np.float64]],
        n: int,
        origin: tuple[float, float, float],
        lengths: tuple[float, float, float],
    ) -> None:
        ox, oy, oz = origin
        lx, ly, lz = lengths
        self.check(
            int(
                self._fn("seams_system_set_ortho")(
                    ptr, _f64(xyz), n, ox, oy, oz, lx, ly, lz, None, None
                )
            ),
            ptr,
        )

    def _out_handle(self, code: int, parent: int, slot: ctypes.c_void_p) -> int:
        if code != 0 or not slot.value:
            self.check(code if code != 0 else 6, parent)
            raise SeamsError(6, "null handle")
        return int(slot.value)

    def neighbours_build(self, ptr: int, cutoff: float, type_: int, graph: int) -> int:
        slot = _P()
        code = int(self._fn("seams_neighbours_build")(ptr, cutoff, type_, graph, ctypes.byref(slot)))
        return self._out_handle(code, ptr, slot)

    def knn(self, ptr: int, k: int, cutoff: float, type_: int, mutual: int) -> int:
        slot = _P()
        code = int(self._fn("seams_neighbours_knn")(ptr, k, cutoff, type_, mutual, ctypes.byref(slot)))
        return self._out_handle(code, ptr, slot)

    def neighbour_count(self, ptr: int) -> int:
        return int(self._fn("seams_neighbours_count")(ptr))

    def neighbour_copy(self, ptr: int) -> tuple[np.ndarray[Any, np.dtype[np.int32]], np.ndarray[Any, np.dtype[np.int32]]]:
        n_off = _S()
        nnz = _S()
        self.check(
            int(self._fn("seams_neighbours_copy")(ptr, None, 0, ctypes.byref(n_off), None, 0, ctypes.byref(nnz))),
            ptr,
        )
        offsets = np.empty(int(n_off.value), dtype=np.int32)
        cols = np.empty(int(nnz.value), dtype=np.int32)
        self.check(
            int(
                self._fn("seams_neighbours_copy")(
                    ptr, _i32(offsets), int(n_off.value), ctypes.byref(n_off), _i32(cols), int(nnz.value), ctypes.byref(nnz)
                )
            ),
            ptr,
        )
        return offsets, cols

    def chill(
        self, ptr: int, cutoff: float, type_: int, graph: int, rule: int, n: int
    ) -> tuple[np.ndarray[Any, np.dtype[np.int32]], np.ndarray[Any, np.dtype[np.int32]]]:
        labels = np.empty(n, dtype=np.int32)
        counts = np.empty(9, dtype=np.int32)
        self.check(
            int(self._fn("seams_chill")(ptr, cutoff, type_, graph, rule, _i32(labels), n, _i32(counts))),
            ptr,
        )
        return labels, counts

    def chill_plus(
        self, ptr: int, cutoff: float, type_: int, graph: int, rule: int, n: int
    ) -> tuple[np.ndarray[Any, np.dtype[np.int32]], np.ndarray[Any, np.dtype[np.int32]]]:
        labels = np.empty(n, dtype=np.int32)
        counts = np.empty(9, dtype=np.int32)
        self.check(
            int(self._fn("seams_chill_plus")(ptr, cutoff, type_, graph, rule, _i32(labels), _i32(counts))),
            ptr,
        )
        return labels, counts

    def ql(
        self,
        ptr: int,
        cutoff: float,
        type_: int,
        order: int,
        graph: int,
        n: int,
        average: bool,
    ) -> tuple[np.ndarray[Any, np.dtype[np.float64]], np.ndarray[Any, np.dtype[np.float64]] | None]:
        ql = np.empty(n, dtype=np.float64)
        bar = np.empty(n, dtype=np.float64) if average else None
        self.check(
            int(
                self._fn("seams_ql")(
                    ptr,
                    cutoff,
                    type_,
                    order,
                    graph,
                    _f64(ql),
                    n,
                    None if bar is None else _f64(bar),
                    n,
                )
            ),
            ptr,
        )
        return ql, bar

    def steinhardt(
        self, ptr: int, cutoff: float, type_: int, order: int, n: int, average: bool
    ) -> tuple[np.ndarray[Any, np.dtype[np.float64]], np.ndarray[Any, np.dtype[np.float64]] | None]:
        ql = np.empty(n, dtype=np.float64)
        bar = np.empty(n, dtype=np.float64) if average else None
        self.check(
            int(
                self._fn("seams_steinhardt")(
                    ptr, cutoff, type_, order, _f64(ql), n, None if bar is None else _f64(bar), n
                )
            ),
            ptr,
        )
        return ql, bar

    def q4q6(
        self, ptr: int, cutoff: float, type_: int, n: int
    ) -> tuple[np.ndarray[Any, np.dtype[np.float64]], np.ndarray[Any, np.dtype[np.float64]]]:
        q4 = np.empty(n, dtype=np.float64)
        q6 = np.empty(n, dtype=np.float64)
        self.check(int(self._fn("seams_q4q6")(ptr, cutoff, type_, _f64(q4), _f64(q6))), ptr)
        return q4, q6

    def q6bar(self, ptr: int, cutoff: float, type_: int, ice_cut: float, n: int, store: bool) -> tuple[np.ndarray[Any, np.dtype[np.float64]] | None, float, float, int]:
        values = np.empty(n, dtype=np.float64) if store else None
        mean = _D()
        median = _D()
        n_ice = _C()
        self.check(
            int(
                self._fn("seams_q6bar")(
                    ptr,
                    cutoff,
                    type_,
                    ice_cut,
                    None if values is None else _f64(values),
                    n,
                    ctypes.byref(mean),
                    ctypes.byref(median),
                    ctypes.byref(n_ice),
                )
            ),
            ptr,
        )
        return values, float(mean.value), float(median.value), int(n_ice.value)

    def rings_find(self, ptr: int, cutoff: float, type_: int, graph: int, max_depth: int) -> int:
        slot = _P()
        code = int(self._fn("seams_rings_find")(ptr, cutoff, type_, graph, max_depth, ctypes.byref(slot)))
        return self._out_handle(code, ptr, slot)

    def rings_count(self, ptr: int) -> int:
        return int(self._fn("seams_rings_count")(ptr))

    def rings_histogram(self, ptr: int, length: int) -> np.ndarray[Any, np.dtype[np.int32]]:
        counts = np.zeros(length, dtype=np.int32)
        self.check(int(self._fn("seams_rings_histogram")(ptr, _i32(counts), length)), ptr)
        return counts

    def rings_labels(self, ptr: int, n: int) -> tuple[np.ndarray[Any, np.dtype[np.int32]], np.ndarray[Any, np.dtype[np.int32]]]:
        sizes = np.empty(n, dtype=np.int32)
        symmetry = np.empty(n, dtype=np.int32)
        self.check(int(self._fn("seams_rings_labels")(ptr, _i32(sizes), _i32(symmetry))), ptr)
        return sizes, symmetry

    def rings_vertices(self, ptr: int, index: int) -> np.ndarray[Any, np.dtype[np.int32]]:
        n_out = _S()
        self.check(int(self._fn("seams_rings_vertices")(ptr, index, None, 0, ctypes.byref(n_out))), ptr)
        atoms = np.empty(int(n_out.value), dtype=np.int32)
        if int(n_out.value) == 0:
            return atoms
        self.check(
            int(self._fn("seams_rings_vertices")(ptr, index, _i32(atoms), int(n_out.value), ctypes.byref(n_out))),
            ptr,
        )
        return atoms

    def rings_fraction(self, ptr: int, size: int) -> float:
        fraction = _D()
        self.check(int(self._fn("seams_rings_fraction")(ptr, size, ctypes.byref(fraction))), ptr)
        return float(fraction.value)

    def cages(self, ptr: int, cutoff: float, type_: int, graph: int, k: int, complete: int) -> int:
        slot = _P()
        code = int(self._fn("seams_cages")(ptr, cutoff, type_, graph, k, complete, ctypes.byref(slot)))
        return self._out_handle(code, ptr, slot)

    def cage_counts(self, ptr: int) -> tuple[int, int, int, int]:
        n_hex, n_cubic, n_water, n_both = _C(), _C(), _C(), _C()
        self.check(
            int(
                self._fn("seams_cages_counts")(
                    ptr, ctypes.byref(n_hex), ctypes.byref(n_cubic), ctypes.byref(n_water), ctypes.byref(n_both)
                )
            ),
            ptr,
        )
        return int(n_hex.value), int(n_cubic.value), int(n_water.value), int(n_both.value)

    def cage_labels(self, ptr: int, n: int) -> np.ndarray[Any, np.dtype[np.int32]]:
        labels = np.empty(n, dtype=np.int32)
        self.check(int(self._fn("seams_cages_labels")(ptr, _i32(labels), n)), ptr)
        return labels

    def signature(self, ptr: int, cutoff: float, type_: int, graph: int, spec: bytes) -> tuple[int, int]:
        n_cages, n_atoms = _S(), _S()
        self.check(
            int(self._fn("seams_signature_cages")(ptr, cutoff, type_, graph, spec, ctypes.byref(n_cages), ctypes.byref(n_atoms))),
            ptr,
        )
        return int(n_cages.value), int(n_atoms.value)

    def symmetry(self, ptr: int, shell_rmsd: float, n: int) -> np.ndarray[Any, np.dtype[np.int32]]:
        labels = np.empty(n, dtype=np.int32)
        self.check(int(self._fn("seams_system_symmetry")(ptr, shell_rmsd, _i32(labels), n)), ptr)
        return labels

    def brine(self, ptr: int, ice: np.ndarray[Any, np.dtype[np.int32]], water: int, ion: int, bins: int, axis: int) -> int:
        slot = _P()
        code = int(
            self._fn("seams_brine_analyse")(ptr, _i32(ice), int(ice.shape[0]), water, ion, bins, axis, ctypes.byref(slot))
        )
        return self._out_handle(code, ptr, slot)

    def brine_fronts(self, ptr: int) -> tuple[float, float, float, bool]:
        enter, exit_, thick, found = _D(), _D(), _D(), _C()
        self.check(
            int(
                self._fn("seams_brine_fronts")(
                    ptr, ctypes.byref(enter), ctypes.byref(exit_), ctypes.byref(thick), ctypes.byref(found)
                )
            ),
            ptr,
        )
        return float(enter.value), float(exit_.value), float(thick.value), bool(found.value)

    def brine_rejection(self, ptr: int) -> tuple[float, float, float, bool]:
        c_ice, c_liquid, k, defined = _D(), _D(), _D(), _C()
        self.check(
            int(
                self._fn("seams_brine_rejection")(
                    ptr, ctypes.byref(c_ice), ctypes.byref(c_liquid), ctypes.byref(k), ctypes.byref(defined)
                )
            ),
            ptr,
        )
        return float(c_ice.value), float(c_liquid.value), float(k.value), bool(defined.value)

    def brine_profile(self, ptr: int) -> tuple[np.ndarray[Any, np.dtype[np.float64]], np.ndarray[Any, np.dtype[np.float64]], np.ndarray[Any, np.dtype[np.float64]]]:
        n_out = _S()
        self.check(int(self._fn("seams_brine_profile")(ptr, None, None, None, 0, ctypes.byref(n_out))), ptr)
        n = int(n_out.value)
        z = np.empty(n, dtype=np.float64)
        mol = np.empty(n, dtype=np.float64)
        ice = np.empty(n, dtype=np.float64)
        self.check(int(self._fn("seams_brine_profile")(ptr, _f64(z), _f64(mol), _f64(ice), n, ctypes.byref(n_out))), ptr)
        return z, mol, ice

    def brine_pockets(self, ptr: int) -> tuple[int, int]:
        pockets, waters = _C(), _C()
        self.check(int(self._fn("seams_brine_pockets")(ptr, ctypes.byref(pockets), ctypes.byref(waters))), ptr)
        return int(pockets.value), int(waters.value)

    def rdf(self, ptr: int, rmax: float, bins: int, type_i: int, type_j: int) -> tuple[np.ndarray[Any, np.dtype[np.float64]], np.ndarray[Any, np.dtype[np.float64]], np.ndarray[Any, np.dtype[np.int32]], float]:
        n_out = _S()
        volume = _D()
        self.check(
            int(self._fn("seams_rdf")(ptr, rmax, bins, type_i, type_j, None, None, None, 0, ctypes.byref(n_out), ctypes.byref(volume))),
            ptr,
        )
        n = int(n_out.value)
        radius = np.empty(n, dtype=np.float64)
        g = np.empty(n, dtype=np.float64)
        count = np.empty(n, dtype=np.int32)
        self.check(
            int(
                self._fn("seams_rdf")(
                    ptr, rmax, bins, type_i, type_j, _f64(radius), _f64(g), _i32(count), n, ctypes.byref(n_out), None
                )
            ),
            ptr,
        )
        return radius, g, count, float(volume.value)

    def running_cn(self, ptr: int, rmax: float, bins: int, type_i: int, type_j: int) -> tuple[np.ndarray[Any, np.dtype[np.float64]], np.ndarray[Any, np.dtype[np.float64]]]:
        n_out = _S()
        self.check(
            int(self._fn("seams_running_cn")(ptr, rmax, bins, type_i, type_j, None, None, 0, ctypes.byref(n_out))),
            ptr,
        )
        n = int(n_out.value)
        radius = np.empty(n, dtype=np.float64)
        cn = np.empty(n, dtype=np.float64)
        self.check(
            int(self._fn("seams_running_cn")(ptr, rmax, bins, type_i, type_j, _f64(radius), _f64(cn), n, ctypes.byref(n_out))),
            ptr,
        )
        return radius, cn

    def coordination(self, ptr: int, rmax: float, bins: int, type_i: int, type_j: int) -> float:
        cn = _D()
        self.check(int(self._fn("seams_coordination")(ptr, rmax, bins, type_i, type_j, ctypes.byref(cn))), ptr)
        return float(cn.value)

    def hbonds(self, ptr: int, cutoff: float, oxygen: int, hydrogen: int, dist: float, angle: float) -> int:
        n_bonds = _C()
        self.check(
            int(self._fn("seams_hbonds")(ptr, cutoff, oxygen, hydrogen, dist, angle, ctypes.byref(n_bonds))),
            ptr,
        )
        return int(n_bonds.value)

    def density(self, ptr: int, type_: int, bins: int, axis: int) -> tuple[np.ndarray[Any, np.dtype[np.float64]], np.ndarray[Any, np.dtype[np.float64]]]:
        n_out = _S()
        self.check(int(self._fn("seams_density_axis")(ptr, type_, bins, axis, None, None, 0, ctypes.byref(n_out))), ptr)
        n = int(n_out.value)
        z = np.empty(n, dtype=np.float64)
        rho = np.empty(n, dtype=np.float64)
        self.check(
            int(self._fn("seams_density_axis")(ptr, type_, bins, axis, _f64(z), _f64(rho), n, ctypes.byref(n_out))),
            ptr,
        )
        return z, rho

    def fingerprint(self, ptr: int, cutoff: float, type_: int, k: int, graph: int, hops: int, colour: int) -> int:
        slot = _P()
        code = int(self._fn("seams_fingerprint")(ptr, cutoff, type_, k, graph, hops, colour, ctypes.byref(slot)))
        return self._out_handle(code, ptr, slot)

    def _key(self, ptr: int, fn: str, index: int | None) -> str:
        n_out = _S()
        if index is None:
            self.check(int(self._fn(fn)(ptr, None, 0, ctypes.byref(n_out))), ptr)
        else:
            self.check(int(self._fn(fn)(ptr, index, None, 0, ctypes.byref(n_out))), ptr)
        cap = int(n_out.value) + 1
        buf = ctypes.create_string_buffer(cap)
        if index is None:
            self.check(int(self._fn(fn)(ptr, buf, cap, ctypes.byref(n_out))), ptr)
        else:
            self.check(int(self._fn(fn)(ptr, index, buf, cap, ctypes.byref(n_out))), ptr)
        return buf.value.decode()

    def fingerprint_key(self, ptr: int) -> str:
        return self._key(ptr, "seams_fingerprint_key", None)

    def fingerprint_atom(self, ptr: int, index: int) -> str:
        return self._key(ptr, "seams_fingerprint_atom", index)

    def fingerprint_rings(self, ptr: int, length: int) -> np.ndarray[Any, np.dtype[np.int32]]:
        counts = np.zeros(length, dtype=np.int32)
        self.check(int(self._fn("seams_fingerprint_rings")(ptr, _i32(counts), length)), ptr)
        return counts

    def pairs(self, ptr: int, spec: bytes) -> tuple[int, int, int]:
        n_pairs, n_cation, n_anion = _S(), _C(), _C()
        self.check(
            int(self._fn("seams_contact_pairs")(ptr, spec, ctypes.byref(n_pairs), ctypes.byref(n_cation), ctypes.byref(n_anion))),
            ptr,
        )
        return int(n_pairs.value), int(n_cation.value), int(n_anion.value)

    def domains(self, ptr: int, cutoff: float, spec: bytes, subset: int) -> tuple[int, int, float]:
        n_subset, largest, percolation = _C(), _C(), _D()
        self.check(
            int(
                self._fn("seams_domains")(
                    ptr, cutoff, spec, subset, ctypes.byref(n_subset), ctypes.byref(largest), ctypes.byref(percolation)
                )
            ),
            ptr,
        )
        return int(n_subset.value), int(largest.value), float(percolation.value)

    def prisms(self, ptr: int, cutoff: float, type_: int) -> tuple[np.ndarray[Any, np.dtype[np.int32]], np.ndarray[Any, np.dtype[np.int32]]]:
        perfect = np.zeros(7, dtype=np.int32)
        imperfect = np.zeros(7, dtype=np.int32)
        self.check(int(self._fn("seams_prisms")(ptr, cutoff, type_, _i32(perfect), _i32(imperfect), 7)), ptr)
        return perfect, imperfect

    def monolayer(self, ptr: int, cutoff: float, type_: int, max_depth: int, area: float) -> np.ndarray[Any, np.dtype[np.int32]]:
        counts = np.zeros(max_depth + 1, dtype=np.int32)
        self.check(
            int(self._fn("seams_monolayer")(ptr, cutoff, type_, max_depth, area, _i32(counts), max_depth + 1)),
            ptr,
        )
        return counts

    def site_table(self, ptr: int, spec: bytes, n: int) -> np.ndarray[Any, np.dtype[np.int32]]:
        kinds = np.empty(n, dtype=np.int32)
        self.check(int(self._fn("seams_site_table")(ptr, spec, _i32(kinds), n)), ptr)
        return kinds


def _f64(arr: np.ndarray[Any, np.dtype[np.float64]]) -> Any:
    return arr.ctypes.data_as(ctypes.POINTER(ctypes.c_double))


def _i32(arr: np.ndarray[Any, np.dtype[np.int32]]) -> Any:
    return arr.ctypes.data_as(ctypes.POINTER(ctypes.c_int))


_API: Api | None = None


def api() -> Api:
    """Return the process-wide :class:`Api`, loading it on first use.

    The path is ``SEAMS_LIB`` or the loader name ``libseams``.
    """
    global _API
    if _API is None:
        _API = Api(os.environ.get("SEAMS_LIB", "libseams"))
    return _API
