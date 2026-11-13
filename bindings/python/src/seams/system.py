# Copyright (c) d-SEAMS developers
"""Frame handles. Each handle is a context manager and is safe to drop."""

from __future__ import annotations

import os
from collections.abc import Callable
from typing import Any

import numpy as np
import numpy.typing as npt

from seams import _lib
from seams.errors import ERR_ARG, ERR_STATE, SeamsError
from seams.results import (
    BrineFronts,
    BrinePockets,
    BrineProfile,
    BrineRejection,
    CageCounts,
    ChillResult,
    Contacts,
    Density,
    Domains,
    NeighbourCSR,
    OrderResult,
    Prisms,
    Q4Q6,
    Q6Bar,
    Radial,
    RingLabels,
    RunningCN,
    Signature,
)

ArrayLike = npt.ArrayLike

CUTOFF = 0
MUTUAL = 1
SEEDED = 2
PUBLISHED = 0
OVITO = 1

LABEL_CUBIC = 0
LABEL_HEXAGONAL = 1
LABEL_WATER = 2
LABEL_INTERFACIAL = 3
LABEL_CLATHRATE = 4
LABEL_INTERCLATHRATE = 5
LABEL_UNCLASSIFIED = 6
LABEL_RECUBIC = 7
LABEL_REHEX = 8

CAGE_WATER = 0
CAGE_HEX = 1
CAGE_CUBIC = 2
CAGE_BOTH = 3


def as_box9(box: ArrayLike) -> npt.NDArray[np.float64]:
    """Nine LAMMPS cell numbers.

    Parameters
    ----------
    box : array_like
        Three edge lengths, or ``ox, oy, oz, lx, ly, lz, xy, xz, yz``.
        Lengths are angstrom. Three lengths place the origin at zero
        and set the tilt to zero.

    Returns
    -------
    numpy.ndarray
        Shape ``(9,)``, dtype ``float64``.

    Raises
    ------
    ValueError
        If ``box`` does not have 3 or 9 numbers.

    Examples
    --------
    >>> as_box9((4.0, 5.0, 6.0)).tolist()
    [0.0, 0.0, 0.0, 4.0, 5.0, 6.0, 0.0, 0.0, 0.0]
    """
    raw = np.asarray(box, dtype=np.float64).reshape(-1)
    out = np.zeros(9, dtype=np.float64)
    if raw.size == 3:
        out[3:6] = raw
    elif raw.size == 9:
        out[:] = raw
    else:
        raise ValueError("box is 3 lengths or 9 LAMMPS cell numbers")
    return out


def _positions(positions: ArrayLike) -> tuple[npt.NDArray[np.float64], int]:
    pos = np.ascontiguousarray(positions, dtype=np.float64)
    if pos.ndim != 2 or pos.shape[1] != 3:
        raise ValueError("positions are (n, 3) angstrom coordinates")
    if pos.shape[0] <= 0:
        raise ValueError("positions are empty")
    flat = np.ascontiguousarray(pos.reshape(-1))
    return flat, int(pos.shape[0])


def _optional_i32(values: ArrayLike | None, n: int, name: str) -> npt.NDArray[np.int32] | None:
    if values is None:
        return None
    arr = np.ascontiguousarray(values, dtype=np.int32).reshape(-1)
    if arr.shape[0] != n:
        raise ValueError(f"{name} has length {arr.shape[0]}, expected {n}")
    return arr


class _Owned:
    """A C handle freed by ``close``, by the context manager, or by ``__del__``."""

    def __init__(self, ptr: int, destroy: Callable[[int], None]) -> None:
        self._ptr: int | None = ptr
        self._destroy = destroy

    def close(self) -> None:
        """Free the handle. A second call does nothing."""
        ptr = self._ptr
        if ptr is None:
            return
        self._ptr = None
        self._destroy(ptr)

    def __enter__(self) -> _Owned:
        return self

    def __exit__(self, *_exc: object) -> None:
        self.close()

    def __del__(self) -> None:
        try:
            self.close()
        except Exception:
            return

    def _live(self) -> int:
        if self._ptr is None:
            raise SeamsError(ERR_STATE, "handle is closed")
        return self._ptr

    @property
    def closed(self) -> bool:
        """True after :meth:`close`."""
        return self._ptr is None


class Neighbours(_Owned):
    """A neighbour list copied out of a frame.

    The list does not borrow the :class:`System`. Closing the system
    leaves this object valid. :meth:`copy` returns NumPy arrays the
    caller owns. The deprecated ``seams_neighbours_csr`` borrow is not
    exposed.

    Examples
    --------
    Use it as a context manager so the handle is freed on the way out.
    """

    def __enter__(self) -> Neighbours:
        return self

    def count(self) -> int:
        """Number of neighbour rows.

        Returns
        -------
        int
        """
        return _lib.api().neighbour_count(self._live())

    def copy(self) -> NeighbourCSR:
        """Copy the CSR into new arrays.

        Returns
        -------
        NeighbourCSR
            Offsets of length ``n + 1`` and columns that omit self.
        """
        offsets, columns = _lib.api().neighbour_copy(self._live())
        return NeighbourCSR(offsets, columns)


class Rings(_Owned):
    """Primitive rings (Franzblau) for one frame.

    The caller owns the handle. Copied arrays do not borrow it.
    """

    def __enter__(self) -> Rings:
        return self

    def count(self) -> int:
        """Number of primitive rings.

        Returns
        -------
        int
        """
        return _lib.api().rings_count(self._live())

    def histogram(self, length: int = 13) -> npt.NDArray[np.int32]:
        """Count rings by size.

        Parameters
        ----------
        length : int
            ``counts[k]`` is the number of rings of size ``k``.

        Returns
        -------
        numpy.ndarray
        """
        if length <= 0:
            raise ValueError("length is positive")
        return _lib.api().rings_histogram(self._live(), length)

    def labels(self) -> RingLabels:
        """Size and symmetry of every ring.

        Returns
        -------
        RingLabels
        """
        n = self.count()
        if n == 0:
            empty = np.empty(0, dtype=np.int32)
            return RingLabels(empty, empty)
        sizes, symmetry = _lib.api().rings_labels(self._live(), n)
        return RingLabels(sizes, symmetry)

    def vertices(self, index: int = 0) -> npt.NDArray[np.int32]:
        """Atom indices of one ring.

        Parameters
        ----------
        index : int
            Ring index, from 0.

        Returns
        -------
        numpy.ndarray
            Positions in the frame.
        """
        return _lib.api().rings_vertices(self._live(), index)

    def fraction(self, size: int = 0) -> float:
        """Fraction of atoms that sit on a ring.

        Parameters
        ----------
        size : int
            Ring size. Zero counts atoms that sit on any ring.

        Returns
        -------
        float
            A value in ``[0, 1]``.
        """
        return _lib.api().rings_fraction(self._live(), size)


class CageSet(_Owned):
    """Hexagonal and double-diamond membership.

    The four population counts are disjoint.
    """

    def __enter__(self) -> CageSet:
        return self

    def counts(self) -> CageCounts:
        """Disjoint population counts.

        Returns
        -------
        CageCounts
        """
        n_hex, n_cubic, n_water, n_both = _lib.api().cage_counts(self._live())
        return CageCounts(n_hex, n_cubic, n_water, n_both)

    def labels(self, n: int) -> npt.NDArray[np.int32]:
        """One cage code per atom.

        Parameters
        ----------
        n : int
            Frame count stored when the set was built. Atoms the species
            filter dropped are ``CAGE_WATER``.

        Returns
        -------
        numpy.ndarray
            ``seams_cage_label`` codes.
        """
        return _lib.api().cage_labels(self._live(), n)


class Fingerprint(_Owned):
    """Topological fingerprint of one bonded frame."""

    def __enter__(self) -> Fingerprint:
        return self

    def key(self) -> str:
        """Frame key text.

        Returns
        -------
        str
        """
        return _lib.api().fingerprint_key(self._live())

    def atom(self, index: int) -> str:
        """Local key of one atom.

        Parameters
        ----------
        index : int
            Atom index, from 0.

        Returns
        -------
        str
        """
        return _lib.api().fingerprint_atom(self._live(), index)

    def rings(self, length: int = 13) -> npt.NDArray[np.int32]:
        """Primitive-ring census stored on the fingerprint.

        Parameters
        ----------
        length : int
            ``counts[k]`` is the number of rings of size ``k``.

        Returns
        -------
        numpy.ndarray
        """
        return _lib.api().fingerprint_rings(self._live(), length)


class Brine(_Owned):
    """One planar ice/brine front.

    Lengths are angstrom in the frame cell. Molality is formula units
    per kilogram of water.
    """

    def __enter__(self) -> Brine:
        return self

    def fronts(self) -> BrineFronts:
        """The two planar crossings of the ice slab.

        Returns
        -------
        BrineFronts
        """
        enter, exit_, thick, found = _lib.api().brine_fronts(self._live())
        return BrineFronts(enter, exit_, thick, found)

    def rejection(self) -> BrineRejection:
        """Molal rejection coefficient.

        Returns
        -------
        BrineRejection
        """
        c_ice, c_liquid, k, defined = _lib.api().brine_rejection(self._live())
        return BrineRejection(c_ice, c_liquid, k, defined)

    def profile(self) -> BrineProfile:
        """Concentration profile along the front axis.

        Returns
        -------
        BrineProfile
        """
        z, mol, ice = _lib.api().brine_profile(self._live())
        return BrineProfile(z, mol, ice)

    def pockets(self) -> BrinePockets:
        """Brine pockets, excluding the liquid reservoir.

        Returns
        -------
        BrinePockets
        """
        n_pockets, waters = _lib.api().brine_pockets(self._live())
        return BrinePockets(n_pockets, waters)


class System(_Owned):
    """One frame.

    The constructor calls ``seams_system_create``. :meth:`close` calls
    ``seams_system_destroy``. Neighbour, ring, cage, fingerprint, and
    brine handles keep their own copies, so they stay valid after this
    frame is closed.

    A failing call raises :class:`seams.errors.SeamsError`. ``status``
    is the ``seams_status`` code and ``message`` is the string from
    ``seams_handle_error`` on this handle.

    Parameters
    ----------
    path : str, optional
        When set, the frame is read before the constructor returns.

    frame : int
        1-based frame index. Used only with ``path``.

    type_filter : int
        Negative keeps every atom. Zero guesses oxygen. A positive
        value keeps that species. Used only with ``path``.

    Examples
    --------
    >>> err = SeamsError(5, "no frame")
    >>> err.status == 5
    True
    """

    def __init__(self, path: str | os.PathLike[str] | None = None, *, frame: int = 1, type_filter: int = -1) -> None:
        ptr = _lib.api().create()
        super().__init__(ptr, _lib.api().destroy_system)
        if path is not None:
            try:
                self.read(path, frame=frame, type_filter=type_filter)
            except Exception:
                self.close()
                raise

    def __enter__(self) -> System:
        return self

    def __repr__(self) -> str:
        if self.closed:
            return "System(closed)"
        return f"System(n={self.count()})"

    def count(self) -> int:
        """Number of points in the frame.

        Returns
        -------
        int
            Zero when the handle holds no frame.
        """
        return _lib.api().count(self._live())

    def read(self, path: str | os.PathLike[str], *, frame: int = 1, type_filter: int = -1) -> int:
        """Read one frame from a trajectory.

        Parameters
        ----------
        path : path-like
            CON, LAMMPS dump, or another format the linked readers accept.
        frame : int
            1-based index, in file order.
        type_filter : int
            Negative keeps every atom. Zero guesses oxygen. A positive
            value keeps that species.

        Returns
        -------
        int
            Atom count after the read.

        Raises
        ------
        seams.SeamsError
            On a rejected path, index, or reader failure. The previous
            frame is kept.
        """
        _lib.api().read(self._live(), os.fsencode(path), frame, type_filter)
        return self.count()

    def set_frame(
        self,
        positions: ArrayLike,
        box: ArrayLike,
        *,
        types: ArrayLike | None = None,
        ids: ArrayLike | None = None,
    ) -> None:
        """Store points and a cell.

        Parameters
        ----------
        positions : array_like
            Shape ``(n, 3)``, angstrom. The library copies them.
        box : array_like
            Three lengths or nine LAMMPS numbers.
        types, ids : array_like, optional
            Length ``n``. Omitted types are 1. Omitted ids are ``1..n``.
            Ids must be unique.

        Raises
        ------
        ValueError
            If the arrays have the wrong shape.
        seams.SeamsError
            If the cell has no volume or an id is repeated.
        """
        flat, n = _positions(positions)
        box9 = np.ascontiguousarray(as_box9(box))
        _lib.api().set_atoms(
            self._live(),
            flat,
            n,
            box9,
            _optional_i32(types, n, "types"),
            _optional_i32(ids, n, "ids"),
        )

    def set_ortho(
        self,
        positions: ArrayLike,
        lengths: tuple[float, float, float],
        *,
        origin: tuple[float, float, float] = (0.0, 0.0, 0.0),
    ) -> None:
        """Store points in an orthorhombic cell.

        Parameters
        ----------
        positions : array_like
            Shape ``(n, 3)``, angstrom.
        lengths : tuple of float
            Edge lengths ``lx, ly, lz``, angstrom.
        origin : tuple of float
            Origin ``ox, oy, oz``, angstrom. The tilt is zero.
        """
        flat, n = _positions(positions)
        if n > np.iinfo(np.int32).max:
            raise ValueError("n does not fit in a C int; use set_frame")
        _lib.api().set_ortho(self._live(), flat, n, origin, lengths)

    def neighbours(self, cutoff: float = 3.5, *, type: int = 0, graph: int = CUTOFF) -> Neighbours:
        """Build a cutoff or mutual neighbour list.

        Parameters
        ----------
        cutoff : float
            Positive cutoff, angstrom. The mutual graph keeps the four
            nearest neighbours inside ``cutoff + 1.5``.
        type : int
            Species to keep. Zero and negative keep every atom.
        graph : int
            ``CUTOFF`` or ``MUTUAL``. ``SEEDED`` is rejected.

        Returns
        -------
        Neighbours
        """
        if graph == SEEDED:
            raise SeamsError(ERR_ARG, "seeded graph is not a neighbour list")
        ptr = _lib.api().neighbours_build(self._live(), float(cutoff), type, graph)
        return Neighbours(ptr, _lib.api().destroy_neighbours)

    def knn(self, k: int = 4, cutoff: float = 3.5, *, type: int = 0, mutual: bool = True) -> Neighbours:
        """Build a k-nearest neighbour list.

        Parameters
        ----------
        k : int
            Positive neighbour count.
        cutoff : float
            Nominal cutoff, angstrom. Candidates lie inside ``cutoff + 1.5``.
        type : int
            Species filter. Zero keeps every atom.
        mutual : bool
            True keeps the mutual graph.

        Returns
        -------
        Neighbours
        """
        ptr = _lib.api().knn(self._live(), k, float(cutoff), type, 1 if mutual else 0)
        return Neighbours(ptr, _lib.api().destroy_neighbours)

    def chill(
        self,
        cutoff: float = 3.5,
        *,
        type: int = 0,
        graph: int = CUTOFF,
        rule: int = PUBLISHED,
    ) -> ChillResult:
        """Original CHILL labels.

        Parameters
        ----------
        cutoff : float
            Positive cutoff, angstrom.
        type : int
            Species filter. Zero and negative keep every atom.
        graph : int
            ``CUTOFF`` or ``MUTUAL``.
        rule : int
            ``PUBLISHED`` (inclusive edges) or ``OVITO`` (open edges).

        Returns
        -------
        ChillResult
            Labels and nine population counts.
        """
        n = self.count()
        labels, counts = _lib.api().chill(self._live(), float(cutoff), type, graph, rule, n)
        return ChillResult(labels, counts)

    def chill_plus(
        self,
        cutoff: float = 3.5,
        *,
        type: int = 0,
        graph: int = CUTOFF,
        rule: int = PUBLISHED,
    ) -> ChillResult:
        """CHILL+ labels.

        Parameters
        ----------
        cutoff : float
            Positive cutoff, angstrom.
        type : int
            Species filter.
        graph : int
            ``CUTOFF`` or ``MUTUAL``.
        rule : int
            ``PUBLISHED`` is the Nguyen and Molinero table. ``OVITO`` is
            the Chill+ modifier.

        Returns
        -------
        ChillResult
        """
        n = self.count()
        labels, counts = _lib.api().chill_plus(self._live(), float(cutoff), type, graph, rule, n)
        return ChillResult(labels, counts)

    def ql(
        self,
        order: int,
        cutoff: float = 3.5,
        *,
        type: int = 0,
        graph: int = CUTOFF,
        average: bool = False,
    ) -> OrderResult:
        """Local ``q_l``, and the Lechner-Dellago average on request.

        Parameters
        ----------
        order : int
            Spherical-harmonic degree. 3, 4, 6, 8, or 12.
        cutoff : float
            Positive cutoff, angstrom.
        type : int
            Species filter.
        graph : int
            ``CUTOFF`` or ``MUTUAL``.
        average : bool
            When true, also return the neighbour mean of ``q_lm``.

        Returns
        -------
        OrderResult

        Notes
        -----
        This call shares a process-global scratch buffer with
        :meth:`q4q6`, :meth:`q6bar`, and :meth:`steinhardt`.
        """
        n = self.count()
        ql, bar = _lib.api().ql(self._live(), float(cutoff), type, order, graph, n, average)
        return OrderResult(ql, bar)

    def steinhardt(
        self,
        order: int,
        cutoff: float = 3.5,
        *,
        type: int = 0,
        average: bool = False,
    ) -> OrderResult:
        """Local ``q_l`` on the cutoff shell.

        Parameters
        ----------
        order : int
            3, 4, 6, 8, or 12.
        cutoff : float
            Positive cutoff, angstrom.
        type : int
            Species filter.
        average : bool
            When true, include the Lechner-Dellago average.

        Returns
        -------
        OrderResult
        """
        n = self.count()
        ql, bar = _lib.api().steinhardt(self._live(), float(cutoff), type, order, n, average)
        return OrderResult(ql, bar)

    def q4q6(self, cutoff: float = 3.5, *, type: int = 0) -> Q4Q6:
        """Local q4 and q6 on the cutoff shell.

        Parameters
        ----------
        cutoff : float
            Positive cutoff, angstrom.
        type : int
            Species filter.

        Returns
        -------
        Q4Q6
        """
        n = self.count()
        q4, q6 = _lib.api().q4q6(self._live(), float(cutoff), type, n)
        return Q4Q6(q4, q6)

    def q6bar(self, cutoff: float = 3.5, *, type: int = 0, ice_cut: float = 0.36, values: bool = True) -> Q6Bar:
        """Lechner-Dellago q6-bar and the ice count above a cut.

        Parameters
        ----------
        cutoff : float
            Positive cutoff, angstrom.
        type : int
            Species filter.
        ice_cut : float
            Finite threshold. The nucleation census uses 0.36.
        values : bool
            When true, keep the per-atom array.

        Returns
        -------
        Q6Bar
        """
        n = self.count()
        stored, mean, median, n_ice = _lib.api().q6bar(self._live(), float(cutoff), type, float(ice_cut), n, values)
        return Q6Bar(stored, mean, median, n_ice)

    def rings(self, cutoff: float = 3.5, *, type: int = 0, graph: int = CUTOFF, max_depth: int = 6) -> Rings:
        """Find primitive rings up to a maximum size.

        Parameters
        ----------
        cutoff : float
            Positive cutoff, angstrom.
        type : int
            Species filter.
        graph : int
            ``CUTOFF`` or ``MUTUAL``.
        max_depth : int
            Largest ring size, from 3 through 12.

        Returns
        -------
        Rings
        """
        ptr = _lib.api().rings_find(self._live(), float(cutoff), type, graph, max_depth)
        return Rings(ptr, _lib.api().destroy_rings)

    def cages(
        self,
        cutoff: float = 3.5,
        *,
        type: int = 0,
        graph: int = SEEDED,
        k: int = 4,
        complete: bool = False,
    ) -> CageSet:
        """Hexagonal and double-diamond cage membership.

        Parameters
        ----------
        cutoff : float
            Positive cutoff, angstrom. The seeded graph uses
            ``cutoff + 1.5`` as the candidate shell.
        type : int
            Species filter.
        graph : int
            ``CUTOFF``, ``MUTUAL``, or ``SEEDED``.
        k : int
            Neighbour count for the seeded graph. Zero selects 4.
        complete : bool
            Run ring-adjacent completion on the seeded labels.

        Returns
        -------
        CageSet
        """
        ptr = _lib.api().cages(self._live(), float(cutoff), type, graph, k, 1 if complete else 0)
        return CageSet(ptr, _lib.api().destroy_cages)

    def signature(self, spec: str = "512", cutoff: float = 3.5, *, type: int = 0, graph: int = CUTOFF) -> Signature:
        """Count closed cages whose ring census matches ``spec``.

        Parameters
        ----------
        spec : str
            A census such as ``"512"``, ``"51262"``, or the names ``"hc"``
            and ``"ddc"``.
        cutoff : float
            Positive cutoff, angstrom.
        type : int
            Species filter.
        graph : int
            ``CUTOFF`` or ``MUTUAL``.

        Returns
        -------
        Signature
        """
        n_cages, n_atoms = _lib.api().signature(self._live(), float(cutoff), type, graph, spec.encode())
        return Signature(n_cages, n_atoms)

    def symmetry(self, shell_rmsd: float = 0.0) -> npt.NDArray[np.int32]:
        """Per-atom ice shell label.

        Parameters
        ----------
        shell_rmsd : float
            Positive scale-free residual ceiling. A non-positive value
            selects the library ceiling.

        Returns
        -------
        numpy.ndarray
            0 ice Ih, 1 ice Ic, 2 stacking disorder, or -1.
        """
        return _lib.api().symmetry(self._live(), float(shell_rmsd), self.count())

    def brine(
        self,
        ice: ArrayLike,
        *,
        water_type: int = 0,
        ion_type: int = 0,
        n_bins: int = 0,
        axis: int = 2,
    ) -> Brine:
        """Planar ice/brine front from an ice assignment.

        Parameters
        ----------
        ice : array_like
            Per-atom ice flags, length equal to the frame count.
            Non-zero means ice.
        water_type : int
            Species kept as water, or 0 to keep every atom that is not
            ``ion_type``.
        ion_type : int
            Species kept as ions, or 0 when the frame has no ions.
        n_bins : int
            Bin count along ``axis``. Zero selects 20.
        axis : int
            0, 1, or 2.

        Returns
        -------
        Brine
        """
        flags = np.ascontiguousarray(ice, dtype=np.int32).reshape(-1)
        if flags.shape[0] != self.count():
            raise ValueError("ice flags must match the frame count")
        ptr = _lib.api().brine(self._live(), flags, water_type, ion_type, n_bins, axis)
        return Brine(ptr, _lib.api().destroy_brine)

    def rdf(self, rmax: float = 6.0, bins: int = 0, *, type_i: int = 0, type_j: int = 0) -> Radial:
        """Partial radial distribution.

        Parameters
        ----------
        rmax : float
            Maximum radius, angstrom.
        bins : int
            Bin count. Zero selects about ``rmax / 0.1``.
        type_i, type_j : int
            Species pair. Zero keeps every atom.

        Returns
        -------
        Radial
        """
        r, g, count, volume = _lib.api().rdf(self._live(), float(rmax), bins, type_i, type_j)
        return Radial(r, g, count, volume)

    def running_cn(self, rmax: float = 6.0, bins: int = 0, *, type_i: int = 0, type_j: int = 0) -> RunningCN:
        """Running coordination number along the partial RDF.

        Parameters
        ----------
        rmax : float
            Maximum radius, angstrom.
        bins : int
            Bin count. Zero selects the library default.
        type_i, type_j : int
            Species pair.

        Returns
        -------
        RunningCN
        """
        r, cn = _lib.api().running_cn(self._live(), float(rmax), bins, type_i, type_j)
        return RunningCN(r, cn)

    def coordination(self, rmax: float = 3.5, bins: int = 8, *, type_i: int = 0, type_j: int = 0) -> float:
        """Coordination number integrated to ``rmax``.

        Parameters
        ----------
        rmax : float
            Cutoff, angstrom.
        bins : int
            Bin count used by the integral.
        type_i, type_j : int
            Species pair.

        Returns
        -------
        float
        """
        return _lib.api().coordination(self._live(), float(rmax), bins, type_i, type_j)

    def hbonds(
        self,
        cutoff: float = 3.5,
        *,
        oxygen_type: int = 1,
        hydrogen_type: int = 2,
        dist: float = 2.42,
        angle: float = 30.0,
    ) -> int:
        """Hydrogen-bond count.

        Parameters
        ----------
        cutoff : float
            Heavy-atom neighbour cutoff, angstrom.
        oxygen_type : int
            Heavy-atom species. Zero keeps every atom that is not
            ``hydrogen_type``.
        hydrogen_type : int
            Positive hydrogen species.
        dist : float
            Acceptor-H distance, angstrom.
        angle : float
            O-O-H angle, degrees.

        Returns
        -------
        int
        """
        return _lib.api().hbonds(self._live(), float(cutoff), oxygen_type, hydrogen_type, float(dist), float(angle))

    def density(self, *, type: int = 0, bins: int = 0, axis: int = 2) -> Density:
        """Number density along one cell axis.

        Parameters
        ----------
        type : int
            Species filter. Zero keeps every atom.
        bins : int
            Bin count. Zero selects about one bin per 0.1 angstrom.
        axis : int
            0, 1, or 2.

        Returns
        -------
        Density
        """
        z, rho = _lib.api().density(self._live(), type, bins, axis)
        return Density(z, rho)

    def fingerprint(
        self,
        cutoff: float = 3.5,
        *,
        type: int = 0,
        k: int = 4,
        graph: int = CUTOFF,
        hops: int = 2,
        colour_types: bool = False,
    ) -> Fingerprint:
        """Topological fingerprint of the bonded graph.

        Parameters
        ----------
        cutoff : float
            Positive cutoff, angstrom.
        type : int
            Species filter.
        k : int
            Neighbour count.
        graph : int
            ``CUTOFF`` or ``MUTUAL``.
        hops : int
            Neighbourhood depth. Zero selects 2.
        colour_types : bool
            Colour vertices by LAMMPS type.

        Returns
        -------
        Fingerprint
        """
        ptr = _lib.api().fingerprint(
            self._live(), float(cutoff), type, k, graph, hops, 1 if colour_types else 0
        )
        return Fingerprint(ptr, _lib.api().destroy_fingerprint)

    def contact_pairs(self, spec: str) -> Contacts:
        """Mutual nearest unlike contacts of a site table.

        Parameters
        ----------
        spec : str
            Text such as ``"1=cationHead,2=anion"``.

        Returns
        -------
        Contacts
        """
        n_pairs, n_cation, n_anion = _lib.api().pairs(self._live(), spec.encode())
        return Contacts(n_pairs, n_cation, n_anion)

    def domains(self, spec: str, cutoff: float = 3.5, *, subset: int = 0) -> Domains:
        """Largest domain of a polar or apolar site subset.

        Parameters
        ----------
        spec : str
            Site-table text.
        cutoff : float
            Positive cutoff, angstrom.
        subset : int
            0 is apolar, 1 is polar.

        Returns
        -------
        Domains
        """
        n_subset, largest, percolation = _lib.api().domains(self._live(), float(cutoff), spec.encode(), subset)
        return Domains(n_subset, largest, percolation)

    def prisms(self, cutoff: float = 3.5, *, type: int = 0) -> Prisms:
        """Perfect and imperfect prism counts by ring size.

        Parameters
        ----------
        cutoff : float
            Positive cutoff, angstrom.
        type : int
            Species filter.

        Returns
        -------
        Prisms
        """
        perfect, imperfect = _lib.api().prisms(self._live(), float(cutoff), type)
        return Prisms(perfect, imperfect)

    def monolayer(self, cutoff: float = 3.5, *, type: int = 0, max_depth: int = 6, sheet_area: float = 100.0) -> npt.NDArray[np.int32]:
        """Primitive-ring histogram of a monolayer.

        Parameters
        ----------
        cutoff : float
            Positive cutoff, angstrom.
        type : int
            Species filter.
        max_depth : int
            Largest ring size.
        sheet_area : float
            Sheet area, square angstrom. Must be positive.

        Returns
        -------
        numpy.ndarray
            ``counts[k]`` is the number of rings of size ``k``.
        """
        return _lib.api().monolayer(self._live(), float(cutoff), type, max_depth, float(sheet_area))

    def site_table(self, spec: str) -> npt.NDArray[np.int32]:
        """Per-atom site kind.

        Parameters
        ----------
        spec : str
            Site-table text such as ``"1=waterO"``.

        Returns
        -------
        numpy.ndarray
            One ``seams_site_kind`` per atom.
        """
        return _lib.api().site_table(self._live(), spec.encode(), self.count())


def version() -> str:
    """Library version text.

    Returns
    -------
    str
        For example ``"2.14.0"``. The caller does not free the C string.
        This function returns a Python copy.
    """
    return _lib.api().version()


def version_tuple() -> tuple[int, int, int]:
    """Major, minor, and patch.

    Returns
    -------
    tuple of int
    """
    return _lib.api().version_parts()


def abi() -> int:
    """ABI number compiled into this library.

    Returns
    -------
    int
        ``SEAMS_ABI``, which is 1.
    """
    return _lib.api().abi()


def set_num_threads(n: int) -> None:
    """Set the OpenMP team used by later analysis calls.

    Parameters
    ----------
    n : int
        Positive width, or 0 to restore the width captured at first use.

    Raises
    ------
    seams.SeamsError
        If ``n`` is negative.

    Notes
    -----
    The width is process-global. Do not call this concurrently with analysis.
    """
    _lib.api().set_threads(n)


def num_threads() -> int:
    """Team width the next analysis call will request.

    Returns
    -------
    int
    """
    return _lib.api().threads()


def nframes(path: str | os.PathLike[str]) -> int:
    """Number of frames in a trajectory path.

    Parameters
    ----------
    path : path-like

    Returns
    -------
    int
    """
    return _lib.api().nframes(os.fsencode(path))


def label_name(code: int) -> str:
    """Name of a CHILL label.

    Parameters
    ----------
    code : int
        A ``seams_label`` value.

    Returns
    -------
    str
        ``"unknown"`` when the code is not a shipped label.
    """
    return _lib.api().label_name(code)


def status_string(code: int) -> str:
    """Short name of a status code.

    Parameters
    ----------
    code : int

    Returns
    -------
    str
    """
    return _lib.api().status_string(code)


def from_pydseams(frame: Any) -> System:
    """Copy a pydseamslib frame into a :class:`System`.

    ``pydseamslib`` (``import pydseams``) is the trajectory package. A
    frame exposes ``positions`` as ``(x, y, z)`` rows and ``box`` as
    three lengths or the LAMMPS restricted-triclinic numbers. This
    package is the ctypes binding of ``libseams``. The in-tree nanobind
    module ``dseams`` covers ``chill_plus``, ``q4q6``, neighbours, and
    rings. Use this function when a pydseams frame should run an
    analysis that module does not wrap.

    Parameters
    ----------
    frame : object
        Any object with ``positions`` and ``box``.

    Returns
    -------
    System
        A new handle. The caller closes it.

    Raises
    ------
    TypeError
        If ``positions`` or ``box`` is missing.
    """
    positions = getattr(frame, "positions", None)
    box = getattr(frame, "box", None)
    if positions is None or box is None:
        raise TypeError("pydseams frame needs positions and box")
    sys = System()
    try:
        sys.set_frame(positions, box)
    except Exception:
        sys.close()
        raise
    return sys
