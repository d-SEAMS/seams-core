# Copyright (c) d-SEAMS developers
"""Owned NumPy results. The library does not keep these buffers."""

from __future__ import annotations

from dataclasses import dataclass

import numpy as np
import numpy.typing as npt

F64 = npt.NDArray[np.float64]
I32 = npt.NDArray[np.int32]


@dataclass(frozen=True, slots=True)
class ChillResult:
    """Per-atom labels and the nine population counts.

    Attributes
    ----------
    labels : numpy.ndarray
        One ``seams_label`` code per atom.
    counts : numpy.ndarray
        Length 9, in label order. ``counts[i]`` is the population of
        label ``i``.
    """

    labels: I32
    counts: I32


@dataclass(frozen=True, slots=True)
class OrderResult:
    """Local Steinhardt parameters.

    Attributes
    ----------
    ql : numpy.ndarray
        One value per atom. Dimensionless.
    qlbar : numpy.ndarray or None
        Lechner-Dellago neighbour average, or ``None`` when it was not
        requested.
    """

    ql: F64
    qlbar: F64 | None


@dataclass(frozen=True, slots=True)
class Q4Q6:
    """Local q4 and q6 on the cutoff shell.

    Attributes
    ----------
    q4, q6 : numpy.ndarray
        One dimensionless value per atom.
    """

    q4: F64
    q6: F64


@dataclass(frozen=True, slots=True)
class Q6Bar:
    """Lechner-Dellago q6-bar and the summary above a cut.

    Attributes
    ----------
    values : numpy.ndarray or None
        Per-atom q6-bar when it was requested.
    mean, median : float
        Over finite values. An empty finite set is 0.
    n_ice : int
        Count of finite values strictly above the cut.
    """

    values: F64 | None
    mean: float
    median: float
    n_ice: int


@dataclass(frozen=True, slots=True)
class Signature:
    """Closed cages of one ring census.

    Attributes
    ----------
    n_cages, n_atoms : int
        Closed cages, and the distinct vertices in those cages.
    """

    n_cages: int
    n_atoms: int


@dataclass(frozen=True, slots=True)
class CageCounts:
    """Disjoint cage populations. The four counts sum to the kept atoms.

    Attributes
    ----------
    n_hex, n_cubic, n_water, n_both : int
        Hexagonal only, double-diamond only, neither, both.
    """

    n_hex: int
    n_cubic: int
    n_water: int
    n_both: int


@dataclass(frozen=True, slots=True)
class NeighbourCSR:
    """Caller-owned CSR. Column indices omit the self entry.

    Attributes
    ----------
    offsets : numpy.ndarray
        Length ``n + 1``.
    columns : numpy.ndarray
        Neighbour indices, length ``nnz``.
    """

    offsets: I32
    columns: I32


@dataclass(frozen=True, slots=True)
class RingLabels:
    """Size and symmetry of each primitive ring.

    Attributes
    ----------
    sizes, symmetry : numpy.ndarray
        ``symmetry`` holds ``seams_ring_symmetry`` codes.
    """

    sizes: I32
    symmetry: I32


@dataclass(frozen=True, slots=True)
class Radial:
    """Partial radial distribution.

    Attributes
    ----------
    r : numpy.ndarray
        Bin centres, angstrom.
    g : numpy.ndarray
        Dimensionless g(r).
    count : numpy.ndarray
        Pair counts.
    volume : float
        Cell volume, cubic angstrom.
    """

    r: F64
    g: F64
    count: I32
    volume: float


@dataclass(frozen=True, slots=True)
class RunningCN:
    """Running coordination number.

    Attributes
    ----------
    r, cn : numpy.ndarray
        ``r`` is angstrom. ``cn`` is dimensionless.
    """

    r: F64
    cn: F64


@dataclass(frozen=True, slots=True)
class Density:
    """Number density along one cell axis.

    Attributes
    ----------
    z, rho : numpy.ndarray
        ``z`` is angstrom. ``rho`` is number density.
    """

    z: F64
    rho: F64


@dataclass(frozen=True, slots=True)
class Contacts:
    """Mutual nearest unlike contacts.

    Attributes
    ----------
    n_pairs, n_cation, n_anion : int
        Pair count, and molecule counts in the ion cloud.
    """

    n_pairs: int
    n_cation: int
    n_anion: int


@dataclass(frozen=True, slots=True)
class Domains:
    """Largest polar or apolar domain.

    Attributes
    ----------
    n_subset, largest : int
        Atoms in the subset, and the largest domain.
    percolation : float
        Dimensionless percolation score.
    """

    n_subset: int
    largest: int
    percolation: float


@dataclass(frozen=True, slots=True)
class Prisms:
    """Perfect and imperfect prism counts. Index is the ring size.

    Attributes
    ----------
    perfect, imperfect : numpy.ndarray
        Length 7. Sizes 3 through 6 are filled.
    """

    perfect: I32
    imperfect: I32


@dataclass(frozen=True, slots=True)
class BrineFronts:
    """The two planar crossings of an ice slab.

    Attributes
    ----------
    z_enter, z_exit, thickness : float
        Angstrom, in the frame cell.
    found : bool
        True when both crossings exist.
    """

    z_enter: float
    z_exit: float
    thickness: float
    found: bool


@dataclass(frozen=True, slots=True)
class BrineRejection:
    """Molal rejection coefficient ``k = c_ice / c_liquid``.

    Attributes
    ----------
    c_ice, c_liquid, k : float
        Molality is formula units per kilogram of water.
    defined : bool
        True when the ratio is defined.
    """

    c_ice: float
    c_liquid: float
    k: float
    defined: bool


@dataclass(frozen=True, slots=True)
class BrineProfile:
    """Concentration profile along the front axis.

    Attributes
    ----------
    z_mid : numpy.ndarray
        Bin centres, angstrom.
    molality : numpy.ndarray
        mol/kg.
    ice_fraction : numpy.ndarray
        Dimensionless.
    """

    z_mid: F64
    molality: F64
    ice_fraction: F64


@dataclass(frozen=True, slots=True)
class BrinePockets:
    """Brine pockets trapped in the ice.

    Attributes
    ----------
    n_pockets, reservoir_waters : int
        Pocket count, and waters in the liquid reservoir.
    """

    n_pockets: int
    reservoir_waters: int
