# Copyright (c) d-SEAMS developers
"""Multi-frame TrajectorySequence, PBC unwrapping, molecular dipoles, and TCFs."""

from __future__ import annotations

import ctypes
from collections.abc import Iterable, Iterator, Sequence
from contextlib import contextmanager
from dataclasses import dataclass
from pathlib import Path

import numpy as np

from . import _capi

_CHEMFILES_EXTS = {
    ".pdb",
    ".gro",
    ".dcd",
    ".xtc",
    ".trr",
    ".nc",
    ".netcdf",
    ".tng",
    ".mmtf",
}
_LAMMPS_EXTS = {".lammpstrj", ".dump", ".lammps", ".trj"}

@dataclass(frozen=True)
class MoleculeCentresResult:
    """Per-molecule intact periodic centres for a single frame."""

    mol_ids: np.ndarray
    centres: np.ndarray


@dataclass(frozen=True)
class UnwrappedCentresResult:
    """Multi-frame unwrapped periodic centres of molecules."""

    mol_ids: np.ndarray
    centres: np.ndarray


@dataclass(frozen=True)
class WaterVectors:
    """Per-water orientation vectors under minimum image."""

    mol_ids: np.ndarray
    oh_mid: np.ndarray
    bisector: np.ndarray
    hh: np.ndarray
    hh_unit: np.ndarray


@dataclass(frozen=True)
class MoleculeDipoles:
    """Per-molecule dipole or geometric orientation vectors."""

    mol_ids: np.ndarray
    dipoles: np.ndarray
    unit_vectors: np.ndarray
    magnitudes: np.ndarray


@dataclass(frozen=True)
class OrientationOrder:
    """Axis-projected orientational order parameters P1 and P2."""

    cos_theta: np.ndarray
    p2: np.ndarray
    mean_p1: float
    mean_p2: float


@dataclass(frozen=True)
class DiffusionFit:
    """Linear Einstein diffusion fit over total, axial, and lateral MSD."""

    d_total: float
    d_axial: float
    d_lateral: float
    alpha: float
    """Log-log exponent; NaN when the window has fewer than two lags with
    t > 0 and MSD > 1e-20, so the log-log fit was skipped."""

    @property
    def alpha_fitted(self) -> bool:
        """False when the log-log fit was skipped and ``alpha`` is NaN."""
        return not np.isnan(self.alpha)

    @property
    def d_total_cm2_s(self) -> float:
        """Convert Angstrom^2 / ps to cm^2 / s (1 A^2/ps = 1e-4 cm^2/s)."""
        return self.d_total * 1.0e-4

    @property
    def d_axial_cm2_s(self) -> float:
        return self.d_axial * 1.0e-4

    @property
    def d_lateral_cm2_s(self) -> float:
        return self.d_lateral * 1.0e-4


@dataclass(frozen=True)
class MsdResult:
    """Multi-origin total, axial, and lateral MSD curves."""

    lag_time: np.ndarray
    total: np.ndarray
    axial: np.ndarray
    lateral: np.ndarray
    n_origins: np.ndarray
    diffusion: DiffusionFit | None = None

    def fit(self, t_start: float = 0.0, t_end: float = 0.0) -> DiffusionFit:
        return fit_diffusion(
            self.lag_time,
            self.total,
            self.axial,
            self.lateral,
            t_start=t_start,
            t_end=t_end,
        )


@dataclass(frozen=True)
class HbondLifetimeFit:
    """Integrated and single-exponential lifetimes for C_HB and S_HB."""

    tau_c_int: float
    tau_c_exp: float
    tau_s_int: float
    tau_s_exp: float


@dataclass(frozen=True)
class HbondLifetimeResult:
    """Intermittent C_HB(t) and continuous S_HB(t) lifetime correlation curves."""

    lag_time: np.ndarray
    c_intermittent: np.ndarray
    s_continuous: np.ndarray
    n_origins: np.ndarray
    lifetimes: HbondLifetimeFit | None = None

    @property
    def c_hb(self) -> np.ndarray:
        return self.c_intermittent

    @property
    def s_hb(self) -> np.ndarray:
        return self.s_continuous

    def fit(self, t_start: float = 0.0, t_end: float = 0.0) -> HbondLifetimeFit:
        c_int, c_exp = fit_exponential_tau(
            self.lag_time, self.c_intermittent, t_start=t_start, t_end=t_end
        )
        s_int, s_exp = fit_exponential_tau(
            self.lag_time, self.s_continuous, t_start=t_start, t_end=t_end
        )
        return HbondLifetimeFit(
            tau_c_int=c_int, tau_c_exp=c_exp, tau_s_int=s_int, tau_s_exp=s_exp
        )


@dataclass(frozen=True)
class ReorientationFit:
    """Integrated and single-exponential reorientation times tau_1 and tau_2."""

    tau1_int: float
    tau1_exp: float
    tau2_int: float
    tau2_exp: float


@dataclass(frozen=True)
class ReorientationResult:
    """First and second Legendre orientational correlation curves C1(t) and C2(t)."""

    lag_time: np.ndarray
    c1: np.ndarray
    c2: np.ndarray
    n_origins: np.ndarray
    relaxation: ReorientationFit | None = None

    def fit(self, t_start: float = 0.0, t_end: float = 0.0) -> ReorientationFit:
        t1_int, t1_exp = fit_exponential_tau(
            self.lag_time, self.c1, t_start=t_start, t_end=t_end
        )
        t2_int, t2_exp = fit_exponential_tau(
            self.lag_time, self.c2, t_start=t_start, t_end=t_end
        )
        return ReorientationFit(
            tau1_int=t1_int, tau1_exp=t1_exp, tau2_int=t2_int, tau2_exp=t2_exp
        )


class CtypesFrame:
    """Lightweight single-frame container backed directly by libseams."""

    def __init__(
        self,
        positions: np.ndarray,
        box9: np.ndarray,
        types: np.ndarray | None = None,
        atom_ids: np.ndarray | None = None,
        mol_ids: np.ndarray | None = None,
        cutoff: float = 3.5,
    ):
        xyz = np.ascontiguousarray(positions, dtype=np.float64)
        if xyz.ndim != 2 or xyz.shape[1] != 3 or xyz.shape[0] == 0:
            raise ValueError("positions must have shape (N, 3) with N >= 1")
        n = xyz.shape[0]
        b9 = _capi.normalize_box9(box9)
        t_arr = (
            np.ones(n, dtype=np.int32)
            if types is None
            else np.ascontiguousarray(types, dtype=np.int32)
        )
        id_arr = (
            np.arange(1, n + 1, dtype=np.int32)
            if atom_ids is None
            else np.ascontiguousarray(atom_ids, dtype=np.int32)
        )
        m_arr = (
            id_arr.copy()
            if mol_ids is None
            else np.ascontiguousarray(mol_ids, dtype=np.int32)
        )
        if t_arr.shape != (n,) or id_arr.shape != (n,) or m_arr.shape != (n,):
            raise ValueError("types, atom_ids, and mol_ids must have shape (N,)")
        self._xyz = xyz
        self._box9 = b9
        self._types = t_arr
        self._atom_ids = id_arr
        self._mol_ids = m_arr
        self.cutoff = float(cutoff)
        self._seams_raw_data = (xyz, b9, t_arr, id_arr, m_arr)

    @property
    def n_atoms(self) -> int:
        return int(self._xyz.shape[0])

    @property
    def positions(self) -> np.ndarray:
        return self._xyz.copy()

    @property
    def box9(self) -> np.ndarray:
        return self._box9.copy()

    @property
    def box(self) -> list[float]:
        dump_box, _ = _capi.box9_to_dump_bounds(self._box9)
        return dump_box

    @property
    def types(self) -> np.ndarray:
        return self._types.copy()

    @property
    def mol_ids(self) -> np.ndarray:
        return self._mol_ids.copy()

    def chill_plus(self, cutoff: float | None = None, species: int = 0) -> dict[str, int]:
        cut = self.cutoff if cutoff is None else float(cutoff)
        names = (
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
        with _system_handle_from_arrays(
            self._xyz, self._box9, self._types, self._atom_ids, self._mol_ids
        ) as sys_h:
            labels = np.zeros(self.n_atoms, dtype=np.int32)
            counts = np.zeros(9, dtype=np.int32)
            lib = _capi.get_lib()
            st = lib.seams_chill_plus_n(
                sys_h,
                cut,
                int(species),
                0,
                0,
                _capi.int_ptr(labels),
                self.n_atoms,
                _capi.int_ptr(counts),
            )
            _capi.check_status(st, sys_h)
        return {names[i]: int(counts[i]) for i in range(9)}


def _extract_raw_from_system(sys_h) -> tuple[np.ndarray, np.ndarray, np.ndarray, np.ndarray, np.ndarray]:
    lib = _capi.get_lib()
    n_c = ctypes.c_size_t(0)
    _capi.check_status(lib.seams_system_natoms(sys_h, ctypes.byref(n_c)), sys_h)
    n = int(n_c.value)
    xyz = np.empty((n, 3), dtype=np.float64)
    box9 = np.empty(9, dtype=np.float64)
    types = np.empty(n, dtype=np.int32)
    ids = np.empty(n, dtype=np.int32)
    mols = np.empty(n, dtype=np.int32)
    _capi.check_status(
        lib.seams_system_positions(sys_h, _capi.double_ptr(xyz), 3 * n), sys_h
    )
    _capi.check_status(lib.seams_system_box(sys_h, _capi.double_ptr(box9), 9), sys_h)
    _capi.check_status(lib.seams_system_types(sys_h, _capi.int_ptr(types), n), sys_h)
    _capi.check_status(lib.seams_system_ids(sys_h, _capi.int_ptr(ids), n), sys_h)
    _capi.check_status(lib.seams_mols(sys_h, _capi.int_ptr(mols), n), sys_h)
    return xyz, box9, types, ids, mols


def _extract_frame_data(
    frame,
) -> tuple[np.ndarray, np.ndarray, np.ndarray, np.ndarray, np.ndarray]:
    """Extract (xyz (N,3), box9 (9,), types (N,), atom_ids (N,), mol_ids (N,))."""
    raw = getattr(frame, "_seams_raw_data", None)
    if raw is not None:
        return raw

    cloud = getattr(frame, "cloud", None)
    if cloud is not None:
        n = int(cloud.nop)
        pts = cloud.pts
        xyz = np.empty((n, 3), dtype=np.float64)
        types = np.empty(n, dtype=np.int32)
        ids = np.empty(n, dtype=np.int32)
        mols = np.empty(n, dtype=np.int32)
        for i in range(n):
            pt = pts[i]
            xyz[i, 0] = pt.x
            xyz[i, 1] = pt.y
            xyz[i, 2] = pt.z
            types[i] = pt.c_type
            ids[i] = pt.atomID
            mols[i] = pt.molID
        box9 = _capi.normalize_box9(list(cloud.box), list(cloud.boxLow))
        raw = (xyz, box9, types, ids, mols)
        try:
            frame._seams_raw_data = raw
        except AttributeError:
            pass
        return raw

    raise TypeError(f"unsupported frame object: {type(frame)!r}")


@contextmanager
def _system_handle_from_arrays(
    xyz: np.ndarray,
    box9: np.ndarray,
    types: np.ndarray | None = None,
    atom_ids: np.ndarray | None = None,
    mol_ids: np.ndarray | None = None,
):
    lib = _capi.get_lib()
    xyz_c = np.ascontiguousarray(xyz, dtype=np.float64)
    if xyz_c.ndim != 2 or xyz_c.shape[1] != 3 or xyz_c.shape[0] == 0:
        raise ValueError("positions must have shape (N, 3) with N >= 1")
    n = xyz_c.shape[0]
    b9_c = np.ascontiguousarray(_capi.normalize_box9(box9), dtype=np.float64)
    t_c = (
        None
        if types is None
        else np.ascontiguousarray(types, dtype=np.int32).reshape(n)
    )
    id_c = (
        None
        if atom_ids is None
        else np.ascontiguousarray(atom_ids, dtype=np.int32).reshape(n)
    )
    sys_h = lib.seams_system_create()
    if not sys_h:
        raise MemoryError("seams_system_create failed")
    try:
        st = lib.seams_system_set_atoms(
            sys_h,
            _capi.double_ptr(xyz_c),
            n,
            _capi.double_ptr(b9_c),
            9,
            None if t_c is None else _capi.int_ptr(t_c),
            None if id_c is None else _capi.int_ptr(id_c),
        )
        _capi.check_status(st, sys_h)
        if mol_ids is not None:
            m_c = np.ascontiguousarray(mol_ids, dtype=np.int32).reshape(n)
            st = lib.seams_set_mols(sys_h, _capi.int_ptr(m_c), n)
            _capi.check_status(st, sys_h)
        yield sys_h
    finally:
        lib.seams_system_destroy(sys_h)


@contextmanager
def _system_handle_from_frame(frame):
    xyz, box9, types, ids, mols = _extract_frame_data(frame)
    with _system_handle_from_arrays(xyz, box9, types, ids, mols) as sys_h:
        yield sys_h


def _make_frame_from_raw(
    xyz: np.ndarray,
    box9: np.ndarray,
    types: np.ndarray,
    atom_ids: np.ndarray,
    mol_ids: np.ndarray,
    *,
    atom_type: int | None = None,
    cutoff: float | None = None,
    bonded: str = "cutoff",
):
    import pydseams as _pkg

    yoda_mod = getattr(_pkg, "yoda", None)
    frame_cls = getattr(_pkg, "Frame", CtypesFrame)
    if yoda_mod is not None and frame_cls is not CtypesFrame:
        cloud = yoda_mod.PointCloudDouble()
        n = int(xyz.shape[0])
        pts = []
        for i in range(n):
            pt = yoda_mod.PointDouble()
            pt.x = float(xyz[i, 0])
            pt.y = float(xyz[i, 1])
            pt.z = float(xyz[i, 2])
            pt.c_type = int(types[i])
            pt.atomID = int(atom_ids[i])
            pt.molID = int(mol_ids[i])
            pt.inSlice = True
            pts.append(pt)
        cloud.pts = pts
        cloud.nop = n
        dump_box, dump_low = _capi.box9_to_dump_bounds(box9)
        cloud.box = dump_box
        cloud.boxLow = dump_low
        id_map = {int(atom_ids[i]): i for i in range(n)}
        cloud.idIndexMap = id_map
        chosen_type = int(atom_type) if atom_type is not None else int(types[0])
        fr = frame_cls(
            atom_type=chosen_type,
            cutoff=cutoff,
            bonded=bonded,
            cloud=cloud,
        )
        fr._seams_raw_data = (
            xyz.copy(),
            box9.copy(),
            types.copy(),
            atom_ids.copy(),
            mol_ids.copy(),
        )
        return fr

    return CtypesFrame(
        xyz,
        box9,
        types=types,
        atom_ids=atom_ids,
        mol_ids=mol_ids,
        cutoff=3.5 if cutoff is None else float(cutoff),
    )


def _count_path_frames(path: Path) -> int:
    if not path.is_file():
        raise FileNotFoundError(str(path))
    lib = _capi.get_lib()
    n_c = ctypes.c_size_t(0)
    st = lib.seams_path_nframes(str(path).encode("utf-8"), ctypes.byref(n_c))
    _capi.check_status(st)
    return int(n_c.value)


def _load_frame_at_step(
    path: Path,
    frame_1based: int,
    *,
    atom_type: int | None = None,
    all_atoms: bool = False,
    cutoff: float | None = None,
    bonded: str = "auto",
):
    import pydseams as _pkg

    frame_cls = getattr(_pkg, "Frame", CtypesFrame)
    yoda_mod = getattr(_pkg, "yoda", None)
    suffix = path.suffix.lower()

    if yoda_mod is not None and frame_cls is not CtypesFrame:
        if suffix in _LAMMPS_EXTS or (
            suffix not in _CHEMFILES_EXTS and suffix not in {".xyz", ".con"}
        ):
            kw = {"frame": frame_1based, "bonded": bonded, "all_atoms": all_atoms}
            if atom_type is not None:
                kw["atom_type"] = atom_type
            if cutoff is not None:
                kw["cutoff"] = cutoff
            return frame_cls.from_file(str(path), **kw)
        if suffix == ".con":
            kw = {"frame": frame_1based}
            if bonded != "auto":
                kw["bonded"] = bonded
            if atom_type is not None:
                kw["atom_type"] = atom_type
            if cutoff is not None:
                kw["cutoff"] = cutoff
            return frame_cls.from_con(str(path), **kw)
        if suffix == ".xyz" and frame_1based == 1 and not all_atoms:
            kw = {}
            if bonded != "auto":
                kw["bonded"] = bonded
            if atom_type is not None:
                kw["atom_type"] = atom_type
            if cutoff is not None:
                kw["cutoff"] = cutoff
            return frame_cls.from_xyz(str(path), **kw)

    type_filter = -1 if (all_atoms or atom_type is None) else int(atom_type)
    lib = _capi.get_lib()
    sys_h = lib.seams_system_create()
    if not sys_h:
        raise MemoryError("seams_system_create failed")
    try:
        st = lib.seams_system_read(
            sys_h, str(path).encode("utf-8"), int(frame_1based), int(type_filter)
        )
        _capi.check_status(st, sys_h)
        xyz, box9, types, ids, mols = _extract_raw_from_system(sys_h)
    finally:
        lib.seams_system_destroy(sys_h)

    use_bonded = "cutoff" if bonded == "auto" else bonded
    return _make_frame_from_raw(
        xyz,
        box9,
        types,
        ids,
        mols,
        atom_type=atom_type,
        cutoff=cutoff,
        bonded=use_bonded,
    )


def molecule_centres(
    frame_or_positions,
    cell=None,
    *,
    mol_ids=None,
    masses=None,
) -> MoleculeCentresResult:
    """Compute intact periodic centres of mass or geometry for one frame."""
    lib = _capi.get_lib()
    if cell is None and not isinstance(frame_or_positions, np.ndarray):
        xyz, box9, types, ids, stored_mols = _extract_frame_data(frame_or_positions)
        use_mols = stored_mols if mol_ids is None else np.asarray(mol_ids, dtype=np.int32)
    else:
        if cell is None:
            raise ValueError("cell is required when frame_or_positions is an array")
        xyz = np.ascontiguousarray(frame_or_positions, dtype=np.float64)
        if xyz.ndim != 2 or xyz.shape[1] != 3:
            raise ValueError("positions must have shape (N, 3)")
        box9 = _capi.normalize_box9(cell)
        types = None
        ids = None
        use_mols = (
            np.arange(1, xyz.shape[0] + 1, dtype=np.int32)
            if mol_ids is None
            else np.ascontiguousarray(mol_ids, dtype=np.int32)
        )

    n = int(xyz.shape[0])
    m_c = np.ascontiguousarray(use_mols, dtype=np.int32).reshape(n)
    w_c = (
        None
        if masses is None
        else np.ascontiguousarray(masses, dtype=np.float64).reshape(n)
    )

    with _system_handle_from_arrays(xyz, box9, types, ids, m_c) as sys_h:
        nm_c = ctypes.c_size_t(0)
        st = lib.seams_molecule_centres(
            sys_h,
            _capi.int_ptr(m_c),
            n,
            None if w_c is None else _capi.double_ptr(w_c),
            0 if w_c is None else n,
            None,
            None,
            0,
            ctypes.byref(nm_c),
        )
        _capi.check_status(st, sys_h)
        nm = int(nm_c.value)
        out_ids = np.empty(nm, dtype=np.int32)
        out_centres = np.empty((nm, 3), dtype=np.float64)
        st = lib.seams_molecule_centres(
            sys_h,
            _capi.int_ptr(m_c),
            n,
            None if w_c is None else _capi.double_ptr(w_c),
            0 if w_c is None else n,
            _capi.int_ptr(out_ids),
            _capi.double_ptr(out_centres),
            nm,
            ctypes.byref(nm_c),
        )
        _capi.check_status(st, sys_h)
    return MoleculeCentresResult(mol_ids=out_ids, centres=out_centres)


def _normalize_boxes_sequence(cell, nframes: int) -> np.ndarray:
    arr = np.asarray(cell, dtype=np.float64)
    if arr.shape in {(3,), (6,), (9,), (3, 3)}:
        return np.ascontiguousarray(_capi.normalize_box9(arr)[None, :], dtype=np.float64)
    if arr.ndim >= 2 and arr.shape[0] == nframes:
        rows = [_capi.normalize_box9(arr[i]) for i in range(nframes)]
        return np.ascontiguousarray(np.stack(rows, axis=0), dtype=np.float64)
    raise ValueError(
        f"cell must specify 1 box or {nframes} per-frame boxes; got shape {arr.shape}"
    )


def unwrap_positions(positions_or_traj, cell=None) -> np.ndarray:
    """Unwrap periodic coordinates across frames (orthorhombic or triclinic)."""
    if isinstance(positions_or_traj, TrajectorySequence):
        return positions_or_traj.unwrap_positions()
    if cell is None:
        raise ValueError("cell is required when unwrapping a coordinate array")
    xyz = np.ascontiguousarray(positions_or_traj, dtype=np.float64)
    if xyz.ndim != 3 or xyz.shape[2] != 3:
        raise ValueError("positions must have shape (n_frames, n_atoms, 3)")
    nframes, natoms, _ = xyz.shape
    boxes = _normalize_boxes_sequence(cell, nframes)
    out = np.empty_like(xyz)
    lib = _capi.get_lib()
    st = lib.seams_traj_unwrap_positions(
        _capi.double_ptr(xyz),
        nframes,
        natoms,
        _capi.double_ptr(boxes),
        boxes.shape[0],
        _capi.double_ptr(out),
        out.size,
    )
    _capi.check_status(st)
    return out


def unwrap_centres(
    positions_or_traj,
    cell=None,
    *,
    mol_ids=None,
    masses=None,
) -> UnwrappedCentresResult:
    """Assemble intact periodic molecule centres per frame and unwrap across time."""
    if isinstance(positions_or_traj, TrajectorySequence):
        return positions_or_traj.unwrap_centres(mol_ids=mol_ids, masses=masses)
    if cell is None:
        raise ValueError("cell is required when unwrapping a coordinate array")
    xyz = np.ascontiguousarray(positions_or_traj, dtype=np.float64)
    if xyz.ndim != 3 or xyz.shape[2] != 3:
        raise ValueError("positions must have shape (n_frames, n_atoms, 3)")
    nframes, natoms, _ = xyz.shape
    boxes = _normalize_boxes_sequence(cell, nframes)
    m_c = (
        np.arange(1, natoms + 1, dtype=np.int32)
        if mol_ids is None
        else np.ascontiguousarray(mol_ids, dtype=np.int32).reshape(natoms)
    )
    w_c = (
        None
        if masses is None
        else np.ascontiguousarray(masses, dtype=np.float64).reshape(natoms)
    )

    lib = _capi.get_lib()
    nm_c = ctypes.c_size_t(0)
    st = lib.seams_traj_unwrap_centres(
        _capi.double_ptr(xyz),
        nframes,
        natoms,
        _capi.double_ptr(boxes),
        boxes.shape[0],
        _capi.int_ptr(m_c),
        natoms,
        None if w_c is None else _capi.double_ptr(w_c),
        0 if w_c is None else natoms,
        None,
        0,
        None,
        0,
        ctypes.byref(nm_c),
    )
    _capi.check_status(st)
    nm = int(nm_c.value)
    out_ids = np.empty(nm, dtype=np.int32)
    out_centres = np.empty((nframes, nm, 3), dtype=np.float64)
    st = lib.seams_traj_unwrap_centres(
        _capi.double_ptr(xyz),
        nframes,
        natoms,
        _capi.double_ptr(boxes),
        boxes.shape[0],
        _capi.int_ptr(m_c),
        natoms,
        None if w_c is None else _capi.double_ptr(w_c),
        0 if w_c is None else natoms,
        _capi.int_ptr(out_ids),
        nm,
        _capi.double_ptr(out_centres),
        out_centres.size,
        ctypes.byref(nm_c),
    )
    _capi.check_status(st)
    return UnwrappedCentresResult(mol_ids=out_ids, centres=out_centres)


def _infer_water_types(types: np.ndarray, oxygen_type: int | None, hydrogen_type: int | None) -> tuple[int, int]:
    if oxygen_type is not None and hydrogen_type is not None:
        return int(oxygen_type), int(hydrogen_type)
    vals, counts = np.unique(types, return_counts=True)
    count_map = dict(zip(vals.tolist(), counts.tolist()))
    if oxygen_type is not None:
        ot = int(oxygen_type)
        for t, c in count_map.items():
            if t != ot and count_map.get(ot, 0) > 0 and c == 2 * count_map[ot]:
                return ot, t
        others = [t for t in vals.tolist() if t != ot]
        return ot, (others[0] if others else (1 if ot != 1 else 2))
    if 2 in count_map and 1 in count_map:
        if count_map[1] == 2 * count_map[2]:
            return 2, 1
        if count_map[2] == 2 * count_map[1]:
            return 1, 2
    if 8 in count_map and 1 in count_map:
        return 8, 1
    for t1 in vals.tolist():
        for t2 in vals.tolist():
            if t1 != t2 and count_map[t2] == 2 * count_map[t1]:
                return t1, t2
    return 2, 1


def water_vectors(
    frame,
    *,
    oxygen_type: int | None = None,
    hydrogen_type: int | None = None,
) -> WaterVectors:
    """Return O->H-midpoint and H1->H2 vectors (raw and unit) for water molecules."""
    xyz, box9, types, ids, mols = _extract_frame_data(frame)
    ot, ht = _infer_water_types(types, oxygen_type, hydrogen_type)
    lib = _capi.get_lib()
    with _system_handle_from_arrays(xyz, box9, types, ids, mols) as sys_h:
        nw_c = ctypes.c_size_t(0)
        st = lib.seams_water_vectors(
            sys_h, ot, ht, None, None, None, None, None, 0, ctypes.byref(nw_c)
        )
        _capi.check_status(st, sys_h)
        nw = int(nw_c.value)
        m_ids = np.empty(nw, dtype=np.int32)
        oh_mid = np.empty((nw, 3), dtype=np.float64)
        bisector = np.empty((nw, 3), dtype=np.float64)
        hh = np.empty((nw, 3), dtype=np.float64)
        hh_unit = np.empty((nw, 3), dtype=np.float64)
        if nw > 0:
            st = lib.seams_water_vectors(
                sys_h,
                ot,
                ht,
                _capi.int_ptr(m_ids),
                _capi.double_ptr(oh_mid),
                _capi.double_ptr(bisector),
                _capi.double_ptr(hh),
                _capi.double_ptr(hh_unit),
                nw,
                ctypes.byref(nw_c),
            )
            _capi.check_status(st, sys_h)
    return WaterVectors(
        mol_ids=m_ids,
        oh_mid=oh_mid,
        bisector=bisector,
        hh=hh,
        hh_unit=hh_unit,
    )


def molecule_dipoles(
    frame_or_positions,
    cell=None,
    *,
    mol_ids=None,
    charges=None,
) -> MoleculeDipoles:
    """Compute per-molecule charge-weighted dipoles or geometric orientation vectors."""
    if cell is None and not isinstance(frame_or_positions, np.ndarray):
        xyz, box9, types, ids, stored_mols = _extract_frame_data(frame_or_positions)
        use_mols = stored_mols if mol_ids is None else np.asarray(mol_ids, dtype=np.int32)
    else:
        if cell is None:
            raise ValueError("cell is required when frame_or_positions is an array")
        xyz = np.ascontiguousarray(frame_or_positions, dtype=np.float64)
        if xyz.ndim != 2 or xyz.shape[1] != 3:
            raise ValueError("positions must have shape (N, 3)")
        box9 = _capi.normalize_box9(cell)
        types = None
        ids = None
        use_mols = (
            np.arange(1, xyz.shape[0] + 1, dtype=np.int32)
            if mol_ids is None
            else np.ascontiguousarray(mol_ids, dtype=np.int32)
        )

    n = int(xyz.shape[0])
    m_c = np.ascontiguousarray(use_mols, dtype=np.int32).reshape(n)
    q_c = (
        None
        if charges is None
        else np.ascontiguousarray(charges, dtype=np.float64).reshape(n)
    )
    lib = _capi.get_lib()
    with _system_handle_from_arrays(xyz, box9, types, ids, m_c) as sys_h:
        nm_c = ctypes.c_size_t(0)
        st = lib.seams_molecule_dipoles(
            sys_h,
            _capi.int_ptr(m_c),
            n,
            None if q_c is None else _capi.double_ptr(q_c),
            0 if q_c is None else n,
            None,
            None,
            None,
            0,
            ctypes.byref(nm_c),
        )
        _capi.check_status(st, sys_h)
        nm = int(nm_c.value)
        out_ids = np.empty(nm, dtype=np.int32)
        dips = np.empty((nm, 3), dtype=np.float64)
        units = np.empty((nm, 3), dtype=np.float64)
        if nm > 0:
            st = lib.seams_molecule_dipoles(
                sys_h,
                _capi.int_ptr(m_c),
                n,
                None if q_c is None else _capi.double_ptr(q_c),
                0 if q_c is None else n,
                _capi.int_ptr(out_ids),
                _capi.double_ptr(dips),
                _capi.double_ptr(units),
                nm,
                ctypes.byref(nm_c),
            )
            _capi.check_status(st, sys_h)
    mags = np.linalg.norm(dips, axis=1) if nm > 0 else np.empty(0, dtype=np.float64)
    return MoleculeDipoles(
        mol_ids=out_ids,
        dipoles=dips,
        unit_vectors=units,
        magnitudes=mags,
    )


def orientation_order(vectors, axis=(0.0, 0.0, 1.0)) -> OrientationOrder:
    """Compute per-molecule cos(theta), P2(cos(theta)), and ensemble <P1>, <P2>."""
    vecs = np.ascontiguousarray(vectors, dtype=np.float64)
    if vecs.ndim != 2 or vecs.shape[1] != 3:
        raise ValueError("vectors must have shape (N, 3)")
    ax = np.ascontiguousarray(axis, dtype=np.float64).reshape(3)
    nv = int(vecs.shape[0])
    cos_t = np.empty(nv, dtype=np.float64)
    p2_arr = np.empty(nv, dtype=np.float64)
    mp1 = ctypes.c_double(0.0)
    mp2 = ctypes.c_double(0.0)
    lib = _capi.get_lib()
    st = lib.seams_orientation_order(
        None if nv == 0 else _capi.double_ptr(vecs),
        nv,
        _capi.double_ptr(ax),
        None if nv == 0 else _capi.double_ptr(cos_t),
        None if nv == 0 else _capi.double_ptr(p2_arr),
        nv,
        ctypes.byref(mp1),
        ctypes.byref(mp2),
    )
    _capi.check_status(st)
    return OrientationOrder(
        cos_theta=cos_t,
        p2=p2_arr,
        mean_p1=float(mp1.value),
        mean_p2=float(mp2.value),
    )


def fit_diffusion(
    lag_time,
    msd_total,
    msd_axial=None,
    msd_lateral=None,
    *,
    t_start: float = 0.0,
    t_end: float = 0.0,
) -> DiffusionFit:
    """Ordinary least-squares Einstein diffusion fit over [t_start, t_end]."""
    lt = np.asarray(lag_time, dtype=np.float64)
    tot = np.asarray(msd_total, dtype=np.float64)
    ax = np.zeros_like(tot) if msd_axial is None else np.asarray(msd_axial, dtype=np.float64)
    lat = tot if msd_lateral is None else np.asarray(msd_lateral, dtype=np.float64)
    if lt.ndim != 1 or lt.size < 2 or tot.shape != lt.shape or ax.shape != lt.shape or lat.shape != lt.shape:
        raise ValueError("lag_time and MSD arrays must be 1D of equal length >= 2")
    t_max = float(lt[-1])
    if t_start == 0.0 and t_end == 0.0:
        lo, hi = 0.1 * t_max, 0.5 * t_max
        if np.count_nonzero((lt >= lo) & (lt <= hi)) < 2:
            lo, hi = float(lt[1]), t_max
    else:
        if t_start < 0.0 or t_end <= t_start:
            raise ValueError("fit window must satisfy 0 <= t_start < t_end")
        lo, hi = float(t_start), float(t_end)
    mask = (
        (lt >= lo - 1e-12 * max(1.0, hi))
        & (lt <= hi + 1e-12 * max(1.0, hi))
        & np.isfinite(lt)
        & np.isfinite(tot)
        & np.isfinite(ax)
        & np.isfinite(lat)
    )
    if np.count_nonzero(mask) < 2:
        raise ValueError("fewer than two valid lag points in fit window")
    tx = lt[mask]
    s_tot = float(np.polyfit(tx, tot[mask], 1)[0])
    s_ax = float(np.polyfit(tx, ax[mask], 1)[0])
    s_lat = float(np.polyfit(tx, lat[mask], 1)[0])
    pos = mask & (lt > 0.0) & (tot > 1e-20)
    alpha = float(np.polyfit(np.log(lt[pos]), np.log(tot[pos]), 1)[0]) if np.count_nonzero(pos) >= 2 else float("nan")
    return DiffusionFit(
        d_total=s_tot / 6.0,
        d_axial=s_ax / 2.0,
        d_lateral=s_lat / 4.0,
        alpha=alpha,
    )


def fit_exponential_tau(
    lag_time,
    corr,
    *,
    t_start: float = 0.0,
    t_end: float = 0.0,
) -> tuple[float, float]:
    """Return (tau_integral, tau_exponential) for a normalized correlation curve."""
    lt = np.asarray(lag_time, dtype=np.float64)
    c = np.asarray(corr, dtype=np.float64)
    if lt.ndim != 1 or lt.size < 2 or c.shape != lt.shape:
        raise ValueError("lag_time and corr must be 1D arrays of length >= 2")

    integral = 0.0
    for i in range(1, lt.size):
        if not (np.isfinite(lt[i - 1]) and np.isfinite(lt[i]) and np.isfinite(c[i - 1]) and np.isfinite(c[i])):
            break
        dt_step = float(lt[i] - lt[i - 1])
        if dt_step <= 0.0 or c[i - 1] <= 0.0:
            break
        if c[i] <= 0.0:
            frac = float(c[i - 1] / (c[i - 1] - c[i]))
            integral += 0.5 * float(c[i - 1]) * dt_step * frac
            break
        integral += 0.5 * float(c[i - 1] + c[i]) * dt_step

    lo = 0.0 if (t_start == 0.0 and t_end == 0.0) else float(t_start)
    hi = float(lt[-1]) if (t_start == 0.0 and t_end == 0.0) else float(t_end)
    if lo < 0.0 or hi <= lo:
        raise ValueError("fit window must satisfy 0 <= t_start < t_end")

    tx_list: list[float] = []
    ly_list: list[float] = []
    w_list: list[float] = []
    for t_val, c_val in zip(lt.tolist(), c.tolist()):
        if t_val < lo - 1e-12 * max(1.0, hi) or t_val > hi + 1e-12 * max(1.0, hi):
            continue
        if not (np.isfinite(t_val) and np.isfinite(c_val)):
            continue
        if c_val <= 1e-6:
            if t_val > lo:
                break
            continue
        tx_list.append(t_val)
        ly_list.append(float(np.log(c_val)))
        w_list.append(c_val * c_val)

    tau_exp = float("nan")
    if len(tx_list) >= 2:
        tx_a = np.asarray(tx_list, dtype=np.float64)
        ly_a = np.asarray(ly_list, dtype=np.float64)
        w_a = np.asarray(w_list, dtype=np.float64)
        sw = float(np.sum(w_a))
        mx = float(np.sum(w_a * tx_a) / sw)
        my = float(np.sum(w_a * ly_a) / sw)
        sxx = float(np.sum(w_a * (tx_a - mx) ** 2))
        sxy = float(np.sum(w_a * (tx_a - mx) * (ly_a - my)))
        if sxx > 0.0:
            slope = sxy / sxx
            if slope < -1e-15:
                tau_exp = -1.0 / slope
    return integral, tau_exp


class MsdAccumulator:
    """Streaming multi-origin MSD accumulator backed by libseams."""

    def __init__(
        self,
        n_particles: int,
        *,
        max_lag: int = 0,
        origin_stride: int = 1,
        dt: float = 1.0,
        axis=(0.0, 0.0, 1.0),
    ):
        ax = np.ascontiguousarray(axis, dtype=np.float64).reshape(3)
        lib = _capi.get_lib()
        raw = ctypes.c_void_p()
        st = lib.seams_msd_create(
            int(n_particles),
            int(max_lag),
            int(origin_stride),
            float(dt),
            _capi.double_ptr(ax),
            ctypes.byref(raw),
        )
        self._handle: ctypes.c_void_p | None = None
        _capi.check_status(st)
        self._handle = raw
        self.n_particles = int(n_particles)

    def close(self) -> None:
        if self._handle is not None and self._handle.value:
            _capi.get_lib().seams_msd_destroy(self._handle)
            self._handle = None

    def __del__(self):
        try:
            self.close()
        except Exception:
            pass

    def __enter__(self) -> MsdAccumulator:
        return self

    def __exit__(self, exc_type, exc, tb) -> None:
        self.close()

    def push(self, positions, cell=None) -> None:
        if self._handle is None:
            raise RuntimeError("MsdAccumulator is closed")
        xyz = np.ascontiguousarray(positions, dtype=np.float64)
        if xyz.shape != (self.n_particles, 3):
            raise ValueError(
                f"positions must have shape ({self.n_particles}, 3); got {xyz.shape}"
            )
        lib = _capi.get_lib()
        if cell is None:
            st = lib.seams_msd_push(
                self._handle, _capi.double_ptr(xyz), self.n_particles, None, 0
            )
        else:
            b9 = _capi.normalize_box9(cell)
            st = lib.seams_msd_push(
                self._handle,
                _capi.double_ptr(xyz),
                self.n_particles,
                _capi.double_ptr(b9),
                9,
            )
        _capi.check_status(st, self._handle)

    def curve(self) -> MsdResult:
        if self._handle is None:
            raise RuntimeError("MsdAccumulator is closed")
        lib = _capi.get_lib()
        nl_c = ctypes.c_size_t(0)
        st = lib.seams_msd_curve(
            self._handle, None, None, None, None, None, 0, ctypes.byref(nl_c)
        )
        _capi.check_status(st, self._handle)
        nl = int(nl_c.value)
        lt = np.empty(nl, dtype=np.float64)
        tot = np.empty(nl, dtype=np.float64)
        ax = np.empty(nl, dtype=np.float64)
        lat = np.empty(nl, dtype=np.float64)
        orig = np.empty(nl, dtype=np.uintp)
        st = lib.seams_msd_curve(
            self._handle,
            _capi.double_ptr(lt),
            _capi.double_ptr(tot),
            _capi.double_ptr(ax),
            _capi.double_ptr(lat),
            _capi.sizet_ptr(orig),
            nl,
            ctypes.byref(nl_c),
        )
        _capi.check_status(st, self._handle)
        return MsdResult(
            lag_time=lt,
            total=tot,
            axial=ax,
            lateral=lat,
            n_origins=orig.astype(np.int64),
        )

    def fit(self, t_start: float = 0.0, t_end: float = 0.0) -> DiffusionFit:
        if self._handle is None:
            raise RuntimeError("MsdAccumulator is closed")
        lib = _capi.get_lib()
        dtot = ctypes.c_double(0.0)
        dax = ctypes.c_double(0.0)
        dlat = ctypes.c_double(0.0)
        alpha = ctypes.c_double(0.0)
        st = lib.seams_msd_fit(
            self._handle,
            float(t_start),
            float(t_end),
            ctypes.byref(dtot),
            ctypes.byref(dax),
            ctypes.byref(dlat),
            ctypes.byref(alpha),
        )
        _capi.check_status(st, self._handle)
        return DiffusionFit(
            d_total=float(dtot.value),
            d_axial=float(dax.value),
            d_lateral=float(dlat.value),
            alpha=float(alpha.value),
        )


class HbondLifetimeAccumulator:
    """Streaming H-bond lifetime accumulator (intermittent C_HB and continuous S_HB)."""

    def __init__(
        self,
        *,
        max_lag: int = 0,
        origin_stride: int = 1,
        dt: float = 1.0,
        directed: bool = True,
    ):
        lib = _capi.get_lib()
        raw = ctypes.c_void_p()
        st = lib.seams_hbond_tcf_create(
            int(max_lag),
            int(origin_stride),
            float(dt),
            1 if directed else 0,
            ctypes.byref(raw),
        )
        self._handle: ctypes.c_void_p | None = None
        _capi.check_status(st)
        self._handle = raw

    def close(self) -> None:
        if self._handle is not None and self._handle.value:
            _capi.get_lib().seams_hbond_tcf_destroy(self._handle)
            self._handle = None

    def __del__(self):
        try:
            self.close()
        except Exception:
            pass

    def __enter__(self) -> HbondLifetimeAccumulator:
        return self

    def __exit__(self, exc_type, exc, tb) -> None:
        self.close()

    def push(self, edges, acceptors=None) -> None:
        if self._handle is None:
            raise RuntimeError("HbondLifetimeAccumulator is closed")
        if acceptors is not None:
            d_arr = np.ascontiguousarray(edges, dtype=np.int32).reshape(-1)
            a_arr = np.ascontiguousarray(acceptors, dtype=np.int32).reshape(-1)
            if d_arr.size != a_arr.size:
                raise ValueError("donors and acceptors must have the same length")
        elif isinstance(edges, np.ndarray):
            if edges.size == 0:
                d_arr = np.empty(0, dtype=np.int32)
                a_arr = np.empty(0, dtype=np.int32)
            else:
                arr = np.ascontiguousarray(edges, dtype=np.int32)
                if arr.ndim != 2 or arr.shape[1] != 2:
                    raise ValueError("edge array must have shape (n_edges, 2)")
                d_arr = np.ascontiguousarray(arr[:, 0], dtype=np.int32)
                a_arr = np.ascontiguousarray(arr[:, 1], dtype=np.int32)
        else:
            pairs = list(edges)
            if not pairs:
                d_arr = np.empty(0, dtype=np.int32)
                a_arr = np.empty(0, dtype=np.int32)
            else:
                arr = np.asarray(pairs, dtype=np.int32)
                if arr.ndim != 2 or arr.shape[1] != 2:
                    raise ValueError("edges must be an iterable of (donor, acceptor) pairs")
                d_arr = np.ascontiguousarray(arr[:, 0], dtype=np.int32)
                a_arr = np.ascontiguousarray(arr[:, 1], dtype=np.int32)

        lib = _capi.get_lib()
        n_edges = int(d_arr.size)
        st = lib.seams_hbond_tcf_push(
            self._handle,
            None if n_edges == 0 else _capi.int_ptr(d_arr),
            None if n_edges == 0 else _capi.int_ptr(a_arr),
            n_edges,
        )
        _capi.check_status(st, self._handle)

    def curve(self) -> HbondLifetimeResult:
        if self._handle is None:
            raise RuntimeError("HbondLifetimeAccumulator is closed")
        lib = _capi.get_lib()
        nl_c = ctypes.c_size_t(0)
        st = lib.seams_hbond_tcf_curve(
            self._handle, None, None, None, None, 0, ctypes.byref(nl_c)
        )
        _capi.check_status(st, self._handle)
        nl = int(nl_c.value)
        lt = np.empty(nl, dtype=np.float64)
        chb = np.empty(nl, dtype=np.float64)
        shb = np.empty(nl, dtype=np.float64)
        orig = np.empty(nl, dtype=np.uintp)
        st = lib.seams_hbond_tcf_curve(
            self._handle,
            _capi.double_ptr(lt),
            _capi.double_ptr(chb),
            _capi.double_ptr(shb),
            _capi.sizet_ptr(orig),
            nl,
            ctypes.byref(nl_c),
        )
        _capi.check_status(st, self._handle)
        return HbondLifetimeResult(
            lag_time=lt,
            c_intermittent=chb,
            s_continuous=shb,
            n_origins=orig.astype(np.int64),
        )

    def fit(self, t_start: float = 0.0, t_end: float = 0.0) -> HbondLifetimeFit:
        if self._handle is None:
            raise RuntimeError("HbondLifetimeAccumulator is closed")
        lib = _capi.get_lib()
        tci = ctypes.c_double(0.0)
        tce = ctypes.c_double(0.0)
        tsi = ctypes.c_double(0.0)
        tse = ctypes.c_double(0.0)
        st = lib.seams_hbond_tcf_fit(
            self._handle,
            float(t_start),
            float(t_end),
            ctypes.byref(tci),
            ctypes.byref(tce),
            ctypes.byref(tsi),
            ctypes.byref(tse),
        )
        _capi.check_status(st, self._handle)
        return HbondLifetimeFit(
            tau_c_int=float(tci.value),
            tau_c_exp=float(tce.value),
            tau_s_int=float(tsi.value),
            tau_s_exp=float(tse.value),
        )


class ReorientationAccumulator:
    """Streaming orientational correlation accumulator for C1(t) and C2(t)."""

    def __init__(
        self,
        n_vectors: int,
        *,
        max_lag: int = 0,
        origin_stride: int = 1,
        dt: float = 1.0,
    ):
        lib = _capi.get_lib()
        raw = ctypes.c_void_p()
        st = lib.seams_reorient_create(
            int(n_vectors),
            int(max_lag),
            int(origin_stride),
            float(dt),
            ctypes.byref(raw),
        )
        self._handle: ctypes.c_void_p | None = None
        _capi.check_status(st)
        self._handle = raw
        self.n_vectors = int(n_vectors)

    def close(self) -> None:
        if self._handle is not None and self._handle.value:
            _capi.get_lib().seams_reorient_destroy(self._handle)
            self._handle = None

    def __del__(self):
        try:
            self.close()
        except Exception:
            pass

    def __enter__(self) -> ReorientationAccumulator:
        return self

    def __exit__(self, exc_type, exc, tb) -> None:
        self.close()

    def push(self, vectors) -> None:
        if self._handle is None:
            raise RuntimeError("ReorientationAccumulator is closed")
        vecs = np.ascontiguousarray(vectors, dtype=np.float64)
        if vecs.shape != (self.n_vectors, 3):
            raise ValueError(
                f"vectors must have shape ({self.n_vectors}, 3); got {vecs.shape}"
            )
        st = _capi.get_lib().seams_reorient_push(
            self._handle, _capi.double_ptr(vecs), self.n_vectors
        )
        _capi.check_status(st, self._handle)

    def curve(self) -> ReorientationResult:
        if self._handle is None:
            raise RuntimeError("ReorientationAccumulator is closed")
        lib = _capi.get_lib()
        nl_c = ctypes.c_size_t(0)
        st = lib.seams_reorient_curve(
            self._handle, None, None, None, None, 0, ctypes.byref(nl_c)
        )
        _capi.check_status(st, self._handle)
        nl = int(nl_c.value)
        lt = np.empty(nl, dtype=np.float64)
        c1 = np.empty(nl, dtype=np.float64)
        c2 = np.empty(nl, dtype=np.float64)
        orig = np.empty(nl, dtype=np.uintp)
        st = lib.seams_reorient_curve(
            self._handle,
            _capi.double_ptr(lt),
            _capi.double_ptr(c1),
            _capi.double_ptr(c2),
            _capi.sizet_ptr(orig),
            nl,
            ctypes.byref(nl_c),
        )
        _capi.check_status(st, self._handle)
        return ReorientationResult(
            lag_time=lt,
            c1=c1,
            c2=c2,
            n_origins=orig.astype(np.int64),
        )

    def fit(self, t_start: float = 0.0, t_end: float = 0.0) -> ReorientationFit:
        if self._handle is None:
            raise RuntimeError("ReorientationAccumulator is closed")
        lib = _capi.get_lib()
        t1i = ctypes.c_double(0.0)
        t1e = ctypes.c_double(0.0)
        t2i = ctypes.c_double(0.0)
        t2e = ctypes.c_double(0.0)
        st = lib.seams_reorient_fit(
            self._handle,
            float(t_start),
            float(t_end),
            ctypes.byref(t1i),
            ctypes.byref(t1e),
            ctypes.byref(t2i),
            ctypes.byref(t2e),
        )
        _capi.check_status(st, self._handle)
        return ReorientationFit(
            tau1_int=float(t1i.value),
            tau1_exp=float(t1e.value),
            tau2_int=float(t2i.value),
            tau2_exp=float(t2e.value),
        )


def hbond_lifetime(
    edge_sets: Iterable,
    *,
    max_lag: int = 0,
    origin_stride: int = 1,
    dt: float = 1.0,
    directed: bool = True,
    fit_window: tuple[float, float] = (0.0, 0.0),
) -> HbondLifetimeResult:
    """Compute intermittent C_HB(t) and continuous S_HB(t) from per-frame edge sets."""
    with HbondLifetimeAccumulator(
        max_lag=max_lag, origin_stride=origin_stride, dt=dt, directed=directed
    ) as acc:
        for edges in edge_sets:
            acc.push(edges)
        res = acc.curve()
        fit_res = (
            acc.fit(fit_window[0], fit_window[1]) if res.lag_time.size >= 2 else None
        )
    return HbondLifetimeResult(
        lag_time=res.lag_time,
        c_intermittent=res.c_intermittent,
        s_continuous=res.s_continuous,
        n_origins=res.n_origins,
        lifetimes=fit_res,
    )


def reorientation_tcf(
    vector_frames: Iterable[np.ndarray],
    *,
    max_lag: int = 0,
    origin_stride: int = 1,
    dt: float = 1.0,
    fit_window: tuple[float, float] = (0.0, 0.0),
) -> ReorientationResult:
    """Stream per-frame (N, 3) molecular vectors into C1(t), C2(t) and tau_1, tau_2."""
    iterator = iter(vector_frames)
    first = next(iterator, None)
    if first is None:
        raise ValueError("vector_frames is empty")
    first_arr = np.ascontiguousarray(first, dtype=np.float64)
    if first_arr.ndim != 2 or first_arr.shape[1] != 3 or first_arr.shape[0] == 0:
        raise ValueError("each vector frame must have shape (N, 3) with N >= 1")
    with ReorientationAccumulator(
        first_arr.shape[0], max_lag=max_lag, origin_stride=origin_stride, dt=dt
    ) as acc:
        acc.push(first_arr)
        for vecs in iterator:
            acc.push(vecs)
        res = acc.curve()
        fit_res = (
            acc.fit(fit_window[0], fit_window[1]) if res.lag_time.size >= 2 else None
        )
    return ReorientationResult(
        lag_time=res.lag_time,
        c1=res.c1,
        c2=res.c2,
        n_origins=res.n_origins,
        relaxation=fit_res,
    )


def _mic_vectors_np(delta: np.ndarray, box9: np.ndarray) -> np.ndarray:
    """Vectorized restricted-triclinic minimum image for shape (..., 3)."""
    lx, ly, lz, xy, xz, yz = map(float, box9[3:9])
    if abs(xy) == 0.0 and abs(xz) == 0.0 and abs(yz) == 0.0:
        lengths = np.array([lx, ly, lz], dtype=np.float64)
        return delta - lengths * np.round(delta / lengths)
    dz = delta[..., 2]
    dy = delta[..., 1]
    dx = delta[..., 0]
    sz = dz / lz
    sy = (dy - yz * sz) / ly
    sx = (dx - xy * sy - xz * sz) / lx
    dsx = sx - np.round(sx)
    dsy = sy - np.round(sy)
    dsz = sz - np.round(sz)
    h = np.array([[lx, xy, xz], [0.0, ly, yz], [0.0, 0.0, lz]], dtype=np.float64)
    ds = np.stack([dsx, dsy, dsz], axis=-1)
    base = ds @ h.T
    out = base.copy()
    if out.size == 0:
        return out
    # Rounding each fractional coordinate is not the Euclidean minimum image
    # on a tilted cell. The minimum image is no longer than the rounded one,
    # and |t_i| <= |v| * |row i of H^-1| bounds the lattice shifts to try.
    best2 = np.einsum("...i,...i->...", out, out)
    rmax = float(np.sqrt(best2.max())) * (1.0 + 1e-12)
    spans = rmax * np.linalg.norm(np.linalg.inv(h), axis=1)
    ranges = [np.arange(int(np.ceil(-0.5 - s)), int(np.floor(0.5 + s)) + 1) for s in spans]
    for ix in ranges[0]:
        for iy in ranges[1]:
            for iz in ranges[2]:
                if ix == 0 and iy == 0 and iz == 0:
                    continue
                cand = base - h @ np.array([ix, iy, iz], dtype=np.float64)
                c2 = np.einsum("...i,...i->...", cand, cand)
                better = c2 < best2
                if np.any(better):
                    out[better] = cand[better]
                    best2[better] = c2[better]
    return out


def water_hbond_edges(
    frame,
    *,
    oxygen_type: int | None = None,
    hydrogen_type: int | None = None,
    oo_cutoff: float = 3.5,
    oh_cutoff: float = 2.42,
    angle_cutoff: float = 30.0,
) -> set[tuple[int, int]]:
    """Directed (donor_mol_id, acceptor_mol_id) H-bond edges matching bond::donatedHydrogenBond."""
    xyz, box9, types, _ids, mols = _extract_frame_data(frame)
    ot, ht = _infer_water_types(types, oxygen_type, hydrogen_type)
    o_mask = types == ot
    h_mask = types == ht
    o_xyz = xyz[o_mask]
    o_mols = mols[o_mask]
    h_xyz = xyz[h_mask]
    h_mols = mols[h_mask]
    n_o = int(o_xyz.shape[0])
    if n_o < 2 or h_xyz.shape[0] == 0:
        return set()

    mol_to_h: dict[int, list[int]] = {}
    for hi, mid in enumerate(h_mols.tolist()):
        mol_to_h.setdefault(mid, []).append(hi)

    cos_cut = float(np.cos(np.deg2rad(angle_cutoff)))
    oo_cut2 = float(oo_cutoff * oo_cutoff)
    oh_cut2 = float(oh_cutoff * oh_cutoff)
    edges: set[tuple[int, int]] = set()

    for d_idx in range(n_o):
        d_mol = int(o_mols[d_idx])
        h_indices = mol_to_h.get(d_mol)
        if not h_indices:
            continue
        oo_vec = _mic_vectors_np(o_xyz - o_xyz[d_idx], box9)
        oo_r2 = np.einsum("ij,ij->i", oo_vec, oo_vec)
        cand = np.flatnonzero((oo_r2 > 1e-12) & (oo_r2 <= oo_cut2))
        if cand.size == 0:
            continue
        oo_c = oo_vec[cand]
        oo_norm = np.sqrt(oo_r2[cand])
        for hi in h_indices:
            oh_c = _mic_vectors_np(o_xyz[cand] - h_xyz[hi], box9)
            oh_r2 = np.einsum("ij,ij->i", oh_c, oh_c)
            within = oh_r2 < oh_cut2
            if not np.any(within):
                continue
            sub_idx = np.flatnonzero(within)
            oh_sub = oh_c[sub_idx]
            oo_sub = oo_c[sub_idx]
            denom = oo_norm[sub_idx] * np.sqrt(oh_r2[sub_idx])
            valid = denom > 1e-15
            cos_ang = np.zeros_like(denom)
            cos_ang[valid] = np.einsum("ij,ij->i", oo_sub[valid], oh_sub[valid]) / denom[valid]
            for k in np.flatnonzero(cos_ang >= cos_cut).tolist():
                a_idx = int(cand[sub_idx[k]])
                edges.add((d_mol, int(o_mols[a_idx])))
    return edges


class TrajectorySequence(Sequence):
    """Multi-frame container with lazy/streaming I/O, slicing, and dynamics.

    Supports ``len(traj)``, ``for frame in traj``, integer indexing
    ``traj[i]`` (including negative indices), and lazy slicing
    ``traj[start:stop:step]``. Open a file with
    :func:`pydseams.open_trajectory` (or :meth:`TrajectorySequence.open`).

    ``pydseams.Trajectory`` is not this class: throughout 2.14.x it stays an
    alias of the single-frame :class:`~pydseams.frame.Frame`.
    """

    def __init__(
        self,
        source=None,
        *,
        atom_type=None,
        cutoff=None,
        bonded="auto",
        all_atoms=False,
        dt: float = 1.0,
        frames=None,
        _indices: range | tuple[int, ...] | None = None,
        _array_backing: tuple[np.ndarray, np.ndarray, np.ndarray | None, np.ndarray | None] | None = None,
    ):
        if not (float(dt) > 0.0) or not np.isfinite(dt):
            raise ValueError("dt must be finite and positive")
        self.dt = float(dt)
        self._atom_type = atom_type
        self._all_atoms = bool(all_atoms)
        self._cutoff = cutoff
        self._bonded = bonded
        self._path: Path | None = None
        self._frames: tuple | None = None
        self._array_backing = _array_backing

        if _array_backing is not None:
            xyz_seq, _boxes, _nums, _mols = _array_backing
            total = int(xyz_seq.shape[0])
            self._indices = range(total) if _indices is None else _indices
            return

        if frames is not None:
            materialized = tuple(frames)
            self._frames = materialized
            self._indices = range(len(materialized)) if _indices is None else _indices
            return

        if source is None:
            raise TypeError(
                "TrajectorySequence requires a path, frames=..., or TrajectorySequence.from_arrays(...)"
            )
        if isinstance(source, (str, Path)):
            p = Path(source).resolve()
            self._path = p
            if _indices is None:
                total = _count_path_frames(p)
                self._indices = range(1, total + 1)
            else:
                self._indices = _indices
            return

        if isinstance(source, Iterable):
            materialized = tuple(source)
            self._frames = materialized
            self._indices = range(len(materialized)) if _indices is None else _indices
            return

        raise TypeError(f"unsupported TrajectorySequence source: {type(source)!r}")

    @classmethod
    def open(
        cls,
        path,
        *,
        atom_type: int | None = None,
        all_atoms: bool = False,
        cutoff: float | None = None,
        bonded: str = "auto",
        dt: float = 1.0,
    ) -> TrajectorySequence:
        """Open a multi-frame trajectory file lazily without reading all frames."""
        return cls(
            source=path,
            atom_type=atom_type,
            all_atoms=all_atoms,
            cutoff=cutoff,
            bonded=bonded,
            dt=dt,
        )

    from_file = open

    @classmethod
    def from_frames(cls, frames: Iterable, *, dt: float = 1.0) -> TrajectorySequence:
        """Construct a TrajectorySequence from an iterable of Frame objects."""
        return cls(frames=frames, dt=dt)

    @classmethod
    def from_arrays(
        cls,
        positions,
        cell,
        *,
        numbers=None,
        mol_ids=None,
        dt: float = 1.0,
        cutoff: float | None = None,
        bonded: str = "cutoff",
    ) -> TrajectorySequence:
        """Build a multi-frame TrajectorySequence from (T, N, 3) positions and cell(s)."""
        arr = np.asarray(positions, dtype=np.float64)
        if arr.ndim != 3 or arr.shape[2] != 3 or arr.shape[0] == 0 or arr.shape[1] == 0:
            raise ValueError("positions must have shape (n_frames, n_atoms, 3)")
        xyz_seq = np.ascontiguousarray(arr, dtype=np.float64)
        nframes, natoms, _ = xyz_seq.shape
        boxes = _normalize_boxes_sequence(cell, nframes)
        nums = (
            None
            if numbers is None
            else np.ascontiguousarray(numbers, dtype=np.int32).reshape(natoms)
        )
        mols = (
            None
            if mol_ids is None
            else np.ascontiguousarray(mol_ids, dtype=np.int32).reshape(natoms)
        )
        return cls(
            cutoff=cutoff,
            bonded=bonded,
            dt=dt,
            _array_backing=(xyz_seq, boxes, nums, mols),
        )

    @property
    def path(self) -> Path | None:
        return self._path

    def __len__(self) -> int:
        return len(self._indices)

    def _frame_at_pos(self, pos: int):
        idx = self._indices[pos]
        if self._array_backing is not None:
            xyz_seq, boxes, nums, mols = self._array_backing
            b9 = boxes[0 if boxes.shape[0] == 1 else idx]
            n = xyz_seq.shape[1]
            t_arr = np.ones(n, dtype=np.int32) if nums is None else nums
            id_arr = np.arange(1, n + 1, dtype=np.int32)
            m_arr = id_arr if mols is None else mols
            return _make_frame_from_raw(
                xyz_seq[idx],
                b9,
                t_arr,
                id_arr,
                m_arr,
                atom_type=self._atom_type,
                cutoff=self._cutoff,
                bonded=self._bonded,
            )
        if self._frames is not None:
            return self._frames[idx]
        assert self._path is not None
        return _load_frame_at_step(
            self._path,
            int(idx),
            atom_type=self._atom_type,
            all_atoms=self._all_atoms,
            cutoff=self._cutoff,
            bonded=self._bonded,
        )

    def _raw_stream(self) -> Iterator[tuple[np.ndarray, np.ndarray, np.ndarray, np.ndarray, np.ndarray]]:
        """Stream (xyz, box9, types, ids, mols) per frame with O(1) frame memory."""
        if self._array_backing is not None:
            xyz_seq, boxes, nums, mols = self._array_backing
            n = xyz_seq.shape[1]
            t_arr = np.ones(n, dtype=np.int32) if nums is None else nums
            id_arr = np.arange(1, n + 1, dtype=np.int32)
            m_arr = id_arr if mols is None else mols
            for idx in self._indices:
                b9 = boxes[0 if boxes.shape[0] == 1 else idx]
                yield xyz_seq[idx], b9, t_arr, id_arr, m_arr
            return

        if self._frames is not None:
            for idx in self._indices:
                yield _extract_frame_data(self._frames[idx])
            return

        assert self._path is not None
        type_filter = (
            -1
            if (self._all_atoms or self._atom_type is None)
            else int(self._atom_type)
        )
        lib = _capi.get_lib()
        sys_h = lib.seams_system_create()
        if not sys_h:
            raise MemoryError("seams_system_create failed")
        encoded = str(self._path).encode("utf-8")
        try:
            for step_1based in self._indices:
                st = lib.seams_system_read(sys_h, encoded, int(step_1based), type_filter)
                _capi.check_status(st, sys_h)
                yield _extract_raw_from_system(sys_h)
        finally:
            lib.seams_system_destroy(sys_h)

    def __iter__(self) -> Iterator:
        for pos in range(len(self._indices)):
            yield self._frame_at_pos(pos)

    def __getitem__(self, item):
        if isinstance(item, slice):
            sub_indices = self._indices[item]
            step = item.step if item.step is not None else 1
            sub_dt = self.dt * abs(step)
            if self._array_backing is not None:
                return TrajectorySequence(
                    atom_type=self._atom_type,
                    all_atoms=self._all_atoms,
                    cutoff=self._cutoff,
                    bonded=self._bonded,
                    dt=sub_dt,
                    _indices=sub_indices,
                    _array_backing=self._array_backing,
                )
            if self._frames is not None:
                return TrajectorySequence(
                    frames=self._frames,
                    dt=sub_dt,
                    _indices=sub_indices,
                )
            return TrajectorySequence(
                source=self._path,
                atom_type=self._atom_type,
                all_atoms=self._all_atoms,
                cutoff=self._cutoff,
                bonded=self._bonded,
                dt=sub_dt,
                _indices=sub_indices,
            )

        if isinstance(item, (int, np.integer)):
            n = len(self._indices)
            idx = int(item)
            if idx < 0:
                idx += n
            if idx < 0 or idx >= n:
                raise IndexError(f"TrajectorySequence index {item} out of range for {n} frames")
            return self._frame_at_pos(idx)

        raise TypeError(f"TrajectorySequence indices must be integers or slices, got {type(item).__name__}")

    def __repr__(self) -> str:
        src = str(self._path.name) if self._path is not None else "in-memory"
        return f"TrajectorySequence(source={src!r}, n_frames={len(self)}, dt={self.dt})"

    def unwrap_positions(self) -> np.ndarray:
        """Stream frames and return PBC-unwrapped coordinates of shape (T, N, 3)."""
        frames_xyz: list[np.ndarray] = []
        boxes: list[np.ndarray] = []
        for xyz, box9, _t, _i, _m in self._raw_stream():
            frames_xyz.append(xyz)
            boxes.append(box9)
        if not frames_xyz:
            raise ValueError("TrajectorySequence has no frames")
        stacked = np.stack(frames_xyz, axis=0)
        box_arr = np.stack(boxes, axis=0)
        return unwrap_positions(stacked, box_arr)

    def unwrap_centres(self, *, mol_ids=None, masses=None) -> UnwrappedCentresResult:
        """Stream frames, assemble intact molecules under PBC, and unwrap centres."""
        frames_xyz: list[np.ndarray] = []
        boxes: list[np.ndarray] = []
        first_mols: np.ndarray | None = None
        for xyz, box9, _t, _i, m_arr in self._raw_stream():
            frames_xyz.append(xyz)
            boxes.append(box9)
            if first_mols is None:
                first_mols = m_arr
        if not frames_xyz:
            raise ValueError("TrajectorySequence has no frames")
        stacked = np.stack(frames_xyz, axis=0)
        box_arr = np.stack(boxes, axis=0)
        use_mols = first_mols if mol_ids is None else mol_ids
        return unwrap_centres(stacked, box_arr, mol_ids=use_mols, masses=masses)

    def water_vectors(
        self,
        *,
        oxygen_type: int | None = None,
        hydrogen_type: int | None = None,
    ) -> list[WaterVectors]:
        """Extract WaterVectors for each frame in the trajectory."""
        return [
            water_vectors(fr, oxygen_type=oxygen_type, hydrogen_type=hydrogen_type)
            for fr in self
        ]

    def molecule_dipoles(self, *, mol_ids=None, charges=None) -> list[MoleculeDipoles]:
        """Extract MoleculeDipoles for each frame in the trajectory."""
        return [
            molecule_dipoles(fr, mol_ids=mol_ids, charges=charges)
            for fr in self
        ]

    def orientation_order(
        self,
        axis=(0.0, 0.0, 1.0),
        *,
        kind: str = "bisector",
        oxygen_type: int | None = None,
        hydrogen_type: int | None = None,
        mol_ids=None,
        charges=None,
    ) -> OrientationOrder:
        """Aggregate axis-projected P1 and P2 order parameters across all frames."""
        per_frame_vecs: list[np.ndarray] = []
        for fr in self:
            if kind == "bisector":
                per_frame_vecs.append(
                    water_vectors(fr, oxygen_type=oxygen_type, hydrogen_type=hydrogen_type).bisector
                )
            elif kind == "hh":
                per_frame_vecs.append(
                    water_vectors(fr, oxygen_type=oxygen_type, hydrogen_type=hydrogen_type).hh_unit
                )
            elif kind == "dipole":
                per_frame_vecs.append(
                    molecule_dipoles(fr, mol_ids=mol_ids, charges=charges).unit_vectors
                )
            else:
                raise ValueError("kind must be 'bisector', 'hh', or 'dipole'")
        if not per_frame_vecs:
            raise ValueError("TrajectorySequence has no frames")
        all_vecs = np.concatenate(per_frame_vecs, axis=0)
        return orientation_order(all_vecs, axis=axis)

    def msd(
        self,
        *,
        axis=(0.0, 0.0, 1.0),
        max_lag: int = 0,
        origin_stride: int = 1,
        dt: float | None = None,
        by_molecule: bool = False,
        mol_ids=None,
        masses=None,
        fit_window: tuple[float, float] = (0.0, 0.0),
    ) -> MsdResult:
        """Stream frames through MsdAccumulator with axial/lateral split and fit."""
        step_dt = self.dt if dt is None else float(dt)
        acc: MsdAccumulator | None = None
        try:
            for xyz, box9, types, ids, stored_mols in self._raw_stream():
                if by_molecule:
                    use_mols = stored_mols if mol_ids is None else mol_ids
                    mc = molecule_centres(
                        xyz, box9, mol_ids=use_mols, masses=masses
                    )
                    pts = mc.centres
                else:
                    pts = xyz
                if acc is None:
                    acc = MsdAccumulator(
                        pts.shape[0],
                        max_lag=max_lag,
                        origin_stride=origin_stride,
                        dt=step_dt,
                        axis=axis,
                    )
                acc.push(pts, cell=box9)
            if acc is None:
                raise ValueError("TrajectorySequence has no frames")
            res = acc.curve()
            fit_res = (
                acc.fit(fit_window[0], fit_window[1])
                if res.lag_time.size >= 2
                else None
            )
            return MsdResult(
                lag_time=res.lag_time,
                total=res.total,
                axial=res.axial,
                lateral=res.lateral,
                n_origins=res.n_origins,
                diffusion=fit_res,
            )
        finally:
            if acc is not None:
                acc.close()

    def hbond_lifetime(
        self,
        edge_sets: Iterable | None = None,
        *,
        max_lag: int = 0,
        origin_stride: int = 1,
        dt: float | None = None,
        directed: bool = True,
        oxygen_type: int | None = None,
        hydrogen_type: int | None = None,
        oo_cutoff: float = 3.5,
        oh_cutoff: float = 2.42,
        angle_cutoff: float = 30.0,
        fit_window: tuple[float, float] = (0.0, 0.0),
    ) -> HbondLifetimeResult:
        """Compute continuous S_HB(t) and intermittent C_HB(t) lifetime curves."""
        step_dt = self.dt if dt is None else float(dt)
        if edge_sets is None:
            edge_iter = (
                water_hbond_edges(
                    fr,
                    oxygen_type=oxygen_type,
                    hydrogen_type=hydrogen_type,
                    oo_cutoff=oo_cutoff,
                    oh_cutoff=oh_cutoff,
                    angle_cutoff=angle_cutoff,
                )
                for fr in self
            )
        else:
            edge_iter = edge_sets
        return hbond_lifetime(
            edge_iter,
            max_lag=max_lag,
            origin_stride=origin_stride,
            dt=step_dt,
            directed=directed,
            fit_window=fit_window,
        )

    def reorientation(
        self,
        vectors: Iterable[np.ndarray] | None = None,
        *,
        kind: str = "bisector",
        max_lag: int = 0,
        origin_stride: int = 1,
        dt: float | None = None,
        oxygen_type: int | None = None,
        hydrogen_type: int | None = None,
        mol_ids=None,
        charges=None,
        fit_window: tuple[float, float] = (0.0, 0.0),
    ) -> ReorientationResult:
        """Compute dipole/orientation C1(t) and C2(t) correlations and tau fits."""
        step_dt = self.dt if dt is None else float(dt)
        if vectors is None:
            def _vec_gen():
                for fr in self:
                    if kind == "bisector":
                        yield water_vectors(
                            fr, oxygen_type=oxygen_type, hydrogen_type=hydrogen_type
                        ).bisector
                    elif kind == "hh":
                        yield water_vectors(
                            fr, oxygen_type=oxygen_type, hydrogen_type=hydrogen_type
                        ).hh_unit
                    elif kind == "dipole":
                        yield molecule_dipoles(
                            fr, mol_ids=mol_ids, charges=charges
                        ).unit_vectors
                    else:
                        raise ValueError("kind must be 'bisector', 'hh', or 'dipole'")

            vec_iter = _vec_gen()
        else:
            vec_iter = vectors
        return reorientation_tcf(
            vec_iter,
            max_lag=max_lag,
            origin_stride=origin_stride,
            dt=step_dt,
            fit_window=fit_window,
        )
