# Copyright (c) d-SEAMS developers
"""Analytic-oracle and negative-control tests for pydseams TrajectorySequence and dynamics."""

from __future__ import annotations

from pathlib import Path

import numpy as np
import pytest

import pydseams as ds


def _wrap_triclinic(xyz: np.ndarray, box9: np.ndarray) -> np.ndarray:
    """Wrap (..., 3) Cartesian coordinates into primary restricted-triclinic cell."""
    ox, oy, oz, lx, ly, lz, xy, xz, yz = map(float, box9)
    sz = (xyz[..., 2] - oz) / lz
    sy = (xyz[..., 1] - oy - yz * sz) / ly
    sx = (xyz[..., 0] - ox - xy * sy - xz * sz) / lx
    sx -= np.floor(sx)
    sy -= np.floor(sy)
    sz -= np.floor(sz)
    out = np.empty_like(xyz)
    out[..., 0] = ox + lx * sx + xy * sy + xz * sz
    out[..., 1] = oy + ly * sy + yz * sz
    out[..., 2] = oz + lz * sz
    return out


def _write_multiframe_lammps_dump(
    path: Path,
    frames_xyz: np.ndarray,
    box_len: float,
    types: np.ndarray,
    mols: np.ndarray,
    orders: list[np.ndarray] | None = None,
    ids: np.ndarray | None = None,
) -> None:
    """Atom i has ID ids[i] (default i + 1); frame f lists atoms in orders[f].

    mols and ids may also be (n_frames, n_atoms), one row per frame.
    """
    lines: list[str] = []
    n_frames, n_atoms, _ = frames_xyz.shape
    atom_ids = np.arange(1, n_atoms + 1) if ids is None else np.asarray(ids)
    for f in range(n_frames):
        frame_mols = mols[f] if np.ndim(mols) == 2 else mols
        frame_ids = atom_ids[f] if atom_ids.ndim == 2 else atom_ids
        lines.extend(
            [
                "ITEM: TIMESTEP",
                str(f * 100),
                "ITEM: NUMBER OF ATOMS",
                str(n_atoms),
                "ITEM: BOX BOUNDS pp pp pp",
                f"0.0 {box_len:.6f}",
                f"0.0 {box_len:.6f}",
                f"0.0 {box_len:.6f}",
                "ITEM: ATOMS id mol type x y z",
            ]
        )
        for i in (range(n_atoms) if orders is None else orders[f]):
            x, y, z = frames_xyz[f, i]
            lines.append(
                f"{int(frame_ids[i])} {int(frame_mols[i])} {int(types[i])} {x:.8f} {y:.8f} {z:.8f}"
            )
    path.write_text("\n".join(lines) + "\n", encoding="utf-8")


def _write_multiframe_pdb(path: Path, frames_xyz: np.ndarray, box_len: float) -> None:
    lines: list[str] = []
    n_frames, n_atoms, _ = frames_xyz.shape
    for f in range(n_frames):
        lines.append(f"MODEL     {f + 1:4d}")
        lines.append(
            f"CRYST1{box_len:9.3f}{box_len:9.3f}{box_len:9.3f}  90.00  90.00  90.00 P 1           1"
        )
        for i in range(n_atoms):
            x, y, z = frames_xyz[f, i]
            elem = "O" if (i % 3 == 0) else "H"
            name = "OW" if elem == "O" else ("HW1" if i % 3 == 1 else "HW2")
            mol_id = (i // 3) + 1
            lines.append(
                f"ATOM  {i + 1:5d} {name:<4s} HOH A{mol_id:4d}    {x:8.3f}{y:8.3f}{z:8.3f}  1.00  0.00          {elem:>2s}"
            )
        lines.append("ENDMDL")
    path.write_text("\n".join(lines) + "\n", encoding="utf-8")


def test_trajectory_sequence_container(tmp_path: Path):
    assert ds.TrajectorySequence is not ds.Frame
    assert not issubclass(ds.TrajectorySequence, ds.Frame)

    rng = np.random.default_rng(42)
    n_frames = 6
    n_waters = 4
    types = np.tile([2, 1, 1], n_waters)
    mols = np.repeat(np.arange(1, n_waters + 1), 3)
    xyz = np.zeros((n_frames, 3 * n_waters, 3), dtype=np.float64)
    for f in range(n_frames):
        for w in range(n_waters):
            base = np.array([2.0 + 3.0 * (w % 2), 2.0 + 3.0 * (w // 2), 1.5 + 0.4 * f]) + 0.02 * rng.normal(size=3)
            xyz[f, 3 * w] = base
            xyz[f, 3 * w + 1] = base + [0.757, 0.0, 0.586]
            xyz[f, 3 * w + 2] = base + [-0.757, 0.0, 0.586]

    dump_path = tmp_path / "water6.lammpstrj"
    _write_multiframe_lammps_dump(dump_path, xyz, 12.0, types, mols)

    traj = ds.open_trajectory(dump_path, all_atoms=True, dt=0.2)
    assert len(traj) == 6
    assert traj[0].n_atoms == 12
    assert traj[-1].n_atoms == 12
    with pytest.raises(IndexError):
        _ = traj[6]
    with pytest.raises(IndexError):
        _ = traj[-7]

    # Slicing is lazy and scales dt by step.
    sub = traj[1:6:2]
    assert isinstance(sub, ds.TrajectorySequence)
    assert len(sub) == 3
    assert sub.dt == pytest.approx(0.4)
    streamed = [fr for fr in sub]
    assert len(streamed) == 3
    assert np.allclose(ds.water_vectors(streamed[0]).mol_ids, [1, 2, 3, 4])

    # Chemfiles multi-frame PDB streaming check.
    pdb_path = tmp_path / "water6.pdb"
    _write_multiframe_pdb(pdb_path, xyz, 12.0)
    pdb_traj = ds.open_trajectory(pdb_path, all_atoms=True, dt=0.5)
    assert len(pdb_traj) == 6
    assert len(pdb_traj[::3]) == 2
    assert pdb_traj[-1].n_atoms == 12


_EXAMPLE_TRAJ = Path(__file__).resolve().parents[1] / "input" / "traj" / "exampleTraj.lammpstrj"


def test_trajectory_stays_the_frame_alias():
    assert ds.Trajectory is ds.Frame
    assert {"Frame", "Trajectory", "TrajectorySequence", "open_trajectory"} <= set(ds.__all__)
    traj = ds.open_trajectory(_EXAMPLE_TRAJ, atom_type=2)
    assert type(traj) is ds.TrajectorySequence
    assert not isinstance(traj, ds.Trajectory)
    assert isinstance(traj[0], ds.Trajectory)


def test_trajectory_sequence_api_does_not_warn():
    import warnings

    with warnings.catch_warnings():
        warnings.simplefilter("error")
        traj = ds.open_trajectory(_EXAMPLE_TRAJ, atom_type=2)
        assert len(traj) == 1
        assert traj[0].n_atoms == 250
        assert isinstance(traj[:1], ds.TrajectorySequence)
        assert isinstance(ds.TrajectorySequence(_EXAMPLE_TRAJ), ds.TrajectorySequence)
        assert isinstance(ds.TrajectorySequence.open(_EXAMPLE_TRAJ), ds.TrajectorySequence)
        assert isinstance(ds.TrajectorySequence.from_file(_EXAMPLE_TRAJ), ds.TrajectorySequence)
        assert isinstance(ds.TrajectorySequence.from_frames([traj[0]]), ds.TrajectorySequence)


def test_trajectory_sequence_has_no_frame_constructor_or_attributes():
    import warnings

    with pytest.raises(TypeError):
        ds.TrajectorySequence(_EXAMPLE_TRAJ, 1)
    with pytest.raises(TypeError):
        ds.TrajectorySequence(_EXAMPLE_TRAJ, frame=1)
    with pytest.raises(ValueError, match="n_frames, n_atoms, 3"):
        ds.TrajectorySequence.from_arrays(np.zeros((4, 3)), [10.0, 10.0, 10.0])
    traj = ds.open_trajectory(_EXAMPLE_TRAJ)
    with warnings.catch_warnings():
        warnings.simplefilter("error")
        for name in ("n_atoms", "classify_chill_plus", "filename", "load_frame", "positions"):
            with pytest.raises(AttributeError):
                getattr(traj, name)


_HAS_PYDSEAMSLIB_FRAME = hasattr(ds.Frame, "classify_chill_plus")


@pytest.mark.skipif(not _HAS_PYDSEAMSLIB_FRAME, reason="pydseamslib Frame is not installed")
def test_pydseamslib_trajectory_alias_workflow():
    import warnings

    # PydSEAMSlib tests/python/test_user_api.py::test_trajectory_alias, run
    # through this adapter, plus the other Frame attributes the alias carries.
    with warnings.catch_warnings():
        warnings.simplefilter("error", DeprecationWarning)
        from pydseams import Frame, Trajectory

        assert Trajectory is Frame
        traj = Trajectory(str(_EXAMPLE_TRAJ.resolve()))
        counts = traj.classify_chill_plus()
        assert counts["hexagonal"] == 0
        assert counts["water"] == 241
        assert counts["interClathrate"] == 9
        assert traj.filename == str(_EXAMPLE_TRAJ.resolve())
        traj.load_frame(1)
        assert traj.frame == 1
        typed = Trajectory(str(_EXAMPLE_TRAJ), atom_type=2)
        assert typed.atom_type == 2
        assert typed.classify_chill_plus()["water"] == 241


def test_triclinic_and_time_varying_box_unwrapping_oracle():
    n_frames = 40
    n_mols = 5
    mols = np.repeat(np.arange(1, n_mols + 1, dtype=np.int32), 3)
    masses = np.tile([15.9994, 1.00794, 1.00794], n_mols)
    m_tot = 15.9994 + 2.0 * 1.00794

    rel = np.array(
        [[0.0, 0.0, 0.0], [0.75695, 0.0, 0.58588], [-0.75695, 0.0, 0.58588]],
        dtype=np.float64,
    )
    com_offset = (masses[:3, None] * rel).sum(axis=0) / m_tot

    rng = np.random.default_rng(2026)
    boxes = np.zeros((n_frames, 9), dtype=np.float64)
    for f in range(n_frames):
        # Slowly breathing triclinic box (NPT-like) where step << halfMin.
        scale = 1.0 + 0.01 * np.sin(0.2 * f)
        boxes[f] = [0.0, 0.0, 0.0, 11.0 * scale, 12.0 * scale, 13.0 * scale, 2.5, -1.5, 2.0]

    s0 = rng.uniform(0.15, 0.85, size=(n_mols, 3))
    r0 = np.column_stack(
        [
            boxes[0, 3] * s0[:, 0] + boxes[0, 6] * s0[:, 1] + boxes[0, 7] * s0[:, 2],
            boxes[0, 4] * s0[:, 1] + boxes[0, 8] * s0[:, 2],
            boxes[0, 5] * s0[:, 2],
        ]
    )
    vel = np.array(
        [
            [1.6, -1.4, 1.8],
            [-1.9, 1.7, -1.5],
            [1.2, 2.0, -1.3],
            [-1.1, -1.6, 2.1],
            [2.2, 0.8, 1.5],
        ],
        dtype=np.float64,
    )

    # Construct true trajectory by integrating displacements and wrapping each frame.
    wrapped_atoms = np.zeros((n_frames, 3 * n_mols, 3), dtype=np.float64)
    true_com = np.zeros((n_frames, n_mols, 3), dtype=np.float64)
    for f in range(n_frames):
        for m in range(n_mols):
            o_pos = r0[m] + vel[m] * f
            true_com[f, m] = o_pos + com_offset
            for a in range(3):
                wrapped_atoms[f, 3 * m + a] = _wrap_triclinic(o_pos + rel[a], boxes[f])

    traj = ds.TrajectorySequence.from_arrays(wrapped_atoms, boxes, mol_ids=mols, dt=0.1)
    res = traj.unwrap_centres(masses=masses)
    assert np.array_equal(res.mol_ids, np.arange(1, n_mols + 1))
    assert np.allclose(res.centres, true_com, atol=1e-10)


_SKEWED_CELL = np.array(
    [0.0, 0.0, 0.0, 18.03850526, 10.63623724, 11.40838967, -8.77049545, 6.89505882, 2.90161048]
)
_SKEWED_EUCLID_STEP = np.array([3.42737868, 3.53897334, -3.29810805])


def _rounded_fractional_image(delta: np.ndarray, box9: np.ndarray) -> np.ndarray:
    _ox, _oy, _oz, lx, ly, lz, xy, xz, yz = map(float, box9)
    h = np.array([[lx, xy, xz], [0.0, ly, yz], [0.0, 0.0, lz]])
    s = np.linalg.solve(h, delta)
    return h @ (s - np.round(s))


def test_skewed_triclinic_unwrap_takes_euclidean_image_across_cell_change():
    box0 = _SKEWED_CELL
    lx, ly = box0[3], box0[4]
    assert abs(box0[6]) <= lx / 2 and abs(box0[7]) <= lx / 2 and abs(box0[8]) <= ly / 2
    raw = np.array([-16.98407757, -4.19565342, 8.11028162])
    frac = _rounded_fractional_image(raw, box0)
    assert np.linalg.norm(frac - _SKEWED_EUCLID_STEP) == pytest.approx(lx, abs=1e-7)

    scale = 1.03
    box1 = box0.copy()
    box1[3:] *= scale
    boxes = np.stack([box0, box0, box1])
    u0 = np.array([4.1, 2.7, 5.3])
    u1 = u0 + _SKEWED_EUCLID_STEP
    u2 = u1 + scale * _SKEWED_EUCLID_STEP
    truth = np.stack([u0, u1, u2])[:, None, :]
    wrapped = np.stack([_wrap_triclinic(truth[f], boxes[f]) for f in range(3)])
    raw2 = wrapped[2, 0] - u1
    assert np.linalg.norm(
        _rounded_fractional_image(raw2, box1) - scale * _SKEWED_EUCLID_STEP
    ) == pytest.approx(scale * lx, abs=1e-7)

    out = ds.unwrap_positions(wrapped, boxes)
    assert np.allclose(out, truth, atol=1e-9)


def test_python_minimum_image_matches_brute_force_on_skewed_cells():
    from pydseams.trajectory import _mic_vectors_np

    rng = np.random.default_rng(20261120)
    shifts = np.stack(
        np.meshgrid(*(np.arange(-16, 17),) * 3, indexing="ij"), axis=-1
    ).reshape(-1, 3)
    rounded_misses = 0
    for _ in range(200):
        lx, ly, lz = rng.uniform(1.0, 30.0, 3)
        box9 = np.array(
            [0.0, 0.0, 0.0, lx, ly, lz, rng.uniform(-0.5, 0.5) * lx,
             rng.uniform(-0.5, 0.5) * lx, rng.uniform(-0.5, 0.5) * ly]
        )
        h = np.array([[lx, box9[6], box9[7]], [0.0, ly, box9[8]], [0.0, 0.0, lz]])
        deltas = rng.uniform(-1.5, 1.5, size=(8, 3)) @ h.T
        got = _mic_vectors_np(deltas, box9)
        for d, g in zip(deltas, got):
            cand = d[None, :] - shifts @ h.T
            want = np.sqrt(np.min(np.einsum("ij,ij->i", cand, cand)))
            assert np.linalg.norm(g) == pytest.approx(want, abs=1e-9)
            if np.linalg.norm(_rounded_fractional_image(d, box9)) > want + 1e-9:
                rounded_misses += 1
    assert rounded_misses > 0


def test_ballistic_and_anisotropic_brownian_msd_oracles():
    # 1. Ballistic oracle along arbitrary channel axis e = (1, 1, 2) / sqrt(6)
    axis = np.array([1.0, 1.0, 2.0])
    e_hat = axis / np.linalg.norm(axis)
    e_perp1 = np.array([1.0, -1.0, 0.0]) / np.sqrt(2.0)
    e_perp2 = np.cross(e_hat, e_perp1)

    dt = 0.2
    n_frames = 30
    v_par = 1.5
    v_lat = 0.9
    vel = v_par * e_hat + v_lat * e_perp1
    box = np.array([0.0, 0.0, 0.0, 8.0, 9.0, 10.0, 1.2, -0.8, 0.6])
    t_grid = np.arange(n_frames) * dt
    raw = (vel[None, :] * t_grid[:, None])[:, None, :]
    wrapped = _wrap_triclinic(raw, box)

    traj_bal = ds.TrajectorySequence.from_arrays(wrapped, box, dt=dt)
    msd_bal = traj_bal.msd(axis=axis, fit_window=(0.2, 4.0))
    assert np.allclose(msd_bal.axial, (v_par**2) * (msd_bal.lag_time**2), atol=1e-11)
    assert np.allclose(msd_bal.lateral, (v_lat**2) * (msd_bal.lag_time**2), atol=1e-11)
    assert np.allclose(msd_bal.total, msd_bal.axial + msd_bal.lateral, atol=1e-12)
    assert msd_bal.diffusion is not None
    assert msd_bal.diffusion.alpha == pytest.approx(2.0, abs=1e-10)

    # 2. Anisotropic Brownian random walk with known D_axial = 0.90, D_lateral = 0.20
    rng = np.random.default_rng(12345)
    n_particles = 500
    n_steps = 240
    dt_br = 0.05
    d_ax_true = 0.90
    d_lat_true = 0.20
    d_tot_true = (d_ax_true + 2.0 * d_lat_true) / 3.0

    sig_ax = np.sqrt(2.0 * d_ax_true * dt_br)
    sig_lat = np.sqrt(2.0 * d_lat_true * dt_br)
    steps = (
        rng.normal(0.0, sig_ax, size=(n_steps, n_particles, 1)) * e_hat
        + rng.normal(0.0, sig_lat, size=(n_steps, n_particles, 1)) * e_perp1
        + rng.normal(0.0, sig_lat, size=(n_steps, n_particles, 1)) * e_perp2
    )
    steps[0] = 0.0
    unwrapped = np.cumsum(steps, axis=0)
    wrapped_br = _wrap_triclinic(unwrapped, [0, 0, 0, 12.0, 12.0, 12.0, 2.0, 0.0, 0.0])

    traj_br = ds.TrajectorySequence.from_arrays(
        wrapped_br, [0, 0, 0, 12.0, 12.0, 12.0, 2.0, 0.0, 0.0], dt=dt_br
    )
    msd_br = traj_br.msd(axis=axis, max_lag=80, fit_window=(0.25, 2.0))
    assert msd_br.diffusion is not None
    assert msd_br.diffusion.d_axial == pytest.approx(d_ax_true, rel=0.05)
    assert msd_br.diffusion.d_lateral == pytest.approx(d_lat_true, rel=0.05)
    assert msd_br.diffusion.d_total == pytest.approx(d_tot_true, rel=0.05)
    assert msd_br.diffusion.alpha == pytest.approx(1.0, abs=0.05)

    # 3. Negative control: purely axial motion has zero lateral MSD
    axial_only = (e_hat[None, :] * np.arange(12, dtype=np.float64)[:, None])[:, None, :]
    res_ax = ds.TrajectorySequence.from_arrays(axial_only, [50.0, 50.0, 50.0], dt=1.0).msd(axis=axis)
    assert np.allclose(res_ax.lateral, 0.0, atol=1e-24)
    assert np.allclose(res_ax.axial, res_ax.total, atol=1e-12)


def test_skipped_loglog_fit_reports_nan_alpha():
    # One origin: lag 1 moves 1e-11 (MSD 1e-22 < 1e-20), lag 2 moves 2.
    acc = ds.MsdAccumulator(1, origin_stride=10, dt=1.0, axis=(1.0, 0.0, 0.0))
    for x in (0.0, 1e-11, 2.0):
        acc.push(np.array([[x, 0.0, 0.0]]))
    fit = acc.fit(0.0, 2.0)
    assert np.isnan(fit.alpha)
    assert not fit.alpha_fitted
    assert fit.d_total == pytest.approx(1.0 / 3.0, abs=1e-12)
    assert fit.d_axial == pytest.approx(1.0, abs=1e-12)

    # The pure-Python fitter agrees, including a two-point line.
    py_fit = ds.fit_diffusion([0.0, 1.0], [0.0, 1.0], t_start=0.0, t_end=1.0)
    assert np.isnan(py_fit.alpha)
    assert not py_fit.alpha_fitted
    assert py_fit.d_total == pytest.approx(1.0 / 6.0, abs=1e-14)

    ballistic = ds.fit_diffusion([0.0, 1.0, 2.0, 3.0], [0.0, 1.0, 4.0, 9.0], t_start=1.0, t_end=3.0)
    assert ballistic.alpha_fitted
    assert ballistic.alpha == pytest.approx(2.0, abs=1e-12)


def test_dipole_orientation_and_order_parameter_oracles():
    # Fibonacci sphere of 600 nearly uniform directions -> <P1> ~ 0, <P2> ~ 0
    n_pts = 600
    idx = np.arange(n_pts, dtype=np.float64) + 0.5
    phi = np.pi * (1.0 + np.sqrt(5.0)) * idx
    cos_t = 1.0 - 2.0 * idx / n_pts
    sin_t = np.sqrt(np.maximum(0.0, 1.0 - cos_t**2))
    iso_vecs = np.column_stack([sin_t * np.cos(phi), sin_t * np.sin(phi), cos_t])

    ord_iso = ds.orientation_order(iso_vecs, axis=(0.3, -0.4, 0.866))
    assert abs(ord_iso.mean_p1) < 2e-3
    assert abs(ord_iso.mean_p2) < 2e-3

    # Magic-angle cone cos(theta) = 1/sqrt(3) -> P2 = 0 identically, P1 = 1/sqrt(3)
    c_magic = 1.0 / np.sqrt(3.0)
    s_magic = np.sqrt(2.0 / 3.0)
    angles = np.linspace(0.0, 2.0 * np.pi, 16, endpoint=False)
    cone = np.column_stack([s_magic * np.cos(angles), s_magic * np.sin(angles), np.full(16, c_magic)])
    ord_cone = ds.orientation_order(cone, axis=(0.0, 0.0, 1.0))
    assert ord_cone.mean_p1 == pytest.approx(c_magic, abs=1e-13)
    assert ord_cone.mean_p2 == pytest.approx(0.0, abs=1e-13)

    # Zero-norm centrosymmetric CO2-like molecule is skipped cleanly without NaN pollution
    xyz = np.array(
        [
            [5.0, 5.0, 5.0],  # C
            [5.0, 5.0, 6.16],  # O1
            [5.0, 5.0, 3.84],  # O2
            [1.0, 1.0, 1.0],  # Polar diatomic A
            [1.0, 1.0, 2.0],  # Polar diatomic B
        ],
        dtype=np.float64,
    )
    mols = np.array([1, 1, 1, 2, 2], dtype=np.int32)
    charges = np.array([0.70, -0.35, -0.35, -0.50, 0.50], dtype=np.float64)
    dips = ds.molecule_dipoles(xyz, [10.0, 10.0, 10.0], mol_ids=mols, charges=charges)
    assert dips.magnitudes[0] == pytest.approx(0.0, abs=1e-14)
    assert dips.magnitudes[1] == pytest.approx(0.50, abs=1e-14)
    ord_mix = ds.orientation_order(dips.dipoles, axis=(0.0, 0.0, 1.0))
    assert np.isnan(ord_mix.cos_theta[0])
    assert ord_mix.mean_p1 == pytest.approx(1.0, abs=1e-14)
    assert ord_mix.mean_p2 == pytest.approx(1.0, abs=1e-14)


def test_hbond_lifetime_and_reorientation_oracles():
    # 1. Exact analytical 4-frame H-bond survival & intermittent table
    # e1=(1,2): [1, 1, 1, 1]
    # e2=(3,4): [1, 0, 1, 1]
    # e3=(5,6): [1, 1, 0, 0]
    edge_sets = [
        {(1, 2), (3, 4), (5, 6)},
        {(1, 2), (5, 6)},
        {(1, 2), (3, 4)},
        {(1, 2), (3, 4)},
    ]
    hb = ds.hbond_lifetime(edge_sets, dt=0.5, directed=True)
    # Lag 0: denom = 3+2+2+2 = 9, num = 9 -> 1.0
    # Lag 1: origins 0(3->2: e1,e3=2), 1(2->2: e1=1), 2(2->2: e1,e2=2); denom = 3+2+2 = 7, C=5/7, S=5/7
    # Lag 2: origins 0(e1,e2=2 intermittent; only e1=1 continuous!), 1(e1=1); denom = 3+2 = 5, C=3/5, S=2/5
    # Lag 3: origin 0(e1,e2=2 intermittent; only e1=1 continuous!); denom = 3, C=2/3, S=1/3
    assert np.allclose(hb.c_hb, [1.0, 5.0 / 7.0, 3.0 / 5.0, 2.0 / 3.0], atol=1e-14)
    assert np.allclose(hb.s_hb, [1.0, 5.0 / 7.0, 2.0 / 5.0, 1.0 / 3.0], atol=1e-14)
    assert np.all(hb.s_hb <= hb.c_hb + 1e-15)

    # Directed vs undirected negative control: (1, 2) -> (2, 1)
    assert ds.hbond_lifetime([{(1, 2)}, {(2, 1)}], directed=True).c_hb[1] == pytest.approx(0.0)
    assert ds.hbond_lifetime([{(1, 2)}, {(2, 1)}], directed=False).c_hb[1] == pytest.approx(1.0)

    # 2. Debye rotational diffusion stochastic ensemble: tau1 = 1/(2 Dr), tau2 = 1/(6 Dr), tau1/tau2 = 3
    rng = np.random.default_rng(777)
    n_rotors = 1200
    n_steps = 160
    dt = 0.02
    d_rot = 0.5  # tau1 = 1.0, tau2 = 1/3
    u = rng.normal(size=(n_rotors, 3))
    u /= np.linalg.norm(u, axis=1, keepdims=True)
    frames_u = [u.copy()]
    sigma = np.sqrt(2.0 * d_rot * dt)
    for _ in range(1, n_steps):
        omega_vec = rng.normal(0.0, sigma, size=(n_rotors, 3))
        # Exact Rodrigues rotation on S^2 so |u| = 1 identically and <u(0).u(t)> = exp(-2 Dr t)
        theta = np.linalg.norm(omega_vec, axis=1, keepdims=True)
        k = np.divide(omega_vec, theta, out=np.zeros_like(omega_vec), where=theta > 1e-15)
        cos_th = np.cos(theta)
        sin_th = np.sin(theta)
        k_dot_u = np.sum(k * u, axis=1, keepdims=True)
        u = u * cos_th + np.cross(k, u) * sin_th + k * k_dot_u * (1.0 - cos_th)
        u /= np.linalg.norm(u, axis=1, keepdims=True)
        frames_u.append(u.copy())

    reor = ds.reorientation_tcf(frames_u, max_lag=80, dt=dt, fit_window=(0.02, 0.8))
    assert reor.relaxation is not None
    tau1_true = 1.0 / (2.0 * d_rot)
    tau2_true = 1.0 / (6.0 * d_rot)
    assert reor.relaxation.tau1_exp == pytest.approx(tau1_true, rel=0.06)
    assert reor.relaxation.tau2_exp == pytest.approx(tau2_true, rel=0.06)
    assert (reor.relaxation.tau1_exp / reor.relaxation.tau2_exp) == pytest.approx(3.0, rel=0.08)


def _moving_waters(n_frames: int, box_len: float) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    """Unwrapped (T, 12, 3) positions of 4 waters that drift and turn, plus types and mols."""
    n_waters = 4
    types = np.tile([2, 1, 1], n_waters)
    mols = np.repeat(np.arange(1, n_waters + 1), 3)
    oh1 = np.array([0.757, 0.0, 0.586])
    oh2 = np.array([-0.757, 0.0, 0.586])
    velocity = np.array([[0.9, 0.1, -0.3], [-0.7, 0.5, 0.2], [0.2, -0.8, 0.6], [0.4, 0.4, -0.9]])
    spin = np.array([0.35, -0.2, 0.5, 0.1])
    xyz = np.zeros((n_frames, 3 * n_waters, 3), dtype=np.float64)
    for f in range(n_frames):
        for w in range(n_waters):
            centre = np.array([2.0 + 5.0 * (w % 2), 2.0 + 5.0 * (w // 2), 6.0]) + f * velocity[w]
            a = spin[w] * f
            rot = np.array([[np.cos(a), -np.sin(a), 0.0], [np.sin(a), np.cos(a), 0.0], [0.0, 0.0, 1.0]])
            xyz[f, 3 * w] = centre
            xyz[f, 3 * w + 1] = centre + rot @ oh1
            xyz[f, 3 * w + 2] = centre + rot @ oh2
    return xyz, types, mols


@pytest.mark.parametrize("atom_ids", [None, 1000 + 7 * np.arange(12)], ids=["ids-1-to-n", "ids-with-gaps"])
def test_reordered_dump_lines_keep_atom_and_molecule_identity(tmp_path: Path, atom_ids):
    # A dump may list atoms in any order each frame. Every trajectory result
    # must match the same dump written in one fixed order. IDs other than
    # 1..N take the native reader, which keeps the dump's line order.
    box_len = 12.0
    n_frames = 8
    truth, types, mols = _moving_waters(n_frames, box_len)
    wrapped = truth - box_len * np.floor(truth / box_len)
    rng = np.random.default_rng(2026)
    orders = [np.arange(12)] + [rng.permutation(12) for _ in range(n_frames - 1)]
    assert any(not np.array_equal(o, np.arange(12)) for o in orders)
    fixed_path = tmp_path / "fixed.lammpstrj"
    shuffled_path = tmp_path / "shuffled.lammpstrj"
    _write_multiframe_lammps_dump(fixed_path, wrapped, box_len, types, mols, ids=atom_ids)
    _write_multiframe_lammps_dump(shuffled_path, wrapped, box_len, types, mols, orders=orders, ids=atom_ids)
    fixed = ds.open_trajectory(fixed_path, all_atoms=True, dt=0.5)
    shuffled = ds.open_trajectory(shuffled_path, all_atoms=True, dt=0.5)

    up_fixed = fixed.unwrap_positions()
    up_shuffled = shuffled.unwrap_positions()
    assert np.allclose(up_fixed - up_fixed[0], truth - truth[0], atol=1e-6)
    assert np.allclose(up_shuffled, up_fixed, atol=1e-6)

    uc_fixed = fixed.unwrap_centres()
    uc_shuffled = shuffled.unwrap_centres()
    assert np.array_equal(uc_shuffled.mol_ids, uc_fixed.mol_ids)
    assert np.allclose(uc_shuffled.centres, uc_fixed.centres, atol=1e-6)

    for by_molecule in (False, True):
        m_fixed = fixed.msd(by_molecule=by_molecule)
        m_shuffled = shuffled.msd(by_molecule=by_molecule)
        assert m_fixed.total[1] > 0.1
        assert np.allclose(m_shuffled.total, m_fixed.total, atol=1e-8)
        assert np.allclose(m_shuffled.axial, m_fixed.axial, atol=1e-8)
        assert np.allclose(m_shuffled.lateral, m_fixed.lateral, atol=1e-8)

    for kind in ("bisector", "hh", "dipole"):
        r_fixed = fixed.reorientation(kind=kind, oxygen_type=2, hydrogen_type=1)
        r_shuffled = shuffled.reorientation(kind=kind, oxygen_type=2, hydrogen_type=1)
        assert np.allclose(r_shuffled.c1, r_fixed.c1, atol=1e-10), kind
        assert np.allclose(r_shuffled.c2, r_fixed.c2, atol=1e-10), kind
    # The turning waters decorrelate, so a row mix-up would show.
    assert fixed.reorientation(kind="hh", oxygen_type=2, hydrogen_type=1).c1[4] < 0.9

    charges = np.tile([-0.8476, 0.4238, 0.4238], 4)
    for kind, kwargs in (
        ("bisector", {"oxygen_type": 2, "hydrogen_type": 1}),
        ("hh", {"oxygen_type": 2, "hydrogen_type": 1}),
        ("dipole", {}),
        ("dipole", {"charges": charges}),
    ):
        o_fixed = fixed.orientation_order(kind=kind, **kwargs)
        o_shuffled = shuffled.orientation_order(kind=kind, **kwargs)
        assert o_shuffled.mean_p1 == pytest.approx(o_fixed.mean_p1, abs=1e-12), kind
        assert o_shuffled.mean_p2 == pytest.approx(o_fixed.mean_p2, abs=1e-12), kind
        assert np.allclose(o_shuffled.cos_theta, o_fixed.cos_theta, atol=1e-12), kind

    for getter in (
        lambda traj: traj.molecule_dipoles(),
        lambda traj: traj.molecule_dipoles(charges=charges),
    ):
        for a, b in zip(getter(fixed), getter(shuffled)):
            assert np.array_equal(a.mol_ids, b.mol_ids)
            assert np.allclose(a.dipoles, b.dipoles, atol=1e-10)
            assert np.allclose(a.unit_vectors, b.unit_vectors, atol=1e-10)

    for a, b in zip(
        fixed.water_vectors(oxygen_type=2, hydrogen_type=1),
        shuffled.water_vectors(oxygen_type=2, hydrogen_type=1),
    ):
        assert np.array_equal(a.mol_ids, b.mol_ids)
        assert np.allclose(a.bisector, b.bisector, atol=1e-10)
        assert np.allclose(a.hh, b.hh, atol=1e-10)
        assert np.allclose(a.hh_unit, b.hh_unit, atol=1e-10)

    h_fixed = fixed.hbond_lifetime(oxygen_type=2, hydrogen_type=1)
    h_shuffled = shuffled.hbond_lifetime(oxygen_type=2, hydrogen_type=1)
    # Long lags are NaN when no origin still has a bond. NaN matches NaN.
    assert np.allclose(h_shuffled.c_hb, h_fixed.c_hb, atol=1e-12, equal_nan=True)
    assert np.allclose(h_shuffled.s_hb, h_fixed.s_hb, atol=1e-12, equal_nan=True)


def test_trajectory_identity_mismatches_are_refused(tmp_path: Path):
    box_len = 12.0
    truth, types, mols = _moving_waters(3, box_len)
    wrapped = truth - box_len * np.floor(truth / box_len)

    # An atom that changes molecule between frames.
    moved = np.stack([mols, mols, mols])
    moved[1, 3], moved[1, 0] = 1, 2
    moved_path = tmp_path / "moved.lammpstrj"
    _write_multiframe_lammps_dump(moved_path, wrapped, box_len, types, moved)
    traj = ds.open_trajectory(moved_path, all_atoms=True)
    for call in (
        traj.unwrap_centres,
        lambda: traj.msd(by_molecule=True),
        lambda: traj.reorientation(kind="dipole"),
        lambda: traj.reorientation(kind="bisector", oxygen_type=2, hydrogen_type=1),
        lambda: traj.reorientation(kind="hh", oxygen_type=2, hydrogen_type=1),
        lambda: traj.orientation_order(kind="dipole"),
        lambda: traj.orientation_order(kind="bisector", oxygen_type=2, hydrogen_type=1),
        lambda: traj.orientation_order(kind="hh", oxygen_type=2, hydrogen_type=1),
        traj.molecule_dipoles,
        lambda: traj.water_vectors(oxygen_type=2, hydrogen_type=1),
        lambda: traj.hbond_lifetime(oxygen_type=2, hydrogen_type=1),
    ):
        with pytest.raises(ValueError, match="molecule"):
            call()
    # Atom-level results do not depend on molecule membership.
    assert traj.unwrap_positions().shape == (3, 12, 3)

    # The same atom count, but different atom IDs in a later frame.
    ids_path = tmp_path / "ids.lammpstrj"
    frame_ids = np.stack([np.arange(1, 13), np.arange(1, 13), np.arange(101, 113)])
    _write_multiframe_lammps_dump(ids_path, wrapped, box_len, types, mols, ids=frame_ids)
    with pytest.raises(ValueError, match="atom IDs differ"):
        ds.open_trajectory(ids_path, all_atoms=True).unwrap_positions()


def test_gapped_permuted_two_water_order_and_molecule_swap(tmp_path: Path):
    """Two waters, three frames, gapped atom IDs, O/H1 line swaps, then a mol swap.

    Frames 1 and 2 list each water's oxygen and first hydrogen in the other
    order. Both orders must keep a +z dipole (mean P1 of 1). Swapping two
    atoms' molecule IDs in frame 1, while leaving the set {1, 2} unchanged,
    must be refused.
    """
    box_len = 20.0
    n_frames = 3
    n_atoms = 6
    types = np.tile([2, 1, 1], 2)
    mols = np.repeat([1, 2], 3)
    ids = np.array([100, 105, 110, 200, 205, 210])
    oh1 = np.array([0.757, 0.0, 0.586])
    oh2 = np.array([-0.757, 0.0, 0.586])
    xyz = np.zeros((n_frames, n_atoms, 3), dtype=np.float64)
    for f in range(n_frames):
        for w, base in enumerate((np.array([2.0, 2.0, 2.0]), np.array([8.0, 3.0, 2.0]))):
            centre = base + np.array([0.0, 0.0, 0.1 * f])
            xyz[f, 3 * w] = centre
            xyz[f, 3 * w + 1] = centre + oh1
            xyz[f, 3 * w + 2] = centre + oh2
    swap_lines = np.array([1, 0, 2, 4, 3, 5])
    orders = [np.arange(n_atoms), swap_lines, swap_lines]
    charges = np.tile([-0.8476, 0.4238, 0.4238], 2)
    fixed_path = tmp_path / "fixed.lammpstrj"
    shuffled_path = tmp_path / "shuffled.lammpstrj"
    _write_multiframe_lammps_dump(fixed_path, xyz, box_len, types, mols, ids=ids)
    _write_multiframe_lammps_dump(shuffled_path, xyz, box_len, types, mols, orders=orders, ids=ids)
    fixed = ds.open_trajectory(fixed_path, all_atoms=True, dt=1.0)
    shuffled = ds.open_trajectory(shuffled_path, all_atoms=True, dt=1.0)
    z_up = np.array([[0.0, 0.0, 1.0], [0.0, 0.0, 1.0]])

    for kwargs in ({}, {"charges": charges}):
        for traj in (fixed, shuffled):
            order = traj.orientation_order(kind="dipole", **kwargs)
            assert order.mean_p1 == pytest.approx(1.0, abs=1e-12)
            assert order.mean_p2 == pytest.approx(1.0, abs=1e-12)
        assert np.allclose(
            shuffled.orientation_order(kind="dipole", **kwargs).cos_theta,
            fixed.orientation_order(kind="dipole", **kwargs).cos_theta,
            atol=1e-12,
        )
        for a, b in zip(fixed.molecule_dipoles(**kwargs), shuffled.molecule_dipoles(**kwargs)):
            assert np.array_equal(a.mol_ids, b.mol_ids)
            assert np.allclose(a.unit_vectors, z_up, atol=1e-12)
            assert np.allclose(a.dipoles, b.dipoles, atol=1e-12)

    swapped = np.stack([mols, mols, mols])
    swapped[1, 0], swapped[1, 3] = int(swapped[1, 3]), int(swapped[1, 0])
    assert set(int(v) for v in swapped[1]) == {1, 2}
    assert not np.array_equal(swapped[1], mols)
    swap_path = tmp_path / "swap.lammpstrj"
    _write_multiframe_lammps_dump(swap_path, xyz, box_len, types, swapped, ids=ids)
    bad = ds.open_trajectory(swap_path, all_atoms=True, dt=1.0)
    for call in (
        bad.unwrap_centres,
        lambda: bad.msd(by_molecule=True),
        lambda: bad.reorientation(kind="dipole"),
        lambda: bad.reorientation(kind="bisector", oxygen_type=2, hydrogen_type=1),
        lambda: bad.reorientation(kind="hh", oxygen_type=2, hydrogen_type=1),
        lambda: bad.orientation_order(kind="dipole"),
        lambda: bad.orientation_order(kind="dipole", charges=charges),
        lambda: bad.orientation_order(kind="bisector", oxygen_type=2, hydrogen_type=1),
        lambda: bad.orientation_order(kind="hh", oxygen_type=2, hydrogen_type=1),
        bad.molecule_dipoles,
        lambda: bad.molecule_dipoles(charges=charges),
        lambda: bad.water_vectors(oxygen_type=2, hydrogen_type=1),
        lambda: bad.hbond_lifetime(oxygen_type=2, hydrogen_type=1),
    ):
        with pytest.raises(ValueError, match="molecule"):
            call()


def test_size_arguments_out_of_size_t_range_are_refused():
    size_max = 2**64 - 1
    for bad in (-1, -(2**63), 2**64):
        with pytest.raises(ValueError, match="n_particles"):
            ds.MsdAccumulator(bad)
        with pytest.raises(ValueError, match="max_lag"):
            ds.MsdAccumulator(1, max_lag=bad)
        with pytest.raises(ValueError, match="origin_stride"):
            ds.MsdAccumulator(1, origin_stride=bad)
        with pytest.raises(ValueError, match="max_lag"):
            ds.HbondLifetimeAccumulator(max_lag=bad)
        with pytest.raises(ValueError, match="n_vectors"):
            ds.ReorientationAccumulator(bad)
        with pytest.raises(ValueError, match="max_lag"):
            ds.ReorientationAccumulator(1, max_lag=bad)
    # In size_t range but too large to allocate: libseams refuses at create.
    for make in (
        lambda: ds.MsdAccumulator(size_max),
        lambda: ds.MsdAccumulator(size_max // 3 + 1),
        lambda: ds.MsdAccumulator(1, max_lag=size_max),
        lambda: ds.HbondLifetimeAccumulator(max_lag=size_max),
        lambda: ds.ReorientationAccumulator(size_max // 3 + 1),
        lambda: ds.ReorientationAccumulator(1, max_lag=size_max),
    ):
        with pytest.raises(ValueError, match="too large"):
            make()


def test_empty_and_one_frame_trajectories():
    empty = ds.TrajectorySequence.from_frames([])
    assert len(empty) == 0
    with pytest.raises(IndexError):
        _ = empty[0]
    for call in (
        empty.unwrap_positions,
        empty.unwrap_centres,
        empty.msd,
        lambda: empty.msd(by_molecule=True),
        empty.orientation_order,
        empty.reorientation,
        empty.hbond_lifetime,
    ):
        with pytest.raises(ValueError):
            call()
    with pytest.raises(ValueError, match="empty"):
        ds.hbond_lifetime([])
    with pytest.raises(ValueError, match="empty"):
        ds.reorientation_tcf([])

    pos = np.array([[[1.0, 2.0, 3.0], [4.0, 5.0, 6.0]]])
    one = ds.TrajectorySequence.from_arrays(pos, [10.0, 10.0, 10.0], dt=0.25)
    assert len(one) == 1
    assert np.allclose(one.unwrap_positions(), pos)
    centres = one.unwrap_centres(mol_ids=[7, 7])
    assert centres.mol_ids.tolist() == [7]
    assert np.allclose(centres.centres, [[[2.5, 3.5, 4.5]]])
    msd = one.msd()
    assert msd.lag_time.tolist() == [0.0]
    assert msd.total.tolist() == [0.0]
    assert msd.n_origins.tolist() == [1]
    assert msd.diffusion is None
    reo = one.reorientation(vectors=[np.array([[0.0, 0.0, 1.0], [1.0, 0.0, 0.0]])])
    assert reo.c1.tolist() == [1.0] and reo.c2.tolist() == [1.0]
    assert reo.relaxation is None
    hb = one.hbond_lifetime(edge_sets=[{(1, 2)}])
    assert hb.c_hb.tolist() == [1.0] and hb.s_hb.tolist() == [1.0]
    assert hb.lifetimes is None
