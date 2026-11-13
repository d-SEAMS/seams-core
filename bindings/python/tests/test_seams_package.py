# Copyright (c) d-SEAMS developers
"""The ctypes package on the cubic cell, and the documented error path."""

from __future__ import annotations

import math
import os
import subprocess
import sys
from pathlib import Path

import numpy as np
import pytest

import seams
from seams import CUTOFF, OVITO, PUBLISHED, SEEDED, SeamsError, System


def repo_root() -> Path:
    here = Path(__file__).resolve()
    for parent in here.parents:
        if (parent / "include" / "seams.h").is_file():
            return parent
    raise RuntimeError("include/seams.h was not found")


def _library() -> None:
    if os.environ.get("SEAMS_LIB"):
        return
    root = repo_root()
    for rel in ("bbdir/src/libseams.so", "build/src/libseams.so", "bbdir-py/src/libseams.so"):
        path = root / rel
        if path.is_file():
            os.environ["SEAMS_LIB"] = str(path)
            return


def diamond(nrep: int = 2) -> tuple[np.ndarray, float]:
    a = 2.75 * 4.0 / math.sqrt(3.0)
    fcc = ((0.0, 0.0, 0.0), (0.5, 0.5, 0.0), (0.5, 0.0, 0.5), (0.0, 0.5, 0.5))
    basis = ((0.0, 0.0, 0.0), (0.25, 0.25, 0.25))
    points = [
        (
            (i + shift[0] + site[0]) * a,
            (j + shift[1] + site[1]) * a,
            (k + shift[2] + site[2]) * a,
        )
        for i in range(nrep)
        for j in range(nrep)
        for k in range(nrep)
        for shift in fcc
        for site in basis
    ]
    return np.asarray(points, dtype=np.float64), nrep * a


@pytest.fixture(scope="module", autouse=True)
def _loaded() -> None:
    _library()


def test_threads_version_and_cubic_ice() -> None:
    seams.set_num_threads(1)
    assert seams.num_threads() == 1
    assert seams.version()
    assert seams.abi() == 1
    assert seams.version_tuple()[0] >= 2
    pos, edge = diamond()
    with System() as sys:
        sys.set_ortho(pos, (edge, edge, edge))
        assert sys.count() == pos.shape[0]
        published = sys.chill_plus(rule=PUBLISHED)
        assert published.labels.shape == (pos.shape[0],)
        assert int(published.counts[0]) == pos.shape[0]
        assert np.all(published.labels == 0)
        ovito = sys.chill_plus(graph=seams.MUTUAL, rule=OVITO)
        assert int(ovito.counts[0]) == pos.shape[0]
        order = sys.q4q6()
        assert float(order.q6.mean()) > 0.3
        assert np.all(np.isfinite(order.q4))
        with sys.neighbours(graph=CUTOFF) as nb:
            csr = nb.copy()
        assert csr.offsets.shape[0] == pos.shape[0] + 1
        assert np.all(np.diff(csr.offsets) == 4)
        with sys.rings() as rings:
            hist = rings.histogram(7)
            assert int(hist[6]) > 0
            assert rings.fraction(6) > 0.5
            assert rings.labels().sizes.shape == (rings.count(),)
        bar = sys.q6bar()
        assert bar.n_ice == pos.shape[0]
        assert bar.mean > 0.3
        assert bar.values is not None


def test_missing_frame_uses_handle_error() -> None:
    with System() as sys:
        with pytest.raises(SeamsError) as caught:
            sys.chill_plus()
    assert caught.value.status != 0
    assert caught.value.message
    assert "no frame" in caught.value.message or caught.value.status == seams.ERR_STATE


def test_closed_handle_and_bad_box() -> None:
    sys = System()
    sys.close()
    assert sys.closed
    sys.close()
    with pytest.raises(SeamsError):
        sys.count()
    with System() as fresh:
        pos, _edge = diamond(1)
        with pytest.raises(ValueError):
            fresh.set_frame(pos, [1.0, 2.0])


def test_seeded_graph_is_rejected_for_neighbours() -> None:
    pos, edge = diamond(1)
    with System() as sys:
        sys.set_frame(pos, (edge, edge, edge))
        with pytest.raises(SeamsError) as caught:
            sys.neighbours(graph=SEEDED)
    assert caught.value.status == seams.ERR_ARG


def test_read_known_trajectories() -> None:
    root = repo_root()
    expect = {
        "input/traj/ice1h.lammpstrj": 128,
        "input/traj/genice_sI.lammpstrj": 46,
        "input/xyz/configuration.xyz": 256,
    }
    for rel, count in expect.items():
        path = root / rel
        if not path.is_file():
            pytest.fail(f"missing {rel}")
        with System(path, type_filter=-1) as sys:
            assert sys.count() == count
        if path.suffix == ".lammpstrj":
            assert seams.nframes(path) >= 1


def test_signature_on_si() -> None:
    path = repo_root() / "input/traj/genice_sI.lammpstrj"
    with System(path, type_filter=1) as sys:
        found = sys.signature("512")
    assert found.n_cages == 2
    assert found.n_atoms == 40


def test_documented_quickstart_output() -> None:
    root = repo_root()
    script = root / "docs" / "examples" / "seams_ctypes_quick.py"
    env = os.environ.copy()
    env["PYTHONPATH"] = str(root / "bindings" / "python" / "src")
    out = subprocess.check_output(
        [sys.executable, str(script), env["SEAMS_LIB"]],
        env=env,
        text=True,
    )
    page = (root / "docs" / "api" / "python.md").read_text()
    assert out.strip() in page


def test_pydseams_duck_type() -> None:
    class Frame:
        def __init__(self, positions: np.ndarray, box: np.ndarray) -> None:
            self.positions = positions
            self.box = box

    pos, edge = diamond(1)
    copied = seams.from_pydseams(Frame(pos, np.array([edge, edge, edge])))
    try:
        labels = copied.chill_plus().labels
        assert np.all(labels == 0)
    finally:
        copied.close()
    with pytest.raises(TypeError):
        seams.from_pydseams(object())
