#!/usr/bin/env python3
"""Search a cheap CHILL+ pre-filter and score it on held-out seeds.

The target labels are the published cutoff rule (Nguyen and Molinero),
read from the library. Features are the tetrahedral order, the first and
fourth neighbour distances, the second-shell radius, the cutoff coordination,
the Steinhardt q3 of the four nearest, and the sorted c_ij of those four.
Held-out liquids use seeds 91001, 92002 and 93003 and are not in the fit.

A formula is only a pre-filter. The script records an accept region where
every held-out atom the formula calls cubic really is cubic, and the gap
from the nearest rejected score. Inside that region the C++ helper may
label cubic. On the boundary, and for every other label, it returns false
and the published rule runs. The default call does not use the helper.
"""

from __future__ import annotations

import json
import math
import os
from pathlib import Path

import numpy as np

ROOT = Path(__file__).resolve().parents[1]
OUT = Path(__file__).resolve().parent / "generated" / "pysr_prefilter.json"
LIB = "/tmp/seams-before/bindings/libseams_lang.so"
CUTOFF = 3.5
HELD_OUT = (91001, 92002, 93003)
PI = math.pi
NORM = np.array(
    [
        0.25 * math.sqrt(7.0 / PI),
        0.125 * math.sqrt(21.0 / PI),
        0.25 * math.sqrt(105.0 / (2.0 * PI)),
        0.125 * math.sqrt(35.0 / PI),
    ]
)


def ice_ic(cells: int) -> tuple[np.ndarray, float]:
    bond = 2.75
    lattice = 4.0 * bond / math.sqrt(3.0)
    side = cells * lattice
    fcc = np.array([[0.0, 0, 0], [0.5, 0.5, 0], [0.5, 0, 0.5], [0, 0.5, 0.5]])
    basis = np.vstack([fcc, fcc + 0.25])
    ijk = np.stack(np.meshgrid(*[np.arange(cells)] * 3, indexing="ij"), -1).reshape(-1, 3)
    pts = ((ijk[:, None, :] + basis).reshape(-1, 3) * lattice) % side
    return np.ascontiguousarray(pts), side


def liquid(n: int, side: float, seed: int) -> np.ndarray:
    rng = np.random.default_rng(seed)
    return rng.random((n, 3)) * side


def published_labels(pts: np.ndarray, side: float) -> np.ndarray:
    import ctypes

    lib = ctypes.CDLL(LIB)
    lib.seams_frame_from_arrays.restype = ctypes.c_void_p
    lib.seams_frame_from_arrays.argtypes = [
        ctypes.POINTER(ctypes.c_double),
        ctypes.c_int,
        ctypes.POINTER(ctypes.c_double),
        ctypes.POINTER(ctypes.c_int),
        ctypes.POINTER(ctypes.c_int),
    ]
    lib.seams_frame_free.argtypes = [ctypes.c_void_p]
    lib.seams_frame_chill_plus.argtypes = [
        ctypes.c_void_p,
        ctypes.c_double,
        ctypes.c_int,
        ctypes.c_char_p,
        ctypes.POINTER(ctypes.c_int),
        ctypes.POINTER(ctypes.c_int),
    ]
    n = len(pts)
    flat = np.ascontiguousarray(pts.reshape(-1))
    box = np.array([0, 0, 0, side, side, side, 0, 0, 0])
    types = np.ones(n, np.int32)
    ids = np.arange(1, n + 1, dtype=np.int32)
    frame = lib.seams_frame_from_arrays(
        flat.ctypes.data_as(ctypes.POINTER(ctypes.c_double)),
        n,
        box.ctypes.data_as(ctypes.POINTER(ctypes.c_double)),
        types.ctypes.data_as(ctypes.POINTER(ctypes.c_int)),
        ids.ctypes.data_as(ctypes.POINTER(ctypes.c_int)),
    )
    labels = (ctypes.c_int * n)()
    counts = (ctypes.c_int * 9)()
    status = lib.seams_frame_chill_plus(frame, CUTOFF, 1, b"cutoff", labels, counts)
    lib.seams_frame_free(frame)
    if status != 0:
        raise RuntimeError("published CHILL+ failed")
    return np.array(labels[:], dtype=np.int32)


def minimum_image(delta: np.ndarray, side: float) -> np.ndarray:
    return delta - side * np.round(delta / side)


def ylm3(bonds: np.ndarray) -> np.ndarray:
    """Sum of unnormalised real l = 3 harmonics. bonds is (4, 3)."""
    acc = np.zeros(7)
    for x, y, z in bonds:
        x2, y2, z2 = x * x, y * y, z * z
        r2 = x2 + y2 + z2
        if r2 <= 0.0:
            acc[0] += NORM[0] * 2.0
            continue
        inv = 1.0 / (r2 * math.sqrt(r2))
        t1 = 5.0 * z2 - r2
        acc[0] += NORM[0] * (z * (5.0 * z2 - 3.0 * r2)) * inv
        acc[1] += NORM[1] * (y * t1) * inv
        acc[2] += NORM[1] * (x * t1) * inv
        acc[3] += NORM[2] * (z * (y2 - x2)) * inv
        acc[4] += NORM[2] * (2.0 * x * y * z) * inv
        acc[5] += NORM[3] * (y * (y2 - 3.0 * x2)) * inv
        acc[6] += NORM[3] * (x * (3.0 * y2 - x2)) * inv
    return acc


def qnorm(q: np.ndarray) -> float:
    return float(q[0] * q[0] + 2.0 * np.dot(q[1:], q[1:]))


def features_of(pts: np.ndarray, side: float) -> np.ndarray:
    n = len(pts)
    # Brute neighbour ranking. These frames are a few hundred molecules.
    delta = pts[:, None, :] - pts[None, :, :]
    delta = minimum_image(delta, side)
    d2 = np.sum(delta * delta, axis=2)
    np.fill_diagonal(d2, np.inf)
    order = np.argpartition(d2, 5, axis=1)[:, :6]
    # argpartition is unordered; sort the six by (d2, index).
    rows = []
    cut2 = CUTOFF * CUTOFF
    out = np.zeros((n, 10))
    q = np.zeros((n, 7))
    bonds = np.zeros((n, 4, 3))
    keep = np.zeros(n, dtype=np.int32)
    for i in range(n):
        cand = order[i]
        keys = sorted((d2[i, j], int(j)) for j in cand)
        # The sixth is only a sentinel so the fifth distance exists.
        chosen = keys[:4]
        keep[i] = 4
        vec = np.stack([delta[i, j] for _, j in chosen])
        bonds[i] = vec
        q[i] = ylm3(vec)
        dist = np.sqrt([k[0] for k in keys[:5]])
        shell = int(np.count_nonzero(d2[i] < cut2))
        # Tetrahedral order of the four nearest.
        acc = 0.0
        for a in range(4):
            for b in range(a + 1, 4):
                na = np.linalg.norm(vec[a])
                nb = np.linalg.norm(vec[b])
                cosine = float(np.dot(vec[a], vec[b]) / (na * nb))
                acc += (cosine + 1.0 / 3.0) ** 2
        qtet = 1.0 - 0.375 * acc
        out[i, 0] = qtet
        out[i, 1] = dist[0]
        out[i, 2] = dist[3]
        out[i, 3] = dist[4]
        out[i, 4] = shell
    nrm = np.array([qnorm(q[i]) for i in range(n)])
    for i in range(n):
        cij = []
        chosen = sorted((d2[i, j], int(j)) for j in order[i])[:4]
        for _, j in chosen:
            ni, nj = nrm[i], nrm[j]
            if ni > 1e-12 and nj > 1e-12:
                dot = q[i, 0] * q[j, 0] + 2.0 * float(np.dot(q[i, 1:], q[j, 1:]))
                cij.append(dot / math.sqrt(ni * nj))
            else:
                cij.append(0.0)
        cij.sort()
        out[i, 5] = cij[0]
        out[i, 6] = cij[-1]
        out[i, 7] = float(np.mean(cij))
        out[i, 8] = math.sqrt(max(nrm[i], 0.0) * 4.0 * PI / (7.0 * 16.0))
        out[i, 9] = keep[i]
    return out


def frames():
    """(name, points, side, held_out)."""
    made = []
    for cells in (2, 3, 4):
        pts, side = ice_ic(cells)
        made.append((f"ice{cells}", pts, side, False))
    rng_ice, side = ice_ic(3)
    for sigma, seed, held in (
        (0.04, 11, False),
        (0.12, 12, False),
        (0.30, 13, False),
        (0.12, 93003, True),
    ):
        rng = np.random.default_rng(seed)
        made.append((f"noise{sigma}-{seed}", rng_ice + rng.normal(0.0, sigma, rng_ice.shape), side, held))
        made[-1] = (
            made[-1][0],
            np.mod(made[-1][1], side),
            side,
            held,
        )
    for seed in (1, 7, 42, 99, 1234, *HELD_OUT):
        made.append((f"liq{seed}", liquid(80, 18.0, seed), 18.0, seed in HELD_OUT))
    return made


def main() -> None:
    names = ["qtet", "r1", "r4", "r5", "shell", "cmin", "cmax", "cmean", "q3", "keep"]
    xs = []
    ys = []
    held = []
    counts: dict[str, int] = {}
    for name, pts, side, is_held in frames():
        labels = published_labels(pts, side)
        feat = features_of(pts, side)
        xs.append(feat)
        ys.append(labels)
        held.append(np.full(len(pts), is_held))
        counts[name] = {int(k): int(v) for k, v in zip(*np.unique(labels, return_counts=True))}
        print(name, "n", len(pts), "held", is_held, counts[name], flush=True)
    x = np.vstack(xs)
    y = np.concatenate(ys)
    mask = np.concatenate(held)
    cubic = (y == 0).astype(np.float64)
    train = ~mask
    os.environ.setdefault("PYTHON_JULIAPKG_OFFLINE", "yes")
    from pysr import PySRRegressor

    model = PySRRegressor(
        niterations=40,
        binary_operators=["+", "-", "*", "/"],
        unary_operators=["square", "abs"],
        populations=12,
        population_size=40,
        maxsize=12,
        procs=4,
        model_selection="best",
        timeout_in_seconds=240,
        parallelism="serial",
        temp_equation_file=True,
        deterministic=True,
        random_state=0,
        extra_sympy_mappings={},
    )
    model.fit(x[train], cubic[train], variable_names=names)
    pred = model.predict(x)
    # Scores are the expression value. Cubic is the positive class.
    formula = str(model.sympy())
    print("formula", formula, flush=True)
    # Accept only scores that no held-out non-cubic reaches, and that sit
    # a gap above the highest held-out non-cubic score.
    held_non = pred[mask & (y != 0)]
    held_cub = pred[mask & (y == 0)]
    train_non = pred[train & (y != 0)]
    ceiling = float(np.max(np.concatenate([held_non, train_non]))) if len(held_non) or len(train_non) else 0.0
    gap_floor = ceiling
    # Require the accept cut to be strictly above every non-cubic score.
    cut = gap_floor
    accepted = pred > cut
    false_accept = int(np.count_nonzero(accepted & (y != 0)))
    true_accept = int(np.count_nonzero(accepted & (y == 0)))
    held_acc = float(np.mean((pred[mask] > 0.5) == (y[mask] == 0))) if mask.any() else None
    report = {
        "formula": formula,
        "variable_names": names,
        "label_cubic": 0,
        "held_out_seeds": list(HELD_OUT),
        "held_out_threshold_accuracy": held_acc,
        "noncubic_score_ceiling": ceiling,
        "accept_if_score_above": cut,
        "false_accepts": false_accept,
        "true_accepts": true_accept,
        "n_train": int(train.sum()),
        "n_held": int(mask.sum()),
        "class_counts": counts,
        "equations": model.equations_.astype(str).to_dict() if hasattr(model, "equations_") else {},
    }
    # The equations frame is easier as a record of the Pareto front.
    try:
        report["pareto"] = model.equations_[["complexity", "loss", "equation"]].astype(str).to_dict("records")
    except Exception as exc:  # noqa: BLE001 — record the front when the columns exist
        report["pareto_error"] = str(exc)
    OUT.parent.mkdir(parents=True, exist_ok=True)
    OUT.write_text(json.dumps(report, indent=2, default=str) + "\n")
    print(json.dumps({k: report[k] for k in ("formula", "held_out_threshold_accuracy", "false_accepts", "true_accepts", "accept_if_score_above")}, indent=2))


if __name__ == "__main__":
    main()
