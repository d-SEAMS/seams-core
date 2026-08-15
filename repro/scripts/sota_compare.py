#!/usr/bin/env python3
"""Measured comparison of per-molecule ice classifiers.

Systems with ground truth by construction: a cubic-diamond oxygen lattice
(ice Ic topology), a lonsdaleite lattice (ice Ih topology), a Gaussian
jitter sweep over both, and a dense disordered null. Both generators
assert four neighbours at the bond length under periodic boundaries, so a
wrong lattice fails loudly instead of poisoning the ground truth.

Methods:
  dseams-topo   per-atom labels from HC/DDC cage affiliation (this work)
  chill+        d-SEAMS CHILL+ (four nearest neighbours)
  freud-q6      Steinhardt q6 threshold, ice versus not (cannot split Ic/Ih)
  freud-ld      nearest-centroid on Lechner-Dellago (q4bar, q6bar)
  ovito-ptm     polyhedral template matching, when the ovito wheel imports
  twist-op      Matsumoto/Yagasaki/Tanaka torsional order parameter,
                nearest-centroid on the ideal lattices, when installed

Scalar/vector methods get their decision fit on the ideal (sigma = 0)
lattices and applied unchanged across the sweep; label-producing methods
use their native labels. Metrics: accuracy on Ic and Ih, Ic/Ih confusion,
false-crystal rate on the null. Output is key=value lines per
(method, system, sigma).
"""

from __future__ import annotations

import importlib.util
import sys
import time

import numpy as np

BOND = 2.75  # O-O first-shell distance in ice, Angstrom
CUTOFF = 3.5
SIGMAS = [0.0, 0.05, 0.10, 0.15, 0.20, 0.25, 0.30]
SEED = 88172645463325252


def cubic_diamond(reps: int):
    a = 4.0 * BOND / np.sqrt(3.0)
    fcc = [(0, 0, 0), (0.5, 0.5, 0), (0.5, 0, 0.5), (0, 0.5, 0.5)]
    basis = [(x, y, z) for (x, y, z) in fcc] + [
        (x + 0.25, y + 0.25, z + 0.25) for (x, y, z) in fcc
    ]
    pos = []
    for i in range(reps):
        for j in range(reps):
            for k in range(reps):
                for b in basis:
                    pos.append(((i + b[0]) * a, (j + b[1]) * a, (k + b[2]) * a))
    box = np.array([reps * a] * 3)
    return np.asarray(pos) % box, box


def lonsdaleite(reps_x: int, reps_y: int, reps_z: int):
    # Wurtzite motif with the ideal parameter u = 3/8: hexagonal cell with
    # c/a = sqrt(8/3), bond length (3/8)c along c and equal laterally. The
    # periodic box is the orthohexagonal double cell (a, a*sqrt(3), c) with
    # the centring copy at (a/2, a*sqrt(3)/2).
    a = BOND / ((3.0 / 8.0) * np.sqrt(8.0 / 3.0))
    c = a * np.sqrt(8.0 / 3.0)
    hexagonal = [
        (1.0 / 3.0, 2.0 / 3.0, 0.0),
        (2.0 / 3.0, 1.0 / 3.0, 0.5),
        (1.0 / 3.0, 2.0 / 3.0, 3.0 / 8.0),
        (2.0 / 3.0, 1.0 / 3.0, 7.0 / 8.0),
    ]
    X, Y = a, a * np.sqrt(3.0)
    ortho = []
    for (hx, hy, hz) in hexagonal:
        # Hexagonal lattice vectors a1 = (a, 0), a2 = (-a/2, a*sqrt(3)/2)
        x = (hx - 0.5 * hy) * a
        y = hy * (np.sqrt(3.0) / 2.0) * a
        for (sx, sy) in ((0.0, 0.0), (X / 2.0, Y / 2.0)):
            ortho.append(((x + sx) % X, (y + sy) % Y, hz * c))
    box = np.array([reps_x * X, reps_y * Y, reps_z * c])
    pos = []
    for i in range(reps_x):
        for j in range(reps_y):
            for k in range(reps_z):
                for (ox, oy, oz) in ortho:
                    pos.append((i * X + ox, j * Y + oy, k * c + oz))
    return np.asarray(pos) % box, box


def dense_null(n: int, box: np.ndarray, rng):
    # Rejection-sampled dense disorder: no pair closer than 2.3 A, no
    # crystalline order by construction
    pos = []
    cell = 2.3
    while len(pos) < n:
        cand = rng.random(3) * box
        ok = True
        for p in pos[-400:] if len(pos) > 400 else pos:
            d = np.abs(cand - p)
            d = np.minimum(d, box - d)
            if (d * d).sum() < cell * cell:
                ok = False
                break
        if ok:
            pos.append(cand)
    return np.asarray(pos)


def neighbour_check(pos: np.ndarray, box: np.ndarray):
    # Every atom must have exactly four neighbours at BOND +- 0.1 under PBC
    n = len(pos)
    counts = np.zeros(n, dtype=int)
    for i in range(n):
        d = np.abs(pos - pos[i])
        d = np.minimum(d, box - d)
        r = np.sqrt((d * d).sum(axis=1))
        counts[i] = int(((r > 0.1) & (np.abs(r - BOND) < 0.1)).sum())
    assert (counts == 4).all(), f"lattice broken: counts {np.unique(counts)}"


def jitter(pos: np.ndarray, box: np.ndarray, sigma: float, rng):
    return (pos + rng.normal(0.0, sigma, pos.shape)) % box if sigma else pos


# ---------------------------------------------------------------- d-SEAMS
import _core  # noqa: E402


def make_cloud(pos: np.ndarray, box: np.ndarray):
    cloud = _core.PointCloudDouble()
    pts = []
    id_map = {}
    for i, (x, y, z) in enumerate(pos):
        p = _core.PointDouble()
        p.c_type = 1
        p.atomID = i + 1
        p.molID = i + 1
        p.x, p.y, p.z = float(x), float(y), float(z)
        pts.append(p)
        id_map[i + 1] = i
    cloud.pts = pts
    cloud.idIndexMap = id_map
    cloud.box = [float(b) for b in box]
    cloud.boxLow = [0.0, 0.0, 0.0]
    cloud.nop = len(pos)
    cloud.currentFrame = 1
    return cloud


def dseams_lists(pos, box):
    cloud = make_cloud(pos, box)
    nl = _core.neighListO(CUTOFF, cloud, 1)
    idx = _core.neighbourListByIndex(cloud, nl)
    return cloud, nl, idx


def _affiliation_labels(cloud, idx, n):
    rings = [r for r in _core.ringNetwork(idx, 7) if len(r) == 6]
    if not rings:
        return ["water"] * n
    hc, ddc = _core.cageAffiliation(rings, idx)
    in_hc = np.zeros(n, dtype=bool)
    in_ddc = np.zeros(n, dtype=bool)
    for ring, h, d in zip(rings, hc, ddc):
        for atom in ring:
            in_hc[atom] |= h
            in_ddc[atom] |= d
    labels = []
    for i in range(n):
        if in_hc[i] and in_ddc[i]:
            labels.append("mixed")
        elif in_hc[i]:
            labels.append("hexagonal")
        elif in_ddc[i]:
            labels.append("cubic")
        else:
            labels.append("water")
    return labels


def dseams_topo_4nn(pos, box, mutual=True):
    """TUM v2: the same cage predicates on the mutual 4-nearest bonded
    graph, which keeps neighbour identity where a cutoff loses it while a
    one-sided nomination never creates a bond."""
    cloud = make_cloud(pos, box)
    knn = _core.kNearestNeighbourList(cloud, 4, 5.0, 1, mutual)
    idx = _core.neighbourListByIndex(cloud, knn)
    return _affiliation_labels(cloud, idx, len(pos))


def dseams_topo(pos, box):
    """Per-atom label from cage affiliation: an atom in any HC-affiliated
    six-ring is hexagonal, in a DDC-affiliated ring cubic, in both mixed,
    in neither water. Six-rings come from the primitive ring network."""
    cloud, nl, idx = dseams_lists(pos, box)
    rings = [r for r in _core.ringNetwork(idx, 7) if len(r) == 6]
    n = len(pos)
    if not rings:
        return ["water"] * n
    hc, ddc = _core.cageAffiliation(rings, idx)
    in_hc = np.zeros(n, dtype=bool)
    in_ddc = np.zeros(n, dtype=bool)
    for ring, h, d in zip(rings, hc, ddc):
        for atom in ring:
            in_hc[atom] |= h
            in_ddc[atom] |= d
    labels = []
    for i in range(n):
        if in_hc[i] and in_ddc[i]:
            labels.append("mixed")
        elif in_hc[i]:
            labels.append("hexagonal")
        elif in_ddc[i]:
            labels.append("cubic")
        else:
            labels.append("water")
    return labels


def chill_plus(pos, box, scratch):
    cloud, nl, idx = dseams_lists(pos, box)
    _core.getCorrelPlus(cloud, nl, False)
    _core.getIceTypePlus(cloud, nl, scratch, 1, False, "sota.txt")
    state = {
        _core.AtomStateType.cubic: "cubic",
        _core.AtomStateType.reCubic: "cubic",
        _core.AtomStateType.hexagonal: "hexagonal",
        _core.AtomStateType.reHex: "hexagonal",
        _core.AtomStateType.interfacial: "interfacial",
        _core.AtomStateType.clathrate: "clathrate",
        _core.AtomStateType.interClathrate: "clathrate",
    }
    return [state.get(pt.iceType, "water") for pt in cloud.pts]


# ------------------------------------------------------------------ freud
def freud_features(pos, box):
    import freud

    fbox = freud.box.Box(Lx=box[0], Ly=box[1], Lz=box[2])
    shifted = pos - box / 2.0
    feats = {}
    for ell, average in ((4, True), (6, True), (6, False)):
        op = freud.order.Steinhardt(l=ell, average=average)
        # Four nearest neighbours, as the water literature uses; a cutoff
        # list inflates ql on low-coordination disordered particles
        op.compute((fbox, shifted), neighbors={"num_neighbors": 4})
        feats[(ell, average)] = np.asarray(op.particle_order)
    return feats


# --------------------------------------------------------------- optional
def have(mod):
    return importlib.util.find_spec(mod) is not None


def ovito_ptm(pos, box):
    from ovito.data import DataCollection, ParticleType, SimulationCell
    from ovito.modifiers import PolyhedralTemplateMatchingModifier
    from ovito.pipeline import Pipeline, StaticSource

    data = DataCollection()
    cell = data.create_cell(
        matrix=[
            [box[0], 0, 0, 0],
            [0, box[1], 0, 0],
            [0, 0, box[2], 0],
        ],
        pbc=(True, True, True),
    )
    particles = data.create_particles(count=len(pos))
    particles.create_property("Position", data=pos)
    pipeline = Pipeline(source=StaticSource(data=data))
    mod = PolyhedralTemplateMatchingModifier(rmsd_cutoff=0.15)
    mod.structures[
        PolyhedralTemplateMatchingModifier.Type.CUBIC_DIAMOND
    ].enabled = True
    mod.structures[
        PolyhedralTemplateMatchingModifier.Type.HEX_DIAMOND
    ].enabled = True
    pipeline.modifiers.append(mod)
    out = pipeline.compute()
    types = np.asarray(out.particles["Structure Type"])
    cubic = int(PolyhedralTemplateMatchingModifier.Type.CUBIC_DIAMOND)
    hexd = int(PolyhedralTemplateMatchingModifier.Type.HEX_DIAMOND)
    labels = []
    for t in types:
        if t == cubic:
            labels.append("cubic")
        elif t == hexd:
            labels.append("hexagonal")
        else:
            labels.append("water")
    return labels


def twist_features(pos, box):
    # twist-op predates the numpy 1.24 alias removals; restore the alias the
    # way its own era provided it
    if not hasattr(np, "float"):
        np.float = float  # noqa: NPY001
    from twist_op import twist_iter

    # Bond list under minimum image, each bond once; normalized vectors
    n = len(pos)
    pairs = []
    vecs = []
    for i in range(n):
        d = pos[i + 1 :] - pos[i]
        d -= box * np.round(d / box)
        r = np.sqrt((d * d).sum(axis=1))
        for off in np.nonzero(r < CUTOFF)[0]:
            j = i + 1 + int(off)
            pairs.append((i, j))
            vecs.append(d[off] / r[off])
    # Per-molecule mean of the real part of the twist: staggered and
    # eclipsed conformations sit at opposite signs of Re exp(3iA)
    per_mol = np.zeros(n)
    counts = np.zeros(n)
    for (i, j), tw in twist_iter(np.asarray(pairs), np.asarray(vecs)):
        per_mol[i] += tw.real
        per_mol[j] += tw.real
        counts[i] += 1
        counts[j] += 1
    counts[counts == 0] = 1
    return per_mol / counts


# ---------------------------------------------------------------- scoring
def score(labels, truth):
    labels = np.asarray(labels)
    n = len(labels)
    return {
        "acc": float((labels == truth).sum() / n),
        "cubic": float((labels == "cubic").sum() / n),
        "hex": float((labels == "hexagonal").sum() / n),
        "crystal": float(
            ((labels == "cubic") | (labels == "hexagonal")).sum() / n
        ),
    }


def centroid_classifier(train):
    # train: {label: feature matrix}; returns a nearest-centroid function
    cents = {k: v.mean(axis=0) for k, v in train.items()}

    def classify(feats):
        labels = []
        for row in feats:
            best = min(cents, key=lambda k: np.linalg.norm(row - cents[k]))
            labels.append(best)
        return labels

    return classify


def main():
    import os
    import tempfile

    rng = np.random.default_rng(SEED)
    scratch = tempfile.mkdtemp() + os.sep

    ic_pos, ic_box = cubic_diamond(4)
    ih_pos, ih_box = lonsdaleite(5, 3, 3)
    neighbour_check(ic_pos, ic_box)
    neighbour_check(ih_pos, ih_box)
    null_pos = dense_null(360, ic_box, rng)
    print(f"# ic n={len(ic_pos)} ih n={len(ih_pos)} null n={len(null_pos)}")

    # Fit the scalar/vector decisions on the ideal lattices plus the null
    ld_train = {
        "cubic": np.column_stack(
            [
                freud_features(ic_pos, ic_box)[(4, True)],
                freud_features(ic_pos, ic_box)[(6, True)],
            ]
        ),
        "hexagonal": np.column_stack(
            [
                freud_features(ih_pos, ih_box)[(4, True)],
                freud_features(ih_pos, ih_box)[(6, True)],
            ]
        ),
        "water": np.column_stack(
            [
                freud_features(null_pos, ic_box)[(4, True)],
                freud_features(null_pos, ic_box)[(6, True)],
            ]
        ),
    }
    ld_classify = centroid_classifier(ld_train)
    q6_ice_threshold = 0.5 * (
        freud_features(ic_pos, ic_box)[(6, False)].mean()
        + freud_features(null_pos, ic_box)[(6, False)].mean()
    )

    twist_ok = have("twist_op")
    twist_classify = None
    if twist_ok:
        try:
            tw_train = {
                "cubic": twist_features(ic_pos, ic_box).reshape(-1, 1),
                "hexagonal": twist_features(ih_pos, ih_box).reshape(-1, 1),
                "water": twist_features(null_pos, ic_box).reshape(-1, 1),
            }
            twist_classify = centroid_classifier(tw_train)
        except Exception as exc:  # noqa: BLE001
            print(f"# twist-op failed: {exc}")
            twist_ok = False

    ptm_ok = have("ovito")

    systems = [("ic", ic_pos, ic_box, "cubic"), ("ih", ih_pos, ih_box, "hexagonal")]
    for name, pos0, box, truth in systems:
        for sigma in SIGMAS:
            pos = jitter(pos0, box, sigma, np.random.default_rng(SEED + int(sigma * 1000)))

            for method, run in (
                ("dseams-topo", lambda: dseams_topo(pos, box)),
                ("dseams-topo-4nn", lambda: dseams_topo_4nn(pos, box)),
                ("chill+", lambda: chill_plus(pos, box, scratch)),
            ):
                t0 = time.perf_counter()
                labels = run()
                dt = time.perf_counter() - t0
                s = score(labels, truth)
                print(
                    f"method={method} system={name} sigma={sigma:.2f} "
                    f"acc={s['acc']:.3f} cubic={s['cubic']:.3f} "
                    f"hex={s['hex']:.3f} crystal={s['crystal']:.3f} "
                    f"ms={dt*1e3:.1f}"
                )

            feats = freud_features(pos, box)
            ld = ld_classify(
                np.column_stack([feats[(4, True)], feats[(6, True)]])
            )
            s = score(ld, truth)
            print(
                f"method=freud-ld system={name} sigma={sigma:.2f} "
                f"acc={s['acc']:.3f} cubic={s['cubic']:.3f} "
                f"hex={s['hex']:.3f} crystal={s['crystal']:.3f} ms=-1"
            )
            q6 = feats[(6, False)]
            ice = q6 > q6_ice_threshold
            print(
                f"method=freud-q6 system={name} sigma={sigma:.2f} "
                f"acc=nan cubic=nan hex=nan "
                f"crystal={float(ice.sum())/len(q6):.3f} ms=-1 "
                f"note=ice-vs-not-only"
            )

            if ptm_ok:
                try:
                    labels = ovito_ptm(pos, box)
                    s = score(labels, truth)
                    print(
                        f"method=ovito-ptm system={name} sigma={sigma:.2f} "
                        f"acc={s['acc']:.3f} cubic={s['cubic']:.3f} "
                        f"hex={s['hex']:.3f} crystal={s['crystal']:.3f} ms=-1"
                    )
                except Exception as exc:  # noqa: BLE001
                    print(f"# ovito-ptm failed: {exc}")
                    ptm_ok = False

            if twist_ok and twist_classify is not None:
                tw = twist_classify(twist_features(pos, box).reshape(-1, 1))
                s = score(tw, truth)
                print(
                    f"method=twist-op system={name} sigma={sigma:.2f} "
                    f"acc={s['acc']:.3f} cubic={s['cubic']:.3f} "
                    f"hex={s['hex']:.3f} crystal={s['crystal']:.3f} ms=-1"
                )

    # False-crystal rate on the null, per label-producing method
    for method, run in (
        ("dseams-topo", lambda: dseams_topo(null_pos, ic_box)),
        ("dseams-topo-4nn", lambda: dseams_topo_4nn(null_pos, ic_box)),
        ("chill+", lambda: chill_plus(null_pos, ic_box, scratch)),
    ):
        labels = run()
        s = score(labels, "water")
        print(
            f"method={method} system=null sigma=0.00 acc={s['acc']:.3f} "
            f"false_crystal={s['crystal']:.3f}"
        )
    feats = freud_features(null_pos, ic_box)
    ld = ld_classify(np.column_stack([feats[(4, True)], feats[(6, True)]]))
    s = score(ld, "water")
    print(
        f"method=freud-ld system=null sigma=0.00 acc={s['acc']:.3f} "
        f"false_crystal={s['crystal']:.3f}"
    )
    q6 = feats[(6, False)]
    print(
        f"method=freud-q6 system=null sigma=0.00 acc=nan "
        f"false_crystal={float((q6 > q6_ice_threshold).sum())/len(q6):.3f}"
    )
    if ptm_ok:
        labels = ovito_ptm(null_pos, ic_box)
        s = score(labels, "water")
        print(
            f"method=ovito-ptm system=null sigma=0.00 acc={s['acc']:.3f} "
            f"false_crystal={s['crystal']:.3f}"
        )
    if twist_ok and twist_classify is not None:
        tw = twist_classify(twist_features(null_pos, ic_box).reshape(-1, 1))
        s = score(tw, "water")
        print(
            f"method=twist-op system=null sigma=0.00 acc={s['acc']:.3f} "
            f"false_crystal={s['crystal']:.3f}"
        )

    print("# availability deepice=no-public-code icecoder=no-public-code "
          "zeron2024=no-public-code (searched GitHub and package indexes)")
    return 0


if __name__ == "__main__":
    sys.exit(main())
