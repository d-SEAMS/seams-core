#!/usr/bin/env python3
"""OVITO polyhedral template matching against d-SEAMS symmetry.

No remote branch carries bc-035d2d78's PTM/IRA tree. This session is
OVITO 3.17.2 PTM, with the cubic and hexagonal diamond templates
enabled, against seams_system_symmetry. Ideal ice Ih and ice Ic, the
stacking frame, and the mixed interface and melted frames.
"""

from __future__ import annotations

import ctypes
import os
import sys
import time
from pathlib import Path

import numpy as np

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "scripts"))

from gen_ice_ic import ice_ic  # noqa: E402
from h2h_freud_ovito import (  # noqa: E402
    PTM,
    SYM,
    core_count,
    frame_system,
    ice_positions,
    library_path,
    load_seams,
    mixed_frames,
    require,
    set_ovito_threads,
)

N = 20


def pin(width: int) -> None:
    try:
        os.sched_setaffinity(0, set(range(width if width > 0 else 1)))
    except OSError:
        pass


def stats(samples):
    arr = np.asarray(samples, dtype=np.float64)
    q1, med, q3 = np.percentile(arr, [25, 50, 75])
    return float(med), float(q3 - q1), float(arr.min())


def confusion(left, right) -> str:
    counts = {}
    for a, b in zip(left, right):
        counts[(a, b)] = counts.get((a, b), 0) + 1
    return " ".join(f"{a}/{b}={counts[a, b]}" for a, b in sorted(counts))


def dseams_sym(lib, pos, lengths):
    sys = frame_system(lib, pos, lengths)
    n = len(pos)
    labels = (ctypes.c_int * n)()
    require(lib, lib.seams_system_symmetry(sys, 3.5, labels, n), "symmetry")
    out = np.ctypeslib.as_array(labels).copy()
    lib.seams_system_destroy(sys)
    return out


def main():
    width = int(sys.argv[1]) if len(sys.argv) > 1 else 1
    os.environ.setdefault("QT_QPA_PLATFORM", "offscreen")
    os.environ["OVITO_THREAD_COUNT"] = str(width)
    pin(width)
    from ovito.data import DataCollection
    from ovito.modifiers import PolyhedralTemplateMatchingModifier
    from ovito.pipeline import Pipeline, StaticSource

    pin(width)
    set_ovito_threads(width)
    lib = load_seams(library_path())
    require(lib, lib.seams_set_num_threads(width), "threads")
    pin(width)
    print(f"cores {os.cpu_count()} affinity {core_count()}", flush=True)
    print(f"pinned_threads {width}", flush=True)
    print(f"samples {N}", flush=True)
    print(
        "note no remote branch holds bc-035d2d78 PTM/IRA work; "
        "the standalone PTM library is not installed, OVITO PTM is the reference",
        flush=True,
    )
    frames = []
    for ncell, label in ((8, "Ih4096"), (16, "Ih32768")):
        pos, lengths = ice_positions(ncell)
        frames.append((label, pos, lengths))
    for rep, label in ((8, "Ic4096"), (16, "Ic32768")):
        pts, side = ice_ic(rep)
        pos = np.ascontiguousarray(pts, dtype=np.float64)
        frames.append((label, pos, (side, side, side)))
    for name, (pos, lengths) in mixed_frames():
        frames.append((name, pos, lengths))

    wins = True
    for name, pos, lengths in frames:
        pos = np.ascontiguousarray(pos, dtype=np.float64)
        n = len(pos)
        d_sym = dseams_sym(lib, pos, lengths)
        data = DataCollection()
        data.create_cell(np.diag(lengths), pbc=(True, True, True))
        parts = data.create_particles(count=n)
        parts.create_property("Position", data=pos)
        ptm = PolyhedralTemplateMatchingModifier(rmsd_cutoff=0.1)
        ptm.structures[PolyhedralTemplateMatchingModifier.Type.CUBIC_DIAMOND].enabled = True
        ptm.structures[PolyhedralTemplateMatchingModifier.Type.HEX_DIAMOND].enabled = True
        pipe = Pipeline(source=StaticSource(data=data))
        pipe.modifiers.append(ptm)
        ov = np.asarray(pipe.compute().particles["Structure Type"])
        sym_name = [SYM.get(int(v), "other") for v in d_sym]
        ptm_name = [PTM.get(int(v), "other") for v in ov]
        same = sum(a == b for a, b in zip(sym_name, ptm_name))
        print(
            f"check {name} molecules {n} agree {same} of {n} "
            + confusion(sym_name, ptm_name),
            flush=True,
        )

        def ours():
            pin(width)
            dseams_sym(lib, pos, lengths)

        def other():
            pin(width)
            ptm.notify_object_changed()
            pipe.compute()

        ours(); other()
        a, b = [], []
        for _ in range(N):
            t0 = time.perf_counter_ns(); ours(); a.append((time.perf_counter_ns() - t0) / 1e6)
            t0 = time.perf_counter_ns(); other(); b.append((time.perf_counter_ns() - t0) / 1e6)
        am, ai, amin = stats(a)
        bm, bi, bmin = stats(b)
        win = am < bm and amin < bmin
        print(
            f"robust op {name} symmetry threads {width} molecules {n} "
            f"dseams_median {am:.4f} dseams_iqr {ai:.4f} dseams_min {amin:.4f} "
            f"ovito_ptm_median {bm:.4f} ovito_ptm_iqr {bi:.4f} ovito_ptm_min {bmin:.4f} "
            f"{'WIN' if win else 'LOSE'}",
            flush=True,
        )
        if not win:
            wins = False
    if not wins:
        raise SystemExit(1)


if __name__ == "__main__":
    main()
