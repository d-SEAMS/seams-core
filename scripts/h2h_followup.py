#!/usr/bin/env python3
"""Recorded-loss rerun, carried rings, and an independent CHILL+ reference.

One session. Serial calls use core 0. Eight-thread calls use cores 0-7.
Eight alternating samples. The core count is printed. No host name.

The ring loss 865 ms vs 511 ms is the type-7 median in
paper/tex/figures/data/ring-head2head.json (git c5da2d9f), one thread,
fresh bench_stages process, the neighbour list built before the timer.
The frames are genice_sH (136), mW cubic (4096) and that cell replicated
3x3x3 (110592). This script times that boundary again.

The 40.6 ms vs 36.1 ms symmetry pair is not a labelled median in the
tree. The nearest identify-diamond record is benchmarks/logs/symmetry-v20-session.txt
on ideal ice of 110592 molecules: the d-SEAMS kernel is 21.743 ms and
identify-diamond is 36.564 ms. The same three frames are timed here.
"""

from __future__ import annotations

import os
import subprocess
import sys
import time
from pathlib import Path

import numpy as np

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "scripts"))

from h2h_freud_ovito import (  # noqa: E402
    CUTOFF,
    SEAMS_CHILL_OVITO,
    SEAMS_GRAPH_CUTOFF,
    core_count,
    frame_system,
    ice_positions,
    library_path,
    load_seams,
    require,
    set_ovito_threads,
)

N = 8
SEAMS_CHILL_PUBLISHED = 0
LABELS = ("cubic", "hexagonal", "water", "interfacial", "clathrate", "interclathrate")


def pin(width: int) -> None:
    try:
        os.sched_setaffinity(0, set(range(width)))
    except OSError:
        pass


def stats(samples):
    arr = np.asarray(samples, dtype=np.float64)
    return float(np.median(arr)), float(arr.min())


def verdict(a, b) -> str:
    return "WIN" if a[0] < b[0] and a[1] < b[1] else "LOSE"


def confusion(left, right) -> str:
    counts = {}
    for a, b in zip(left, right):
        counts[(a, b)] = counts.get((a, b), 0) + 1
    return " ".join(f"{a}/{b}={counts[a, b]}" for a, b in sorted(counts))


def recall_ice(truth_ice, predicted) -> tuple[int, int, int, int, float]:
    """Ice is cubic or hexagonal. Truth is the constructed lattice."""
    tp = fp = fn = tn = 0
    for ice, label in zip(truth_ice, predicted):
        called = label in ("cubic", "hexagonal")
        if ice and called:
            tp += 1
        elif ice and not called:
            fn += 1
        elif called:
            fp += 1
        else:
            tn += 1
    denom = tp + fn
    return tp, fp, fn, tn, (tp / denom if denom else float("nan"))


def chill_labels(lib, sys, n, rule):
    import ctypes

    labels = (ctypes.c_int * n)()
    counts = (ctypes.c_int * 9)()
    require(
        lib,
        lib.seams_chill_plus(sys, CUTOFF, 1, SEAMS_GRAPH_CUTOFF, rule, labels, counts),
        "chill+",
    )
    return [LABELS[int(v)] if 0 <= int(v) < len(LABELS) else "other" for v in labels]


def liquid_shell(seed, lengths, count, gap):
    """Random points outside the seed box, in a longer cell."""
    rng = np.random.default_rng(252)
    lx, ly, lz = lengths
    box = (lx, ly, lz + gap + 12.0)
    pts = []
    while len(pts) < count:
        trial = rng.random(3) * box
        if trial[2] <= lz + gap:
            continue
        pts.append(trial)
    return np.vstack([seed, np.asarray(pts)]), box


def main() -> None:
    os.environ.setdefault("QT_QPA_PLATFORM", "offscreen")
    pin(os.cpu_count() or 1)
    print(f"cores {os.cpu_count()} affinity {core_count()}", flush=True)
    print(f"samples {N}", flush=True)
    print(
        "recorded rings source paper/tex/figures/data/ring-head2head.json "
        "git c5da2d9f one thread fresh bench_stages graph outside the timer "
        "type-7 median orth-1e5 864.701 vs 511.024 "
        "orth-4k 31.617 vs 15.831 tri-136 1.078 vs 0.697",
        flush=True,
    )
    print(
        "recorded symmetry pair 40.6 vs 36.1 is not a labelled median; "
        "nearest is symmetry-v20 ideal-110592 kernel 21.743 "
        "identify-diamond 36.564",
        flush=True,
    )

    bench = ROOT / "bbdir-release" / "tests" / "bench_rings_h2h"
    if not bench.is_file():
        raise SystemExit(f"missing {bench.name}")
    ring_json = ROOT / "benchmarks" / "h2h-recorded-rings.json"
    subprocess.run(
        [sys.executable, "scripts/rings_one_thread.py",
         "--bench", "bbdir-release/tests/bench_rings_h2h",
         "--python", sys.executable,
         "--json", "benchmarks/h2h-recorded-rings.json"],
        cwd=ROOT,
        check=True,
    )
    import json
    record = json.loads(ring_json.read_text())
    wins = True
    for threads, arm in record["threads"].items():
        for name, row in arm.items():
            da = (row["dseams_median_ms"], min(row["dseams_ms"]))
            db = (row["ovito_median_ms"], min(row["ovito_ms"]))
            tag = verdict(da, db)
            if tag != "WIN":
                wins = False
            print(
                f"robust op {name} rings threads {threads} molecules {row['molecules']} "
                f"dseams_median {da[0]:.4f} dseams_min {da[1]:.4f} "
                f"ovito_median {db[0]:.4f} ovito_min {db[1]:.4f} "
                f"affinity {threads} {tag} census {row['census']}",
                flush=True,
            )

    # Symmetry on the same three frames. The neighbour build is inside
    # both calls. Ideal cubic ice is the 4096 and 110592 frames.
    from ovito.data import DataCollection
    from ovito.modifiers import ChillPlusModifier, IdentifyDiamondModifier
    from ovito.pipeline import Pipeline, StaticSource

    import ovito

    print(f"ovito {'.'.join(str(v) for v in ovito.version)}", flush=True)
    pin(os.cpu_count() or 1)
    lib = load_seams(library_path())
    sys.path.insert(0, str(ROOT / "tests"))
    import bench_h2h_rings as h2h

    def lengths_of(frame):
        return (float(frame.H[0, 0]), float(frame.H[1, 1]), float(frame.H[2, 2]))

    cubic = h2h.read_lammps_frame(str(ROOT / "input" / "traj" / "mW_cubic.lammpstrj"), 1)
    replica = h2h.replicate_frame(cubic, 3, 3, 3)
    frames = [
        ("orth-4k", cubic.positions, lengths_of(cubic)),
        ("orth-1e5", replica.positions, lengths_of(replica)),
    ]
    import ctypes

    lib.seams_system_symmetry.argtypes = [
        ctypes.c_void_p, ctypes.c_double, ctypes.POINTER(ctypes.c_int), ctypes.c_int,
    ]
    lib.seams_system_symmetry.restype = ctypes.c_int
    for width in (1, 8):
        pin(width)
        require(lib, lib.seams_set_num_threads(width), "threads")
        set_ovito_threads(width)
        for name, pos, lengths in frames:
            pos = np.ascontiguousarray(pos, dtype=np.float64)
            n = len(pos)
            held = frame_system(lib, pos, lengths)
            buf = (ctypes.c_int * n)()

            def ours(held=held, buf=buf, n=n, width=width):
                pin(width)
                require(lib, lib.seams_system_symmetry(held, 0.0, buf, n), "symmetry")

            data = DataCollection()
            data.create_cell(np.diag(np.asarray(lengths, dtype=np.float64)), pbc=(True, True, True))
            parts = data.create_particles(count=n)
            parts.create_property("Position", data=pos)
            pipe = Pipeline(source=StaticSource(data=data))
            pipe.modifiers.append(IdentifyDiamondModifier())

            def other(pipe=pipe, width=width):
                pin(width)
                pipe.modifiers[0].notify_object_changed()
                pipe.compute()

            ours()
            other()
            a, b = [], []
            for _ in range(N):
                t0 = time.perf_counter_ns()
                ours()
                a.append((time.perf_counter_ns() - t0) / 1e6)
                t0 = time.perf_counter_ns()
                other()
                b.append((time.perf_counter_ns() - t0) / 1e6)
            da, db = stats(a), stats(b)
            tag = verdict(da, db)
            if tag != "WIN":
                wins = False
            print(
                f"robust op {name} symmetry threads {width} molecules {n} "
                f"dseams_median {da[0]:.4f} dseams_min {da[1]:.4f} "
                f"diamond_median {db[0]:.4f} diamond_min {db[1]:.4f} "
                f"affinity {core_count()} {tag}",
                flush=True,
            )
            lib.seams_system_destroy(held)

    carry = ROOT / "bbdir-release" / "tests" / "bench_ring_carry"
    for width in (1, 8):
        env = os.environ.copy()
        # The system lib directory pulls in a libc this binary cannot use.
        env["LD_LIBRARY_PATH"] = str(ROOT / "bbdir-release" / "src")
        env["OMP_NUM_THREADS"] = str(width)
        env["OMP_PROC_BIND"] = "close"
        env["OMP_PLACES"] = "cores"
        env["OMP_WAIT_POLICY"] = "PASSIVE"
        for kind in ("churn", "local"):
            cmd = ["taskset", "-c", "0" if width == 1 else "0-7", str(carry), "4096", "8"]
            if kind == "local":
                cmd.append("local")
            proc = subprocess.run(
                cmd, cwd=ROOT, env=env, check=True, capture_output=True, text=True,
            )
            print(f"ring_carry {kind} threads {width} affinity {width}", flush=True)
            sys.stdout.write(proc.stdout)
            if "identical no" in proc.stdout:
                wins = False
            if kind == "local" and "LOSE" in proc.stdout:
                wins = False

    # Independent ice reference: a constructed Ih lattice, then the same
    # lattice with a liquid slab outside it. The lattice atoms are ice
    # because they were placed on the ice lattice, not because CHILL+ or
    # PTM said so. Published thresholds only.
    print(
        "chill_plus published staggered<=-0.35 eclipsed -0.8..0.25 four nearest",
        flush=True,
    )
    print("chill_plus ovito open edges every neighbour inside 3.5 A", flush=True)
    built = []
    for name, bond in (("mW", 2.76), ("TIP4P", 2.75)):
        # ice_positions is the mW Ih lattice at 2.76. TIP4P uses the same
        # generator scaled so the first-neighbour distance is 2.75.
        pos, lengths = ice_positions(4)
        if bond != 2.76:
            pos = pos * (bond / 2.76)
            lengths = tuple(v * (bond / 2.76) for v in lengths)
        pure = np.ascontiguousarray(pos, dtype=np.float64)
        mixed, box = liquid_shell(pure, lengths, 256, 2.0)
        built.append((f"{name}-lattice", pure, lengths, np.ones(len(pure), dtype=bool)))
        truth = np.zeros(len(mixed), dtype=bool)
        truth[: len(pure)] = True
        built.append((f"{name}-seed", mixed, box, truth))

    for name, pos, lengths, truth in built:
        held = frame_system(lib, pos, lengths)
        paper = chill_labels(lib, held, len(pos), SEAMS_CHILL_PUBLISHED)
        ovito_rule = chill_labels(lib, held, len(pos), SEAMS_CHILL_OVITO)
        data = DataCollection()
        data.create_cell(np.diag(np.asarray(lengths, dtype=np.float64)), pbc=(True, True, True))
        parts = data.create_particles(count=len(pos))
        parts.create_property("Position", data=np.ascontiguousarray(pos, dtype=np.float64))
        pipe = Pipeline(source=StaticSource(data=data))
        pipe.modifiers.append(ChillPlusModifier(cutoff=CUTOFF))
        modifier = np.asarray(pipe.compute().particles["Structure Type"])
        # OVITO ChillPlus: 0 other, 1 hex, 2 cubic, 3 interfacial, 4 clathrate, 5 interclathrate
        ov_names = []
        for value in modifier:
            code = int(value)
            ov_names.append(
                {0: "water", 1: "hexagonal", 2: "cubic", 3: "interfacial",
                 4: "clathrate", 5: "interclathrate"}.get(code, "other")
            )
        truth_name = ["ice" if flag else "other" for flag in truth]
        def ice_name(labels):
            return ["ice" if lab in ("cubic", "hexagonal") else "other" for lab in labels]
        for kind, labels in (
            ("paper", paper),
            ("ovito_rule", ovito_rule),
            ("modifier", ov_names),
        ):
            tp, fp, fn, tn, rec = recall_ice(truth, labels)
            print(
                f"recall {name} {kind} tp {tp} fp {fp} fn {fn} tn {tn} "
                f"recall {rec:.6f} matrix {confusion(truth_name, ice_name(labels))}",
                flush=True,
            )
        print(
            f"matrix {name} paper_vs_ovito_rule {confusion(paper, ovito_rule)}",
            flush=True,
        )
        lib.seams_system_destroy(held)

    sqq = subprocess.run(
        [sys.executable, str(ROOT / "scripts" / "h2h_sqq_box.py")],
        cwd=ROOT,
        check=False,
        capture_output=True,
        text=True,
    )
    sys.stdout.write(sqq.stdout)
    if sqq.returncode != 0:
        sys.stderr.write(sqq.stderr)
        wins = False
    if wins:
        return
    raise SystemExit(1)


if __name__ == "__main__":
    main()
