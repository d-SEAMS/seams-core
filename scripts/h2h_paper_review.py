#!/usr/bin/env python3
"""Paper-review losses: rings, symmetry, library cages, CHILL+ matrices.

Pinned ABAB, 8 samples. Serial is one core. Eight threads are cores
0-7. The machine core count is printed. No host name.
"""

from __future__ import annotations

import ctypes
import os
import subprocess
import sys
import time
from pathlib import Path

import numpy as np

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "scripts"))

from gen_ice_ic import ice_ic  # noqa: E402
from h2h_freud_ovito import (  # noqa: E402
    CUTOFF,
    SEAMS_CHILL_OVITO,
    SEAMS_GRAPH_CUTOFF,
    core_count,
    frame_system,
    ice_positions,
    library_path,
    load_seams,
    mixed_frames,
    require,
    set_ovito_threads,
)
from h2h_rings import bind_rings, load_sii  # noqa: E402

N = 8
SEAMS_CHILL_PUBLISHED = 0
TO_OVITO = np.array([2, 1, 0, 3, 4, 5, 0, 0, 0], dtype=np.int64)
LABELS = ("cubic", "hexagonal", "water", "interfacial", "clathrate", "interclathrate")


def pin(width: int) -> None:
    try:
        os.sched_setaffinity(0, set(range(width)))
    except OSError:
        pass


def stats(samples):
    arr = np.asarray(samples, dtype=np.float64)
    return float(np.median(arr)), float(arr.min())


def abab(ours, other):
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
    return stats(a), stats(b)


def verdict(a, b) -> str:
    return "WIN" if a[0] < b[0] and a[1] < b[1] else "LOSE"


def confusion(left, right) -> str:
    counts = {}
    for a, b in zip(left, right):
        counts[(a, b)] = counts.get((a, b), 0) + 1
    return " ".join(f"{a}/{b}={counts[a, b]}" for a, b in sorted(counts))


def chill(lib, sys, n, rule) -> np.ndarray:
    labels = (ctypes.c_int * n)()
    counts = (ctypes.c_int * 9)()
    require(
        lib,
        lib.seams_chill_plus(sys, CUTOFF, 1, SEAMS_GRAPH_CUTOFF, rule, labels, counts),
        "chill+",
    )
    return np.ctypeslib.as_array(labels).copy()


def main():
    os.environ.setdefault("QT_QPA_PLATFORM", "offscreen")
    # OpenMP places the team on the mask it first sees. Start on every
    # core, then narrow only for a timed serial call.
    pin(os.cpu_count() or 1)
    lib = load_seams(library_path())
    bind_rings(lib)
    lib.seams_cages_signature.argtypes = [
        ctypes.c_void_p, ctypes.c_double, ctypes.c_int, ctypes.c_int, ctypes.c_char_p,
        ctypes.POINTER(ctypes.c_size_t), ctypes.POINTER(ctypes.c_size_t),
    ]
    lib.seams_cages_signature.restype = ctypes.c_int
    # Loading the library can shrink the mask to one core. Restore it
    # before the header, and again inside each timed call.
    pin(os.cpu_count() or 1)
    print(f"cores {os.cpu_count()} affinity {core_count()}", flush=True)
    print(f"samples {N}", flush=True)
    print(
        "recorded rings OVITO faster at 136, 4096, 110592: 865 ms vs 511 ms",
        flush=True,
    )
    print(
        "recorded symmetry d-SEAMS 40.6 ms vs identify-diamond 36.1 ms",
        flush=True,
    )
    print(
        "recorded library cages 282 ms vs command 213 ms",
        flush=True,
    )

    def rings_pair(pos, lengths, name, width):
        pin(width)
        require(lib, lib.seams_set_num_threads(width), "threads")
        set_ovito_threads(width)
        pos = np.ascontiguousarray(pos, dtype=np.float64)
        held = frame_system(lib, pos, lengths)

        def ours():
            pin(width)
            rings = ctypes.c_void_p()
            require(lib, lib.seams_rings_find(held, CUTOFF, 1, 0, 6, ctypes.byref(rings)), "rings")
            hist = (ctypes.c_int * 7)()
            require(lib, lib.seams_rings_histogram(rings, hist, 7), "histogram")
            lib.seams_rings_destroy(rings)

        data = DataCollection()
        data.create_cell(np.diag(lengths), pbc=(True, True, True))
        parts = data.create_particles(count=len(pos))
        parts.create_property("Position", data=pos)
        pipe = Pipeline(source=StaticSource(data=data))
        pipe.modifiers.append(CreateBondsModifier(cutoff=CUTOFF))
        pipe.modifiers.append(FindRingsModifier(maximum_ring_size=6))

        def other():
            pin(width)
            for modifier in pipe.modifiers:
                modifier.notify_object_changed()
            pipe.compute()

        da, db = abab(ours, other)
        print(
            f"robust op {name} rings threads {width} molecules {len(pos)} "
            f"dseams_median {da[0]:.4f} dseams_min {da[1]:.4f} "
            f"ovito_median {db[0]:.4f} ovito_min {db[1]:.4f} "
            f"affinity {core_count()} {verdict(da, db)}",
            flush=True,
        )
        lib.seams_system_destroy(held)
        return verdict(da, db) == "WIN"

    def sym_pair(pos, lengths, name, width):
        pin(width)
        require(lib, lib.seams_set_num_threads(width), "threads")
        set_ovito_threads(width)
        pos = np.ascontiguousarray(pos, dtype=np.float64)
        n = len(pos)
        held = frame_system(lib, pos, lengths)
        buf = (ctypes.c_int * n)()

        def ours():
            pin(width)
            require(lib, lib.seams_system_symmetry(held, 0.0, buf, n), "symmetry")

        data = DataCollection()
        data.create_cell(np.diag(lengths), pbc=(True, True, True))
        parts = data.create_particles(count=n)
        parts.create_property("Position", data=pos)
        pipe = Pipeline(source=StaticSource(data=data))
        pipe.modifiers.append(IdentifyDiamondModifier())

        def other():
            pin(width)
            pipe.modifiers[0].notify_object_changed()
            pipe.compute()

        # Labels once, before the timer, so the speed path is the same call.
        require(lib, lib.seams_system_symmetry(held, 0.0, buf, n), "symmetry")
        ov = np.asarray(pipe.compute().particles["Structure Type"])
        # identify-diamond: 0 other, 1 cubic, 2 cubic 1st, 3 cubic 2nd,
        # 4 hex, 5 hex 1st, 6 hex 2nd. d-SEAMS symmetry: 0 hex, 1 cubic.
        dia = []
        for value in ov:
            code = int(value)
            if code in (1, 2, 3):
                dia.append("cubic")
            elif code in (4, 5, 6):
                dia.append("hex")
            else:
                dia.append("other")
        sym = []
        for value in np.ctypeslib.as_array(buf):
            code = int(value)
            sym.append({0: "hex", 1: "cubic", 2: "stack"}.get(code, "other"))
        same = sum(a == b for a, b in zip(sym, dia))
        print(
            f"check {name} symmetry molecules {n} agree {same} of {n} "
            + confusion(sym, dia),
            flush=True,
        )
        da, db = abab(ours, other)
        print(
            f"robust op {name} symmetry threads {width} molecules {n} "
            f"dseams_median {da[0]:.4f} dseams_min {da[1]:.4f} "
            f"diamond_median {db[0]:.4f} diamond_min {db[1]:.4f} "
            f"affinity {core_count()} {verdict(da, db)}",
            flush=True,
        )
        lib.seams_system_destroy(held)
        return verdict(da, db) == "WIN"

    wins = True
    # Cages before the large ring census. The ring search keeps a
    # per-source pool; timing cages after that pool is a different question.
    # Library signature cages against the command, orthorhombic sII 8x8x8.
    cell = ROOT / "bbdir-release" / "h2h-paper" / "sII8.lammpstrj"
    if not cell.is_file():
        from h2h_cages import replicate_ortho

        cell.parent.mkdir(parents=True, exist_ok=True)
        replicate_ortho(ROOT / "input" / "traj" / "genice_sII.lammpstrj", 8, 8, 8, cell)
    from h2h_cages import read_dump

    pos, bounds, _tilts = read_dump(cell)
    lengths = tuple(hi - lo for lo, hi in bounds)
    print(f"frame cages molecules {len(pos)}", flush=True)
    for width in (1, 8):
        pin(width)
        require(lib, lib.seams_set_num_threads(width), "threads")
        held = frame_system(lib, pos, lengths)
        nc = ctypes.c_size_t()
        na = ctypes.c_size_t()

        def library():
            pin(width)
            require(
                lib,
                lib.seams_cages_signature(
                    held, 3.5, 1, 0, b"51264", ctypes.byref(nc), ctypes.byref(na)
                ),
                "cages",
            )

        env = os.environ.copy()
        env["LD_LIBRARY_PATH"] = str(ROOT / "bbdir-release" / "src")
        env["OMP_NUM_THREADS"] = str(width)
        cmd = [
            "taskset", "-c", "0" if width == 1 else "0-7",
            str(ROOT / "bbdir-release" / "seams"),
            "--jobs", str(width), "--cutoff", "3.5", "--graph", "cutoff",
            "--signature", "51264", "cages", str(cell),
        ]

        def command():
            subprocess.run(cmd, cwd=ROOT, env=env, check=True, capture_output=True, text=True)

        da, db = abab(library, command)
        print(
            f"robust op sII69632 cages threads {width} molecules {len(pos)} "
            f"library_median {da[0]:.4f} library_min {da[1]:.4f} "
            f"command_median {db[0]:.4f} command_min {db[1]:.4f} "
            f"affinity {core_count()} {verdict(da, db)} "
            f"cages {nc.value}",
            flush=True,
        )
        lib.seams_system_destroy(held)
        if verdict(da, db) != "WIN":
            wins = False

    from ovito.data import DataCollection
    from ovito.modifiers import (
        ChillPlusModifier,
        CreateBondsModifier,
        FindRingsModifier,
        IdentifyDiamondModifier,
    )
    from ovito.pipeline import Pipeline, StaticSource

    import ovito

    print(f"ovito {'.'.join(str(v) for v in ovito.version)}", flush=True)
    frames = [("sII136", *load_sii())]
    pos, lengths = ice_positions(8)
    frames.append(("ice4096", pos, lengths))
    pos, lengths = ice_positions(24)
    frames.append(("ice110592", pos, lengths))
    for width in (1, 8):
        for name, pos, lengths in frames:
            if not rings_pair(pos, lengths, name, width):
                wins = False
            if name != "sII136":
                if not sym_pair(pos, lengths, name, width):
                    wins = False

    # CHILL+ confusion. Thresholds are unchanged.
    print("chill_plus thresholds published staggered<=-0.35 eclipsed -0.8..0.25", flush=True)
    print("chill_plus ovito mode uses the open edges and every neighbour inside 3.5 A", flush=True)
    chill_frames = []
    pos, lengths = ice_positions(4)
    chill_frames.append(("Ih512", pos, lengths))
    pts, side = ice_ic(4)
    chill_frames.append(("Ic512", np.ascontiguousarray(pts), (side, side, side)))
    for name, pair in mixed_frames():
        chill_frames.append((name, pair[0], pair[1]))
    require(lib, lib.seams_set_num_threads(1), "threads")
    for name, pos, lengths in chill_frames:
        pos = np.ascontiguousarray(pos, dtype=np.float64)
        n = len(pos)
        sys = frame_system(lib, pos, lengths)
        paper = chill(lib, sys, n, SEAMS_CHILL_PUBLISHED)
        ovito_rule = chill(lib, sys, n, SEAMS_CHILL_OVITO)
        data = DataCollection()
        data.create_cell(np.diag(lengths), pbc=(True, True, True))
        parts = data.create_particles(count=n)
        parts.create_property("Position", data=pos)
        pipe = Pipeline(source=StaticSource(data=data))
        pipe.modifiers.append(ChillPlusModifier(cutoff=CUTOFF))
        modifier = np.asarray(pipe.compute().particles["Structure Type"])
        mapped = TO_OVITO[ovito_rule]
        miss = int(np.count_nonzero(mapped != modifier))
        paper_name = [LABELS[int(v)] if 0 <= int(v) < len(LABELS) else "other" for v in paper]
        rule_name = [LABELS[int(v)] if 0 <= int(v) < len(LABELS) else "other" for v in ovito_rule]
        mod_name = []
        back = {0: "water", 1: "hexagonal", 2: "cubic", 3: "interfacial", 4: "clathrate", 5: "interclathrate"}
        for value in modifier:
            mod_name.append(back.get(int(value), "other"))
        print(
            f"check {name} chill_ovito_vs_modifier mismatches {miss} of {n}",
            flush=True,
        )
        print(
            f"matrix {name} paper_vs_ovito_rule " + confusion(paper_name, rule_name),
            flush=True,
        )
        print(
            f"matrix {name} ovito_rule_vs_modifier " + confusion(rule_name, mod_name),
            flush=True,
        )
        lib.seams_system_destroy(sys)
        if miss:
            wins = False
    if not wins:
        raise SystemExit(1)


if __name__ == "__main__":
    main()
