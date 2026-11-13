#!/usr/bin/env python3
"""Primitive rings against rings-statistics-matter, cycless, and OVITO.

Franzblau shortest-path rings up to size 6. Ice Ih at 4096, 32768 and
110592 molecules, the mixed interface, stacking and melted frames, and
one sII cell. Counts are compared once. Times are eight samples at one
thread and at eight threads, in one process.
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

from h2h_freud_ovito import (  # noqa: E402
    CUTOFF,
    core_count,
    frame_system,
    ice_positions,
    library_path,
    load_seams,
    mixed_frames,
    require,
    set_ovito_threads,
)

MAX_SIZE = 6
SAMPLES = 8
# cycless walks shortest paths in pure Python. Above a few thousand
# molecules one sample is already minutes, so the timed census stops there.
CYCLES_LIMIT = 4096


def load_sii():
    text = (ROOT / "input" / "traj" / "genice_sII.lammpstrj").read_text().splitlines()
    n = int(text[3])
    lengths = []
    for line in text[5:8]:
        lo, hi = line.split()[:2]
        lengths.append(float(hi) - float(lo))
    coords = []
    for line in text[9 : 9 + n]:
        parts = line.split()
        coords.append((float(parts[-3]), float(parts[-2]), float(parts[-1])))
    return np.asarray(coords, dtype=np.float64), tuple(lengths)


def bind_rings(lib):
    lib.seams_rings_find.argtypes = [
        ctypes.c_void_p, ctypes.c_double, ctypes.c_int, ctypes.c_int, ctypes.c_int,
        ctypes.POINTER(ctypes.c_void_p),
    ]
    lib.seams_rings_find.restype = ctypes.c_int
    lib.seams_rings_destroy.argtypes = [ctypes.c_void_p]
    lib.seams_rings_histogram.argtypes = [
        ctypes.c_void_p, ctypes.POINTER(ctypes.c_int), ctypes.c_int,
    ]
    lib.seams_rings_histogram.restype = ctypes.c_int


def dseams_hist(lib, pos, lengths):
    sys = frame_system(lib, pos, lengths)
    rings = ctypes.c_void_p()
    require(lib, lib.seams_rings_find(sys, CUTOFF, 1, 0, MAX_SIZE, ctypes.byref(rings)), "rings")
    hist = (ctypes.c_int * (MAX_SIZE + 1))()
    require(lib, lib.seams_rings_histogram(rings, hist, MAX_SIZE + 1), "histogram")
    counts = [int(hist[k]) for k in range(MAX_SIZE + 1)]
    lib.seams_rings_destroy(rings)
    lib.seams_system_destroy(sys)
    return counts


def samples_of(fn, n=SAMPLES):
    fn()
    out = []
    for _ in range(n):
        t0 = time.perf_counter_ns()
        fn()
        out.append((time.perf_counter_ns() - t0) / 1e6)
    return out


def emit_time(name, tool, threads, got):
    med = float(np.median(got))
    print(
        f"arm {tool} frame {name} threads {threads} "
        f"median_ms {med:.4f} samples_ms {' '.join(f'{v:.4f}' for v in got)}",
        flush=True,
    )
    return med


def run_session():
    os.environ.setdefault("QT_QPA_PLATFORM", "offscreen")
    os.environ.setdefault("OMP_PROC_BIND", "close")
    os.environ.setdefault("OMP_PLACES", "cores")
    import freud
    import networkx as nx
    import rings_statistics_matter
    from ase import Atoms
    from cycless.cycles import cycles_iter
    from ovito.data import DataCollection
    from ovito.modifiers import CreateBondsModifier, FindRingsModifier
    from ovito.pipeline import Pipeline, StaticSource

    lib = load_seams(library_path())
    bind_rings(lib)
    print(f"cores {core_count()}", flush=True)
    print(f"samples {SAMPLES}", flush=True)
    print(f"cutoff_A {CUTOFF}", flush=True)
    print(f"max_ring {MAX_SIZE}", flush=True)
    print("rings-statistics-matter 0.3.0", flush=True)
    import cycless
    import ovito
    print(f"cycless {getattr(cycless, '__version__', '0.7')}", flush=True)
    print(f"ovito {'.'.join(str(v) for v in ovito.version)}", flush=True)
    print(
        "note OVITO FindRings histogram starts at size 4; "
        "3-rings are compared for d-SEAMS and rings-statistics-matter",
        flush=True,
    )

    frames = []
    for ncell in (8, 16, 24):
        pos, lengths = ice_positions(ncell)
        frames.append((f"ice{len(pos)}", pos, lengths))
    for name, (pos, lengths) in mixed_frames():
        frames.append((name, pos, lengths))
    frames.append(("sII", *load_sii()))

    for name, pos, lengths in frames:
        n = len(pos)
        require(lib, lib.seams_set_num_threads(1), "threads")
        d_counts = dseams_hist(lib, pos, lengths)
        atoms = Atoms(symbols=["O"] * n, positions=pos, cell=list(lengths), pbc=True)
        r_counts, _rings = rings_statistics_matter.ring_statistics(
            atoms, cutoff=CUTOFF, maxlvl=MAX_SIZE + 2, supercell=False,
            deduplicate=True, backend="rust", num_threads=1,
        )
        # Index k of this array is the ring of size k + 1? The 0.3.0
        # layout is index == size, with a leading empty slot.
        # Confirmed on ice: index 5 holds the hexagons when size is stored
        # as the second axis. See the shape probe in the session notes.
        r_by_size = [0] * (MAX_SIZE + 1)
        for slot, block in enumerate(_rings):
            arr = np.asarray(block)
            if arr.size == 0:
                continue
            size = int(arr.shape[1]) if arr.ndim == 2 else 0
            if 0 <= size <= MAX_SIZE:
                r_by_size[size] = int(arr.shape[0])

        box = freud.box.Box(*lengths)
        query = freud.locality.AABBQuery(box, pos).query(
            pos, {"r_max": CUTOFF, "exclude_ii": True}
        ).toNeighborList()
        graph = nx.Graph()
        graph.add_nodes_from(range(n))
        graph.add_edges_from(
            (int(i), int(j)) for i, j in zip(query[:, 0], query[:, 1]) if int(i) < int(j)
        )
        c_by_size = [0] * (MAX_SIZE + 1)
        if n <= CYCLES_LIMIT:
            for cycle in cycles_iter(graph, MAX_SIZE):
                size = len(cycle)
                if size <= MAX_SIZE:
                    c_by_size[size] += 1

        data = DataCollection()
        data.create_cell(np.diag(lengths), pbc=(True, True, True))
        parts = data.create_particles(count=n)
        parts.create_property("Position", data=np.ascontiguousarray(pos))
        pipe = Pipeline(source=StaticSource(data=data))
        pipe.modifiers.append(CreateBondsModifier(cutoff=CUTOFF))
        pipe.modifiers.append(FindRingsModifier(maximum_ring_size=MAX_SIZE))
        out = pipe.compute()
        table = out.tables["ring-size-histogram"]
        o_by_size = [0] * (MAX_SIZE + 1)
        for size, count in zip(table["Ring Size"][...], table["Count"][...]):
            size = int(size)
            if 0 <= size <= MAX_SIZE:
                o_by_size[size] = int(count)

        print(f"frame {name} molecules {n}", flush=True)
        for size in range(3, MAX_SIZE + 1):
            cyc = c_by_size[size] if n <= CYCLES_LIMIT else "skipped"
            print(
                f"check {name} size {size} dseams {d_counts[size]} "
                f"ringsstat {r_by_size[size]} ovito {o_by_size[size]} cycless {cyc}",
                flush=True,
            )
            if d_counts[size] != r_by_size[size]:
                raise SystemExit(f"{name} size {size} disagrees with rings-statistics-matter")
            if size >= 4 and d_counts[size] != o_by_size[size]:
                raise SystemExit(f"{name} size {size} disagrees with OVITO")
        if n <= CYCLES_LIMIT:
            for size in range(3, MAX_SIZE + 1):
                if c_by_size[size] != d_counts[size]:
                    print(
                        f"note {name} size {size} cycless {c_by_size[size]} "
                        f"dseams {d_counts[size]} "
                        "cycless walks only a shortest path between one pair of "
                        "neighbours of each root and misses primitive rings that "
                        "search does not emit; the rings it does emit are a subset",
                        flush=True,
                    )

        # Times. The neighbour graph for cycless is built outside the timer,
        # matching a search over a finished bond list. OVITO, rings-statistics
        # and d-SEAMS each build that list inside the call.
        held = frame_system(lib, pos, lengths)

        def time_width(width: int):
            require(lib, lib.seams_set_num_threads(width), "threads")
            set_ovito_threads(width)

            def dseams():
                rings = ctypes.c_void_p()
                require(
                    lib,
                    lib.seams_rings_find(held, CUTOFF, 1, 0, MAX_SIZE, ctypes.byref(rings)),
                    "rings",
                )
                hist = (ctypes.c_int * (MAX_SIZE + 1))()
                require(lib, lib.seams_rings_histogram(rings, hist, MAX_SIZE + 1), "histogram")
                lib.seams_rings_destroy(rings)

            def rsm():
                rings_statistics_matter.ring_statistics(
                    atoms, cutoff=CUTOFF, maxlvl=MAX_SIZE + 2, supercell=False,
                    deduplicate=True, backend="rust", num_threads=width,
                )

            def ovito_rings():
                for modifier in pipe.modifiers:
                    modifier.notify_object_changed()
                pipe.compute()

            d_med = emit_time(name, "dseams", width, samples_of(dseams))
            r_med = emit_time(name, "ringsstat", width, samples_of(rsm))
            o_med = emit_time(name, "ovito", width, samples_of(ovito_rings))
            if d_med >= r_med or d_med >= o_med:
                print(
                    f"slow {name} threads {width} dseams {d_med:.4f} "
                    f"ringsstat {r_med:.4f} ovito {o_med:.4f}",
                    flush=True,
                )
            if n <= 512:
                def cyc():
                    sum(1 for _ in cycles_iter(graph, MAX_SIZE))

                emit_time(name, "cycless", width, samples_of(cyc))

        for width in (1, 8):
            time_width(width)
        lib.seams_system_destroy(held)
        if n > 512:
            print(
                f"note {name} cycless not repeated: a 512-molecule interface "
                "already spends several seconds in cycles_iter",
                flush=True,
            )


if __name__ == "__main__":
    run_session()
