#!/usr/bin/env python3
"""End-to-end neighbours plus rings, and disordered-frame counts.

The timer covers the neighbour build and the ring search. Ice Ic is the
cell in scripts/gen_ice_ic.py. SQQ 0.5.8 counts chordless rings. OVITO
counts the rings FindRingsModifier reports after CreateBondsModifier.
"""

from __future__ import annotations

import os
import sys
import time

import numpy as np

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from gen_ice_ic import ice_ic

CUTOFF = 3.5


def cell_graph(points, side, cutoff=CUTOFF):
    n = len(points)
    ncell = max(1, int(side / cutoff))
    cell = side / ncell
    buckets = [[] for _ in range(ncell ** 3)]

    def cell_of(x):
        c = int(x / cell)
        if c < 0:
            c = 0
        if c >= ncell:
            c = ncell - 1
        return c

    for i, p in enumerate(points):
        buckets[(cell_of(p[0]) * ncell + cell_of(p[1])) * ncell + cell_of(p[2])].append(i)
    adj = [[] for _ in range(n)]
    cut2 = cutoff * cutoff
    for i, pi in enumerate(points):
        cx, cy, cz = cell_of(pi[0]), cell_of(pi[1]), cell_of(pi[2])
        for dx in (-1, 0, 1):
            for dy in (-1, 0, 1):
                for dz in (-1, 0, 1):
                    nx = (cx + dx + ncell) % ncell
                    ny = (cy + dy + ncell) % ncell
                    nz = (cz + dz + ncell) % ncell
                    for j in buckets[(nx * ncell + ny) * ncell + nz]:
                        if j <= i:
                            continue
                        d2 = 0.0
                        for a in range(3):
                            d = pi[a] - points[j][a]
                            if d > side / 2:
                                d -= side
                            elif d < -side / 2:
                                d += side
                            d2 += d * d
                        if 0.0 < d2 < cut2:
                            adj[i].append(j)
                            adj[j].append(i)
    return adj


def sqq_counts(adj, sizes):
    from sqq.core.sqq_py.ring import find_rings

    adjacency = {i: set(row) for i, row in enumerate(adj)}
    found = find_rings(adjacency, sizes, chordless=True, definition="chordless")
    by = {s: len(found.get(s, [])) for s in sizes}
    return by, sum(by.values())


def ovito_counts(points, side, maxsize):
    os.environ["QT_QPA_PLATFORM"] = "offscreen"
    from ovito.data import DataCollection
    from ovito.modifiers import CreateBondsModifier, FindRingsModifier
    from ovito.pipeline import Pipeline, StaticSource

    data = DataCollection()
    data.create_cell(np.diag([side, side, side]), pbc=(True, True, True))
    parts = data.create_particles(count=len(points))
    parts.create_property("Position", data=np.asarray(points, dtype=np.float64))
    pipe = Pipeline(source=StaticSource(data=data))
    pipe.modifiers.append(CreateBondsModifier(cutoff=CUTOFF))
    pipe.modifiers.append(FindRingsModifier(maximum_ring_size=maxsize))
    out = pipe.compute()
    hist = out.tables["ring-size-histogram"]
    sizes = hist["Ring Size"][...]
    counts = hist["Count"][...]
    by = {int(s): int(c) for s, c in zip(sizes, counts)}
    return by, int(np.sum(counts))


def time_once(fn):
    fn()
    t0 = time.perf_counter()
    result = fn()
    return (time.perf_counter() - t0) * 1e3, result


def ice_session(threads):
    os.environ["OVITO_THREAD_COUNT"] = str(threads)
    os.environ["OMP_NUM_THREADS"] = str(threads)
    for rep in (11, 24, 48):
        points, side = ice_ic(rep)
        print(f"foreign molecules {len(points)} threads {threads} building", flush=True)

        def sqq():
            adj = cell_graph(points, side)
            return sqq_counts(adj, [6])

        ms, (by, total) = time_once(sqq)
        print(
            f"sample tool sqq-0.5.8 threads {threads} molecules {len(points)} "
            f"depth 6 ms {ms:.6f} rings {total} n6 {by.get(6, 0)}",
            flush=True,
        )

        def ovi():
            return ovito_counts(points, side, 6)

        ms, (by, total) = time_once(ovi)
        print(
            f"sample tool ovito threads {threads} molecules {len(points)} "
            f"depth 6 ms {ms:.6f} rings {total} n6 {by.get(6, 0)}",
            flush=True,
        )


def main():
    print("session rings v3 foreign", flush=True)
    for threads in (1, 8):
        ice_session(threads)


if __name__ == "__main__":
    main()
