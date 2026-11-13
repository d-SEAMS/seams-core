#!/usr/bin/env python3
"""One session: d-SEAMS, SQQ 0.5.8, OVITO FindRings, and cycless.

Ice Ic is built in memory (the same cell as scripts/gen_ice_ic.py). The
neighbour graph is built before each foreign timer. d-SEAMS times only
ringNetwork, via tests/bench_ring_ic. SQQ is the 0.5.8 chordless search.
OVITO times FindRingsModifier after the bonds exist. cycless times
cycles_iter after the networkx graph exists.
"""

from __future__ import annotations

import os
import subprocess
import sys
import time

import numpy as np

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from gen_ice_ic import ice_ic

CUTOFF = 3.5
ROOT = os.path.abspath(os.path.join(os.path.dirname(__file__), ".."))
BENCH = os.path.join(ROOT, "build-rel", "tests", "bench_ring_ic")


def cell_graph(points, side, cutoff=CUTOFF):
    """Undirected neighbour lists, vertex 0..n-1, no self."""
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
                        pj = points[j]
                        r2 = 0.0
                        for a in range(3):
                            d = pi[a] - pj[a]
                            if d > side / 2:
                                d -= side
                            elif d < -side / 2:
                                d += side
                            r2 += d * d
                        if 1e-12 < r2 <= cut2:
                            adj[i].append(j)
                            adj[j].append(i)
    return adj


def time_calls(fn, repeats, warmup=1):
    samples = []
    last = None
    for i in range(warmup + repeats):
        t0 = time.perf_counter()
        last = fn()
        dt = (time.perf_counter() - t0) * 1e3
        if i >= warmup:
            samples.append(dt)
    return samples, last


def run_dseams(repetitions, depth, threads_note="1 8"):
    env = os.environ.copy()
    env["LD_LIBRARY_PATH"] = os.path.join(ROOT, "build-rel", "src") + os.pathsep + env.get(
        "LD_LIBRARY_PATH", ""
    )
    env["OMP_WAIT_POLICY"] = "PASSIVE"
    cmd = [BENCH, "--depth", str(depth), str(repetitions)]
    print(f"tool dseams depth {depth} repetitions {repetitions}", flush=True)
    proc = subprocess.run(cmd, env=env, check=True, capture_output=True, text=True)
    sys.stdout.write(proc.stdout)
    if proc.stderr:
        sys.stderr.write(proc.stderr)


def sqq_count(adj, sizes):
    from sqq.core.sqq_py.ring import find_rings

    adjacency = {i: set(row) for i, row in enumerate(adj)}
    found = find_rings(adjacency, sizes, chordless=True, definition="chordless")
    return sum(len(v) for v in found.values())


def cycless_count(adj, maxsize):
    import networkx as nx
    from cycless.cycles import cycles_iter

    graph = nx.Graph()
    graph.add_nodes_from(range(len(adj)))
    graph.add_edges_from((i, j) for i, row in enumerate(adj) for j in row if i < j)
    return sum(1 for _ in cycles_iter(graph, maxsize))


def ovito_prepare(points, side):
    """Particles and cutoff bonds. The ring timer starts after this returns.

    libGL lives next to a libc this process must not prefer, so the caller
    starts OVITO with that directory already on the library path.
    """
    os.environ["QT_QPA_PLATFORM"] = "offscreen"
    os.environ["OVITO_THREAD_COUNT"] = "1"
    from ovito.data import DataCollection
    from ovito.modifiers import CreateBondsModifier
    from ovito.pipeline import Pipeline, StaticSource

    data = DataCollection()
    data.create_cell(np.diag([side, side, side]), pbc=(True, True, True))
    parts = data.create_particles(count=len(points))
    parts.create_property("Position", data=np.asarray(points, dtype=np.float64))
    bond_pipe = Pipeline(source=StaticSource(data=data))
    bond_pipe.modifiers.append(CreateBondsModifier(cutoff=CUTOFF))
    return bond_pipe.compute()


def ovito_rings(bonded, maxsize):
    from ovito.modifiers import FindRingsModifier
    from ovito.pipeline import Pipeline, StaticSource

    pipe = Pipeline(source=StaticSource(data=bonded))
    pipe.modifiers.append(FindRingsModifier(maximum_ring_size=maxsize))
    out = pipe.compute()
    hist = out.tables["ring-size-histogram"]
    return int(np.sum(hist["Count"][...]))


def main():
    foreign_only = "--foreign-only" in sys.argv
    print("session rings v2", flush=True)
    if not foreign_only:
        # d-SEAMS at 6, 8 and 10, both sizes, 1 and 8 threads.
        for depth in (6, 8, 10):
            for rep in (11, 24):
                run_dseams(rep, depth)

    for rep in (11, 24):
        points, side = ice_ic(rep)
        print(f"foreign molecules {len(points)} building graph", flush=True)
        adj = cell_graph(points, side)
        edges = sum(len(r) for r in adj) // 2
        print(f"foreign molecules {len(points)} edges {edges}", flush=True)

        def sqq():
            return sqq_count(adj, [6])

        samples, count = time_calls(sqq, repeats=2 if rep == 11 else 1, warmup=1)
        for s in samples:
            print(f"sample tool sqq-0.5.8 molecules {len(points)} ms {s:.6f} rings {count}", flush=True)

        def cyc():
            return cycless_count(adj, 6)

        samples, count = time_calls(cyc, repeats=2 if rep == 11 else 1, warmup=1)
        for s in samples:
            print(f"sample tool cycless molecules {len(points)} ms {s:.6f} rings {count}", flush=True)

        bonded = ovito_prepare(points, side)

        def ovi():
            return ovito_rings(bonded, 6)

        samples, count = time_calls(ovi, repeats=2 if rep == 11 else 1, warmup=1)
        for s in samples:
            print(f"sample tool ovito molecules {len(points)} ms {s:.6f} rings {count}", flush=True)


if __name__ == "__main__":
    main()
