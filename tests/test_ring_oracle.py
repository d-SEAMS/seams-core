# Copyright (c) d-SEAMS developers
"""PydSEAMSlib ring entries against a brute-force primitive ring oracle.

The oracle lists every simple cycle up to the depth from its lowest
vertex and keeps the ones no pair of members can shortcut through the
graph (Franzblau's shortest-path rings). It is plain Python and shares
nothing with the engine. The checks cover pydseams.yoda.ringNetwork,
RingUpdater and cageAffiliation on index rows, and dseams.neighbours,
dseams.rings, dseams.cages and dseams.signature_cages on coordinates.

Run with the source root as the only argument.
"""

import math
import random
import re
import sys
from collections import deque
from pathlib import Path

import numpy as np

import dseams
import pydseams.yoda as yoda


def simple_graph(rows, n):
    adj = [set() for _ in range(n)]
    for i, row in enumerate(rows):
        for j in row:
            if 0 <= j < n and j != i:
                adj[i].add(j)
                adj[j].add(i)
    return [sorted(a) for a in adj]


def simple_cycles(adj, depth):
    out = []
    for s in range(len(adj)):
        path = [s]
        on = {s}
        stack = [iter(adj[s])]
        while stack:
            v = next(stack[-1], None)
            if v is None:
                stack.pop()
                on.discard(path.pop())
                continue
            if v == s and len(path) >= 3 and path[1] < path[-1]:
                out.append(list(path))
            elif v > s and v not in on and len(path) < depth:
                path.append(v)
                on.add(v)
                stack.append(iter(adj[v]))
    return sorted(out)


def ball(adj, src, radius):
    dist = {src: 0}
    queue = deque([src])
    while queue:
        u = queue.popleft()
        if dist[u] == radius:
            continue
        for v in adj[u]:
            if v not in dist:
                dist[v] = dist[u] + 1
                queue.append(v)
    return dist


def primitive_rings(adj, depth):
    radius = depth // 2
    balls = {}
    keep = []
    for cyc in simple_cycles(adj, depth):
        k = len(cyc)
        ok = True
        for a in range(k):
            ba = balls.get(cyc[a])
            if ba is None:
                ba = balls[cyc[a]] = ball(adj, cyc[a], radius)
            for b in range(a + 1, k):
                around = min(b - a, k - (b - a))
                if ba.get(cyc[b], radius + 1) < around:
                    ok = False
                    break
            if not ok:
                break
        if ok:
            keep.append(cyc)
    return keep


def rows_of(adj, rng=None):
    rows = []
    for i, nbrs in enumerate(adj):
        row = [i] + list(nbrs)
        if rng is not None:
            if nbrs and rng.random() < 0.3:
                row.append(rng.choice(nbrs))
            tail = row[1:]
            rng.shuffle(tail)
            row = [i] + tail
        rows.append(row)
    return rows


def lattice(box9):
    return np.array([[box9[3], 0, 0], [box9[6], box9[4], 0], [box9[7], box9[8], box9[5]]])


def cutoff_graph(xyz, box9, cutoff):
    h = lattice(box9)
    shifts = np.array([a * h[0] + b * h[1] + c * h[2]
                       for a in (-2, -1, 0, 1, 2) for b in (-2, -1, 0, 1, 2)
                       for c in (-2, -1, 0, 1, 2)])
    n = len(xyz)
    rows = [[] for _ in range(n)]
    gap = math.inf
    for i in range(n):
        d = xyz[None, :, :] - xyz[i][None, None, :] + shifts[:, None, :]
        r = np.sqrt((d * d).sum(axis=2)).min(axis=0)
        gap = min(gap, float(np.abs(r[np.arange(n) != i] - cutoff).min()) if n > 1 else gap)
        rows[i] = [j for j in range(n) if j != i and r[j] < cutoff]
    assert gap > 1e-3, f"a pair sits within 1e-3 of the cutoff ({gap})"
    return simple_graph(rows, n)


def ice_ih(nx, ny, nz):
    bond = 2.75
    z = 1.0 / 16.0
    frac = [(1 / 3, 2 / 3, z), (2 / 3, 1 / 3, 0.5 + z), (2 / 3, 1 / 3, 1 - z), (1 / 3, 2 / 3, 0.5 - z)]
    c = bond / 0.375
    a = c / math.sqrt(8.0 / 3.0)
    hv = np.array([[a, 0, 0], [-0.5 * a, a * math.sqrt(3) / 2, 0], [0, 0, c]])
    lengths = np.array([a, a * math.sqrt(3), c])
    basis = [np.array(f) @ hv for f in frac]
    basis += [b + hv[1] for b in basis]
    basis = [b - np.floor(b / lengths + 1e-8) * lengths for b in basis]
    xyz = [np.array([i, j, k]) * lengths + b
           for i in range(nx) for j in range(ny) for k in range(nz) for b in basis]
    box9 = [0, 0, 0, nx * lengths[0], ny * lengths[1], nz * lengths[2], 0, 0, 0]
    return np.array(xyz), box9


def ice_ic(side):
    a = 4 * 2.75 / math.sqrt(3)
    basis = [(0, 0, 0), (.5, .5, 0), (.5, 0, .5), (0, .5, .5),
             (.25, .25, .25), (.75, .75, .25), (.75, .25, .75), (.25, .75, .75)]
    xyz = [(np.array([i, j, k]) + np.array(b)) * a
           for i in range(side) for j in range(side) for k in range(side) for b in basis]
    return np.array(xyz), [0, 0, 0, side * a, side * a, side * a, 0, 0, 0]


def read_table(path, name):
    text = Path(path).read_text()
    body = re.search(r"%s\[\]\[3\] = \{(.*?)\n\};" % re.escape(name), text, re.S).group(1)
    nums = [float(v) for v in re.findall(r"-?\d+\.\d+(?:e-?\d+)?", body)]
    return np.array(nums).reshape(-1, 3)


def check(cond, what):
    if not cond:
        raise AssertionError(what)


def check_rows(adj, depth, rng, what):
    want = primitive_rings(adj, depth)
    rows = rows_of(adj, rng)
    got = [list(r) for r in yoda.ringNetwork(rows, depth)]
    check(got == want, f"{what}: ringNetwork depth {depth}\nwant {want}\ngot  {got}")
    return want


def check_frame(xyz, box9, cutoff, depth, what):
    adj = cutoff_graph(xyz, box9, cutoff)
    n = len(adj)
    off, cols = dseams.neighbours(xyz, box9, cutoff, dseams.CUTOFF)
    got_adj = simple_graph([list(cols[off[i]:off[i + 1]]) for i in range(n)], n)
    check(got_adj == adj, f"{what}: dseams.neighbours differs from the cutoff graph")
    want = primitive_rings(adj, depth)
    sizes, symmetry, hist, six_fraction = dseams.rings(xyz, box9, cutoff, depth)
    check(sorted(int(s) for s in sizes) == sorted(len(r) for r in want),
          f"{what}: dseams.rings sizes depth {depth}")
    for s in range(len(hist)):
        check(int(hist[s]) == sum(1 for r in want if len(r) == s), f"{what}: histogram size {s}")
    on_six = {v for r in want if len(r) == 6 for v in r}
    check(abs(six_fraction - len(on_six) / n) < 1e-12, f"{what}: six-ring fraction")
    return adj, want, list(symmetry)


def main():
    root = Path(sys.argv[1])
    rng = random.Random(20261121)

    for trial in range(300):
        depth = rng.randint(3, 10)
        n = rng.randint(0, 18 if depth <= 8 else 14)
        p = min(1.0, rng.uniform(0.5, 3.5) / max(n - 1, 1))
        rows = [[j for j in range(i + 1, n) if rng.random() < p] for i in range(n)]
        check_rows(simple_graph(rows, n), depth, rng, f"random trial {trial}")

    for k in range(3, 12):
        cycle = simple_graph([[(i + 1) % k] for i in range(k)], k)
        for depth in range(3, 12):
            want = check_rows(cycle, depth, rng, f"{k}-cycle")
            check(len(want) == (1 if depth >= k else 0), f"{k}-cycle at depth {depth}")

    for run in range(4):
        n = 24
        depth = rng.randint(4, 9)
        rows = [[j for j in range(i + 1, n) if rng.random() < 3.0 / (n - 1)] for i in range(n)]
        adj = simple_graph(rows, n)
        updater = yoda.RingUpdater(depth)
        for frame in range(6):
            got = [list(r) for r in updater.update(rows_of(adj))]
            check(got == primitive_rings(adj, depth), f"RingUpdater run {run} frame {frame}")
            nxt = [list(a) for a in adj]
            for _ in range(3):
                i, j = rng.randrange(n), rng.randrange(n)
                if i == j:
                    continue
                if j in nxt[i] or i in nxt[j]:
                    nxt[i] = [v for v in nxt[i] if v != j]
                    nxt[j] = [v for v in nxt[j] if v != i]
                else:
                    nxt[i].append(j)
            adj = simple_graph(nxt, n)

    ih, ih_box = ice_ih(4, 2, 2)
    ic, ic_box = ice_ic(2)
    for xyz, box9, name in ((ih, ih_box, "ice Ih"), (ic, ic_box, "ice Ic")):
        for depth in (6, 8, 10):
            adj, want, _ = check_frame(xyz, box9, 3.2, depth, f"{name} depth {depth}")
            check([list(r) for r in yoda.ringNetwork(rows_of(adj), depth)] == want,
                  f"{name}: yoda.ringNetwork on the cutoff rows")
        n = len(xyz)
        six = [r for r in want if len(r) == 6]
        hc, ddc = yoda.cageAffiliation(six, rows_of(adj))
        _, symmetry, _, _ = dseams.rings(xyz, box9, 3.2, 6)
        codes = [int(c) for c in symmetry]
        _, counts = dseams.cages(xyz, box9, 3.2, "cutoff")
        counts = [int(c) for c in counts]
        if name == "ice Ih":
            check(all(hc) and not any(ddc), "ice Ih: every six-ring is hexagonal-cage only")
            check(counts == [n, 0, 0, 0], f"ice Ih: dseams.cages counts {counts}")
            check(codes.count(2) == n // 2 and codes.count(3) == 3 * n // 2,
                  f"ice Ih: dseams.rings symmetry {sorted(set(codes))} "
                  f"basal {codes.count(2)} prismatic {codes.count(3)} hexagon {codes.count(10)}")
        else:
            check(all(ddc) and not any(hc), "ice Ic: every six-ring is double-diamond only")
            check(counts == [0, n, 0, 0], f"ice Ic: dseams.cages counts {counts}")
            check(codes.count(1) == 2 * n, f"ice Ic: dseams.rings symmetry {sorted(set(codes))}")

    si = read_table(root / "tests" / "hydrate_sites.inc", "kSIxyz")
    si_cell = float(re.search(r"kSIcell = ([0-9.]+);",
                              (root / "tests" / "hydrate_sites.inc").read_text()).group(1))
    si_box = [0, 0, 0, si_cell, si_cell, si_cell, 0, 0, 0]
    for depth in (6, 8):
        check_frame(si, si_box, 3.2, depth, f"sI depth {depth}")
    check(tuple(dseams.signature_cages(si, si_box, "512", 3.2, "cutoff")) == (2, 40), "sI 512 cages")
    check(tuple(dseams.signature_cages(si, si_box, "51262", 3.2, "cutoff")) == (6, 46),
          "sI 51262 cages")

    cells = root / "tests" / "ring_oracle_cells.inc"
    text = cells.read_text()
    for name in ("IceII", "IceIII"):
        xyz = read_table(cells, f"k{name}Xyz")
        box9 = [float(v) for v in re.search(r"k%sBox9\[9\] = \{([^}]*)\}" % name, text)
                .group(1).split(",")]
        adj = cutoff_graph(xyz, box9, 3.2)
        for depth in (6, 8, 10):
            check_rows(adj, depth, rng, f"{name} unit cell")

    print("ring oracle: PydSEAMSlib agrees")
    return 0


if __name__ == "__main__":
    sys.exit(main())
