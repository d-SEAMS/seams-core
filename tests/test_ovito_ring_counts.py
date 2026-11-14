#!/usr/bin/env python3
"""Primitive rings against OVITO FindRings.

Both sides see the same cutoff graph. d-SEAMS counts Franzblau primitive
rings. OVITO's shortest-path census is a strict subset: on K2,3, which
has three chordless 4-rings and no shorter chord, FindRings reports 2.
Exit 77 only when OVITO cannot be imported.
"""

from __future__ import annotations

import math
import os
import subprocess
import sys
import tempfile
from collections import deque

CUTOFF = 3.5
DEPTH = 10
SIZES = range(4, 11)
DUMP = os.path.join("input", "traj", "interface_f0.lammpstrj")


def prepare_ovito() -> bool:
    site = os.environ.get("SEAMS_OVITO_SITE", "")
    if site and os.path.isdir(site) and site not in sys.path:
        sys.path.insert(0, site)
    extra = []
    gl = "/usr/lib/x86_64-linux-gnu"
    qt = os.path.join(site, "PySide6", "Qt", "lib")
    for path in (gl, qt):
        if os.path.isdir(path):
            extra.append(path)
    if extra and os.environ.get("SEAMS_OVITO_LIBS") != "1":
        current = os.environ.get("LD_LIBRARY_PATH", "")
        os.environ["LD_LIBRARY_PATH"] = ":".join(extra + ([current] if current else []))
        os.environ["SEAMS_OVITO_LIBS"] = "1"
        os.environ.setdefault("QT_QPA_PLATFORM", "offscreen")
        os.execv(sys.executable, [sys.executable, *sys.argv])
    os.environ.setdefault("QT_QPA_PLATFORM", "offscreen")
    try:
        import ovito  # noqa: F401

        return True
    except Exception as exc:
        print("ovito is not usable: {}".format(exc), file=sys.stderr)
        return False


def load_frame(path: str):
    import numpy as np

    lines = open(path).read().splitlines()
    i = 0
    while "NUMBER OF ATOMS" not in lines[i]:
        i += 1
    n = int(lines[i + 1])
    i += 2
    while not lines[i].startswith("ITEM: BOX"):
        i += 1
    i += 1
    lo, hi = [], []
    for _ in range(3):
        a, b = map(float, lines[i].split()[:2])
        lo.append(a)
        hi.append(b)
        i += 1
    header = lines[i].split()[2:]
    i += 1
    cols = {name: k for k, name in enumerate(header)}
    pts = np.empty((n, 3))
    for k in range(n):
        tok = lines[i].split()
        i += 1
        pts[k] = [float(tok[cols[a]]) for a in "xyz"]
    return pts, lo, hi


def bench(bench_bin: str, dump: str, depth: int, extra: list[str]) -> str:
    env = os.environ.copy()
    root = os.path.abspath(os.path.join(os.path.dirname(bench_bin), ".."))
    env["LD_LIBRARY_PATH"] = os.path.join(root, "src")
    env["OMP_NUM_THREADS"] = "1"
    env.pop("SEAMS_RING_PROFILE", None)
    proc = subprocess.run(
        [bench_bin, "--dump", dump, "--frame", "1", "--depth", str(depth), "--threads", "1", *extra],
        env=env,
        check=True,
        capture_output=True,
        text=True,
    )
    return proc.stdout


def summary_counts(text: str) -> dict[int, int]:
    counts = {s: 0 for s in SIZES}
    for line in text.splitlines():
        if not line.startswith("dump "):
            continue
        parts = line.split()
        for s in SIZES:
            key = "n{}".format(s)
            if key in parts:
                counts[s] = int(parts[parts.index(key) + 1])
    return counts


def summary_edges(text: str) -> int:
    for line in text.splitlines():
        if line.startswith("dump ") and " edges " in line:
            parts = line.split()
            return int(parts[parts.index("edges") + 1])
    raise SystemExit("bench did not report an edge count")


def listed_rings(text: str) -> list[list[int]]:
    rings = []
    for line in text.splitlines():
        parts = line.split()
        if parts and parts[0].isdigit():
            rings.append([int(item) for item in parts])
    return rings


def ovito_histogram(pts, lo, hi, depth: int) -> tuple[dict[int, int], int]:
    import numpy as np
    from ovito.data import DataCollection
    from ovito.modifiers import CreateBondsModifier, FindRingsModifier
    from ovito.pipeline import Pipeline, StaticSource

    os.environ["QT_QPA_PLATFORM"] = "offscreen"
    os.environ.setdefault("OVITO_THREAD_COUNT", "1")
    side = [hi[a] - lo[a] for a in range(3)]
    mat = np.zeros((3, 4))
    mat[:, :3] = np.diag(side)
    mat[:, 3] = lo
    data = DataCollection()
    data.create_cell(mat, pbc=(True, True, True))
    parts = data.create_particles(count=len(pts))
    parts.create_property("Position", data=np.asarray(pts))
    pipe = Pipeline(source=StaticSource(data=data))
    pipe.modifiers.append(CreateBondsModifier(cutoff=CUTOFF))
    pipe.modifiers.append(FindRingsModifier(minimum_ring_size=4, maximum_ring_size=depth))
    out = pipe.compute()
    hist = out.tables["ring-size-histogram"]
    found = {int(s): int(c) for s, c in zip(hist["Ring Size"][...], hist["Count"][...])}
    bonds = out.particles.bonds.count if out.particles.bonds is not None else 0
    return {s: found.get(s, 0) for s in range(4, depth + 1)}, bonds


def cutoff_graph(pts, side):
    import numpy as np

    n = len(pts)
    side = np.asarray(side, dtype=float)
    cut2 = CUTOFF * CUTOFF
    adj = [set() for _ in range(n)]
    edges = 0
    for a in range(n):
        delta = pts[a + 1 :] - pts[a]
        delta -= side * np.round(delta / side)
        r2 = np.einsum("ij,ij->i", delta, delta)
        for j in np.where(r2 < cut2)[0]:
            b = a + 1 + int(j)
            adj[a].add(b)
            adj[b].add(a)
            edges += 1
    return adj, edges


def primitive(rings, adj) -> int:
    """Rings whose shorter arc is a shortest path. Returns the non-primitive count."""

    def dist(src, limit):
        seen = {src: 0}
        queue = deque([src])
        while queue:
            u = queue.popleft()
            if seen[u] == limit:
                continue
            for v in adj[u]:
                if v not in seen:
                    seen[v] = seen[u] + 1
                    queue.append(v)
        return seen

    bad = 0
    for ring in rings:
        length = len(ring)
        ok = True
        for i, a in enumerate(ring):
            reached = dist(a, length // 2)
            for j in range(i + 1, length):
                arc = min(j - i, length - (j - i))
                if reached.get(ring[j], 10**9) < arc:
                    ok = False
                    break
            if not ok:
                break
        if not ok:
            bad += 1
    return bad


def write_k23(path: str) -> None:
    """K2,3 embedded so a 3.5 A cutoff reproduces exactly its six edges.

    Vertices 0 and 3 are one part. Vertices 1, 2 and 4 sit on a circle
    between them, each within the cutoff of both and outside it of each
    other. The three 4-cycles are chordless, so each is primitive.
    """
    radius = 2.5
    pts = [
        (0.0, 0.0, 0.0),
        (2.0, radius, 0.0),
        (2.0, -radius / 2.0, radius * math.sqrt(3.0) / 2.0),
        (4.0, 0.0, 0.0),
        (2.0, -radius / 2.0, -radius * math.sqrt(3.0) / 2.0),
    ]
    lines = [
        "ITEM: TIMESTEP",
        "0",
        "ITEM: NUMBER OF ATOMS",
        "5",
        "ITEM: BOX BOUNDS pp pp pp",
        "0 30",
        "0 30",
        "0 30",
        "ITEM: ATOMS id type x y z",
    ]
    for i, (x, y, z) in enumerate(pts, start=1):
        lines.append("{} 1 {:.8f} {:.8f} {:.8f}".format(i, x, y, z))
    with open(path, "w", encoding="utf-8") as handle:
        handle.write("\n".join(lines) + "\n")


def main() -> None:
    if not prepare_ovito():
        print("ovito is not importable", file=sys.stderr)
        raise SystemExit(77)
    if len(sys.argv) != 2:
        print("usage: test_ovito_ring_counts.py BENCH", file=sys.stderr)
        raise SystemExit(2)
    bench_bin = sys.argv[1]

    with tempfile.TemporaryDirectory(prefix="seams-k23-") as work:
        dump = os.path.join(work, "k23.lammpstrj")
        write_k23(dump)
        pts, lo, hi = load_frame(dump)
        side = [hi[a] - lo[a] for a in range(3)]
        _adj, edges = cutoff_graph(pts, side)
        if edges != 6:
            print("K2,3 cutoff graph has {} edges, expected 6".format(edges), file=sys.stderr)
            raise SystemExit(1)
        ours = summary_counts(bench(bench_bin, dump, 4, []))
        theirs, obonds = ovito_histogram(pts, lo, hi, 4)
        if ours[4] != 3 or obonds != 6 or theirs[4] != 2:
            print("K2,3 dseams n4 {} ovito n4 {} bonds {}".format(ours[4], theirs[4], obonds), file=sys.stderr)
            raise SystemExit(1)
        print("K2,3 primitive 4-rings 3, OVITO FindRings 2, bonds 6")

    dump = os.path.abspath(DUMP)
    pts, lo, hi = load_frame(dump)
    side = [hi[a] - lo[a] for a in range(3)]
    adj, brute = cutoff_graph(pts, side)
    head = bench(bench_bin, dump, DEPTH, [])
    listed = bench(bench_bin, dump, DEPTH, ["--list"])
    rings = listed_rings(listed)
    ours = summary_counts(head)
    edges = summary_edges(head)
    theirs, obonds = ovito_histogram(pts, lo, hi, DEPTH)
    if edges != brute or obonds != brute:
        print("edges bench {} brute {} ovito {}".format(edges, brute, obonds), file=sys.stderr)
        raise SystemExit(1)
    bad = primitive(rings, adj)
    listed_sizes = {}
    for ring in rings:
        listed_sizes[len(ring)] = listed_sizes.get(len(ring), 0) + 1
    # The bench list includes 3-rings. The comparison histogram starts at 4.
    mismatch = [s for s in SIZES if listed_sizes.get(s, 0) != ours[s]]
    if bad or mismatch:
        print("non-primitive {} size mismatch {}".format(bad, mismatch), file=sys.stderr)
        raise SystemExit(1)
    # FindRings is a subset of the primitive census, not the same histogram.
    over = [s for s in SIZES if theirs[s] > ours[s]]
    if over:
        print("OVITO exceeds the primitive census at {}".format(over), file=sys.stderr)
        print("dseams {}".format(ours), file=sys.stderr)
        print("ovito  {}".format(theirs), file=sys.stderr)
        raise SystemExit(1)
    print(
        "same bonds",
        brute,
        "primitive",
        " ".join("n{}={}".format(s, ours[s]) for s in SIZES),
        "ovito",
        " ".join("n{}={}".format(s, theirs[s]) for s in SIZES),
    )


if __name__ == "__main__":
    main()
