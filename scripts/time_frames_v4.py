#!/usr/bin/env python3
"""Serial timings. One thread, eight samples after one warm-up.

The header records the online core count and not a host name.
"""

from __future__ import annotations

import os
import subprocess
import time

import numpy as np

ROOT = os.path.abspath(os.path.join(os.path.dirname(__file__), ".."))
BENCH = os.path.join(ROOT, "build-rel", "tests", "bench_ring_ic")
CUTOFF = 3.5
CORES = int(os.popen("getconf _NPROCESSORS_ONLN").read().strip() or "0")


def header():
    print(f"cores {CORES}", flush=True)
    print("threads 1", flush=True)
    print("samples 8", flush=True)
    print("serial 1", flush=True)


def load_frame(path):
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
    header_cols = lines[i].split()[2:]
    i += 1
    cols = {name: k for k, name in enumerate(header_cols)}
    pts = np.empty((n, 3))
    for k in range(n):
        tok = lines[i].split()
        i += 1
        pts[k] = [float(tok[cols[a]]) for a in "xyz"]
    return pts, np.array(lo), np.array(hi)


def tile(pts, lo, hi, nx, ny, nz):
    side = hi - lo
    copies = []
    for ix in range(nx):
        for iy in range(ny):
            for iz in range(nz):
                shift = np.array([ix * side[0], iy * side[1], iz * side[2]])
                copies.append(pts + shift)
    out = np.vstack(copies)
    nlo = lo.copy()
    nhi = lo + np.array([nx * side[0], ny * side[1], nz * side[2]])
    return out, nlo, nhi


def graph(pts, lo, hi, cutoff=CUTOFF):
    side = hi - lo
    rel = pts - lo
    n = len(pts)
    nc = [max(1, int(side[a] / cutoff)) for a in range(3)]
    cell = side / np.array(nc)
    buckets = {}

    def key(p):
        return tuple(min(nc[a] - 1, max(0, int(p[a] / cell[a]))) for a in range(3))

    for i, p in enumerate(rel):
        buckets.setdefault(key(p), []).append(i)
    adj = [[] for _ in range(n)]
    cut2 = cutoff * cutoff
    for i, p in enumerate(rel):
        c = key(p)
        for dx in (-1, 0, 1):
            for dy in (-1, 0, 1):
                for dz in (-1, 0, 1):
                    nb = ((c[0] + dx) % nc[0], (c[1] + dy) % nc[1], (c[2] + dz) % nc[2])
                    for j in buckets.get(nb, []):
                        if j <= i:
                            continue
                        d = p - rel[j]
                        for a in range(3):
                            if d[a] > side[a] / 2:
                                d[a] -= side[a]
                            elif d[a] < -side[a] / 2:
                                d[a] += side[a]
                        if 0 < float(d @ d) < cut2:
                            adj[i].append(int(j))
                            adj[j].append(int(i))
    return adj


def write_dump(path, pts, lo, hi):
    with open(path, "w") as out:
        out.write("ITEM: TIMESTEP\n0\nITEM: NUMBER OF ATOMS\n")
        out.write(f"{len(pts)}\nITEM: BOX BOUNDS pp pp pp\n")
        for a in range(3):
            out.write(f"{lo[a]} {hi[a]}\n")
        out.write("ITEM: ATOMS id type x y z\n")
        for i, p in enumerate(pts, start=1):
            out.write(f"{i} 1 {p[0]} {p[1]} {p[2]}\n")


def dseams_counts(dump, depth):
    env = os.environ.copy()
    env["LD_LIBRARY_PATH"] = os.path.join(ROOT, "build-rel", "src")
    env["OMP_NUM_THREADS"] = "1"
    env.pop("SEAMS_RING_PROFILE", None)
    proc = subprocess.run(
        [BENCH, "--dump", dump, "--frame", "1", "--depth", str(depth), "--threads", "1"],
        env=env,
        check=True,
        capture_output=True,
        text=True,
    )
    for line in proc.stdout.splitlines():
        if line.startswith("dump "):
            return line
    return proc.stdout


def sqq_count(adj, depth):
    from sqq.core.sqq_py.ring import find_rings

    found = find_rings(
        {i: set(row) for i, row in enumerate(adj)},
        list(range(3, depth + 1)),
        chordless=True,
        definition="chordless",
    )
    return {s: len(found.get(s, [])) for s in range(3, depth + 1)}


def ovito_count(pts, lo, hi, depth):
    os.environ["QT_QPA_PLATFORM"] = "offscreen"
    os.environ["OVITO_THREAD_COUNT"] = "1"
    from ovito.data import DataCollection
    from ovito.modifiers import CreateBondsModifier, FindRingsModifier
    from ovito.pipeline import Pipeline, StaticSource

    side = hi - lo
    mat = np.zeros((3, 4))
    mat[:, :3] = np.diag(side)
    mat[:, 3] = lo
    data = DataCollection()
    data.create_cell(mat, pbc=(True, True, True))
    parts = data.create_particles(count=len(pts))
    parts.create_property("Position", data=np.asarray(pts))
    pipe = Pipeline(source=StaticSource(data=data))
    pipe.modifiers.append(CreateBondsModifier(cutoff=CUTOFF))
    pipe.modifiers.append(FindRingsModifier(maximum_ring_size=depth))
    out = pipe.compute()
    hist = out.tables["ring-size-histogram"]
    return {int(s): int(c) for s, c in zip(hist["Ring Size"][...], hist["Count"][...])}


def time_call(fn, repeats=8, warmup=1):
    for _ in range(warmup):
        fn()
    samples = []
    last = None
    for _ in range(repeats):
        t0 = time.perf_counter()
        last = fn()
        samples.append((time.perf_counter() - t0) * 1e3)
    return samples, last


def median(xs):
    ys = sorted(xs)
    return ys[len(ys) // 2]


def main():
    header()
    env = os.environ.copy()
    env["LD_LIBRARY_PATH"] = os.path.join(ROOT, "build-rel", "src")
    env["OMP_NUM_THREADS"] = "1"
    env["OMP_WAIT_POLICY"] = "PASSIVE"
    env.pop("SEAMS_RING_PROFILE", None)
    print("tool dseams-e2e depth 10 molecules 110592", flush=True)
    proc = subprocess.run(
        [BENCH, "--e2e", "--threads", "1", "--depth", "10", "24"],
        env=env,
        check=True,
        capture_output=True,
        text=True,
    )
    for line in proc.stdout.splitlines():
        if line.startswith("sample ") or line.startswith("summary "):
            print(line, flush=True)

    jobs = [
        ("nucleation", "/tmp/frames/nucleation.lammpstrj", 1, (1, 1, 1), 6),
        ("nucleation", "/tmp/frames/nucleation.lammpstrj", 1, (6, 6, 6), 6),
        ("interface", os.path.join(ROOT, "input/traj/interface_f0.lammpstrj"), 1, (1, 1, 1), 6),
        ("interface", os.path.join(ROOT, "input/traj/interface_f0.lammpstrj"), 1, (16, 8, 8), 6),
    ]
    for name, path, frame, reps, depth in jobs:
        pts, lo, hi = load_frame(path)
        pts, lo, hi = tile(pts, lo, hi, *reps)
        tag = f"{name} tile {reps[0]}x{reps[1]}x{reps[2]} molecules {len(pts)} depth {depth}"
        dump = f"/tmp/v4-{name}-{len(pts)}.lammpstrj"
        if not os.path.exists(dump) or os.path.getsize(dump) < 1000:
            write_dump(dump, pts, lo, hi)
        print(f"system {tag}", flush=True)

        dseams_ms = []
        line = ""
        for i in range(9):
            line = dseams_counts(dump, depth)
            ms = float(line.split("total_ms")[1].split()[0])
            if i:
                dseams_ms.append(ms)
        print(
            f"sample tool dseams {tag} median_ms {median(dseams_ms):.3f} {line}",
            flush=True,
        )
        adj = graph(pts, lo, hi)

        def run_sqq(adj=adj, depth=depth):
            return sqq_count(adj, depth)

        samples, counts = time_call(run_sqq)
        print(
            f"sample tool sqq {tag} median_ms {median(samples):.3f} {counts}",
            flush=True,
        )

        def run_ovito(pts=pts, lo=lo, hi=hi, depth=depth):
            return ovito_count(pts, lo, hi, depth)

        samples, counts = time_call(run_ovito)
        print(
            f"sample tool ovito {tag} median_ms {median(samples):.3f} {counts}",
            flush=True,
        )


if __name__ == "__main__":
    main()
