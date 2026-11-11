#!/usr/bin/env python3
# Copyright (c) d-SEAMS developers
"""Time the cutoff neighbour list against vesin, freud and matscipy.

Ice Ic, oxygen lattice constant 6.35 A, cutoff 3.5 A. The list is the
full both-direction shell (four neighbours, 4N pairs). Each sample is
one build in a fresh process, pinned to the first ``threads`` cores.
"""

import argparse
import os
import subprocess
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
PY = os.environ.get("NEIGH_PYTHON", sys.executable)
BENCH = os.environ.get(
    "NEIGH_BENCH", str(ROOT / "build-release/tests/bench_cutoff_list"))
SIZES = (10648, 110592, 884736)
THREADS = (1, 2, 4, 8)


def ice_code():
    return r"""
import os, sys, time
import numpy as np
kind, molecules, threads, reps = sys.argv[1], int(sys.argv[2]), int(sys.argv[3]), int(sys.argv[4])
n = 1
while n * n * n * 8 < molecules:
    n += 1
a = 6.35
basis = np.array([
    [0, 0, 0], [0, 0.5, 0.5], [0.5, 0, 0.5], [0.5, 0.5, 0],
    [0.25, 0.25, 0.25], [0.25, 0.75, 0.75], [0.75, 0.25, 0.75],
    [0.75, 0.75, 0.25],
])
ix, iy, iz = np.mgrid[0:n, 0:n, 0:n]
origin = np.stack([ix, iy, iz], axis=-1).reshape(-1, 3)[:, None, :]
pts = ((origin + basis[None, :, :]) * a).reshape(-1, 3)
L = n * a
cutoff = 3.5

if kind == "vesin":
    import vesin
    nl = vesin.NeighborList(cutoff=cutoff, full_list=True, n_threads=threads)
    def once():
        i, j = nl.compute(points=pts, box=np.eye(3) * L, periodic=True,
                          quantities="ij", copy=False)
        return len(i)
elif kind == "freud":
    import freud
    freud.parallel.set_num_threads(threads)
    fbox = freud.box.Box.cube(L)
    def once():
        aq = freud.locality.AABBQuery(fbox, pts)
        neigh = aq.query(pts, dict(r_max=cutoff, exclude_ii=True)).toNeighborList()
        return len(neigh.point_indices)
elif kind == "matscipy":
    from ase import Atoms
    from matscipy.neighbours import neighbour_list
    atoms = Atoms(positions=pts, cell=[L, L, L], pbc=True)
    def once():
        i, j = neighbour_list("ij", atoms, cutoff)
        return len(i)
else:
    raise SystemExit("unknown kind")

for _ in range(2):
    pairs = once()
samples = []
for _ in range(reps):
    t0 = time.perf_counter()
    pairs = once()
    samples.append((time.perf_counter() - t0) * 1e3)
samples.sort()
med = samples[len(samples) // 2]
print(f"kind={kind} molecules={molecules} threads={threads} pairs={pairs} "
      f"median_ms={med:.4f} samples_ms={','.join(f'{s:.4f}' for s in samples)}")
"""


def run(cmd, env, cores):
    pin = ["taskset", "-c", cores]
    proc = subprocess.run(pin + cmd, env=env, text=True, capture_output=True)
    if proc.returncode != 0:
        sys.stderr.write(proc.stdout)
        sys.stderr.write(proc.stderr)
        raise SystemExit(f"failed: {cmd} ({proc.returncode})")
    return proc.stdout.strip()


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--reps", type=int, default=8)
    parser.add_argument("--out", type=Path, default=None)
    args = parser.parse_args()
    lines = []
    for molecules in SIZES:
        for threads in THREADS:
            cores = "0" if threads == 1 else f"0-{threads - 1}"
            env = os.environ.copy()
            env.update({
                "OMP_NUM_THREADS": str(threads),
                "OMP_DYNAMIC": "FALSE",
                "OMP_PROC_BIND": "close",
                "OMP_PLACES": "cores",
                "OMP_WAIT_POLICY": "PASSIVE",
                "OPENBLAS_NUM_THREADS": str(threads),
                "MKL_NUM_THREADS": str(threads),
                "RAYON_NUM_THREADS": str(threads),
            })
            text = run([BENCH, str(molecules), str(threads), str(args.reps)], env, cores)
            lines.append(text)
            print(text, flush=True)
            for kind in ("vesin", "freud", "matscipy"):
                text = run([PY, "-c", ice_code(), kind, str(molecules),
                            str(threads), str(args.reps)], env, cores)
                lines.append(text)
                print(text, flush=True)
    if args.out is not None:
        args.out.write_text("\n".join(lines) + "\n")


if __name__ == "__main__":
    main()
