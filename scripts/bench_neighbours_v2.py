#!/usr/bin/env python3
# Copyright (c) d-SEAMS developers
"""Triclinic neighbour list and end-to-end CHILL+, rings, diamond CNA.

Ice Ic, a = 6.35 A. The triclinic cell has xy = 0.25 L, xz = 0.05 L,
yz = 0.12 L. Every cutoff neighbour list uses 3.5 A. Symmetry is the
diamond/CNA kernel on that list, timed against PTM and identify-diamond.
Eight timed calls after one warmup, pinned to cores 0..T-1.
"""

import os
import subprocess
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
PY = os.environ.get("NEIGH_PYTHON", sys.executable)
BENCH = os.environ.get("NEIGH_BENCH", str(ROOT / "bbdir-release/tests/bench_v2"))
GL = "/lib/x86_64-linux-gnu:/usr/lib/x86_64-linux-gnu"


def ice_src():
    return r"""
import os, sys, time
import numpy as np
kind, molecules, threads, reps = sys.argv[1:]
molecules, threads, reps = int(molecules), int(threads), int(reps)
n = 1
while 8 * n * n * n < molecules:
    n += 1
a = 6.35
L = n * a
basis = np.array([
    [0, 0, 0], [0, 0.5, 0.5], [0.5, 0, 0.5], [0.5, 0.5, 0],
    [0.25, 0.25, 0.25], [0.25, 0.75, 0.75], [0.75, 0.25, 0.75],
    [0.75, 0.75, 0.25],
], float)
ix, iy, iz = np.mgrid[0:n, 0:n, 0:n]
origin = np.stack([ix, iy, iz], axis=-1).reshape(-1, 3)[:, None, :]
s = ((origin + basis) / n).reshape(-1, 3)
tilt = kind.startswith("tilt")
xy = 0.25 * L if tilt else 0.0
xz = 0.05 * L if tilt else 0.0
yz = 0.12 * L if tilt else 0.0
pts = np.column_stack([
    L * s[:, 0] + xy * s[:, 1] + xz * s[:, 2],
    L * s[:, 1] + yz * s[:, 2],
    L * s[:, 2],
])
# Row vectors a, b, c matching the positions above.
box = np.array([[L, 0.0, 0.0], [xy, L, 0.0], [xz, yz, L]], float)
cutoff = 3.5

def time_calls(fn):
    fn()
    samples = []
    sink = fn()
    for _ in range(reps):
        t0 = time.perf_counter()
        sink = fn()
        samples.append((time.perf_counter() - t0) * 1e3)
    samples.sort()
    med = samples[len(samples) // 2]
    print(f"kind={kind} molecules={molecules} threads={threads} sink={sink} "
          f"median_ms={med:.4f} samples_ms={','.join(f'{v:.4f}' for v in samples)}")

if kind in ("tilt-vesin", "tilt-freud"):
    if kind == "tilt-vesin":
        import vesin
        nl = vesin.NeighborList(cutoff=cutoff, full_list=True, n_threads=threads)
        def once():
            i, j = nl.compute(points=pts, box=box, periodic=True, quantities="ij", copy=False)
            return len(i)
    else:
        import freud
        freud.parallel.set_num_threads(threads)
        # Column-stacked cell vectors. The other orientation changes the pair count.
        fbox = freud.box.Box.from_matrix(box.T.copy())
        def once():
            aq = freud.locality.AABBQuery(fbox, pts)
            return len(aq.query(pts, dict(r_max=cutoff, exclude_ii=True)).toNeighborList().point_indices)
    time_calls(once)
elif kind.startswith("mdapy-"):
    from mdapy import System, set_num_threads
    set_num_threads(threads)
    system = System(positions=pts, box=box, elements="O")
    op = kind.split("-", 1)[1]
    if op == "chill":
        from mdapy.identify import ChillPlus
        tool = ChillPlus(cutoff)
        def once():
            result = tool.compute(system)
            return int(np.sum(result.structure == 2)) if hasattr(result, "structure") else 0
    elif op == "rings":
        from mdapy.analysis import RingStatistics
        tool = RingStatistics(max_size=6, cutoff=cutoff)
        def once():
            result = tool.compute(system)
            return int(np.sum(result.counts)) if hasattr(result, "counts") else 0
    elif op == "diamond":
        from mdapy.identify import DiamondStructure
        tool = DiamondStructure()
        def once():
            result = tool.compute(system)
            return int(np.sum(np.asarray(result.structure) > 0))
    else:
        from mdapy.identify import PolyhedralTemplateMatching
        # Cubic diamond is DCUB and hexagonal diamond is DHEX. The default
        # set is the close-packed templates, which do not match this ice.
        try:
            tool = PolyhedralTemplateMatching(structure=["dcub", "dhex"])
        except TypeError:
            tool = PolyhedralTemplateMatching()
        def once():
            result = tool.compute(system)
            return int(np.sum(np.asarray(result.structure) > 0))
    time_calls(once)
elif kind.startswith("ovito-"):
    from ovito.data import DataCollection, Particles, SimulationCell
    from ovito.pipeline import Pipeline, StaticSource
    from ovito.modifiers import (
        ChillPlusModifier, CreateBondsModifier, FindRingsModifier,
        IdentifyDiamondModifier, PolyhedralTemplateMatchingModifier)
    op = kind.split("-", 1)[1]
    def once():
        # A pipeline caches its output. Build a new one so each sample
        # runs the modifier on the coordinates.
        data = DataCollection()
        cell = SimulationCell(pbc=(True, True, True))
        matrix = np.zeros((3, 4))
        matrix[:, 0] = box[0]
        matrix[:, 1] = box[1]
        matrix[:, 2] = box[2]
        cell[:, :] = matrix
        parts = Particles()
        parts.create_property("Position", data=pts)
        data.objects.append(cell)
        data.objects.append(parts)
        pipe = Pipeline(source=StaticSource(data=data))
        if op == "chill":
            pipe.modifiers.append(ChillPlusModifier(cutoff=cutoff))
        elif op == "rings":
            pipe.modifiers.append(CreateBondsModifier(cutoff=cutoff))
            pipe.modifiers.append(
                FindRingsModifier(minimum_ring_size=3, maximum_ring_size=6))
        elif op == "diamond":
            pipe.modifiers.append(IdentifyDiamondModifier())
        else:
            mod = PolyhedralTemplateMatchingModifier(rmsd_cutoff=0.1)
            # PTM_MATCH_DCUB is 6 and PTM_MATCH_DHEX is 7. OVITO leaves both off.
            for structure in mod.structures:
                structure.enabled = structure.id in (6, 7)
            pipe.modifiers.append(mod)
        out = pipe.compute()
        if op in ("symmetry", "diamond", "chill") and "Structure Type" in out.particles:
            return int(np.sum(np.asarray(out.particles["Structure Type"]) > 0))
        return int(out.particles.count)
    time_calls(once)
elif kind == "freud-none":
    print(f"kind=freud-none molecules={molecules} threads={threads} "
          "sink=0 median_ms=nan samples_ms=no-chill-rings-or-ptm")
else:
    raise SystemExit("unknown " + kind)
"""


def run(cmd, env, cores):
    proc = subprocess.run(
        ["taskset", "-c", cores, *cmd], env=env, text=True, capture_output=True)
    if proc.returncode != 0:
        sys.stderr.write(proc.stdout)
        sys.stderr.write(proc.stderr)
        raise SystemExit(f"failed {cmd[0]} ({proc.returncode})")
    line = proc.stdout.strip().splitlines()[-1]
    print(line, flush=True)
    return line


def env_for(threads, gl):
    env = os.environ.copy()
    env.pop("LD_LIBRARY_PATH", None)
    env.update({
        "OMP_NUM_THREADS": str(threads),
        "OMP_DYNAMIC": "FALSE",
        "OMP_PROC_BIND": "close",
        "OMP_PLACES": "cores",
        "OMP_WAIT_POLICY": "PASSIVE",
        "OPENBLAS_NUM_THREADS": str(threads),
        "MKL_NUM_THREADS": str(threads),
        "OVITO_THREAD_COUNT": str(threads),
    })
    if gl:
        env["LD_LIBRARY_PATH"] = GL
    return env


def main():
    reps = "8"
    lines = []
    jobs = []
    for molecules in (10648, 110592, 884736):
        for threads in (1, 8):
            jobs.append(("dseams", "tilt", molecules, threads))
            jobs.append(("py", "tilt-vesin", molecules, threads))
            jobs.append(("py", "tilt-freud", molecules, threads))
    for mode in ("chill", "rings", "symmetry"):
        for molecules in (110592, 884736):
            for threads in (1, 8):
                jobs.append(("dseams", mode, molecules, threads))
                jobs.append(("py", f"mdapy-{mode}", molecules, threads))
                jobs.append(("py", f"ovito-{mode}", molecules, threads))
                if mode == "chill":
                    jobs.append(("py", "freud-none", molecules, threads))
                if mode == "symmetry":
                    jobs.append(("py", "mdapy-diamond", molecules, threads))
                    jobs.append(("py", "ovito-diamond", molecules, threads))
    for who, kind, molecules, threads in jobs:
        cores = "0" if threads == 1 else f"0-{threads - 1}"
        if who == "dseams":
            line = run([BENCH, kind, str(molecules), str(threads), reps],
                       env_for(threads, False), cores)
        else:
            line = run([PY, "-c", ice_src(), kind, str(molecules), str(threads), reps],
                       env_for(threads, kind.startswith("ovito")), cores)
        lines.append(line)
    out = Path(sys.argv[1]) if len(sys.argv) > 1 else None
    if out is not None:
        out.write_text("\n".join(lines) + "\n")


if __name__ == "__main__":
    main()
