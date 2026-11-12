"""Head-to-head local q4 and q6 against OVITO and freud.

Ideal ice Ic and one liquid made by shifting that crystal with a normal
draw (seed 91001, sigma 0.9 angstrom). The frames stay in a temporary
directory. One warm-up, then eight samples. The median of an even count
is the mean of the two central values. Quartiles are Hyndman-Fan type 7.
Taskset pins one core or eight. The neighbour build is inside every sample.
"""

from __future__ import annotations

import argparse
import math
import os
import struct
import subprocess
import sys
import tempfile
from pathlib import Path

import numpy as np

ROOT = Path(__file__).resolve().parents[1]
CUTOFF = 3.5
WARMUP = 1
REPEATS = 8
ORDERS = (4, 6)
CELLS = (11, 24, 48)


def ice_ic(cells, bond=2.75):
    lattice = 4.0 * bond / math.sqrt(3.0)
    side = cells * lattice
    fcc = np.array([[0.0, 0.0, 0.0], [0.5, 0.5, 0.0], [0.5, 0.0, 0.5], [0.0, 0.5, 0.5]])
    basis = np.vstack([fcc, fcc + 0.25])
    ijk = np.stack(np.meshgrid(*[np.arange(cells)] * 3, indexing="ij"), -1).reshape(-1, 3)
    pts = (ijk[:, None, :] + basis).reshape(-1, 3) * lattice
    return np.ascontiguousarray(pts % side), float(side)


def liquid(cells, seed=91001, sigma=0.9):
    pts, side = ice_ic(cells)
    rng = np.random.default_rng(seed)
    return np.ascontiguousarray(np.mod(pts + rng.normal(0.0, sigma, pts.shape), side)), side


def write_frame(path, pts, side):
    with open(path, "wb") as handle:
        handle.write(struct.pack("<I", len(pts)))
        handle.write(struct.pack("<d", float(side)))
        handle.write(np.ascontiguousarray(pts, dtype=np.float64).tobytes())


def quantile(sorted_vals, q):
    pos = q * (len(sorted_vals) - 1)
    lo = int(pos)
    frac = pos - lo
    if lo + 1 >= len(sorted_vals):
        return sorted_vals[-1]
    return sorted_vals[lo] * (1.0 - frac) + sorted_vals[lo + 1] * frac


def summarize(samples):
    ordered = sorted(float(v) for v in samples)
    mid = len(ordered) // 2
    median = ordered[mid] if len(ordered) % 2 else 0.5 * (ordered[mid - 1] + ordered[mid])
    return median, quantile(ordered, 0.75) - quantile(ordered, 0.25)


def bench_path():
    env = os.environ.get("SEAMS_QL_BENCH")
    if env:
        return env
    for rel in ("bbdir-release/tests/bench_ql", "bbdir/tests/bench_ql"):
        path = ROOT / rel
        if path.is_file():
            return str(path)
    raise SystemExit("bench_ql is not built")


def run_seams(binary, frame, order, threads, ql_path):
    cpus = "0" if threads == 1 else "0-7"
    env = os.environ.copy()
    env.update({
        "OMP_NUM_THREADS": str(threads),
        "OMP_DYNAMIC": "FALSE",
        "OMP_PROC_BIND": "false",
        "OMP_WAIT_POLICY": "active",
    })
    proc = subprocess.run(
        ["taskset", "-c", cpus, binary, frame, str(order), str(REPEATS), ql_path],
        cwd=ROOT, env=env, capture_output=True, text=True, check=True,
    )
    samples = [float(line) for line in proc.stdout.split()]
    if len(samples) != REPEATS:
        raise RuntimeError(proc.stderr)
    return samples


def run_other(code, frame, order, threads, work):
    cpus = "0" if threads == 1 else "0-7"
    script = r"""
import os, struct, sys, time
import numpy as np
path, order_s, threads_s, code, out = sys.argv[1:]
order = int(order_s)
threads = int(threads_s)
with open(path, "rb") as handle:
    n = struct.unpack("<I", handle.read(4))[0]
    side = struct.unpack("<d", handle.read(8))[0]
    pts = np.frombuffer(handle.read(), dtype=np.float64).reshape(n, 3).copy()
print(f"cores={len(os.sched_getaffinity(0))}", flush=True)
samples = []
last = None
if code == "ovito":
    from ovito.data import Particles, SimulationCell, DataCollection
    from ovito.modifiers import SteinhardtParametersModifier
    from ovito.pipeline import StaticSource, Pipeline
    data = DataCollection()
    cell = SimulationCell(pbc=(True, True, True))
    cell[:, :3] = np.diag([side, side, side])
    data.objects.append(cell)
    part = Particles()
    part.create_property("Position", data=pts)
    data.objects.append(part)
    for i in range(1 + 8):
        pipe = Pipeline(source=StaticSource(data=data))
        mod = SteinhardtParametersModifier(cutoff=3.5, max_l=order, compute_averaged=False)
        mod.neighbor_mode = SteinhardtParametersModifier.NeighborMode.Cutoff
        pipe.modifiers.append(mod)
        t0 = time.perf_counter()
        result = pipe.compute()
        dt = (time.perf_counter() - t0) * 1e3
        last = np.ascontiguousarray(np.asarray(result.particles["Steinhardt"], dtype=np.float64)[:, -1])
        if i >= 1:
            samples.append(dt)
else:
    import freud
    freud.parallel.set_num_threads(threads)
    box = freud.box.Box.cube(side)
    for i in range(1 + 8):
        op = freud.order.Steinhardt(l=order, average=False)
        t0 = time.perf_counter()
        op.compute((box, pts), neighbors={"r_max": 3.5})
        dt = (time.perf_counter() - t0) * 1e3
        last = np.asarray(op.ql, dtype=np.float64).copy()
        if i >= 1:
            samples.append(dt)
np.save(out, last)
for value in samples:
    print(f"{value:.6f}")
"""
    env = os.environ.copy()
    env["OMP_NUM_THREADS"] = str(threads)
    env["OVITO_THREAD_COUNT"] = str(threads)
    env["OMP_DYNAMIC"] = "FALSE"
    out = work / f"{code}.npy"
    proc = subprocess.run(
        ["taskset", "-c", cpus, sys.executable, "-c", script, frame, str(order), str(threads), code, str(out)],
        env=env, capture_output=True, text=True,
    )
    if proc.returncode != 0:
        sys.stderr.write(proc.stderr[-2000:])
        raise RuntimeError(code)
    lines = [line for line in proc.stdout.splitlines() if line]
    cores = lines[0].split("=", 1)[1]
    samples = [float(line) for line in lines[1:]]
    return samples, np.load(out), cores


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--out", required=True)
    args = parser.parse_args()
    binary = bench_path()
    rows = ["code,frame,molecules,l,threads,cores,sample_ms"]
    with tempfile.TemporaryDirectory(prefix="ql-") as folder:
        work = Path(folder)
        for kind in ("ice", "liquid"):
            for cells in CELLS:
                pts, side = ice_ic(cells) if kind == "ice" else liquid(cells)
                frame = work / f"{kind}-{cells}.bin"
                write_frame(frame, pts, side)
                n = len(pts)
                for order in ORDERS:
                    for threads in (1, 8):
                        ql_path = work / "seams.ql"
                        seams = run_seams(binary, str(frame), order, threads, str(ql_path))
                        seams_ql = np.fromfile(ql_path, dtype=np.float64)
                        ov, ov_ql, ov_cores = run_other("ovito", str(frame), order, threads, work)
                        fr, fr_ql, fr_cores = run_other("freud", str(frame), order, threads, work)
                        d_ov = float(np.max(np.abs(seams_ql - ov_ql)))
                        d_fr = float(np.max(np.abs(seams_ql - fr_ql)))
                        for code, samples, cores in (
                            ("seams", seams, str(threads)),
                            ("ovito", ov, ov_cores),
                            ("freud", fr, fr_cores),
                        ):
                            for value in samples:
                                rows.append(
                                    f"{code},{kind},{n},{order},{threads},{cores},{value:.6f}")
                        sm, si = summarize(seams)
                        om, oi = summarize(ov)
                        fm, fi = summarize(fr)
                        print(
                            f"{kind} n={n} l={order} t={threads} "
                            f"seams {sm:.3f} ({si:.3f}) ovito {om:.3f} ({oi:.3f}) "
                            f"freud {fm:.3f} ({fi:.3f}) dOv {d_ov:.3e} dFr {d_fr:.3e}",
                            flush=True,
                        )
    Path(args.out).write_text("\n".join(rows) + "\n", encoding="utf-8")


if __name__ == "__main__":
    main()
