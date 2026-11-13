#!/usr/bin/env python3
"""Time 1 mol/kg NaCl, 50k oxygens, at 1 and 8 threads.

Release -O3. One warm-up, then six timed calls. The median of an even
sample is the mean of the two central values. Quartiles are Hyndman-Fan
type 7. The timed region is CHILL+ or the neighbour list on the water
oxygens (LAMMPS type 1). Ions are not in that set.

One thread is pinned to core 0. Eight threads use cores 0-7. d-SEAMS
takes its width from OMP_NUM_THREADS and gives that width to Rayon.
mdapy and freud set their own pools. OVITO reads OVITO_THREAD_COUNT
before import, so each width is a new process.
"""

from __future__ import annotations

import argparse
import csv
import json
import os
import platform
import subprocess
import sys
import tempfile
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
FRAME = ROOT / "input" / "brine" / "work" / "nacl-n50000-m1.0" / "dump.lammpstrj"
LIB = ROOT / "bbdir" / "src" / "libyodaLib.so"
CUTOFF = 3.5
WARMUP = 1
REPEATS = 6
OXYGENS = 50000
WIDTHS = (1, 8)

DSEAMS_COARSE = {
    0: "cubic",
    1: "hexagonal",
    2: "other",
    3: "interfacial",
    4: "hydrate",
    5: "interfacial_hydrate",
    6: "other",
    7: "cubic",
    8: "hexagonal",
}
MDAPY_COARSE = {
    0: "other",
    1: "hexagonal",
    2: "cubic",
    3: "interfacial",
    4: "hydrate",
    5: "interfacial_hydrate",
}


def quantile_linear(sorted_vals, q):
    pos = q * (len(sorted_vals) - 1)
    lo = int(pos)
    frac = pos - lo
    if lo + 1 >= len(sorted_vals):
        return sorted_vals[-1]
    return sorted_vals[lo] * (1.0 - frac) + sorted_vals[lo + 1] * frac


def summarize(samples):
    if len(samples) < REPEATS:
        raise ValueError(f"need n >= {REPEATS}, got {len(samples)}")
    ordered = sorted(float(v) for v in samples)
    mid = len(ordered) // 2
    if len(ordered) % 2 == 0:
        median = 0.5 * (ordered[mid - 1] + ordered[mid])
    else:
        median = ordered[mid]
    q1 = quantile_linear(ordered, 0.25)
    q3 = quantile_linear(ordered, 0.75)
    return median, q1, q3, q3 - q1


def hardware():
    cpu = platform.processor() or platform.machine()
    try:
        for line in Path("/proc/cpuinfo").read_text().splitlines():
            if line.startswith("model name"):
                cpu = line.split(":", 1)[1].strip()
                break
    except OSError:
        pass
    nproc = subprocess.check_output(["nproc"], text=True).strip()
    nproc_all = subprocess.check_output(["nproc", "--all"], text=True).strip()
    allowed = "unknown"
    try:
        for line in Path("/proc/self/status").read_text().splitlines():
            if line.startswith("Cpus_allowed_list:"):
                allowed = line.split(":", 1)[1].strip()
    except OSError:
        pass
    return f"{cpu}; nproc={nproc}; nproc_all={nproc_all}; affinity={allowed}"


def read_dump(path, oxygen_type=1):
    oxy = []
    box = None
    with open(path) as handle:
        mode = None
        bounds = []
        for line in handle:
            text = line.strip()
            if text.startswith("ITEM: BOX"):
                mode = "box"
                bounds = []
                continue
            if text.startswith("ITEM: ATOMS"):
                mode = "atoms"
                continue
            if text.startswith("ITEM:"):
                mode = None
                continue
            if mode == "box":
                parts = text.split()
                bounds.append(float(parts[1]) - float(parts[0]))
                if len(bounds) == 3:
                    box = bounds
            elif mode == "atoms":
                parts = text.split()
                if int(float(parts[2])) == oxygen_type:
                    oxy.append((float(parts[3]), float(parts[4]), float(parts[5])))
    if not oxy or box is None:
        raise RuntimeError(f"no oxygens in {path}")
    import numpy as np

    return np.asarray(oxy, dtype=np.float64), box


def time_calls(run):
    import time

    samples = []
    last = None
    for i in range(WARMUP + REPEATS):
        t0 = time.perf_counter_ns()
        last = run()
        dt = (time.perf_counter_ns() - t0) / 1e6
        if i >= WARMUP:
            samples.append(dt)
    return samples, last


def mask_for(threads):
    return "0" if threads == 1 else "0-7"


def base_env():
    env = os.environ.copy()
    env.pop("LD_PRELOAD", None)
    env["OMP_DYNAMIC"] = "FALSE"
    env["OMP_PROC_BIND"] = "close"
    env["OPENBLAS_NUM_THREADS"] = "1"
    env["MKL_NUM_THREADS"] = "1"
    env["NUMEXPR_NUM_THREADS"] = "1"
    return env


def worker_dseams(threads):
    import ctypes

    import numpy as np

    points, box = read_dump(FRAME)
    if points.shape[0] != OXYGENS:
        raise SystemExit(f"expected {OXYGENS} oxygens, got {points.shape[0]}")
    n = points.shape[0]
    lib = ctypes.CDLL(str(LIB))
    lib.seams_chill_plus_species.argtypes = [
        ctypes.POINTER(ctypes.c_double),
        ctypes.c_int,
        ctypes.POINTER(ctypes.c_double),
        ctypes.POINTER(ctypes.c_int),
        ctypes.POINTER(ctypes.c_int),
        ctypes.c_int,
        ctypes.c_int,
        ctypes.POINTER(ctypes.c_int),
        ctypes.c_int,
        ctypes.c_double,
        ctypes.c_char_p,
        ctypes.POINTER(ctypes.c_int),
        ctypes.POINTER(ctypes.c_int),
        ctypes.POINTER(ctypes.c_int),
        ctypes.POINTER(ctypes.c_double),
    ]
    lib.seams_chill_plus_species.restype = ctypes.c_int
    lib.seams_oxygen_cutoff_degree.argtypes = [
        ctypes.POINTER(ctypes.c_double),
        ctypes.c_int,
        ctypes.POINTER(ctypes.c_double),
        ctypes.POINTER(ctypes.c_int),
        ctypes.POINTER(ctypes.c_int),
        ctypes.c_int,
        ctypes.c_int,
        ctypes.POINTER(ctypes.c_int),
        ctypes.c_int,
        ctypes.c_double,
        ctypes.POINTER(ctypes.c_int),
    ]
    lib.seams_oxygen_cutoff_degree.restype = ctypes.c_int
    lib.seams_c_error.restype = ctypes.c_char_p
    lib.lc_num_threads.restype = ctypes.c_int
    xyz = np.ascontiguousarray(points.reshape(-1))
    box9 = np.array([0, 0, 0, box[0], box[1], box[2], 0, 0, 0], dtype=np.float64)
    types = np.ones(n, dtype=np.int32)
    labels = np.empty(n, dtype=np.int32)
    degree = np.empty(n, dtype=np.int32)

    def run_chill():
        rc = lib.seams_chill_plus_species(
            xyz.ctypes.data_as(ctypes.POINTER(ctypes.c_double)),
            n,
            box9.ctypes.data_as(ctypes.POINTER(ctypes.c_double)),
            types.ctypes.data_as(ctypes.POINTER(ctypes.c_int)),
            None,
            1,
            0,
            None,
            0,
            CUTOFF,
            b"cutoff",
            labels.ctypes.data_as(ctypes.POINTER(ctypes.c_int)),
            None,
            None,
            None,
        )
        if rc != 0:
            raise RuntimeError(ctypes.string_at(lib.seams_c_error()).decode())
        return labels.copy()

    def run_neigh():
        rc = lib.seams_oxygen_cutoff_degree(
            xyz.ctypes.data_as(ctypes.POINTER(ctypes.c_double)),
            n,
            box9.ctypes.data_as(ctypes.POINTER(ctypes.c_double)),
            types.ctypes.data_as(ctypes.POINTER(ctypes.c_int)),
            None,
            1,
            0,
            None,
            0,
            CUTOFF,
            degree.ctypes.data_as(ctypes.POINTER(ctypes.c_int)),
        )
        if rc != 0:
            raise RuntimeError(ctypes.string_at(lib.seams_c_error()).decode())
        return int(degree.sum())

    chill, chill_labels = time_calls(run_chill)
    neigh, degree_sum = time_calls(run_neigh)
    reported = int(lib.lc_num_threads())
    if reported != threads:
        raise SystemExit(f"linkcell threads {reported} != {threads}")
    out = label_path("dseams", threads)
    np.save(out, chill_labels)
    emit("d-SEAMS", "2.13.0", threads, reported, [
        {"operation": "chill+", "samples": chill, "labels": str(out)},
        {"operation": "neighbours", "samples": neigh, "degree_sum": degree_sum},
    ])


def worker_mdapy(threads):
    import mdapy
    import numpy as np
    from mdapy.identify import ChillPlus
    from mdapy.neighbor import NeighborList

    mdapy.set_num_threads(threads)
    reported = int(mdapy.get_num_threads())
    if reported != threads:
        raise SystemExit(f"mdapy threads {reported} != {threads}")
    points, box = read_dump(FRAME)
    if points.shape[0] != OXYGENS:
        raise SystemExit(f"expected {OXYGENS} oxygens, got {points.shape[0]}")
    system = mdapy.System(points, box)

    def run_chill():
        return np.asarray(ChillPlus(CUTOFF).compute(system).structure).copy()

    def run_neigh():
        return int(NeighborList.build(system, CUTOFF).j.size)

    chill, labels = time_calls(run_chill)
    neigh, pairs = time_calls(run_neigh)
    out = label_path("mdapy", threads)
    np.save(out, labels.astype(np.int32))
    emit("mdapy", mdapy.__version__, threads, reported, [
        {"operation": "chill+", "samples": chill, "labels": str(out)},
        {"operation": "neighbours", "samples": neigh, "pairs": pairs},
    ])


def worker_ovito(threads):
    import ovito
    import numpy as np
    from ovito.data import DataCollection, Particles, SimulationCell
    from ovito.modifiers import ChillPlusModifier
    from ovito.pipeline import Pipeline, StaticSource

    if tuple(ovito.version)[:3] != (3, 17, 2):
        raise SystemExit(f"ovito {ovito.version}")
    points, box = read_dump(FRAME)
    if points.shape[0] != OXYGENS:
        raise SystemExit(f"expected {OXYGENS} oxygens, got {points.shape[0]}")

    def run_chill():
        data = DataCollection()
        particles = Particles()
        particles.create_property("Position", data=points)
        data.objects.append(particles)
        cell = SimulationCell(pbc=(True, True, True))
        cell[:, 0] = (box[0], 0.0, 0.0)
        cell[:, 1] = (0.0, box[1], 0.0)
        cell[:, 2] = (0.0, 0.0, box[2])
        data.objects.append(cell)
        pipe = Pipeline(source=StaticSource(data=data))
        pipe.modifiers.append(ChillPlusModifier(cutoff=CUTOFF))
        return np.asarray(pipe.compute().particles["Structure Type"]).copy()

    chill, labels = time_calls(run_chill)
    out = label_path("ovito", threads)
    np.save(out, labels.astype(np.int32))
    emit("OVITO", ".".join(str(v) for v in tuple(ovito.version)[:3]), threads, threads, [
        {"operation": "chill+", "samples": chill, "labels": str(out)},
    ])


def worker_freud(threads):
    import freud

    freud.parallel.set_num_threads(threads)
    reported = int(freud.parallel.get_num_threads())
    if reported != threads:
        raise SystemExit(f"freud threads {reported} != {threads}")
    points, box = read_dump(FRAME)
    if points.shape[0] != OXYGENS:
        raise SystemExit(f"expected {OXYGENS} oxygens, got {points.shape[0]}")
    fbox = freud.box.Box(box[0], box[1], box[2])

    def run_neigh():
        query = freud.locality.AABBQuery(fbox, points)
        return int(
            query.query(points, {"r_max": CUTOFF, "exclude_ii": True})
            .toNeighborList()
            .neighbor_counts.sum()
        )

    neigh, counts = time_calls(run_neigh)
    emit("freud", freud.__version__, threads, reported, [
        {"operation": "neighbours", "samples": neigh, "degree_sum": counts},
    ])


def emit(code, version, threads, reported, results):
    payload = {
        "code": code,
        "version": version,
        "threads": threads,
        "reported_threads": reported,
        "results": results,
    }
    sys.stdout.write(json.dumps(payload) + "\n")


def coarse(labels, table):
    import numpy as np

    return np.array([table[int(v)] for v in labels], dtype=object)


def agreement(a, b):
    import numpy as np

    if len(a) != len(b) or len(a) == 0:
        return float("nan")
    return float(np.mean(a == b))


def label_path(tool, threads):
    root = os.environ.get("BRINE_LABEL_DIR")
    if not root:
        root = tempfile.mkdtemp(prefix="brine-labels-")
        os.environ["BRINE_LABEL_DIR"] = root
    return Path(root) / f"brine-labels-{tool}-{threads}.npy"


def run_worker(role, threads, python, labels):
    env = base_env()
    env["OMP_NUM_THREADS"] = str(threads)
    env["RAYON_NUM_THREADS"] = str(threads)
    env["BRINE_LABEL_DIR"] = str(labels)
    if role == "dseams":
        env.pop("LD_LIBRARY_PATH", None)
    elif role == "ovito":
        env["LD_LIBRARY_PATH"] = "/usr/lib/x86_64-linux-gnu"
        env["QT_QPA_PLATFORM"] = "offscreen"
        env["OVITO_THREAD_COUNT"] = str(threads)
    else:
        env.pop("LD_LIBRARY_PATH", None)
    cmd = [
        "taskset", "-c", mask_for(threads),
        str(python), str(Path(__file__).resolve()),
        "--worker", role, "--threads", str(threads),
    ]
    done = subprocess.run(cmd, env=env, text=True, capture_output=True, check=False)
    if done.returncode != 0:
        raise RuntimeError(f"{role} threads={threads} failed\n{done.stderr}\n{done.stdout}")
    lines = [ln for ln in done.stdout.splitlines() if ln.strip().startswith("{")]
    if len(lines) != 1:
        raise RuntimeError(f"{role} threads={threads} output\n{done.stdout}\n{done.stderr}")
    return json.loads(lines[0])


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--worker", choices=("dseams", "mdapy", "ovito", "freud"))
    parser.add_argument("--threads", type=int)
    parser.add_argument(
        "--python",
        default=os.environ.get("BRINE_PYTHON", sys.executable),
        help="Interpreter for the worker processes. BRINE_PYTHON overrides the default.",
    )
    parser.add_argument(
        "--label-dir",
        default=os.environ.get("BRINE_LABEL_DIR"),
        help="Directory for the label arrays. BRINE_LABEL_DIR overrides the default.",
    )
    args = parser.parse_args()
    if args.worker:
        os.environ["OMP_NUM_THREADS"] = str(args.threads)
        os.environ["RAYON_NUM_THREADS"] = str(args.threads)
        {
            "dseams": worker_dseams,
            "mdapy": worker_mdapy,
            "ovito": worker_ovito,
            "freud": worker_freud,
        }[args.worker](args.threads)
        return 0

    import numpy as np

    labels = Path(args.label_dir) if args.label_dir else Path(tempfile.mkdtemp(prefix="brine-labels-"))
    labels.mkdir(parents=True, exist_ok=True)
    hw = hardware()
    csv_path = ROOT / "benchmarks" / "h2h-brine-50k-threads.csv"
    log_path = ROOT / "benchmarks" / "h2h-brine-50k-threads.txt"
    fields = [
        "hardware", "build", "frame", "code", "version", "threads", "affinity",
        "operation", "oxygens", "warmup", "n", "median_ms", "q1_ms", "q3_ms",
        "iqr_ms", "agreement", "samples_ms",
    ]
    rows = []
    log_lines = [
        f"frame nacl-n50000-m1.0",
        f"oxygens {OXYGENS}",
        f"build release -O3",
        f"hardware {hw}",
        f"warmup {WARMUP}",
        f"n {REPEATS}",
        f"cutoff {CUTOFF}",
        "one thread pinned to core 0; eight threads on cores 0-7",
    ]
    ref = None
    ref_raw = None
    for threads in WIDTHS:
        for role in ("dseams", "mdapy", "ovito", "freud"):
            print(f"{role} threads={threads}", flush=True)
            payload = run_worker(role, threads, args.python, labels)
            for result in payload["results"]:
                samples = [float(v) for v in result["samples"]]
                med, q1, q3, iqr = summarize(samples)
                agree = ""
                if "labels" in result:
                    labels = np.load(result["labels"])
                    if payload["code"] == "d-SEAMS":
                        if ref_raw is None:
                            ref_raw = labels
                            ref = coarse(labels, DSEAMS_COARSE)
                        else:
                            if not np.array_equal(labels, ref_raw):
                                raise SystemExit(
                                    f"d-SEAMS labels at {threads} threads differ from 1 thread"
                                )
                        agree = "1.0000"
                    else:
                        other = coarse(labels, MDAPY_COARSE)
                        agree = f"{agreement(ref, other):.4f}"
                affinity = mask_for(threads)
                rows.append({
                    "hardware": hw,
                    "build": "release -O3",
                    "frame": "nacl-n50000-m1.0",
                    "code": payload["code"],
                    "version": payload["version"],
                    "threads": threads,
                    "affinity": affinity,
                    "operation": result["operation"],
                    "oxygens": OXYGENS,
                    "warmup": WARMUP,
                    "n": len(samples),
                    "median_ms": f"{med:.6f}",
                    "q1_ms": f"{q1:.6f}",
                    "q3_ms": f"{q3:.6f}",
                    "iqr_ms": f"{iqr:.6f}",
                    "agreement": agree,
                    "samples_ms": ",".join(f"{v:.6f}" for v in samples),
                })
                log_lines.append(
                    f"threads {threads} affinity {affinity} code {payload['code']} "
                    f"version {payload['version']} operation {result['operation']} "
                    f"reported_threads {payload['reported_threads']}"
                )
                for i, sample in enumerate(samples, start=1):
                    log_lines.append(f"sample {i} {sample:.6f}")
                extra = ""
                if "degree_sum" in result:
                    extra = f" degree_sum {result['degree_sum']}"
                if "pairs" in result:
                    extra = f" pairs {result['pairs']}"
                log_lines.append(
                    f"median {med:.6f} q1 {q1:.6f} q3 {q3:.6f} iqr {iqr:.6f} "
                    f"agreement {agree or 'na'}{extra}"
                )
                print(
                    f"  {payload['code']} {result['operation']} "
                    f"{med:.3f} ms iqr {iqr:.3f} agree {agree or 'na'}",
                    flush=True,
                )
    with csv_path.open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields)
        writer.writeheader()
        writer.writerows(rows)
    log_path.write_text("\n".join(log_lines) + "\n")
    print(f"wrote {csv_path}")
    print(f"wrote {log_path}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
