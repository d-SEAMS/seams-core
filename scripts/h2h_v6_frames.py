#!/usr/bin/env python3
"""Labels and times on the 200 ps mW nucleation frames.

Kernel classifyIce, OVITO PTM with the diamond templates, OVITO identify
diamond, and pyscal3 diamond_structure. One warm-up, six calls.
"""

from __future__ import annotations

import csv
import json
import os
import subprocess
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))

import h2h_protocol as proto

ROOT = Path(__file__).resolve().parents[1]
HERE = Path(__file__).resolve()
HOST_LIB = "/usr/lib/x86_64-linux-gnu"
KERNEL = os.environ.get("SYMMETRY_BENCH", "bench_symmetry_ice")
SHA = os.environ.get("SYMMETRY_SHA", "59d85155350cfa3c53da3df11d6ddf45b6b31ede")


def read_xyz(path):
    import numpy as np
    lines = Path(path).read_text().splitlines()
    n = int(lines[0])
    side = float(lines[1].split()[1])
    pos = np.zeros((n, 3))
    for i, line in enumerate(lines[2:2 + n]):
        _el, x, y, z = line.split()
        pos[i] = (float(x), float(y), float(z))
    return pos, side


def kernel(path, threads):
    cores = "0" if threads == 1 else "0-7"
    env = os.environ.copy()
    env["OMP_NUM_THREADS"] = str(threads)
    env["OMP_DYNAMIC"] = "FALSE"
    env.pop("LD_LIBRARY_PATH", None)
    done = subprocess.run(
        ["taskset", "-c", cores, KERNEL, path, "1", "6"],
        env=env, text=True, capture_output=True, check=False)
    if done.returncode != 0:
        tail = (done.stderr or done.stdout or "").strip().splitlines()
        raise SystemExit(tail[-1] if tail else f"kernel exit {done.returncode}")
    text = done.stdout.strip().splitlines()[-1]
    fields = {}
    for part in text.split():
        if "=" in part and not part.startswith("samples_ms"):
            key, value = part.split("=", 1)
            fields[key] = value
    samples = [float(v) for v in text.split("samples_ms=", 1)[1].split()]
    return samples, fields


def python_worker(kind, threads, frame_npz):
    code = r'''
import json, os, sys
from pathlib import Path
sys.path.insert(0, sys.argv[1])
import numpy as np
import h2h_protocol as proto
kind, threads, frame = sys.argv[2], int(sys.argv[3]), sys.argv[4]
data = np.load(frame)
pos, side = data["positions"], float(data["side"])
n = len(pos)
cpus = "unknown"
for line in open("/proc/self/status", encoding="utf-8"):
    if line.startswith("Cpus_allowed_list:"):
        cpus = line.split(":", 1)[1].strip()
if kind == "pyscal":
    import pyscal3
    from ase import Atoms
    pyscal3.set_num_threads(threads)
    atoms = Atoms(symbols=["O"] * n, positions=pos, cell=[side, side, side], pbc=True)
    def call(_a):
        return dict(pyscal3.diamond_structure(atoms))
    samples, result = proto.time_calls(lambda: atoms, call)
    json.dump({"samples_ms": samples, "labels": {str(k): int(v) for k, v in result.items()}, "cpus": cpus}, sys.stdout)
else:
    os.environ["OVITO_THREAD_COUNT"] = str(threads)
    os.environ["QT_QPA_PLATFORM"] = "offscreen"
    from ovito.data import DataCollection
    from ovito.modifiers import IdentifyDiamondModifier, PolyhedralTemplateMatchingModifier
    def prep():
        data = DataCollection()
        data.create_cell(matrix=[[side, 0, 0, 0], [0, side, 0, 0], [0, 0, side, 0]], pbc=(True, True, True))
        particles = data.create_particles(count=n)
        particles.create_property("Position", data=np.ascontiguousarray(pos))
        return data
    if kind == "ptm":
        def call(data):
            mod = PolyhedralTemplateMatchingModifier()
            for structure in mod.structures:
                structure.enabled = structure.id in (6, 7)
            data.apply(mod)
            labels = np.asarray(data.particles["Structure Type"], dtype=np.int32)
            values, counts = np.unique(labels, return_counts=True)
            return {str(int(v)): int(c) for v, c in zip(values, counts)}
    else:
        def call(data):
            data.apply(IdentifyDiamondModifier())
            labels = np.asarray(data.particles["Structure Type"], dtype=np.int32)
            values, counts = np.unique(labels, return_counts=True)
            return {str(int(v)): int(c) for v, c in zip(values, counts)}
    samples, result = proto.time_calls(prep, call)
    json.dump({"samples_ms": samples, "labels": result, "cpus": cpus}, sys.stdout)
'''
    cores = "0" if threads == 1 else "0-7"
    env = os.environ.copy()
    env["OMP_NUM_THREADS"] = str(threads)
    env["OMP_DYNAMIC"] = "FALSE"
    if kind == "pyscal":
        env.pop("LD_LIBRARY_PATH", None)
    else:
        env["LD_LIBRARY_PATH"] = HOST_LIB
        env["OVITO_THREAD_COUNT"] = str(threads)
        env["QT_QPA_PLATFORM"] = "offscreen"
    done = subprocess.run(
        ["taskset", "-c", cores, sys.executable, "-c", code, str(HERE.parent), kind, str(threads), frame_npz],
        env=env, text=True, capture_output=True, check=False)
    if done.returncode != 0:
        tail = (done.stderr or done.stdout or "").strip().splitlines()
        raise SystemExit(tail[-1] if tail else f"{kind} exit {done.returncode}")
    line = [row for row in done.stdout.splitlines() if row.startswith("{")][-1]
    return json.loads(line)


def main():
    import numpy as np
    frames = sys.argv[1:]
    if not frames:
        raise SystemExit("usage: h2h_v6_frames.py seed-91001-200ps.xyz ...")
    out = ROOT / "benchmarks" / "h2h-watch-v6-frames.csv"
    fields = ["frame", "threads", "code", "operation", "molecules", "n",
              "median_ms", "iqr_ms", "detail", "samples_ms"]
    with out.open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields)
        writer.writeheader()
        for path in frames:
            pos, side = read_xyz(path)
            npz = Path(path).with_suffix(".npz")
            np.savez(npz, positions=pos, side=side)
            name = Path(path).stem
            for threads in (1, 8):
                samples, fields_k = kernel(path, threads)
                summary = proto.summarize(samples)
                writer.writerow({
                    "frame": name, "threads": threads, "code": "d-SEAMS",
                    "operation": "classify-ice", "molecules": len(pos),
                    "n": summary["n"],
                    "median_ms": f"{summary['median_ms']:.3f}",
                    "iqr_ms": f"{summary['iqr_ms']:.3f}",
                    "detail": f"sha={SHA} Ih={fields_k.get('Ih')} Ic={fields_k.get('Ic')} Isd={fields_k.get('Isd')} other={fields_k.get('other')}",
                    "samples_ms": " ".join(f"{v:.3f}" for v in samples),
                })
                handle.flush()
                print(name, threads, "kernel", summary["median_ms"], flush=True)
                for kind, code, operation in (
                    ("ptm", "OVITO", "ptm"),
                    ("diamond", "OVITO", "diamond-cna"),
                    ("pyscal", "pyscal3", "diamond-cna"),
                ):
                    payload = python_worker(kind, threads, str(npz))
                    summary = proto.summarize(payload["samples_ms"])
                    writer.writerow({
                        "frame": name, "threads": threads, "code": code,
                        "operation": operation, "molecules": len(pos),
                        "n": summary["n"],
                        "median_ms": f"{summary['median_ms']:.3f}",
                        "iqr_ms": f"{summary['iqr_ms']:.3f}",
                        "detail": json.dumps(payload["labels"]),
                        "samples_ms": " ".join(f"{v:.3f}" for v in payload["samples_ms"]),
                    })
                    handle.flush()
                    print(name, threads, operation, summary["median_ms"], payload["labels"], flush=True)
            npz.unlink(missing_ok=True)


if __name__ == "__main__":
    main()
