#!/usr/bin/env python3
"""Second head-to-head: OVITO 3.17.2, matscipy 1.3.1 and cycless 0.7.

Same protocol as scripts/h2h_watch.py. Ideal ice Ic, one core, one
thread, one warm-up, six timed calls, median and type-7 IQR. pyscal3
has neither CHILL+ nor a ring census, and the Palmer CHILL+ sources at
ee8a3cb do not link, so those two are not timed.

    H2H_PYTHON=.../python scripts/h2h_v2.py \
        --lib bbdir-eessi/bindings/libseams_lang.so \
        --out benchmarks/h2h-watch-v2.csv
"""

from __future__ import annotations

import argparse
import csv
import json
import os
import subprocess
import sys
import tempfile
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))

import h2h_protocol as proto

proto.apply_pin_env()
os.environ["OVITO_THREAD_COUNT"] = "1"
os.environ.setdefault(
    "LD_LIBRARY_PATH",
    "/usr/lib/x86_64-linux-gnu" + os.pathsep + os.environ.get("LD_LIBRARY_PATH", ""),
)

import numpy as np

ROOT = Path(__file__).resolve().parents[1]
HERE = Path(__file__).resolve()

# OVITO ChillPlusModifier.Type on this build.
OVITO_CLASS = {
    0: "other",
    1: "hexagonal",
    2: "cubic",
    3: "interfacial",
    4: "hydrate",
    5: "interfacial_hydrate",
}


def launch_env():
    env = proto.pin_env()
    env["OVITO_THREAD_COUNT"] = "1"
    lib = "/usr/lib/x86_64-linux-gnu"
    env["LD_LIBRARY_PATH"] = lib + os.pathsep + env.get("LD_LIBRARY_PATH", "")
    return env


def pinned(python, args):
    cmd = ["taskset", "-c", proto.CORE, python, str(HERE), *args]
    done = subprocess.run(
        cmd, cwd=ROOT, env=launch_env(), text=True, capture_output=True)
    if done.returncode != 0:
        tail = (done.stderr or done.stdout or "").strip().splitlines()
        detail = tail[-1] if tail else f"exit {done.returncode}"
        return {"error": detail, "returncode": done.returncode}
    line = done.stdout.strip().splitlines()[-1]
    payload = json.loads(line)
    if payload.get("cpus") != proto.CORE:
        return {"error": f"not pinned: {payload.get('cpus')}"}
    return payload


def emit(payload):
    json.dump(payload, sys.stdout)
    sys.stdout.write("\n")


def load_frame(path):
    data = np.load(path)
    return data["positions"], float(data["side"])


def worker_ovito(frame, operation, warmup, repeats):
    from ovito.data import DataCollection
    from ovito.modifiers import (
        ChillPlusModifier, CreateBondsModifier, FindRingsModifier)
    import ovito

    pos, side = load_frame(frame)

    def prep():
        data = DataCollection()
        data.create_cell(
            matrix=[[side, 0, 0, 0], [0, side, 0, 0], [0, 0, side, 0]],
            pbc=(True, True, True))
        particles = data.create_particles(count=len(pos))
        particles.create_property("Position", data=np.ascontiguousarray(pos))
        return data

    if operation == "chill+":
        def call(data):
            data.apply(ChillPlusModifier(cutoff=proto.CUTOFF_A))
            labels = np.asarray(data.particles["Structure Type"], dtype=np.int32)
            return [OVITO_CLASS[int(v)] for v in labels]
    else:
        def call(data):
            data.apply(CreateBondsModifier(cutoff=proto.CUTOFF_A))
            data.apply(FindRingsModifier(
                minimum_ring_size=3, maximum_ring_size=proto.MAX_RING))
            hist = data.tables["ring-size-histogram"]
            sizes = [int(v) for v in np.asarray(hist["Ring Size"])]
            counts = [int(v) for v in np.asarray(hist["Count"])]
            found = dict(zip(sizes, counts))
            return [found.get(size, 0) for size in range(3, proto.MAX_RING + 1)]

    samples, result = proto.time_calls(prep, call, warmup, repeats)
    payload = {
        "code": "OVITO",
        "version": ovito.version_string,
        "operation": operation,
        "molecules": int(len(pos)),
        "samples_ms": samples,
        "cpus": proto.cpus_allowed(),
    }
    if operation == "chill+":
        payload["labels"] = result
    else:
        payload["census"] = result
    emit(payload)


def worker_matscipy(frame, warmup, repeats):
    import matscipy
    from ase import Atoms
    from matscipy.rings import ring_statistics

    pos, side = load_frame(frame)
    atoms = Atoms(
        symbols=["O"] * len(pos), positions=pos,
        cell=[side, side, side], pbc=True)

    def prep():
        return atoms

    def call(frame_atoms):
        # maxlength is one past the largest ring, as in the paper's campaign.
        counts = np.asarray(ring_statistics(
            frame_atoms, proto.CUTOFF_A, maxlength=proto.MAX_RING + 1))
        return [int(counts[size]) if size < len(counts) else 0
                for size in range(3, proto.MAX_RING + 1)]

    samples, census = proto.time_calls(prep, call, warmup, repeats)
    emit({
        "code": "matscipy",
        "version": matscipy.__version__,
        "operation": "rings",
        "molecules": int(len(pos)),
        "samples_ms": samples,
        "census": census,
        "cpus": proto.cpus_allowed(),
    })


def worker_cycless(frame, warmup, repeats):
    import importlib.metadata
    import networkx as nx
    import pairlist
    from cycless.cycles import cycles_iter

    pos, side = load_frame(frame)
    version = importlib.metadata.version("cycless")

    def prep():
        return pos

    def call(points):
        fractional = points / side
        graph = nx.Graph()
        graph.add_nodes_from(range(len(points)))
        graph.add_edges_from(pairlist.pairs_iter(
            fractional, proto.CUTOFF_A, np.diag([side, side, side]),
            distance=False))
        census = [0, 0, 0, 0]
        for ring in cycles_iter(graph, proto.MAX_RING, pos=fractional):
            size = len(ring)
            if 3 <= size <= proto.MAX_RING:
                census[size - 3] += 1
        return census

    samples, census = proto.time_calls(prep, call, warmup, repeats)
    emit({
        "code": "cycless",
        "version": version,
        "operation": "rings",
        "molecules": int(len(pos)),
        "samples_ms": samples,
        "census": census,
        "cpus": proto.cpus_allowed(),
    })


def agreement(reference, other, operation):
    if "error" in other or "error" in reference:
        return "", other.get("error") or reference.get("error")
    if operation == "chill+":
        left, right = reference["labels"], other["labels"]
        matched = sum(a == b for a, b in zip(left, right))
        return f"{matched / len(left):.6f}", f"matched={matched}/{len(left)}"
    same = reference["census"] == other["census"]
    return f"{float(same):.6f}", f"dseams={reference['census']};other={other['census']}"


def write_pair(writer, reference, other, operation, repetitions):
    fraction, detail = agreement(reference, other, operation)
    molecules = int(reference.get("molecules") or other.get("molecules") or 0)
    failed = "error" in reference or "error" in other
    for payload in (reference, other):
        if "error" in payload:
            writer.writerow({
                "code": payload.get("code", "unknown"),
                "version": payload.get("version", ""),
                "operation": operation,
                "molecules": molecules,
                "repetitions": repetitions,
                "warmup": proto.WARMUP,
                "n": 0,
                "median_ms": "",
                "q1_ms": "",
                "q3_ms": "",
                "iqr_ms": "",
                "agreement": "",
                "detail": payload["error"],
                "samples_ms": "",
            })
            continue
        summary = proto.summarize(payload["samples_ms"])
        row_detail = detail
        row_fraction = fraction
        if failed:
            row_fraction = ""
            row_detail = "timed; the other code aborted: " + detail
        writer.writerow({
            "code": payload["code"],
            "version": payload["version"],
            "operation": operation,
            "molecules": payload["molecules"],
            "repetitions": repetitions,
            "warmup": proto.WARMUP,
            "n": summary["n"],
            "median_ms": f"{summary['median_ms']:.6f}",
            "q1_ms": f"{summary['q1_ms']:.6f}",
            "q3_ms": f"{summary['q3_ms']:.6f}",
            "iqr_ms": f"{summary['iqr_ms']:.6f}",
            "agreement": row_fraction,
            "detail": row_detail,
            "samples_ms": ",".join(f"{v:.6f}" for v in payload["samples_ms"]),
        })


def dseams(python, lib, frame, operation):
    cmd = [
        "--worker", "--lib", lib, "--operation", operation,
        "--warmup", str(proto.WARMUP), "--repeats", str(proto.REPEATS),
        "--frame", frame,
    ]
    # Reuse the v1 worker so the d-SEAMS call is the same binary path.
    full = ["taskset", "-c", proto.CORE, python,
            str(ROOT / "scripts/h2h_watch.py"), *cmd]
    done = subprocess.run(
        full, cwd=ROOT, env=launch_env(), text=True, capture_output=True)
    if done.returncode != 0:
        tail = (done.stderr or done.stdout or "").strip().splitlines()
        return {"error": tail[-1] if tail else f"exit {done.returncode}",
                "code": "d-SEAMS"}
    payload = json.loads(done.stdout.strip().splitlines()[-1])
    if payload.get("cpus") != proto.CORE:
        return {"error": f"not pinned: {payload.get('cpus')}", "code": "d-SEAMS"}
    return payload


def campaign(args):
    python = args.python or sys.executable
    out = Path(args.out)
    out.parent.mkdir(parents=True, exist_ok=True)
    fields = [
        "code", "version", "operation", "molecules", "repetitions",
        "warmup", "n", "median_ms", "q1_ms", "q3_ms", "iqr_ms",
        "agreement", "detail", "samples_ms",
    ]
    sizes = [int(v) for v in args.sizes.split(",")] if args.sizes else list(proto.REPETITIONS)
    with out.open("w", newline="") as handle, tempfile.TemporaryDirectory(prefix="h2h-v2-") as tmp:
        writer = csv.DictWriter(handle, fieldnames=fields)
        writer.writeheader()
        for repetitions in sizes:
            pos, side = proto.ice_ic(repetitions)
            frame = str(Path(tmp) / f"ic-{len(pos)}.npz")
            np.savez(frame, positions=pos, side=np.float64(side))
            print(f"frame molecules={len(pos)}", flush=True)
            jobs = [
                ("chill+", "ovito", ["--worker", "ovito", "--operation", "chill+"]),
                ("rings", "ovito", ["--worker", "ovito", "--operation", "rings"]),
                ("rings", "matscipy", ["--worker", "matscipy"]),
                ("rings", "cycless", ["--worker", "cycless"]),
            ]
            # One d-SEAMS sample set per operation, shared by the ring codes.
            reference = {}
            for operation, _name, _extra in jobs:
                if operation not in reference:
                    reference[operation] = dseams(python, args.lib, frame, operation)
            for operation, name, extra in jobs:
                other = pinned(python, [
                    *extra, "--frame", frame,
                    "--warmup", str(proto.WARMUP),
                    "--repeats", str(proto.REPEATS),
                ])
                if "error" not in other:
                    other.setdefault("code", name)
                else:
                    other["code"] = name
                write_pair(writer, reference[operation], other, operation, repetitions)
                handle.flush()
                fraction, detail = agreement(reference[operation], other, operation)
                if "samples_ms" in other and "samples_ms" in reference[operation]:
                    left = proto.summarize(reference[operation]["samples_ms"])
                    right = proto.summarize(other["samples_ms"])
                    print(
                        f"  {operation} {name} d-SEAMS {left['median_ms']:.3f} "
                        f"other {right['median_ms']:.3f} iqr {right['iqr_ms']:.3f} "
                        f"agreement {fraction} {detail}",
                        flush=True)
                else:
                    print(f"  {operation} {name} FAILED {detail}", flush=True)
    print(f"wrote {out}", flush=True)


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--lib", default=str(ROOT / "bbdir-eessi/bindings/libseams_lang.so"))
    parser.add_argument("--out", default=str(ROOT / "benchmarks/h2h-watch-v2.csv"))
    parser.add_argument("--python", default="")
    parser.add_argument("--sizes", default="")
    parser.add_argument("--worker", default="")
    parser.add_argument("--frame", default="")
    parser.add_argument("--operation", default="rings")
    parser.add_argument("--warmup", type=int, default=proto.WARMUP)
    parser.add_argument("--repeats", type=int, default=proto.REPEATS)
    args = parser.parse_args()
    if args.worker == "ovito":
        worker_ovito(args.frame, args.operation, args.warmup, args.repeats)
    elif args.worker == "matscipy":
        worker_matscipy(args.frame, args.warmup, args.repeats)
    elif args.worker == "cycless":
        worker_cycless(args.frame, args.warmup, args.repeats)
    elif args.worker:
        raise SystemExit(args.worker)
    else:
        campaign(args)


if __name__ == "__main__":
    main()
