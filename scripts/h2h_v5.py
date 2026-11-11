#!/usr/bin/env python3
"""Symmetry head-to-head on ideal ice Ic, 110592 molecules.

The d-SEAMS column is classifyIce from the feat/symmetry kernel
(scale-free 17-atom shell, ceiling 0.18). Comparators are OVITO
polyhedral template matching, OVITO's diamond identification (extended
common-neighbour analysis), OVITO adaptive common-neighbour analysis,
and pyscal3's diamond, adaptive CNA and Ackland-Jones classifiers.
One warm-up, six calls, median and type-7 IQR. One core is taskset -c 0.
Eight cores are taskset -c 0-7. Each code uses its own thread count.
"""

from __future__ import annotations

import argparse
import csv
import json
import os
import subprocess
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
HERE = Path(__file__).resolve()
HOST_LIB = "/usr/lib/x86_64-linux-gnu"
# Tip of haozeke/feat/symmetry. 59d85155 is not an object on that branch.
KERNEL_SHA = "59d85155350cfa3c53da3df11d6ddf45b6b31ede"

FIELDS = [
    "threads", "code", "version", "operation", "molecules",
    "warmup", "n", "median_ms", "q1_ms", "q3_ms", "iqr_ms",
    "agreement", "detail", "samples_ms", "cpus",
]


def parse_args(argv):
    parser = argparse.ArgumentParser()
    parser.add_argument("--out", default="benchmarks/h2h-watch-v5.csv")
    parser.add_argument("--python")
    parser.add_argument("--kernel-bin", default=os.environ.get("SYMMETRY_BENCH", "/tmp/bench_symmetry_ice"))
    parser.add_argument("--kernel-sha", default=os.environ.get("SYMMETRY_SHA", KERNEL_SHA))
    parser.add_argument("--threads", type=int, default=1)
    parser.add_argument("--worker", choices=(
        "ovito-ptm", "ovito-ptm-all", "ovito-diamond", "ovito-cna", "pyscal"))
    parser.add_argument("--frame")
    parser.add_argument("--repetitions", type=int, default=24)
    return parser.parse_args(argv)


def child_env(kind, threads):
    env = os.environ.copy()
    width = str(threads)
    env["OMP_NUM_THREADS"] = width
    env["OMP_DYNAMIC"] = "FALSE"
    env["MKL_NUM_THREADS"] = "1"
    env["OPENBLAS_NUM_THREADS"] = "1"
    env["NUMEXPR_NUM_THREADS"] = "1"
    env.pop("LD_LIBRARY_PATH", None)
    if kind == "ovito":
        env["LD_LIBRARY_PATH"] = HOST_LIB
        env["OVITO_THREAD_COUNT"] = width
        env["QT_QPA_PLATFORM"] = "offscreen"
    else:
        env.pop("OVITO_THREAD_COUNT", None)
    return env


def emit(payload):
    json.dump(payload, sys.stdout)
    sys.stdout.write("\n")


def run_worker(args):
    pinned = open("/proc/self/status", encoding="utf-8").read()
    cpus = "unknown"
    for line in pinned.splitlines():
        if line.startswith("Cpus_allowed_list:"):
            cpus = line.split(":", 1)[1].strip()
    os.environ["OMP_NUM_THREADS"] = str(args.threads)
    sys.path.insert(0, str(HERE.parent))
    import numpy as np
    import h2h_protocol as proto

    data = np.load(args.frame)
    pos = data["positions"]
    side = float(data["side"])
    n = len(pos)

    if args.worker == "pyscal":
        import pyscal3
        from ase import Atoms
        pyscal3.set_num_threads(args.threads)
        atoms = Atoms(symbols=["O"] * n, positions=pos, cell=[side, side, side], pbc=True)
        version = getattr(pyscal3, "__version__", "")
        jobs = (
            ("diamond-cna", lambda: counts_of(pyscal3.diamond_structure(atoms))),
            ("adaptive-cna", lambda: counts_of(pyscal3.common_neighbor_analysis(atoms))),
            ("ackland-jones", lambda: counts_of(pyscal3.identify_ackland_jones(atoms))),
        )
        for operation, call in jobs:
            samples, result = proto.time_calls(lambda: atoms, lambda _a, call=call: call())
            emit({
                "code": "pyscal3", "version": version, "operation": operation,
                "molecules": int(n), "samples_ms": samples, "result": result,
                "cpus": cpus, "threads": args.threads,
            })
        return

    os.environ["OVITO_THREAD_COUNT"] = str(args.threads)
    os.environ["QT_QPA_PLATFORM"] = "offscreen"
    from ovito.data import DataCollection
    from ovito.modifiers import (
        CommonNeighborAnalysisModifier, IdentifyDiamondModifier,
        PolyhedralTemplateMatchingModifier)
    import ovito

    def prep():
        data = DataCollection()
        data.create_cell(
            matrix=[[side, 0, 0, 0], [0, side, 0, 0], [0, 0, side, 0]],
            pbc=(True, True, True))
        particles = data.create_particles(count=n)
        particles.create_property("Position", data=np.ascontiguousarray(pos))
        return data

    if args.worker in ("ovito-ptm", "ovito-ptm-all"):
        def call(data):
            mod = PolyhedralTemplateMatchingModifier()
            # The ice comparison enables the two diamond templates. The
            # all-template pass is a separate row.
            # Default OVITO leaves cubic and hexagonal diamond off.
            # The ice row turns those two on and the others off.
            keep = None if args.worker == "ovito-ptm-all" else {6, 7}
            for structure in mod.structures:
                structure.enabled = keep is None or structure.id in keep
            data.apply(mod)
            labels = np.asarray(data.particles["Structure Type"], dtype=np.int32)
            return histogram(labels)
        operation = "ptm" if args.worker == "ovito-ptm" else "ptm-all"
    elif args.worker == "ovito-diamond":
        def call(data):
            data.apply(IdentifyDiamondModifier())
            labels = np.asarray(data.particles["Structure Type"], dtype=np.int32)
            return histogram(labels)
        operation = "diamond-cna"
    else:
        def call(data):
            data.apply(CommonNeighborAnalysisModifier())
            labels = np.asarray(data.particles["Structure Type"], dtype=np.int32)
            return histogram(labels)
        operation = "adaptive-cna"
    samples, result = proto.time_calls(prep, call)
    emit({
        "code": "OVITO", "version": ovito.version_string, "operation": operation,
        "molecules": int(n), "samples_ms": samples, "histogram": result,
        "cpus": cpus, "threads": args.threads,
    })


def histogram(labels):
    import numpy as np
    values, counts = np.unique(labels, return_counts=True)
    return {str(int(v)): int(c) for v, c in zip(values, counts)}


def counts_of(value):
    import numpy as np
    if isinstance(value, dict):
        return {str(k): int(v) for k, v in value.items()}
    array = np.asarray(value)
    if array.ndim == 0:
        return {"value": str(value)}
    values, counts = np.unique(array, return_counts=True)
    return {str(v): int(c) for v, c in zip(values, counts)}


def kernel_run(binary, threads, repetitions):
    cores = "0" if threads == 1 else "0-7"
    cmd = ["taskset", "-c", cores, binary, str(repetitions), "1", "6"]
    done = subprocess.run(
        cmd, cwd=ROOT, env=child_env("kernel", threads), text=True, capture_output=True)
    if done.returncode not in (0, 2):
        tail = (done.stderr or done.stdout or "").strip().splitlines()
        return {"error": tail[-1] if tail else f"exit {done.returncode}", "code": "d-SEAMS"}
    text = done.stdout.strip().splitlines()[-1]
    fields = {}
    samples = []
    for part in text.split():
        if part.startswith("samples_ms="):
            continue
        if "=" in part:
            key, value = part.split("=", 1)
            fields[key] = value
        else:
            samples.append(float(part))
    # samples follow the samples_ms= token and are bare numbers already collected
    # when the token itself holds nothing. Re-parse the tail.
    tail = text.split("samples_ms=", 1)[1]
    samples = [float(v) for v in tail.split()]
    return {
        "code": "d-SEAMS", "version": fields.get("cpus", ""),
        "samples_ms": samples, "fields": fields, "cpus": fields.get("cpus", ""),
        "threads": threads, "molecules": int(fields["molecules"]),
    }


def pinned_python(python, threads, worker, frame):
    cores = "0" if threads == 1 else "0-7"
    kind = "ovito" if worker.startswith("ovito") else "pyscal"
    cmd = ["taskset", "-c", cores, python, str(HERE),
           "--worker", worker, "--threads", str(threads), "--frame", frame]
    done = subprocess.run(
        cmd, cwd=ROOT, env=child_env(kind, threads), text=True, capture_output=True)
    if done.returncode != 0:
        tail = (done.stderr or done.stdout or "").strip().splitlines()
        return {"error": tail[-1] if tail else f"exit {done.returncode}",
                "code": worker, "threads": threads}
    lines = [row for row in done.stdout.splitlines() if row.startswith("{")]
    if not lines:
        return [{"error": "no json", "code": worker, "threads": threads}]
    payloads = [json.loads(row) for row in lines]
    expect = cores
    for payload in payloads:
        if payload.get("cpus") != expect:
            return [{"error": f"not pinned: {payload.get('cpus')}", "code": worker,
                     "threads": threads}]
    return payloads


def write_row(writer, payload, operation, fraction, detail):
    sys.path.insert(0, str(HERE.parent))
    import h2h_protocol as proto
    if "error" in payload:
        writer.writerow({
            "threads": payload.get("threads", ""),
            "code": payload.get("code", ""),
            "version": "",
            "operation": operation,
            "molecules": payload.get("molecules", ""),
            "warmup": "",
            "n": 0,
            "median_ms": "", "q1_ms": "", "q3_ms": "", "iqr_ms": "",
            "agreement": "",
            "detail": payload["error"],
            "samples_ms": "",
            "cpus": "",
        })
        return
    summary = proto.summarize(payload["samples_ms"])
    writer.writerow({
        "threads": payload["threads"],
        "code": payload["code"],
        "version": payload.get("version", ""),
        "operation": operation,
        "molecules": payload["molecules"],
        "warmup": proto.WARMUP,
        "n": summary["n"],
        "median_ms": f"{summary['median_ms']:.3f}",
        "q1_ms": f"{summary['q1_ms']:.3f}",
        "q3_ms": f"{summary['q3_ms']:.3f}",
        "iqr_ms": f"{summary['iqr_ms']:.3f}",
        "agreement": fraction,
        "detail": detail,
        "samples_ms": " ".join(f"{v:.3f}" for v in payload["samples_ms"]),
        "cpus": payload.get("cpus", ""),
    })


def cubic_fraction(payload, molecules):
    """Share of atoms named cubic diamond / ice Ic."""
    if "error" in payload:
        return ""
    if payload["code"] == "d-SEAMS":
        ic = int(payload["fields"].get("Ic", "0"))
        return f"{ic / molecules:.6f}"
    hist = payload.get("histogram") or payload.get("result") or {}
    # PTM cubic diamond is 6. Identify diamond's cubic diamond is 1.
    # pyscal3 names the same lattice "cubic diamond".
    for key in ("6", "1", "cubic diamond"):
        if key in hist and hist[key]:
            return f"{int(hist[key]) / molecules:.6f}"
    return "0.000000"


def campaign(args):
    sys.path.insert(0, str(HERE.parent))
    os.environ.pop("LD_LIBRARY_PATH", None)
    import numpy as np
    import h2h_protocol as proto
    python = args.python or os.environ.get("H2H_PYTHON") or sys.executable
    out = ROOT / args.out
    out.parent.mkdir(parents=True, exist_ok=True)
    pos, side = proto.ice_ic(args.repetitions)
    frame = ROOT / "benchmarks" / f".ice-ic-{len(pos)}.npz"
    np.savez(frame, positions=pos, side=side)
    workers = ("ovito-ptm", "ovito-ptm-all", "ovito-diamond", "ovito-cna", "pyscal")
    with out.open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=FIELDS)
        writer.writeheader()
        for threads in (1, 8):
            kernel = kernel_run(args.kernel_bin, threads, args.repetitions)
            kernel["threads"] = threads
            kernel["version"] = args.kernel_sha
            kernel["molecules"] = len(pos)
            detail = ""
            if "fields" in kernel:
                detail = (
                    f"Ih={kernel['fields'].get('Ih')} Ic={kernel['fields'].get('Ic')} "
                    f"Isd={kernel['fields'].get('Isd')} other={kernel['fields'].get('other')}"
                )
            write_row(writer, kernel, "classify-ice",
                      cubic_fraction(kernel, len(pos)) if "error" not in kernel else "",
                      detail)
            handle.flush()
            print(f"threads={threads} kernel {kernel.get('samples_ms', kernel.get('error'))}",
                  flush=True)
            for worker in workers:
                rivals = pinned_python(python, threads, worker, str(frame))
                for rival in rivals:
                    rival["molecules"] = len(pos)
                    rival["threads"] = threads
                    operation = rival.get("operation", worker)
                    if "result" in rival:
                        detail = json.dumps(rival["result"])
                        fraction = ""
                    else:
                        detail = json.dumps(rival.get("histogram", rival.get("error")))
                        fraction = cubic_fraction(rival, len(pos)) if "error" not in rival else ""
                    write_row(writer, rival, operation, fraction, detail)
                    handle.flush()
                    print(
                        f"threads={threads} {operation} "
                        f"{rival.get('samples_ms', rival.get('error'))}",
                        flush=True)
    frame.unlink(missing_ok=True)


def main(argv=None):
    args = parse_args(sys.argv[1:] if argv is None else argv)
    if args.worker:
        run_worker(args)
        return
    campaign(args)


if __name__ == "__main__":
    main()
