#!/usr/bin/env python3
"""fastatomstruct King rings and sharp-cutoff q6 on one ice Ic frame.

King's shortest-path rings are the R1 census of RINGS. q_l uses the
near-hard radial weight p1=48, which on this crystal is the four
neighbours inside 3.5 A. Prints one JSON object.
"""

from __future__ import annotations

import argparse
import json
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))

import h2h_protocol as proto

proto.apply_pin_env()

import numpy as np
import fastatomstruct as fs


def load_frame(path):
    data = np.load(path)
    return data["positions"], float(data["side"])


def atoms_of(pos, side):
    elements = [fs.Element.O] * len(pos)
    cell = np.eye(3, dtype=np.float64) * side
    return fs.Atoms(pos, cell, elements)


def run(path, operation, warmup, repeats, once):
    fs.disable_telemetry()
    fs.suppress_telemetry_notice()
    fs.set_num_threads(1)
    if fs.get_num_threads() != 1:
        raise RuntimeError(f"fastatomstruct threads {fs.get_num_threads()}")
    pos, side = load_frame(path)
    atoms = atoms_of(pos, side)
    cutoff = fs.analysis.distance.CutoffMode.Fixed(proto.CUTOFF_A)

    def prep():
        return atoms

    if operation == "rings":
        method = fs.analysis.rings.RingMethod.King

        def call(frame):
            sizes, counts = fs.analysis.rings.ring_size_distribution(
                frame, cutoff, proto.MAX_RING, method)
            return {
                "sizes": [int(v) for v in sizes],
                "counts": [int(round(float(v))) for v in counts],
            }
    elif operation == "q6":
        def call(frame):
            values = fs.analysis.bond_order.q_l(frame, 6, cutoff, 48)
            return np.asarray(values, dtype=np.float64)
    else:
        raise SystemExit(f"unknown operation {operation}")

    if once:
        call(prep())
        return
    samples, result = proto.time_calls(prep, call, warmup, repeats)
    payload = {
        "code": "fastatomstruct",
        "version": fs.__version__,
        "operation": operation,
        "molecules": int(len(pos)),
        "samples_ms": samples,
        "cpus": proto.cpus_allowed(),
        "threads": int(fs.get_num_threads()),
    }
    if operation == "rings":
        payload["census"] = result
    else:
        payload["q6"] = result.tolist()
    json.dump(payload, sys.stdout)
    sys.stdout.write("\n")


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--frame", required=True)
    parser.add_argument("--operation", choices=("rings", "q6"), required=True)
    parser.add_argument("--warmup", type=int, default=proto.WARMUP)
    parser.add_argument("--repeats", type=int, default=proto.REPEATS)
    parser.add_argument("--once", action="store_true")
    args = parser.parse_args()
    run(args.frame, args.operation, args.warmup, args.repeats, args.once)


if __name__ == "__main__":
    main()
