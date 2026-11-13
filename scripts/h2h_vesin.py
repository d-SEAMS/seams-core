#!/usr/bin/env python3
# Copyright (c) d-SEAMS developers
"""vesin cutoff neighbour list on one generated ice Ic frame.

The timer covers NeighborList.compute with a full (directed) list at
3.5 A. The calculator is built outside the timer. n_threads is the
vesin pool; zero would follow OMP_NUM_THREADS instead. Prints one JSON
object.
"""

from __future__ import annotations

import argparse
import importlib.metadata
import json
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))

import h2h_protocol as proto

proto.apply_pin_env()

import numpy as np
from vesin import NeighborList


def load_frame(path):
    data = np.load(path)
    return data["positions"], float(data["side"])


def run(path, warmup, repeats, once, threads=1):
    pos, side = load_frame(path)
    cell = np.eye(3, dtype=np.float64) * side
    calc = NeighborList(cutoff=proto.CUTOFF_A, full_list=True, n_threads=threads)

    def prep():
        return pos

    def call(points):
        left, right, _shifts = calc.compute(
            points=points, box=cell, periodic=True, quantities="ijS"
        )
        return int(len(left)), int(len(right))

    if once:
        call(prep())
        return
    samples, result = proto.time_calls(prep, call, warmup, repeats)
    pairs, other = result
    if pairs != other:
        raise RuntimeError(f"vesin pair ends differ: {pairs} {other}")
    payload = {
        "code": "vesin",
        "version": importlib.metadata.version("vesin"),
        "operation": "neighbours",
        "molecules": int(len(pos)),
        "samples_ms": samples,
        "pairs": pairs,
        "cpus": proto.cpus_allowed(),
        "threads": threads,
        "full_list": True,
        "cutoff_A": proto.CUTOFF_A,
    }
    json.dump(payload, sys.stdout)
    sys.stdout.write("\n")


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--frame", required=True)
    parser.add_argument("--warmup", type=int, default=proto.WARMUP)
    parser.add_argument("--repeats", type=int, default=proto.REPEATS)
    parser.add_argument("--threads", type=int, default=1)
    parser.add_argument("--once", action="store_true")
    args = parser.parse_args()
    run(args.frame, args.warmup, args.repeats, args.once, args.threads)


if __name__ == "__main__":
    main()
