#!/usr/bin/env python3
# Copyright (c) d-SEAMS developers
"""mdapy 2.0 cutoff q6 on one generated ice Ic frame.

SteinhardtParameters with a 3.5 A cutoff builds the neighbour list inside
compute. The System is built outside the timer. The local q6 is not the
Lechner-Dellago average.
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
import mdapy
from mdapy import System
from mdapy.analysis import SteinhardtParameters


def load_frame(path):
    data = np.load(path)
    return data["positions"], float(data["side"])


def run(path, warmup, repeats, once, threads=1):
    mdapy.set_num_threads(threads)
    if mdapy.get_num_threads() != threads:
        raise RuntimeError(f"mdapy threads {mdapy.get_num_threads()}")
    pos, side = load_frame(path)
    order = SteinhardtParameters(degrees=(6,), cutoff=proto.CUTOFF_A, average=False)

    def prep():
        return System(pos, box=side)

    def call(system):
        result = order.compute(system)
        return np.asarray(result.q, dtype=np.float64).reshape(-1)

    if once:
        call(prep())
        return
    samples, result = proto.time_calls(prep, call, warmup, repeats)
    payload = {
        "code": "mdapy",
        "version": mdapy.__version__,
        "operation": "q6",
        "molecules": int(len(pos)),
        "samples_ms": samples,
        "q6": result.tolist(),
        "cpus": proto.cpus_allowed(),
        "threads": int(mdapy.get_num_threads()),
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
