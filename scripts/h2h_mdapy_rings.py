#!/usr/bin/env python3
# Copyright (c) d-SEAMS developers
"""mdapy 2.0 shortest-path rings on one generated ice Ic frame.

RingStatistics uses Franzblau's test: no two atoms of a ring are closer
through the network than along the ring. The timer covers compute, which
builds the cutoff bonds. The System is built outside the timer.
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
from mdapy.analysis import RingStatistics


def load_frame(path):
    data = np.load(path)
    return data["positions"], float(data["side"])


def census_of(result):
    counts = [int(result.counts[k]) for k in range(3, proto.MAX_RING + 1)]
    return {"sizes": list(range(3, proto.MAX_RING + 1)), "counts": counts}


def run(path, warmup, repeats, once, threads=1):
    mdapy.set_num_threads(threads)
    if mdapy.get_num_threads() != threads:
        raise RuntimeError(f"mdapy threads {mdapy.get_num_threads()}")
    pos, side = load_frame(path)

    def prep():
        return System(pos, box=side)

    def call(system):
        found = RingStatistics(
            max_size=proto.MAX_RING, min_size=3, cutoff=proto.CUTOFF_A
        ).compute(system)
        return census_of(found)

    if once:
        call(prep())
        return
    samples, result = proto.time_calls(prep, call, warmup, repeats)
    payload = {
        "code": "mdapy",
        "version": mdapy.__version__,
        "operation": "rings",
        "molecules": int(len(pos)),
        "samples_ms": samples,
        "census": result,
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
