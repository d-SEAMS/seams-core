#!/usr/bin/env python3
"""mdapy 2.0 ChillPlus on one generated ice Ic frame.

The timer covers ChillPlus.compute, which builds the cutoff neighbour
list. The System is built outside the timer. Prints one JSON object.
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
from mdapy.identify import ChillPlus


def load_frame(path):
    data = np.load(path)
    return data["positions"], float(data["side"])


def run(path, warmup, repeats, once):
    # mdapy 2.0 takes its pool from the cpus this process may use. Under
    # taskset -c 0 that is one thread. OMP_NUM_THREADS does not control it.
    mdapy.set_num_threads(1)
    if mdapy.get_num_threads() != 1:
        raise RuntimeError(f"mdapy threads {mdapy.get_num_threads()}")
    pos, side = load_frame(path)

    def prep():
        return System(pos, box=side)

    def call(system):
        return ChillPlus(cutoff=proto.CUTOFF_A).compute(system)

    if once:
        system = prep()
        call(system)
        return
    samples, result = proto.time_calls(prep, call, warmup, repeats)
    labels = np.asarray(result.structure, dtype=np.int32)
    mapped = [proto.MDAPY_CLASS[int(v)] for v in labels]
    payload = {
        "code": "mdapy",
        "version": mdapy.__version__,
        "operation": "chill+",
        "molecules": int(len(pos)),
        "samples_ms": samples,
        "labels": mapped,
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
    parser.add_argument("--once", action="store_true")
    args = parser.parse_args()
    run(args.frame, args.warmup, args.repeats, args.once)


if __name__ == "__main__":
    main()
