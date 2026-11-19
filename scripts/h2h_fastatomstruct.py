#!/usr/bin/env python3
"""fastatomstruct King rings and sharp-cutoff q6 on one ice Ic frame.

King's shortest-path rings are the R1 census of RINGS. q_l uses the
near-hard radial weight p1=48, which on this crystal is the four
neighbours inside 3.5 A. The frame is built with the cubic cell and
periodic boundaries on every axis. In 2.3.0 a frame with no cell is
non-periodic, and the Rust ring and neighbour calls return Result; the
Python entry points used here surface an error instead of a value.
Prints one JSON object.
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


def _matrix(value):
    magnitude = getattr(value, "magnitude", None)
    if magnitude is not None and not hasattr(value, "shape"):
        value = magnitude
    return np.asarray(value, dtype=np.float64)


def atoms_of(pos, side):
    """Periodic cubic cell. A missing cell is a different problem on 2.3.0."""
    elements = [fs.Element.O] * len(pos)
    cell = np.eye(3, dtype=np.float64) * float(side)
    atoms = fs.Atoms(np.ascontiguousarray(pos, dtype=np.float64), cell, elements)
    # 2.2.0 has no set_pbc. 2.3.0 treats a cell without it as non-periodic,
    # so the sitting checks run only on that API.
    if hasattr(atoms, "set_pbc"):
        atoms.set_pbc([True, True, True])
        pbc = [bool(flag) for flag in atoms.get_pbc()]
        if pbc != [True, True, True]:
            raise RuntimeError(f"frame pbc {pbc}")
        got = _matrix(atoms.get_cell())
        if got.shape != (3, 3) or not np.allclose(got, cell, rtol=0.0, atol=1e-8):
            raise RuntimeError("frame cell is not the periodic ice Ic cube")
    return atoms


def unwrap_result(value):
    """Return the success value if a binding hands back a Result object."""
    current = value
    for _ in range(3):
        is_ok = getattr(current, "is_ok", None)
        unwrap = getattr(current, "unwrap", None)
        if not callable(is_ok) or not callable(unwrap):
            return current
        if not is_ok():
            reveal = getattr(current, "unwrap_err", None)
            err = reveal() if callable(reveal) else current
            raise RuntimeError(f"fastatomstruct Result: {err}")
        current = unwrap()
    return current


def run(path, operation, warmup, repeats, once, threads=1):
    if hasattr(fs, "disable_telemetry"):
        fs.disable_telemetry()
    if hasattr(fs, "suppress_telemetry_notice"):
        fs.suppress_telemetry_notice()
    fs.set_num_threads(threads)
    if fs.get_num_threads() != threads:
        raise RuntimeError(f"fastatomstruct threads {fs.get_num_threads()}")
    pos, side = load_frame(path)
    atoms = atoms_of(pos, side)
    cutoff = fs.analysis.distance.CutoffMode.Fixed(proto.CUTOFF_A)

    def prep():
        return atoms

    if operation == "rings":
        method = fs.analysis.rings.RingMethod.King

        def call(frame):
            got = unwrap_result(fs.analysis.rings.ring_size_distribution(
                frame, cutoff, proto.MAX_RING, method))
            if not isinstance(got, tuple) or len(got) != 2:
                raise RuntimeError("ring result is not a pair")
            sizes, counts = got
            return {
                "sizes": [int(v) for v in sizes],
                "counts": [float(v) for v in counts],
            }
    elif operation == "q6":
        def call(frame):
            values = unwrap_result(fs.analysis.bond_order.q_l(frame, 6, cutoff, 48))
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
        "warmup": int(warmup),
        "periodic": True,
        "cell_side_A": side,
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
    parser.add_argument("--threads", type=int, default=1)
    parser.add_argument("--once", action="store_true")
    args = parser.parse_args()
    run(args.frame, args.operation, args.warmup, args.repeats, args.once, args.threads)


if __name__ == "__main__":
    main()
