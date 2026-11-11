#!/usr/bin/env python3
"""SQQ 0.5.8 shortest-path rings on ideal ice Ic.

Oxygen-only SOL residues, O-O cutoff 0.35 nm, ring sizes 4-6.
SQQ does not search three-rings. Cages and order parameters are off.
One warm-up, six calls, median and type-7 IQR. The ring walk is one
Python process; the eight-core row sets parallel.workers to 8.
"""

from __future__ import annotations

import argparse
import csv
import os
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))

import h2h_protocol as proto


def frame_of(pos, side):
    import numpy as np
    from sqq.models import Atom, Frame

    pos_nm = np.asarray(pos, dtype=float) / 10.0
    side_nm = float(side) / 10.0
    atoms = [
        Atom(index=i, resid=i + 1, resname="SOL", atomname="O", atomid=i + 1,
             xyz=pos_nm[i], molecule_id=i + 1)
        for i in range(len(pos_nm))
    ]
    return Frame(name=f"ice-{len(atoms)}", atoms=atoms,
                 box=np.array([side_nm, side_nm, side_nm]))


def config_for(workers):
    from sqq.api import load_config
    return load_config({
        "engine": "py",
        "graph": {"mode": "oo", "oo_cutoff_nm": 0.35},
        "ring": {"sizes": [4, 5, 6], "definition": "shortest_path", "chordless": True},
        "half_cage": {"enabled": False},
        "quasi_cage": {"enabled": False},
        "cage": {"enabled": False},
        "order_parameter": {"enabled": []},
        "parallel": {"workers": workers, "math_threads": 1},
    }, engine="py")


def census(result):
    return [len(result.rings.get(size, [])) for size in (4, 5, 6)]


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--out", default="benchmarks/h2h-watch-v6-sqq.csv")
    parser.add_argument("--repetitions", type=int, default=24)
    parser.add_argument("--threads", type=int, default=0)
    args = parser.parse_args()
    from sqq.api import analyze_frame
    import sqq

    pos, side = proto.ice_ic(args.repetitions)
    out = Path(__file__).resolve().parents[1] / args.out
    out.parent.mkdir(parents=True, exist_ok=True)
    version = getattr(sqq, "__version__", "0.5.8")
    fields = [
        "threads", "code", "version", "operation", "molecules",
        "warmup", "n", "median_ms", "q1_ms", "q3_ms", "iqr_ms",
        "agreement", "detail", "samples_ms", "cpus",
    ]
    with out.open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields)
        writer.writeheader()
        widths = (args.threads,) if args.threads else (1, 8)
        for threads in widths:
            cfg = config_for(threads)
            built = frame_of(pos, side)

            def call(_unused, cfg=cfg, built=built):
                return census(analyze_frame(built, cfg, engine="py"))

            samples, result = proto.time_calls(lambda: built, call)
            summary = proto.summarize(samples)
            expect = [0, 0, 2 * len(pos)]
            writer.writerow({
                "threads": threads,
                "code": "SQQ",
                "version": version,
                "operation": "rings",
                "molecules": len(pos),
                "warmup": proto.WARMUP,
                "n": summary["n"],
                "median_ms": f"{summary['median_ms']:.3f}",
                "q1_ms": f"{summary['q1_ms']:.3f}",
                "q3_ms": f"{summary['q3_ms']:.3f}",
                "iqr_ms": f"{summary['iqr_ms']:.3f}",
                "agreement": f"{float(result == expect):.6f}",
                "detail": f"sizes 4,5,6 = {result}; expect {expect}",
                "samples_ms": " ".join(f"{v:.3f}" for v in samples),
                "cpus": proto.cpus_allowed(),
            })
            handle.flush()
            print(f"threads={threads} {result} {samples}", flush=True)


if __name__ == "__main__":
    main()
