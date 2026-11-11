#!/usr/bin/env python3
"""OVITO 3.17.2 Chill+ on one generated ice Ic frame.

OVITO_THREAD_COUNT is set before ovito is imported. The timer covers
pipeline.compute, which builds the cutoff bonds and then classifies.
The particle data is built outside the timer. Prints one JSON object.
"""

from __future__ import annotations

import argparse
import json
import os
import sys
from pathlib import Path

os.environ["OVITO_THREAD_COUNT"] = "1"
os.environ.setdefault("QT_QPA_PLATFORM", "offscreen")

sys.path.insert(0, str(Path(__file__).resolve().parent))

import h2h_protocol as proto

proto.apply_pin_env()
os.environ["OVITO_THREAD_COUNT"] = "1"

import numpy as np
import ovito
from ovito.data import DataCollection
from ovito.modifiers import ChillPlusModifier
from ovito.pipeline import Pipeline, StaticSource


def load_frame(path):
    data = np.load(path)
    return data["positions"], float(data["side"])


def run(path, warmup, repeats, once):
    if tuple(ovito.version)[:3] != (3, 17, 2):
        raise RuntimeError(f"need OVITO 3.17.2, got {ovito.version}")
    pos, side = load_frame(path)
    names = {
        int(ChillPlusModifier.Type.OTHER): "other",
        int(ChillPlusModifier.Type.HEXAGONAL_ICE): "hexagonal",
        int(ChillPlusModifier.Type.CUBIC_ICE): "cubic",
        int(ChillPlusModifier.Type.INTERFACIAL_ICE): "interfacial",
        int(ChillPlusModifier.Type.HYDRATE): "hydrate",
        int(ChillPlusModifier.Type.INTERFACIAL_HYDRATE): "interfacial_hydrate",
    }

    def prep():
        data = DataCollection()
        data.create_cell(np.diag([side, side, side]), pbc=(True, True, True))
        parts = data.create_particles(count=len(pos))
        parts.create_property("Position", data=pos)
        return data

    def call(data):
        pipe = Pipeline(source=StaticSource(data=data))
        pipe.modifiers.append(ChillPlusModifier(cutoff=proto.CUTOFF_A))
        out = pipe.compute()
        return np.asarray(out.particles["Structure Type"])

    if once:
        call(prep())
        return
    samples, result = proto.time_calls(prep, call, warmup, repeats)
    mapped = [names[int(v)] for v in result]
    payload = {
        "code": "OVITO",
        "version": ".".join(str(v) for v in tuple(ovito.version)[:3]),
        "operation": "chill+",
        "molecules": int(len(pos)),
        "samples_ms": samples,
        "labels": mapped,
        "cpus": proto.cpus_allowed(),
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
