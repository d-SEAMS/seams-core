#!/usr/bin/env python3
"""Head-to-head on the brine frames.

Serial, one core, one thread. One warm-up, then at least six timed calls.
The median of an even sample is the mean of the two central values.
Quartiles are Hyndman-Fan type 7. The timed region is the classifier or
the neighbour list on the water oxygens (LAMMPS type 1). Ions are not
in that set. Frames come from input/brine/regenerate.py and are not
committed.

    taskset -c 0 python3 scripts/h2h_brine.py
"""

from __future__ import annotations

import ctypes
import math
import os
import sys
from pathlib import Path

os.environ["OMP_NUM_THREADS"] = "1"
os.environ["OMP_DYNAMIC"] = "FALSE"
os.environ["MKL_NUM_THREADS"] = "1"
os.environ["OPENBLAS_NUM_THREADS"] = "1"
os.environ["NUMEXPR_NUM_THREADS"] = "1"
os.environ["RAYON_NUM_THREADS"] = "1"

try:
    os.sched_setaffinity(0, {0})
except (AttributeError, OSError):
    pass

import numpy as np

ROOT = Path(__file__).resolve().parents[1]
CUTOFF = 3.5
WARMUP = 1
REPEATS = 6

# Coarse labels shared with mdapy ChillPlus. d-SEAMS reCubic/reHex fold
# into cubic/hexagonal. Liquid and unclassified fold into other.
DSEAMS_COARSE = {
    0: "cubic",
    1: "hexagonal",
    2: "other",
    3: "interfacial",
    4: "hydrate",
    5: "interfacial_hydrate",
    6: "other",
    7: "cubic",
    8: "hexagonal",
}
MDAPY_COARSE = {
    0: "other",
    1: "hexagonal",
    2: "cubic",
    3: "interfacial",
    4: "hydrate",
    5: "interfacial_hydrate",
}


def quantile_linear(sorted_vals, q):
    if len(sorted_vals) == 1:
        return sorted_vals[0]
    pos = q * (len(sorted_vals) - 1)
    lo = int(pos)
    frac = pos - lo
    if lo + 1 >= len(sorted_vals):
        return sorted_vals[-1]
    return sorted_vals[lo] * (1.0 - frac) + sorted_vals[lo + 1] * frac


def summarize(samples):
    if len(samples) < REPEATS:
        raise ValueError(f"need n >= {REPEATS}, got {len(samples)}")
    ordered = sorted(float(v) for v in samples)
    mid = len(ordered) // 2
    if len(ordered) % 2 == 0:
        median = 0.5 * (ordered[mid - 1] + ordered[mid])
    else:
        median = ordered[mid]
    q1 = quantile_linear(ordered, 0.25)
    q3 = quantile_linear(ordered, 0.75)
    return median, q3 - q1


def time_calls(prep, run):
    samples = []
    last = None
    for i in range(WARMUP + REPEATS):
        ready = prep()
        t0 = time_ns()
        last = run(ready)
        dt = (time_ns() - t0) / 1e6
        if i >= WARMUP:
            samples.append(dt)
    return samples, last


def time_ns():
    return __import__("time").perf_counter_ns()


def read_dump(path, oxygen_type=1):
    oxy = []
    box = None
    with open(path) as handle:
        mode = None
        bounds = []
        for line in handle:
            text = line.strip()
            if text.startswith("ITEM: BOX"):
                mode = "box"
                bounds = []
                continue
            if text.startswith("ITEM: ATOMS"):
                mode = "atoms"
                continue
            if text.startswith("ITEM:"):
                mode = None
                continue
            if mode == "box":
                parts = text.split()
                bounds.append(float(parts[1]) - float(parts[0]))
                if len(bounds) == 3:
                    box = bounds
            elif mode == "atoms":
                parts = text.split()
                if int(float(parts[2])) == oxygen_type:
                    oxy.append((float(parts[3]), float(parts[4]), float(parts[5])))
            if box is not None and mode is None and oxy and text.startswith("ITEM: TIMESTEP") and len(oxy) > 0:
                break
    if not oxy or box is None:
        raise RuntimeError(f"no oxygens in {path}")
    points = np.asarray(oxy, dtype=np.float64)
    return points, box


def load_dseams():
    lib = ctypes.CDLL(str(ROOT / "bbdir" / "src" / "libyodaLib.so"))
    lib.seams_chill_plus_species.argtypes = [
        ctypes.POINTER(ctypes.c_double),
        ctypes.c_int,
        ctypes.POINTER(ctypes.c_double),
        ctypes.POINTER(ctypes.c_int),
        ctypes.POINTER(ctypes.c_int),
        ctypes.c_int,
        ctypes.c_int,
        ctypes.POINTER(ctypes.c_int),
        ctypes.c_int,
        ctypes.c_double,
        ctypes.c_char_p,
        ctypes.POINTER(ctypes.c_int),
        ctypes.POINTER(ctypes.c_int),
        ctypes.POINTER(ctypes.c_int),
        ctypes.POINTER(ctypes.c_double),
    ]
    lib.seams_chill_plus_species.restype = ctypes.c_int
    lib.seams_oxygen_cutoff_degree.argtypes = [
        ctypes.POINTER(ctypes.c_double),
        ctypes.c_int,
        ctypes.POINTER(ctypes.c_double),
        ctypes.POINTER(ctypes.c_int),
        ctypes.POINTER(ctypes.c_int),
        ctypes.c_int,
        ctypes.c_int,
        ctypes.POINTER(ctypes.c_int),
        ctypes.c_int,
        ctypes.c_double,
        ctypes.POINTER(ctypes.c_int),
    ]
    lib.seams_oxygen_cutoff_degree.restype = ctypes.c_int
    lib.seams_c_error.restype = ctypes.c_char_p
    return lib


def frames():
    work = ROOT / "input" / "brine" / "work"
    found = sorted(work.glob("*/dump.lammpstrj"))
    # One bulk concentration and the two ice geometries, at every size
    # that finished. 1 mol/kg is the middle of the set.
    keep = []
    for path in found:
        name = path.parent.name
        if name.startswith("nacl-") and "-m1.0" in name:
            keep.append(path)
        elif name.startswith("slab-") or name.startswith("front-"):
            keep.append(path)
    return keep


def coarse_dseams(labels):
    return np.array([DSEAMS_COARSE[int(v)] for v in labels], dtype=object)


def agreement(a, b):
    if len(a) != len(b) or len(a) == 0:
        return float("nan")
    return float(np.mean(a == b))


def main():
    paths = frames()
    if not paths:
        print("no brine dumps under input/brine/work; run input/brine/regenerate.py")
        return 1
    lib = load_dseams()
    mdapy = None
    freud = None
    try:
        import mdapy as mdapy_mod
        from mdapy.identify import ChillPlus
        from mdapy.neighbor import NeighborList

        mdapy = (mdapy_mod, ChillPlus, NeighborList)
    except ImportError as exc:
        print("mdapy unavailable:", exc)
    try:
        import freud as freud_mod

        freud = freud_mod
    except ImportError as exc:
        print("freud unavailable:", exc)
    os.environ.setdefault("QT_QPA_PLATFORM", "offscreen")
    libgl = "/usr/lib/x86_64-linux-gnu"
    os.environ["LD_LIBRARY_PATH"] = libgl + ":" + os.environ.get("LD_LIBRARY_PATH", "")
    ovito_note = "not imported"
    try:
        from ovito.data import DataCollection, Particles, SimulationCell
        from ovito.modifiers import ChillPlusModifier
        from ovito.pipeline import Pipeline, StaticSource

        ovito_ok = True
    except ImportError as exc:
        ovito_ok = False
        ovito_note = str(exc)

    rows = []
    for path in paths:
        points, box = read_dump(path)
        n = points.shape[0]
        xyz = np.ascontiguousarray(points.reshape(-1))
        box9 = np.array([0, 0, 0, box[0], box[1], box[2], 0, 0, 0], dtype=np.float64)
        types = np.ones(n, dtype=np.int32)
        labels = np.empty(n, dtype=np.int32)
        degree = np.empty(n, dtype=np.int32)

        def prep_ds():
            return None

        def run_chill(_):
            rc = lib.seams_chill_plus_species(
                xyz.ctypes.data_as(ctypes.POINTER(ctypes.c_double)),
                n,
                box9.ctypes.data_as(ctypes.POINTER(ctypes.c_double)),
                types.ctypes.data_as(ctypes.POINTER(ctypes.c_int)),
                None,
                1,
                0,
                None,
                0,
                CUTOFF,
                b"cutoff",
                labels.ctypes.data_as(ctypes.POINTER(ctypes.c_int)),
                None,
                None,
                None,
            )
            if rc != 0:
                raise RuntimeError(lib.seams_c_error())
            return labels.copy()

        def run_neigh(_):
            rc = lib.seams_oxygen_cutoff_degree(
                xyz.ctypes.data_as(ctypes.POINTER(ctypes.c_double)),
                n,
                box9.ctypes.data_as(ctypes.POINTER(ctypes.c_double)),
                types.ctypes.data_as(ctypes.POINTER(ctypes.c_int)),
                None,
                1,
                0,
                None,
                0,
                CUTOFF,
                degree.ctypes.data_as(ctypes.POINTER(ctypes.c_int)),
            )
            if rc != 0:
                raise RuntimeError(lib.seams_c_error())
            return int(degree.sum())

        chill_samples, dseams_labels = time_calls(prep_ds, run_chill)
        neigh_samples, dseams_degree = time_calls(prep_ds, run_neigh)
        dmed, diqr = summarize(chill_samples)
        nmed, niqr = summarize(neigh_samples)
        rows.append((path.parent.name, "d-SEAMS", "chill+", n, dmed, diqr, 1.0))
        rows.append((path.parent.name, "d-SEAMS", "neighbours", n, nmed, niqr, 1.0))
        ref = coarse_dseams(dseams_labels)

        if mdapy is not None:
            mdapy_mod, ChillPlus, NeighborList = mdapy
            system = mdapy_mod.System(points, box)

            def run_mdapy(_):
                return ChillPlus(CUTOFF).compute(system).structure.copy()

            def run_mdapy_n(_):
                return NeighborList.build(system, CUTOFF).j.size

            samples, structure = time_calls(prep_ds, run_mdapy)
            med, iqr = summarize(samples)
            other = np.array([MDAPY_COARSE[int(v)] for v in structure], dtype=object)
            rows.append((path.parent.name, "mdapy", "chill+", n, med, iqr, agreement(ref, other)))
            samples, _ = time_calls(prep_ds, run_mdapy_n)
            med, iqr = summarize(samples)
            rows.append((path.parent.name, "mdapy", "neighbours", n, med, iqr, float("nan")))

        if freud is not None:
            fbox = freud.box.Box(box[0], box[1], box[2])

            def run_freud(_):
                query = freud.locality.AABBQuery(fbox, points)
                return query.query(points, {"r_max": CUTOFF, "exclude_ii": True}).toNeighborList().neighbor_counts.sum()

            samples, _ = time_calls(prep_ds, run_freud)
            med, iqr = summarize(samples)
            rows.append((path.parent.name, "freud", "neighbours", n, med, iqr, float("nan")))

        if ovito_ok:
            def run_ovito(_):
                data = DataCollection()
                particles = Particles()
                particles.create_property("Position", data=points)
                data.objects.append(particles)
                cell = SimulationCell(pbc=(True, True, True))
                cell[:, 0] = (box[0], 0.0, 0.0)
                cell[:, 1] = (0.0, box[1], 0.0)
                cell[:, 2] = (0.0, 0.0, box[2])
                data.objects.append(cell)
                pipe = Pipeline(source=StaticSource(data=data))
                pipe.modifiers.append(ChillPlusModifier(cutoff=CUTOFF))
                return np.asarray(pipe.compute().particles["Structure Type"])

            try:
                samples, structure = time_calls(prep_ds, run_ovito)
                med, iqr = summarize(samples)
                # OVITO 3.17.2 Chill+ uses the same integers as mdapy:
                # 0 other, 1 hexagonal, 2 cubic, 3 interfacial, 4 hydrate,
                # 5 interfacial hydrate. That modifier is the open-edge
                # rule, not Table 1's closed edges and neighbour test.
                other = np.array(
                    [MDAPY_COARSE.get(int(v), "other") for v in structure],
                    dtype=object,
                )
                rows.append(
                    (path.parent.name, "OVITO", "chill+", n, med, iqr, agreement(ref, other))
                )
            except Exception as exc:
                ovito_note = str(exc)

        print(
            f"{path.parent.name} n={n} dseams_chill {dmed:.3f}/{diqr:.3f} "
            f"neigh {nmed:.3f}/{niqr:.3f} degree_sum {dseams_degree}",
            flush=True,
        )

    print("system,code,operation,oxygens,n,median_ms,iqr_ms,agreement")
    for row in rows:
        name, code, op, n, med, iqr, agree = row
        agree_s = "" if math.isnan(agree) else f"{agree:.4f}"
        print(f"{name},{code},{op},{n},{REPEATS},{med:.4f},{iqr:.4f},{agree_s}")
    if not ovito_ok:
        print("ovito: not imported")
    elif ovito_note != "not imported":
        print("ovito:", ovito_note)
    return 0


if __name__ == "__main__":
    sys.exit(main())
