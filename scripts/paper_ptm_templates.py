#!/usr/bin/env python3
"""PTM template sets, and per-atom agreement with the symmetry labels.

The classifier writes ice Ih, ice Ic, stacking-disordered ice, or other
for every atom. Polyhedral template matching with only the two diamond
templates names cubic diamond, hexagonal diamond, or Other. Its default
set names FCC, HCP, BCC, or Other and leaves both diamond templates off.
Those are different searches. This script times each search at one and
eight threads, then checks the labels on ideal ice Ih, ideal ice Ic, and
the 230 K brine start frame.

The brine frame is the starting coordinates of input/brine/grow.py at
230 K and 1 mol/kg (seed 1330, three by three by two ice cells), before
minimisation or dynamics. Oxygens only. Ions are not labelled.
CHILL+ is compared where a class exists on both sides: ice Ih with
hexagonal ice, ice Ic with cubic ice, and other with Other.
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
sys.path.insert(0, str(ROOT / "scripts"))
import paper_stats as stats

MOLECULES = 110592
OURS = {-1: "other", 0: "Ih", 1: "Ic", 2: "Isd"}
# Classes both methods can name. Stacking-disordered ice has no PTM
# template and no CHILL+ class.
PTM_MATCH = {"Ih": 7, "Ic": 6, "other": 0}
CHILL_MATCH = {"Ih": 1, "Ic": 2, "other": 0}
PTM_NAME = {
    0: "Other",
    1: "FCC",
    2: "HCP",
    3: "BCC",
    4: "ICO",
    5: "SC",
    6: "cubic-diamond",
    7: "hex-diamond",
    8: "graphene",
}
CHILL_NAME = {
    0: "Other",
    1: "hexagonal-ice",
    2: "cubic-ice",
    3: "interfacial-ice",
    4: "hydrate",
    5: "interfacial-hydrate",
}


def lib_path() -> str:
    build = ROOT / "bbdir-release"
    parts = [
        build / "src",
        build / "subprojects/readcon-core",
        build / "subprojects/readcon-db",
        build / "subprojects/linkcell",
        build / "subprojects/vesin",
        build / "subprojects/minimage",
        Path("/usr/lib/x86_64-linux-gnu"),
    ]
    return os.pathsep.join(str(path) for path in parts)


def ice_ic(cells: int):
    import h2h_protocol as proto

    pos, side = proto.ice_ic(cells)
    lattice = [
        [side, 0.0, 0.0],
        [0.0, side, 0.0],
        [0.0, 0.0, side],
    ]
    return pos, lattice


def ice_ih():
    """Orthogonal ice Ih from the brine cell, replicated until the box
    is large enough for the orthorhombic shell walk."""
    sys.path.insert(0, str(ROOT / "input" / "brine"))
    import regenerate as regen
    import numpy as np

    nx = ny = nz = 3
    pts = []
    for ix in range(nx):
        for iy in range(ny):
            for iz in range(nz):
                shift = (
                    ix * regen.ICE_BOX[0],
                    iy * regen.ICE_BOX[1],
                    iz * regen.ICE_BOX[2],
                )
                for oxygen in regen.ICE_OXYGEN:
                    pts.append(
                        (
                            oxygen[0] + shift[0],
                            oxygen[1] + shift[1],
                            oxygen[2] + shift[2],
                        )
                    )
    box = (
        nx * regen.ICE_BOX[0],
        ny * regen.ICE_BOX[1],
        nz * regen.ICE_BOX[2],
    )
    lattice = [
        [box[0], 0.0, 0.0],
        [0.0, box[1], 0.0],
        [0.0, 0.0, box[2]],
    ]
    return np.ascontiguousarray(pts, dtype=float), lattice


def brine_start():
    """Starting oxygens of the 230 K, 1 mol/kg growth slab."""
    sys.path.insert(0, str(ROOT / "input" / "brine"))
    import grow
    import numpy as np

    temperature = 230
    molality = 1.0
    seed = 1000 + temperature + int(molality * 100)
    waters, _na, _cl, box, n_ice = grow.build_growth(3, 3, 2, molality, seed)
    oxygens = np.ascontiguousarray([trip[0] for trip in waters], dtype=float)
    lattice = [
        [box[0], 0.0, 0.0],
        [0.0, box[1], 0.0],
        [0.0, 0.0, box[2]],
    ]
    return oxygens, lattice, n_ice, seed


def classify(lib, pos, lattice):
    import ctypes
    import numpy as np

    xyz = np.ascontiguousarray(pos, dtype=np.float64)
    lat = np.ascontiguousarray(lattice, dtype=np.float64).reshape(9)
    n = len(xyz)
    labels = np.empty(n, dtype=np.int32)
    rc = lib.seams_symmetry_classify(
        xyz.ctypes.data_as(ctypes.POINTER(ctypes.c_double)),
        n,
        lat.ctypes.data_as(ctypes.POINTER(ctypes.c_double)),
        labels.ctypes.data_as(ctypes.POINTER(ctypes.c_int)),
    )
    if rc != 0:
        raise RuntimeError(f"seams_symmetry_classify returned {rc}")
    return labels


def open_core(path: str):
    import ctypes

    lib = ctypes.CDLL(path)
    lib.seams_symmetry_classify.restype = ctypes.c_int
    lib.seams_symmetry_classify.argtypes = [
        ctypes.POINTER(ctypes.c_double),
        ctypes.c_int,
        ctypes.POINTER(ctypes.c_double),
        ctypes.POINTER(ctypes.c_int),
    ]
    return lib


def ovito_labels(pos, lattice, kind: str):
    import numpy as np
    from ovito.data import DataCollection, Particles, SimulationCell
    from ovito.modifiers import ChillPlusModifier, PolyhedralTemplateMatchingModifier
    from ovito.pipeline import Pipeline, StaticSource

    data = DataCollection()
    cell = SimulationCell(pbc=(True, True, True))
    matrix = np.zeros((3, 4))
    lat = np.asarray(lattice, dtype=float)
    matrix[:, 0] = lat[0]
    matrix[:, 1] = lat[1]
    matrix[:, 2] = lat[2]
    cell[:, :] = matrix
    parts = Particles()
    parts.create_property("Position", data=np.ascontiguousarray(pos, dtype=float))
    data.objects.append(cell)
    data.objects.append(parts)
    pipe = Pipeline(source=StaticSource(data=data))
    if kind == "chill":
        pipe.modifiers.append(ChillPlusModifier(cutoff=3.5))
    else:
        mod = PolyhedralTemplateMatchingModifier(rmsd_cutoff=0.1)
        if kind == "ice":
            for structure in mod.structures:
                structure.enabled = structure.id in (6, 7)
        elif kind != "default":
            raise RuntimeError(kind)
        pipe.modifiers.append(mod)
    out = pipe.compute()
    return np.asarray(out.particles["Structure Type"], dtype=np.int32)


def time_ptm(pos, lattice, kind: str, warmup: int, repeats: int):
    import time

    samples = []
    for index in range(warmup + repeats):
        t0 = time.perf_counter()
        labels = ovito_labels(pos, lattice, kind)
        elapsed = (time.perf_counter() - t0) * 1000.0
        if labels.shape[0] != len(pos):
            raise RuntimeError("PTM dropped atoms")
        if index >= warmup:
            samples.append(elapsed)
    return samples


def crosstab(ours, theirs, names):
    import numpy as np

    rows = []
    for our_code, our_name in OURS.items():
        mask = ours == our_code
        total = int(np.count_nonzero(mask))
        if total == 0:
            continue
        seen = theirs[mask]
        values, counts = np.unique(seen, return_counts=True)
        for value, count in zip(values, counts):
            their = names.get(int(value), str(int(value)))
            rows.append((our_name, their, int(count), total))
    return rows


def frame_of(name: str):
    if name == "ic":
        pos, lattice = ice_ic(24)
        note = "ideal ice Ic, 24^3 conventional cells"
        return pos, lattice, note
    if name == "ih":
        pos, lattice = ice_ih()
        note = "orthogonal ice Ih, 3x3x3 brine cells"
        return pos, lattice, note
    if name == "brine":
        pos, lattice, n_ice, seed = brine_start()
        note = (
            f"230 K brine start, 1 mol/kg, seed {seed}, "
            f"ice oxygens {n_ice}, liquid oxygens {len(pos) - n_ice}"
        )
        return pos, lattice, note
    raise RuntimeError(name)


def emit_samples(samples):
    print(" ".join(f"{value:.6f}" for value in samples))


def worker(args) -> int:
    os.environ["OVITO_THREAD_COUNT"] = str(args.threads)
    os.environ["QT_QPA_PLATFORM"] = "offscreen"
    os.environ["OMP_NUM_THREADS"] = str(args.threads)
    os.environ["OMP_DYNAMIC"] = "FALSE"
    if args.worker == "time":
        pos, lattice = ice_ic(24)
        if len(pos) != MOLECULES:
            raise SystemExit(f"ice Ic has {len(pos)} molecules")
        emit_samples(time_ptm(pos, lattice, args.template, args.warmup, args.repeats))
        return 0
    pos, lattice, note = frame_of(args.frame)
    lib = open_core(args.lib)
    ours = classify(lib, pos, lattice)
    ptm = ovito_labels(pos, lattice, "ice")
    chill = ovito_labels(pos, lattice, "chill")
    payload = {
        "frame": args.frame,
        "note": note,
        "molecules": int(len(pos)),
        "ptm": crosstab(ours, ptm, PTM_NAME),
        "chill": crosstab(ours, chill, CHILL_NAME),
    }
    print(json.dumps(payload))
    return 0


def launch(args, worker_name, template, threads, frame):
    cmd = [
        "taskset",
        "-c",
        stats.taskset_cores(threads),
        sys.executable,
        str(Path(__file__).resolve()),
        "--worker",
        worker_name,
        "--template",
        template,
        "--threads",
        str(threads),
        "--warmup",
        str(args.warmup),
        "--repeats",
        str(args.repeats),
        "--frame",
        frame,
        "--lib",
        args.lib,
    ]
    env = stats.pin_env(threads)
    env["LD_LIBRARY_PATH"] = lib_path()
    env["OVITO_THREAD_COUNT"] = str(threads)
    env["QT_QPA_PLATFORM"] = "offscreen"
    done = subprocess.run(cmd, cwd=ROOT, env=env, text=True, capture_output=True)
    if done.returncode != 0:
        tail = (done.stderr or done.stdout or "").strip().splitlines()
        raise RuntimeError(tail[-1] if tail else f"exit {done.returncode}")
    return [line for line in done.stdout.splitlines() if line.strip()][-1]


def write_time(path: Path, rows: list[dict]) -> None:
    with path.open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(rows[0]))
        writer.writeheader()
        writer.writerows(rows)


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--out", type=Path, default=ROOT / "benchmarks/paper-cells")
    parser.add_argument("--lib", default=str(ROOT / "bbdir-release/src/libyodaLib.so"))
    parser.add_argument("--worker", choices=("time", "agree"))
    parser.add_argument("--template", choices=("default", "ice"), default="ice")
    parser.add_argument("--frame", choices=("ic", "ih", "brine"), default="ic")
    parser.add_argument("--threads", type=int, default=8)
    parser.add_argument("--warmup", type=int, default=1)
    parser.add_argument("--repeats", type=int, default=8)
    args = parser.parse_args()
    if args.worker:
        return worker(args)
    if args.repeats < 8:
        raise SystemExit("n must be at least 8")
    if not Path(args.lib).is_file():
        raise SystemExit(f"missing {args.lib}")
    args.out.mkdir(parents=True, exist_ok=True)
    cores, affinity = stats.cores_now()
    log = [
        f"cores {cores}",
        f"affinity {affinity}",
        f"warmup {args.warmup}",
        f"n {args.repeats}",
        "rmsd_cutoff 0.1",
        "timing frame ideal ice Ic 110592",
        "default templates FCC HCP BCC",
        "ice templates cubic-diamond hex-diamond",
        "classifier labels every atom; ice-only PTM names only the two diamond structures",
    ]
    for template, filename in (("default", "ptm-default.csv"), ("ice", "ptm-ice.csv")):
        rows = []
        for threads in (1, 8):
            text = launch(args, "time", template, threads, "ic")
            samples = [float(piece) for piece in text.split()]
            if len(samples) != args.repeats:
                raise SystemExit(f"{template} threads {threads} returned {len(samples)} samples")
            summary = stats.summarize(samples)
            rows.append({
                "template_set": template,
                "molecules": MOLECULES,
                "threads": threads,
                "cores": threads,
                "affinity": stats.taskset_cores(threads),
                "warmup": args.warmup,
                "n": summary["n"],
                "median_ms": f"{summary['median_ms']:.6f}",
                "iqr_ms": f"{summary['iqr_ms']:.6f}",
                "q1_ms": f"{summary['q1_ms']:.6f}",
                "q3_ms": f"{summary['q3_ms']:.6f}",
                "rmsd_cutoff": 0.1,
                "samples_ms": stats.fmt_samples(samples),
            })
            log.append(
                f"time {template} threads {threads} median_ms {summary['median_ms']:.6f} "
                f"iqr_ms {summary['iqr_ms']:.6f} samples {stats.fmt_samples(samples)}"
            )
            print(
                f"{template} threads {threads} {summary['median_ms']:.4f} ms",
                flush=True,
            )
        write_time(args.out / filename, rows)
    agree_rows = []
    for frame in ("ih", "ic", "brine"):
        payload = json.loads(launch(args, "agree", "ice", 8, frame))
        log.append(f"frame {frame} {payload['note']} molecules {payload['molecules']}")
        for code, match in (("PTM", PTM_MATCH), ("CHILL+", CHILL_MATCH)):
            key = "ptm" if code == "PTM" else "chill"
            for our_name, their_name, count, total in payload[key]:
                expected = match.get(our_name)
                applies = our_name in match
                agreed = applies and their_name == {
                    7: "hex-diamond",
                    6: "cubic-diamond",
                    0: "Other",
                    1: "hexagonal-ice",
                    2: "cubic-ice",
                }.get(expected, "")
                agree_rows.append({
                    "frame": frame,
                    "molecules": payload["molecules"],
                    "code": code,
                    "our_class": our_name,
                    "their_class": their_name,
                    "count": count,
                    "our_total": total,
                    "applies": "yes" if applies else "no",
                    "agree": "yes" if agreed else ("no" if applies else "na"),
                })
                log.append(
                    f"agree {frame} {code} {our_name} {their_name} {count} of {total}"
                )
                print(
                    f"{frame} {code} {our_name} -> {their_name} {count}/{total}",
                    flush=True,
                )
    write_time(args.out / "ptm-agreement.csv", agree_rows)
    (args.out / "raw-ptm-templates.txt").write_text("\n".join(log) + "\n", encoding="utf-8")
    return 0


if __name__ == "__main__":
    sys.exit(main())
