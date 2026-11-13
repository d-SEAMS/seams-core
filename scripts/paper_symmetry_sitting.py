#!/usr/bin/env python3
"""One sitting of the symmetry classifier against both PTM template sets.

Ideal ice Ic, 110592 molecules, one and eight threads. bench_v2 mode
symmetry counts every atom (sink) and the ice labels (crystal). PTM is
the default set (FCC, HCP, BCC) and the ice-only set (cubic and
hexagonal diamond, RMSD cutoff 0.1). The ratios are medians from this
run.

The interface check reads an equilibrated oxygen dump. There is no
230 K brine trajectory in the tree, so the frame is the TIP4P/2005
ice/water slab after the hop equilibration. CHILL+ atoms that are
hexagonal while the classifier says other are binned by the first
kernel gate they fail, and by how many of their four nearest oxygens
the classifier calls ice.
"""

from __future__ import annotations

import argparse
import csv
import os
import subprocess
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "scripts"))
import paper_ptm_templates as templates
import paper_stats as stats

MOLECULES = 110592
N_ICE = 2048


def bench_env(threads: int) -> dict[str, str]:
    env = stats.pin_env(threads)
    build = ROOT / "bbdir-release"
    env["LD_LIBRARY_PATH"] = os.pathsep.join(
        str(build / name)
        for name in (
            "src",
            "subprojects/readcon-core",
            "subprojects/readcon-db",
            "subprojects/linkcell",
            "subprojects/vesin",
            "subprojects/minimage",
        )
    )
    return env


def time_classifier(binary: Path, threads: int, repeats: int) -> tuple[list[float], int, int]:
    cmd = [
        "taskset",
        "-c",
        stats.taskset_cores(threads),
        str(binary),
        "symmetry",
        str(MOLECULES),
        str(threads),
        str(repeats),
    ]
    done = subprocess.run(cmd, cwd=ROOT, env=bench_env(threads), text=True, capture_output=True)
    if done.returncode != 0:
        tail = (done.stderr or done.stdout or "").strip().splitlines()
        raise RuntimeError(tail[-1] if tail else f"exit {done.returncode}")
    line = [row for row in done.stdout.splitlines() if row.startswith("mode=")][-1]
    fields = {}
    for token in line.split():
        if "=" in token and not token.startswith("samples_ms="):
            key, value = token.split("=", 1)
            fields[key] = value
    samples_token = line.split("samples_ms=", 1)[1]
    samples = [float(piece) for piece in samples_token.split(",") if piece]
    if len(samples) != repeats:
        raise RuntimeError(f"classifier returned {len(samples)} samples")
    return samples, int(fields["sink"]), int(fields["crystal"])


def time_ptm(kind: str, threads: int, warmup: int, repeats: int) -> list[float]:
    cmd = [
        "taskset",
        "-c",
        stats.taskset_cores(threads),
        sys.executable,
        str(Path(templates.__file__).resolve()),
        "--worker",
        "time",
        "--template",
        kind,
        "--threads",
        str(threads),
        "--warmup",
        str(warmup),
        "--repeats",
        str(repeats),
        "--lib",
        str(ROOT / "bbdir-release/src/libyodaLib.so"),
    ]
    env = stats.pin_env(threads)
    env["LD_LIBRARY_PATH"] = templates.lib_path()
    env["PYTHONPATH"] = os.environ.get("PYTHONPATH", "")
    env["OVITO_THREAD_COUNT"] = str(threads)
    env["QT_QPA_PLATFORM"] = "offscreen"
    done = subprocess.run(cmd, cwd=ROOT, env=env, text=True, capture_output=True)
    if done.returncode != 0:
        tail = (done.stderr or done.stdout or "").strip().splitlines()
        raise RuntimeError(tail[-1] if tail else f"exit {done.returncode}")
    text = [line for line in done.stdout.splitlines() if line.strip()][-1]
    samples = [float(piece) for piece in text.split()]
    if len(samples) != repeats:
        raise RuntimeError(f"PTM {kind} returned {len(samples)} samples")
    return samples


def write_csv(path: Path, rows: list[dict]) -> None:
    with path.open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(rows[0]))
        writer.writeheader()
        writer.writerows(rows)


def run_timing(args) -> None:
    binary = Path(args.bench)
    if not binary.is_file():
        raise SystemExit(f"missing {binary}")
    args.out.mkdir(parents=True, exist_ok=True)
    measured = []
    log = [
        *stats.machine_lines(ROOT),
        "sitting one",
        f"warmup {args.warmup}",
        f"n {args.repeats}",
        "frame ideal ice Ic 110592",
        "checksum sink is the atom count, crystal is labels >= 0",
        "PTM default is FCC HCP BCC; ice is cubic and hexagonal diamond, rmsd 0.1",
    ]
    for threads in (1, 8):
        samples, sink, crystal = time_classifier(binary, threads, args.repeats)
        if sink != MOLECULES:
            raise SystemExit(f"classifier sink {sink} is not every atom")
        summary = stats.summarize(samples)
        measured.append(("classifier", "classifier", threads, summary, samples, sink, crystal))
        log.append(
            f"classifier threads {threads} median_ms {summary['median_ms']:.6f} "
            f"iqr_ms {summary['iqr_ms']:.6f} sink {sink} crystal {crystal} "
            f"samples {stats.fmt_samples(samples)}"
        )
        print(
            f"classifier threads {threads} {summary['median_ms']:.4f} ms "
            f"sink {sink} crystal {crystal}",
            flush=True,
        )
        for kind in ("default", "ice"):
            samples = time_ptm(kind, threads, args.warmup, args.repeats)
            summary = stats.summarize(samples)
            measured.append((f"PTM-{kind}", kind, threads, summary, samples, "", ""))
            log.append(
                f"ptm {kind} threads {threads} median_ms {summary['median_ms']:.6f} "
                f"iqr_ms {summary['iqr_ms']:.6f} samples {stats.fmt_samples(samples)}"
            )
            print(f"PTM {kind} threads {threads} {summary['median_ms']:.4f} ms", flush=True)
    by_thread = {}
    for _code, _kind, threads, summary, _samples, _sink, _crystal in measured:
        if _code == "classifier":
            by_thread[threads] = summary["median_ms"]
    rows = []
    for code, kind, threads, summary, samples, sink, crystal in measured:
        classifier_ms = by_thread[threads]
        rows.append({
            "code": code,
            "template_set": kind,
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
            "sink": sink,
            "crystal": crystal,
            "ratio_this_over_classifier": f"{summary['median_ms'] / classifier_ms:.6f}",
            "ratio_classifier_over_this": f"{classifier_ms / summary['median_ms']:.6f}",
            "samples_ms": stats.fmt_samples(samples),
        })
        log.append(
            f"ratio {code} threads {threads} "
            f"this_over_classifier {summary['median_ms'] / classifier_ms:.6f} "
            f"classifier_over_this {classifier_ms / summary['median_ms']:.6f}"
        )
    write_csv(args.out / "symmetry-vs-ptm.csv", rows)
    (args.out / "raw-symmetry-vs-ptm.txt").write_text("\n".join(log) + "\n", encoding="utf-8")
    print(f"wrote {args.out / 'symmetry-vs-ptm.csv'}", flush=True)


def read_oxygens(path: Path):
    import numpy as np

    oxy_id = []
    oxy = []
    bounds = []
    with path.open() as handle:
        mode = None
        for line in handle:
            text = line.strip()
            if text.startswith("ITEM: TIMESTEP"):
                if oxy:
                    break
                mode = None
                continue
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
            if mode == "box" and len(bounds) < 3:
                lo, hi = text.split()[:2]
                bounds.append(float(hi) - float(lo))
            elif mode == "atoms":
                parts = text.split()
                if len(parts) < 5 or int(parts[1]) != 1:
                    continue
                oxy_id.append(int(parts[0]))
                oxy.append((float(parts[2]), float(parts[3]), float(parts[4])))
    pos = np.ascontiguousarray(oxy, dtype=float)
    lattice = [
        [bounds[0], 0.0, 0.0],
        [0.0, bounds[1], 0.0],
        [0.0, 0.0, bounds[2]],
    ]
    return np.asarray(oxy_id, dtype=int), pos, lattice


def mic(a, b, box):
    d = b - a
    d -= box * np_round(d / box)
    return d


def np_round(values):
    import numpy as np

    return np.round(values)


def tetra_q(vectors) -> float:
    import numpy as np

    acc = 0.0
    for i in range(4):
        vi = vectors[i]
        ni = float(np.linalg.norm(vi))
        for j in range(i + 1, 4):
            vj = vectors[j]
            nj = float(np.linalg.norm(vj))
            cosine = float(np.dot(vi, vj) / (ni * nj))
            acc += (cosine + 1.0 / 3.0) ** 2
    return 1.0 - (3.0 / 8.0) * acc


def neighbour_table(pos, box):
    import numpy as np

    n = len(pos)
    box = np.asarray(box, dtype=float)
    nearest = []
    within = []
    for i in range(n):
        delta = pos - pos[i]
        delta -= box * np.round(delta / box)
        dist = np.sqrt(np.sum(delta * delta, axis=1))
        dist[i] = 1e9
        order = np.argpartition(dist, 4)[:4]
        order = order[np.argsort(dist[order])]
        nearest.append([(float(dist[j]), int(j), delta[j]) for j in order])
        within.append(int(np.count_nonzero(dist < 3.5)))
    return nearest, within


def gate_of(nearest_i, within_i, labels) -> str:
    if within_i < 4:
        return "fewer-than-4-within-3.5"
    lengths = [item[0] for item in nearest_i]
    if min(lengths) < 2.30 or max(lengths) > 3.35:
        return "four-nearest-outside-2.30-3.35"
    vectors = [item[2] for item in nearest_i]
    if not (max(lengths) <= min(lengths) * 1.00001) and tetra_q(vectors) < 0.94:
        return "tetrahedral-q-below-0.94"
    crystal = sum(1 for _dist, j, _vec in nearest_i if int(labels[j]) >= 0)
    if crystal == 0:
        return "no-crystal-neighbour"
    return "passes-early-gates"


def run_interface(args) -> None:
    import numpy as np

    path = Path(args.interface)
    if not path.is_file():
        raise SystemExit(f"missing {path}")
    ids, pos, lattice = read_oxygens(path)
    os.environ["OVITO_THREAD_COUNT"] = "8"
    os.environ["QT_QPA_PLATFORM"] = "offscreen"
    os.environ["OMP_NUM_THREADS"] = "8"
    os.environ.setdefault("LD_LIBRARY_PATH", templates.lib_path())
    lib = templates.open_core(args.lib)
    ours = templates.classify(lib, pos, lattice)
    ptm = templates.ovito_labels(pos, lattice, "ice")
    chill = templates.ovito_labels(pos, lattice, "chill")
    box = np.array([lattice[0][0], lattice[1][1], lattice[2][2]], dtype=float)
    nearest, within = neighbour_table(pos, box)
    mol = (ids - 1) // 3 + 1
    ice_origin = mol <= N_ICE
    extra = (ours == -1) & (chill == 1)
    gates = {"fewer-than-4-within-3.5": 0, "four-nearest-outside-2.30-3.35": 0,
             "tetrahedral-q-below-0.94": 0, "no-crystal-neighbour": 0,
             "passes-early-gates": 0}
    origin_ice = 0
    crystal_neighbour_hist = [0, 0, 0, 0, 0]
    for index in np.flatnonzero(extra):
        gates[gate_of(nearest[index], within[index], ours)] += 1
        origin_ice += int(ice_origin[index])
        crystal_n = sum(1 for _d, j, _v in nearest[index] if int(ours[j]) >= 0)
        crystal_neighbour_hist[crystal_n] += 1
    crystal_z = pos[ours >= 0, 2]
    surface = 0
    if len(crystal_z):
        z0 = float(crystal_z.min())
        z1 = float(crystal_z.max())
        band = 3.5
        for index in np.flatnonzero(extra):
            z = float(pos[index, 2])
            if min(abs(z - z0), abs(z - z1)) <= band:
                surface += 1
    rows = []
    log = [
        *stats.machine_lines(ROOT),
        f"frame {path.name}",
        "source TIP4P/2005 ice/water slab after hop equilibration",
        "no 230 K brine trajectory is in the tree",
        f"oxygens {len(pos)}",
        f"original ice molecules {int(np.count_nonzero(ice_origin))}",
        f"original liquid molecules {int(np.count_nonzero(~ice_origin))}",
        "CHILL+ hexagonal is type 1; interfacial ice is type 3",
        f"extra hexagonal (classifier other, CHILL+ hexagonal) {int(np.count_nonzero(extra))}",
        f"of those, original ice {origin_ice}",
        f"of those, within 3.5 A of the crystal z edge {surface}",
        "crystal-neighbour count among the four nearest "
        + " ".join(f"{n}:{crystal_neighbour_hist[n]}" for n in range(5)),
    ]
    for name, count in gates.items():
        log.append(f"gate {name} {count}")
    for code, their, names in (
        ("PTM", ptm, templates.PTM_NAME),
        ("CHILL+", chill, templates.CHILL_NAME),
    ):
        for our_name, their_name, count, total in templates.crosstab(ours, their, names):
            rows.append({
                "frame": "tip4p-equilibrated-interface",
                "molecules": len(pos),
                "code": code,
                "our_class": our_name,
                "their_class": their_name,
                "count": count,
                "our_total": total,
            })
            log.append(f"agree {code} {our_name} {their_name} {count} of {total}")
            print(f"{code} {our_name} -> {their_name} {count}/{total}", flush=True)
    for name, count in gates.items():
        rows.append({
            "frame": "tip4p-equilibrated-interface",
            "molecules": len(pos),
            "code": "CHILL+-extra-hexagonal",
            "our_class": "other",
            "their_class": name,
            "count": count,
            "our_total": int(np.count_nonzero(extra)),
        })
        print(f"extra {name} {count}", flush=True)
    args.out.mkdir(parents=True, exist_ok=True)
    write_csv(args.out / "ptm-interface-agreement.csv", rows)
    (args.out / "raw-interface-labels.txt").write_text("\n".join(log) + "\n", encoding="utf-8")
    print(f"wrote {args.out / 'ptm-interface-agreement.csv'}", flush=True)


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--out", type=Path, default=ROOT / "benchmarks/paper-cells")
    parser.add_argument("--bench", default=str(ROOT / "bbdir-release/tests/bench_v2"))
    parser.add_argument("--lib", default=str(ROOT / "bbdir-release/src/libyodaLib.so"))
    parser.add_argument("--warmup", type=int, default=1)
    parser.add_argument("--repeats", type=int, default=8)
    parser.add_argument("--timing", action="store_true")
    parser.add_argument("--interface", type=Path, default=None)
    args = parser.parse_args()
    if args.repeats < 8:
        raise SystemExit("n must be at least 8")
    if args.timing:
        run_timing(args)
    if args.interface is not None:
        run_interface(args)
    if not args.timing and args.interface is None:
        raise SystemExit("pass --timing and/or --interface")
    return 0


if __name__ == "__main__":
    sys.exit(main())
