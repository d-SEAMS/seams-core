#!/usr/bin/env python3
"""GenIce sI and sII cage census, then d-SEAMS and ClavisCage times.

The census is whatever scripts/rerun_main_text.py prints for the two
unit cells. The times are bench_cage_h2h (ring construction outside the
timer; hexagonal cages, double-diamond cages, and signatures 512,
51262, 51264) and ClavisCage 0.3.18 analyze_waters. One warm-up and
eight calls, eight threads, on each dump. The first of nine printed
bench samples is the warm-up.
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
import paper_stats as stats

DUMPS = {
    "sI": ROOT / "input/traj/genice_sI.lammpstrj",
    "sII": ROOT / "input/traj/genice_sII.lammpstrj",
}


def read_oxygens(path: Path):
    import numpy as np
    lines = path.read_text(encoding="utf-8").splitlines()
    if not lines[0].startswith("ITEM: TIMESTEP"):
        raise RuntimeError(f"{path.name} is not a LAMMPS dump")
    count = int(lines[3])
    bounds = []
    for row in lines[5:8]:
        lo, hi = row.split()[:2]
        bounds.append(float(hi) - float(lo))
    header = lines[8].split()[2:]
    type_at = header.index("type")
    xyz_at = [header.index(axis) for axis in ("x", "y", "z")]
    points = []
    for line in lines[9:9 + count]:
        parts = line.split()
        if int(parts[type_at]) != 1:
            continue
        points.append([float(parts[index]) for index in xyz_at])
    pos = np.ascontiguousarray(points, dtype=np.float64)
    lattice = np.diag(bounds)
    return {"positions": pos, "lattice": lattice, "molecules": len(pos)}


def bench_samples(binary: Path, dump: Path, threads: int, reps: int) -> tuple[list[float], str]:
    env = stats.pin_env(threads)
    build = binary.resolve().parents[1]
    env["LD_LIBRARY_PATH"] = ":".join(
        str(build / name) for name in (
            "src", "subprojects/readcon-core", "subprojects/readcon-db",
            "subprojects/linkcell", "subprojects/vesin", "subprojects/minimage"))
    cmd = [
        "taskset", "-c", stats.taskset_cores(threads),
        str(binary), str(dump.resolve().relative_to(ROOT)), "1", str(reps), "1",
    ]
    done = subprocess.run(cmd, cwd=ROOT, env=env, text=True, capture_output=True)
    if done.returncode != 0:
        tail = (done.stderr or done.stdout or "").strip().splitlines()
        raise RuntimeError(tail[-1] if tail else f"exit {done.returncode}")
    samples = [float(line.split()[-1]) for line in done.stdout.splitlines()
               if line.startswith("both_sample_ms ")]
    if len(samples) != reps:
        raise RuntimeError(f"{dump.name} printed {len(samples)} samples, need {reps}")
    return samples, done.stdout


def clavis_samples(frame: dict, threads: int, warmup: int, repeats: int) -> tuple[list[float], tuple, str]:
    import json
    import time
    import h2h_literature as lit
    import clavis_engine as engine

    version = engine.version()
    if not isinstance(version, str):
        version = str(version)
    if not version.startswith("0.3.18"):
        raise RuntimeError(f"ClavisCage {version}, need 0.3.18")
    text, box_nm, _fourth = lit.hydrogen_waters(frame)
    samples = []
    last = None
    for index in range(warmup + repeats):
        t0 = time.perf_counter()
        raw = engine.analyze_waters(text, "[]", box_nm, "trace", 0.35, 35.0)
        elapsed = (time.perf_counter() - t0) * 1000.0
        parsed = json.loads(raw)
        last = (
            len(parsed["cages_512"]),
            len(parsed["cages_51262"]),
            len(parsed["cages_51264"]),
        )
        if index >= warmup:
            samples.append(elapsed)
    return samples, last, version


def rerun(seams: Path, out_dir: Path) -> str:
    dest = out_dir / "rerun-main"
    dest.mkdir(parents=True, exist_ok=True)
    cmd = [
        sys.executable, str(ROOT / "scripts/rerun_main_text.py"),
        "--seams", str(seams.resolve().relative_to(ROOT)),
        "--input", "input/traj",
        "--archive", "input/traj/genice_sI.lammpstrj",
        "--source", ".",
        "--out", str(dest.resolve().relative_to(ROOT)),
    ]
    done = subprocess.run(cmd, cwd=ROOT, text=True, capture_output=True)
    report = dest / "report.txt"
    if not report.is_file():
        tail = (done.stderr or done.stdout or "").strip().splitlines()
        raise RuntimeError(tail[-1] if tail else "rerun_main_text wrote no report")
    return report.read_text(encoding="utf-8")


def census_rows(report: str) -> list[str]:
    found = []
    for line in report.splitlines():
        parts = line.split("\t")
        if len(parts) < 2:
            continue
        if parts[1].startswith("sI ") or parts[1].startswith("sII "):
            if not line.startswith("match\t"):
                raise RuntimeError(line)
            found.append(line)
    if len(found) < 4:
        raise RuntimeError("rerun did not record the sI and sII signatures")
    return found


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--out", type=Path, default=ROOT / "benchmarks/paper-cells")
    parser.add_argument("--bench", type=Path, default=ROOT / "bbdir-release/tests/bench_cage_h2h")
    parser.add_argument("--seams", type=Path, default=ROOT / "bbdir-release/seams")
    parser.add_argument("--threads", type=int, default=8)
    parser.add_argument("--repeats", type=int, default=8)
    parser.add_argument("--warmup", type=int, default=1)
    args = parser.parse_args()
    if args.repeats < 8:
        raise SystemExit("n must be at least 8")
    if not args.bench.is_file() or not args.seams.is_file():
        raise SystemExit("missing seams or bench_cage_h2h")
    args.out.mkdir(parents=True, exist_ok=True)
    os.environ.update(stats.pin_env(args.threads))
    try:
        os.sched_setaffinity(0, set(range(args.threads)))
    except (AttributeError, OSError):
        pass
    cores, affinity = stats.cores_now()
    report = rerun(args.seams, args.out)
    counted = census_rows(report)
    log = [
        *stats.machine_lines(ROOT),
        f"cores {cores}",
        f"affinity {affinity}",
        f"threads {args.threads}",
        f"warmup {args.warmup}",
        f"n {args.repeats}",
        "census scripts/rerun_main_text.py",
        "archive argument is input/traj/genice_sI.lammpstrj; "
        "the release tarball is not in this checkout",
        *counted,
    ]
    rows = [{
        "cell": "sI sII cage counts",
        "code": "d-SEAMS",
        "phase": "sI/sII",
        "molecules": "",
        "threads": args.threads,
        "cores": cores,
        "affinity": affinity,
        "warmup": "",
        "n": "",
        "median_ms": "",
        "iqr_ms": "",
        "samples_ms": "",
        "detail": "; ".join(line.replace("\t", " ") for line in counted),
    }]
    timed = []
    for phase, dump in DUMPS.items():
        frame = read_oxygens(dump)
        samples, text = bench_samples(
            args.bench, dump, args.threads, args.warmup + args.repeats)
        kept = samples[args.warmup:]
        summary = stats.summarize(kept)
        clath = ""
        for line in text.splitlines():
            if line.startswith("clathrate"):
                clath = line
        log.append(
            f"d-SEAMS {phase} molecules {frame['molecules']} {clath} "
            f"median_ms {summary['median_ms']:.6f} iqr_ms {summary['iqr_ms']:.6f} "
            f"samples {stats.fmt_samples(kept)}"
        )
        timed.append((phase, "d-SEAMS", frame["molecules"], summary, kept, clath))
        print(f"d-SEAMS {phase} {summary['median_ms']:.4f} ms", flush=True)
        samples, cages, version = clavis_samples(
            frame, args.threads, args.warmup, args.repeats)
        summary = stats.summarize(samples)
        detail = f"ClavisCage {version} TRACE 512={cages[0]} 51262={cages[1]} 51264={cages[2]}"
        log.append(
            f"{detail} {phase} median_ms {summary['median_ms']:.6f} "
            f"iqr_ms {summary['iqr_ms']:.6f} samples {stats.fmt_samples(samples)}"
        )
        timed.append((phase, "ClavisCage", frame["molecules"], summary, samples, detail))
        print(f"ClavisCage {phase} {summary['median_ms']:.4f} ms {detail}", flush=True)
    for phase, code, molecules, summary, samples, detail in timed:
        rows.append({
            "cell": "sI sII cages" if code == "d-SEAMS" else "sI sII ClavisCage",
            "code": code,
            "phase": phase,
            "molecules": molecules,
            "threads": args.threads,
            "cores": cores,
            "affinity": affinity,
            "warmup": args.warmup,
            "n": summary["n"],
            "median_ms": f"{summary['median_ms']:.6f}",
            "iqr_ms": f"{summary['iqr_ms']:.6f}",
            "samples_ms": stats.fmt_samples(samples),
            "detail": detail,
        })
    (args.out / "raw-cages.txt").write_text("\n".join(log) + "\n", encoding="utf-8")
    frag = args.out / "frag-lit-cages.csv"
    with frag.open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(rows[0]))
        writer.writeheader()
        writer.writerows(rows)
    print(f"wrote {frag}", flush=True)
    return 0


if __name__ == "__main__":
    import os
    sys.exit(main())
