#!/usr/bin/env python3
# Copyright (c) d-SEAMS developers

"""Regenerate homogeneous mW quenches and census the saved frames.

The dumps are large and stay outside the git work tree. Summaries written
by nucleation-pathway are the artefact that is small enough to keep.

n=46 is 97336 molecules. n=100 is 1000000 molecules. n=16 is 4096.
The deck is input/nucleation/in.homogeneous. The integrator step is
10 fs. Temperature defaults to 180 K.

A seeded run passes --seeds. Each seed is the LAMMPS velocity seed of
an independent 350 K liquid that is then quenched. One picosecond is
100 steps. The coordinate dumps are recorded by md5 in a DVC-style
manifest and are not committed.

Usage:
  python3 input/nucleation/regenerate.py --lmp /path/to/lmp --n 46 --n 100
  python3 input/nucleation/regenerate.py --lmp /path/to/lmp --n 16 \\
      --seeds 20261111-20261130 --melt 10000 --quench 100000 --dump-every 100
"""

from __future__ import annotations

import argparse
import hashlib
import os
import shutil
import subprocess
import time
from concurrent.futures import ThreadPoolExecutor, as_completed
from pathlib import Path

ROOT = Path(__file__).resolve().parents[2]
DECK = Path(__file__).resolve().parent
SW = ROOT / "input" / "mw215" / "mW.sw"


def parse_seeds(text: str) -> list[int]:
    seeds: list[int] = []
    for piece in text.split(","):
        piece = piece.strip()
        if not piece:
            continue
        if "-" in piece:
            lo, hi = piece.split("-", 1)
            seeds.extend(range(int(lo), int(hi) + 1))
        else:
            seeds.append(int(piece))
    if len(seeds) != len(set(seeds)):
        raise SystemExit(f"repeated seed in {text}")
    return seeds


def md5_file(path: Path) -> str:
    digest = hashlib.md5()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1 << 20), b""):
            digest.update(block)
    return digest.hexdigest()


def lammps_command(lmp: str, n: int, melt: int, quench: int, dump_every: int,
                   temperature: float, seed: int, nskip: int) -> list[str]:
    return [
        lmp,
        "-in",
        "in.homogeneous",
        "-var",
        "n",
        str(n),
        "-var",
        "nmelt",
        str(melt),
        "-var",
        "nquench",
        str(quench),
        "-var",
        "ndump",
        str(dump_every),
        "-var",
        "nskip",
        str(nskip),
        "-var",
        "Tq",
        str(temperature),
        "-var",
        "seed",
        str(seed),
        "-log",
        "log.lammps",
        "-screen",
        "none",
    ]


def run_one(lmp: str, n: int, melt: int, quench: int, dump_every: int,
            temperature: float, seed: int, nskip: int, work: Path, seams: str,
            threads: int, skip_md: bool) -> dict:
    work.mkdir(parents=True, exist_ok=True)
    env = os.environ.copy()
    env["OMP_NUM_THREADS"] = str(threads)
    lmp_wall = 0.0
    if not skip_md:
        shutil.copy(SW, work / "mW.sw")
        shutil.copy(DECK / "in.homogeneous", work / "in.homogeneous")
        started = time.perf_counter()
        subprocess.run(
            lammps_command(lmp, n, melt, quench, dump_every, temperature, seed, nskip),
            cwd=work,
            env=env,
            check=True,
        )
        lmp_wall = time.perf_counter() - started
    census_wall = 0.0
    frames = 0
    if seams:
        started = time.perf_counter()
        for dump, out, summary in (
            ("liquid.lammpstrj", "liquid.tsv", "liquid.summary.tsv"),
            ("quench.lammpstrj", "frames.tsv", "summary.tsv"),
        ):
            subprocess.run(
                [
                    seams,
                    str(work / dump),
                    "--out",
                    str(work / out),
                    "--summary",
                    str(work / summary),
                    "--type",
                    "1",
                    "--threads",
                    str(threads),
                ],
                cwd=ROOT,
                env=env,
                check=True,
            )
        census_wall = time.perf_counter() - started
        frames = max(0, (work / "frames.tsv").read_text().count("\n") - 1)
    quench_path = work / "quench.lammpstrj"
    liquid_path = work / "liquid.lammpstrj"
    return {
        "seed": seed,
        "molecules": n ** 3,
        "lmp_wall_s": lmp_wall,
        "census_wall_s": census_wall,
        "frames": frames,
        "quench_md5": md5_file(quench_path) if quench_path.is_file() else "",
        "quench_bytes": quench_path.stat().st_size if quench_path.is_file() else 0,
        "liquid_md5": md5_file(liquid_path) if liquid_path.is_file() else "",
        "liquid_bytes": liquid_path.stat().st_size if liquid_path.is_file() else 0,
        "work": str(work),
    }


def write_manifest(path: Path, rows: list[dict], command: str) -> None:
    lines = [
        "# DVC-style outs. The lammpstrj files are not in git.",
        f"# regenerate: {command}",
        "outs:",
    ]
    for row in rows:
        for kind in ("quench", "liquid"):
            lines.append(f"- path: {row['work']}/{kind}.lammpstrj")
            lines.append(f"  md5: {row[kind + '_md5']}")
            lines.append(f"  size: {row[kind + '_bytes']}")
            lines.append(f"  seed: {row['seed']}")
    path.write_text("\n".join(lines) + "\n")


def write_timing(path: Path, rows: list[dict], threads: int, jobs: int, cores: int,
                 wall: float) -> None:
    lines = [
        f"# cores\t{cores}",
        f"# jobs\t{jobs}",
        f"# threads_per_job\t{threads}",
        f"# wall_s\t{wall:.3f}",
        "seed\tmolecules\tlmp_wall_s\tcensus_wall_s\tframes\tquench_bytes",
    ]
    for row in rows:
        lines.append(
            f"{row['seed']}\t{row['molecules']}\t{row['lmp_wall_s']:.3f}\t"
            f"{row['census_wall_s']:.3f}\t{row['frames']}\t{row['quench_bytes']}"
        )
    path.write_text("\n".join(lines) + "\n")


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--lmp", default=os.environ.get("LMP", "lmp"))
    parser.add_argument("--n", type=int, action="append", default=[])
    parser.add_argument("--melt", type=int, default=2000)
    parser.add_argument("--quench", type=int, default=8000)
    parser.add_argument("--dump-every", type=int, default=2000)
    parser.add_argument("--nskip", type=int, default=0)
    parser.add_argument("--temperature", type=float, default=180.0)
    parser.add_argument("--seed", type=int, default=20261111)
    parser.add_argument("--seeds", default="")
    parser.add_argument(
        "--dest",
        default=os.environ.get("DSEAMS_NUCLEATION_DIR", "nucleation-frames"),
    )
    parser.add_argument("--seams", default=os.environ.get("NUCLEATION_PATHWAY", ""))
    parser.add_argument("--threads", type=int, default=1)
    parser.add_argument("--jobs", type=int, default=1)
    parser.add_argument("--skip-md", action="store_true")
    parser.add_argument("--manifest", type=Path, default=None)
    args = parser.parse_args()
    dest = Path(args.dest)
    dest.mkdir(parents=True, exist_ok=True)
    lmp = shutil.which(args.lmp) or args.lmp
    if not args.skip_md and not Path(lmp).is_file() and shutil.which(lmp) is None:
        raise SystemExit(f"LAMMPS executable not found: {args.lmp}")
    seeds = parse_seeds(args.seeds) if args.seeds else []
    if seeds:
        sizes = args.n or [16]
        if len(sizes) != 1:
            raise SystemExit("--seeds runs one system size")
        n = sizes[0]
        parent = dest / f"mw-{n ** 3}"
        parent.mkdir(parents=True, exist_ok=True)
        started = time.perf_counter()
        rows: list[dict] = []
        with ThreadPoolExecutor(max_workers=max(1, args.jobs)) as pool:
            futures = [
                pool.submit(
                    run_one,
                    lmp,
                    n,
                    args.melt,
                    args.quench,
                    args.dump_every,
                    args.temperature,
                    seed,
                    args.nskip,
                    parent / f"seed-{seed}",
                    args.seams,
                    args.threads,
                    args.skip_md,
                )
                for seed in seeds
            ]
            for future in as_completed(futures):
                row = future.result()
                rows.append(row)
                print(
                    f"seed {row['seed']} lmp {row['lmp_wall_s']:.1f}s "
                    f"census {row['census_wall_s']:.1f}s frames {row['frames']}",
                    flush=True,
                )
        rows.sort(key=lambda row: row["seed"])
        wall = time.perf_counter() - started
        cores = os.cpu_count() or 0
        write_timing(parent / "timing.tsv", rows, args.threads, args.jobs, cores, wall)
        manifest = args.manifest or (parent / "outs.dvc")
        command = " ".join(str(part) for part in [
            "python3 input/nucleation/regenerate.py",
            "--n", str(n),
            "--seeds", args.seeds,
            "--melt", str(args.melt),
            "--quench", str(args.quench),
            "--dump-every", str(args.dump_every),
            "--nskip", str(args.nskip),
            "--temperature", str(args.temperature),
            "--threads", str(args.threads),
            "--jobs", str(args.jobs),
        ])
        write_manifest(manifest, rows, command)
        print(f"wall {wall:.1f}s cores {cores} jobs {args.jobs} threads {args.threads}")
        print(manifest)
        return 0

    sizes = args.n or [46, 100]
    for n in sizes:
        molecules = n ** 3
        work = dest / f"mw-{molecules}"
        row = run_one(
            lmp,
            n,
            args.melt,
            args.quench,
            args.dump_every,
            args.temperature,
            args.seed,
            args.nskip,
            work,
            args.seams,
            args.threads,
            args.skip_md,
        )
        print(f"wrote {work} frames {row['frames']}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
