#!/usr/bin/env python3
# Copyright (c) d-SEAMS developers

"""Regenerate short homogeneous mW quenches and census the saved frames.

The dumps are large and stay outside the git work tree. Summaries written
by nucleation-pathway are the artefact that is small enough to keep.

n=46 is 97336 molecules. n=100 is 1000000 molecules. Both fit the
10^5 to 10^6 window. The deck is input/nucleation/in.homogeneous.
The integrator step is 10 fs. One quench at 180 K is a few tens of
picoseconds, which is the run that fits on 8 cores.

Usage:
  python3 input/nucleation/regenerate.py --lmp /path/to/lmp --n 46 --n 100
"""

from __future__ import annotations

import argparse
import os
import shutil
import subprocess
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[2]
DECK = Path(__file__).resolve().parent
SW = ROOT / "input" / "mw215" / "mW.sw"


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--lmp", default=os.environ.get("LMP", "lmp"))
    parser.add_argument("--n", type=int, action="append", default=[])
    parser.add_argument("--melt", type=int, default=2000)
    parser.add_argument("--quench", type=int, default=8000)
    parser.add_argument("--dump-every", type=int, default=2000)
    parser.add_argument("--temperature", type=float, default=180.0)
    parser.add_argument(
        "--dest",
        default=os.environ.get("DSEAMS_NUCLEATION_DIR", "/tmp/dseams-nucleation"),
    )
    parser.add_argument("--seams", default=os.environ.get("NUCLEATION_PATHWAY", ""))
    parser.add_argument("--threads", type=int, default=8)
    parser.add_argument("--skip-md", action="store_true")
    args = parser.parse_args()
    sizes = args.n or [46, 100]
    dest = Path(args.dest)
    dest.mkdir(parents=True, exist_ok=True)
    lmp = shutil.which(args.lmp) or args.lmp
    if not args.skip_md and not Path(lmp).is_file() and shutil.which(lmp) is None:
        print(f"LAMMPS executable not found: {args.lmp}", file=sys.stderr)
        return 1
    for n in sizes:
        molecules = n ** 3
        work = dest / f"mw-{molecules}"
        work.mkdir(parents=True, exist_ok=True)
        dump = work / "quench.lammpstrj"
        if not args.skip_md:
            shutil.copy(SW, work / "mW.sw")
            shutil.copy(DECK / "in.homogeneous", work / "in.homogeneous")
            env = os.environ.copy()
            env["OMP_NUM_THREADS"] = str(args.threads)
            cmd = [
                lmp,
                "-in",
                "in.homogeneous",
                "-var",
                "n",
                str(n),
                "-var",
                "nmelt",
                str(args.melt),
                "-var",
                "nquench",
                str(args.quench),
                "-var",
                "ndump",
                str(args.dump_every),
                "-var",
                "Tq",
                str(args.temperature),
                "-log",
                "log.lammps",
            ]
            print(" ".join(cmd), flush=True)
            subprocess.run(cmd, cwd=work, env=env, check=True)
        if args.seams:
            summary = work / "summary.tsv"
            table = work / "frames.tsv"
            env = os.environ.copy()
            env["OMP_NUM_THREADS"] = str(args.threads)
            subprocess.run(
                [
                    args.seams,
                    str(dump),
                    "--out",
                    str(table),
                    "--summary",
                    str(summary),
                    "--type",
                    "1",
                    "--threads",
                    str(args.threads),
                ],
                cwd=ROOT,
                env=env,
                check=True,
            )
            print(f"wrote {summary}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
