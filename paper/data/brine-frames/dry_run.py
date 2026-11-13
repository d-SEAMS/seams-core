#!/usr/bin/env python3
"""Dry-run the brine-frame rules and require every labelled frame.

snakemake -n has to resolve all twelve
paper/data/brine-frames/out/T230_m{0,20,40}_r0_{t0,t1,t2,tend}.xyz targets.
Placeholder restarts live under paper/data/brine-frames/.dry-run/, which
is gitignored. Pass a real dataset with --config to run the jobs.
"""

from __future__ import annotations

import shutil
import subprocess
import sys
from pathlib import Path

HERE = Path(__file__).resolve().parent
REPO = HERE.parents[2]
PAIRS = (0, 20, 40)
STEPS = {"t0": 0, "t1": 830000, "t2": 1670000, "tend": 2500000}


def main() -> None:
    snakemake = shutil.which("snakemake")
    if snakemake is None:
        raise SystemExit("snakemake is not on PATH")
    root = Path("paper/data/brine-frames/.dry-run/d-SEAMS/brine-frames")
    restarts = root / "data" / "restarts"
    frames = root / "data" / "frames"
    (REPO / frames).mkdir(parents=True, exist_ok=True)
    for pairs in PAIRS:
        dest = REPO / restarts / f"T230_m{pairs}_r0"
        dest.mkdir(parents=True, exist_ok=True)
        for step in STEPS.values():
            path = dest / f"step-{step}.restart"
            path.write_bytes(b"placeholder\n")
    targets = [
        f"paper/data/brine-frames/out/T230_m{pairs}_r0_{name}.xyz"
        for pairs in PAIRS
        for name in STEPS
    ]
    cmd = [
        snakemake,
        "-s",
        "paper/data/brine-frames/Snakefile",
        "-n",
        *targets,
        "--config",
        f"dataset={root}",
        f"restarts={restarts}",
        f"frames={frames}",
    ]
    proc = subprocess.run(cmd, cwd=REPO, text=True, capture_output=True)
    sys.stdout.write(proc.stdout)
    sys.stderr.write(proc.stderr)
    if proc.returncode != 0:
        raise SystemExit(proc.returncode)
    text = proc.stdout + proc.stderr
    missing = [target for target in targets if target not in text]
    if missing:
        raise SystemExit("dry-run did not name:\n" + "\n".join(missing))
    print(f"resolved {len(targets)} frame targets")


if __name__ == "__main__":
    main()
