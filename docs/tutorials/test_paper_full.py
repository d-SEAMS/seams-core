#!/usr/bin/env python3
"""Run one full paper tutorial through its Snakemake rule and check the report."""

import os
import shutil
import subprocess
import sys
from pathlib import Path

EXPECT = {
    "interface": [
        "interface full input molecules 864 frames 1 file-frames 1 file input/traj/interface_f0.lammpstrj",
        "interface full paper molecules 864 frames 1",
        "interface full published hexagonal 352 interfacial 49 water 434 clathrate 3 interClathrate 25",
        "interface full ovito hexagonal 346 interfacial 57 water 437 clathrate 3 interClathrate 20",
        "interface full phi-c 0.00283286",
        "interface full stack HHHHHHHH....C.",
    ],
    "cages": [
        "cages full input sI molecules 368 frames 1 file sI-full.lammpstrj",
        "cages full input sII molecules 544 frames 1 file sII-full.lammpstrj",
        "cages full input sH molecules 272 frames 1 file sH-full.lammpstrj",
        "cages full sI 16 48",
        "cages full sII 64 32",
        "cages full sH 24 16 8",
        "cages full paper unit sI molecules 46 cages 2 6",
        "cages full paper unit sII molecules 136 cages 16 8",
        "cages full paper unit sH molecules 68 cages 6 4 2",
    ],
    "rings": [
        "rings full input molecules 4096 frames 1 file-frames 11 file input/traj/mW_cubic.lammpstrj",
        "rings full paper molecules 4096 frames 1 rings 8192",
        "rings full nop 4096 rings 8192",
        "rings full histogram 6 8192",
    ],
    "brine": [
        "brine full input molecules 873 frames 1 file-frames 1 file input/traj/nacl_front.lammpstrj",
        "brine full paper waters 864 ions 9 frames 1",
        "brine full waters 864 ions 9 hexagonal 352 water 434",
    ],
    "nucleation": [
        "nucleation full input molecules 4096 frames 148 file-frames 148 file docs/.cache/figshare-11448702/nucleation.lammpstrj",
        "nucleation full paper molecules 4096 frames 148",
        "nucleation full frame 1 q6bar 0.26063",
        "nucleation full frame 148 q6bar 0.287256",
        "nucleation full lines 148",
    ],
}
FIGURE = {
    "interface": "interface-full.png",
    "cages": "cages-full.png",
    "rings": "rings-full.png",
    "brine": "brine-full.png",
    "nucleation": "nucleation-full.png",
}


def snakemake_bin():
    sibling = Path(sys.executable).resolve().parent / "snakemake"
    if sibling.is_file():
        return str(sibling)
    return shutil.which("snakemake")


def main():
    topic, root, seams, libdir, build = sys.argv[1:6]
    # Off unless the build asked for the long Snakemake paper rules.
    if os.environ.get("SEAMS_PAPER_FULL") != "1":
        print(
            "full paper runs are opt-in: meson setup -Dpaper_full=true",
            file=sys.stderr,
        )
        raise SystemExit(77)
    program = snakemake_bin()
    if not program:
        print("snakemake is absent")
        raise SystemExit(77)
    root = Path(root)
    work = Path(build) / "paper-full"
    work.mkdir(parents=True, exist_ok=True)
    out = work / topic
    if out.exists():
        shutil.rmtree(out)
    env = os.environ.copy()
    env["SEAMS"] = seams
    env["SEAMS_LIBDIR"] = libdir
    env["PAPER_OUT"] = str(Path(build) / "paper-full")
    env["PYTHON"] = sys.executable
    env["NO_COLOR"] = "1"
    env["OPENMM_CPU_THREADS"] = "1"
    done = subprocess.run(
        [
            program,
            "-s",
            str(root / "docs" / "tutorials" / "Snakefile"),
            "--cores",
            "1",
            "--force",
            topic,
        ],
        check=False,
        capture_output=True,
        text=True,
        env=env,
        cwd=str(work),
    )
    blob = done.stdout + "\n" + done.stderr
    if done.returncode != 0:
        if "absent" in blob:
            print(blob.strip())
            raise SystemExit(77)
        sys.stderr.write(blob)
        raise SystemExit(done.returncode or 1)
    report = (out / "report.txt").read_text()
    missing = [line for line in EXPECT[topic] if line not in report]
    if missing:
        sys.stderr.write(report)
        raise SystemExit("missing lines:\n" + "\n".join(missing))
    figure = out / FIGURE[topic]
    if not figure.is_file() or figure.stat().st_size > 300 * 1024:
        raise SystemExit("figure {} failed".format(figure.name))
    print(report.strip())


if __name__ == "__main__":
    main()
