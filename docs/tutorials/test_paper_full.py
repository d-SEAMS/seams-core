#!/usr/bin/env python3
"""Run one full paper tutorial through its Snakemake rule and check the report."""

import os
import shutil
import subprocess
import sys
from pathlib import Path

EXPECT = {
    "interface": [
        "interface full waters 2048",
        "interface full published hexagonal 840 interfacial 8 water 1131 clathrate 16 interClathrate 53",
        "interface full ovito hexagonal 807 interfacial 43 water 1140 clathrate 17 interClathrate 41",
        "interface full phi-c 0",
        "interface full stack HHHHHHHHHH..........",
    ],
    "cages": [
        "cages full sI 16 48",
        "cages full sII 64 32",
        "cages full sH 24 16 8",
    ],
    "rings": [
        "rings full nop 4096 rings 8192",
        "rings full histogram 6 8192",
    ],
    "brine": [
        "brine full frame 0 hexagonal 93 water 152 front 18.50 stack HHHHH.......",
        "brine full frame 1 hexagonal 160 water 86 front 25.90 stack HHHHHHH.....",
        "brine full frame 2 hexagonal 218 water 38 front 37.00 stack HHHHHHHHHH..",
    ],
    "nucleation": [
        "nucleation full frame 0 q6bar 0.26828",
        "nucleation full frame 1 q6bar 0.287644",
        "nucleation full frame 2 q6bar 0.56248",
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
