#!/usr/bin/env python3
"""Run one paper tutorial and check the numbers it prints."""

import subprocess
import sys
from pathlib import Path

EXPECT = {
    "interface": [
        "interface published hexagonal 68 interfacial 19 water 41",
        "interface ovito hexagonal 68 interfacial 20 water 40",
    ],
    "cages": [
        "cages sI 2 6",
        "cages sII 16 8",
        "cages sH 6 4 2",
    ],
    "rings": [
        "rings primitive 4 24 6 32",
        "rings shortest 4 24 6 32",
        "rings searched 4 24 6 128",
    ],
    "brine": [
        "brine frame 0 hexagonal 31 water 38 front 7.40",
        "brine frame 1 hexagonal 32 water 36 front 7.40",
        "brine frame 2 hexagonal 32 water 35 front 7.40",
    ],
    "nucleation": [
        "nucleation frame 0 q6bar 0.276023",
        "nucleation frame 1 q6bar 0.320504",
        "nucleation frame 2 q6bar 0.334254",
    ],
}
FIGURE = {
    "interface": "interface.png",
    "cages": "cages.png",
    "rings": "rings.png",
    "brine": "brine.png",
    "nucleation": "nucleation.png",
}


def main():
    topic, root, seams, libdir, build = sys.argv[1:6]
    out = Path(build) / "paper-tutorials" / topic
    script = Path(root) / "docs" / "tutorials" / "paper_science.py"
    done = subprocess.run(
        [sys.executable, str(script), topic, seams, libdir, str(out)],
        check=False,
        capture_output=True,
        text=True,
    )
    if done.returncode == 77:
        print(done.stdout.strip())
        raise SystemExit(77)
    if done.returncode != 0:
        sys.stderr.write(done.stdout)
        sys.stderr.write(done.stderr)
        raise SystemExit(done.returncode)
    missing = [line for line in EXPECT[topic] if line not in done.stdout]
    if missing:
        sys.stderr.write(done.stdout)
        raise SystemExit("missing lines:\n" + "\n".join(missing))
    figure = out / FIGURE[topic]
    if not figure.is_file() or figure.stat().st_size > 300 * 1024:
        raise SystemExit("figure {} failed".format(figure.name))
    print(done.stdout.strip())


if __name__ == "__main__":
    main()
