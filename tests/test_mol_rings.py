# Copyright (c) d-SEAMS developers
"""seams --mol water rings keeps the oxygen and drops the ion."""

import subprocess
import sys
import tempfile
from pathlib import Path


def main():
    seams = sys.argv[1]
    lines = [
        "ITEM: TIMESTEP",
        "0",
        "ITEM: NUMBER OF ATOMS",
        "4",
        "ITEM: BOX BOUNDS pp pp pp",
        "0 20",
        "0 20",
        "0 20",
        "ITEM: ATOMS id mol type x y z",
        "1 1 1 0 0 0",
        "2 1 2 0.96 0 0",
        "3 1 2 -0.24 0.93 0",
        "4 2 3 8 0 0",
    ]
    with tempfile.TemporaryDirectory() as tmp:
        path = Path(tmp) / "ion.lammpstrj"
        path.write_text("\n".join(lines) + "\n")
        done = subprocess.run(
            [seams, "--mol", "water", "rings", str(path)],
            check=False,
            capture_output=True,
            text=True,
        )
        if done.returncode != 0:
            sys.stderr.write(done.stdout + done.stderr)
            return 1
        first = done.stdout.splitlines()[0]
        if not first.startswith("nop 1 "):
            sys.stderr.write(done.stdout)
            return 1
    return 0


if __name__ == "__main__":
    sys.exit(main())
