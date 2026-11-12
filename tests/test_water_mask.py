# Copyright (c) d-SEAMS developers
"""Cutoff rings keep every --water-types pair, including cross types."""

import math
import subprocess
import sys
import tempfile
from pathlib import Path


def main():
    seams = sys.argv[1]
    side = 2.0
    diag = side * math.sqrt(2.0)
    cutoff = (side + diag) / 2.0
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
        "2 2 2 2 0 0",
        "3 3 1 2 2 0",
        "4 4 2 0 2 0",
    ]
    with tempfile.TemporaryDirectory() as tmp:
        path = Path(tmp) / "cross.lammpstrj"
        path.write_text("\n".join(lines) + "\n")

        def rings(types):
            done = subprocess.run(
                [
                    seams,
                    "--water-types",
                    types,
                    "--cutoff",
                    str(cutoff),
                    "--max-ring",
                    "4",
                    "rings",
                    str(path),
                ],
                check=False,
                capture_output=True,
                text=True,
            )
            return done

        both = rings("1,2")
        only = rings("1")
        if both.returncode != 0 or only.returncode != 0:
            sys.stderr.write(both.stdout + both.stderr + only.stdout + only.stderr)
            return 1
        if "4 1" not in both.stdout.splitlines():
            sys.stderr.write(both.stdout)
            return 1
        if "4 1" in only.stdout.splitlines():
            sys.stderr.write(only.stdout)
            return 1
    return 0


if __name__ == "__main__":
    sys.exit(main())
