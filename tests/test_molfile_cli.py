#!/usr/bin/env python3
# Copyright (c) d-SEAMS developers
"""seams reads the VMD molfile fixtures through chemfiles.

    test_molfile_cli.py SEAMS INPUT_DIR

A suffix the CLI does not route to chemfiles is read as a LAMMPS dump and
gives "nop 0".
"""
import pathlib
import subprocess
import sys

EXPECTED = {
    "water.pqr": "nop 3 frame 1 box \n",
    "water.rst7": "nop 3 frame 1 box \n",
    "water.cube": "nop 3 frame 1 box 10 10 10\n",
    "water.xsf": "nop 3 frame 1 box 10 10 10\n",
}


def main():
    if len(sys.argv) != 3:
        sys.exit(__doc__)
    seams = sys.argv[1]
    root = pathlib.Path(sys.argv[2]) / "molfile"
    failed = False
    for name, want in EXPECTED.items():
        for extra in ([], ["--last", "2"]):
            got = subprocess.run(
                [seams, "read", str(root / name), "--type", "-1", *extra],
                capture_output=True, text=True, timeout=60)
            if got.returncode != 0 or got.stdout != want:
                sys.stderr.write("{} {}: exit {}, stdout {!r}, want {!r}\n{}".format(
                    name, " ".join(extra), got.returncode, got.stdout, want, got.stderr))
                failed = True
    return 1 if failed else 0


if __name__ == "__main__":
    sys.exit(main())
