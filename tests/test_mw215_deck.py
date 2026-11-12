# Copyright (c) d-SEAMS developers

"""The 215 K deck is small. The liquid frame is written outside the tree."""

import os
import pathlib
import sys

ROOT = pathlib.Path(__file__).resolve().parents[1]
DECK = ROOT / "input" / "mw215" / "in.equilibrate"
FETCH = ROOT / "input" / "mw215" / "fetch.py"
FRAME = pathlib.Path(os.environ.get("MW215_DIR", "/var/tmp/dseams/mw215")) / "mw-liquid-512.lammpstrj"


def main() -> int:
    if len(sys.argv) > 1 and sys.argv[1] == "frame":
        if not FRAME.is_file():
            sys.stderr.write("skip missing 215 K frame\n")
            return 77
        if FRAME.stat().st_size < 1000:
            sys.stderr.write("frame is empty\n")
            return 1
        return 0
    text = DECK.read_text()
    if "velocity all create 500.0 215 " not in text:
        sys.stderr.write("velocity seed is not 215\n")
        return 1
    if "temp 215.0 215.0" not in text:
        sys.stderr.write("thermostat is not 215 K\n")
        return 1
    fetch = FETCH.read_text()
    if "d-SEAMS/mw-liquid-215K" not in fetch:
        sys.stderr.write("fetch script is missing the dataset id\n")
        return 1
    if "regenerate.py" not in fetch:
        sys.stderr.write("fetch script has no regenerate fallback\n")
        return 1
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
