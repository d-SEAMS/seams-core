# Copyright (c) d-SEAMS developers

"""The 215 K deck is small. The liquid frame is written outside the tree."""

import os
import pathlib
import shutil
import sys

ROOT = pathlib.Path(__file__).resolve().parents[1]
DECK = ROOT / "input" / "mw215" / "in.equilibrate"
FETCH = ROOT / "input" / "mw215" / "fetch.py"
FRAME = pathlib.Path(os.environ.get("MW215_DIR", "mw215-frames")) / "mw-liquid-512.lammpstrj"


def main() -> int:
    if len(sys.argv) > 1 and sys.argv[1] == "frame":
        if not FRAME.is_file():
            # The liquid dump is produced by the deck, not stored in git.
            if shutil.which("lmp") is None:
                sys.stderr.write("lmp is not on PATH\n")
                return 77
            sys.path.insert(0, str(ROOT / "input" / "mw215"))
            from regenerate import generate

            try:
                generate(8)
            except FileNotFoundError as exc:
                sys.stderr.write("lmp is not runnable: {}\n".format(exc))
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
