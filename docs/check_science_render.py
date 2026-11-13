#!/usr/bin/env python3
"""The vendored renderer must match science/interface-nucleation f65657b9."""

import hashlib
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
SOURCE = ROOT / "scripts" / "science_render.py"
MARKER = b"# END-HEADER\n"
# sha256 of scripts/science_render.py at f65657b9b89481e7dd6617b1987bc58b2ec0f01c.
EXPECTED = "c39425bc9a993261389521a7fb3c58117762938cb1b547d6c1d8b68c90105465"


def main():
    text = SOURCE.read_bytes()
    if MARKER not in text:
        raise SystemExit("science_render.py is missing the END-HEADER marker")
    body = text.split(MARKER, 1)[1]
    digest = hashlib.sha256(body).hexdigest()
    if digest != EXPECTED:
        raise SystemExit(f"science_render.py body {digest} is not f65657b9")
    print("science_render.py matches f65657b9")
    return 0


if __name__ == "__main__":
    sys.exit(main())
