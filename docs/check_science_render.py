#!/usr/bin/env python3
"""The vendored renderer stays the f65657b9 body plus the shared palette.

The palette import and the triclinic cell read are the only edits.
docs/check_sh_cell.py and docs/check_palette.py lock those two edits.
"""

import hashlib
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
SOURCE = ROOT / "scripts" / "science_render.py"
MARKER = b"# END-HEADER\n"
# sha256 of the body after the palette import and the triclinic cell read.
# The renderer underneath is still f65657b9b89481e7dd6617b1987bc58b2ec0f01c.
EXPECTED = "95fb812110839cbc6310776f547356acab306b3f081bab13c93fc33139111064"


def main():
    text = SOURCE.read_bytes()
    if MARKER not in text:
        raise SystemExit("science_render.py is missing the END-HEADER marker")
    body = text.split(MARKER, 1)[1]
    digest = hashlib.sha256(body).hexdigest()
    if digest != EXPECTED:
        raise SystemExit(f"science_render.py body {digest} changed")
    print("science_render.py body matches the shared-palette copy")
    return 0


if __name__ == "__main__":
    sys.exit(main())
