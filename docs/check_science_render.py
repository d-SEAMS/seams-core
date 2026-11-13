#!/usr/bin/env python3
"""The vendored renderer stays the bb1dbdbf body plus the shared palette.

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
# The renderer underneath is still bb1dbdbf25553174a47da6fa3f4be8f31238d7f0.
EXPECTED = "122aff01ddd900e6f654675a728c5235157ffc45d733c0bbff9eac21f3dc8a70"


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
