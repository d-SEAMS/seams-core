#!/usr/bin/env python3
# Copyright (c) d-SEAMS developers
"""The vendored renderer stays the shared-palette body.

The palette import, the triclinic cell read, and the American spelling
of that body are the edits. docs/check_sh_cell.py and
docs/check_palette.py lock the first two.
"""

import hashlib
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
SOURCE = ROOT / "scripts" / "science_render.py"
MARKER = b"# END-HEADER\n"
# sha256 of the body after the palette import, the triclinic cell read,
# and the American spelling of that body.
EXPECTED = "3baea759676b85fc7784c1126fac6910914830b504c8d949f3c82d73fc2f7738"


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
