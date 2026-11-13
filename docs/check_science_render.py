#!/usr/bin/env python3
"""The vendored renderer must match science/interface-nucleation bb1dbdbf."""

import hashlib
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
SOURCE = ROOT / "scripts" / "science_render.py"
MARKER = b"# END-HEADER\n"
# sha256 of scripts/science_render.py at bb1dbdbf25553174a47da6fa3f4be8f31238d7f0.
EXPECTED = "0d2c27b5d8db1341e932ab62618b73197a013b4b29c9524ccc6cb8e09971e6bf"


def main():
    text = SOURCE.read_bytes()
    if MARKER not in text:
        raise SystemExit("science_render.py is missing the END-HEADER marker")
    body = text.split(MARKER, 1)[1]
    digest = hashlib.sha256(body).hexdigest()
    if digest != EXPECTED:
        raise SystemExit(f"science_render.py body {digest} is not bb1dbdbf")
    print("science_render.py matches bb1dbdbf")
    return 0


if __name__ == "__main__":
    sys.exit(main())
