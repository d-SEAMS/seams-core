#!/usr/bin/env python3
# Copyright (c) d-SEAMS developers
"""No color in the shared palette may mean two roles."""

import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "scripts"))
from chill_palette import COLOR, assert_distinct


def main():
    distance = assert_distinct()
    assert_distinct(COLOR)
    print(f"palette {len(COLOR)} roles, closest CIEDE2000 {distance:.2f}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
