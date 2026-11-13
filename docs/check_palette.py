#!/usr/bin/env python3
"""No colour in the shared palette may mean two roles."""

import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "scripts"))
from chill_palette import COLOUR, assert_distinct


def main():
    distance = assert_distinct()
    assert_distinct(COLOUR)
    print(f"palette {len(COLOUR)} roles, closest CIEDE2000 {distance:.2f}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
