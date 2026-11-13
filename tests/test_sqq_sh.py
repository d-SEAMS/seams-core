#!/usr/bin/env python3
"""SQQ rejects the triclinic sH dump and a triclinic box mode.

Exit 77 when SQQ is not installed. The orthorhombic sI dump is accepted
by the same reader, so the rejection is the tilt.
"""

from __future__ import annotations

import runpy
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]


def main() -> None:
    try:
        import sqq  # noqa: F401
    except ImportError:
        print("sqq is not installed", file=sys.stderr)
        raise SystemExit(77)
    sys.argv = [str(ROOT / "scripts" / "h2h_sqq_box.py")]
    runpy.run_path(str(ROOT / "scripts" / "h2h_sqq_box.py"), run_name="__main__")


if __name__ == "__main__":
    main()
