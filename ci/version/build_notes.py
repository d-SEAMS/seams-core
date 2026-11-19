#!/usr/bin/env python3
"""Render changelog.d into CHANGELOG.md during a cog bump.

Skip when that version's heading is already in CHANGELOG.md, so a second
bump of the same version does not duplicate the section. Skip when
changelog.d has no fragments.
"""

from __future__ import annotations

import subprocess
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[2]


def main() -> int:
    if len(sys.argv) != 2:
        print("usage: build_notes.py VERSION", file=sys.stderr)
        return 2
    version = sys.argv[1]
    changelog = (ROOT / "CHANGELOG.md").read_text(encoding="utf-8")
    if f"## [{version}]" in changelog:
        print(f"## [{version}] is already in CHANGELOG.md")
        return 0
    fragments = [
        path
        for path in (ROOT / "changelog.d").iterdir()
        if path.suffix == ".md" and path.name != "README.md"
    ]
    if not fragments:
        print("no unreleased fragments")
        return 0
    subprocess.check_call(
        [sys.executable, "-m", "towncrier", "build", "--yes", "--version", version],
        cwd=ROOT,
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
