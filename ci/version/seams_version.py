#!/usr/bin/env python3
"""The release version, read from meson.build.

meson.build is the one hand-written copy. cogapp regenerates every path
in ci/version/files.txt from VERSION. Call cogapp as ``python -m cogapp``.
The cogapp wheel also installs a ``cog`` script, and that script would
shadow cocogitto's ``cog``.
"""

from __future__ import annotations

import re
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[2]
MESON = ROOT / "meson.build"
_VERSION = re.compile(
    r"(project\(\s*'[^']+',\s*\[[^\]]+\],\s*version:\s*')([^']+)'",
    re.DOTALL,
)


def _read() -> str:
    match = _VERSION.search(MESON.read_text(encoding="utf-8"))
    if match is None:
        raise SystemExit("meson.build has no project version")
    return match.group(2)


VERSION = _read()


def set_version(version: str) -> None:
    if re.fullmatch(r"\d+\.\d+\.\d+(?:\.[0-9A-Za-z]+)*", version) is None:
        raise SystemExit(f"not a release version: {version}")
    text = MESON.read_text(encoding="utf-8")
    new, count = _VERSION.subn(lambda m: m.group(1) + version + "'", text, count=1)
    if count != 1:
        raise SystemExit("meson.build version was not updated")
    MESON.write_text(new, encoding="utf-8")


def main(argv: list[str]) -> int:
    if not argv:
        print(VERSION)
        return 0
    if len(argv) == 2 and argv[0] == "set":
        set_version(argv[1])
        print(argv[1])
        return 0
    print("usage: seams_version.py [set VERSION]", file=sys.stderr)
    return 2


if __name__ == "__main__":
    raise SystemExit(main(sys.argv[1:]))
