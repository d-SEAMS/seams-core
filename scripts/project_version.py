#!/usr/bin/env python3
"""Read the project version from the Meson project declaration."""

from __future__ import annotations

import re
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]


def project_version() -> str:
    text = (ROOT / "meson.build").read_text(encoding="utf-8")
    match = re.search(r"^\s*version:\s*'([^']+)'", text, re.MULTILINE)
    if match is None:
        raise RuntimeError("meson.build has no project version")
    return match.group(1)


if __name__ == "__main__":
    print(project_version())
