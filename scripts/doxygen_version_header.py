#!/usr/bin/env python3
# Copyright (c) d-SEAMS developers
"""Write seams_version.h for Doxygen from the Meson project version."""

import re
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]


def main():
    text = (ROOT / "meson.build").read_text(encoding="utf-8")
    match = re.search(r"version:\s*'([0-9]+)\.([0-9]+)\.([0-9]+)'", text)
    if match is None:
        print("meson project version is missing", file=sys.stderr)
        return 1
    major, minor, patch = match.groups()
    version = "{}.{}.{}".format(major, minor, patch)
    dest = ROOT / "docs" / "source" / "doxygen-include"
    dest.mkdir(parents=True, exist_ok=True)
    number = (int(major) << 16) | (int(minor) << 8) | int(patch)
    body = (
        "#ifndef SEAMS_VERSION_H_\n"
        "#define SEAMS_VERSION_H_\n"
        "#define SEAMS_VERSION_MAJOR {major}\n"
        "#define SEAMS_VERSION_MINOR {minor}\n"
        "#define SEAMS_VERSION_PATCH {patch}\n"
        '#define SEAMS_VERSION_STRING "{version}"\n'
        "#define SEAMS_VERSION_NUMBER {number}\n"
        "#endif\n"
    ).format(major=major, minor=minor, patch=patch, version=version, number=number)
    (dest / "seams_version.h").write_text(body, encoding="utf-8")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
