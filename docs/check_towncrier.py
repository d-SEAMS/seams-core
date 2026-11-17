#!/usr/bin/env python3
# Copyright (c) d-SEAMS developers
"""Render the 2.14.0 towncrier draft and require it in CHANGELOG.md.

Fragments under changelog.d stay on disk. The draft is compared to the
changelog, and every fragment file is required to still be there.

Fragments for the next release sit in changelog.d/2.14.1. Towncrier does
not read that directory, so the 2.14.0 draft leaves them out, and each one
is required in the 2.14.1 section of CHANGELOG.md.
"""

import subprocess
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
CHANGELOG = ROOT / "CHANGELOG.md"
FRAGMENTS = ROOT / "changelog.d"
NEXT = FRAGMENTS / "2.14.1"


def draft():
    done = subprocess.run(
        [
            sys.executable,
            "-m",
            "towncrier",
            "build",
            "--draft",
            "--version",
            "2.14.0",
            "--date",
            "2026-11-15",
        ],
        cwd=str(ROOT),
        check=False,
        capture_output=True,
        text=True,
    )
    if done.returncode != 0:
        sys.stderr.write(done.stdout)
        sys.stderr.write(done.stderr)
        return None
    text = done.stdout
    marker = text.find("## [2.14.0]")
    if marker < 0:
        sys.stderr.write(text)
        return None
    return text[marker:].strip() + "\n"


def section(text, version):
    start = text.find("## [{}]".format(version))
    if start < 0:
        return ""
    end = text.find("\n## [", start + 1)
    return text[start:] if end < 0 else text[start:end]


def missing_next(changelog):
    if not NEXT.is_dir():
        return []
    notes = " ".join(section(changelog, "2.14.1").split())
    return [
        path.name
        for path in sorted(NEXT.glob("*.md"))
        if " ".join(path.read_text().split()) not in notes
    ]


def main():
    rendered = draft()
    if not rendered:
        print("towncrier draft failed", file=sys.stderr)
        return 1
    changelog = CHANGELOG.read_text()
    if rendered.strip() not in changelog:
        print("CHANGELOG.md does not contain the towncrier 2.14.0 draft", file=sys.stderr)
        return 1
    absent = missing_next(changelog)
    if absent:
        print("CHANGELOG.md 2.14.1 lacks " + ", ".join(absent), file=sys.stderr)
        return 1
    names = [
        path.name
        for path in FRAGMENTS.iterdir()
        if path.is_file() and path.suffix == ".md" and path.name not in {"README.md"}
    ]
    if len(names) < 100:
        print("changelog.d lost fragments ({})".format(len(names)), file=sys.stderr)
        return 1
    print("towncrier 2.14.0 renders, {} fragments kept".format(len(names)))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
