#!/usr/bin/env python3
# Copyright (c) d-SEAMS developers
"""Render towncrier output from the configured fragment directory."""

import datetime as dt
import re
import subprocess
import sys
import tomllib
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
CHANGELOG = ROOT / "CHANGELOG.md"
CONFIG = ROOT / "towncrier.toml"
RELEASE = re.compile(r"^## \[(\d+\.\d+\.\d+)\] - (\d{4}-\d{2}-\d{2})$", re.M)

def read_config():
    data = tomllib.loads(CONFIG.read_text())
    tool = data.get("tool", {}).get("towncrier", {})
    directory = tool.get("directory")
    if not isinstance(directory, str) or not directory:
        raise ValueError("tool.towncrier.directory is missing")
    return ROOT / directory


def release_map():
    releases = {}
    for version, date in RELEASE.findall(CHANGELOG.read_text()):
        releases[version] = date
    return releases


def draft(version: str, date: str):
    done = subprocess.run(
        [
            sys.executable,
            "-m",
            "towncrier",
            "build",
            "--draft",
            "--version",
            version,
            "--date",
            date,
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
    marker = text.find("## [{}]".format(version))
    if marker < 0:
        sys.stderr.write(text)
        return None
    return text[marker:].strip() + "\n"


def main():
    try:
        fragments = read_config()
    except ValueError as exc:
        print(str(exc), file=sys.stderr)
        return 1
    if not fragments.is_dir():
        print("missing fragment directory: {}".format(fragments), file=sys.stderr)
        return 1
    releases = release_map()
    name = fragments.name
    if re.fullmatch(r"\d+\.\d+\.\d+", name):
        version = name
        date = releases.get(version, dt.date.today().isoformat())
    else:
        if not releases:
            print("CHANGELOG.md has no release headings", file=sys.stderr)
            return 1
        version, date = next(iter(releases.items()))
    rendered = draft(version, date)
    if not rendered:
        print("towncrier draft failed", file=sys.stderr)
        return 1
    changelog = CHANGELOG.read_text()
    present = rendered.strip() in changelog
    if version in releases and name == "changelog.d":
        if not present:
            print("CHANGELOG.md does not contain the towncrier {} draft".format(version),
                  file=sys.stderr)
            return 1
        print("towncrier {} draft matches CHANGELOG.md".format(version))
        return 0
    if present:
        print("CHANGELOG.md already contains the towncrier {} draft".format(version),
              file=sys.stderr)
        return 1
    print("towncrier {} renders from {}".format(version, fragments.relative_to(ROOT)))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
