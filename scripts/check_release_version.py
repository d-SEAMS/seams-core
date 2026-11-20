#!/usr/bin/env python3
# Copyright (c) d-SEAMS developers
"""Require one release version across the build, packaging and citation files.

meson.build is the source of the version. pixi.toml, the conda recipe
and CITATION.cff must carry the same string, and the newest CHANGELOG.md
heading must be that version with the CITATION.cff ``date-released``.
``--tag vX.Y.Z`` also requires the tag to name that version, which is the
guard the release workflow runs before it builds anything. ``--print``
writes the version to stdout.

>>> first(r"^project\\([^)]*?version:\\s*'([^']+)'", "# c\\nproject('x',\\n  version: '1.2.3',\\n)", "m")
'1.2.3'
"""

from __future__ import annotations

import argparse
import json
import re
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]


def first(pattern, text, name):
    match = re.search(pattern, text, re.M)
    if match is None:
        raise SystemExit("{}: no match for {}".format(name, pattern))
    return match.group(1)


def read(name):
    return (ROOT / name).read_text(encoding="utf-8")


def codemeta_version(text):
    try:
        data = json.loads(text)
    except json.JSONDecodeError as exc:
        raise SystemExit("codemeta.json: invalid JSON ({})".format(exc))
    version = data.get("version")
    if not isinstance(version, str) or not version:
        raise SystemExit("codemeta.json: missing version")
    return version


def versions():
    found = {
        "meson.build": first(r"^project\([^)]*?version:\s*'([^']+)'", read("meson.build"), "meson.build"),
        "pixi.toml": first(r'^\[workspace\][^\[]*?^version\s*=\s*"([^"]+)"', read("pixi.toml"), "pixi.toml"),
        "recipe/recipe.yaml": first(r'^context:\s*\n\s+version:\s*"([^"]+)"', read("recipe/recipe.yaml"), "recipe/recipe.yaml"),
        "CITATION.cff": first(r'^version:\s*"?([^"\s]+)"?', read("CITATION.cff"), "CITATION.cff"),
        "codemeta.json": codemeta_version(read("codemeta.json")),
        "skills/d-seams/SKILL.md": first(
            r'^metadata:\s*\n(?:[ \t]+[^\n]*\n)*?[ \t]+version:\s*"?([^"\s]+)"?',
            read("skills/d-seams/SKILL.md"),
            "skills/d-seams/SKILL.md",
        ),
        "packaging/yggdrasil/L/libseams/build_tarballs.jl": first(
            r'^version\s*=\s*v"([^"]+)"',
            read("packaging/yggdrasil/L/libseams/build_tarballs.jl"),
            "packaging/yggdrasil/L/libseams/build_tarballs.jl",
        ),
    }
    return found


def main(argv):
    parser = argparse.ArgumentParser(description=__doc__.split("\n\n")[0])
    parser.add_argument("--tag", help="release tag, such as vX.Y.Z")
    parser.add_argument("--print", action="store_true", help="print the version")
    args = parser.parse_args(argv)

    found = versions()
    version = found["meson.build"]
    errors = [
        "{} has {}, meson.build has {}".format(name, value, version)
        for name, value in found.items()
        if value != version
    ]

    heading = re.search(r"^## \[([^\]]+)\] - (\d{4}-\d{2}-\d{2})$", read("CHANGELOG.md"), re.M)
    released = first(r'^date-released:\s*"?(\d{4}-\d{2}-\d{2})"?', read("CITATION.cff"), "CITATION.cff")
    if heading is None:
        errors.append("CHANGELOG.md has no '## [X.Y.Z] - YYYY-MM-DD' heading")
    else:
        if heading.group(1) != version:
            errors.append("the newest CHANGELOG.md section is {}, not {}".format(heading.group(1), version))
        if heading.group(2) != released:
            errors.append(
                "CHANGELOG.md dates {} {}, CITATION.cff date-released is {}".format(
                    heading.group(1), heading.group(2), released
                )
            )

    if args.tag is not None and args.tag != "v" + version:
        errors.append("tag {} does not name version {}".format(args.tag, version))

    if errors:
        print("\n".join(errors), file=sys.stderr)
        return 1
    if args.print:
        print(version)
    else:
        print("release version {} agrees in {} files".format(version, len(found) + 1))
    return 0


if __name__ == "__main__":
    raise SystemExit(main(sys.argv[1:]))
