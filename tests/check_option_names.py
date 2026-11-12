#!/usr/bin/env python3
"""Parent option names must not be copied onto a subproject, and the CLI
must not list the same long option twice.

A wrap patch that ships meson.build and no options file makes
`meson configure` load the parent meson_options.txt under that
subproject. The copies keep the file defaults, so one name is auto in
the parent and disabled in the subproject.
"""

import re
import subprocess
import sys
from pathlib import Path


def patch_dirs(root: Path):
    wraps = root / "subprojects"
    found = []
    for wrap in sorted(wraps.glob("*.wrap")):
        text = wrap.read_text()
        match = re.search(r"(?m)^patch_directory\s*=\s*(\S+)\s*$", text)
        if not match:
            continue
        patch = wraps / "packagefiles" / match.group(1)
        if (patch / "meson.build").is_file():
            found.append(patch)
    return found


def options_file(patch: Path):
    for name in ("meson.options", "meson_options.txt"):
        path = patch / name
        if path.is_file():
            return path
    return None


def parent_option_names(root: Path):
    text = (root / "meson_options.txt").read_text()
    return set(re.findall(r"option\(\s*'([^']+)'", text))


def long_options(help_text: str):
    names = re.findall(r"(?m)^\s+(--[a-z0-9-]+)", help_text)
    if not names:
        names = re.findall(r"--[a-z0-9-]+", help_text.split("positional arguments:", 1)[0])
    return names


def main():
    root = Path(sys.argv[1]).resolve()
    seams = Path(sys.argv[2])
    failed = False
    parent = parent_option_names(root)
    for patch in patch_dirs(root):
        path = options_file(patch)
        if path is None:
            print(f"missing options file in {patch}", file=sys.stderr)
            failed = True
            continue
        text = path.read_text()
        declared = set(re.findall(r"option\(\s*'([^']+)'", text))
        copied = sorted(declared & parent)
        if copied:
            print(f"{path} repeats parent options: {', '.join(copied)}", file=sys.stderr)
            failed = True
    help_run = subprocess.run([str(seams), "--help"], check=False, capture_output=True, text=True)
    if help_run.returncode != 0:
        print(help_run.stderr, file=sys.stderr)
        print("seams --help failed", file=sys.stderr)
        return 1
    names = long_options(help_run.stdout)
    seen = {}
    for name in names:
        if name in seen:
            print(f"repeated CLI option {name}", file=sys.stderr)
            failed = True
        seen[name] = True
    return 1 if failed else 0


if __name__ == "__main__":
    sys.exit(main())
