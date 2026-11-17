#!/usr/bin/env python3
# Copyright (c) d-SEAMS developers
"""Compile the library sources the language bindings build with extra warnings.

Usage: check_binding_warnings.py BUILDDIR [FLAG ...]

PydSEAMSlib and dseams-r compile these sources under their own, stricter
flags. Each C++ source of the targets below is compiled again with the
command Meson recorded for it, FLAG appended, and -fsyntax-only. Subproject
sources are never compiled, so third-party code keeps its own warning set.
Any diagnostic fails the check. Without FLAG the default is
-Wsign-conversion -Wshorten-64-to-32 -Wunneeded-internal-declaration
-Wunused-function -Werror.
"""

import json
import os
import subprocess
import sys
from concurrent.futures import ThreadPoolExecutor
from pathlib import Path

TARGETS = ("yodaLib", "seams", "seams_lang", "shell_avx", "ylm_avx512")
DEFAULT_FLAGS = [
    "-Wsign-conversion",
    "-Wshorten-64-to-32",
    "-Wunneeded-internal-declaration",
    "-Wunused-function",
    "-Werror",
]


def commands(build):
    intro = json.loads((build / "meson-info" / "intro-targets.json").read_text())
    found = set()
    for target in intro:
        if target["name"] not in TARGETS or target["subproject"]:
            continue
        found.add(target["name"])
        for group in target["target_sources"]:
            if group.get("language") != "cpp":
                continue
            for source in group["sources"]:
                yield target["name"], source, group["compiler"] + group["parameters"]
    missing = set(TARGETS) - found - {"shell_avx", "ylm_avx512"}
    if missing:
        sys.exit(f"targets not configured in {build}: {', '.join(sorted(missing))}")


def check(build, item, flags):
    name, source, argv = item
    argv = [a for a in argv if a != "-fdiagnostics-color=always"]
    result = subprocess.run(
        argv + flags + ["-fno-diagnostics-color", "-fsyntax-only", source],
        cwd=build,
        capture_output=True,
        text=True,
    )
    output = (result.stdout + result.stderr).strip()
    return name, source, result.returncode, output


def main():
    if len(sys.argv) < 2:
        sys.exit(__doc__)
    build = Path(sys.argv[1]).resolve()
    flags = sys.argv[2:] or DEFAULT_FLAGS
    items = list(commands(build))
    failed = 0
    with ThreadPoolExecutor(max_workers=os.cpu_count() or 2) as pool:
        for name, source, code, output in pool.map(lambda i: check(build, i, flags), items):
            if code != 0 or output:
                failed += 1
                print(f"== {name}: {source}\n{output}\n", flush=True)
    print(f"{len(items) - failed} of {len(items)} sources clean under {' '.join(flags)}")
    return 1 if failed else 0


if __name__ == "__main__":
    sys.exit(main())
