#!/usr/bin/env python3
# Copyright (c) d-SEAMS developers
"""Compile c_api_tour.c, run it, and require docs/api/c.md to quote every line."""

import os
import re
import subprocess
import sys
from pathlib import Path


def fail(message):
    print(message, file=sys.stderr)
    raise SystemExit(1)


def main():
    root = Path(sys.argv[1])
    libdir = Path(sys.argv[2])
    build = Path(sys.argv[3])
    page = (root / "docs" / "api" / "c.md").read_text()
    source = root / "docs" / "examples" / "c_api_tour.c"
    header = (root / "include" / "seams.h").read_text()
    functions = re.findall(r"SEAMS_API\s+(?:const\s+char\s+\*|[\w\s\*]+?)\s*(seams_[A-Za-z0-9_]+)\s*\(", header)
    binary = build / "c_api_tour"
    binary.parent.mkdir(parents=True, exist_ok=True)
    compile_cmd = [
        "gcc",
        "-I" + str(root / "include"),
        "-I" + str(libdir.parent / "src") if (libdir.parent / "src" / "seams_version.h").is_file() else "-I" + str(libdir),
        str(source),
        "-L" + str(libdir),
        "-lseams",
        "-lm",
        "-Wl,-rpath," + str(libdir),
        "-o",
        str(binary),
    ]
    # seams_version.h is generated next to the library sources.
    version = None
    for folder in (libdir, libdir.parent / "src", Path(sys.argv[2])):
        candidate = Path(folder) / "seams_version.h"
        if candidate.is_file():
            version = candidate.parent
            break
    if version is None:
        fail("seams_version.h was not found beside the library")
    compile_cmd = [
        "gcc",
        "-I" + str(root / "include"),
        "-I" + str(version),
        str(source),
        "-L" + str(libdir),
        "-lseams",
        "-lm",
        "-Wl,-rpath," + str(libdir),
        "-o",
        str(binary),
    ]
    done = subprocess.run(compile_cmd, text=True, capture_output=True)
    if done.returncode != 0:
        fail(done.stderr[-800:])
    env = os.environ.copy()
    env["LD_LIBRARY_PATH"] = str(libdir) + os.pathsep + env.get("LD_LIBRARY_PATH", "")
    env["OMP_NUM_THREADS"] = "1"
    run = subprocess.run([str(binary), str(root)], text=True, capture_output=True, env=env)
    if run.returncode != 0:
        fail(run.stderr[-800:])
    lines = [row for row in run.stdout.splitlines() if row.startswith("capi ")]
    if len(lines) < 60:
        fail("tour printed {} lines".format(len(lines)))
    missing = [row for row in lines if row not in page]
    if missing:
        fail("docs/api/c.md is missing:\n" + "\n".join(missing))
    documented = re.findall(r"(?m)^capi \S.*$", page)
    extra = [row for row in documented if row not in lines]
    if extra:
        fail("docs/api/c.md quotes lines the tour did not print:\n" + "\n".join(extra[:20]))
    print("c api documented lines {}".format(len(documented)))
    quiet = [name for name in functions if "capi " + name not in run.stdout]
    if quiet:
        fail("tour did not print " + " ".join(quiet))
    print("c api tour {} lines".format(len(lines)))


if __name__ == "__main__":
    main()
