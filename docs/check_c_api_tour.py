#!/usr/bin/env python3
# Copyright (c) d-SEAMS developers
"""Compile c_api_tour.c, run it, and require docs/api/c.md to quote every line."""

import os
import re
import subprocess
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))
from toolchain import c_compiler, link_command  # noqa: E402


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
    # seams_version.h is generated next to the library sources.
    version = None
    for folder in (libdir, libdir.parent / "src", Path(sys.argv[2])):
        candidate = Path(folder) / "seams_version.h"
        if candidate.is_file():
            version = candidate.parent
            break
    if version is None:
        fail("seams_version.h was not found beside the library")
    flags = [
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
    runtime = os.environ.get("TSAN_RUNTIME", "")
    if runtime and os.path.isfile(runtime):
        flags[0:0] = [
            "-fsanitize=thread",
            "-shared-libsan",
            "-Wl,-rpath," + os.path.dirname(runtime),
        ]
    compile_cmd = link_command(c_compiler(), flags)
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
