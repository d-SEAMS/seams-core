#!/usr/bin/env python3
"""Compile and run the C++ language-front Ice Ih census.

Usage: run_lang_ice.py ROOT BUILD
"""

import os
import shutil
import subprocess
import sys
from pathlib import Path


def main():
    if len(sys.argv) != 3:
        raise SystemExit("usage: run_lang_ice.py ROOT BUILD")
    root = Path(sys.argv[1])
    build = Path(sys.argv[2])
    cc = shutil.which("cc") or shutil.which("gcc")
    if cc is None:
        raise SystemExit("cc is required")
    work = build / "docs-lang-ice"
    work.mkdir(parents=True, exist_ok=True)
    binary = work / "lang_ice"
    src_lib = str(build / "src")
    lang_lib = str(build / "bindings")
    subprocess.check_call(
        [
            cc,
            "-I",
            str(root / "bindings" / "cxx"),
            str(root / "docs" / "examples" / "lang_ice.c"),
            "-L" + lang_lib,
            "-L" + src_lib,
            "-lseams_lang",
            "-lseams",
            "-Wl,-rpath," + lang_lib,
            "-Wl,-rpath," + src_lib,
            "-o",
            str(binary),
        ]
    )
    env = os.environ.copy()
    env["LD_LIBRARY_PATH"] = lang_lib + os.pathsep + src_lib + os.pathsep + env.get("LD_LIBRARY_PATH", "")
    subprocess.check_call([str(binary)], cwd=root, env=env)
    print("lang_ice ok")


if __name__ == "__main__":
    main()
