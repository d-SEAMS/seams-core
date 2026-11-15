#!/usr/bin/env python3
"""Compile and run the C++ language-front Ice Ih census.

Usage: run_lang_ice.py ROOT BUILD
"""

import os
import subprocess
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
from toolchain import c_compiler, link_command  # noqa: E402


def main():
    if len(sys.argv) != 3:
        raise SystemExit("usage: run_lang_ice.py ROOT BUILD")
    root = Path(sys.argv[1])
    build = Path(sys.argv[2])
    cc = c_compiler()
    work = build / "docs-lang-ice"
    work.mkdir(parents=True, exist_ok=True)
    binary = work / "lang_ice"
    src_lib = str(build / "src")
    lang_lib = str(build / "bindings")
    subprocess.check_call(
        link_command(
            cc,
            [
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
            ],
        )
    )
    env = os.environ.copy()
    env["LD_LIBRARY_PATH"] = lang_lib + os.pathsep + src_lib + os.pathsep + env.get("LD_LIBRARY_PATH", "")
    subprocess.check_call([str(binary)], cwd=root, env=env)
    print("lang_ice ok")


if __name__ == "__main__":
    main()
