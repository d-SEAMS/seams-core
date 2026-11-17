#!/usr/bin/env python3
# Copyright (c) d-SEAMS developers
"""Compile and run the Ice Ih and sII examples in C, C++ and ctypes Python.

The Lua, Julia and R copies of this census live in the binding
repositories.

Usage: run_core_counts.py ROOT BUILD
"""

import os
import subprocess
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT / "docs"))
from toolchain import c_compiler, cxx_compiler, link_command  # noqa: E402


def library_path(build):
    parts = [str(build / "src")]
    current = os.environ.get("LD_LIBRARY_PATH")
    if current:
        parts.append(current)
    return os.pathsep.join(parts)


def run(cmd, env=None, cwd=None):
    subprocess.check_call(cmd, env=env, cwd=cwd or ROOT)


def compile_c_family(build, env):
    cc = c_compiler()
    cxx = cxx_compiler()
    work = build / "docs-core-counts"
    work.mkdir(parents=True, exist_ok=True)
    c_bin = work / "core_counts_c"
    run(
        link_command(
            cc,
            [
                "-I",
                str(ROOT / "include"),
                "-I",
                str(build / "src"),
                str(ROOT / "docs/examples/core_counts.c"),
                "-L" + str(build / "src"),
                "-lseams",
                "-Wl,-rpath," + str(build / "src"),
                "-o",
                str(c_bin),
            ],
        )
    )
    run([str(c_bin)], env=env)
    cpp_bin = work / "core_counts_cpp"
    run(
        link_command(
            cxx,
            [
                "-I",
                str(ROOT / "include"),
                "-I",
                str(build / "src"),
                str(ROOT / "docs/examples/core_counts.cpp"),
                "-L" + str(build / "src"),
                "-lseams",
                "-Wl,-rpath," + str(build / "src"),
                "-o",
                str(cpp_bin),
            ],
        )
    )
    run([str(cpp_bin)], env=env)


def run_python(env):
    env = env.copy()
    env["SEAMS_LIBDIR"] = str(Path(env["SEAMS_LIBDIR"]))
    run([sys.executable, "-m", "doctest", "-v", str(ROOT / "docs/examples/core_counts.py")], env=env)
    run([sys.executable, str(ROOT / "docs/examples/core_counts.py")], env=env)


def main():
    if len(sys.argv) != 3:
        raise SystemExit("usage: run_core_counts.py ROOT BUILD")
    build = Path(sys.argv[2])
    env = os.environ.copy()
    env["LD_LIBRARY_PATH"] = library_path(build)
    env["SEAMS_LIBDIR"] = str(build / "src")
    compile_c_family(build, env)
    run_python(env)
    print("core counts ok")


if __name__ == "__main__":
    main()
