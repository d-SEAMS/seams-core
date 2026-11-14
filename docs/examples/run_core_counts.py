#!/usr/bin/env python3
"""Compile and run the Ice Ih and sII examples.

Usage: run_core_counts.py ROOT BUILD
"""

import os
import shutil
import subprocess
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[2]


def library_path(build):
    parts = [str(build / "src"), str(build / "bindings")]
    current = os.environ.get("LD_LIBRARY_PATH")
    if current:
        parts.append(current)
    return os.pathsep.join(parts)


def run(cmd, env=None, cwd=None):
    subprocess.check_call(cmd, env=env, cwd=cwd or ROOT)


def compile_c_family(build, env):
    cc = shutil.which("cc") or shutil.which("gcc")
    cxx = shutil.which("c++") or shutil.which("g++")
    if cc is None or cxx is None:
        raise SystemExit("cc and c++ are required")
    work = build / "docs-core-counts"
    work.mkdir(parents=True, exist_ok=True)
    c_bin = work / "core_counts_c"
    run(
        [
            cc,
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
        ]
    )
    run([str(c_bin)], env=env)
    cpp_bin = work / "core_counts_cpp"
    run(
        [
            cxx,
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
        ]
    )
    run([str(cpp_bin)], env=env)


def run_python(env):
    env = env.copy()
    env["SEAMS_LIBDIR"] = str(Path(env["SEAMS_LIBDIR"]))
    run([sys.executable, "-m", "doctest", "-v", str(ROOT / "docs/examples/core_counts.py")], env=env)
    run([sys.executable, str(ROOT / "docs/examples/core_counts.py")], env=env)


def run_lua(build, env):
    lua = shutil.which("lua5.4") or shutil.which("lua")
    if lua is None:
        raise SystemExit("lua is not on PATH")
    modules = list((build / "bindings").rglob("dseams_core.so"))
    if not modules:
        raise SystemExit("dseams_core.so was not built")
    lua_env = env.copy()
    lua_env["LUA_PATH"] = str(ROOT / "bindings" / "lua") + "/?.lua;;"
    lua_env["LUA_CPATH"] = str(modules[0].parent) + "/?.so;;"
    run([lua, str(ROOT / "docs/examples/core_counts.lua")], env=lua_env)


def run_julia(build, env):
    sys.path.insert(0, str(ROOT / "docs/examples"))
    import run_binding

    julia = run_binding.julia_bin()
    julia_env = env.copy()
    julia_env["DSEAMS_ROOT"] = str(ROOT)
    julia_env["SEAMS_LANG_LIB"] = str(build / "bindings" / "libseams_lang.so")
    run(
        [
            julia,
            "--project=" + str(ROOT / "bindings/julia/Dseams"),
            str(ROOT / "docs/examples/core_counts.jl"),
        ],
        env=julia_env,
    )


def run_r(build, env):
    sys.path.insert(0, str(ROOT / "docs/examples"))
    import run_binding

    script, r_env = run_binding.install_r_package(build)
    r_env = {**env, **r_env}
    run([script, str(ROOT / "docs/examples/core_counts.R")], env=r_env)


def main():
    if len(sys.argv) != 3:
        raise SystemExit("usage: run_core_counts.py ROOT BUILD")
    build = Path(sys.argv[2])
    env = os.environ.copy()
    env["LD_LIBRARY_PATH"] = library_path(build)
    env["SEAMS_LIBDIR"] = str(build / "src")
    compile_c_family(build, env)
    run_python(env)
    run_lua(build, env)
    run_julia(build, env)
    run_r(build, env)
    print("core counts ok")


if __name__ == "__main__":
    main()
