#!/usr/bin/env python3
# Copyright (c) d-SEAMS developers
"""Build the six family examples out of tree against an installed libseams.

Usage: check_installed.py BUILD SOURCE

Installs BUILD into a fresh directory under BUILD/install-stages, then
builds every example in docs/examples/c twice, once with the flags
`pkg-config seams` reports and
once with CMake's find_package(seams), and runs each binary against the
installed library. CC names the C compiler.
"""

import importlib.util
import os
import shlex
import subprocess
import sys
from pathlib import Path


def installed_seams(stage):
    for name in ("libseams.so", "libseams.dylib", "seams.dll"):
        found = [path for path in stage.rglob(name) if path.is_file()]
        if found:
            return found
    return []

FAMILIES = ("neighbours", "rings", "chill_plus", "q4q6", "cages", "errors")
WARN = ["-std=c11", "-Wall", "-Wextra", "-Wpedantic", "-Werror"]


def helpers(source):
    path = source / "examples" / "c" / "check_installed.py"
    spec = importlib.util.spec_from_file_location("quickstart_installed", path)
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


def assert_loads_installed(label, family, binary, run_env, stage, libname):
    if sys.platform == "darwin":
        listed = subprocess.check_output(["otool", "-L", str(binary)], env=run_env, text=True)
        load = subprocess.check_output(["otool", "-l", str(binary)], env=run_env, text=True)
        named = [line for line in listed.splitlines() if libname in line]
        if len(named) != 1 or str(stage) not in load:
            raise SystemExit("{} {} does not load the installed libseams:\n{}".format(label, family, listed + load))
        return
    listed = subprocess.check_output(["ldd", str(binary)], env=run_env, text=True)
    resolved = [line for line in listed.splitlines() if line.strip().startswith(libname)]
    if len(resolved) != 1 or str(stage) not in resolved[0]:
        raise SystemExit("{} {} does not load the installed libseams:\n{}".format(label, family, listed))


def run_examples(label, binaries, env, source, search, help_mod, libname):
    stage = search[0]
    for family, binary in binaries:
        args = [str(binary)] if family == "errors" else [str(binary), str(source)]
        run_env = help_mod.library_dirs(binary, env, search)
        assert_loads_installed(label, family, binary, run_env, stage, libname)
        done = subprocess.run(args, env=run_env, capture_output=True, text=True)
        if done.returncode != 0:
            raise SystemExit("{} {} exit {}\n{}{}".format(label, family, done.returncode, done.stdout, done.stderr))
        print("{} {}: {}".format(label, family, done.stdout.strip()))


def main(argv):
    if len(argv) != 3:
        print(__doc__, file=sys.stderr)
        return 2
    build = Path(argv[1]).resolve()
    source = Path(argv[2]).resolve()
    help_mod = helpers(source)
    base = help_mod.scratch(build, "capi-family")
    stage = base / "destdir"
    work = base / "work"
    work.mkdir()
    subprocess.check_call(
        ["meson", "install", "-C", str(build), "--no-rebuild", "--quiet", "--destdir", str(stage)]
    )

    pc_files = list(stage.rglob("seams.pc"))
    configs = list(stage.rglob("seams-config.cmake"))
    libs = installed_seams(stage)
    if len(pc_files) != 1 or len(configs) != 1 or not libs:
        raise SystemExit("install did not produce one seams.pc, seams-config.cmake and libseams")
    env = os.environ.copy()
    env["PKG_CONFIG_PATH"] = str(pc_files[0].parent)
    env["PKG_CONFIG_SYSROOT_DIR"] = str(stage)
    # Only the runs see the staged libdir: ld would otherwise resolve the
    # staged libyodaLib and then want its dependencies, which are not staged.
    run_env = dict(env)
    run_env["LD_LIBRARY_PATH"] = os.pathsep.join([str(libs[0].parent), env.get("LD_LIBRARY_PATH", "")])
    if sys.platform == "darwin":
        run_env["DYLD_LIBRARY_PATH"] = run_env["LD_LIBRARY_PATH"]
        link_line = shlex.split(os.environ.get("SEAMS_OPENMP_LINK", ""))
        link_line = link_line + ["-Wl,-rpath," + str(libs[0].parent)]
    else:
        link_line = shlex.split(os.environ.get("SEAMS_OPENMP_LINK", ""))
    search = [stage, build]
    cc = shlex.split(os.environ.get("CC", "cc"))
    libname = libs[0].name

    pkg_config = help_mod.pkg_config_program()
    flags = subprocess.check_output([pkg_config, "--cflags", "--libs", "seams"], env=env, text=True).split()
    examples = source / "docs" / "examples" / "c"
    binaries = []
    for family in FAMILIES:
        binary = work / ("pkg-" + family)
        subprocess.check_call(
            [*cc, *WARN, str(examples / (family + ".c")), "-o", str(binary), *flags, *link_line],
            env=env,
        )
        binaries.append((family, binary))
    run_examples("pkg-config", binaries, run_env, source, search, help_mod, libname)

    installed = stage / Path(help_mod.prefix_of(build)).relative_to("/")
    cmake_build = work / "cmake"
    subprocess.check_call(
        [
            "cmake", "-S", str(examples), "-B", str(cmake_build),
            "-DCMAKE_C_COMPILER=" + cc[0],
            "-DCMAKE_PREFIX_PATH=" + str(installed),
            "-Dseams_DIR=" + str(configs[0].parent),
            "-DCMAKE_BUILD_TYPE=Release",
            "-DCMAKE_EXE_LINKER_FLAGS=" + " ".join(link_line),
        ],
        env=env,
    )
    subprocess.check_call(["cmake", "--build", str(cmake_build)], env=env)
    run_examples("cmake", [(family, cmake_build / family) for family in FAMILIES], run_env, source, search, help_mod, libname)
    print("installed family examples ok")
    return 0


if __name__ == "__main__":
    raise SystemExit(main(sys.argv))
