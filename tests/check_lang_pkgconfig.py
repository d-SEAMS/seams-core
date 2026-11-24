#!/usr/bin/env python3
# Copyright (c) d-SEAMS developers
"""Build a seams-lang consumer from pkg-config alone, in and out of the tree.

Usage: check_lang_pkgconfig.py BUILD SOURCE

Builds tests/lang_pkgconfig_consumer.c twice with the flags
`pkg-config --cflags --libs seams-lang` reports: once with PKG_CONFIG_PATH
at meson-uninstalled in the top build directory (BUILD, or the parent
project's when seams-core is a subproject), the build-tree flow a binding
uses against a seams-core checkout, and once against a fresh install under
BUILD/install-stages. Each binary records the -L directories as its
runpath and runs with LD_LIBRARY_PATH unset and LD_BIND_NOW=1, so a
symbol the loaded libraries do not define stops it at startup. ldd must
show libseams_lang, libyodaLib and the vendored libseams-chemfiles from
the tree under test, and no unpatched libchemfiles. CC names the C
compiler and SEAMS_OPENMP_LINK the consumer link line of this build.
"""

import importlib.util
import os
import shlex
import subprocess
import sys
from pathlib import Path

WARN = ["-std=c11", "-Wall", "-Wextra", "-Wpedantic", "-Werror"]
CHEMFILES = "libseams-chemfiles.so.0.10"


def helpers(source):
    path = source / "examples" / "c" / "check_installed.py"
    spec = importlib.util.spec_from_file_location("quickstart_installed", path)
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


def pkg_config(program, env, *args):
    return subprocess.check_output([program, *args], env=env, text=True).strip()


def loaded(binary, env):
    """ldd of BINARY as {soname: path}, failing on any library not found."""
    listed = subprocess.check_output(["ldd", str(binary)], env=env, text=True)
    found = {}
    for line in listed.splitlines():
        if "=>" not in line:
            continue
        name, _, rest = line.strip().partition(" => ")
        if rest.startswith("not found"):
            raise SystemExit("{} does not find {}:\n{}".format(binary, name, listed))
        found[name] = rest.split(" (")[0]
    return found, listed


def check(label, pc_env, define, root, source, work, version, vendored, help_mod):
    program = help_mod.pkg_config_program()
    extra = [define] if define else []
    modversion = pkg_config(program, pc_env, *extra, "--modversion", "seams-lang")
    if modversion != version:
        raise SystemExit("{}: pkg-config says seams-lang {}, project is {}".format(label, modversion, version))
    flags = shlex.split(pkg_config(program, pc_env, *extra, "--cflags", "--libs", "seams-lang"))
    libdirs = [flag[2:] for flag in flags if flag.startswith("-L")]
    for libdir in libdirs:
        if not Path(libdir).is_dir():
            raise SystemExit("{}: pkg-config names -L{}, which does not exist".format(label, libdir))
    cc = shlex.split(os.environ.get("CC", "").strip())
    if not cc:
        raise SystemExit("CC is unset. The test environment must name the compiler meson used.")
    link_line = shlex.split(os.environ.get("SEAMS_OPENMP_LINK", ""))
    binary = work / ("consumer-" + label)
    subprocess.check_call(
        [*cc, *WARN, str(source / "tests" / "lang_pkgconfig_consumer.c"), "-o", str(binary),
         *flags, *link_line, *["-Wl,-rpath," + libdir for libdir in libdirs]],
        env=pc_env,
    )

    run_env = dict(os.environ)
    run_env.pop("LD_LIBRARY_PATH", None)
    run_env["LD_BIND_NOW"] = "1"
    libs, listed = loaded(binary, run_env)
    wanted = ["libseams_lang.so", "libyodaLib.so"] + ([CHEMFILES] if vendored else [])
    for name in wanted:
        path = libs.get(name)
        if path is None or not Path(path).resolve().is_relative_to(root.resolve()):
            raise SystemExit("{}: {} does not load {} from {}:\n{}".format(label, binary.name, name, root, listed))
    unpatched = [name for name in libs if name.startswith("libchemfiles.so")]
    if vendored and unpatched:
        raise SystemExit("{}: {} loads the unpatched {}:\n{}".format(label, binary.name, unpatched[0], listed))

    sample = source / "input" / "xyz" / "configuration.xyz"
    done = subprocess.run([str(binary), version, str(sample)], env=run_env, capture_output=True, text=True)
    expected = "seams-lang {} frames 1 atoms 256".format(version)
    if done.returncode != 0 or done.stdout.strip() != expected:
        raise SystemExit("{}: {} exit {}, expected {!r}\n{}{}".format(
            label, binary.name, done.returncode, expected, done.stdout, done.stderr))
    print("{}: {}".format(label, done.stdout.strip()))


def main(argv):
    if len(argv) != 3:
        print(__doc__, file=sys.stderr)
        return 2
    build = Path(argv[1]).resolve()
    source = Path(argv[2]).resolve()
    help_mod = helpers(source)
    version = help_mod.project_version(build)
    vendored = (build / "chemfiles-prefix" / "lib" / CHEMFILES).exists()
    base = help_mod.scratch(build, "lang-pkgconfig")
    work = base / "work"
    work.mkdir()

    env = os.environ.copy()
    env.pop("PKG_CONFIG_SYSROOT_DIR", None)
    previous = env.get("PKG_CONFIG_PATH", "")
    uninstalled = str(help_mod.meson_root(build) / "meson-uninstalled")
    tree_env = dict(env)
    tree_env["PKG_CONFIG_PATH"] = uninstalled + (os.pathsep + previous if previous else "")
    check("uninstalled", tree_env, None, build, source, work, version, vendored, help_mod)

    stage = base / "destdir"
    subprocess.check_call(
        ["meson", "install", "-C", str(help_mod.meson_root(build)), "--no-rebuild", "--quiet", "--destdir", str(stage)]
    )
    pc_files = list(stage.rglob("seams-lang.pc"))
    if len(pc_files) != 1:
        raise SystemExit("install did not produce one seams-lang.pc")
    installed = stage / Path(help_mod.prefix_of(build)).relative_to("/")
    stage_env = dict(env)
    stage_env["PKG_CONFIG_PATH"] = str(pc_files[0].parent) + (os.pathsep + previous if previous else "")
    stage_env["PKG_CONFIG_DISABLE_UNINSTALLED"] = "1"
    # A multiarch libdir makes pkg-config relocate prefix; name the staged one.
    define = "--define-variable=prefix=" + str(installed)
    check("installed", stage_env, define, stage, source, work, version, vendored, help_mod)
    print("seams-lang pkg-config consumers ok")
    return 0


if __name__ == "__main__":
    raise SystemExit(main(sys.argv))
