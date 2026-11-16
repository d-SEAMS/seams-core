#!/usr/bin/env python3
# Copyright (c) d-SEAMS developers
"""Install libseams and build examples/c against that install."""

import atexit
import os
import shutil
import subprocess
import sys
import tempfile
from pathlib import Path


def installed_seams(stage):
    """The installed libseams. .so is tried first, so Linux finds that file."""
    for name in ("libseams.so", "libseams.dylib", "seams.dll"):
        found = [path for path in stage.rglob(name) if path.is_file()]
        if found:
            return found
    return []


def pkg_config_program():
    """The conda pkg-config wrapper applies PKG_CONFIG_SYSROOT_DIR twice."""
    found = shutil.which("pkg-config")
    if found:
        sibling = Path(found).with_name("pkg-config.bin")
        if sibling.is_file():
            return str(sibling)
    return found or "pkg-config"


def run(cmd, env=None, cwd=None):
    subprocess.check_call(cmd, env=env, cwd=cwd)


def prefix_of(build):
    import json

    options = json.loads(
        subprocess.check_output(
            ["meson", "introspect", "--buildoptions", str(build)], text=True
        )
    )
    for option in options:
        if option["name"] == "prefix":
            return option["value"]
    raise SystemExit("meson prefix is missing")


def project_version(build):
    import json

    info = json.loads(
        subprocess.check_output(["meson", "introspect", "--projectinfo", str(build)], text=True)
    )
    return info["version"]


def scratch(build, prefix):
    """A fresh directory under BUILD that is removed at exit. A --repeat copy
    of the same test runs beside this one and must not share its install."""
    parent = build / "install-stages"
    parent.mkdir(exist_ok=True)
    path = Path(tempfile.mkdtemp(prefix=prefix + "-", dir=parent))
    atexit.register(shutil.rmtree, path, True)
    return path



def darwin_library_dirs(binary, env, search):
    """Same as library_dirs, using otool. dyld reads DYLD_LIBRARY_PATH."""
    listed = subprocess.check_output(["otool", "-L", str(binary)], text=True, env=env)
    missing = []
    for line in listed.splitlines()[1:]:
        text = line.strip()
        if not text:
            continue
        name = text.split(" (", 1)[0].strip()
        base = os.path.basename(name)
        if name.startswith("@rpath/") or name.startswith("@loader_path/") or not os.path.isfile(name):
            missing.append(base)
    found = dict(env)
    extra = []
    for base in missing:
        for root in search:
            for candidate in Path(root).rglob(base):
                if candidate.is_file():
                    extra.append(str(candidate.parent))
                    break
    path = found.get("DYLD_LIBRARY_PATH", found.get("LD_LIBRARY_PATH", ""))
    if extra or path:
        found["DYLD_LIBRARY_PATH"] = os.pathsep.join(extra + ([path] if path else []))
        found["LD_LIBRARY_PATH"] = found["DYLD_LIBRARY_PATH"]
    return found


def library_dirs(binary, env, search):
    """Add directories that hold shared libraries ldd cannot yet see."""
    if sys.platform == "darwin":
        return darwin_library_dirs(binary, env, search)
    missing = []
    listed = subprocess.check_output(["ldd", str(binary)], text=True, env=env)
    for line in listed.splitlines():
        if "not found" not in line:
            continue
        missing.append(line.strip().split()[0])
    if not missing:
        return env
    found = dict(env)
    path = found.get("LD_LIBRARY_PATH", "")
    extra = []
    for name in missing:
        for root in search:
            for candidate in root.rglob(name):
                if candidate.is_file():
                    extra.append(str(candidate.parent))
                    break
    if extra:
        found["LD_LIBRARY_PATH"] = os.pathsep.join(extra + ([path] if path else []))
    return found


def main():
    build = Path(sys.argv[1]).resolve()
    source = Path(sys.argv[2]).resolve()
    base = scratch(build, "capi")
    stage = base / "destdir"
    work = base / "work"
    work.mkdir()
    run(["meson", "install", "-C", str(build), "--no-rebuild", "--destdir", str(stage)])

    pc_files = list(stage.rglob("seams.pc"))
    configs = list(stage.rglob("seams-config.cmake"))
    if len(pc_files) != 1 or len(configs) != 1:
        raise SystemExit("install did not produce seams.pc and seams-config.cmake")
    pc = pc_files[0]
    libs_on_disk = installed_seams(stage)
    if not libs_on_disk:
        raise SystemExit("install did not produce libseams")
    libdir = libs_on_disk[0].parent
    prefix = prefix_of(build)
    installed = stage / Path(prefix).relative_to("/")
    env = os.environ.copy()
    env["PKG_CONFIG_PATH"] = str(pc.parent)
    # A multiarch libdir makes pkg-config relocate prefix to .../lib,
    # and PKG_CONFIG_SYSROOT_DIR then prepends the destdir a second time.
    # Name the destdir prefix explicitly and do not set a sysroot.
    env.pop("PKG_CONFIG_SYSROOT_DIR", None)
    pkg_config = pkg_config_program()
    define = "--define-variable=prefix=" + str(installed)
    cflags = subprocess.check_output(
        [pkg_config, define, "--cflags", "seams"], env=env, text=True
    ).split()
    libs = subprocess.check_output(
        [pkg_config, define, "--libs", "seams"], env=env, text=True
    ).split()
    binary = work / "quickstart-pkg"
    cc = os.environ.get("CC", "").strip()
    if not cc:
        raise SystemExit("CC is unset. The test environment must name the compiler meson used.")
    omp = os.environ.get("SEAMS_OPENMP_LINK", "").split()
    run(
        [cc, str(source / "examples" / "c" / "quickstart.c"), "-o", str(binary), *cflags, *libs, "-lm", *omp],
        env=env,
    )
    compile_env = dict(env)
    env["LD_LIBRARY_PATH"] = os.pathsep.join(
        [str(libdir), str(build / "src"), env.get("LD_LIBRARY_PATH", "")]
    )
    if sys.platform == "darwin":
        env["DYLD_LIBRARY_PATH"] = env["LD_LIBRARY_PATH"]
    env = library_dirs(binary, env, [stage, build])
    sample = source / "input" / "traj" / "agi_water_tiny.lammpstrj"
    run([str(binary), str(sample)], env=env)

    version = project_version(build)
    modversion = subprocess.check_output(
        [pkg_config, "--modversion", "seams"], env=env, text=True
    ).strip()
    if modversion != version:
        raise SystemExit("pkg-config says {}, project is {}".format(modversion, version))
    pinned = subprocess.run(
        [pkg_config, "--exists", "seams >= " + version], env=env
    )
    if pinned.returncode != 0:
        raise SystemExit("pkg-config --exists 'seams >= {}' failed".format(version))
    find_check = source / "examples" / "c" / "find_check"
    version_bin = work / "version-pkg"
    run(
        [cc, str(find_check / "version.c"), "-o", str(version_bin), *cflags, *libs, *omp],
        env=compile_env,
    )
    env = library_dirs(version_bin, env, [stage, build])
    run([str(version_bin), version], env=env)

    cmake_build = work / "cmake"
    run(
        [
            "cmake",
            "-S",
            str(source / "examples" / "c"),
            "-B",
            str(cmake_build),
            "-DCMAKE_PREFIX_PATH=" + str(installed),
            "-Dseams_DIR=" + str(configs[0].parent),
            "-DCMAKE_BUILD_TYPE=Release",
        ]
    )
    run(["cmake", "--build", str(cmake_build)])
    cmake_bin = cmake_build / "quickstart"
    env = library_dirs(cmake_bin, env, [stage, build])
    run([str(cmake_bin), str(sample)], env=env)

    find_build = work / "cmake-find"
    run(
        [
            "cmake",
            "-S",
            str(find_check),
            "-B",
            str(find_build),
            "-DCMAKE_PREFIX_PATH=" + str(installed),
            "-Dseams_DIR=" + str(configs[0].parent),
            "-DSEAMS_EXPECT_VERSION=" + version,
        ]
    )
    run(["cmake", "--build", str(find_build)])
    find_bin = find_build / "version"
    env = library_dirs(find_bin, env, [stage, build])
    run([str(find_bin), version], env=env)
    print("installed C quickstart ok")


if __name__ == "__main__":
    main()
