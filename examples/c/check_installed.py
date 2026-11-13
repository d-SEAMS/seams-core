#!/usr/bin/env python3
"""Install libseams and build examples/c against that install."""

import os
import shutil
import subprocess
import sys
from pathlib import Path


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


def library_dirs(binary, env, search):
    """Add directories that hold shared libraries ldd cannot yet see."""
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
    stage = build / "capi-destdir"
    work = build / "capi-work"
    if stage.exists():
        shutil.rmtree(stage)
    work.mkdir(parents=True, exist_ok=True)
    run(["meson", "install", "-C", str(build), "--destdir", str(stage)])

    pc_files = list(stage.rglob("seams.pc"))
    configs = list(stage.rglob("seams-config.cmake"))
    if len(pc_files) != 1 or len(configs) != 1:
        raise SystemExit("install did not produce seams.pc and seams-config.cmake")
    pc = pc_files[0]
    libs_on_disk = list(stage.rglob("libseams.so"))
    if not libs_on_disk:
        raise SystemExit("install did not produce libseams.so")
    libdir = libs_on_disk[0].parent
    env = os.environ.copy()
    env["PKG_CONFIG_PATH"] = str(pc.parent)
    env["PKG_CONFIG_SYSROOT_DIR"] = str(stage)
    cflags = subprocess.check_output(
        ["pkg-config", "--cflags", "seams"], env=env, text=True
    ).split()
    libs = subprocess.check_output(
        ["pkg-config", "--libs", "seams"], env=env, text=True
    ).split()
    binary = work / "quickstart-pkg"
    cc = os.environ.get("CC", "cc")
    run(
        [cc, str(source / "examples" / "c" / "quickstart.c"), "-o", str(binary), *cflags, *libs, "-lm"],
        env=env,
    )
    env["LD_LIBRARY_PATH"] = os.pathsep.join(
        [str(libdir), str(build / "src"), env.get("LD_LIBRARY_PATH", "")]
    )
    env = library_dirs(binary, env, [stage, build])
    sample = source / "input" / "traj" / "agi_water_tiny.lammpstrj"
    run([str(binary), str(sample)], env=env)

    prefix = prefix_of(build)
    installed = stage / Path(prefix).relative_to("/")
    cmake_build = work / "cmake"
    if cmake_build.exists():
        shutil.rmtree(cmake_build)
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
    print("installed C quickstart ok")


if __name__ == "__main__":
    main()
