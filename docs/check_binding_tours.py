#!/usr/bin/env python3
# Copyright (c) d-SEAMS developers
"""Run the Rust, Fortran, and Python call-set examples and require the docs to quote them."""

import os
import shutil
import subprocess
import sys
import tempfile
from pathlib import Path

REPO = Path(__file__).resolve().parents[1]


def fail(message):
    print(message, file=sys.stderr)
    raise SystemExit(1)


def run(cmd, env, cwd):
    done = subprocess.run(cmd, cwd=str(cwd), env=env, text=True, stdout=subprocess.PIPE, stderr=subprocess.STDOUT)
    if done.returncode != 0:
        fail("exit {} for {}\n{}".format(done.returncode, cmd, done.stdout[-1200:]))
    return done.stdout


def tagged(text, prefix):
    found = []
    for line in text.splitlines():
        if prefix in line:
            found.append(line[line.index(prefix) :].strip())
    if not found:
        fail("no lines starting with {}".format(prefix))
    return found


def require(lines, pages):
    bodies = [(page, page.read_text()) for page in pages]
    for line in lines:
        for page, body in bodies:
            if line not in body:
                fail("{} is missing from {}".format(line, page.relative_to(REPO)))


def fortran(build, env):
    gfortran = shutil.which("gfortran")
    if not gfortran:
        fail("gfortran is not on PATH")
    work = Path(tempfile.mkdtemp(prefix="seams-fortran-"))
    lib = build / "src"
    run(
        [gfortran, "-J", str(work), "-c", str(REPO / "bindings/fortran/src/seams_c.f90"), "-o", str(work / "seams.o")],
        env,
        REPO,
    )
    exe = work / "test_seams"
    run(
        [
            gfortran, "-I", str(work),
            str(REPO / "bindings/fortran/test_seams.f90"),
            str(work / "seams.o"),
            "-o", str(exe),
            "-L" + str(lib), "-lseams", "-Wl,-rpath," + str(lib),
        ],
        env,
        REPO,
    )
    child = env.copy()
    child["LD_LIBRARY_PATH"] = os.pathsep.join(
        [str(lib), str(build / "bindings"), child.get("LD_LIBRARY_PATH", "")]
    )
    child["OMP_NUM_THREADS"] = "1"
    return run([str(exe), str(REPO)], child, REPO)


def rust(build, env):
    cargo = shutil.which("cargo")
    if not cargo:
        fail("cargo is not on PATH")
    child = env.copy()
    child["SEAMS_LIB_DIR"] = str(build / "src")
    child["SEAMS_CPP_ORACLE"] = str(build / "bindings" / "binding_oracle")
    child["SEAMS_SOURCE_ROOT"] = str(REPO)
    child["OMP_NUM_THREADS"] = "1"
    child["LD_LIBRARY_PATH"] = os.pathsep.join(
        [str(build / "src"), str(build / "bindings"), child.get("LD_LIBRARY_PATH", "")]
    )
    return run(
        [
            cargo, "test", "--offline", "--manifest-path", str(REPO / "bindings/rust/Cargo.toml"),
            "--test", "cpp_match", "--", "--nocapture", "--test-threads=1",
        ],
        child,
        REPO,
    )


def python(build, env):
    child = env.copy()
    child["LD_LIBRARY_PATH"] = os.pathsep.join(
        [str(build / "src"), child.get("LD_LIBRARY_PATH", "")]
    )
    child["OMP_NUM_THREADS"] = "1"
    return run(
        [sys.executable, str(REPO / "docs/examples/python_api.py"), str(build / "src")],
        child,
        REPO,
    )


def main():
    if len(sys.argv) != 2:
        fail("usage: check_binding_tours.py BUILD")
    build = Path(sys.argv[1])
    env = os.environ.copy()
    require(tagged(python(build, env), "python call "), [
        REPO / "docs/tutorials/python.md",
        REPO / "docs/api/python.md",
    ])
    require(tagged(fortran(build, env), "fortran call "), [
        REPO / "docs/tutorials/fortran.md",
        REPO / "docs/api/fortran.md",
    ])
    require(tagged(rust(build, env), "rustcall "), [
        REPO / "docs/tutorials/rust.md",
        REPO / "docs/api/rust.md",
    ])
    print("binding tours ok")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
