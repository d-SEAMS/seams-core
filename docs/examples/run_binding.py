#!/usr/bin/env python3
"""Run one documented language example against a seams-core build."""

import os
import shutil
import subprocess
import sys
import tempfile
from pathlib import Path

ROOT = Path(__file__).resolve().parents[2]


def library_path(build):
    parts = [str(build / "src"), str(build / "bindings")]
    current = os.environ.get("LD_LIBRARY_PATH")
    if current:
        parts.append(current)
    return os.pathsep.join(parts)


def run_fortran(build):
    gfortran = shutil.which("gfortran")
    if not gfortran:
        raise SystemExit("gfortran is not on PATH")
    lib = build / "src"
    work = Path(tempfile.mkdtemp(prefix="seams-fortran-"))
    subprocess.check_call(
        [gfortran, "-J", str(work), "-c", str(ROOT / "bindings/fortran/seams.f90"), "-o", str(work / "seams.o")],
        cwd=work,
    )
    exe = work / "test_seams"
    subprocess.check_call(
        [
            gfortran,
            "-I",
            str(work),
            str(ROOT / "bindings/fortran/test_seams.f90"),
            str(work / "seams.o"),
            "-o",
            str(exe),
            "-L" + str(lib),
            "-lseams",
            "-ldl",
            "-Wl,-rpath," + str(lib),
        ],
        cwd=work,
    )
    env = os.environ.copy()
    env["LD_LIBRARY_PATH"] = library_path(build)
    env["OMP_NUM_THREADS"] = "1"
    subprocess.check_call([str(exe), str(ROOT)], env=env)


def run_rust(build):
    cargo = shutil.which("cargo")
    if not cargo:
        raise SystemExit("cargo is not on PATH")
    env = os.environ.copy()
    env["SEAMS_LIB_DIR"] = str(build / "src")
    env["LD_LIBRARY_PATH"] = library_path(build)
    env["SEAMS_CPP_ORACLE"] = str(build / "bindings" / "binding_oracle")
    env["SEAMS_SOURCE_ROOT"] = str(ROOT)
    env["OMP_NUM_THREADS"] = "1"
    subprocess.check_call(
        [
            cargo,
            "test",
            "--offline",
            "--manifest-path",
            str(ROOT / "bindings/rust/Cargo.toml"),
            "--",
            "--test-threads=1",
        ],
        env=env,
        cwd=ROOT,
    )


def rscript():
    found = shutil.which("Rscript")
    if found:
        return found
    eessi = Path(
        "/cvmfs/software.eessi.io/versions/2026.06/software/linux/x86_64/"
        "intel/sapphirerapids/software/R/4.6.1-gfbf-2026.1/bin/Rscript"
    )
    if eessi.is_file():
        return str(eessi)
    raise SystemExit("Rscript is not on PATH")


def install_r_package(build):
    """Install bindings/r/dseams once. Parallel tests share the library dir."""
    import fcntl

    script = rscript()
    r_cmd = str(Path(script).with_name("R"))
    dest = Path(build) / "docs-r-lib"
    dest.mkdir(parents=True, exist_ok=True)
    env = os.environ.copy()
    env["DSEAMS_ROOT"] = str(ROOT)
    env["SEAMS_INCLUDE"] = str(ROOT / "bindings/cxx")
    env["SEAMS_LIBDIR"] = str(Path(build) / "bindings")
    env["LD_LIBRARY_PATH"] = library_path(Path(build))
    env["R_LIBS"] = str(dest)
    env["R_LIBS_USER"] = str(dest)
    with open(dest / ".install.lock", "a", encoding="utf-8") as handle:
        fcntl.flock(handle, fcntl.LOCK_EX)
        stale = dest / "00LOCK-dseams"
        if stale.exists():
            shutil.rmtree(stale, ignore_errors=True)
        if not (dest / "dseams" / "DESCRIPTION").is_file():
            subprocess.check_call(
                [r_cmd, "CMD", "INSTALL", "-l", str(dest), str(ROOT / "bindings/r/dseams")],
                env=env,
            )
    return script, env


def run_r(build):
    script, env = install_r_package(build)
    subprocess.check_call([script, str(ROOT / "docs/examples/r_api.R")], env=env)


def julia_bin():
    found = shutil.which("julia")
    if found:
        return found
    home = Path.home() / ".local/julia"
    if home.is_dir():
        matches = sorted(home.glob("julia-*/bin/julia"))
        if matches:
            return str(matches[-1])
    raise SystemExit("julia is not on PATH")


def run_julia(build):
    julia = julia_bin()
    env = os.environ.copy()
    env["DSEAMS_ROOT"] = str(ROOT)
    env["SEAMS_LANG_LIB"] = str(build / "bindings" / "libseams_lang.so")
    env["LD_LIBRARY_PATH"] = library_path(build)
    subprocess.check_call(
        [julia, "--project=" + str(ROOT / "bindings/julia/Dseams"), str(ROOT / "docs/examples/julia_api.jl")],
        env=env,
        cwd=ROOT,
    )


def main():
    if len(sys.argv) != 3:
        raise SystemExit("usage: run_binding.py fortran|rust|r|julia BUILD")
    kind, build = sys.argv[1], Path(sys.argv[2])
    {"fortran": run_fortran, "rust": run_rust, "r": run_r, "julia": run_julia}[kind](build)


if __name__ == "__main__":
    main()
