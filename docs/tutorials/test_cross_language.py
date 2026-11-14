"""One task, three languages.

CHILL+ on frame 1 of ``mW_cubic.lammpstrj``, LAMMPS type 1.
Python (pydseamslib) and R and Julia (``libseams_lang``) default to the
cutoff neighbour list. The count is 4096 cubic in each language.
"""

import os
import shutil
import subprocess
from pathlib import Path

import pytest

pydseams = pytest.importorskip("pydseams")

ROOT = Path(os.environ.get("DSEAMS_ROOT", Path(__file__).resolve().parents[2]))


def _lib_path() -> str:
    if os.environ.get("SEAMS_LANG_LIB"):
        return os.environ["SEAMS_LANG_LIB"]
    candidate = ROOT / "build" / "bindings" / "libseams_lang.so"
    return str(candidate)


def test_python_cubic_count(cubic):
    frame = pydseams.read(cubic, atom_type=1)
    assert frame.chill_plus().cubic == 4096


def _env():
    env = os.environ.copy()
    libdir = str(Path(_lib_path()).resolve().parent)
    extra = ":".join(
        [
            libdir,
            str(ROOT / "build" / "src"),
            str(ROOT / "build" / "subprojects" / "readcon-core"),
            str(ROOT / "build" / "subprojects" / "readcon-db"),
        ]
    )
    previous = env.get("LD_LIBRARY_PATH")
    env["LD_LIBRARY_PATH"] = extra + (":" + previous if previous else "")
    env["DSEAMS_ROOT"] = str(ROOT)
    env["SEAMS_LANG_LIB"] = str(Path(_lib_path()).resolve())
    env["SEAMS_INCLUDE"] = str(ROOT / "bindings" / "cxx")
    env["SEAMS_LIBDIR"] = libdir
    env["R_LIBS_USER"] = env.get("R_LIBS_USER", str(Path.home() / "R" / "library"))
    return env


@pytest.mark.skipif(shutil.which("Rscript") is None, reason="Rscript is not installed")
def test_r_cubic_count():
    probe = subprocess.run(
        ["Rscript", "-e", "library(dseams)"],
        check=False,
        capture_output=True,
        text=True,
        env=_env(),
    )
    if probe.returncode != 0:
        pytest.skip("the dseams R package is not installed")
    script = r"""
root <- Sys.getenv("DSEAMS_ROOT")
suppressPackageStartupMessages(library(dseams))
frame <- read_frame(file.path(root, "input/traj/mW_cubic.lammpstrj"), type = 1L)
ice <- chill_plus(frame)
stopifnot(sum(ice$name == "cubic") == 4096L)
"""
    proc = subprocess.run(
        ["Rscript", "-e", script],
        check=False,
        capture_output=True,
        text=True,
        env=_env(),
    )
    assert proc.returncode == 0, proc.stderr


@pytest.mark.skipif(shutil.which("julia") is None, reason="julia is not installed")
def test_julia_cubic_count():
    julia = shutil.which("julia")
    proc = subprocess.run(
        [
            julia,
            f"--project={ROOT / 'bindings' / 'julia' / 'Dseams'}",
            "-e",
            'include("dev/develop_jll.jl"); include("test/runtests.jl")',
        ],
        check=False,
        capture_output=True,
        text=True,
        env=_env(),
        cwd=str(ROOT / "bindings" / "julia" / "Dseams"),
    )
    assert proc.returncode == 0, proc.stdout + proc.stderr
