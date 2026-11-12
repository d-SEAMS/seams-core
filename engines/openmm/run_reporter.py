# Copyright (c) d-SEAMS developers
"""Run the OpenMM reporter. Skip only when OpenMM itself is absent."""

import importlib.util
import os
import runpy
import shutil
import sys
from pathlib import Path

if importlib.util.find_spec("openmm") is None:
    print("openmm is not installed")
    sys.exit(77)

build = Path(os.environ.get("DSEAMS_BUILD", "build"))
src = Path(os.environ.get("DSEAMS_SRC", Path(__file__).resolve().parents[2]))
cores = list(build.rglob("_core*.so"))
if not cores:
    print("openmm is installed but dseams._core was not built")
    sys.exit(1)

stage = build / "python-stage"
pkg = stage / "dseams"
pkg.mkdir(parents=True, exist_ok=True)
root = src / "bindings" / "python"
for name in ("__init__.py", "__main__.py"):
    shutil.copy(root / "dseams" / name, pkg / name)
shutil.copy(root / "frames.py", pkg / "frames.py")
for so in cores:
    dest = pkg / so.name
    if so.resolve() != dest.resolve():
        shutil.copy(so, dest)

libs = list(build.rglob("libseams_lang.so"))
if not libs:
    print("openmm is installed but libseams_lang.so was not built")
    sys.exit(1)

if os.environ.get("DSEAMS_PY_READY") != "1":
    env = os.environ.copy()
    lib = str(libs[0].parent)
    src_lib = str(build / "src")
    env["LD_LIBRARY_PATH"] = lib + os.pathsep + src_lib + os.pathsep + env.get("LD_LIBRARY_PATH", "")
    env["PYTHONPATH"] = str(stage) + os.pathsep + str(src / "adapters" / "python") + os.pathsep + env.get("PYTHONPATH", "")
    env["SEAMS_LANG_LIB"] = str(libs[0])
    env["DSEAMS_PY_READY"] = "1"
    os.execve(sys.executable, [sys.executable, __file__], env)

import dseams

if not dseams.version():
    sys.exit("dseams.version is empty")

here = Path(__file__).resolve().parent
sys.path.insert(0, str(here))
sys.path.insert(0, str(src / "adapters" / "python"))
runpy.run_path(str(here / "test_reporter.py"), run_name="__main__")
