# Copyright (c) d-SEAMS developers
"""Import the built nanobind module and run its unittest. Exit 77 if it was not built."""

import os
import shutil
import sys
from pathlib import Path

try:
    import numpy  # noqa: F401
except ImportError:
    print("numpy is not installed")
    sys.exit(77)

build = Path(os.environ.get("DSEAMS_BUILD", "build"))
src = Path(os.environ.get("DSEAMS_SRC", Path(__file__).resolve().parents[2]))
# A conda python and the extension must share a stable ABI tag. A module
# built for another interpreter is a skip, not a failed import.
tag = sys.implementation.cache_tag  # cpython-312
cores = [path for path in build.rglob("_core*.so") if tag in path.name]
if not cores:
    built = [path.name for path in build.rglob("_core*.so")]
    print(f"dseams._core was not built for {tag}; found {built or 'nothing'}")
    sys.exit(77)

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
    # The unittest loader puts bindings/python on sys.path ahead of the
    # staged tree. The extension has to sit next to that package too.
    sibling = root / "dseams" / so.name
    if so.resolve() != sibling.resolve():
        shutil.copy(so, sibling)

if os.environ.get("DSEAMS_PY_READY") != "1":
    env = os.environ.copy()
    lib = build / "src"
    env["LD_LIBRARY_PATH"] = str(lib) + os.pathsep + env.get("LD_LIBRARY_PATH", "")
    ctypes_src = root / "src"
    env["PYTHONPATH"] = os.pathsep.join(
        [str(stage), str(ctypes_src), env.get("PYTHONPATH", "")]
    )
    env["DSEAMS_PY_READY"] = "1"
    os.execve(sys.executable, [sys.executable, __file__], env)

import unittest

try:
    suite = unittest.defaultTestLoader.discover(str(root / "tests"), pattern="test_*.py")
    result = unittest.TextTestRunner(verbosity=1).run(suite)
    sys.exit(0 if result.wasSuccessful() else 1)
finally:
    for so in cores:
        sibling = root / "dseams" / so.name
        if sibling.is_file() and so.resolve() != sibling.resolve():
            sibling.unlink()
