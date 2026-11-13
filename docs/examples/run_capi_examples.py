#!/usr/bin/env python3
"""Compile and run the six C API group examples.

Usage: run_capi_examples.py ROOT BUILD
"""

import os
import shutil
import subprocess
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[2]
EXAMPLES = (
    "capi_atoms.c",
    "capi_neighbours.c",
    "capi_order.c",
    "capi_rings.c",
    "capi_cages.c",
    "capi_errors.c",
)


def main():
    if len(sys.argv) != 3:
        raise SystemExit("usage: run_capi_examples.py ROOT BUILD")
    build = Path(sys.argv[2])
    cc = shutil.which("cc") or shutil.which("gcc")
    if cc is None:
        raise SystemExit("cc is required")
    work = build / "docs-capi-examples"
    work.mkdir(parents=True, exist_ok=True)
    env = os.environ.copy()
    lib = str(build / "src")
    env["LD_LIBRARY_PATH"] = lib + os.pathsep + env.get("LD_LIBRARY_PATH", "")
    for name in EXAMPLES:
        binary = work / name.replace(".c", "")
        subprocess.check_call(
            [
                cc,
                "-I",
                str(ROOT / "include"),
                "-I",
                lib,
                str(ROOT / "docs" / "examples" / name),
                "-L" + lib,
                "-lseams",
                "-Wl,-rpath," + lib,
                "-o",
                str(binary),
            ]
        )
        subprocess.check_call([str(binary)], cwd=ROOT, env=env)
        print(name, "ok")
    print("capi examples ok")


if __name__ == "__main__":
    main()
