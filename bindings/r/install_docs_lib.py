#!/usr/bin/env python3
# Copyright (c) d-SEAMS developers
"""Install the dseams R package into BUILD/docs-r-lib once, before the tests.

usage: install_docs_lib.py R SOURCE_ROOT BUILD LANGDIR STAMP

The tests that run R read this library in parallel. Installing it from
each test let one test load a package another was still writing.
"""

import os
import shutil
import subprocess
import sys
import tempfile
from pathlib import Path


def main(argv):
    if len(argv) != 5:
        raise SystemExit(__doc__)
    r_cmd, root, build, langdir, stamp = argv
    root = Path(root)
    dest = Path(build) / "docs-r-lib"
    dest.mkdir(parents=True, exist_ok=True)
    env = os.environ.copy()
    env["DSEAMS_ROOT"] = str(root)
    env["SEAMS_INCLUDE"] = str(root / "bindings" / "cxx")
    env["SEAMS_LIBDIR"] = langdir
    work = Path(tempfile.mkdtemp(prefix="docs-r-lib-", dir=dest))
    try:
        done = subprocess.run(
            [r_cmd, "CMD", "INSTALL", "-l", str(work), str(root / "bindings" / "r" / "dseams")],
            env=env,
            text=True,
            capture_output=True,
        )
        if done.returncode != 0:
            sys.stderr.write(done.stdout + done.stderr)
            return done.returncode
        installed = dest / "dseams"
        if installed.exists():
            shutil.rmtree(installed)
        (work / "dseams").rename(installed)
    finally:
        shutil.rmtree(work, ignore_errors=True)
    Path(stamp).write_text("dseams\n")
    return 0


if __name__ == "__main__":
    raise SystemExit(main(sys.argv[1:]))
