# Copyright (c) d-SEAMS developers
"""A -Dwith_tests=false configure of the source tree must register no tests.

A consumer that builds seams-core as a subproject with -Dwith_tests=false
inherits every test seams-core registers into its own meson test run, so
none may be registered outside get_option('with_tests'). Wraps are not
downloaded, so the configure writes nothing into SOURCE_ROOT/subprojects;
it uses the subprojects already there and fails if one is missing. WORK_DIR is
removed first and afterwards.

Usage: with_tests_false.py MESON SOURCE_ROOT WORK_DIR
"""

import json
import shutil
import subprocess
import sys
from pathlib import Path


def main(argv: list[str]) -> int:
    if len(argv) != 4:
        print(__doc__, file=sys.stderr)
        return 2
    meson, source, work = argv[1], Path(argv[2]), Path(argv[3])
    shutil.rmtree(work, ignore_errors=True)
    try:
        setup = subprocess.run(
            [meson, "setup", str(work), str(source), "-Dwith_tests=false", "--wrap-mode=nodownload"],
            capture_output=True,
            text=True,
        )
        if setup.returncode != 0:
            print(setup.stdout + setup.stderr, file=sys.stderr)
            print(f"meson setup -Dwith_tests=false exited {setup.returncode}", file=sys.stderr)
            return 1
        listed = subprocess.run(
            [meson, "introspect", "--tests", str(work)],
            capture_output=True,
            text=True,
            check=True,
        )
        names = sorted(t["name"] for t in json.loads(listed.stdout))
    finally:
        shutil.rmtree(work, ignore_errors=True)
    for name in names:
        print(f"registered with -Dwith_tests=false: {name}", file=sys.stderr)
    if names:
        return 1
    print("-Dwith_tests=false registers no tests")
    return 0


if __name__ == "__main__":
    raise SystemExit(main(sys.argv))
