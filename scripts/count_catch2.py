#!/usr/bin/env python3
"""Count Catch2 cases from a build, not from a handwritten total.

Each test binary is asked for ``--list-tests``. A case is a line indented
two spaces whose text is not a ``[tag]``. The paper's supplement claimed
41 binaries and 509 cases; this script prints the count the binaries
actually register.

    python3 scripts/count_catch2.py [BUILD_DIR]

BUILD_DIR defaults to ``build`` under the repository root, or ``SEAMS_BUILD``.
"""

from __future__ import annotations

import json
import os
import subprocess
import sys
from pathlib import Path


def repo_root() -> Path:
    return Path(__file__).resolve().parents[1]


def cases_of(binary: Path) -> list[str]:
    proc = subprocess.run(
        [str(binary), "--list-tests"],
        check=False,
        capture_output=True,
        text=True,
    )
    if proc.returncode != 0:
        raise SystemExit(f"{binary} --list-tests exited {proc.returncode}\n{proc.stderr}")
    names = []
    for line in proc.stdout.splitlines():
        if line.startswith("  ") and not line.startswith("   ") and not line.startswith("  ["):
            names.append(line.strip())
    return names


def main() -> int:
    root = repo_root()
    build = Path(sys.argv[1] if len(sys.argv) > 1 else os.environ.get("SEAMS_BUILD", root / "build"))
    if not build.is_absolute():
        build = (Path.cwd() / build).resolve()
    tests = build / "tests"
    if not tests.is_dir():
        print(f"no test directory at {tests}", file=sys.stderr)
        return 1
    binaries = sorted(
        p for p in tests.iterdir() if p.is_file() and os.access(p, os.X_OK) and p.name.startswith("test_")
    )
    rows = []
    total = 0
    for binary in binaries:
        names = cases_of(binary)
        rows.append({"binary": binary.name, "cases": len(names)})
        total += len(names)
    report = {
        "build": str(build),
        "binaries": len(rows),
        "cases": total,
        "per_binary": rows,
    }
    json.dump(report, sys.stdout, indent=2)
    sys.stdout.write("\n")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
