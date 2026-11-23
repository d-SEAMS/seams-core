#!/usr/bin/env python3
# Copyright (c) d-SEAMS developers
"""Run a tracked-file docs check when the source tree has no Git metadata."""

import subprocess
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))
import tracked_files


def main():
    if len(sys.argv) != 2:
        raise SystemExit("usage: run_doc_check.py CHECK")
    root = Path(__file__).resolve().parents[1]
    check = Path(sys.argv[1]).resolve()
    # Fail fast with the same tracked-file resolution the checks use.
    tracked_files.list_tracked(root, "docs")
    return subprocess.run([sys.executable, str(check)], cwd=root).returncode


if __name__ == "__main__":
    raise SystemExit(main())
