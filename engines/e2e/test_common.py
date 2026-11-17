# Copyright (c) d-SEAMS developers
"""Regression checks for the shared engine end-to-end helpers."""

import os
import sys
import tempfile
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
import e2e_common as e2e  # noqa: E402


def main():
    with tempfile.TemporaryDirectory(prefix="dseams-e2e-common-") as tmp:
        root = Path(tmp)
        build = root / "build"
        plugin = build / "engines" / "plugin.so"
        plugin.parent.mkdir(parents=True)
        plugin.touch()
        away = root / "work"
        away.mkdir()
        old = Path.cwd()
        try:
            os.chdir(root)
            found = e2e.find_build_file("build", "plugin.so")
            os.chdir(away)
            if found != plugin.resolve() or not found.is_file():
                return 1
        finally:
            os.chdir(old)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
