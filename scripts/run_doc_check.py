#!/usr/bin/env python3
# Copyright (c) d-SEAMS developers
"""Run a tracked-file docs check when the source tree has no Git metadata."""

import os
import subprocess
import sys
import tempfile
from pathlib import Path


def main():
    if len(sys.argv) != 2:
        raise SystemExit("usage: run_doc_check.py CHECK")
    root = Path(__file__).resolve().parents[1]
    check = Path(sys.argv[1]).resolve()
    if (root / ".git").exists():
        return subprocess.run([sys.executable, str(check)], cwd=root).returncode
    with tempfile.TemporaryDirectory() as temporary:
        git_dir = Path(temporary) / "git"
        subprocess.run(["git", "init", "--quiet", str(git_dir)], check=True)
        env = os.environ.copy()
        env["GIT_DIR"] = str(git_dir / ".git")
        env["GIT_WORK_TREE"] = str(root)
        paths = [root / "README.md"]
        for base in (root / "docs", root / "changelog.d"):
            for path in base.rglob("*"):
                relative = path.relative_to(root)
                if path.suffix not in {".md", ".org", ".rst"}:
                    continue
                if {".cache", "build", "doctrees", "__pycache__"}.intersection(relative.parts):
                    continue
                if path.is_relative_to(root / "docs" / "source" / "api"):
                    continue
                paths.append(path)
        command = ["git", "add", "-f", "--"] + [str(path.relative_to(root)) for path in paths]
        subprocess.run(command, cwd=root, env=env, check=True)
        return subprocess.run([sys.executable, str(check)], cwd=root, env=env).returncode


if __name__ == "__main__":
    raise SystemExit(main())
