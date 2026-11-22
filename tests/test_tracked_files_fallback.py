#!/usr/bin/env python3
# Copyright (c) d-SEAMS developers
"""Checks for tracked-files fallback behavior outside git work trees."""

from __future__ import annotations

import os
import subprocess
import sys
import tempfile
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "scripts"))
import tracked_files  # noqa: E402


def _write_manifest(root: Path, rows: list[str]) -> None:
    path = root / tracked_files.MANIFEST
    path.parent.mkdir(parents=True, exist_ok=True)
    data = b"\0".join(row.encode() for row in rows) + b"\0"
    path.write_bytes(data)


def _without_git() -> dict[str, str]:
    env = os.environ.copy()
    env["PATH"] = ""
    return env


def _assert(condition: bool, message: str) -> None:
    if not condition:
        raise RuntimeError(message)


def _no_git_with_manifest() -> None:
    with tempfile.TemporaryDirectory() as tmp:
        root = Path(tmp)
        _write_manifest(root, ["docs/source/alpha.rst", "docs/guide.md", "README.md"])
        old = os.environ.copy()
        os.environ.clear()
        os.environ.update(_without_git())
        try:
            rows = tracked_files.list_tracked(root, "docs/source")
        finally:
            os.environ.clear()
            os.environ.update(old)
        _assert(rows == ["docs/source/alpha.rst"], "no-git manifest fallback did not match docs/source")


def _inside_unrelated_git_tree() -> None:
    with tempfile.TemporaryDirectory() as tmp:
        outer = Path(tmp) / "outer"
        outer.mkdir()
        subprocess.run(["git", "init", "--quiet"], cwd=outer, check=True)
        (outer / "tracked.txt").write_text("x\n", encoding="utf-8")
        subprocess.run(["git", "add", "tracked.txt"], cwd=outer, check=True)
        inner = outer / "sdist"
        inner.mkdir()
        _write_manifest(inner, ["docs/source/from-manifest.rst"])
        rows = tracked_files.list_tracked(inner, "docs/source")
        _assert(
            rows == ["docs/source/from-manifest.rst"],
            "sdist under unrelated git tree did not use manifest",
        )


def _no_git_and_no_manifest() -> None:
    with tempfile.TemporaryDirectory() as tmp:
        root = Path(tmp)
        old = os.environ.copy()
        os.environ.clear()
        os.environ.update(_without_git())
        try:
            try:
                tracked_files.list_tracked(root, "docs/source")
            except RuntimeError as exc:
                _assert("missing" in str(exc), "missing manifest error is not loud enough")
            else:
                raise RuntimeError("expected RuntimeError when git and manifest are unavailable")
        finally:
            os.environ.clear()
            os.environ.update(old)


def main() -> int:
    _no_git_with_manifest()
    _inside_unrelated_git_tree()
    _no_git_and_no_manifest()
    print("tracked_files fallback checks ok")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
