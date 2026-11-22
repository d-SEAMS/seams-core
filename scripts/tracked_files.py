#!/usr/bin/env python3
# Copyright (c) d-SEAMS developers
"""Resolve tracked files from git, or from a dist-time manifest."""

from __future__ import annotations

import subprocess
from pathlib import Path

MANIFEST = Path("docs/.cache/tracked-files.zlist")


def _matches(path: str, specs: tuple[str, ...]) -> bool:
    if not specs:
        return True
    for spec in specs:
        clean = spec.strip("/")
        if path == clean or path.startswith(clean + "/"):
            return True
    return False


def _from_git(root: Path, specs: tuple[str, ...]) -> list[str] | None:
    cmd = ["git", "-C", str(root), "ls-files", "-z"]
    if specs:
        cmd += ["--", *specs]
    done = subprocess.run(cmd, check=False, capture_output=True)
    if done.returncode != 0:
        return None
    return [row.decode() for row in done.stdout.split(b"\0") if row]


def _from_manifest(root: Path, specs: tuple[str, ...]) -> list[str]:
    path = root / MANIFEST
    rows = [row.decode() for row in path.read_bytes().split(b"\0") if row]
    return [row for row in rows if _matches(row, specs)]


def list_tracked(root: Path, *specs: str) -> list[str]:
    """git-tracked files under ``specs``; falls back to the dist manifest."""
    rows = _from_git(root, specs)
    if rows is not None:
        return rows
    path = root / MANIFEST
    if not path.is_file():
        raise RuntimeError("git is unavailable and {} is missing".format(MANIFEST.as_posix()))
    return _from_manifest(root, specs)
