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


def _git_toplevel(root: Path) -> Path | None:
    cmd = ["git", "-C", str(root), "rev-parse", "--show-toplevel"]
    try:
        done = subprocess.run(cmd, check=False, capture_output=True, text=True)
    except OSError:
        return None
    if done.returncode != 0:
        return None
    return Path(done.stdout.strip()).resolve()


def _from_git(root: Path, specs: tuple[str, ...]) -> list[str]:
    cmd = ["git", "-C", str(root), "ls-files", "-z"]
    if specs:
        cmd += ["--", *specs]
    try:
        done = subprocess.run(cmd, check=False, capture_output=True)
    except OSError as exc:
        raise RuntimeError("git is unavailable: {}".format(exc)) from exc
    if done.returncode != 0:
        err = done.stderr.decode(errors="replace").strip()
        raise RuntimeError("git ls-files failed: {}".format(err or "unknown error"))
    return [row.decode() for row in done.stdout.split(b"\0") if row]


def _from_manifest(root: Path, specs: tuple[str, ...]) -> list[str]:
    path = root / MANIFEST
    rows = [row.decode() for row in path.read_bytes().split(b"\0") if row]
    return [row for row in rows if _matches(row, specs)]


def list_tracked(root: Path, *specs: str) -> list[str]:
    """git-tracked files under ``specs``; falls back to the dist manifest."""
    root = root.resolve()
    if _git_toplevel(root) == root:
        return _from_git(root, specs)
    path = root / MANIFEST
    if not path.is_file():
        raise RuntimeError(
            "source root is not a git toplevel and {} is missing".format(MANIFEST.as_posix())
        )
    return _from_manifest(root, specs)
