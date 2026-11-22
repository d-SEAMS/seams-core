#!/usr/bin/env python3
# Copyright (c) d-SEAMS developers
"""Resolve tracked files from git, or from a dist-time manifest."""

from __future__ import annotations

import os
import subprocess
from pathlib import Path

MANIFEST = Path("docs/.cache/tracked-files.zlist")
_GIT_ENV_PREFIX = "GIT_"


def _git_env() -> dict[str, str]:
    return {key: value for key, value in os.environ.items() if not key.startswith(_GIT_ENV_PREFIX)}


def _matches(path: str, specs: tuple[str, ...]) -> bool:
    if not specs:
        return True
    for spec in specs:
        clean = spec.strip("/")
        if path == clean or path.startswith(clean + "/"):
            return True
    return False


def _git_worktree(root: Path) -> bool:
    # A work tree has a .git directory; a linked work tree has a .git file.
    # rev-parse exits 128 for "not a repository" and for real failures, and
    # git translates that text, so the split is the marker on disk.
    marker = root / ".git"
    return marker.is_dir() or marker.is_file() or marker.is_symlink()


def _git_toplevel(root: Path) -> Path | None:
    cmd = ["git", "-C", str(root), "rev-parse", "--show-toplevel"]
    try:
        done = subprocess.run(cmd, check=False, capture_output=True, text=True, env=_git_env())
    except FileNotFoundError:
        if _git_worktree(root):
            raise
        return None
    if done.returncode == 0:
        return Path(done.stdout.strip()).resolve()
    # Exit 128 with no marker is "not a repository", in any locale.
    # A marker means this is a work tree and the failure is real.
    if done.returncode == 128 and not _git_worktree(root):
        return None
    err = (done.stderr or "").strip()
    raise RuntimeError("git rev-parse failed: {}".format(err or "unknown error"))


def _from_git(root: Path, specs: tuple[str, ...]) -> list[str]:
    cmd = ["git", "-C", str(root), "ls-files", "-z"]
    if specs:
        cmd += ["--", *specs]
    try:
        done = subprocess.run(cmd, check=False, capture_output=True, env=_git_env())
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
    top = _git_toplevel(root)
    if top == root:
        return _from_git(root, specs)
    # Git ignores a marker that is not a repository and keeps walking. An
    # enclosing repository would then supply the file list.
    if _git_worktree(root):
        raise RuntimeError("{} has a .git marker but git resolved {}".format(root, top))
    path = root / MANIFEST
    if not path.is_file():
        raise RuntimeError(
            "source root is not a git toplevel and {} is missing".format(MANIFEST.as_posix())
        )
    return _from_manifest(root, specs)
