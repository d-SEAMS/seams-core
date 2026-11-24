#!/usr/bin/env python3
# Copyright (c) d-SEAMS developers
"""dist_vendor writes the source root's tracked files, whatever git env it inherits."""

from __future__ import annotations

import importlib.util
import os
import subprocess
import tempfile
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]


def _load_dist_vendor():
    path = ROOT / "scripts" / "dist_vendor.py"
    spec = importlib.util.spec_from_file_location("dist_vendor", path)
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


dist_vendor = _load_dist_vendor()

REPO_VARS = (
    "GIT_DIR",
    "GIT_WORK_TREE",
    "GIT_INDEX_FILE",
    "GIT_COMMON_DIR",
    "GIT_OBJECT_DIRECTORY",
    "GIT_ALTERNATE_OBJECT_DIRECTORIES",
    "GIT_NAMESPACE",
    "GIT_PREFIX",
    "GIT_IMPLICIT_WORK_TREE",
)


def _clean_env() -> dict[str, str]:
    env = {key: value for key, value in os.environ.items() if key not in REPO_VARS}
    env["GIT_CONFIG_NOSYSTEM"] = "1"
    env["GIT_CONFIG_GLOBAL"] = os.devnull
    return env


def _make_repo(path: Path, files: list[str]) -> None:
    env = _clean_env()
    path.mkdir(parents=True)
    subprocess.run(["git", "init", "--quiet", str(path)], check=True, env=env)
    for name in files:
        target = path / name
        target.parent.mkdir(parents=True, exist_ok=True)
        target.write_text(name + "\n", encoding="utf-8")
    subprocess.run(["git", "-C", str(path), "add", "--", *files], check=True, env=env)


def _assert(condition: bool, message: str) -> None:
    if not condition:
        raise RuntimeError(message)


def _foreign_git_dir_and_index_are_ignored() -> None:
    with tempfile.TemporaryDirectory() as tmp:
        base = Path(tmp).resolve()
        source = base / "source"
        foreign = base / "foreign"
        dist = base / "dist"
        _make_repo(source, ["README.md", "docs/source/page.rst"])
        _make_repo(foreign, ["foreign-only.txt", "docs/foreign.md"])
        dist.mkdir()

        foreign_index = foreign / ".git" / "index"
        index_bytes = foreign_index.read_bytes()
        saved = {key: os.environ.get(key) for key in ("GIT_DIR", "GIT_INDEX_FILE")}
        os.environ["GIT_DIR"] = str(foreign / ".git")
        os.environ["GIT_INDEX_FILE"] = str(foreign_index)
        try:
            dist_vendor.write_tracked_manifest(source, dist)
        finally:
            for key, old in saved.items():
                if old is None:
                    os.environ.pop(key, None)
                else:
                    os.environ[key] = old

        manifest = dist / "docs" / ".cache" / "tracked-files.zlist"
        rows = sorted(row.decode() for row in manifest.read_bytes().split(b"\0") if row)
        _assert(
            rows == ["README.md", "docs/source/page.rst"],
            "manifest did not list only the source root's tracked files: {}".format(rows),
        )
        _assert(foreign_index.read_bytes() == index_bytes, "foreign GIT_INDEX_FILE was modified")


def main() -> int:
    _foreign_git_dir_and_index_are_ignored()
    print("dist_vendor manifest checks ok")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
