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


LOCK_FIXTURE = """\
[[package]]
name = "cfg-if"
version = "1.0.5"
source = "registry+https://github.com/rust-lang/crates.io-index"

[[package]]
name = "linkcell"
version = "0.6.0"

[[package]]
name = "wasip2"
version = "1.0.4+wasi-0.2.12"
source = "registry+https://github.com/rust-lang/crates.io-index"
"""


def _assert_raises(message_part: str, fn) -> None:
    try:
        fn()
    except SystemExit as exc:
        text = str(exc)
        _assert(message_part in text, "expected {!r} in {!r}".format(message_part, text))
        return
    raise RuntimeError("expected SystemExit containing {!r}".format(message_part))


def _lock_names_come_from_the_lock() -> None:
    got = dist_vendor.lock_sourced_crates(LOCK_FIXTURE)
    _assert(
        got == [("cfg-if", "1.0.5"), ("wasip2", "1.0.4+wasi-0.2.12")],
        "sourced crates were {}".format(got),
    )


def _vendor_must_hold_every_lock_crate() -> None:
    with tempfile.TemporaryDirectory() as tmp:
        project = Path(tmp) / "linkcell"
        project.mkdir()
        (project / "Cargo.lock").write_text(LOCK_FIXTURE, encoding="utf-8")
        _assert_raises("cfg-if 1.0.5", lambda: dist_vendor.require_lock_vendored(project))
        for name in ("cfg-if", "wasip2"):
            (project / "vendor" / name).mkdir(parents=True)
        dist_vendor.require_lock_vendored(project)


def _install_copies_the_vendored_tree() -> None:
    with tempfile.TemporaryDirectory() as tmp:
        base = Path(tmp)
        dist = base / "dist"
        dest = base / "dest"
        src = dist / "subprojects" / "linkcell"
        (src / "vendor" / "cfg-if").mkdir(parents=True)
        (src / "vendor" / "wasip2").mkdir()
        (src / ".cargo").mkdir()
        (src / ".cargo" / "config.toml").write_text(dist_vendor.VENDOR_CONFIG, encoding="utf-8")
        (src / "Cargo.lock").write_text(LOCK_FIXTURE, encoding="utf-8")
        target = dest / "subprojects" / "linkcell"
        target.mkdir(parents=True)
        (target / "Cargo.lock").write_text(LOCK_FIXTURE, encoding="utf-8")
        (target / "vendor").mkdir()
        dist_vendor.install_vendored_cargo(dist, dest)
        _assert((target / "vendor" / "cfg-if").is_dir(), "cfg-if was not installed")
        _assert((target / ".cargo" / "config.toml").is_file(), "cargo config was not installed")
        (target / "Cargo.lock").write_text(LOCK_FIXTURE + "\n", encoding="utf-8")
        _assert_raises(
            "Cargo.lock does not match",
            lambda: dist_vendor.install_vendored_cargo(dist, dest),
        )


def _install_rejects_an_empty_dist() -> None:
    with tempfile.TemporaryDirectory() as tmp:
        base = Path(tmp)
        (base / "dist" / "subprojects").mkdir(parents=True)
        (base / "dest").mkdir()
        _assert_raises(
            "no cargo vendor directories",
            lambda: dist_vendor.install_vendored_cargo(base / "dist", base / "dest"),
        )


def main() -> int:
    _foreign_git_dir_and_index_are_ignored()
    _lock_names_come_from_the_lock()
    _vendor_must_hold_every_lock_crate()
    _install_copies_the_vendored_tree()
    _install_rejects_an_empty_dist()
    print("dist_vendor manifest checks ok")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
