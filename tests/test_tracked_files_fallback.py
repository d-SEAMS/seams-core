#!/usr/bin/env python3
# Copyright (c) d-SEAMS developers
"""Checks for tracked-files fallback behavior outside git work trees."""

from __future__ import annotations

import importlib.util
import os
import subprocess
import sys
import tempfile
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]


def _load_tracked_files():
    impl = os.environ.get("TRACKED_FILES_IMPL")
    path = Path(impl) if impl else (ROOT / "scripts" / "tracked_files.py")
    spec = importlib.util.spec_from_file_location("tracked_files", path)
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


tracked_files = _load_tracked_files()


def _write_manifest(root: Path, rows: list[str]) -> None:
    path = root / tracked_files.MANIFEST
    path.parent.mkdir(parents=True, exist_ok=True)
    data = b"\0".join(row.encode() for row in rows) + b"\0"
    path.write_bytes(data)


def _assert(condition: bool, message: str) -> None:
    if not condition:
        raise RuntimeError(message)


class _PatchRun:
    def __init__(self, fn):
        self._fn = fn
        self._old = None

    def __enter__(self):
        self._old = tracked_files.subprocess.run
        tracked_files.subprocess.run = self._fn
        return self

    def __exit__(self, _exc_type, _exc, _tb):
        tracked_files.subprocess.run = self._old


def _cp(args, code=0, stdout="", stderr=""):
    return subprocess.CompletedProcess(args=args, returncode=code, stdout=stdout, stderr=stderr)


def _no_git_with_manifest() -> None:
    with tempfile.TemporaryDirectory() as tmp:
        root = Path(tmp)
        _write_manifest(root, ["docs/source/alpha.rst", "docs/guide.md", "README.md"])

        def fake_run(*_args, **_kwargs):
            raise FileNotFoundError("git")

        with _PatchRun(fake_run):
            rows = tracked_files.list_tracked(root, "docs/source")
        _assert(rows == ["docs/source/alpha.rst"], "no-git manifest fallback did not match docs/source")


def _inside_unrelated_git_tree() -> None:
    with tempfile.TemporaryDirectory() as tmp:
        outer = Path(tmp) / "outer"
        outer.mkdir()
        (outer / ".git").mkdir()
        (outer / ".git" / "HEAD").write_text("ref: refs/heads/main\n", encoding="utf-8")
        (outer / ".git" / "objects").mkdir()
        (outer / ".git" / "refs").mkdir()
        (outer / ".git" / "refs" / "heads").mkdir(parents=True)
        inner = outer / "sdist"
        inner.mkdir()
        _write_manifest(inner, ["docs/source/from-manifest.rst"])

        def fake_run(args, check=False, capture_output=False, text=False):
            if args[3:] == ["rev-parse", "--show-toplevel"]:
                out = str(outer) + ("\n" if text else "")
                return _cp(args, stdout=out)
            raise RuntimeError("unexpected call {}".format(args))

        with _PatchRun(fake_run):
            rows = tracked_files.list_tracked(inner, "docs/source")
        _assert(
            rows == ["docs/source/from-manifest.rst"],
            "sdist under unrelated git tree did not use manifest",
        )


def _no_git_and_no_manifest() -> None:
    with tempfile.TemporaryDirectory() as tmp:
        root = Path(tmp)

        def fake_run(*_args, **_kwargs):
            raise FileNotFoundError("git")

        with _PatchRun(fake_run):
            try:
                tracked_files.list_tracked(root, "docs/source")
            except RuntimeError as exc:
                _assert("missing" in str(exc), "missing manifest error is not loud enough")
            else:
                raise RuntimeError("expected RuntimeError when git and manifest are unavailable")


def _git_toplevel_uses_git_not_manifest() -> None:
    with tempfile.TemporaryDirectory() as tmp:
        root = Path(tmp)
        _write_manifest(root, ["docs/source/from-manifest.rst"])
        tracked = b"docs/source/from-git.rst\0"

        def fake_run(args, check=False, capture_output=False, text=False):
            if args[3:] == ["rev-parse", "--show-toplevel"]:
                out = str(root) + ("\n" if text else "")
                return _cp(args, stdout=out)
            if args[3] == "ls-files":
                return _cp(args, stdout=tracked)
            raise RuntimeError("unexpected call {}".format(args))

        with _PatchRun(fake_run):
            rows = tracked_files.list_tracked(root, "docs/source")
        _assert(rows == ["docs/source/from-git.rst"], "git-tracked rows were not preferred at git toplevel")


def _git_ls_files_failure_raises() -> None:
    with tempfile.TemporaryDirectory() as tmp:
        root = Path(tmp)
        _write_manifest(root, ["docs/source/from-manifest.rst"])

        def fake_run(args, check=False, capture_output=False, text=False):
            if args[3:] == ["rev-parse", "--show-toplevel"]:
                out = str(root) + ("\n" if text else "")
                return _cp(args, stdout=out)
            if args[3] == "ls-files":
                err = "fatal: index file corrupt"
                return _cp(args, code=128, stderr=err.encode())
            raise RuntimeError("unexpected call {}".format(args))

        with _PatchRun(fake_run):
            try:
                tracked_files.list_tracked(root, "docs/source")
            except RuntimeError as exc:
                _assert("ls-files failed" in str(exc), "ls-files failure was not raised")
            else:
                raise RuntimeError("expected RuntimeError on git ls-files failure")


def _rev_parse_exit_128_raises() -> None:
    with tempfile.TemporaryDirectory() as tmp:
        root = Path(tmp)
        (root / ".git").mkdir()
        _write_manifest(root, ["docs/source/from-manifest.rst"])

        def fake_run(args, check=False, capture_output=False, text=False):
            err = "fatal: detected dubious ownership in repository"
            return _cp(args, code=128, stderr=err)

        with _PatchRun(fake_run):
            try:
                tracked_files.list_tracked(root, "docs/source")
            except RuntimeError as exc:
                _assert("rev-parse failed" in str(exc), "rev-parse exit-128 was not raised")
            else:
                raise RuntimeError("expected RuntimeError on rev-parse exit-128")


def _not_git_repo_uses_manifest() -> None:
    with tempfile.TemporaryDirectory() as tmp:
        root = Path(tmp)
        _write_manifest(root, ["docs/source/from-manifest.rst"])

        def fake_run(args, check=False, capture_output=False, text=False):
            err = "fatal: not a git repository (or any of the parent directories): .git"
            return _cp(args, code=128, stderr=err)

        with _PatchRun(fake_run):
            rows = tracked_files.list_tracked(root, "docs/source")
        _assert(rows == ["docs/source/from-manifest.rst"], "not-a-repository should use manifest")


def _localized_not_a_repo_uses_manifest() -> None:
    messages = (
        "Schwerwiegend: Kein Git-Repository (oder irgendeines der Elternverzeichnisse): .git",
        "fatal\u00a0: ni ceci ni aucun de ses répertoires parents n'est un dépôt git\u00a0: .git",
    )
    for err in messages:
        _assert(
            "not a git repository" not in err.lower(),
            "localized stderr accidentally contains the English phrase",
        )
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            _write_manifest(root, ["docs/source/from-manifest.rst"])

            def fake_run(args, check=False, capture_output=False, text=False, err=err):
                return _cp(args, code=128, stderr=err)

            with _PatchRun(fake_run):
                rows = tracked_files.list_tracked(root, "docs/source")
            _assert(
                rows == ["docs/source/from-manifest.rst"],
                "localized not-a-repository error should use the manifest: {}".format(err),
            )


def _ls_files_oserror_raises() -> None:
    with tempfile.TemporaryDirectory() as tmp:
        root = Path(tmp)
        _write_manifest(root, ["docs/source/from-manifest.rst"])

        def fake_run(args, check=False, capture_output=False, text=False):
            if args[3:] == ["rev-parse", "--show-toplevel"]:
                out = str(root) + ("\n" if text else "")
                return _cp(args, stdout=out)
            if args[3] == "ls-files":
                raise OSError("ls-files")
            raise RuntimeError("unexpected call {}".format(args))

        with _PatchRun(fake_run):
            try:
                tracked_files.list_tracked(root, "docs/source")
            except (OSError, RuntimeError) as exc:
                _assert("ls-files" in str(exc), "ls-files OSError was not raised")
            else:
                raise RuntimeError("ls-files OSError returned rows")


def _rev_parse_probe_does_not_catch_exception() -> None:
    with tempfile.TemporaryDirectory() as tmp:
        root = Path(tmp)
        _write_manifest(root, ["docs/source/from-manifest.rst"])

        def fake_run(*_args, **_kwargs):
            raise RuntimeError("rev-parse probe")

        with _PatchRun(fake_run):
            try:
                tracked_files.list_tracked(root, "docs/source")
            except RuntimeError as exc:
                _assert(str(exc) == "rev-parse probe", "rev-parse probe swallowed Exception")
            else:
                raise RuntimeError("rev-parse probe caught Exception and fell back")


def _unreadable_manifest_raises() -> None:
    with tempfile.TemporaryDirectory() as tmp:
        root = Path(tmp)
        _write_manifest(root, ["docs/source/from-manifest.rst"])
        path = root / tracked_files.MANIFEST
        path.chmod(0)
        try:
            path.read_bytes()
        except PermissionError:
            pass
        else:
            raise RuntimeError("manifest stayed readable; unreadable-manifest check did not run")

        def fake_run(*_args, **_kwargs):
            raise FileNotFoundError("git")

        try:
            with _PatchRun(fake_run):
                try:
                    tracked_files.list_tracked(root, "docs/source")
                except PermissionError:
                    pass
                else:
                    raise RuntimeError("unreadable manifest returned rows")
        finally:
            path.chmod(0o644)


def _permission_error_launching_git_raises() -> None:
    with tempfile.TemporaryDirectory() as tmp:
        root = Path(tmp)
        _write_manifest(root, ["docs/source/from-manifest.rst"])

        def fake_run(*_args, **_kwargs):
            raise PermissionError("git")

        with _PatchRun(fake_run):
            try:
                tracked_files.list_tracked(root, "docs/source")
            except PermissionError:
                pass
            else:
                raise RuntimeError("PermissionError launching git fell back to the manifest")


def _sibling_prefix_rejected() -> None:
    with tempfile.TemporaryDirectory() as tmp:
        root = Path(tmp)
        _write_manifest(root, ["docs/source/page.rst", "docs/sourcebook/page.rst"])

        def fake_run(*_args, **_kwargs):
            raise FileNotFoundError("git")

        with _PatchRun(fake_run):
            rows = tracked_files.list_tracked(root, "docs/source")
        _assert(rows == ["docs/source/page.rst"], "sibling prefix docs/sourcebook should not match docs/source")


def main() -> int:
    _no_git_with_manifest()
    _inside_unrelated_git_tree()
    _no_git_and_no_manifest()
    _git_toplevel_uses_git_not_manifest()
    _git_ls_files_failure_raises()
    _rev_parse_exit_128_raises()
    _not_git_repo_uses_manifest()
    _localized_not_a_repo_uses_manifest()
    _ls_files_oserror_raises()
    _rev_parse_probe_does_not_catch_exception()
    _unreadable_manifest_raises()
    _permission_error_launching_git_raises()
    _sibling_prefix_rejected()
    print("tracked_files fallback checks ok")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
