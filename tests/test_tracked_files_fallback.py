#!/usr/bin/env python3
# Copyright (c) d-SEAMS developers
"""Checks for tracked-files fallback behavior outside git work trees."""

from __future__ import annotations

import importlib.util
import os
import shlex
import shutil
import subprocess
import sys
import tempfile
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]


def _tracked_files_path() -> Path:
    impl = os.environ.get("TRACKED_FILES_IMPL")
    return Path(impl) if impl else (ROOT / "scripts" / "tracked_files.py")


def _load_tracked_files():
    path = _tracked_files_path()
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


class _PatchReadBytes:
    def __init__(self, target: Path, exc: BaseException):
        self._target = target.resolve()
        self._exc = exc
        self._old = None
        self.calls = 0

    def __enter__(self):
        old = tracked_files.Path.read_bytes
        self._old = old
        target = self._target
        exc = self._exc
        patch = self

        def wrapped(path_self: Path):
            if path_self.resolve() == target:
                patch.calls += 1
                raise exc
            return old(path_self)

        tracked_files.Path.read_bytes = wrapped
        return self

    def __exit__(self, _exc_type, _exc, _tb):
        tracked_files.Path.read_bytes = self._old


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

        def fake_run(args, check=False, capture_output=False, text=False, env=None, **_kwargs):
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

        def fake_run(args, check=False, capture_output=False, text=False, env=None, **_kwargs):
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

        def fake_run(args, check=False, capture_output=False, text=False, env=None, **_kwargs):
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

        def fake_run(args, check=False, capture_output=False, text=False, env=None, **_kwargs):
            err = "fatal: detected dubious ownership in repository"
            return _cp(args, code=128, stderr=err)

        with _PatchRun(fake_run):
            try:
                tracked_files.list_tracked(root, "docs/source")
            except RuntimeError as exc:
                _assert("rev-parse failed" in str(exc), "rev-parse exit-128 was not raised")
            else:
                raise RuntimeError("expected RuntimeError on rev-parse exit-128")


def _rev_parse_exit_128_with_git_file_raises() -> None:
    with tempfile.TemporaryDirectory() as tmp:
        root = Path(tmp)
        (root / ".git").write_text("gitdir: /nonexistent/worktree\n", encoding="utf-8")
        _write_manifest(root, ["docs/source/from-manifest.rst"])

        def fake_run(args, check=False, capture_output=False, text=False, env=None, **_kwargs):
            err = "fatal: not a git repository: /nonexistent/worktree"
            return _cp(args, code=128, stderr=err)

        with _PatchRun(fake_run):
            try:
                tracked_files.list_tracked(root, "docs/source")
            except RuntimeError as exc:
                _assert("rev-parse failed" in str(exc), "rev-parse exit-128 with .git file was not raised")
            else:
                raise RuntimeError("expected RuntimeError on rev-parse exit-128 with .git file")


def _rev_parse_non_128_without_marker_raises() -> None:
    for code in (1, 129):
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            _write_manifest(root, ["docs/source/from-manifest.rst"])

            def fake_run(args, check=False, capture_output=False, text=False, env=None, code=code, **_kwargs):
                return _cp(args, code=code, stderr="git failed with code {}".format(code))

            with _PatchRun(fake_run):
                try:
                    tracked_files.list_tracked(root, "docs/source")
                except RuntimeError as exc:
                    _assert(
                        "rev-parse failed" in str(exc),
                        "rev-parse exit-{} was not raised".format(code),
                    )
                else:
                    raise RuntimeError(
                        "expected RuntimeError on rev-parse exit-{} without .git marker".format(code)
                    )


def _not_git_repo_uses_manifest() -> None:
    with tempfile.TemporaryDirectory() as tmp:
        root = Path(tmp)
        _write_manifest(root, ["docs/source/from-manifest.rst"])

        def fake_run(args, check=False, capture_output=False, text=False, env=None, **_kwargs):
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

            def fake_run(args, check=False, capture_output=False, text=False, env=None, err=err, **_kwargs):
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

        def fake_run(args, check=False, capture_output=False, text=False, env=None, **_kwargs):
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

        def fake_run(*_args, **_kwargs):
            raise FileNotFoundError("git")

        err = PermissionError(13, "Permission denied", str(path))
        with _PatchRun(fake_run), _PatchReadBytes(path, err) as patch_read:
            try:
                tracked_files.list_tracked(root, "docs/source")
            except PermissionError:
                pass
            else:
                raise RuntimeError("unreadable manifest returned rows")
            _assert(patch_read.calls == 1, "unreadable-manifest check did not invoke read_bytes")


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


def _git_marker_missing_git_raises() -> None:
    for marker_kind in ("dir", "file"):
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            if marker_kind == "dir":
                (root / ".git").mkdir()
            else:
                (root / ".git").write_text("gitdir: /nonexistent/worktree\n", encoding="utf-8")
            _write_manifest(root, ["docs/source/from-manifest.rst"])

            def fake_run(*_args, **_kwargs):
                raise FileNotFoundError(2, "No such file or directory", "git")

            with _PatchRun(fake_run):
                try:
                    tracked_files.list_tracked(root, "docs/source")
                except FileNotFoundError:
                    pass
                else:
                    raise RuntimeError(
                        "FileNotFoundError with .git {} marker fell back to the manifest".format(
                            marker_kind
                        )
                    )


def _git_subprocess_env_scrubs_git_vars() -> None:
    with tempfile.TemporaryDirectory() as tmp:
        root = Path(tmp)
        (root / ".git").mkdir()
        _write_manifest(root, ["docs/source/from-manifest.rst"])
        old_git_dir = os.environ.get("GIT_DIR")
        old_git_work_tree = os.environ.get("GIT_WORK_TREE")
        old_git_index_file = os.environ.get("GIT_INDEX_FILE")
        os.environ["GIT_DIR"] = "/tmp/outer/.git"
        os.environ["GIT_WORK_TREE"] = "/tmp/outer"
        os.environ["GIT_INDEX_FILE"] = "/tmp/outer/index"

        def fake_run(args, check=False, capture_output=False, text=False, env=None, **_kwargs):
            _assert(env is not None, "git subprocess env was not set")
            _assert("GIT_DIR" not in env, "GIT_DIR leaked into git subprocess env")
            _assert("GIT_WORK_TREE" not in env, "GIT_WORK_TREE leaked into git subprocess env")
            _assert("GIT_INDEX_FILE" not in env, "GIT_INDEX_FILE leaked into git subprocess env")
            if args[3:] == ["rev-parse", "--show-toplevel"]:
                out = str(root) + ("\n" if text else "")
                return _cp(args, stdout=out)
            if args[3] == "ls-files":
                return _cp(args, stdout=b"docs/source/from-git.rst\0")
            raise RuntimeError("unexpected call {}".format(args))

        try:
            with _PatchRun(fake_run):
                rows = tracked_files.list_tracked(root, "docs/source")
        finally:
            if old_git_dir is None:
                del os.environ["GIT_DIR"]
            else:
                os.environ["GIT_DIR"] = old_git_dir
            if old_git_work_tree is None:
                del os.environ["GIT_WORK_TREE"]
            else:
                os.environ["GIT_WORK_TREE"] = old_git_work_tree
            if old_git_index_file is None:
                del os.environ["GIT_INDEX_FILE"]
            else:
                os.environ["GIT_INDEX_FILE"] = old_git_index_file
        _assert(rows == ["docs/source/from-git.rst"], "git rows were not returned with scrubbed env")


def _git_config_reaches_subprocess_and_index_file_is_scrubbed() -> None:
    with tempfile.TemporaryDirectory() as tmp:
        root = Path(tmp)
        (root / ".git").mkdir()
        _write_manifest(root, ["docs/source/from-manifest.rst"])
        saved = {
            key: os.environ.get(key)
            for key in (
                "GIT_DIR",
                "GIT_WORK_TREE",
                "GIT_INDEX_FILE",
                "GIT_CONFIG_COUNT",
                "GIT_CONFIG_KEY_0",
                "GIT_CONFIG_VALUE_0",
                "GIT_CONFIG_PARAMETERS",
                "GIT_CONFIG_GLOBAL",
            )
        }
        os.environ["GIT_DIR"] = "/tmp/outer/.git"
        os.environ["GIT_WORK_TREE"] = "/tmp/outer"
        os.environ["GIT_INDEX_FILE"] = "/tmp/outer/index"
        os.environ["GIT_CONFIG_COUNT"] = "1"
        os.environ["GIT_CONFIG_KEY_0"] = "safe.directory"
        os.environ["GIT_CONFIG_VALUE_0"] = "*"
        os.environ["GIT_CONFIG_PARAMETERS"] = "'safe.directory=*'"
        os.environ["GIT_CONFIG_GLOBAL"] = "/tmp/outer/gitconfig"

        def fake_run(args, check=False, capture_output=False, text=False, env=None, **_kwargs):
            _assert(env is not None, "git subprocess env was not set")
            for leaked in ("GIT_DIR", "GIT_WORK_TREE", "GIT_INDEX_FILE"):
                _assert(leaked not in env, "{} leaked into git subprocess env".format(leaked))
            _assert(env.get("GIT_CONFIG_COUNT") == "1", "GIT_CONFIG_COUNT was not preserved for git")
            _assert(
                env.get("GIT_CONFIG_PARAMETERS") == "'safe.directory=*'",
                "GIT_CONFIG_PARAMETERS was not preserved for git",
            )
            _assert(env.get("GIT_CONFIG_KEY_0") == "safe.directory", "GIT_CONFIG_KEY_0 was not preserved for git")
            _assert(env.get("GIT_CONFIG_VALUE_0") == "*", "GIT_CONFIG_VALUE_0 was not preserved for git")
            _assert(
                env.get("GIT_CONFIG_GLOBAL") == "/tmp/outer/gitconfig",
                "GIT_CONFIG_GLOBAL was not preserved for git",
            )
            if args[3:] == ["rev-parse", "--show-toplevel"]:
                out = str(root) + ("\n" if text else "")
                return _cp(args, stdout=out)
            if args[3] == "ls-files":
                return _cp(args, stdout=b"docs/source/from-git.rst\0")
            raise RuntimeError("unexpected call {}".format(args))

        try:
            with _PatchRun(fake_run):
                rows = tracked_files.list_tracked(root, "docs/source")
        finally:
            for key, old in saved.items():
                if old is None:
                    os.environ.pop(key, None)
                else:
                    os.environ[key] = old
        _assert(
            rows == ["docs/source/from-git.rst"],
            "git rows were not returned when config env was preserved",
        )


def _rev_parse_exit_128_with_dangling_git_symlink_raises() -> None:
    with tempfile.TemporaryDirectory() as tmp:
        root = Path(tmp)
        (root / ".git").symlink_to(root / "missing-git-dir", target_is_directory=True)
        _write_manifest(root, ["docs/source/from-manifest.rst"])

        def fake_run(args, check=False, capture_output=False, text=False, env=None, **_kwargs):
            return _cp(args, code=128, stderr="fatal: not a git repository")

        with _PatchRun(fake_run):
            try:
                tracked_files.list_tracked(root, "docs/source")
            except RuntimeError as exc:
                _assert("rev-parse failed" in str(exc), "dangling .git symlink did not raise rev-parse failure")
            else:
                raise RuntimeError("expected RuntimeError on rev-parse exit-128 with dangling .git symlink")


def _nonrepo_marker_inside_outer_repo_raises() -> None:
    for variant in ("dangling", "emptydir", "symlink-nonrepo"):
        with tempfile.TemporaryDirectory() as tmp:
            outer = Path(tmp) / "outer"
            root = outer / "sdist"
            root.mkdir(parents=True)
            if variant == "dangling":
                (root / ".git").symlink_to(root / "missing-git-dir", target_is_directory=True)
            elif variant == "emptydir":
                (root / ".git").mkdir()
            else:
                plain = root / "plain-dir"
                plain.mkdir()
                (root / ".git").symlink_to(plain, target_is_directory=True)
            _write_manifest(root, ["docs/source/from-manifest.rst"])

            def fake_run(args, check=False, capture_output=False, text=False, env=None, **_kwargs):
                if args[3:] == ["rev-parse", "--show-toplevel"]:
                    out = str(outer.resolve()) + ("\n" if text else "")
                    return _cp(args, stdout=out)
                raise RuntimeError("unexpected call {}".format(args))

            with _PatchRun(fake_run):
                try:
                    tracked_files.list_tracked(root, "docs/source")
                except RuntimeError as exc:
                    _assert(
                        "has a .git marker" in str(exc),
                        "{} marker inside an outer repo did not name the marker: {}".format(variant, exc),
                    )
                    _assert(
                        str(outer.resolve()) in str(exc),
                        "{} marker inside an outer repo did not name the resolved toplevel: {}".format(
                            variant, exc
                        ),
                    )
                else:
                    raise RuntimeError(
                        "{} .git marker inside an outer repo silently used the manifest".format(variant)
                    )


def _docs_prose_broken_git_does_not_pass_stale_manifest() -> None:
    with tempfile.TemporaryDirectory() as tmp:
        base = Path(tmp).resolve()
        root = base / "repo"
        (root / ".git").mkdir(parents=True)
        (root / "scripts").mkdir(parents=True)
        (root / "docs").mkdir(parents=True)
        shutil.copy2(_tracked_files_path(), root / "scripts" / "tracked_files.py")
        shutil.copy2(ROOT / "scripts" / "run_doc_check.py", root / "scripts" / "run_doc_check.py")
        shutil.copy2(ROOT / "docs" / "check_prose.py", root / "docs" / "check_prose.py")
        (root / "docs" / "clean.md").write_text("The clean page is seamless.\n", encoding="utf-8")
        (root / "docs" / "new-page.md").write_text("Supplementary documentation page.\n", encoding="utf-8")
        _write_manifest(root, ["docs/clean.md"])

        goodbin = base / "goodbin"
        goodbin.mkdir()
        helper = base / "fake_git.py"
        helper.write_text(
            "import sys\n"
            "args = sys.argv[1:]\n"
            "sub = args[2:] if len(args) >= 2 and args[0] == '-C' else args\n"
            "if sub == ['rev-parse', '--show-toplevel']:\n"
            "    sys.stdout.write({root!r} + '\\n')\n"
            "    raise SystemExit(0)\n"
            "if sub and sub[0] == 'ls-files':\n"
            "    sys.stdout.buffer.write(b'docs/clean.md\\0docs/new-page.md\\0')\n"
            "    raise SystemExit(0)\n"
            "sys.stderr.write('unexpected git call: ' + repr(args) + '\\n')\n"
            "raise SystemExit(2)\n".format(root=str(root)),
            encoding="utf-8",
        )
        good_git = goodbin / "git"
        good_git.write_text(
            "#!/bin/sh\nexec {} {} \"$@\"\n".format(
                shlex.quote(sys.executable), shlex.quote(str(helper))
            ),
            encoding="utf-8",
        )
        good_git.chmod(0o755)

        badbin = base / "badbin"
        badbin.mkdir()
        bad_git = badbin / "git"
        bad_git.write_text("#!/nonexistent/interpreter\n", encoding="utf-8")
        bad_git.chmod(0o755)

        nobin = base / "nobin"
        nobin.mkdir()

        cmd = [
            sys.executable,
            str(root / "scripts" / "run_doc_check.py"),
            str(root / "docs" / "check_prose.py"),
        ]
        good_env = dict(os.environ, PATH=str(goodbin))
        good_done = subprocess.run(cmd, cwd=root, env=good_env, capture_output=True, text=True)
        _assert(
            good_done.returncode == 1 and "1 prose hits" in good_done.stderr,
            "docs_prose with working git should catch the new rule-breaking page: rc={} stderr={!r}".format(
                good_done.returncode, good_done.stderr
            ),
        )

        for label, bindir in (("broken-shebang", badbin), ("missing-git", nobin)):
            env = dict(os.environ, PATH=str(bindir))
            done = subprocess.run(cmd, cwd=root, env=env, capture_output=True, text=True)
            _assert(
                done.returncode != 0 and "prose check ok" not in done.stdout,
                "docs_prose with {} silently passed on a stale manifest: rc={} stdout={!r} stderr={!r}".format(
                    label, done.returncode, done.stdout, done.stderr
                ),
            )
            _assert(
                "FileNotFoundError" in done.stderr,
                "docs_prose with {} did not fail with FileNotFoundError: stderr={!r}".format(
                    label, done.stderr
                ),
            )


def _docs_prose_inherited_git_dir_outer_repo_fails() -> None:
    with tempfile.TemporaryDirectory() as tmp:
        base = Path(tmp).resolve()
        outer = base / "outer"
        outer.mkdir()
        subprocess.run(["git", "init", str(outer)], check=True, capture_output=True)

        root = outer / "sdist"
        (root / "scripts").mkdir(parents=True)
        (root / "docs").mkdir(parents=True)
        shutil.copy2(_tracked_files_path(), root / "scripts" / "tracked_files.py")
        shutil.copy2(ROOT / "scripts" / "run_doc_check.py", root / "scripts" / "run_doc_check.py")
        shutil.copy2(ROOT / "docs" / "check_prose.py", root / "docs" / "check_prose.py")
        (root / "docs" / "clean.md").write_text("The clean page is seamless.\n", encoding="utf-8")
        (root / "docs" / "new-page.md").write_text("Supplementary documentation page.\n", encoding="utf-8")
        _write_manifest(root, ["docs/clean.md"])

        subprocess.run(["git", "-C", str(outer), "add", "-A"], check=True, capture_output=True)
        cmd = [
            sys.executable,
            str(root / "scripts" / "run_doc_check.py"),
            str(root / "docs" / "check_prose.py"),
        ]
        env = dict(os.environ, GIT_DIR=str(outer / ".git"), LC_ALL="C.UTF-8", LANG="C.UTF-8")
        done = subprocess.run(cmd, cwd=root, env=env, capture_output=True, text=True)
        _assert(
            done.returncode == 1 and "1 prose hits" in done.stderr,
            "docs_prose with inherited outer GIT_DIR did not fail on staged prose violation: rc={} stdout={!r} stderr={!r}".format(
                done.returncode, done.stdout, done.stderr
            ),
        )


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
    _rev_parse_exit_128_with_git_file_raises()
    _rev_parse_non_128_without_marker_raises()
    _not_git_repo_uses_manifest()
    _localized_not_a_repo_uses_manifest()
    _ls_files_oserror_raises()
    _rev_parse_probe_does_not_catch_exception()
    _unreadable_manifest_raises()
    _permission_error_launching_git_raises()
    _git_marker_missing_git_raises()
    _git_subprocess_env_scrubs_git_vars()
    _git_config_reaches_subprocess_and_index_file_is_scrubbed()
    _rev_parse_exit_128_with_dangling_git_symlink_raises()
    _nonrepo_marker_inside_outer_repo_raises()
    _docs_prose_broken_git_does_not_pass_stale_manifest()
    _docs_prose_inherited_git_dir_outer_repo_fails()
    _sibling_prefix_rejected()
    print("tracked_files fallback checks ok")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
