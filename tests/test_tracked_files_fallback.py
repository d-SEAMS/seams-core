#!/usr/bin/env python3
# Copyright (c) d-SEAMS developers
"""Checks for tracked-files fallback behavior outside git work trees."""

from __future__ import annotations

import importlib.util
import os
import shlex
import shutil
import socket
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



def _make_special_git_marker(root: Path, kind: str) -> None:
    marker = root / ".git"
    if kind == "fifo":
        os.mkfifo(marker)
    elif kind == "socket":
        with socket.socket(socket.AF_UNIX, socket.SOCK_STREAM) as sock:
            sock.bind(str(marker))
    else:
        raise RuntimeError("unknown special marker kind {}".format(kind))


_MANIFEST_ROW = ["docs/source/from-manifest.rst"]
_DE_NOT_REPO = (
    "Schwerwiegend: Kein Git-Repository (oder irgendeines der Elternverzeichnisse): .git"
)
_FR_NOT_REPO = (
    "fatal\u00a0: ni ceci ni aucun de ses répertoires parents n'est un dépôt git\u00a0: .git"
)


def _git_outcome_cases():
    """Git result to the rows or the error the caller must see.

    The GIT_DIR scrub cases are a different contract and stay beside this table.
    """
    rows = _MANIFEST_ROW
    plain = "plain"
    return [
        {
            "name": "no-git manifest fallback",
            "layout": plain,
            "manifest": ["docs/source/alpha.rst", "docs/guide.md", "README.md"],
            "git": "missing",
            "expect": ("rows", ["docs/source/alpha.rst"]),
        },
        {
            "name": "sdist under unrelated git tree",
            "layout": "unrelated",
            "manifest": rows,
            "git": "outer-top",
            "top": "str",
            "expect": ("rows", ["docs/source/from-manifest.rst"]),
        },
        {
            "name": "git and manifest missing",
            "layout": plain,
            "manifest": None,
            "git": "missing",
            "expect": ("error", RuntimeError, "missing"),
        },
        {
            "name": "git toplevel prefers git rows",
            "layout": plain,
            "manifest": rows,
            "git": "toplevel-git",
            "expect": ("rows", ["docs/source/from-git.rst"]),
        },
        {
            "name": "ls-files failure",
            "layout": plain,
            "manifest": rows,
            "git": "toplevel-ls-fail",
            "expect": ("error", RuntimeError, "ls-files failed"),
        },
        {
            "name": "rev-parse exit 128 with .git directory",
            "layout": "gitdir",
            "manifest": rows,
            "git": "code",
            "code": 128,
            "stderr": "fatal: detected dubious ownership in repository",
            "expect": ("error", RuntimeError, "rev-parse failed"),
        },
        {
            "name": "rev-parse exit 128 with .git file",
            "layout": "gitfile",
            "manifest": rows,
            "git": "code",
            "code": 128,
            "stderr": "fatal: not a git repository: /nonexistent/worktree",
            "expect": ("error", RuntimeError, "rev-parse failed"),
        },
        {
            "name": "rev-parse exit 1 without marker",
            "layout": plain,
            "manifest": rows,
            "git": "code",
            "code": 1,
            "stderr": "git failed with code 1",
            "expect": ("error", RuntimeError, "rev-parse failed"),
        },
        {
            "name": "rev-parse exit 129 without marker",
            "layout": plain,
            "manifest": rows,
            "git": "code",
            "code": 129,
            "stderr": "git failed with code 129",
            "expect": ("error", RuntimeError, "rev-parse failed"),
        },
        {
            "name": "not a repository uses the manifest",
            "layout": plain,
            "manifest": rows,
            "git": "code",
            "code": 128,
            "stderr": "fatal: not a git repository (or any of the parent directories): .git",
            "expect": ("rows", ["docs/source/from-manifest.rst"]),
        },
        {
            "name": "German not-a-repository uses the manifest",
            "layout": plain,
            "manifest": rows,
            "git": "code",
            "code": 128,
            "stderr": _DE_NOT_REPO,
            "stderr_lacks": "not a git repository",
            "expect": ("rows", ["docs/source/from-manifest.rst"]),
        },
        {
            "name": "French not-a-repository uses the manifest",
            "layout": plain,
            "manifest": rows,
            "git": "code",
            "code": 128,
            "stderr": _FR_NOT_REPO,
            "stderr_lacks": "not a git repository",
            "expect": ("rows", ["docs/source/from-manifest.rst"]),
        },
        {
            "name": "ls-files OSError",
            "layout": plain,
            "manifest": rows,
            "git": "toplevel-oserror",
            "expect": ("error", (OSError, RuntimeError), "ls-files"),
        },
        {
            "name": "rev-parse probe is not caught",
            "layout": plain,
            "manifest": rows,
            "git": "probe",
            "expect": ("exact", RuntimeError, "rev-parse probe"),
        },
        {
            "name": "unreadable manifest",
            "layout": plain,
            "manifest": rows,
            "git": "missing",
            "expect": ("unreadable",),
        },
        {
            "name": "PermissionError launching git",
            "layout": plain,
            "manifest": rows,
            "git": "permission",
            "expect": ("error", PermissionError, None),
        },
        {
            "name": "missing git with .git directory",
            "layout": "gitdir",
            "manifest": rows,
            "git": "missing-errno",
            "expect": ("error", FileNotFoundError, None),
        },
        {
            "name": "missing git with .git file",
            "layout": "gitfile",
            "manifest": rows,
            "git": "missing-errno",
            "expect": ("error", FileNotFoundError, None),
        },
        {
            "name": "dangling .git symlink",
            "layout": "dangling",
            "manifest": rows,
            "git": "code",
            "code": 128,
            "stderr": "fatal: not a git repository",
            "expect": ("error", RuntimeError, "rev-parse failed"),
        },
        {
            "name": "dangling marker inside an outer repo",
            "layout": "outer-dangling",
            "manifest": rows,
            "git": "outer-top",
            "expect": ("marker",),
        },
        {
            "name": "empty .git inside an outer repo",
            "layout": "outer-empty",
            "manifest": rows,
            "git": "outer-top",
            "expect": ("marker",),
        },
        {
            "name": "symlink marker inside an outer repo",
            "layout": "outer-symlink",
            "manifest": rows,
            "git": "outer-top",
            "expect": ("marker",),
        },
        {
            "name": "fifo marker inside an outer repo",
            "layout": "outer-fifo",
            "manifest": rows,
            "git": "outer-top",
            "expect": ("marker",),
        },
        {
            "name": "socket marker inside an outer repo",
            "layout": "outer-socket",
            "manifest": rows,
            "git": "outer-top",
            "expect": ("marker",),
        },
        {
            "name": "fifo marker rev-parse 128",
            "layout": "standalone",
            "marker": "fifo",
            "manifest": rows,
            "git": "code",
            "code": 128,
            "stderr": "fatal: not a git repository",
            "expect": ("error", RuntimeError, "rev-parse failed"),
        },
        {
            "name": "socket marker rev-parse 128",
            "layout": "standalone",
            "marker": "socket",
            "manifest": rows,
            "git": "code",
            "code": 128,
            "stderr": "fatal: not a git repository",
            "expect": ("error", RuntimeError, "rev-parse failed"),
        },
        {
            "name": "fifo marker missing git",
            "layout": "standalone",
            "marker": "fifo",
            "manifest": rows,
            "git": "missing-errno",
            "expect": ("error", FileNotFoundError, None),
        },
        {
            "name": "socket marker missing git",
            "layout": "standalone",
            "marker": "socket",
            "manifest": rows,
            "git": "missing-errno",
            "expect": ("error", FileNotFoundError, None),
        },
        {
            "name": "sibling prefix",
            "layout": plain,
            "manifest": ["docs/source/page.rst", "docs/sourcebook/page.rst"],
            "git": "missing",
            "expect": ("rows", ["docs/source/page.rst"]),
        },
    ]


def _prepare_outcome(tmp: str, case: dict):
    base = Path(tmp)
    layout = case["layout"]
    outer = None
    if layout == "unrelated":
        outer = base / "outer"
        outer.mkdir()
        (outer / ".git").mkdir()
        (outer / ".git" / "HEAD").write_text("ref: refs/heads/main\n", encoding="utf-8")
        (outer / ".git" / "objects").mkdir()
        (outer / ".git" / "refs" / "heads").mkdir(parents=True)
        root = outer / "sdist"
        root.mkdir()
    elif layout.startswith("outer-"):
        outer = base / "outer"
        root = outer / "sdist"
        root.mkdir(parents=True)
        kind = layout.split("-", 1)[1]
        if kind == "dangling":
            (root / ".git").symlink_to(root / "missing-git-dir", target_is_directory=True)
        elif kind == "empty":
            (root / ".git").mkdir()
        elif kind == "symlink":
            plain = root / "plain-dir"
            plain.mkdir()
            (root / ".git").symlink_to(plain, target_is_directory=True)
        elif kind in ("fifo", "socket"):
            _make_special_git_marker(root, kind)
        else:
            raise RuntimeError("unknown layout {}".format(layout))
    elif layout == "standalone":
        root = base / "standalone"
        root.mkdir()
        _make_special_git_marker(root, case["marker"])
    else:
        root = base
        if layout == "gitdir":
            (root / ".git").mkdir()
        elif layout == "gitfile":
            (root / ".git").write_text("gitdir: /nonexistent/worktree\n", encoding="utf-8")
        elif layout == "dangling":
            (root / ".git").symlink_to(root / "missing-git-dir", target_is_directory=True)
        elif layout != "plain":
            raise RuntimeError("unknown layout {}".format(layout))
    if case["manifest"] is not None:
        _write_manifest(root, case["manifest"])
    return root, outer


def _fake_git(case: dict, root: Path, outer):
    def fake_run(args, check=False, capture_output=False, text=False, env=None, **_kwargs):
        kind = case["git"]
        if kind == "missing":
            raise FileNotFoundError("git")
        if kind == "missing-errno":
            raise FileNotFoundError(2, "No such file or directory", "git")
        if kind == "permission":
            raise PermissionError("git")
        if kind == "probe":
            raise RuntimeError("rev-parse probe")
        if kind == "code":
            err = case["stderr"]
            stderr = err if text else err.encode()
            return _cp(args, code=case["code"], stderr=stderr)
        if kind == "toplevel-git":
            if args[3:] == ["rev-parse", "--show-toplevel"]:
                return _cp(args, stdout=str(root) + ("\n" if text else ""))
            if args[3] == "ls-files":
                return _cp(args, stdout=b"docs/source/from-git.rst\0")
            raise RuntimeError("unexpected call {}".format(args))
        if kind == "toplevel-ls-fail":
            if args[3:] == ["rev-parse", "--show-toplevel"]:
                return _cp(args, stdout=str(root) + ("\n" if text else ""))
            if args[3] == "ls-files":
                return _cp(args, code=128, stderr=b"fatal: index file corrupt")
            raise RuntimeError("unexpected call {}".format(args))
        if kind == "toplevel-oserror":
            if args[3:] == ["rev-parse", "--show-toplevel"]:
                return _cp(args, stdout=str(root) + ("\n" if text else ""))
            if args[3] == "ls-files":
                raise OSError("ls-files")
            raise RuntimeError("unexpected call {}".format(args))
        if kind == "outer-top":
            if args[3:] == ["rev-parse", "--show-toplevel"]:
                shown = outer if case.get("top") == "str" else outer.resolve()
                return _cp(args, stdout=str(shown) + ("\n" if text else ""))
            raise RuntimeError("unexpected call {}".format(args))
        raise RuntimeError("unknown git kind {}".format(kind))

    return fake_run


def _one_git_outcome(case: dict) -> None:
    lacks = case.get("stderr_lacks")
    if lacks is not None:
        _assert(
            lacks not in case["stderr"].lower(),
            "localized stderr accidentally contains the English phrase",
        )
    with tempfile.TemporaryDirectory() as tmp:
        root, outer = _prepare_outcome(tmp, case)
        fake = _fake_git(case, root, outer)
        expect = case["expect"]
        kind = expect[0]
        if kind == "rows":
            with _PatchRun(fake):
                got = tracked_files.list_tracked(root, "docs/source")
            _assert(got == expect[1], "{}: {}".format(case["name"], got))
            return
        if kind == "unreadable":
            path = root / tracked_files.MANIFEST
            err = PermissionError(13, "Permission denied", str(path))
            with _PatchRun(fake), _PatchReadBytes(path, err) as patch_read:
                try:
                    tracked_files.list_tracked(root, "docs/source")
                except PermissionError:
                    pass
                else:
                    raise RuntimeError("unreadable manifest returned rows")
                _assert(patch_read.calls == 1, "unreadable-manifest check did not invoke read_bytes")
            return
        if kind == "marker":
            with _PatchRun(fake):
                try:
                    tracked_files.list_tracked(root, "docs/source")
                except RuntimeError as exc:
                    text = str(exc)
                    _assert(
                        "has a .git marker" in text,
                        "{} did not name the marker: {}".format(case["name"], text),
                    )
                    _assert(
                        str(outer.resolve()) in text,
                        "{} did not name the resolved toplevel: {}".format(case["name"], text),
                    )
                else:
                    raise RuntimeError("{} silently used the manifest".format(case["name"]))
            return
        exc_type, needle = expect[1], expect[2]
        with _PatchRun(fake):
            try:
                tracked_files.list_tracked(root, "docs/source")
            except exc_type as exc:
                if kind == "exact":
                    _assert(str(exc) == needle, "{}: {}".format(case["name"], exc))
                elif needle is not None:
                    _assert(needle in str(exc), "{}: {}".format(case["name"], exc))
            else:
                raise RuntimeError("{} did not raise".format(case["name"]))


def _git_result_outcomes() -> None:
    for case in _git_outcome_cases():
        _one_git_outcome(case)


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


def _git_repo_env_scrubs_every_repo_selecting_var() -> None:
    repo_vars = (
        ("GIT_DIR", "/tmp/outer/.git"),
        ("GIT_WORK_TREE", "/tmp/outer"),
        ("GIT_INDEX_FILE", "/tmp/outer/index"),
        ("GIT_COMMON_DIR", "/tmp/outer/.git"),
        ("GIT_OBJECT_DIRECTORY", "/tmp/outer/.git/objects"),
        ("GIT_ALTERNATE_OBJECT_DIRECTORIES", "/tmp/outer/.git/objects"),
        ("GIT_NAMESPACE", "outer-ns"),
        ("GIT_PREFIX", "subdir/"),
        ("GIT_IMPLICIT_WORK_TREE", "1"),
    )
    config_vars = (
        ("GIT_CONFIG_COUNT", "1"),
        ("GIT_CONFIG_KEY_0", "safe.directory"),
        ("GIT_CONFIG_VALUE_0", "*"),
        ("GIT_CONFIG_PARAMETERS", "'safe.directory=*'"),
        ("GIT_CONFIG_GLOBAL", "/tmp/outer/gitconfig"),
    )
    with tempfile.TemporaryDirectory() as tmp:
        root = Path(tmp)
        (root / ".git").mkdir()
        _write_manifest(root, ["docs/source/from-manifest.rst"])
        for var_name, var_value in repo_vars:
            touched = [var_name] + [name for name, _ in config_vars]
            saved = {name: os.environ.get(name) for name in touched}
            os.environ[var_name] = var_value
            for cfg_name, cfg_value in config_vars:
                os.environ[cfg_name] = cfg_value
            seen_subcommands: list[str] = []

            def fake_run(
                args,
                check=False,
                capture_output=False,
                text=False,
                env=None,
                var_name=var_name,
                seen_subcommands=seen_subcommands,
                **_kwargs,
            ):
                _assert(env is not None, "git subprocess env was not set for {}".format(var_name))
                _assert(
                    var_name not in env,
                    "{} leaked into git subprocess env".format(var_name),
                )
                for cfg_name, cfg_value in config_vars:
                    _assert(
                        env.get(cfg_name) == cfg_value,
                        "{} was not preserved for git when scrubbing {}".format(cfg_name, var_name),
                    )
                if args[3:] == ["rev-parse", "--show-toplevel"]:
                    seen_subcommands.append("rev-parse")
                    out = str(root) + ("\n" if text else "")
                    return _cp(args, stdout=out)
                if args[3] == "ls-files":
                    seen_subcommands.append("ls-files")
                    return _cp(args, stdout=b"docs/source/from-git.rst\0")
                raise RuntimeError("unexpected call {}".format(args))

            try:
                with _PatchRun(fake_run):
                    rows = tracked_files.list_tracked(root, "docs/source")
            finally:
                for name, old in saved.items():
                    if old is None:
                        os.environ.pop(name, None)
                    else:
                        os.environ[name] = old
            _assert(
                seen_subcommands == ["rev-parse", "ls-files"],
                "expected rev-parse and ls-files calls for {}, saw {}".format(
                    var_name, seen_subcommands
                ),
            )
            _assert(
                rows == ["docs/source/from-git.rst"],
                "git rows were not returned when scrubbing {}".format(var_name),
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
        (root / "docs" / "clean.md").write_text("Clean documentation page.\n", encoding="utf-8")
        (root / "docs" / "new-page.md").write_text("The new page is seamless.\n", encoding="utf-8")
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
    # Same fake-git pattern as the broken-git check. There is no real git
    # process, so nothing can stage into an inherited GIT_DIR or rewrite an
    # inherited GIT_INDEX_FILE, and the check still runs when git is not on PATH.
    with tempfile.TemporaryDirectory() as tmp:
        base = Path(tmp).resolve()
        outer = base / "outer"
        root = outer / "sdist"
        (root / "scripts").mkdir(parents=True)
        (root / "docs").mkdir(parents=True)
        shutil.copy2(_tracked_files_path(), root / "scripts" / "tracked_files.py")
        shutil.copy2(ROOT / "scripts" / "run_doc_check.py", root / "scripts" / "run_doc_check.py")
        shutil.copy2(ROOT / "docs" / "check_prose.py", root / "docs" / "check_prose.py")
        (root / "docs" / "clean.md").write_text("The clean page is seamless.\n", encoding="utf-8")
        (root / "docs" / "new-page.md").write_text("Supplementary documentation page.\n", encoding="utf-8")
        _write_manifest(root, ["docs/clean.md"])

        index_file = base / "inherited-index"
        index_bytes = b"sentinel-index-do-not-touch\n"
        index_file.write_bytes(index_bytes)

        shimbin = base / "shimbin"
        shimbin.mkdir()
        helper = base / "fake_git.py"
        repo_vars = (
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
        helper.write_text(
            "import os\n"
            "import sys\n"
            "repo_vars = {repo_vars!r}\n"
            "leaked = [name for name in repo_vars if name in os.environ]\n"
            "args = sys.argv[1:]\n"
            "sub = args[2:] if len(args) >= 2 and args[0] == '-C' else args\n"
            "if leaked:\n"
            "    if sub == ['rev-parse', '--show-toplevel']:\n"
            "        sys.stdout.write({root!r} + '\\n')\n"
            "        raise SystemExit(0)\n"
            "    if sub and sub[0] == 'ls-files':\n"
            "        raise SystemExit(0)\n"
            "    sys.stderr.write('repo git var leaked: ' + ','.join(leaked) + '\\n')\n"
            "    raise SystemExit(2)\n"
            "if sub == ['rev-parse', '--show-toplevel']:\n"
            "    sys.stdout.write({outer!r} + '\\n')\n"
            "    raise SystemExit(0)\n"
            "if sub and sub[0] == 'ls-files':\n"
            "    sys.stdout.buffer.write(b'docs/clean.md\\0docs/new-page.md\\0')\n"
            "    raise SystemExit(0)\n"
            "sys.stderr.write('unexpected git call: ' + repr(args) + '\\n')\n"
            "raise SystemExit(2)\n".format(repo_vars=repo_vars, root=str(root), outer=str(outer)),
            encoding="utf-8",
        )
        shim_git = shimbin / "git"
        shim_git.write_text(
            "#!/bin/sh\nexec {} {} \"$@\"\n".format(
                shlex.quote(sys.executable), shlex.quote(str(helper))
            ),
            encoding="utf-8",
        )
        shim_git.chmod(0o755)

        cmd = [
            sys.executable,
            str(root / "scripts" / "run_doc_check.py"),
            str(root / "docs" / "check_prose.py"),
        ]
        env = dict(os.environ)
        env["PATH"] = str(shimbin)
        env["GIT_DIR"] = str(outer / ".git")
        env["GIT_WORK_TREE"] = str(outer)
        env["GIT_INDEX_FILE"] = str(index_file)
        env["LC_ALL"] = "C.UTF-8"
        env["LANG"] = "C.UTF-8"
        done = subprocess.run(cmd, cwd=root, env=env, capture_output=True, text=True)
        _assert(index_file.read_bytes() == index_bytes, "inherited GIT_INDEX_FILE was modified")
        _assert(
            done.returncode == 1 and "1 prose hits" in done.stderr,
            "docs_prose with inherited outer GIT_DIR did not fail on staged prose violation: rc={} stdout={!r} stderr={!r}".format(
                done.returncode, done.stdout, done.stderr
            ),
        )



def main() -> int:
    _git_result_outcomes()
    _git_subprocess_env_scrubs_git_vars()
    _git_config_reaches_subprocess_and_index_file_is_scrubbed()
    _git_repo_env_scrubs_every_repo_selecting_var()
    _docs_prose_broken_git_does_not_pass_stale_manifest()
    _docs_prose_inherited_git_dir_outer_repo_fails()
    print("tracked_files fallback checks ok")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
