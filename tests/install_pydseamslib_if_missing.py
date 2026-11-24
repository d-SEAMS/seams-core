#!/usr/bin/env python3
"""Regression tests for scripts/install_pydseamslib.sh --if-missing."""

from __future__ import annotations

import json
import os
import shlex
import shutil
import subprocess
import sys
import tempfile
from pathlib import Path

sys.dont_write_bytecode = True
sys.path.insert(0, str(Path(__file__).resolve().parent))
from host_tool_env import without_asan_preload  # noqa: E402

ROOT = Path(__file__).resolve().parents[1]
SCRIPT = ROOT / "scripts" / "install_pydseamslib.sh"


def _run(
    cmd: list[str],
    *,
    env: dict[str, str] | None = None,
    cwd: Path | None = None,
) -> subprocess.CompletedProcess[str]:
    """Run a host tool; ENV defaults to os.environ without the sanitizer preload."""
    try:
        return subprocess.run(
            cmd,
            check=True,
            text=True,
            capture_output=True,
            env=without_asan_preload() if env is None else env,
            cwd=cwd,
        )
    except subprocess.CalledProcessError as exc:
        print(f"command failed with exit {exc.returncode}: {cmd}", file=sys.stderr)
        print(f"--- stdout ---\n{exc.stdout}", file=sys.stderr)
        print(f"--- stderr ---\n{exc.stderr}", file=sys.stderr)
        raise


def _write(path: Path, text: str) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(text, encoding="utf-8")


def _make_stub_module(path: Path) -> None:
    _write(path / "dseams.py", "value = 'stub'\n")
    _write(
        path / "pydseams.py",
        "class Frame:\n"
        "    @staticmethod\n"
        "    def classify_chill_plus(*_args, **_kwargs):\n"
        "        return 'stub'\n",
    )


# pip's build isolation installs nothing for an empty requires list, so the
# fake package builds with --no-index and without a network.
FAKE_BACKEND = '''\
import base64
import hashlib
import zipfile
from pathlib import Path

NAME = "pydseamslib"
VERSION = "0.0.0"
MODULES = ("dseams.py", "pydseams.py")
DIST_INFO = f"{NAME}-{VERSION}.dist-info"


def get_requires_for_build_wheel(config_settings=None):
    return []


def _record_line(arcname, data):
    digest = base64.urlsafe_b64encode(hashlib.sha256(data).digest()).rstrip(b"=").decode()
    return f"{arcname},sha256={digest},{len(data)}\\n"


def build_wheel(wheel_directory, config_settings=None, metadata_directory=None):
    wheel_name = f"{NAME}-{VERSION}-py3-none-any.whl"
    files = {name: Path(name).read_bytes() for name in MODULES}
    files[f"{DIST_INFO}/METADATA"] = (
        f"Metadata-Version: 2.1\\nName: {NAME}\\nVersion: {VERSION}\\n"
    ).encode()
    files[f"{DIST_INFO}/WHEEL"] = (
        b"Wheel-Version: 1.0\\nGenerator: fake_backend\\n"
        b"Root-Is-Purelib: true\\nTag: py3-none-any\\n"
    )
    record = "".join(_record_line(arcname, data) for arcname, data in files.items())
    record += f"{DIST_INFO}/RECORD,,\\n"
    files[f"{DIST_INFO}/RECORD"] = record.encode()
    with zipfile.ZipFile(Path(wheel_directory) / wheel_name, "w") as whl:
        for arcname, data in files.items():
            whl.writestr(arcname, data)
    return wheel_name
'''


def _make_fake_pydseamslib(path: Path) -> None:
    _write(
        path / "pyproject.toml",
        "[build-system]\n"
        "requires = []\n"
        'build-backend = "fake_backend"\n'
        'backend-path = ["."]\n',
    )
    _write(path / "fake_backend.py", FAKE_BACKEND)
    _write(path / "dseams.py", "value = 'installed'\n")
    _write(
        path / "pydseams.py",
        "class Frame:\n"
        "    @staticmethod\n"
        "    def classify_chill_plus(*_args, **_kwargs):\n"
        "        return 'installed'\n",
    )
    (path / "subprojects" / "packagefiles").mkdir(parents=True, exist_ok=True)


def _make_fake_git(path: Path, real_git: str) -> None:
    _write(
        path / "git",
        "#!/bin/sh\n"
        "set -eu\n"
        "if [ \"$1\" = \"clone\" ]; then\n"
        "  for last; do :; done\n"
        "  cp -R \"$FAKE_PYDSEAMSLIB_SOURCE\" \"$last\"\n"
        "  : >\"$FAKE_GIT_CLONE_LOG\"\n"
        "  exit 0\n"
        "fi\n"
        f"exec {shlex.quote(real_git)} \"$@\"\n",
    )
    os.chmod(path / "git", 0o755)


def _script_env(
    fake_bin: Path,
    venv: Path,
    fake_src: Path,
    clone_log: Path,
    *,
    pythonpath: Path | None = None,
) -> dict[str, str]:
    env = without_asan_preload()
    env.pop("PYTHONPATH", None)
    if pythonpath is not None:
        env["PYTHONPATH"] = str(pythonpath)
    env["PATH"] = os.pathsep.join([str(fake_bin), str(venv / "bin"), env.get("PATH", "")])
    env["FAKE_PYDSEAMSLIB_SOURCE"] = str(fake_src)
    env["FAKE_GIT_CLONE_LOG"] = str(clone_log)
    return env


def _clean_inspect_env(venv: Path) -> dict[str, str]:
    env = without_asan_preload()
    env.pop("PYTHONPATH", None)
    env["PATH"] = os.pathsep.join([str(venv / "bin"), env.get("PATH", "")])
    return env


def _venv_purelib(venv: Path) -> Path:
    res = _run(
        [
            str(venv / "bin" / "python"),
            "-I",
            "-c",
            "import sysconfig; from pathlib import Path; "
            "print(Path(sysconfig.get_path('purelib')).resolve())",
        ],
        env=_clean_inspect_env(venv),
    )
    return Path(res.stdout.strip())


def _inspect_target(venv: Path) -> dict[str, object]:
    """Inspect the installed pydseams distribution in venv with PYTHONPATH cleared."""
    res = _run(
        [
            str(venv / "bin" / "python"),
            "-I",
            "-c",
            (
                "import json, sysconfig, dseams, pydseams; "
                "from pathlib import Path; "
                "targets = sorted({str(Path(p).resolve()) for p in "
                "(sysconfig.get_path('purelib'), sysconfig.get_path('platlib')) if p}); "
                "p_path = Path(pydseams.__file__).resolve(); "
                "d_path = Path(dseams.__file__).resolve(); "
                "stamp = p_path.parent / '.seams-core-build.json'; "
                "data = json.loads(stamp.read_text(encoding='utf-8')) if stamp.is_file() else None; "
                "print(json.dumps({"
                "'targets': targets, "
                "'pydseams_file': str(p_path), "
                "'dseams_file': str(d_path), "
                "'classify': pydseams.Frame.classify_chill_plus(), "
                "'stamp_path': str(stamp), "
                "'stamp': data"
                "}))"
            ),
        ],
        env=_clean_inspect_env(venv),
    )
    return json.loads(res.stdout)


def _assert_valid_target_install(venv: Path, expected_root: str, expected_tree: str | None) -> Path:
    info = _inspect_target(venv)
    if info["classify"] != "installed":
        raise SystemExit(f"expected installed pydseams in venv, got {info['classify']!r}")
    targets = [Path(t) for t in info["targets"]]
    for key in ("pydseams_file", "dseams_file"):
        p = Path(info[key])
        if not any(p.is_relative_to(t) for t in targets):
            raise SystemExit(f"{key} resolved outside target site-packages: {p} not in {targets}")
    stamp = info["stamp"]
    if not isinstance(stamp, dict):
        raise SystemExit(f"missing install stamp in target site-packages at {info['stamp_path']}")
    if stamp.get("root") != expected_root:
        raise SystemExit(f"install stamp root mismatch: {stamp.get('root')!r} != {expected_root!r}")
    if stamp.get("tree") != expected_tree:
        raise SystemExit(f"install stamp tree mismatch: {stamp.get('tree')!r} != {expected_tree!r}")
    return Path(info["stamp_path"])


def _git_tree(root: Path, rev: str) -> str:
    return _run(["git", "-C", str(root), "rev-parse", f"{rev}^{{tree}}"]).stdout.strip()


def _stage_source_tree(root: Path) -> None:
    """Copy the script into a fresh git repository with two commits whose trees differ.

    The test then does not depend on the source tree under test being a git
    checkout, which an unpacked meson dist tarball is not.
    """
    _write(root / "scripts" / SCRIPT.name, SCRIPT.read_text(encoding="utf-8"))
    _write(root / "subprojects" / "packagefiles" / "README", "staged packagefiles\n")
    env = without_asan_preload(extra={"GIT_CONFIG_GLOBAL": os.devnull, "GIT_CONFIG_NOSYSTEM": "1"})
    ident = ["-c", "user.name=seams-test", "-c", "user.email=seams-test@invalid"]
    _run(["git", "init", "-q", str(root)], env=env)
    for label in ("parent", "head"):
        _write(root / "revision.txt", f"{label}\n")
        _run(["git", "-C", str(root), "add", "-A"], env=env)
        _run(["git", "-C", str(root), *ident, "commit", "-q", "-m", label], env=env)


def _install_if_missing(root: str, env: dict[str, str]) -> subprocess.CompletedProcess[str]:
    script = Path(root) / "scripts" / SCRIPT.name
    return _run(["/bin/sh", str(script), "--if-missing"], env=env, cwd=Path(root))


def test_stub_no_stamp_stamps_target_not_stub(
    tdir: Path,
    fake_bin: Path,
    fake_src: Path,
    expected_root: str,
    expected_tree: str,
) -> Path:
    venv = tdir / "venv-stub"
    stub_dir = tdir / "stub-no-stamp"
    clone_log = tdir / "clone-stub.log"
    _make_stub_module(stub_dir)
    _run([sys.executable, "-m", "venv", str(venv)])

    env = _script_env(fake_bin, venv, fake_src, clone_log, pythonpath=stub_dir)
    _install_if_missing(expected_root, env)
    if not clone_log.is_file():
        raise SystemExit("stub without stamp: --if-missing exited without installing")
    if (stub_dir / ".seams-core-build.json").exists():
        raise SystemExit("stub without stamp: install stamp was written onto PYTHONPATH stub")
    _assert_valid_target_install(venv, expected_root, expected_tree)

    clone_log.unlink()
    _install_if_missing(expected_root, env)
    if clone_log.is_file():
        raise SystemExit("valid target install: --if-missing reinstalled unnecessarily")
    if (stub_dir / ".seams-core-build.json").exists():
        raise SystemExit("valid target install: stamp unexpectedly appeared on PYTHONPATH stub")
    return venv


def test_forged_stub_rejected(
    tdir: Path,
    fake_bin: Path,
    fake_src: Path,
    expected_root: str,
    expected_tree: str,
) -> None:
    venv = tdir / "venv-forged-stub"
    stub_dir = tdir / "stub-forged"
    clone_log = tdir / "clone-forged.log"
    _make_stub_module(stub_dir)
    _write(
        stub_dir / ".seams-core-build.json",
        json.dumps({"root": expected_root, "tree": expected_tree}, sort_keys=True) + "\n",
    )
    _run([sys.executable, "-m", "venv", str(venv)])

    env = _script_env(fake_bin, venv, fake_src, clone_log, pythonpath=stub_dir)
    _install_if_missing(expected_root, env)
    if not clone_log.is_file():
        raise SystemExit("forged PYTHONPATH stub: --if-missing accepted stub without installing")
    _assert_valid_target_install(venv, expected_root, expected_tree)


def test_relative_dot_root_rejected(
    venv: Path,
    tdir: Path,
    fake_bin: Path,
    fake_src: Path,
    expected_root: str,
    expected_tree: str,
) -> None:
    stamp_path = _assert_valid_target_install(venv, expected_root, expected_tree)
    _write(
        stamp_path,
        json.dumps({"root": ".", "tree": expected_tree}, sort_keys=True) + "\n",
    )
    clone_log = tdir / "clone-dot-root.log"
    env = _script_env(fake_bin, venv, fake_src, clone_log)
    _install_if_missing(expected_root, env)
    if not clone_log.is_file():
        raise SystemExit("relative root '.': --if-missing accepted relative root against cwd")
    _assert_valid_target_install(venv, expected_root, expected_tree)


def test_parent_tree_stamp_rejected(
    venv: Path,
    tdir: Path,
    fake_bin: Path,
    fake_src: Path,
    expected_root: str,
    expected_tree: str,
) -> None:
    parent_tree = _git_tree(Path(expected_root), "HEAD~1")
    if parent_tree == expected_tree:
        raise SystemExit("staged source tree: HEAD~1 has the same tree as HEAD")
    stamp_path = _assert_valid_target_install(venv, expected_root, expected_tree)
    _write(
        stamp_path,
        json.dumps({"root": expected_root, "tree": parent_tree}, sort_keys=True) + "\n",
    )
    clone_log = tdir / "clone-parent-tree.log"
    env = _script_env(fake_bin, venv, fake_src, clone_log)
    _install_if_missing(expected_root, env)
    if not clone_log.is_file():
        raise SystemExit("parent tree stamp: --if-missing accepted stale tree without reinstalling")
    _assert_valid_target_install(venv, expected_root, expected_tree)


def test_other_checkout_root_stamp_rejected(
    venv: Path,
    tdir: Path,
    fake_bin: Path,
    fake_src: Path,
    expected_root: str,
    expected_tree: str,
) -> None:
    other_root = (tdir / "other-checkout").resolve()
    other_root.mkdir(parents=True, exist_ok=True)
    stamp_path = _assert_valid_target_install(venv, expected_root, expected_tree)
    _write(
        stamp_path,
        json.dumps({"root": str(other_root), "tree": expected_tree}, sort_keys=True) + "\n",
    )
    clone_log = tdir / "clone-other-root.log"
    env = _script_env(fake_bin, venv, fake_src, clone_log)
    _install_if_missing(expected_root, env)
    if not clone_log.is_file():
        raise SystemExit("other checkout root: --if-missing accepted foreign root without reinstalling")
    _assert_valid_target_install(venv, expected_root, expected_tree)


def test_pth_outside_target_rejected(
    tdir: Path,
    fake_bin: Path,
    fake_src: Path,
    expected_root: str,
    expected_tree: str,
) -> None:
    venv = tdir / "venv-pth"
    outside_dir = tdir / "outside-pth-stub"
    clone_log = tdir / "clone-pth.log"
    _make_stub_module(outside_dir)
    _write(
        outside_dir / ".seams-core-build.json",
        json.dumps({"root": expected_root, "tree": expected_tree}, sort_keys=True) + "\n",
    )
    _run([sys.executable, "-m", "venv", str(venv)])
    _write(_venv_purelib(venv) / "shadow_pydseams.pth", f"{outside_dir}\n")

    env = _script_env(fake_bin, venv, fake_src, clone_log)
    _install_if_missing(expected_root, env)
    if not clone_log.is_file():
        raise SystemExit("outside .pth module: --if-missing accepted module outside target site-packages")
    _assert_valid_target_install(venv, expected_root, expected_tree)


def test_non_git_tree_always_reinstalls(
    staged: Path,
    tdir: Path,
    fake_bin: Path,
    fake_src: Path,
) -> None:
    plain = tdir / "plain-source"
    shutil.copytree(staged, plain, ignore=shutil.ignore_patterns(".git"))
    plain_root = str(plain.resolve(strict=True))
    venv = tdir / "venv-plain"
    _run([sys.executable, "-m", "venv", str(venv)])
    for attempt in ("first", "second"):
        clone_log = tdir / f"clone-plain-{attempt}.log"
        env = _script_env(fake_bin, venv, fake_src, clone_log)
        res = _install_if_missing(plain_root, env)
        if not clone_log.is_file():
            raise SystemExit(f"non-git tree, {attempt} run: --if-missing exited without installing")
        if "is not a git checkout" not in res.stderr:
            raise SystemExit(f"non-git tree, {attempt} run: no notice on stderr: {res.stderr!r}")
        _assert_valid_target_install(venv, plain_root, None)


def test_null_tree_stamp_rejected_in_checkout(
    venv: Path,
    tdir: Path,
    fake_bin: Path,
    fake_src: Path,
    expected_root: str,
    expected_tree: str,
) -> None:
    stamp_path = _assert_valid_target_install(venv, expected_root, expected_tree)
    _write(
        stamp_path,
        json.dumps({"root": expected_root, "tree": None}, sort_keys=True) + "\n",
    )
    clone_log = tdir / "clone-null-tree.log"
    env = _script_env(fake_bin, venv, fake_src, clone_log)
    _install_if_missing(expected_root, env)
    if not clone_log.is_file():
        raise SystemExit("null tree stamp: --if-missing in a git checkout accepted it without reinstalling")
    _assert_valid_target_install(venv, expected_root, expected_tree)


def test_unresolvable_git_head_aborts(
    staged: Path,
    venv: Path,
    tdir: Path,
    fake_bin: Path,
    fake_src: Path,
) -> None:
    git_env = without_asan_preload(extra={"GIT_CONFIG_GLOBAL": os.devnull, "GIT_CONFIG_NOSYSTEM": "1"})
    unborn = tdir / "unborn-source"
    shutil.copytree(staged, unborn, ignore=shutil.ignore_patterns(".git"))
    _run(["git", "init", "-q", str(unborn)], env=git_env)
    corrupt = tdir / "corrupt-source"
    shutil.copytree(staged, corrupt, ignore=shutil.ignore_patterns(".git"))
    _write(corrupt / ".git", f"gitdir: {tdir / 'missing-gitdir'}\n")
    for label, root in (("unborn .git", unborn), ("corrupt .git", corrupt)):
        clone_log = tdir / f"clone-{root.name}.log"
        env = _script_env(fake_bin, venv, fake_src, clone_log)
        res = subprocess.run(
            ["/bin/sh", str(root / "scripts" / SCRIPT.name), "--if-missing"],
            text=True,
            capture_output=True,
            env=env,
            cwd=root,
        )
        if res.returncode == 0:
            raise SystemExit(f"{label}: --if-missing exited 0: {res.stderr!r}")
        if clone_log.is_file():
            raise SystemExit(f"{label}: --if-missing installed despite an unresolvable HEAD")
        if "is not a git checkout" in res.stderr:
            raise SystemExit(f"{label}: treated as a non-git tree: {res.stderr!r}")


def main() -> int:
    real_git = shutil.which("git")
    if real_git is None:
        raise SystemExit("git is not on PATH")

    with tempfile.TemporaryDirectory(prefix="seams-pydseamslib-") as tmp:
        tdir = Path(tmp)
        staged = tdir / "source"
        _stage_source_tree(staged)
        expected_root = str(staged.resolve(strict=True))
        expected_tree = _git_tree(staged, "HEAD")
        fake_src = tdir / "fake-source"
        fake_bin = tdir / "fake-bin"
        _make_fake_pydseamslib(fake_src)
        _make_fake_git(fake_bin, real_git)

        venv = test_stub_no_stamp_stamps_target_not_stub(
            tdir, fake_bin, fake_src, expected_root, expected_tree
        )
        test_forged_stub_rejected(
            tdir, fake_bin, fake_src, expected_root, expected_tree
        )
        test_relative_dot_root_rejected(
            venv, tdir, fake_bin, fake_src, expected_root, expected_tree
        )
        test_parent_tree_stamp_rejected(
            venv, tdir, fake_bin, fake_src, expected_root, expected_tree
        )
        test_other_checkout_root_stamp_rejected(
            venv, tdir, fake_bin, fake_src, expected_root, expected_tree
        )
        test_pth_outside_target_rejected(
            tdir, fake_bin, fake_src, expected_root, expected_tree
        )
        test_non_git_tree_always_reinstalls(staged, tdir, fake_bin, fake_src)
        test_null_tree_stamp_rejected_in_checkout(
            venv, tdir, fake_bin, fake_src, expected_root, expected_tree
        )
        test_unresolvable_git_head_aborts(staged, venv, tdir, fake_bin, fake_src)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
