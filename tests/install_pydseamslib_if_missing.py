#!/usr/bin/env python3
"""Regression tests for scripts/install_pydseamslib.sh --if-missing."""

from __future__ import annotations

import json
import os
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


def _make_fake_pydseamslib(path: Path) -> None:
    _write(
        path / "pyproject.toml",
        "[build-system]\n"
        'requires = ["setuptools>=68"]\n'
        'build-backend = "setuptools.build_meta"\n',
    )
    _write(
        path / "setup.cfg",
        "[metadata]\n"
        "name = pydseamslib\n"
        "version = 0.0.0\n\n"
        "[options]\n"
        "py_modules = dseams, pydseams\n",
    )
    _write(path / "dseams.py", "value = 'installed'\n")
    _write(
        path / "pydseams.py",
        "class Frame:\n"
        "    @staticmethod\n"
        "    def classify_chill_plus(*_args, **_kwargs):\n"
        "        return 'installed'\n",
    )
    (path / "subprojects" / "packagefiles").mkdir(parents=True, exist_ok=True)


def _make_fake_git(path: Path) -> None:
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
        "exec /usr/bin/git \"$@\"\n",
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


def _assert_valid_target_install(venv: Path, expected_root: str, expected_tree: str) -> Path:
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


def _git_tree(rev: str) -> str | None:
    res = subprocess.run(
        ["git", "-C", str(ROOT), "rev-parse", f"{rev}^{{tree}}"],
        text=True,
        capture_output=True,
        env=without_asan_preload(),
    )
    if res.returncode != 0:
        return None
    return res.stdout.strip()


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
    _run(["/bin/sh", str(SCRIPT), "--if-missing"], env=env, cwd=ROOT)
    if not clone_log.is_file():
        raise SystemExit("stub without stamp: --if-missing exited without installing")
    if (stub_dir / ".seams-core-build.json").exists():
        raise SystemExit("stub without stamp: install stamp was written onto PYTHONPATH stub")
    _assert_valid_target_install(venv, expected_root, expected_tree)

    clone_log.unlink()
    _run(["/bin/sh", str(SCRIPT), "--if-missing"], env=env, cwd=ROOT)
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
    _run(["/bin/sh", str(SCRIPT), "--if-missing"], env=env, cwd=ROOT)
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
    _run(["/bin/sh", str(SCRIPT), "--if-missing"], env=env, cwd=ROOT)
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
    parent_tree = _git_tree("HEAD~1")
    if not parent_tree or parent_tree == expected_tree:
        parent_tree = "0" * 40 if expected_tree != "0" * 40 else "1" * 40
    stamp_path = _assert_valid_target_install(venv, expected_root, expected_tree)
    _write(
        stamp_path,
        json.dumps({"root": expected_root, "tree": parent_tree}, sort_keys=True) + "\n",
    )
    clone_log = tdir / "clone-parent-tree.log"
    env = _script_env(fake_bin, venv, fake_src, clone_log)
    _run(["/bin/sh", str(SCRIPT), "--if-missing"], env=env, cwd=ROOT)
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
    _run(["/bin/sh", str(SCRIPT), "--if-missing"], env=env, cwd=ROOT)
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
    _run(["/bin/sh", str(SCRIPT), "--if-missing"], env=env, cwd=ROOT)
    if not clone_log.is_file():
        raise SystemExit("outside .pth module: --if-missing accepted module outside target site-packages")
    _assert_valid_target_install(venv, expected_root, expected_tree)


def main() -> int:
    expected_root = str(ROOT.resolve(strict=True))
    expected_tree = _git_tree("HEAD")
    if not expected_tree:
        raise SystemExit("could not resolve HEAD^{tree}")

    with tempfile.TemporaryDirectory(prefix="seams-pydseamslib-") as tmp:
        tdir = Path(tmp)
        fake_src = tdir / "fake-source"
        fake_bin = tdir / "fake-bin"
        _make_fake_pydseamslib(fake_src)
        _make_fake_git(fake_bin)

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
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
