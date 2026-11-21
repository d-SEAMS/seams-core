#!/usr/bin/env python3
"""Regression test for scripts/install_pydseamslib.sh --if-missing."""

from __future__ import annotations

import os
import subprocess
import tempfile
from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]
SCRIPT = ROOT / "scripts" / "install_pydseamslib.sh"


def _run(cmd: list[str], *, env: dict[str, str] | None = None) -> subprocess.CompletedProcess[str]:
    return subprocess.run(cmd, check=True, text=True, capture_output=True, env=env)


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


def main() -> int:
    with tempfile.TemporaryDirectory(prefix="seams-pydseamslib-") as tmp:
        tdir = Path(tmp)
        stub_dir = tdir / "stub"
        fake_src = tdir / "fake-source"
        fake_bin = tdir / "fake-bin"
        clone_log = tdir / "clone.log"
        venv = tdir / "venv"
        _make_stub_module(stub_dir)
        _make_fake_pydseamslib(fake_src)
        _make_fake_git(fake_bin)
        _run(["python3", "-m", "venv", str(venv)])

        env = dict(os.environ)
        env["PATH"] = os.pathsep.join([str(fake_bin), str(venv / "bin"), env.get("PATH", "")])
        env["PYTHONPATH"] = str(stub_dir)
        env["FAKE_PYDSEAMSLIB_SOURCE"] = str(fake_src)
        env["FAKE_GIT_CLONE_LOG"] = str(clone_log)

        _run(["/bin/sh", str(SCRIPT), "--if-missing"], env=env)
        if not clone_log.is_file():
            raise SystemExit("--if-missing exited without installing from this checkout")

        check = _run(
            [
                str(venv / "bin" / "python"),
                "-c",
                (
                    "import json, pydseams; from pathlib import Path; "
                    "stamp = Path(pydseams.__file__).resolve().parent / '.seams-core-build.json'; "
                    "print(json.loads(stamp.read_text(encoding='utf-8'))['root'])"
                ),
            ],
            env=env,
        )
        got = check.stdout.strip()
        if got != str(ROOT):
            raise SystemExit(f"install stamp root mismatch: {got!r} != {str(ROOT)!r}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
