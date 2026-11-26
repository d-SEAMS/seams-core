#!/usr/bin/env python3
# Copyright (c) d-SEAMS developers
"""cargo-from-source uses the manifest directory only for a vendored config."""

from __future__ import annotations

import os
import stat
import subprocess
import sys
import tempfile
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
WRAPPER = ROOT / "scripts" / "cargo-from-source"


def _assert(condition: bool, message: str) -> None:
    if not condition:
        raise RuntimeError(message)


def _write_project(path: Path, vendored: bool) -> None:
    (path / "src").mkdir(parents=True)
    (path / "Cargo.toml").write_text(
        '[package]\nname = "probe"\nversion = "0.0.0"\nedition = "2021"\n',
        encoding="utf-8",
    )
    if vendored:
        cargo = path / ".cargo"
        cargo.mkdir()
        (cargo / "config.toml").write_text(
            '[source.crates-io]\nreplace-with = "vendored-sources"\n'
            '\n[source.vendored-sources]\ndirectory = "vendor"\n',
            encoding="utf-8",
        )


def _fake_cargo(directory: Path) -> None:
    directory.mkdir()
    program = directory / "cargo"
    program.write_text(
        "#!/bin/sh\nprintf '%s\\n' \"$PWD\" > \"$CARGO_FROM_SOURCE_RECORD\"\n",
        encoding="utf-8",
    )
    program.chmod(program.stat().st_mode | stat.S_IEXEC)


def _recorded(base: Path, project: Path, args: list[str]) -> str:
    record = base / "cwd.txt"
    if record.exists():
        record.unlink()
    other = base / "other"
    other.mkdir(exist_ok=True)
    env = {key: value for key, value in os.environ.items() if key != "CARGO_FROM_SOURCE_RECORD"}
    env["PATH"] = str(base / "bin")
    env["CARGO_FROM_SOURCE_RECORD"] = str(record)
    subprocess.run(
        [sys.executable, str(WRAPPER), *args],
        cwd=other,
        env=env,
        check=True,
        capture_output=True,
        text=True,
    )
    return record.read_text(encoding="utf-8").strip()


def _vendored_config_selects_the_manifest_directory() -> None:
    with tempfile.TemporaryDirectory() as tmp:
        base = Path(tmp)
        project = base / "proj"
        _write_project(project, vendored=True)
        _fake_cargo(base / "bin")
        manifest = str(project / "Cargo.toml")
        got = _recorded(base, project, ["metadata", "--manifest-path", manifest])
        _assert(got == str(project), "vendored cargo cwd was {}".format(got))
        got = _recorded(base, project, ["metadata", "--manifest-path=" + manifest])
        _assert(got == str(project), "equals-form cwd was {}".format(got))


def _plain_config_keeps_the_working_directory() -> None:
    with tempfile.TemporaryDirectory() as tmp:
        base = Path(tmp)
        project = base / "proj"
        _write_project(project, vendored=False)
        _fake_cargo(base / "bin")
        got = _recorded(
            base,
            project,
            ["metadata", "--manifest-path", str(project / "Cargo.toml")],
        )
        _assert(got == str(base / "other"), "plain cargo cwd was {}".format(got))


def main() -> int:
    _vendored_config_selects_the_manifest_directory()
    _plain_config_keeps_the_working_directory()
    print("cargo-from-source checks ok")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
