#!/usr/bin/env python3
"""Build libseams_lang as meson release -O3, or refuse a debug library.

The related-work watch used to time a meson debug -O0 tree. That build
is not a timing result. This script is the asv build command and the
check the harness runs before the first sample.
"""

from __future__ import annotations

import argparse
import json
import shutil
import subprocess
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
BUILD = ROOT / "bbdir-release"


def _option(build_dir: Path, name: str):
    info = Path(build_dir) / "meson-info" / "intro-buildoptions.json"
    if not info.is_file():
        return None
    for row in json.loads(info.read_text(encoding="utf-8")):
        if row.get("name") == name:
            return row.get("value")
    return None


def describe(build_dir: Path) -> str:
    kind = _option(build_dir, "buildtype")
    opt = _option(build_dir, "optimization")
    if kind is None or opt is None:
        raise SystemExit(f"refusing to time {build_dir}: no meson build options")
    return f"{kind} -O{opt}"


def assert_release(build_dir: Path) -> str:
    label = describe(build_dir)
    kind, _, opt = label.partition(" -O")
    if kind != "release" or opt != "3":
        raise SystemExit(
            f"refusing to time a {label} build; the harness requires release -O3")
    commands = Path(build_dir) / "compile_commands.json"
    if commands.is_file():
        text = commands.read_text(encoding="utf-8", errors="replace")
        if " -O0 " in text or " -O0\n" in text or text.rstrip().endswith(" -O0"):
            raise SystemExit("refusing to time a build whose compile commands use -O0")
    return label


def build_dir_of(lib: Path) -> Path | None:
    for parent in [lib.parent, *lib.parents]:
        if (parent / "meson-info" / "intro-buildoptions.json").is_file():
            return parent
        if parent == ROOT:
            break
    return None


def setup_release() -> Path:
    meson = shutil.which("meson")
    if meson is None:
        raise SystemExit("meson is not on PATH; cannot build release -O3")
    if not (BUILD / "meson-info" / "intro-buildoptions.json").is_file():
        subprocess.run(
            [meson, "setup", str(BUILD), str(ROOT),
             "-Dbuildtype=release", "-Doptimization=3",
             "-Dwith_python=false", "-Dwith_tests=true"],
            check=True,
        )
    else:
        kind = _option(BUILD, "buildtype")
        opt = str(_option(BUILD, "optimization"))
        if kind != "release" or opt != "3":
            subprocess.run(
                [meson, "setup", str(BUILD), "--reconfigure",
                 "-Dbuildtype=release", "-Doptimization=3"],
                check=True,
            )
    ninja = shutil.which("ninja") or "ninja"
    subprocess.run(
        [ninja, "-C", str(BUILD), "bindings/libseams_lang.so", "src/libyodaLib.so",
         "tests/bench_rings_h2h", "tests/bench_rings_ih"],
        check=True,
    )
    lib = BUILD / "bindings" / "libseams_lang.so"
    if not lib.is_file():
        raise SystemExit(f"release build did not produce {lib}")
    assert_release(BUILD)
    return lib


def ensure_release(lib: str | None = None) -> tuple[str, str]:
    """Return (library path, 'release -O3'). A debug library is an error."""
    if lib:
        path = Path(lib)
        if not path.is_file():
            raise SystemExit(f"library not found: {path}")
        build = build_dir_of(path.resolve())
        if build is None:
            raise SystemExit(f"refusing to time {path}: not inside a meson build")
        return str(path), assert_release(build)
    built = setup_release()
    return str(built), assert_release(BUILD)


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--lib", default="", help="check this library instead of building")
    args = parser.parse_args()
    path, label = ensure_release(args.lib or None)
    print(f"{label} {path}")


if __name__ == "__main__":
    main()
