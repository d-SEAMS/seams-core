#!/usr/bin/env python3
"""Fetch, wire, and inventory the locked d-SEAMS source set."""

from __future__ import annotations

import argparse
import json
import os
from pathlib import Path
import re
import subprocess
from typing import Any


LOCK_SCHEMA = "dseams.ecosystem-lock/v1"
MANIFEST_SCHEMA = "dseams.source-manifest/v1"
REQUIRED_COMPONENTS = {"linkcell", "pydseamslib", "yodastruct"}
FULL_REVISION = re.compile(r"[0-9a-f]{40}")


def load_lock(path: Path) -> dict[str, Any]:
    """Load and validate a source lock with immutable revisions."""
    data = json.loads(path.read_text(encoding="utf-8"))
    if data.get("schema") != LOCK_SCHEMA:
        raise ValueError(f"lock schema must be {LOCK_SCHEMA!r}")
    components = data.get("components")
    if not isinstance(components, dict):
        raise ValueError("lock components must be an object")
    missing = REQUIRED_COMPONENTS - components.keys()
    if missing:
        raise ValueError(f"lock is missing components: {', '.join(sorted(missing))}")
    for name, component in components.items():
        if not isinstance(component, dict):
            raise ValueError(f"component {name!r} must be an object")
        directory = component.get("directory")
        if (
            not isinstance(directory, str)
            or not directory
            or Path(directory).name != directory
            or directory in {".", ".."}
        ):
            raise ValueError(f"component {name!r} has an unsafe directory")
        repository = component.get("repository")
        if not isinstance(repository, str) or not repository:
            raise ValueError(f"component {name!r} has no repository")
        revision = component.get("revision")
        if not isinstance(revision, str) or FULL_REVISION.fullmatch(revision) is None:
            raise ValueError(
                f"component {name!r} revision must be a 40-character lowercase SHA"
            )
    return data


def _git(*args: str, cwd: Path | None = None) -> str:
    process = subprocess.run(
        ["git", *args],
        cwd=cwd,
        capture_output=True,
        text=True,
    )
    if process.returncode != 0:
        rendered = " ".join(("git", *args))
        raise RuntimeError(
            f"{rendered} failed ({process.returncode})\n"
            f"stdout:\n{process.stdout}\nstderr:\n{process.stderr}"
        )
    return process.stdout.strip()


def _component_path(lock: dict[str, Any], root: Path, name: str) -> Path:
    return root / lock["components"][name]["directory"]


def _assert_clean_repository(path: Path, *, include_untracked: bool = True) -> None:
    args = ["status", "--porcelain"]
    if not include_untracked:
        args.append("--untracked-files=no")
    status = _git(*args, cwd=path)
    if status:
        raise RuntimeError(f"locked source has local changes: {path}")


def fetch(lock: dict[str, Any], root: Path) -> None:
    """Materialize every component at its locked revision."""
    root.mkdir(parents=True, exist_ok=True)
    for name in sorted(lock["components"]):
        component = lock["components"][name]
        destination = _component_path(lock, root, name)
        if destination.exists():
            if not (destination / ".git").exists():
                raise FileExistsError(
                    f"source path is not a Git checkout: {destination}"
                )
            _assert_clean_repository(destination)
            configured = _git("remote", "get-url", "origin", cwd=destination)
            if configured != component["repository"]:
                raise RuntimeError(
                    f"origin mismatch for {name}: {configured!r} != "
                    f"{component['repository']!r}"
                )
        else:
            _git(
                "clone",
                "--filter=blob:none",
                "--no-checkout",
                component["repository"],
                str(destination),
            )
        _git(
            "fetch",
            "--depth=1",
            "origin",
            component["revision"],
            cwd=destination,
        )
        _git("checkout", "--detach", component["revision"], cwd=destination)
        actual = _git("rev-parse", "HEAD", cwd=destination)
        if actual != component["revision"]:
            raise RuntimeError(
                f"revision mismatch for {name}: {actual} != {component['revision']}"
            )


def _ensure_symlink(path: Path, target: Path) -> None:
    expected = target.resolve()
    if path.is_symlink():
        if path.resolve() == expected:
            return
        raise FileExistsError(f"refusing to replace symlink {path}")
    if path.exists():
        raise FileExistsError(f"refusing to replace real path {path}")
    path.parent.mkdir(parents=True, exist_ok=True)
    os.symlink(expected, path, target_is_directory=True)


def wire(lock: dict[str, Any], root: Path, core: Path) -> None:
    """Point every frontend and the engine at one locked source graph."""
    core = core.resolve()
    root = root.resolve()
    linkcell = _component_path(lock, root, "linkcell")
    python = _component_path(lock, root, "pydseamslib")
    lua = _component_path(lock, root, "yodastruct")
    for source in (linkcell, python, lua, core):
        if not source.is_dir():
            raise FileNotFoundError(f"source directory is missing: {source}")
    _ensure_symlink(core / "subprojects" / "linkcell", linkcell)
    _ensure_symlink(python / "subprojects" / "seams-core", core)
    _ensure_symlink(lua / "subprojects" / "seams-core", core)


def source_manifest(lock: dict[str, Any], root: Path, core: Path) -> dict[str, Any]:
    """Verify and describe the exact source graph used by a run."""
    components: dict[str, dict[str, str]] = {}
    for name, locked in sorted(lock["components"].items()):
        source = _component_path(lock, root, name)
        _assert_clean_repository(source, include_untracked=False)
        actual = _git("rev-parse", "HEAD", cwd=source)
        if actual != locked["revision"]:
            raise RuntimeError(
                f"revision mismatch for {name}: {actual} != {locked['revision']}"
            )
        components[name] = {
            "repository": locked["repository"],
            "revision": actual,
        }

    core = core.resolve()
    tracked = _git("status", "--porcelain", "--untracked-files=no", cwd=core)
    if tracked:
        raise RuntimeError(f"seams-core has tracked local changes: {core}")
    components["seams-core"] = {
        "repository": _git("remote", "get-url", "origin", cwd=core),
        "revision": _git("rev-parse", "HEAD", cwd=core),
    }
    return {"schema": MANIFEST_SCHEMA, "components": components}


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument(
        "--lock",
        type=Path,
        default=Path("repro/ecosystem-lock.json"),
    )
    subparsers = parser.add_subparsers(dest="command", required=True)

    fetch_parser = subparsers.add_parser("fetch")
    fetch_parser.add_argument("--root", type=Path, required=True)

    wire_parser = subparsers.add_parser("wire")
    wire_parser.add_argument("--root", type=Path, required=True)
    wire_parser.add_argument("--core", type=Path, default=Path.cwd())

    manifest_parser = subparsers.add_parser("manifest")
    manifest_parser.add_argument("--root", type=Path, required=True)
    manifest_parser.add_argument("--core", type=Path, default=Path.cwd())
    manifest_parser.add_argument("--output", type=Path, required=True)

    args = parser.parse_args()
    lock = load_lock(args.lock.resolve())
    if args.command == "fetch":
        fetch(lock, args.root.resolve())
    elif args.command == "wire":
        wire(lock, args.root.resolve(), args.core.resolve())
    elif args.command == "manifest":
        manifest = source_manifest(lock, args.root.resolve(), args.core.resolve())
        output = args.output.resolve()
        output.parent.mkdir(parents=True, exist_ok=True)
        output.write_text(
            json.dumps(manifest, indent=2, sort_keys=True) + "\n",
            encoding="utf-8",
        )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
