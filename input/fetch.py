#!/usr/bin/env python3
# Copyright (c) d-SEAMS developers

"""Fetch the large public trajectory inputs and verify their checksums."""

from __future__ import annotations

import argparse
import hashlib
import shutil
import socket
import urllib.error
import urllib.request
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
DATASET = "d-SEAMS/seams-core-inputs"
REVISION = "163ca2d8dfcb7ef21d88aea1c67381d841f63641"
BASE_URL = f"https://huggingface.co/datasets/{DATASET}/resolve/{REVISION}/traj"
# A hung connection must not stall a docs or test build indefinitely.
SOCKET_TIMEOUT_S = 30
TRAJECTORIES = {
    "dump-bcc.lammpstrj": "cd77663fa08e7446bc0a9e0367b7f1832244d843e94b2e9b9c99fbb31eb8ef47",
    "dump-l.lammpstrj": "700b7ed79c78f94cbb9d923887e96845275ba4271514371fbe385540b04a5cfc",
    "mW_cubic.lammpstrj": "66d854880bb9431abc9eb25d68b0622d94a0160dd03e1ac7d636343d9830485d",
}


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as stream:
        for chunk in iter(lambda: stream.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def fetch(name: str, destination: Path) -> Path:
    expected = TRAJECTORIES[name]
    output = destination / name
    if output.is_file() and sha256(output) == expected:
        print(f"verified {output}")
        return output

    destination.mkdir(parents=True, exist_ok=True)
    partial = output.with_name(output.name + ".partial")
    partial.unlink(missing_ok=True)
    url = f"{BASE_URL}/{name}"
    try:
        request = urllib.request.Request(
            url, headers={"User-Agent": "seams-core-input-fetch/1"}
        )
        try:
            with urllib.request.urlopen(
                request, timeout=SOCKET_TIMEOUT_S
            ) as source, partial.open("wb") as target:
                shutil.copyfileobj(source, target)
        except (urllib.error.URLError, TimeoutError, socket.timeout) as exc:
            raise RuntimeError(
                f"{name}: could not fetch {url} (timeout {SOCKET_TIMEOUT_S}s): {exc}"
            ) from exc
        actual = sha256(partial)
        if actual != expected:
            raise RuntimeError(
                f"{name}: expected sha256 {expected}, downloaded {actual}"
            )
        partial.replace(output)
    except Exception:
        partial.unlink(missing_ok=True)
        raise
    print(f"fetched and verified {output}")
    return output


def ensure(name: str, override: str | None, destination: Path | None = None) -> Path:
    """The local path for a pinned trajectory, honoring an explicit override.

    With no override this fetches and verifies the trajectory (a no-op
    network-wise when the local copy already matches the pinned digest),
    so callers always get a real, checksummed file or a loud failure.
    """
    if override:
        return Path(override)
    return fetch(name, destination if destination is not None else ROOT / "input" / "traj")


def main() -> int:
    parser = argparse.ArgumentParser(
        description=(
            f"fetch public trajectories from {DATASET} at revision {REVISION}"
        )
    )
    parser.add_argument(
        "names",
        nargs="*",
        choices=tuple(TRAJECTORIES),
        metavar="TRAJECTORY",
    )
    parser.add_argument(
        "--output-dir",
        type=Path,
        default=ROOT / "input" / "traj",
        help="destination directory (default: input/traj)",
    )
    args = parser.parse_args()
    for name in args.names or TRAJECTORIES:
        fetch(name, args.output_dir)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
