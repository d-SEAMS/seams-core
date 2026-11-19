#!/usr/bin/env python3
# Copyright (c) d-SEAMS developers

"""Require and exercise the complete, checksum-pinned trajectory inputs."""

from __future__ import annotations

import hashlib
import subprocess
import sys
from pathlib import Path

EXPECTED = {
    "dump-bcc.lammpstrj": (
        "cd77663fa08e7446bc0a9e0367b7f1832244d843e94b2e9b9c99fbb31eb8ef47",
        201,
        432,
        1_000_000,
    ),
    "dump-l.lammpstrj": (
        "700b7ed79c78f94cbb9d923887e96845275ba4271514371fbe385540b04a5cfc",
        201,
        512,
        1_000_000,
    ),
    "mW_cubic.lammpstrj": (
        "66d854880bb9431abc9eb25d68b0622d94a0160dd03e1ac7d636343d9830485d",
        11,
        4096,
        1000,
    ),
}
MARKER = b"ITEM: TIMESTEP\n"


def frames(data: bytes) -> list[bytes]:
    starts: list[int] = []
    offset = 0
    while True:
        offset = data.find(MARKER, offset)
        if offset < 0:
            break
        starts.append(offset)
        offset += len(MARKER)
    ends = starts[1:] + [len(data)]
    return [data[start:end] for start, end in zip(starts, ends)]


def frame_metadata(frame: bytes) -> tuple[int, int, int]:
    lines = frame.splitlines()
    step = int(lines[1])
    atoms = int(lines[3])
    atom_id = int(lines[9 + atoms - 1].split()[0])
    return step, atoms, atom_id


def main() -> int:
    seams = Path(sys.argv[1]).resolve()
    root = Path(sys.argv[2]).resolve()
    trajectory_dir = root / "input" / "traj"
    for name, (
        expected_hash,
        expected_frames,
        expected_atoms,
        expected_last_step,
    ) in EXPECTED.items():
        path = trajectory_dir / name
        data = path.read_bytes()
        actual_hash = hashlib.sha256(data).hexdigest()
        if actual_hash != expected_hash:
            raise RuntimeError(f"{name}: expected {expected_hash}, got {actual_hash}")
        snapshots = frames(data)
        if len(snapshots) != expected_frames:
            raise RuntimeError(
                f"{name}: expected {expected_frames} frames, got {len(snapshots)}"
            )
        first = frame_metadata(snapshots[0])
        last = frame_metadata(snapshots[-1])
        if first[1] != expected_atoms or last[1] != expected_atoms:
            raise RuntimeError(f"{name}: atom counts {first[1]} and {last[1]}")
        if first[0] != 0 or last[0] != expected_last_step:
            raise RuntimeError(f"{name}: timesteps {first[0]}..{last[0]}")
        if name == "mW_cubic.lammpstrj" and last[2] != 4096:
            raise RuntimeError(f"{name}: last-frame final atom ID is {last[2]}")
        if first[0] == last[0] or snapshots[0] == snapshots[-1]:
            raise RuntimeError(f"{name}: first and last frames are not distinct")
        print(
            f"{name}: sha256={actual_hash} frames={len(snapshots)} "
            f"steps={first[0]}..{last[0]} atoms={expected_atoms} "
            f"last_atom_id={last[2]}"
        )

    subprocess.run(
        [
            sys.executable,
            str(root / "tests" / "cli_carry.py"),
            str(seams),
            str(trajectory_dir),
            "full",
        ],
        check=True,
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
