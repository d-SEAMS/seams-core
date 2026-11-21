#!/usr/bin/env python3
# Copyright (c) d-SEAMS developers

"""Extract byte-for-byte real-frame fixtures from the fetched trajectories."""

from __future__ import annotations

import hashlib
from pathlib import Path

HERE = Path(__file__).resolve().parent
MARKER = b"ITEM: TIMESTEP\n"
SOURCES = {
    "dump-bcc.lammpstrj": {
        "sha256": "cd77663fa08e7446bc0a9e0367b7f1832244d843e94b2e9b9c99fbb31eb8ef47",
        "indices": tuple(range(30)),
        "output": "bcc-real-30.fixture.lammpstrj",
    },
    "dump-l.lammpstrj": {
        "sha256": "700b7ed79c78f94cbb9d923887e96845275ba4271514371fbe385540b04a5cfc",
        "indices": tuple(range(30)),
        "output": "liquid-real-30.fixture.lammpstrj",
    },
    "mW_cubic.lammpstrj": {
        "sha256": "66d854880bb9431abc9eb25d68b0622d94a0160dd03e1ac7d636343d9830485d",
        "indices": (0, 2, 4, 6, 8, 9, 10),
        "output": "cubic-real-7.fixture.lammpstrj",
    },
}


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


def timestep(frame: bytes) -> int:
    return int(frame.splitlines()[1])


def main() -> int:
    for name, spec in SOURCES.items():
        source = HERE / name
        data = source.read_bytes()
        actual = hashlib.sha256(data).hexdigest()
        if actual != spec["sha256"]:
            raise RuntimeError(f"{name}: expected {spec['sha256']}, got {actual}")
        source_frames = frames(data)
        indices = spec["indices"]
        selected = [source_frames[index] for index in indices]
        output = HERE / spec["output"]
        output.write_bytes(b"".join(selected))
        steps = ",".join(str(timestep(frame)) for frame in selected)
        print(f"{output.name}: frames={len(selected)} timesteps={steps}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
