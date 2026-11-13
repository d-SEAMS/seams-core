#!/usr/bin/env python3
# Copyright (c) d-SEAMS developers

"""Fetch the deposited mW nucleation dump.

The file is frame-for-frame the trajectory named in
input/traj/nucleation_f1.source.txt. The MD5 is the one recorded there.
The dump is not committed. When it exceeds 1 MB its relative path is
listed in paper/tex/figures/data/science/large-frames.csv.
"""

from __future__ import annotations

import hashlib
import urllib.request
from pathlib import Path

URL = "https://ndownloader.figshare.com/files/20418309"
MD5 = "c28fc209a2854d2ba49e8f9b55accbe8"
EXPECTED_FRAMES = 148


def repo_root() -> Path:
    return Path(__file__).resolve().parents[1]


def md5_file(path: Path) -> str:
    digest = hashlib.md5()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(1 << 20), b""):
            digest.update(chunk)
    return digest.hexdigest()


def record_large(path: Path, root: Path) -> None:
    csv_path = root / "paper" / "tex" / "figures" / "data" / "science" / "large-frames.csv"
    csv_path.parent.mkdir(parents=True, exist_ok=True)
    rel = path.relative_to(root).as_posix()
    kept = []
    if csv_path.is_file():
        for line in csv_path.read_text().splitlines():
            if not line or line.startswith("path ") or line.startswith(rel + " "):
                continue
            kept.append(line)
    lines = ["path bytes", *kept]
    size = path.stat().st_size
    if size > 1_000_000:
        lines.append(f"{rel} {size}")
    csv_path.write_text("\n".join(lines) + "\n")


def main() -> int:
    root = repo_root()
    dest = root / "data" / "science-frames" / "nucleation.lammpstrj"
    dest.parent.mkdir(parents=True, exist_ok=True)
    if not dest.is_file() or md5_file(dest) != MD5:
        print(f"fetch {URL}", flush=True)
        tmp = dest.with_suffix(".lammpstrj.part")
        with urllib.request.urlopen(URL) as response, tmp.open("wb") as out:
            while True:
                chunk = response.read(1 << 20)
                if not chunk:
                    break
                out.write(chunk)
        tmp.replace(dest)
    got = md5_file(dest)
    if got != MD5:
        raise SystemExit(f"md5 {got} != {MD5}")
    frames = 0
    with dest.open() as handle:
        for line in handle:
            if line.startswith("ITEM: TIMESTEP"):
                frames += 1
    if frames != EXPECTED_FRAMES:
        raise SystemExit(f"frames {frames} != {EXPECTED_FRAMES}")
    record_large(dest, root)
    print(f"ok frames={frames} md5={got} path={dest.relative_to(root).as_posix()}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
