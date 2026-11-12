#!/usr/bin/env python3
"""Copy coexistence frames off /tmp every 0.25 ns.

Dumps and restart files land in data/brine-front/, which is gitignored.
A VM restart can resume from the newest dump: positions and velocities
are in the file. This script does not upload anywhere.
"""

from __future__ import annotations

import argparse
import shutil
import time
from pathlib import Path

ROOT = Path(__file__).resolve().parents[2]
DEST = ROOT / "data" / "brine-front"
# Production dumps are every 20 ps. 13 frames is 0.26 ns.
STRIDE = 13
MARKER = b"ITEM: TIMESTEP"


class FrameCursor:
    """Count new frames without reading the dump from the start.

    A full reread every poll holds the growing trajectory in memory and
    competes with the integrator. Only the bytes appended since the last
    poll are scanned. The marker can straddle that boundary.
    """

    def __init__(self) -> None:
        self.offset = 0
        self.count = 0
        self.inode = None
        self.carry = b""

    def update(self, path: Path) -> int:
        if not path.is_file():
            return 0
        st = path.stat()
        if self.inode != st.st_ino or st.st_size < self.offset:
            self.offset = 0
            self.count = 0
            self.inode = st.st_ino
            self.carry = b""
        if st.st_size > self.offset:
            with path.open("rb") as handle:
                handle.seek(self.offset)
                chunk = handle.read()
            self.offset = st.st_size
            # carry is shorter than the marker, so it holds no complete hit.
            # A marker split across the previous poll is completed here.
            data = self.carry + chunk
            self.count += data.count(MARKER)
            self.carry = data[-(len(MARKER) - 1) :]
        return self.count


def snapshot(work: Path, dest: Path, count: int) -> None:
    dest.mkdir(parents=True, exist_ok=True)
    stamp = dest / f"dump-{count:04d}.lammpstrj"
    shutil.copy2(work / "dump.lammpstrj", stamp)
    for name in ("restart.a", "restart.b", "meta.txt", "in.coexist"):
        src = work / name
        if src.is_file():
            shutil.copy2(src, dest / name)
    (dest / "latest.txt").write_text(f"{count} {stamp.name}\n")
    print(f"checkpoint {work.name} frames {count} -> {stamp}", flush=True)


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--src", type=Path, default=Path("/tmp/brine-science-v2"))
    parser.add_argument("--poll", type=float, default=30.0)
    args = parser.parse_args()
    seen: dict[str, int] = {}
    cursors: dict[str, FrameCursor] = {}
    while True:
        if not args.src.is_dir():
            time.sleep(args.poll)
            continue
        for dump in sorted(args.src.glob("*/dump.lammpstrj")):
            tag = dump.parent.name
            cursor = cursors.setdefault(tag, FrameCursor())
            count = cursor.update(dump)
            previous = seen.get(tag, 0)
            # Stride is measured from the last copy. Advancing the mark
            # on every extra frame never reaches STRIDE.
            if count >= previous + STRIDE or (count > 0 and previous == 0):
                snapshot(dump.parent, DEST / tag, count)
                seen[tag] = count
        time.sleep(args.poll)


if __name__ == "__main__":
    main()
