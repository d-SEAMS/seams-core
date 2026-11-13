#!/usr/bin/env python3
"""GenIce sI, sII and sH cage census.

sI is 2 of 5^12 and 6 of 5^12 6^2. sII is 16 and 8. The bundled sH
cell is four unit cells, so its census is four times 3, 2 and 1.
"""

import os
import re
import subprocess
import sys
from pathlib import Path


def cages(seams, path, signature, libdir, edge=None):
    cmd = [
        seams,
        "cages",
        str(path),
        "--type",
        "1",
        "--cutoff",
        "3.5",
        "--graph",
        "cutoff",
        "--signature",
        signature,
    ]
    if edge is not None:
        cmd.extend(["--cage-edge", edge])
    env = os.environ.copy()
    env["LD_LIBRARY_PATH"] = str(libdir) + os.pathsep + env.get("LD_LIBRARY_PATH", "")
    text = subprocess.check_output(cmd, env=env, text=True)
    found = {}
    for line in text.splitlines():
        match = re.search(r"signature (\S+) cages (\d+)", line)
        if match:
            found[match.group(1)] = int(match.group(2))
    return found


def main():
    root = Path(sys.argv[1])
    seams = sys.argv[2]
    libdir = sys.argv[3]
    traj = root / "input" / "traj"
    si = cages(seams, traj / "genice_sI.lammpstrj", "512+51262", libdir)
    sii = cages(seams, traj / "genice_sII.lammpstrj", "512+51264", libdir)
    sh = cages(seams, traj / "genice_sH.lammpstrj", "512+sh+51268", libdir, edge="none")
    if si != {"5:12": 2, "5:12,6:2": 6}:
        raise SystemExit(f"sI {si}")
    if sii != {"5:12": 16, "5:12,6:4": 8}:
        raise SystemExit(f"sII {sii}")
    unit = (3, 2, 1)
    got = (sh.get("5:12"), sh.get("4:3,5:6,6:3"), sh.get("5:12,6:8"))
    if got != tuple(4 * n for n in unit):
        raise SystemExit(f"sH {sh}")
    print("clathrate ok", si, sii, got)


if __name__ == "__main__":
    main()
