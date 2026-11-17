# Copyright (c) d-SEAMS developers
"""seams ions prints the capped and broken ring counts."""

import math
import subprocess
import sys
import tempfile
from pathlib import Path


def hexagon_dump(path: Path) -> None:
    lines = [
        "ITEM: TIMESTEP",
        "0",
        "ITEM: NUMBER OF ATOMS",
        "7",
        "ITEM: BOX BOUNDS pp pp pp",
        "0 30",
        "0 30",
        "0 30",
        "ITEM: ATOMS id type x y z",
        "1 2 15 15 15",
    ]
    for i in range(6):
        angle = i * math.pi / 3.0
        x = 15.0 + 2.8 * math.cos(angle)
        y = 15.0 + 2.8 * math.sin(angle)
        lines.append(f"{i + 2} 1 {x:.6f} {y:.6f} 15")
    path.write_text("\n".join(lines) + "\n", encoding="utf-8")


def field(text: str, name: str) -> int:
    parts = text.split()
    for i, token in enumerate(parts):
        if token == name and i + 1 < len(parts):
            return int(parts[i + 1])
    raise ValueError(name)


def main() -> int:
    if len(sys.argv) != 2:
        print("usage: test_ion_shell_rings.py SEAMS", file=sys.stderr)
        return 2
    seams = sys.argv[1]
    with tempfile.TemporaryDirectory(prefix="dseams-ions-") as work:
        dump = Path(work) / "hex.lammpstrj"
        hexagon_dump(dump)
        proc = subprocess.run(
            [
                seams,
                "-k",
                "2",
                "--cutoff",
                "3.5",
                "--type",
                "1",
                "ions",
                str(dump),
                "--ion-types",
                "2",
                "--ion-cutoff",
                "3.5",
            ],
            check=False,
            capture_output=True,
            text=True,
        )
    if proc.returncode != 0:
        sys.stderr.write(proc.stdout + proc.stderr)
        return 1
    try:
        capped = field(proc.stdout, "capped")
        broken = field(proc.stdout, "broken")
    except ValueError as exc:
        sys.stderr.write(proc.stdout + proc.stderr)
        print(exc, file=sys.stderr)
        return 1
    if capped != 1 or broken < 1:
        print(f"capped {capped} broken {broken}", file=sys.stderr)
        print(proc.stdout, file=sys.stderr)
        return 1
    print(f"capped {capped} broken {broken}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
