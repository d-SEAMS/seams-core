# Copyright (c) d-SEAMS developers
"""The head-to-head timer runs on one small frame and writes a CSV."""

import subprocess
import sys
import tempfile
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]


def main() -> int:
    if len(sys.argv) != 2:
        print("usage: test_h2h_harness.py BENCH", file=sys.stderr)
        return 2
    bench = Path(sys.argv[1])
    traj = ROOT / "input" / "traj" / "mW_cubic.lammpstrj"
    if not bench.is_file() or not traj.is_file():
        print("missing harness or frame", file=sys.stderr)
        return 1
    with tempfile.TemporaryDirectory(prefix="dseams-h2h-") as work:
        csv = Path(work) / "neighbours.csv"
        proc = subprocess.run(
            [str(bench), "neighbours", str(traj), "3.5", "6", str(csv)],
            check=False,
            capture_output=True,
            text=True,
        )
        if proc.returncode != 0:
            sys.stderr.write(proc.stdout + proc.stderr)
            return 1
        lines = csv.read_text(encoding="utf-8").splitlines()
        if len(lines) < 2 or not lines[0].strip():
            print("empty neighbour csv", file=sys.stderr)
            return 1
    summary = ROOT / "scripts" / "h2h_cpc_summary.py"
    if not summary.is_file():
        print("missing summary script", file=sys.stderr)
        return 1
    print(f"h2h neighbours {len(lines) - 1}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
