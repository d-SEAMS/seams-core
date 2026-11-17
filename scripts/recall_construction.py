#!/usr/bin/env python3
"""Score d-SEAMS and CHILL+ against construction labels.

The binary builds an ideal cubic block and a spaced liquid block. The
truth is that construction. A cubic site is interior only when the cell
holds its 4 first and 12 second diamond neighbours. The all-site score
and the interior score are both written. Per-site rows, including the
distance to the open face, go to the sites table.
"""

from __future__ import annotations

import argparse
import csv
import os
import subprocess
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]


def child_env(binary: Path) -> dict[str, str]:
    env = os.environ.copy()
    env["OMP_NUM_THREADS"] = "1"
    env["OMP_DYNAMIC"] = "FALSE"
    origin = binary.resolve().parent
    libs = [
        Path("/cvmfs/software.eessi.io/versions/2026.06/compat/linux/x86_64/lib64"),
        origin.parent / "src",
        origin.parent / "chemfiles-prefix" / "lib",
        origin.parent / "subprojects" / "readcon-core",
        origin.parent / "subprojects" / "readcon-db",
        Path("/lib/x86_64-linux-gnu"),
    ]
    parts = [str(path) for path in libs if path.is_dir()]
    old = env.get("LD_LIBRARY_PATH", "")
    if old:
        parts.append(old)
    if parts:
        env["LD_LIBRARY_PATH"] = ":".join(parts)
    return env


def head_commit() -> str:
    done = subprocess.run(
        ["git", "rev-parse", "HEAD"], cwd=ROOT, text=True, capture_output=True, check=True
    )
    return done.stdout.strip()


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--bench", type=Path, required=True)
    parser.add_argument("--out", type=Path, required=True)
    parser.add_argument("--sites", type=Path, required=True)
    parser.add_argument("--commit", default="")
    args = parser.parse_args()
    if not args.bench.is_file():
        raise SystemExit(f"missing {args.bench}")
    args.sites.parent.mkdir(parents=True, exist_ok=True)
    done = subprocess.run(
        [str(args.bench), str(args.sites)],
        cwd=ROOT,
        env=child_env(args.bench),
        text=True,
        capture_output=True,
    )
    if done.returncode != 0:
        sys.stderr.write(done.stderr)
        raise SystemExit(f"recall exited {done.returncode}")
    commit = args.commit or head_commit()
    lines = [line for line in done.stdout.splitlines() if line.strip()]
    if not lines or not lines[0].startswith("scope,method,"):
        raise SystemExit("recall wrote no table")
    lines = [line for line in lines if line.startswith("scope,") or line.startswith("all,") or line.startswith("interior,")]
    header = lines[0].split(",") + ["commit"]
    args.out.parent.mkdir(parents=True, exist_ok=True)
    with args.out.open("w", newline="") as handle:
        writer = csv.writer(handle, lineterminator="\n")
        writer.writerow(header)
        for line in lines[1:]:
            writer.writerow(line.split(",") + [commit])
    with args.sites.open() as handle:
        site_lines = handle.read().splitlines()
    if not site_lines or not site_lines[0].startswith("site,"):
        raise SystemExit("recall wrote no per-site table")
    stamped = [site_lines[0] + ",commit"]
    stamped.extend(line + "," + commit for line in site_lines[1:] if line.strip())
    args.sites.write_text("\n".join(stamped) + "\n")
    print(done.stdout, end="")
    print(f"wrote {args.out}")
    print(f"wrote {args.sites}")


if __name__ == "__main__":
    main()
