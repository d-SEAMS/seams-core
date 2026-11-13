#!/usr/bin/env python3
"""Rerun the engine counts behind the main-text cage and label numbers.

Writes one report that names the archive it hashed and the seams binary it
ran. A row is ``match`` when the fresh count equals the deposited record,
``mismatch`` when it does not, and ``not-rerun`` when the input is not in
this checkout (the full nucleation trajectory, the comparator installs, the
GPU node). The process exits 1 when any rerun mismatches.

Usage:
  python3 scripts/rerun_main_text.py \\
      --seams build/seams --input input/traj \\
      --archive program-files/seams-core-v2.12.0.tar.gz \\
      --out rerun-main
"""

from __future__ import annotations

import argparse
import hashlib
import os
import re
import subprocess
import sys
from pathlib import Path


def affinity_label() -> str:
    try:
        allowed = sorted(os.sched_getaffinity(0))
    except (AttributeError, OSError):
        return "not recorded"
    if not allowed:
        return "not recorded"
    if allowed == list(range(allowed[0], allowed[-1] + 1)):
        return f"{allowed[0]}-{allowed[-1]}" if len(allowed) > 1 else str(allowed[0])
    return ",".join(str(cpu) for cpu in allowed)


def lscpu_field(label: str) -> str:
    proc = subprocess.run(["lscpu"], check=True, capture_output=True, text=True)
    prefix = label + ":"
    for line in proc.stdout.splitlines():
        if line.startswith(prefix):
            return line.split(":", 1)[1].strip()
    raise RuntimeError(f"lscpu has no {label}")


def sha256(path: Path) -> str:
    h = hashlib.sha256()
    with path.open("rb") as fh:
        for chunk in iter(lambda: fh.read(1 << 20), b""):
            h.update(chunk)
    return h.hexdigest()


def run(cmd: list[str]) -> str:
    proc = subprocess.run(cmd, check=False, capture_output=True, text=True)
    if proc.returncode != 0:
        raise RuntimeError(
            f"{' '.join(cmd)} exited {proc.returncode}\n{proc.stderr}"
        )
    return proc.stdout


def cage_line(seams: Path, traj: Path, *extra: str) -> str:
    return run([str(seams), "cages", str(traj), "--cutoff", "3.5", *extra]).strip()


def field(line: str, name: str) -> int:
    match = re.search(rf"(?:^|\s){name}\s+(-?\d+)", line)
    if match is None:
        raise RuntimeError(f"no {name} in {line!r}")
    return int(match.group(1))


def signature_cages(line: str) -> int:
    return field(line, "cages")


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--seams", type=Path, required=True)
    parser.add_argument("--input", type=Path, required=True)
    parser.add_argument("--archive", type=Path, required=True)
    parser.add_argument(
        "--source",
        type=Path,
        default=None,
        help="git checkout the seams binary was built from; its tree hash is recorded",
    )
    parser.add_argument("--out", type=Path, required=True)
    parser.add_argument(
        "--nucleation",
        type=Path,
        default=None,
        help="full nucleation.lammpstrj if it is available; frame 148 is checked",
    )
    args = parser.parse_args()
    if not args.seams.is_file():
        sys.exit(f"no seams binary at {args.seams}")
    if not args.archive.is_file():
        sys.exit(f"no archive at {args.archive}")
    args.out.mkdir(parents=True, exist_ok=True)

    archive_hash = sha256(args.archive)
    version = run([str(args.seams), "--version"]).strip()
    source_tree = ""
    source_commit = ""
    if args.source is not None:
        source_commit = run(
            ["git", "-C", str(args.source), "rev-parse", "HEAD"]
        ).strip()
        source_tree = run(
            ["git", "-C", str(args.source), "write-tree"]
        ).strip()
    rows: list[str] = []
    mismatches: list[str] = []

    def check(name: str, got: int, expected: int, source: str) -> None:
        status = "match" if got == expected else "mismatch"
        rows.append(f"{status}\t{name}\texpected {expected}\tgot {got}\t{source}")
        if status == "mismatch":
            mismatches.append(name)

    def absent(name: str, reason: str) -> None:
        rows.append(f"not-rerun\t{name}\t{reason}")

    traj = args.input
    seams = args.seams

    # Cubic mW, frame 1, 4096 molecules. One double-diamond cage for every
    # two molecules: the conflict graph is the line graph of a 4-regular
    # bipartite graph, so the maximum packing is n/2. bench cages_v2
    # deposited 2040, the count from the unpatched library.
    mw = traj / "mW_cubic.lammpstrj"
    if mw.is_file():
        ddc = cage_line(seams, mw, "--type", "1", "--graph", "cutoff", "--signature", "ddc")
        check("mW DDC cages", signature_cages(ddc), 2048, "maximum packing n/2")
        seeded = cage_line(seams, mw, "--type", "1", "--graph", "seeded")
        check("mW seeded hexagonal", field(seeded, "hexagonal"), 0, "mW cubic")
        check("mW seeded cubic", field(seeded, "cubic"), 4096, "mW cubic")
        check("mW seeded both", field(seeded, "both"), 0, "mW cubic")
        check("mW seeded water", field(seeded, "water"), 0, "mW cubic")
    else:
        absent("mW", f"missing {mw}")

    # GenIce unit cells that ship with the engine. sI is the 12 A cell.
    crystals = {
        "sI": ("genice_sI.lammpstrj", {"512": 2, "51262": 6, "51264": 0, "51268": 0}),
        "sII": ("genice_sII.lammpstrj", {"512": 16, "51264": 8, "51262": 0, "51268": 0}),
        "sH": ("genice_sH.lammpstrj", {"512": 12, "51262": 0, "51264": 0, "51268": 4, "sh": 8}),
    }
    for name, (filename, expect) in crystals.items():
        path = traj / filename
        if not path.is_file():
            absent(name, f"missing {path}")
            continue
        for spec, n in expect.items():
            line = cage_line(seams, path, "--type", "1", "--graph", "cutoff", "--signature", spec)
            check(f"{name} {spec}", signature_cages(line), n, filename)

    # Nucleation frame 1 is in the engine tree. Frame 148 needs the full dump.
    # Deposited disjoint seeded counts (walks/figshare-nucleation.txt):
    # frame 1 is 141 / 98 / 41, frame 148 is 255 / 363 / 377.
    frame1 = traj / "nucleation_f1.lammpstrj"
    if frame1.is_file():
        line = cage_line(seams, frame1, "--type", "1", "--graph", "seeded")
        check("nucleation frame 1 ih", field(line, "hexagonal"), 141, "walks frame 1")
        check("nucleation frame 1 ic", field(line, "cubic"), 98, "walks frame 1")
        check("nucleation frame 1 both", field(line, "both"), 41, "walks frame 1")
    else:
        absent("nucleation frame 1", f"missing {frame1}")

    full = args.nucleation
    if full is None:
        candidate = traj / "nucleation.lammpstrj"
        full = candidate if candidate.is_file() else None
    if full is not None and full.is_file():
        line = cage_line(
            seams, full, "--type", "1", "--graph", "seeded", "--frame", "148"
        )
        check("nucleation frame 148 ih", field(line, "hexagonal"), 255, "walks frame 148")
        check("nucleation frame 148 ic", field(line, "cubic"), 363, "walks frame 148")
        check("nucleation frame 148 both", field(line, "both"), 377, "walks frame 148")
    else:
        absent(
            "nucleation frame 148",
            "full nucleation.lammpstrj is not in this checkout "
            "(figshare 10.6084/m9.figshare.11448702); frame 1 was rerun",
        )

    absent(
        "Table 2 head-to-head",
        "needs OVITO, freud, GenIce2, cycless and R.I.N.G.S. "
        "from prepare_h2h.sh; not installed here",
    )
    absent(
        "GPU batch",
        "no NVIDIA GPU on this machine; the deposited record names no engine commit",
    )
    absent(
        "legacy seeding and brine campaign",
        "starting liquids were not archived; the campaign cannot be regenerated",
    )

    dirty = ""
    if args.source is not None:
        dirty = run(
            ["git", "-C", str(args.source), "status", "--porcelain", "--untracked-files=no"]
        ).strip()
    report = args.out / "report.txt"
    body = [
        f"archive_sha256 {archive_hash}",
        f"archive {args.archive}",
        f"seams {args.seams}",
        f"seams_version {version}",
        f"source_commit {source_commit or 'not recorded'}",
        f"cpu_model {lscpu_field('Model name')}",
        f"cpu_cores {lscpu_field('CPU(s)')}",
        f"affinity {affinity_label()}",
        f"source_tree {source_tree or 'not recorded'}",
        "source_dirty "
        + ("yes" if dirty else "no" if args.source is not None else "not recorded"),
        "note the archive hash is the file passed to --archive; "
        "source_commit is the checkout the binary was built from",
        "",
        *rows,
        "",
        f"mismatches {len(mismatches)}",
    ]
    report.write_text("\n".join(body) + "\n", encoding="utf-8")
    print(report.read_text(), end="")
    return 1 if mismatches else 0


if __name__ == "__main__":
    sys.exit(main())
