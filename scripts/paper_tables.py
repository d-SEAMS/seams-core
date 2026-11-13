#!/usr/bin/env python3
"""Write one CSV per paper table from benchmarks/paper-cells/sitting.csv.

The sitting is the only timing source. These files are selections of
its rows. A table that is not a 110592-molecule or hydrate-cage timing
is left to its own rule.
"""

from __future__ import annotations

import csv
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]


def load(path: Path) -> list[dict]:
    with path.open(newline="") as handle:
        return list(csv.DictReader(handle))


def write(path: Path, rows: list[dict]) -> None:
    if not rows:
        raise SystemExit(f"no rows for {path.name}")
    with path.open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(rows[0]))
        writer.writeheader()
        writer.writerows(rows)
    print(f"wrote {path} ({len(rows)} rows)", flush=True)


def select(rows: list[dict], **want) -> list[dict]:
    found = []
    for row in rows:
        if all(row.get(key) == value for key, value in want.items()):
            found.append(row)
    return found


def main() -> int:
    out = Path(sys.argv[1] if len(sys.argv) > 1 else ROOT / "benchmarks/paper-cells")
    rows = load(out / "sitting.csv")
    ice = [row for row in rows if row["phase"] == "ice"]
    write(out / "tab-rings-214.csv", [
        row for row in ice
        if row["operation"] == "rings"
        and row["code"] in ("d-SEAMS", "rings-statistics-matter")
    ])
    write(out / "tab-ql-eq.csv", [
        row for row in ice
        if row["operation"] in ("q4", "q6", "q4-bar", "q6-bar", "q4+q6", "q4-core", "q6-core")
        and row["code"] == "d-SEAMS"
    ])
    write(out / "tab-neighbours.csv", [
        row for row in ice if row["operation"] == "neighbours"
    ])
    write(out / "tab-chill.csv", [
        row for row in ice if row["operation"] == "chill+"
    ])
    write(out / "frag-cli-q.csv", [
        row for row in ice if row["operation"] == "q4+q6" and row["code"] == "d-SEAMS"
    ])
    write(out / "frag-cli-rsm.csv", [
        row for row in ice
        if row["operation"] == "rings" and row["code"] == "rings-statistics-matter"
    ])
    write(out / "frag-lit-q6.csv", [
        row for row in ice if row["operation"] in ("q4", "q6", "q4-core", "q6-core")
    ])
    cages = [
        row for row in rows
        if row["operation"] == "cages" and row["phase"] != "ice"
    ]
    write(out / "frag-lit-cages.csv", cages)
    write(out / "tab-cages.csv", cages)
    return 0


if __name__ == "__main__":
    sys.exit(main())
