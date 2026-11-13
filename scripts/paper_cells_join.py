#!/usr/bin/env python3
"""Assemble one CSV per shared paper table from the producer fragments."""

from __future__ import annotations

import csv
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]


def load(path: Path) -> tuple[list[str], list[dict]]:
    with path.open(newline="") as handle:
        reader = csv.DictReader(handle)
        rows = list(reader)
        return list(reader.fieldnames or []), rows


def write_union(dest: Path, fragments: list[Path]) -> None:
    header: list[str] = []
    rows: list[dict] = []
    for path in fragments:
        fields, got = load(path)
        for name in fields:
            if name not in header:
                header.append(name)
        rows.extend(got)
    with dest.open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=header, extrasaction="raise")
        writer.writeheader()
        writer.writerows(rows)
    print(f"wrote {dest} ({len(rows)} rows)", flush=True)


def main() -> int:
    out = Path(sys.argv[1] if len(sys.argv) > 1 else ROOT / "benchmarks/paper-cells")
    write_union(out / "tab-cli.csv", [
        out / "frag-cli-ptm.csv",
        out / "frag-cli-q.csv",
        out / "frag-cli-rsm.csv",
    ])
    write_union(out / "tab-lit.csv", [
        out / "frag-lit-q6.csv",
        out / "frag-lit-cages.csv",
    ])
    return 0


if __name__ == "__main__":
    sys.exit(main())
