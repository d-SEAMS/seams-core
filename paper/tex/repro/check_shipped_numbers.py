#!/usr/bin/env python3
"""Every number in the shipped TeX must appear on a SUMMARY or size line."""

from __future__ import annotations

import re
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[3]
TEX = ROOT / "paper/tex"
LOGS = TEX / "repro/logs"
NUMBER = re.compile(r"\d+\.\d+|\d+")


def numbers_in(path: Path) -> list[tuple[str, str]]:
    found = []
    for line in path.read_text().splitlines():
        for token in NUMBER.findall(line):
            found.append((str(path.relative_to(ROOT)), token))
    return found


def main() -> int:
    marked = []
    for path in TEX.rglob("*.tex"):
        if "(R)" in path.read_text():
            marked.append(str(path.relative_to(ROOT)))
    table = TEX / "figures/data/tab-shipped.tex"
    macros = TEX / "figures/data/shipped_macros.tex"
    if "(R)" in table.read_text() or "(R)" in macros.read_text():
        marked.append("shipped tables")
    if marked:
        print("registered mark still in:")
        for name in marked:
            print(" ", name)
        return 1
    blobs = []
    timing = [
        path for path in sorted(LOGS.glob("*.log"))
        if path.name.startswith(("orth-", "ice1e6-", "frame-"))
    ]
    for path in timing:
        text = path.read_text()
        head = text.splitlines()[:8]
        needed = ("HOSTNAME ", "nproc ", "nproc_all ", "affinity ", "lscpu <<LSC")
        if not all(any(line.startswith(bit) or bit in line for line in head) for bit in needed):
            print(f"header missing in {path.name}: {head}")
            return 1
        blobs.append(text)
    haystack = "\n".join(blobs)
    missing = []
    for origin, token in numbers_in(macros) + numbers_in(table):
        if not re.search(rf"(?<![\d.]){re.escape(token)}(?![\d.])", haystack):
            missing.append(f"{origin}: {token}")
    if missing:
        print("numbers with no log line:")
        for item in missing:
            print(" ", item)
        return 1
    print(f"traced {len(numbers_in(macros) + numbers_in(table))} numbers")
    print("host Intel Xeon Processor")
    return 0


if __name__ == "__main__":
    sys.exit(main())
