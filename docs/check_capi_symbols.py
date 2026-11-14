#!/usr/bin/env python3
"""Fail if a SEAMS_API declaration is missing from the generated reference."""

import re
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "docs"))

from render_capi import HEADER, OUT, SYMBOLS, parse_header  # noqa: E402


def names_in(text):
    return set(re.findall(r"\bseams_[A-Za-z0-9_]+\b", text))


def main():
    declared = [record["name"] for record in parse_header(HEADER.read_text(encoding="utf-8"))]
    if len(declared) != len(set(declared)):
        raise SystemExit("include/seams.h declares a SEAMS_API name twice")
    listed = [line.strip() for line in SYMBOLS.read_text(encoding="utf-8").splitlines() if line.strip()]
    if sorted(listed) != sorted(declared) or len(listed) != len(set(listed)):
        raise SystemExit(
            "docs/capi/symbols.txt is not the SEAMS_API list in include/seams.h. "
            "Run python3 docs/render_capi.py."
        )
    book = "\n".join(path.read_text(encoding="utf-8") for path in OUT.glob("*.rst"))
    present = names_in(book)
    missing = [name for name in declared if name not in present]
    if missing:
        raise SystemExit("C API reference is missing:\n" + "\n".join(missing))
    stale = sorted(
        name
        for name in (
            "seams_chill_plus_n",
            "seams_last_error_message",
            "seams_neighbours_nrows",
            "seams_q4q6_n",
            "seams_rings_histogram_n",
            "seams_rings_labels_n",
            "seams_rings_nrings",
            "seams_system_natoms",
        )
        if name in present
    )
    if stale:
        raise SystemExit("C API reference still names functions this header does not declare:\n" + "\n".join(stale))
    print("capi reference covers", len(declared), "exports from include/seams.h")
    return 0


if __name__ == "__main__":
    sys.exit(main() or 0)
