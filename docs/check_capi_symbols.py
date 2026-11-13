#!/usr/bin/env python3
"""Fail if an exported C API symbol is missing from the reference.

``docs/capi/symbols.txt`` is the expansion of ``seams_*`` in
``src/libseams.map`` on capi/polish-2.14 (0127f41e) against
``include/seams.h``. When this tree grows an explicit map, those names
are checked too.
"""

import re
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
PAGES = ROOT / "docs" / "source" / "capi"
SYMBOLS = ROOT / "docs" / "capi" / "symbols.txt"
MAP = ROOT / "src" / "libseams.map"


def names_in(text):
    return set(re.findall(r"\bseams_[A-Za-z0-9_]+\b", text))


def main():
    wanted = [line.strip() for line in SYMBOLS.read_text().splitlines() if line.strip()]
    if len(wanted) != 77 or len(set(wanted)) != 77:
        raise SystemExit("symbols.txt must list 77 unique seams_* functions")
    book = "\n".join(path.read_text(encoding="utf-8") for path in PAGES.glob("*.rst"))
    present = names_in(book)
    missing = [name for name in wanted if name not in present]
    if missing:
        raise SystemExit("C API reference is missing:\n" + "\n".join(missing))
    if MAP.is_file():
        explicit = [
            name
            for name in names_in(MAP.read_text(encoding="utf-8"))
            if name != "seams_"
        ]
        # A wildcard map ("seams_*") yields the token seams_ only after the
        # star is stripped by the pattern above, so an explicit list is
        # anything that is a full function name.
        explicit = [name for name in explicit if not name.endswith("_")]
        absent = [name for name in explicit if name not in present]
        if absent:
            raise SystemExit("libseams.map symbols missing from the reference:\n" + "\n".join(absent))
    print("capi reference covers", len(wanted), "exports")


if __name__ == "__main__":
    sys.exit(main() or 0)
