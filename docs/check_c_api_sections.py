#!/usr/bin/env python3
# Copyright (c) d-SEAMS developers
"""Fail if docs/api/c.md drops or gains a heading that tests/capi_sections.txt does not list.

docs/api/c.md is hand-written around the generated symbol reference, and a
generator that rewrites it can lose a hand-written section or a reference
entry without failing anything else. tests/capi_sections.txt lists every
Markdown heading of the page, one per line, outside fenced code. The
generators never write it: after an intended change to the page's headings,
run this script with --update and review the diff.
"""

import collections
import re
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
PAGE = ROOT / "docs" / "api" / "c.md"
MANIFEST = ROOT / "tests" / "capi_sections.txt"
HEADING = re.compile(r"#{1,6} \S.*")
FENCE = re.compile(r"(```|~~~)")


def headings(text):
    out = []
    fence = None
    for line in text.splitlines():
        stripped = line.strip()
        m = FENCE.match(stripped)
        if m:
            if fence is None:
                fence = m.group(1)
            elif stripped.startswith(fence):
                fence = None
            continue
        if fence is None and HEADING.fullmatch(line.rstrip()):
            out.append(line.rstrip())
    return out


def main(argv):
    found = headings(PAGE.read_text(encoding="utf-8"))
    if argv[1:] == ["--update"]:
        MANIFEST.write_text("\n".join(found) + "\n", encoding="utf-8")
        print(f"wrote {len(found)} headings to {MANIFEST.relative_to(ROOT)}")
        return 0
    if argv[1:]:
        raise SystemExit("usage: check_c_api_sections.py [--update]")
    listed = [line for line in MANIFEST.read_text(encoding="utf-8").splitlines() if line]
    have = collections.Counter(found)
    want = collections.Counter(listed)
    dropped = sorted((want - have).elements())
    unlisted = sorted((have - want).elements())
    if dropped or unlisted:
        lines = [f"{PAGE.relative_to(ROOT)} does not match {MANIFEST.relative_to(ROOT)}."]
        if dropped:
            lines.append("Headings the page no longer has:")
            lines += ["  " + h for h in dropped]
        if unlisted:
            lines.append("Headings the manifest does not list:")
            lines += ["  " + h for h in unlisted]
        lines.append("If the change is intended, run python3 docs/check_c_api_sections.py --update.")
        raise SystemExit("\n".join(lines))
    print(f"c.md keeps all {len(listed)} listed headings")
    return 0


if __name__ == "__main__":
    sys.exit(main(sys.argv))
