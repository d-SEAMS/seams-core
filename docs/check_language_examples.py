#!/usr/bin/env python3
# Copyright (c) d-SEAMS developers
"""Count the code examples on the eight language pages and require each to be tested.

A literalinclude is compiled by docs/run_snippets.py. A language fence must
assert its result. A ``capi`` line on the C page is an expected-output
example checked by docs/check_c_api_tour.py.
"""

import re
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "docs"))
import run_snippets  # noqa: E402

PAGES = ("c", "cpp", "python", "fortran", "julia", "r", "rust", "lua")
FENCE = re.compile(r"```([^\n]*)\n(.*?)```", re.S)


def fail(message):
    print(message, file=sys.stderr)
    raise SystemExit(1)


def main():
    programs = 0
    asserted = 0
    documented = 0
    for name in PAGES:
        path = ROOT / "docs" / "api" / (name + ".md")
        text = path.read_text(encoding="utf-8")
        picked = run_snippets.snippets_of(path)
        if not picked:
            fail(name + " has no executed example")
        for kind, lang, body in picked:
            if kind == "include":
                target = run_snippets.resolve_include(path, lang)
                if not target.is_file():
                    fail("missing " + str(target))
                programs += 1
                continue
            if lang in {"bash", "sh"}:
                continue
            if lang == "python" and "assert " not in body and not body.lstrip().startswith(">>>"):
                fail("Python example on " + name + " does not assert")
            if lang == "rust" and ("fn main" not in body or "assert" not in body):
                fail("Rust example is not a complete asserted program")
            asserted += 1
        for match in FENCE.finditer(text):
            lang = match.group(1).strip()
            body = match.group(2).strip()
            if lang or not body.startswith("capi "):
                continue
            documented += 1
    total = programs + asserted + documented
    print(
        "language examples {} ({} compiled sources, {} asserted snippets, {} C output lines)".format(
            total, programs, asserted, documented
        )
    )
    if total < 1:
        fail("no language examples")
    return 0


if __name__ == "__main__":
    sys.exit(main())
