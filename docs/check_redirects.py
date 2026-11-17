#!/usr/bin/env python3
# Copyright (c) d-SEAMS developers
"""Every moved docname points at a page that exists, and the old path is gone.

docs/page_redirects.json is the map the HTML build publishes. A target
that is not an Org file, or an old path that still is, fails this check.
"""

import json
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
ORG = ROOT / "docs" / "orgmode"
MAP = ROOT / "docs" / "page_redirects.json"


def problems():
    if not MAP.is_file():
        return ["missing {}".format(MAP)]
    mapping = json.loads(MAP.read_text(encoding="utf-8"))
    errors = []
    if not isinstance(mapping, dict) or not mapping:
        return ["redirect map is empty"]
    for old, new in mapping.items():
        if not isinstance(old, str) or not isinstance(new, str) or not old or not new:
            errors.append("bad redirect entry {!r} -> {!r}".format(old, new))
            continue
        if old == new or old.startswith("/") or new.startswith("/"):
            errors.append("redirect {} does not move the page".format(old))
        if (ORG / (old + ".org")).is_file():
            errors.append("old page still exists: {}".format(old))
        if not (ORG / (new + ".org")).is_file():
            errors.append("redirect target is missing: {} -> {}".format(old, new))
    return errors


def main():
    errors = problems()
    if errors:
        print("\n".join(errors), file=sys.stderr)
        return 1
    print("redirects ok ({})".format(len(json.loads(MAP.read_text(encoding="utf-8")))))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
