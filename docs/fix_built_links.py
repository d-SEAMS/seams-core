#!/usr/bin/env python3
# Copyright (c) d-SEAMS developers
"""Repair links in the generated API HTML.

Doxygen names the tilde index ``index__7E`` and links it as ``#index_~``.
A letter index that links to ``#index_<letter>`` sometimes omits that id.
FORD links ``#moduleprocedure-<name>`` where the anchor id is ``proc-<name>``,
and it can emit the build machine's ``docs/ford`` path.
"""

import os
import re
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
BUILD = ROOT / "docs" / "build"
FORD_MARK = "/docs/ford/"
INDEX_HREF = re.compile(r'href="#(index_[^"]+)"')


def rewrite(tree, old, new):
    if not tree.is_dir():
        return 0
    changed = 0
    for path in tree.rglob("*.html"):
        text = path.read_text(encoding="utf-8", errors="replace")
        updated = text.replace(old, new)
        if updated != text:
            path.write_text(updated, encoding="utf-8")
            changed += 1
    return changed


def add_missing_index_ids(tree):
    if not tree.is_dir():
        return 0
    changed = 0
    for path in tree.rglob("*.html"):
        text = path.read_text(encoding="utf-8", errors="replace")
        ids = set(re.findall(r'\bid="([^"]+)"', text))

        def repl(match):
            frag = match.group(1)
            if frag in ids:
                return match.group(0)
            ids.add(frag)
            return 'id="{0}" href="#{0}"'.format(frag)

        updated = INDEX_HREF.sub(repl, text)
        if updated != text:
            path.write_text(updated, encoding="utf-8")
            changed += 1
    return changed


def relativize_ford(tree):
    """Point FORD's absolute docs/ford links at the published fortran tree."""
    if not tree.is_dir():
        return 0
    changed = 0
    for path in tree.rglob("*.html"):
        text = path.read_text(encoding="utf-8", errors="replace")
        if FORD_MARK not in text:
            continue

        def repl(match):
            quote, rest = match.group(1), match.group(2)
            target = (BUILD / "fortran" / rest).resolve()
            rel = os.path.relpath(target, path.parent.resolve())
            return "href={}{}{}".format(quote, rel, quote)

        updated = re.sub(
            r"href=(['\"])(?:file://)?[^'\"]*?/docs/ford/([^'\"]+)\1",
            repl,
            text,
        )
        if updated != text:
            path.write_text(updated, encoding="utf-8")
            changed += 1
    return changed


def main():
    api = rewrite(BUILD / "api-cpp", "#index_~", "#index__7E")
    anchors = add_missing_index_ids(BUILD / "api-cpp")
    fortran = rewrite(BUILD / "fortran", "#moduleprocedure-", "#proc-")
    paths = relativize_ford(BUILD / "fortran")
    print(
        "built links: api-cpp {} anchors {} fortran {} paths {}".format(
            api, anchors, fortran, paths
        )
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
