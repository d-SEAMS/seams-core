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
# LDoc writes `<a name = "symbol">`. The link checker matches fragments to id.
NAME_ANCHOR = re.compile(
    r"<a\b([^>]*?)\bname\s*=\s*([\"'])([^\"']+)\2([^>]*)>",
    re.IGNORECASE,
)
# Documenter appends #Lstart-Lend. GitHub's static HTML has no such id, and a
# commit that is not on the remote 404s. The file on main is the published source.
# Documenter records the clone it was built from and a #L range. GitHub's
# file page has no such id, and a commit that was never pushed 404s.
# A citation of a historical commit has no #L range and is left alone.
SOURCE_LINE = re.compile(
    r"https://github\.com/(?:HaoZeke|d-SEAMS)/seams-core/blob/"
    r"(?:[0-9a-f]{7,40}|master|HEAD)"
    r"(/[^\"'#\s<>]+)#L\d+(?:-L\d+)?"
)
OBFUSCATED_MAILTO = re.compile(r'href="(mai(?:[^"]*?))"')


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


def add_name_ids(tree):
    """Give each named anchor an id so a fragment link resolves."""
    if not tree.is_dir():
        return 0
    changed = 0
    for path in tree.rglob("*.html"):
        text = path.read_text(encoding="utf-8", errors="replace")

        def repl(match):
            pre, quote, name, post = match.groups()
            if re.search(r"\bid\s*=", pre + post, re.IGNORECASE):
                return match.group(0)
            return "<a{0} name={1}{2}{1} id={1}{2}{1}{3}>".format(
                pre, quote, name, post
            )

        updated = NAME_ANCHOR.sub(repl, text)
        if updated != text:
            path.write_text(updated, encoding="utf-8")
            changed += 1
    return changed


def publish_source_lines(tree):
    """Point Documenter line links at the file on main, without a fragment."""
    if not tree.is_dir():
        return 0
    changed = 0
    for path in tree.rglob("*.html"):
        text = path.read_text(encoding="utf-8", errors="replace")
        updated = SOURCE_LINE.sub(
            r"https://github.com/d-SEAMS/seams-core/blob/main\1", text
        )
        if updated != text:
            path.write_text(updated, encoding="utf-8")
            changed += 1
    return changed


def decode_doxygen_mail(tree):
    """Turn Doxygen's concatenated mail attributes into mailto links."""
    if not tree.is_dir():
        return 0
    changed = 0
    for path in tree.rglob("*.html"):
        text = path.read_text(encoding="utf-8", errors="replace")

        def repl(match):
            address = match.group(1).replace("'+'", "").replace("+'", "")
            address = address.replace("'", "")
            return 'href="{}"'.format(address)

        updated = OBFUSCATED_MAILTO.sub(repl, text)
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
    lua = add_name_ids(BUILD / "lua")
    sources = publish_source_lines(BUILD)
    mail = decode_doxygen_mail(BUILD / "api-cpp")
    print(
        "built links: api-cpp {} anchors {} fortran {} paths {} lua {} sources {} mail {}".format(
            api, anchors, fortran, paths, lua, sources, mail
        )
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
