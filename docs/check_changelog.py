#!/usr/bin/env python3
# Copyright (c) d-SEAMS developers
"""Render the 2.14 release notes from towncrier fragments and validate them.

changelog.d is the only source. Top-level fragments are the 2.14.0 notes,
and changelog.d/<X.Y.Z> holds each 2.14 patch release. The rendering has
one section per release, newest first, and the newest must be the project
version in meson.build. The docs build writes the rendering to
docs/source/reference/release-2.14-notes.rst. ``--check`` writes it to a
temporary file instead and validates that file.
"""

from __future__ import annotations

import importlib.util
import pathlib
import re
import sys
import tempfile
import tomllib

ROOT = pathlib.Path(__file__).resolve().parents[1]
CONFIG = ROOT / "towncrier.toml"
NOTES = ROOT / "docs" / "source" / "reference" / "release-2.14-notes.rst"
PAGE = ROOT / "docs" / "orgmode" / "reference" / "release-2.14.org"
SERIES = "2.14"
SERIES_BASE = SERIES + ".0"
MINIMUM_FRAGMENTS = 100
IGNORED = {"README.md", ".gitkeep"}
CODE = re.compile(r"`([^`]+)`")
VERSION_DIR = re.compile(r"^\d+\.\d+\.\d+$")
VERSION_RULE = "-"
TYPE_RULE = "~"


def read_config():
    data = tomllib.loads(CONFIG.read_text())
    towncrier = data.get("tool", {}).get("towncrier", {})
    directory = towncrier.get("directory")
    if not isinstance(directory, str) or not directory:
        raise ValueError("tool.towncrier.directory is missing")
    types = tuple(
        item.get("directory")
        for item in data.get("tool", {}).get("towncrier", {}).get("type", [])
        if isinstance(item.get("directory"), str)
    )
    if not types:
        raise ValueError("tool.towncrier.type is missing")
    prefix = towncrier.get("orphan_prefix")
    if not isinstance(prefix, str) or not prefix:
        raise ValueError("tool.towncrier.orphan_prefix is missing")
    headings = {
        item["directory"]: item.get("name", item["directory"].title())
        for item in data.get("tool", {}).get("towncrier", {}).get("type", [])
        if isinstance(item.get("directory"), str)
    }
    return ROOT / directory, types, prefix, headings


def project_version():
    """The release version, read the way docs/source/conf.py reads it."""
    script = ROOT / "scripts" / "project_version.py"
    spec = importlib.util.spec_from_file_location("project_version", script)
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module.project_version()


def version_key(version):
    return tuple(int(part) for part in version.split("."))


def load_directory(directory, version, errors, found):
    for path in sorted(directory.iterdir()):
        if path.is_dir() and directory == FRAGMENTS:
            continue
        if not path.is_file():
            errors.append("{} is not a news fragment".format(path.relative_to(FRAGMENTS)))
            continue
        if path.name in IGNORED:
            continue
        name = str(path.relative_to(FRAGMENTS))
        if path.suffix != ".md":
            errors.append("{} is not a news fragment".format(name))
            continue
        match = NAME.fullmatch(path.name)
        if match is None:
            errors.append("{} does not match fragment naming rules".format(name))
            continue
        body = path.read_text().strip()
        if not body:
            errors.append("{} is empty".format(name))
            continue
        found.append((version, match.group(2), name, body))


def load():
    """Return (errors, fragments) with fragments as (version, type, name, body)."""
    errors = []
    found = []
    if not FRAGMENTS.is_dir():
        return ["missing {}".format(FRAGMENTS.relative_to(ROOT))], []
    load_directory(FRAGMENTS, SERIES_BASE, errors, found)
    for path in sorted(p for p in FRAGMENTS.iterdir() if p.is_dir()):
        version = path.name
        if not VERSION_DIR.fullmatch(version):
            errors.append("{} is not a release directory".format(path.relative_to(ROOT)))
            continue
        if not version.startswith(SERIES + ".") or version_key(version) <= version_key(
            SERIES_BASE
        ):
            errors.append(
                "{} is not a {} patch release".format(path.relative_to(ROOT), SERIES)
            )
            continue
        before = len(found)
        load_directory(path, version, errors, found)
        if len(found) == before:
            errors.append("{} has no fragments".format(path.relative_to(ROOT)))
    seen = {}
    for _version, _typ, name, body in found:
        key = " ".join(body.split())
        if key in seen:
            errors.append("{} repeats {}".format(name, seen[key]))
        else:
            seen[key] = name
    return errors, found


def ordered(fragments):
    """Fragments in rendering order: newest release, then towncrier type."""
    releases = sorted({item[0] for item in fragments}, key=version_key, reverse=True)
    result = []
    for version in releases:
        for typ in TYPES:
            result.extend(
                item for item in fragments if item[0] == version and item[1] == typ
            )
    return releases, result


def rst_body(body):
    """Markdown code spans become reStructuredText literals."""
    return CODE.sub(lambda match: "``{}``".format(match.group(1)), body)


def render(fragments):
    lines = [
        ".. Generated by docs/check_changelog.py from changelog.d.",
        ".. Edit the fragment, not this file.",
        "",
    ]
    releases, items = ordered(fragments)
    for version in releases:
        lines.append(version)
        lines.append(VERSION_RULE * len(version))
        lines.append("")
        for typ in TYPES:
            bodies = [body for ver, kind, _name, body in items if ver == version and kind == typ]
            if not bodies:
                continue
            title = HEADINGS[typ]
            lines.append(title)
            lines.append(TYPE_RULE * len(title))
            lines.append("")
            for body in bodies:
                text = rst_body(body).splitlines()
                lines.append("- {}".format(text[0]))
                for extra in text[1:]:
                    lines.append("  {}".format(extra) if extra.strip() else "")
                lines.append("")
    return "\n".join(lines).rstrip() + "\n"


def parse(text):
    """Read a rendering back as (errors, releases, entries).

    entries are (version, type, body) with literals restored to Markdown.
    """
    errors = []
    releases = []
    entries = []
    kinds = {title: typ for typ, title in HEADINGS.items()}
    version = None
    kind = None
    current = None

    def close():
        nonlocal current
        if current is None:
            return
        body = "\n".join(current).strip()
        entries.append((version, kind, re.sub(r"``([^`]+)``", r"`\1`", body)))
        current = None

    lines = text.splitlines()
    index = 0
    while index < len(lines):
        line = lines[index]
        rule = lines[index + 1] if index + 1 < len(lines) else ""
        number = index + 1
        index += 1
        if line.startswith(".. "):
            if version is not None:
                errors.append("line {}: comment inside the notes".format(number))
            continue
        if line.strip() == "":
            if current is not None and rule.startswith("  ") and rule.strip():
                current.append("")
            else:
                close()
            continue
        if line.startswith("- "):
            close()
            if kind is None:
                errors.append("line {}: note outside a type section".format(number))
            current = [line[2:]]
            continue
        if line.startswith("  ") and current is not None:
            current.append(line[2:])
            continue
        if rule and len(rule) == len(line) and set(rule) == {VERSION_RULE}:
            close()
            index += 1
            if not VERSION_DIR.fullmatch(line):
                errors.append("line {}: release heading {!r} is not X.Y.Z".format(number, line))
            releases.append(line)
            version = line
            kind = None
            continue
        if rule and len(rule) == len(line) and set(rule) == {TYPE_RULE}:
            close()
            index += 1
            if line not in kinds:
                errors.append("line {}: {!r} is not a towncrier type".format(number, line))
            if version is None:
                errors.append("line {}: type heading before a release heading".format(number))
            kind = kinds.get(line)
            continue
        close()
        errors.append("line {}: unexpected text {!r}".format(number, line[:60]))
    close()
    return errors, releases, entries


def check(fragments, text, version):
    errors = []
    base = [item for item in fragments if item[0] == SERIES_BASE]
    if len(base) < MINIMUM_FRAGMENTS:
        errors.append(
            "expected at least {} towncrier fragments for {}, found {}".format(
                MINIMUM_FRAGMENTS, SERIES_BASE, len(base)
            )
        )
    syntax, releases, parsed = parse(text)
    errors.extend(syntax)
    wanted_releases, items = ordered(fragments)
    if releases != wanted_releases:
        errors.append(
            "release headings are {}, fragments give {}".format(releases, wanted_releases)
        )
    if not releases or releases[0] != version:
        errors.append(
            "newest release heading is {}, meson.build version is {}".format(
                releases[0] if releases else None, version
            )
        )
    if len(set(releases)) != len(releases):
        errors.append("a release heading repeats: {}".format(releases))
    wanted = [(ver, typ, body) for ver, typ, _name, body in items]
    names = {(ver, typ, body): name for ver, typ, name, body in items}
    if parsed != wanted:
        missing = [entry for entry in wanted if entry not in parsed]
        extra = [entry for entry in parsed if entry not in wanted]
        repeated = [entry for entry in wanted if parsed.count(entry) > 1]
        for entry in missing:
            errors.append("orphan fragment: {}".format(names[entry]))
        for entry in extra:
            errors.append(
                "note is not a fragment of {} {}: {}".format(
                    entry[0], entry[1], entry[2].splitlines()[0][:80]
                )
            )
        for entry in repeated:
            errors.append("duplicate note: {}".format(names[entry]))
        if not missing and not extra and not repeated:
            errors.append("rendered notes are not in release and type order")
    page = PAGE.read_text()
    if ".. include:: {}".format(NOTES.name) not in page:
        errors.append("{} does not include {}".format(PAGE.name, NOTES.name))
    title = re.search(r"^#\+TITLE:\s*(.+?)\s*$", page, re.M)
    expected = "{} release notes".format(version)
    if title is None or title.group(1) != expected:
        errors.append(
            "{} title is {!r}, expected {!r}".format(
                PAGE.name, title.group(1) if title else None, expected
            )
        )
    return errors


def main(argv):
    global FRAGMENTS, TYPES, HEADINGS, NAME
    if argv not in ([], ["--check"]):
        print("usage: docs/check_changelog.py [--check]", file=sys.stderr)
        return 2
    try:
        FRAGMENTS, TYPES, orphan_prefix, HEADINGS = read_config()
    except ValueError as exc:
        print(str(exc), file=sys.stderr)
        return 1
    NAME = re.compile(
        r"^(" + re.escape(orphan_prefix) + r"[A-Za-z0-9][A-Za-z0-9.-]*|[0-9]+)\.("
        + "|".join(TYPES)
        + r")\.md$"
    )
    version = project_version()
    errors, fragments = load()
    if errors:
        print("\n".join(errors), file=sys.stderr)
        return 1
    rendered = render(fragments)
    with tempfile.TemporaryDirectory() as scratch:
        target = pathlib.Path(scratch) / NOTES.name if argv else NOTES
        target.parent.mkdir(parents=True, exist_ok=True)
        target.write_text(rendered)
        errors = check(fragments, target.read_text(), version)
    if errors:
        print("\n".join(errors), file=sys.stderr)
        return 1
    releases, _items = ordered(fragments)
    counts = ", ".join(
        "{} {}".format(ver, sum(1 for item in fragments if item[0] == ver)) for ver in releases
    )
    print(
        "{} release notes render {} fragments from {} ({})".format(
            version, len(fragments), FRAGMENTS.relative_to(ROOT), counts
        )
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main(sys.argv[1:]))
