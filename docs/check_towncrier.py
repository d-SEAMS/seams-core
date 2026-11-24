#!/usr/bin/env python3
# Copyright (c) d-SEAMS developers
"""Require every sourced CHANGELOG.md release section to be its towncrier draft.

Top-level changelog.d fragments are the 2.14.0 release, and changelog.d/<X.Y.Z>
holds each 2.14 patch release, as docs/check_changelog.py reads them. Each of
those releases is rendered by towncrier from its own fragments, and its
CHANGELOG.md section has to equal that draft: the same entries, under the same
type headings, in towncrier order. A 2.14 patch-release section without a
fragment directory fails, and so does a fragment directory without a section.
"""

from __future__ import annotations

import difflib
import re
import shutil
import subprocess
import sys
import tempfile
import tomllib
from collections import Counter
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
CHANGELOG = ROOT / "CHANGELOG.md"
CONFIG = ROOT / "towncrier.toml"
RELEASE = re.compile(r"^## \[(\d+\.\d+\.\d+)\] - (\d{4}-\d{2}-\d{2})$", re.M)

sys.path.insert(0, str(ROOT / "docs"))
from check_changelog import (  # noqa: E402
    IGNORED,
    MINIMUM_FRAGMENTS,
    SERIES,
    SERIES_BASE,
    VERSION_DIR,
    version_key,
)


def read_config():
    data = tomllib.loads(CONFIG.read_text())
    tool = data.get("tool", {}).get("towncrier", {})
    directory = tool.get("directory")
    if not isinstance(directory, str) or not directory:
        raise ValueError("tool.towncrier.directory is missing")
    template = tool.get("template")
    if not isinstance(template, str) or not template:
        raise ValueError("tool.towncrier.template is missing")
    return directory, template


def sources(fragments):
    """Return (errors, {version: [fragment paths]}) for each sourced release."""
    errors = []
    found = {
        SERIES_BASE: sorted(
            path for path in fragments.iterdir() if path.is_file() and path.name not in IGNORED
        )
    }
    for path in sorted(p for p in fragments.iterdir() if p.is_dir()):
        if not VERSION_DIR.fullmatch(path.name):
            errors.append("{} is not a release directory".format(path.relative_to(ROOT)))
            continue
        found[path.name] = sorted(
            item for item in path.iterdir() if item.is_file() and item.name not in IGNORED
        )
    return errors, found


def draft(version, date, paths, directory, template):
    """Render one release from its own fragments, or return (None, error)."""
    with tempfile.TemporaryDirectory(prefix="seams-towncrier-") as scratch:
        tree = Path(scratch)
        shutil.copy(CONFIG, tree / CONFIG.name)
        (tree / template).parent.mkdir(parents=True, exist_ok=True)
        shutil.copy(ROOT / template, tree / template)
        (tree / directory).mkdir(parents=True)
        for path in paths:
            shutil.copy(path, tree / directory / path.name)
        done = subprocess.run(
            [
                sys.executable,
                "-m",
                "towncrier",
                "build",
                "--draft",
                "--version",
                version,
                "--date",
                date,
            ],
            cwd=str(tree),
            check=False,
            capture_output=True,
            text=True,
        )
    if done.returncode != 0:
        return None, (done.stdout + done.stderr).strip() or "towncrier draft failed"
    marker = done.stdout.find("## [{}]".format(version))
    if marker < 0:
        return None, "towncrier draft has no {} heading".format(version)
    return done.stdout[marker:].strip() + "\n", None


def section_spans(text):
    """(version, start, end) for each release section, end exclusive."""
    matches = list(RELEASE.finditer(text))
    return [
        (
            match.group(1),
            match.start(),
            matches[index + 1].start() if index + 1 < len(matches) else len(text),
        )
        for index, match in enumerate(matches)
    ]


def section_text(text, version):
    for name, start, end in section_spans(text):
        if name == version:
            return text[start:end]
    return None


def section_entries(section):
    """Ordered (type heading, body) entries under one release heading."""
    found = []
    category = None
    current = None

    def close():
        nonlocal current
        if current is not None and category is not None:
            body = "\n".join(current).strip()
            if body:
                found.append((category, body))
        current = None

    for line in section.splitlines():
        if line.startswith("## "):
            close()
            continue
        if line.startswith("### "):
            close()
            category = line[4:].strip()
            continue
        if line.startswith("- "):
            close()
            current = [line[2:]]
            continue
        if current is not None and (line.startswith(" ") or line.strip() == ""):
            current.append(line[2:] if line.startswith("  ") else line)
            continue
        close()
    close()
    return found


def preview(entry):
    return entry.splitlines()[0].strip()[:80]


def edited(old, new):
    if old == new:
        return False
    shorter = min(len(old), len(new))
    if shorter >= 12 and (new.startswith(old) or old.startswith(new)):
        return True
    if shorter < 24:
        return False
    return difflib.SequenceMatcher(a=old, b=new, autojunk=False).ratio() >= 0.8


def grouped(items):
    order = []
    groups = {}
    for category, entry in items:
        if category not in groups:
            order.append(category)
            groups[category] = []
        groups[category].append(entry)
    return order, groups


def wrong_section(found, expected, entry):
    return "entry is in the wrong section (found in {} {}, expected in {} {}): {}".format(
        found[0], found[1], expected[0], expected[1], preview(entry)
    )


def entry_errors(version, expected, actual, expected_at, found_at):
    """Name what differs between one rendered draft and its CHANGELOG.md section."""
    if Counter(expected) == Counter(actual):
        if expected != actual:
            return ["{} entries are not in towncrier order".format(version)]
        return ["{} section does not exactly match the towncrier draft".format(version)]
    exp_order, exp_groups = grouped(expected)
    act_order, act_groups = grouped(actual)
    errors = []
    for category in exp_order + [cat for cat in act_order if cat not in exp_groups]:
        exp_items = exp_groups.get(category, [])
        act_items = act_groups.get(category, [])
        if exp_items == act_items:
            continue
        here = (version, category)
        missing = list((Counter(exp_items) - Counter(act_items)).elements())
        extra = list((Counter(act_items) - Counter(exp_items)).elements())
        for entry in missing:
            hits = [spot for spot in found_at.get(entry, []) if spot != here]
            if hits:
                errors.append(wrong_section(hits[0], here, entry))
                continue
            partner = next((item for item in extra if edited(entry, item)), None)
            if partner is not None:
                extra.remove(partner)
                errors.append(
                    "{} has a stale or edited entry under {}: {}".format(
                        version, category, preview(entry)
                    )
                )
            else:
                errors.append(
                    "{} is missing an entry under {}: {}".format(version, category, preview(entry))
                )
        for entry in extra:
            owners = [spot for spot in expected_at.get(entry, []) if spot != here]
            if owners:
                errors.append(wrong_section(here, owners[0], entry))
            else:
                errors.append(
                    "{} has an extra entry under {}: {}".format(version, category, preview(entry))
                )
    if not errors:
        errors.append("{} section does not exactly match the towncrier draft".format(version))
    return errors


def changelog_errors(drafts, changelog):
    """Errors for CHANGELOG.md text against {version: rendered draft}."""
    errors = []
    headings = [match.group(1) for match in RELEASE.finditer(changelog)]
    for version, count in Counter(headings).items():
        if count > 1:
            errors.append("CHANGELOG.md repeats the [{}] heading".format(version))
    if headings != sorted(headings, key=version_key, reverse=True):
        errors.append("CHANGELOG.md releases are not newest first: {}".format(headings))
    for version in headings:
        if (
            version.startswith(SERIES + ".")
            and version_key(version) > version_key(SERIES_BASE)
            and version not in drafts
        ):
            errors.append(
                "CHANGELOG.md [{}] has no changelog.d/{} fragments".format(version, version)
            )
    expected_at = {}
    for version, text in drafts.items():
        for category, entry in section_entries(text):
            expected_at.setdefault(entry, []).append((version, category))
    found_at = {}
    for name, start, end in section_spans(changelog):
        for category, entry in section_entries(changelog[start:end]):
            found_at.setdefault(entry, []).append((name, category))
    for version in sorted(drafts, key=version_key, reverse=True):
        text = drafts[version]
        section = section_text(changelog, version)
        if section is None:
            errors.append("CHANGELOG.md has no [{}] section for its fragments".format(version))
            continue
        if section.strip() == text.strip():
            continue
        errors.extend(
            entry_errors(
                version,
                section_entries(text),
                section_entries(section),
                expected_at,
                found_at,
            )
        )
    unique = []
    for error in errors:
        if error not in unique:
            unique.append(error)
    return unique


def main():
    try:
        directory, template = read_config()
    except ValueError as exc:
        print(str(exc), file=sys.stderr)
        return 1
    fragments = ROOT / directory
    if not fragments.is_dir():
        print("missing fragment directory: {}".format(directory), file=sys.stderr)
        return 1
    errors, found = sources(fragments)
    if errors:
        print("\n".join(errors), file=sys.stderr)
        return 1
    if len(found[SERIES_BASE]) < MINIMUM_FRAGMENTS:
        print(
            "{} lost fragments ({})".format(directory, len(found[SERIES_BASE])),
            file=sys.stderr,
        )
        return 1
    changelog = CHANGELOG.read_text()
    dates = {version: date for version, date in RELEASE.findall(changelog)}
    drafts = {}
    errors = []
    for version, paths in found.items():
        if version not in dates:
            errors.append("CHANGELOG.md has no [{}] section for its fragments".format(version))
            continue
        text, error = draft(version, dates[version], paths, directory, template)
        if error is not None:
            errors.append("{}: {}".format(version, error))
            continue
        drafts[version] = text
    errors.extend(changelog_errors(drafts, changelog))
    if errors:
        print("\n".join(dict.fromkeys(errors)), file=sys.stderr)
        return 1
    counts = ", ".join(
        "{} {}".format(version, len(found[version]))
        for version in sorted(found, key=version_key, reverse=True)
    )
    print("towncrier drafts match CHANGELOG.md ({})".format(counts))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
