#!/usr/bin/env python3
# Copyright (c) d-SEAMS developers
"""Public C symbols in the installed headers must appear in the docs."""

import pathlib
import re

import check_examples

ROOT = pathlib.Path(__file__).resolve().parents[1]
HEADERS = (
    ROOT / "include" / "seams.h",
    ROOT / "bindings" / "cxx" / "seams_lang.h",
)


def public_c_symbols():
    """Names declared for callers in the installed C headers."""
    found = set()
    for header in HEADERS:
        text = header.read_text()
        found.update(re.findall(r"\bseams_[A-Za-z0-9_]+\b", text))
        found.update(re.findall(r"\bSEAMS_[A-Z0-9_]+\b", text))
    found.discard("SEAMS_API")
    return sorted(name for name in found if not name.endswith("_H_"))


def docs_text():
    """Prose plus the examples the docs pages include."""
    parts = []
    for folder in ("api", "engines", "tutorials"):
        for page in sorted((ROOT / "docs" / folder).glob("*.md")):
            parts.append(page.read_text())
    for page in sorted((ROOT / "docs" / "orgmode").rglob("*.org")):
        parts.append(page.read_text())
    errors = []
    corpus = check_examples.harvest(errors)
    if errors:
        raise RuntimeError("\n".join(errors))
    for chunks in corpus.by_doc.values():
        parts.extend(chunks)
    return "\n".join(parts)


def missing_public_c_symbols():
    """Public C symbols that no docs page names.

    >>> missing_public_c_symbols()
    []
    """
    blob = docs_text()
    missing = []
    for name in public_c_symbols():
        if re.search(r"(?<!\w)" + re.escape(name) + r"(?!\w)", blob) is None:
            missing.append(name)
    return missing
