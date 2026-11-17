#!/usr/bin/env python3
# Copyright (c) d-SEAMS developers
"""Public C symbols in the installed headers must appear in the docs."""

import pathlib
import re

import check_examples

ROOT = pathlib.Path(__file__).resolve().parents[1]
HEADERS = (
    ROOT / "include" / "seams.h",
    ROOT / "include" / "seams_lammps_box.h",
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
    for header in HEADERS:
        found.difference_update(re.findall(r"#ifndef ([A-Z0-9_]+)", header.read_text()))
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
    return _missing(public_c_symbols(), docs_text())


def _missing(names, blob):
    missing = []
    for name in names:
        if re.search(r"(?<!\w)" + re.escape(name) + r"(?!\w)", blob) is None:
            missing.append(name)
    return missing


def cxx_symbols():
    """Functions and constants in the C++ language header."""
    text = (ROOT / "bindings" / "cxx" / "seams_lang.h").read_text()
    # @file seams_lang.h is a filename, not a callable.
    text = re.sub(r"\bseams_[A-Za-z0-9_]+\.h\b", " ", text)
    found = set(re.findall(r"\bseams_[A-Za-z0-9_]+\b", text))
    found.update(re.findall(r"\bSEAMS_[A-Z0-9_]+\b", text))
    return sorted(name for name in found if not name.endswith("_H_"))


def _header_names(text):
    found = set(re.findall(r"\bseams_[A-Za-z0-9_]+\b", text))
    found.update(re.findall(r"\bSEAMS_[A-Z0-9_]+\b", text))
    found.discard("SEAMS_API")
    found.difference_update(re.findall(r"#ifndef ([A-Z0-9_]+)", text))
    return {name for name in found if not name.endswith("_H_")}


def c_header_symbols():
    """Functions and constants in the installed C headers."""
    found = set()
    for name in ("seams.h", "seams_lammps_box.h"):
        found.update(_header_names((ROOT / "include" / name).read_text()))
    return sorted(found)


def missing_public_symbols():
    """Undocumented public names in the C API and the C++ language header.

    The binding repositories check their own languages.

    >>> missing_public_symbols()
    []
    """
    blob = docs_text()
    missing = []
    missing.extend("c:" + name for name in _missing(c_header_symbols(), blob))
    missing.extend("cxx:" + name for name in _missing(cxx_symbols(), blob))
    return missing
