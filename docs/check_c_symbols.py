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
    found = set(re.findall(r"\bseams_[A-Za-z0-9_]+\b", text))
    found.update(re.findall(r"\bSEAMS_[A-Z0-9_]+\b", text))
    return sorted(name for name in found if not name.endswith("_H_"))


def c_header_symbols():
    """Functions and constants in include/seams.h."""
    text = (ROOT / "include" / "seams.h").read_text()
    found = set(re.findall(r"\bseams_[A-Za-z0-9_]+\b", text))
    found.update(re.findall(r"\bSEAMS_[A-Z0-9_]+\b", text))
    found.discard("SEAMS_API")
    return sorted(name for name in found if not name.endswith("_H_"))


def python_symbols():
    """Names a caller gets from ``import dseams``."""
    text = (ROOT / "bindings" / "python" / "dseams" / "__init__.py").read_text()
    names = set(re.findall(r"^def ([A-Za-z][A-Za-z0-9_]*)\(", text, re.M))
    names.update(re.findall(r"^([A-Za-z][A-Za-z0-9_]*) = ", text, re.M))
    core = (ROOT / "bindings" / "python" / "dseams" / "_core.cpp").read_text()
    names.update(re.findall(r'\.def\(\s*"([A-Za-z][A-Za-z0-9_]*)"', core))
    names.update(re.findall(r'\.def\(\s*"([A-Za-z][A-Za-z0-9_]*)"', core))
    names.update(re.findall(r'm\.def\(\s*"([A-Za-z][A-Za-z0-9_]*)"', core))
    return sorted(names)


def rust_symbols():
    """Public items of the safe ``seams`` crate."""
    text = (ROOT / "bindings" / "rust" / "seams" / "src" / "lib.rs").read_text()
    names = set(re.findall(r"^pub (?:struct|fn|enum|const|type) ([A-Za-z0-9_]+)", text, re.M))
    names.update(re.findall(r"^    pub fn ([A-Za-z0-9_]+)", text, re.M))
    return sorted(names)


def fortran_symbols():
    """Parameters and procedures exported by the Fortran module."""
    text = (ROOT / "bindings" / "fortran" / "seams.f90").read_text()
    names = set(re.findall(r"parameter :: ([A-Z0-9_]+)", text))
    names.update(re.findall(r"(?:function|subroutine) ([a-z0-9_]+)\(", text))
    return sorted(names)


def julia_symbols():
    """Names on the Dseams export line."""
    text = (ROOT / "bindings" / "julia" / "Dseams" / "src" / "Dseams.jl").read_text()
    match = re.search(r"^export (.+)$", text, re.M)
    if match is None:
        return []
    return sorted(part.strip() for part in match.group(1).split(",") if part.strip())


def r_symbols():
    """Functions marked @export in the R package."""
    text = (ROOT / "bindings" / "r" / "dseams" / "R" / "dseams.R").read_text()
    return sorted(re.findall(r"#' @export\n([A-Za-z0-9_.]+) <- function", text))


def missing_public_symbols():
    """Undocumented public names in each shipped language front.

    >>> missing_public_symbols()
    []
    """
    blob = docs_text()
    missing = []
    missing.extend("c:" + name for name in _missing(c_header_symbols(), blob))
    missing.extend("cxx:" + name for name in _missing(cxx_symbols(), blob))
    missing.extend("python:" + name for name in _missing(python_symbols(), blob))
    missing.extend("rust:" + name for name in _missing(rust_symbols(), blob))
    missing.extend("fortran:" + name for name in _missing(fortran_symbols(), blob))
    missing.extend("julia:" + name for name in _missing(julia_symbols(), blob))
    missing.extend("r:" + name for name in _missing(r_symbols(), blob))
    return missing
