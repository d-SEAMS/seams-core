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


def python_symbols():
    """Names a caller gets from ``import dseams``, including Frame methods."""
    text = (ROOT / "bindings" / "python" / "dseams" / "__init__.py").read_text()
    names = set(re.findall(r"^def ([A-Za-z][A-Za-z0-9_]*)\(", text, re.M))
    names.update(re.findall(r"^([A-Za-z][A-Za-z0-9_]*) = ", text, re.M))
    names.update(re.findall(r"^class ([A-Za-z][A-Za-z0-9_]*)", text, re.M))
    core = (ROOT / "bindings" / "python" / "dseams" / "_core.cpp").read_text()
    head, _, tail = core.partition("nb::class_")
    names.update(re.findall(r'm\.def\(\s*"([A-Za-z][A-Za-z0-9_]*)"', head))
    for method in re.findall(r'\.def\(\s*"([A-Za-z][A-Za-z0-9_]*)"', tail):
        names.add(method)
        names.add("Frame." + method)
    user = (ROOT / "bindings" / "python" / "dseams" / "userframe.py").read_text()
    current = None
    for line in user.splitlines():
        matched = re.match(r"^class ([A-Za-z][A-Za-z0-9_]*)", line)
        if matched:
            current = matched.group(1)
            names.add("userframe." + current)
            continue
        matched = re.match(r"^def ([A-Za-z][A-Za-z0-9_]*)\(", line)
        if matched:
            current = None
            if not matched.group(1).startswith("_"):
                names.add("userframe." + matched.group(1))
            continue
        matched = re.match(r"^    def ([A-Za-z][A-Za-z0-9_]*)\(", line)
        if matched and current and not matched.group(1).startswith("_"):
            names.add("userframe.{}.{}".format(current, matched.group(1)))
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


def _julia_export(path):
    text = path.read_text()
    match = re.search(r"^export (.+)$", text, re.M)
    if match is None:
        return []
    return [part.strip() for part in match.group(1).split(",") if part.strip()]


def julia_symbols():
    """Exported names of both Julia packages, plus the Ice Ih path helper."""
    dseams = ROOT / "bindings" / "julia" / "Dseams" / "src" / "Dseams.jl"
    seams = ROOT / "bindings" / "julia" / "Seams" / "src" / "Seams.jl"
    names = set(_julia_export(dseams))
    if "ice_ih_path()" in dseams.read_text():
        names.add("ice_ih_path")
    seams_text = seams.read_text()
    for name in _julia_export(seams):
        names.add("Seams." + name)
    names.update("Seams." + item for item in re.findall(r"^const ([A-Z][A-Z0-9_]*) =", seams_text, re.M))
    return sorted(names)


def lua_symbols():
    """Names on the table returned by bindings/lua/dseams.lua."""
    text = (ROOT / "bindings" / "lua" / "dseams.lua").read_text()
    names = set(re.findall(r"^function dseams\.([A-Za-z_][A-Za-z0-9_]*)", text, re.M))
    names.update(re.findall(r"^dseams\.([A-Za-z_][A-Za-z0-9_]*)\s*=", text, re.M))
    names.update(re.findall(r'^dseams\["([^"]+)"\]\s*=', text, re.M))
    names.update(re.findall(r"^\s+([A-Z][A-Z0-9_]*) =", text, re.M))
    names.update({"core", "cloud.close", "cloud.box", "cloud.nop"})
    return sorted(names)


def rust_sys_symbols():
    """Public items of the generated seams-sys crate. They mirror include/seams.h."""
    text = (ROOT / "bindings" / "rust" / "seams-sys" / "src" / "lib.rs").read_text()
    names = set(re.findall(r"^pub (?:const|fn|struct|enum|type|union) ([A-Za-z0-9_]+)", text, re.M))
    return sorted(names)


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
    missing.extend("lua:" + name for name in _missing(lua_symbols(), blob))
    return missing
