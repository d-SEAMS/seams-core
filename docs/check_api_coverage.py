#!/usr/bin/env python3
# Copyright (c) d-SEAMS developers
"""Each public symbol needs a purpose, parameters, units, and an example."""

import re
import sys
from pathlib import Path

import check_c_symbols

ROOT = Path(__file__).resolve().parents[1]
API = ROOT / "docs" / "api"
BEGIN = "<!-- API_SYMBOLS_BEGIN -->"
END = "<!-- API_SYMBOLS_END -->"
PAGES = {
    "c": "c.md",
    "cxx": "cpp.md",
    "python": "python.md",
    "rust": "rust.md",
    "fortran": "fortran.md",
    "julia": "julia.md",
    "r": "r.md",
}
EXAMPLES = {
    "c": "docs/examples/minimal_c.c",
    "cxx": "docs/examples/minimal_cpp.cpp",
    "python": "docs/examples/minimal_py.py",
    "rust": "docs/examples/minimal_rust.rs",
    "fortran": "docs/examples/minimal_fortran.f90",
    "julia": "docs/examples/minimal_julia.jl",
    "r": "docs/examples/minimal_r.R",
}


def groups():
    return {
        "c": check_c_symbols.c_header_symbols(),
        "cxx": check_c_symbols.cxx_symbols(),
        "python": check_c_symbols.python_symbols(),
        "rust": check_c_symbols.rust_symbols(),
        "fortran": check_c_symbols.fortran_symbols(),
        "julia": check_c_symbols.julia_symbols(),
        "r": check_c_symbols.r_symbols(),
    }


def reference_body(path):
    text = path.read_text()
    start = text.find(BEGIN)
    stop = text.find(END)
    if start < 0 or stop < start:
        return ""
    return text[start:stop]


def sections(body):
    found = {}
    chunks = re.split(r"^### ", body, flags=re.M)
    for chunk in chunks[1:]:
        line, _, rest = chunk.partition("\n")
        found[line.strip()] = rest
    return found


def problems():
    errors = []
    for lang, names in groups().items():
        page = API / PAGES[lang]
        if not page.is_file():
            errors.append("missing {}".format(page))
            continue
        body = reference_body(page)
        if not body:
            errors.append("{} has no symbol reference block".format(page.name))
            continue
        have = sections(body)
        example = EXAMPLES[lang]
        if not (ROOT / example).is_file():
            errors.append("missing example {}".format(example))
        for name in names:
            block = have.get(name)
            if block is None:
                errors.append("{}: {} has no section".format(lang, name))
                continue
            purpose = re.search(r"^Purpose: \S.*$", block, re.M)
            parameters = re.search(r"^Parameters: \S.*$", block, re.M)
            units = re.search(r"^Units: \S.*$", block, re.M)
            if purpose is None:
                errors.append("{}: {} needs a one-line purpose".format(lang, name))
            if parameters is None:
                errors.append("{}: {} needs parameters".format(lang, name))
            if units is None:
                errors.append("{}: {} needs units".format(lang, name))
            if example not in block:
                errors.append("{}: {} needs a link to {}".format(lang, name, example))
    return errors


def main():
    errors = problems()
    if errors:
        print("\n".join(errors), file=sys.stderr)
        return 1
    count = sum(len(names) for names in groups().values())
    print("api coverage ok ({} symbols)".format(count))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
