#!/usr/bin/env python3
# Copyright (c) d-SEAMS developers
"""Each public symbol needs a purpose, parameters, units, and an example."""

import re
import sys
from pathlib import Path

import check_c_symbols
import run_snippets

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
    "lua": "lua.md",
}
EXAMPLES = {
    "c": "docs/examples/c_api_tour.c",
    "cxx": "docs/examples/minimal_cpp.cpp",
    "python": "docs/examples/minimal_py.py",
    "rust": "docs/examples/minimal_rust.rs",
    "fortran": "docs/examples/minimal_fortran.f90",
    "julia": "docs/examples/minimal_julia.jl",
    "r": "docs/examples/minimal_r.R",
    "lua": "docs/examples/minimal_lua.lua",
}
FIELDS = ("Purpose", "Parameters", "Returns", "Ownership", "Errors", "Units")


def groups():
    return {
        "c": check_c_symbols.c_header_symbols(),
        "cxx": check_c_symbols.cxx_symbols(),
        "python": check_c_symbols.python_symbols(),
        "rust": check_c_symbols.rust_symbols(),
        "fortran": check_c_symbols.fortran_symbols(),
        "julia": check_c_symbols.julia_symbols(),
        "r": check_c_symbols.r_symbols(),
        "lua": check_c_symbols.lua_symbols(),
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
        if not run_snippets.snippets_of(page):
            errors.append("{} has no executed example".format(page.name))
        for name in names:
            block = have.get(name)
            if block is None:
                errors.append("undocumented {} {}".format(lang, name))
                continue
            for label in FIELDS:
                if re.search(r"^" + label + r": \S.*$", block, re.M) is None:
                    errors.append("{}: {} needs {}".format(lang, name, label))
            example = re.search(r"^Example: `([^`]+)`", block, re.M)
            if example is None or not (ROOT / example.group(1)).is_file():
                errors.append("{}: {} needs an example file that exists".format(lang, name))
        for name in have:
            if name not in names:
                errors.append("documented symbol does not exist: {}: {}".format(lang, name))
    c_page = sections(reference_body(API / PAGES["c"]))
    for name in check_c_symbols.rust_sys_symbols():
        block = c_page.get(name)
        if block is None:
            errors.append("undocumented rust-sys {}".format(name))
            continue
        for label in FIELDS:
            if re.search(r"^" + label + r": \S.*$", block, re.M) is None:
                errors.append("rust-sys {} needs {} on the C page".format(name, label))
    return errors


def main():
    errors = problems()
    if errors:
        print("\n".join(errors), file=sys.stderr)
        return 1
    count = sum(len(names) for names in groups().values())
    count += len(check_c_symbols.rust_sys_symbols())
    undocumented = [item for item in errors if item.startswith("undocumented")]
    print("api coverage ok ({} symbols, {} undocumented)".format(count, len(undocumented)))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
