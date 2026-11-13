#!/usr/bin/env python3
"""Fail when an exported binding symbol has no docstring.

Julia exports, R NAMESPACE exports, Lua ``dseams.*`` functions, Python
public functions, and ``seams_lang.h`` declarations are all checked.
Each language's docs must also state the Ice Ih census: 128 hexagonal
waters and 256 six-rings.
"""

import ast
import re
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
MISSING = []


def fail(message):
    MISSING.append(message)


def _doc_blocks(text):
    return re.findall(r'"""(.*?)"""', text, flags=re.S)


def check_julia(path):
    text = path.read_text()
    export_line = next(
        (line for line in text.splitlines() if line.startswith("export ")),
        "",
    )
    names = [part.strip() for part in export_line[len("export ") :].split(",") if part.strip()]
    blocks = _doc_blocks(text)
    for name in names:
        bare = name.rstrip("!")
        if not any(re.search(r"\b" + re.escape(bare) + r"\b", block) for block in blocks):
            fail("{} export {} has no docstring".format(path, name))
    if "128" not in text or "256" not in text:
        fail("{} is missing the Ice Ih census (128 hexagonal, 256 six-rings)".format(path))


def check_r():
    namespace = (ROOT / "bindings/r/dseams/NAMESPACE").read_text()
    names = re.findall(r"^export\(([^)]+)\)", namespace, flags=re.M)
    pages = list((ROOT / "bindings/r/dseams/man").glob("*.Rd"))
    aliases = []
    for page in pages:
        aliases.extend(re.findall(r"\\alias\{([^}]+)\}", page.read_text()))
    for name in names:
        if name not in aliases:
            fail("R export {} has no Rd alias".format(name))
    ice = (ROOT / "bindings/r/dseams/man/chill_plus.Rd").read_text()
    if "128" not in ice or "256" not in ice:
        fail("R chill_plus.Rd is missing the Ice Ih census")


def check_lua():
    path = ROOT / "bindings/lua/dseams.lua"
    lines = path.read_text().splitlines()
    for index, line in enumerate(lines):
        stripped = line.strip()
        is_fun = stripped.startswith("function dseams.")
        is_alias = stripped.startswith("dseams.") and "=" in stripped.split("--")[0]
        if not is_fun and not is_alias:
            continue
        prev = ""
        for earlier in range(index - 1, -1, -1):
            if lines[earlier].strip():
                prev = lines[earlier].strip()
                break
        if not prev.startswith("---"):
            fail("Lua {} has no LDoc comment".format(stripped.split("(")[0]))
    text = "\n".join(lines)
    if "128" not in text or "256" not in text:
        fail("Lua module doc is missing the Ice Ih census")


def _public_functions(path):
    tree = ast.parse(path.read_text(), filename=str(path))
    found = []
    for node in tree.body:
        if isinstance(node, ast.FunctionDef) and not node.name.startswith("_"):
            found.append(node)
        elif isinstance(node, ast.ClassDef) and not node.name.startswith("_"):
            if ast.get_docstring(node) is None:
                fail("{} class {} has no docstring".format(path.name, node.name))
            for item in node.body:
                if isinstance(item, ast.FunctionDef) and not item.name.startswith("_"):
                    found.append(item)
    for node in found:
        if ast.get_docstring(node) is None:
            fail("{} function {} has no docstring".format(path.name, node.name))


def check_python():
    init = ROOT / "bindings/python/dseams/__init__.py"
    frames = ROOT / "bindings/python/frames.py"
    _public_functions(init)
    _public_functions(frames)
    text = init.read_text()
    if "128" not in text or "256" not in text:
        fail("Python module doc is missing the Ice Ih census")
    page = ROOT / "docs/source/lang/python.rst"
    page_text = page.read_text() if page.is_file() else ""
    if "core_counts" not in page_text or "128" not in page_text or "256" not in page_text:
        fail("Python API page is missing the Ice Ih doctest")


def check_cxx():
    path = ROOT / "bindings/cxx/seams_lang.h"
    lines = path.read_text().splitlines()
    hidden = False
    for index, line in enumerate(lines):
        if "DOXYGEN_SHOULD_SKIP_THIS" in line and line.strip().startswith("#ifndef"):
            hidden = True
            continue
        if hidden and line.strip().startswith("#endif"):
            hidden = False
            continue
        if hidden or line.strip().startswith(("*", "/", "#")):
            continue
        if not re.search(r"\bseams_[A-Za-z0-9_]+\s*\(", line):
            continue
        window = "\n".join(lines[max(0, index - 6) : index])
        if "/**" not in window:
            fail("C++ {} has no Doxygen comment".format(line.strip()))
    text = path.read_text() + (ROOT / "docs/examples/lang_ice.c").read_text()
    if "128" not in text or "256" not in text:
        fail("C++ docs are missing the Ice Ih census")


def main():
    check_julia(ROOT / "bindings/julia/Dseams/src/Dseams.jl")
    check_julia(ROOT / "bindings/julia/Seams/src/Seams.jl")
    check_r()
    check_lua()
    check_python()
    check_cxx()
    if MISSING:
        sys.stderr.write("\n".join(MISSING) + "\n")
        return 1
    print("binding docstrings cover every export")
    return 0


if __name__ == "__main__":
    sys.exit(main())
