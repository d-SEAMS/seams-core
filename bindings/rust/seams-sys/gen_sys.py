#!/usr/bin/env python3
# Copyright (c) d-SEAMS developers
# SPDX-License-Identifier: MIT
"""Write seams-sys/src/lib.rs from include/seams.h.

LuaJIT is not involved. The Rust cdef is the C declarations with the
export macro and the comments removed. Run it again when seams.h grows:

    python3 bindings/rust/seams-sys/gen_sys.py include/seams.h bindings/rust/seams-sys/src/lib.rs
"""

import pathlib
import re
import shutil
import subprocess
import sys

RUST_KEYWORDS = {
    "type",
    "box",
    "ref",
    "move",
    "async",
    "await",
    "dyn",
    "impl",
    "fn",
    "mod",
    "pub",
    "use",
    "self",
    "super",
    "crate",
    "where",
    "as",
    "if",
    "else",
    "for",
    "in",
    "loop",
    "while",
    "return",
    "match",
    "const",
    "static",
    "mut",
    "let",
    "trait",
    "struct",
    "enum",
    "union",
    "extern",
    "unsafe",
}


def rust_name(name):
    if name in RUST_KEYWORDS:
        return "r#{}".format(name)
    return name


ENUM_VALUE = {
    "int": "c_int",
    "double": "f64",
    "size_t": "usize",
    "seams_status": "c_int",
}


def strip_comments(text):
    text = re.sub(r"/\*.*?\*/", "", text, flags=re.S)
    text = re.sub(r"//.*?$", "", text, flags=re.M)
    return text


def body_of(text):
    start = text.find("typedef enum seams_status")
    end = text.rfind("#ifdef __cplusplus")
    if start < 0 or end < start:
        raise SystemExit("seams.h has no declaration span the generator knows")
    return text[start:end]


def rust_base(ctype):
    ctype = " ".join(ctype.replace("const", "").split())
    if ctype in ("int", "seams_status"):
        return "c_int"
    if ctype == "double":
        return "f64"
    if ctype == "size_t":
        return "usize"
    if ctype == "char":
        return "c_char"
    if ctype == "void":
        return "c_void"
    if re.fullmatch(r"seams_[A-Za-z0-9_]+", ctype):
        return ctype
    raise SystemExit("unknown C type: {}".format(ctype))


def rust_type(raw):
    text = " ".join(raw.split())
    if text in ("void", ""):
        return None
    array = re.fullmatch(r"(.+?)\s*\[\d+\]", text)
    if array:
        text = array.group(1).strip() + " *"
    stars = text.count("*")
    base_src = text.replace("*", " ")
    const = "const" in base_src
    base = rust_base(base_src)
    if stars == 0:
        return base
    # A const before a single star qualifies that pointer. Two stars with
    # const qualify the inner pointer; the outer pointer is written.
    if stars == 1:
        return "*const {}".format(base) if const else "*mut {}".format(base)
    if stars == 2:
        inner = "*const {}".format(base) if const else "*mut {}".format(base)
        return "*mut {}".format(inner)
    raise SystemExit("too many stars: {}".format(raw))


def split_params(args):
    args = args.strip()
    if args == "" or args == "void":
        return []
    parts = []
    depth = 0
    current = []
    for ch in args:
        if ch == "(":
            depth += 1
        elif ch == ")":
            depth -= 1
        if ch == "," and depth == 0:
            parts.append("".join(current).strip())
            current = []
        else:
            current.append(ch)
    tail = "".join(current).strip()
    if tail:
        parts.append(tail)
    parsed = []
    for part in parts:
        if "[" in part:
            name = re.search(r"(\w+)\s*\[", part).group(1)
            ctype = part[: part.rfind(name)].strip() + " *"
        elif "*" in part:
            name = part.split("*")[-1].strip()
            ctype = part[: part.rfind("*") + 1].strip()
        else:
            name = part.split()[-1]
            ctype = part[: part.rfind(name)].strip()
        rust = rust_type(ctype)
        if rust is None:
            continue
        parsed.append((name, rust))
    return parsed


def enums_of(text):
    found = []
    for match in re.finditer(r"typedef\s+enum\s+(\w+)?\s*\{(.*?)\}\s*(\w+)\s*;", text, re.S):
        body = match.group(2)
        value = None
        for item in body.split(","):
            item = item.strip()
            if not item:
                continue
            if "=" in item:
                name, raw = item.split("=", 1)
                value = int(raw.strip(), 0)
                name = name.strip()
            else:
                name = item.split()[0]
                value = 0 if value is None else value + 1
            found.append((name, value))
    return found


def structs_of(text):
    return re.findall(r"typedef\s+struct\s+(\w+)\s+\1\s*;", text)


def functions_of(text):
    found = []
    for match in re.finditer(r"SEAMS_API\s+(.*?)\s*(\w+)\s*\((.*?)\)\s*;", text, re.S):
        ret = " ".join(match.group(1).split())
        name = match.group(2)
        args = " ".join(match.group(3).split())
        found.append((rust_type(ret), name, split_params(args)))
    return found


def render(header):
    text = body_of(strip_comments(header))
    lines = [
        "/* Copyright (c) d-SEAMS developers */",
        "/* Generated from include/seams.h by gen_sys.py. Do not edit. */",
        "#![allow(non_camel_case_types, non_snake_case)]",
        "",
        "use std::os::raw::{c_char, c_int, c_void};",
        "",
    ]
    for name, value in enums_of(text):
        lines.append("pub const {}: c_int = {};".format(name, value))
    lines.append("")
    for name in structs_of(text):
        lines.append("#[repr(C)]")
        lines.append("pub struct {} {{".format(name))
        lines.append("    _private: [u8; 0],")
        lines.append("}")
    lines.append("")
    lines.append('extern "C" {')
    for ret, name, params in functions_of(text):
        if params:
            lines.append("    pub fn {}(".format(name))
            for pname, ptype in params:
                lines.append("        {}: {},".format(rust_name(pname), ptype))
            lines.append("    ){};".format("" if ret is None else " -> {}".format(ret)))
        else:
            lines.append("    pub fn {}(){};".format(name, "" if ret is None else " -> {}".format(ret)))
    lines.append("}")
    lines.append("")
    return "\n".join(lines)


def main(argv):
    header = pathlib.Path(argv[1]).read_text()
    out = pathlib.Path(argv[2])
    out.write_text(render(header))
    rustfmt = shutil.which("rustfmt")
    if rustfmt:
        subprocess.check_call([rustfmt, "--edition", "2021", str(out)])
    return 0


if __name__ == "__main__":
    raise SystemExit(main(sys.argv))
