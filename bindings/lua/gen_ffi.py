#!/usr/bin/env python3
# Copyright (c) d-SEAMS developers
# SPDX-License-Identifier: MIT
"""Write a LuaJIT FFI module from include/seams.h.

LuaJIT's cdef parser does not run the C preprocessor, so the generated
text is the declarations with the include, the export macro, and the
windows branch removed.
"""

import pathlib
import sys


def declarations(header):
    text = header
    start = text.find("typedef enum seams_status")
    end = text.rfind("#ifdef __cplusplus")
    if start < 0 or end < 0 or end < start:
        raise SystemExit("seams.h has no declaration span the FFI generator knows")
    body = text[start:end]
    body = body.replace("SEAMS_API ", "")
    lines = []
    for line in body.splitlines():
        stripped = line.strip()
        if stripped.startswith("#"):
            continue
        if stripped.startswith("/*") or stripped.startswith("*"):
            continue
        lines.append(line)
    return "\n".join(lines).strip() + "\n"


def render(header_text):
    cdef = declarations(header_text)
    return (
        "-- Generated from include/seams.h. Do not edit.\n"
        "-- Copyright (c) d-SEAMS developers\n"
        "local ffi = require('ffi')\n"
        "ffi.cdef[[\n"
        + cdef
        + "]]\n"
        "local name = os.getenv('SEAMS_LIB') or 'seams'\n"
        "return ffi.load(name)\n"
    )


def main(argv):
    header = pathlib.Path(argv[1]).read_text()
    pathlib.Path(argv[2]).write_text(render(header))
    return 0


if __name__ == "__main__":
    raise SystemExit(main(sys.argv))
