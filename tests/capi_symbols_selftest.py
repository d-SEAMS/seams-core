#!/usr/bin/env python3
# Copyright (c) d-SEAMS developers
"""capi_symbols.py must reject an unlisted export and a listed name the
library lost, and must skip rather than pass when it cannot read exports.

Usage: capi_symbols_selftest.py <seams.h> <capi_symbols.txt> <libseams.so>
"""

import shutil
import sys
from pathlib import Path

sys.dont_write_bytecode = True
sys.path.insert(0, str(Path(__file__).resolve().parent))
import capi_symbols  # noqa: E402

WHICH = shutil.which


def only_tool(tool):
    return lambda name: WHICH(name) if name == tool else None


def main(argv):
    if len(argv) != 4:
        print(__doc__, file=sys.stderr)
        return 2
    header_path, list_path, lib = argv[1:]
    header = capi_symbols.header_names(header_path)
    listed = capi_symbols.listed_names(list_path)
    failures = []

    readers = {}
    for tool in ("nm", "readelf"):
        if WHICH(tool) is None:
            continue
        shutil.which = only_tool(tool)
        readers[tool] = capi_symbols.exported(lib)
    shutil.which = WHICH
    if not readers:
        print("neither nm nor readelf found", file=sys.stderr)
        return capi_symbols.SKIP
    if len(set(map(tuple, readers.values()))) != 1:
        failures.append("nm and readelf disagree on the exports of " + lib)
    exports = next(iter(readers.values()))

    if capi_symbols.compare(header, listed, exports):
        failures.append("the checked-in list does not match the library to start with")

    added = "seams_unannounced_addition"
    errors = capi_symbols.compare(header, listed, exports + [added])
    if "exported by libseams, not listed: " + added not in errors:
        failures.append("an unlisted export passed")

    removed = listed[len(listed) // 2]
    errors = capi_symbols.compare(header, listed, [n for n in exports if n != removed])
    if "listed, not exported by libseams: " + removed not in errors:
        failures.append("a listed name missing from the library passed")

    shutil.which = lambda name: None
    status = capi_symbols.main(["capi_symbols.py", header_path, list_path, lib])
    shutil.which = WHICH
    if status != capi_symbols.SKIP:
        failures.append("without nm or readelf the check returned {}, not a skip".format(status))

    if failures:
        print("\n".join(failures), file=sys.stderr)
        return 1
    print("capi_symbols rejects additions and removals ({} readers)".format(len(readers)))
    return 0


if __name__ == "__main__":
    raise SystemExit(main(sys.argv))
