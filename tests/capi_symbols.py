#!/usr/bin/env python3
# Copyright (c) d-SEAMS developers
"""Header declarations, tests/capi_symbols.txt, and libseams exports must agree.

Usage: capi_symbols.py <seams.h> <capi_symbols.txt> <libseams.so>

Exit 0 when they agree, 1 on any difference, and 77 (a meson skip) when
neither nm nor readelf can list the library's exports.
"""

import re
import shutil
import subprocess
import sys
from pathlib import Path

DECL = re.compile(r"^SEAMS_API\b[^(;]*?\b(seams_[A-Za-z0-9_]+)\s*\(", re.M)
SKIP = 77


def header_names(path):
    text = Path(path).read_text()
    text = re.sub(r"/\*.*?\*/", "", text, flags=re.S)
    return sorted(set(DECL.findall(text)))


def listed_names(path):
    rows = Path(path).read_text().splitlines()
    return [row.strip() for row in rows if row.strip() and not row.startswith("#")]


def exported(lib):
    nm = shutil.which("nm")
    if nm is not None:
        done = subprocess.run(
            [nm, "-D", "--defined-only", lib], text=True, capture_output=True, check=True
        )
        names = []
        for row in done.stdout.splitlines():
            parts = row.split()
            if len(parts) >= 3 and parts[1] in "TtWwDdBbRrVvi":
                names.append(parts[2].split("@")[0])
        return sorted(set(names))
    readelf = shutil.which("readelf")
    if readelf is not None:
        done = subprocess.run(
            [readelf, "--dyn-syms", "-W", lib], text=True, capture_output=True, check=True
        )
        names = []
        for row in done.stdout.splitlines():
            parts = row.split()
            # Num: Value Size Type Bind Vis Ndx Name; UND entries are imports.
            if (len(parts) >= 8 and parts[0].endswith(":") and parts[6] != "UND"
                    and parts[3] in ("FUNC", "OBJECT", "IFUNC", "TLS")
                    and parts[4] in ("GLOBAL", "WEAK")):
                names.append(parts[7].split("@")[0])
        return sorted(set(names))
    return None


def compare(header, listed, exports):
    """Every difference between the three name sets, as messages."""
    errors = []
    if listed != sorted(listed) or len(listed) != len(set(listed)):
        errors.append("capi_symbols.txt is not sorted and unique")
    for name in sorted(set(header) - set(listed)):
        errors.append("declared in seams.h, missing from capi_symbols.txt: " + name)
    for name in sorted(set(listed) - set(header)):
        errors.append("listed, not declared in seams.h: " + name)
    for name in sorted(set(listed) - set(exports)):
        errors.append("listed, not exported by libseams: " + name)
    for name in sorted(set(exports) - set(listed)):
        errors.append("exported by libseams, not listed: " + name)
    return errors


def main(argv):
    if len(argv) != 4:
        print(__doc__, file=sys.stderr)
        return 2
    exports = exported(argv[3])
    if exports is None:
        print("neither nm nor readelf found; libseams exports not checked", file=sys.stderr)
        return SKIP
    listed = listed_names(argv[2])
    errors = compare(header_names(argv[1]), listed, exports)
    if errors:
        print("\n".join(errors), file=sys.stderr)
        return 1
    print("capi symbols ok ({} functions)".format(len(listed)))
    return 0


if __name__ == "__main__":
    raise SystemExit(main(sys.argv))
