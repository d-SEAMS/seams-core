#!/usr/bin/env python3
# Copyright (c) d-SEAMS developers
"""Header declarations, tests/capi_symbols.txt, and libseams exports must agree.

Usage: capi_symbols.py <seams.h> <capi_symbols.txt> <libseams.so>
"""

import re
import shutil
import subprocess
import sys
from pathlib import Path

DECL = re.compile(r"^SEAMS_API\b[^(;]*?\b(seams_[A-Za-z0-9_]+)\s*\(", re.M)


def header_names(path):
    text = Path(path).read_text()
    text = re.sub(r"/\*.*?\*/", "", text, flags=re.S)
    return sorted(set(DECL.findall(text)))


def listed_names(path):
    rows = Path(path).read_text().splitlines()
    return [row.strip() for row in rows if row.strip() and not row.startswith("#")]


def exported(lib):
    nm = shutil.which("nm")
    if nm is None:
        return None
    done = subprocess.run(
        [nm, "-D", "--defined-only", lib], text=True, capture_output=True, check=True
    )
    names = []
    for row in done.stdout.splitlines():
        parts = row.split()
        if len(parts) >= 3 and parts[1] in "TtWwDdBbRrVvi":
            names.append(parts[2].split("@")[0])
    return sorted(set(names))


def main(argv):
    if len(argv) != 4:
        print(__doc__, file=sys.stderr)
        return 2
    header = header_names(argv[1])
    listed = listed_names(argv[2])
    errors = []
    if listed != sorted(listed) or len(listed) != len(set(listed)):
        errors.append("capi_symbols.txt is not sorted and unique")
    for name in sorted(set(header) - set(listed)):
        errors.append("declared in seams.h, missing from capi_symbols.txt: " + name)
    for name in sorted(set(listed) - set(header)):
        errors.append("listed, not declared in seams.h: " + name)
    exports = exported(argv[3])
    if exports is None:
        print("nm not found; library exports not checked", file=sys.stderr)
    else:
        for name in sorted(set(listed) - set(exports)):
            errors.append("listed, not exported by libseams: " + name)
        for name in sorted(set(exports) - set(listed)):
            errors.append("exported by libseams, not listed: " + name)
    if errors:
        print("\n".join(errors), file=sys.stderr)
        return 1
    print("capi symbols ok ({} functions)".format(len(listed)))
    return 0


if __name__ == "__main__":
    raise SystemExit(main(sys.argv))
