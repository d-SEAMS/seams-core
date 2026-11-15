#!/usr/bin/env python3
# Copyright (c) d-SEAMS developers
"""capi_docs.py must catch each kind of omission in seams.h and on the page.

Usage: capi_docs_selftest.py <seams.h> <capi_symbols.txt> <docs/api/c.md>

Each case removes one piece from an in-memory copy of the real header or
page and requires the matching message.
"""

import re
import sys
from pathlib import Path

sys.dont_write_bytecode = True
sys.path.insert(0, str(Path(__file__).resolve().parent))
import capi_docs  # noqa: E402

VICTIM = "seams_neighbours_copy"


def doc_comment(header, name):
    fn = re.search(r"/\*\*(?:(?!\*/).)*?\*/\s*SEAMS_API[^;]*?\b" + name + r"\s*\(", header, re.S)
    return fn.group(0)


def page_section(page, name):
    body = page.split(capi_docs.BEGIN, 1)[1]
    return re.search(r"^### " + name + r"\n.*?(?=^### |" + re.escape(capi_docs.END) + ")",
                     body, re.S | re.M).group(0)


def main(argv):
    if len(argv) != 4:
        print(__doc__, file=sys.stderr)
        return 2
    header = Path(argv[1]).read_text()
    names = capi_docs.listed_names(argv[2])
    page = Path(argv[3]).read_text()
    failures = []

    start = capi_docs.problems(header, names, page)
    if start:
        failures.append("the real header and page fail to start with:\n" + "\n".join(start))

    comment = doc_comment(header, VICTIM)
    section = page_section(page, VICTIM)
    cases = [
        ("no doc comment",
         header.replace(comment, re.sub(r"/\*\*.*?\*/\s*", "", comment, flags=re.S)), page, names,
         VICTIM + ": no doc comment above the declaration"),
        ("no @brief",
         header.replace(comment, re.sub(r" \* @brief .*\n", "", comment)), page, names,
         VICTIM + ": no @brief"),
        ("a missing @param",
         header.replace(comment, re.sub(r" \* @param cols_cap .*\n", "", comment)), page, names,
         VICTIM + ": no @param for cols_cap"),
        ("a stale @param",
         header.replace(comment, comment.replace("@param nnz ", "@param n_cols ")), page, names,
         VICTIM + ": @param n_cols names no parameter"),
        ("no @return",
         header.replace(comment, re.sub(r" \* @return .*\n(?: \*  .*\n)*", "", comment)), page, names,
         VICTIM + ": no @return"),
        ("no thread safety",
         header.replace(comment, comment.replace("@par Thread safety", "@par Threads")), page, names,
         VICTIM + ": no @par Thread safety"),
        ("no ownership",
         header.replace(comment, comment.replace("@par Ownership", "@par Notes")), page, names,
         VICTIM + ": takes or returns a pointer and has no @par Ownership"),
        ("no page section",
         header, page.replace(section, ""), names,
         VICTIM + ": no section on the C API page"),
        ("a stale page parameter",
         header, page.replace(section, section.replace("- `nnz`: ", "- `nnz`: Old text. ")), names,
         VICTIM + ": the C API page lacks - `nnz` as documented in seams.h"),
        ("an undocumented export",
         header + "\nSEAMS_API int seams_unannounced(void);\n", page, sorted(names + ["seams_unannounced"]),
         "seams_unannounced: no doc comment above the declaration"),
    ]
    for label, case_header, case_page, case_names, message in cases:
        if case_header == header and case_page == page:
            failures.append("the {} case changed nothing".format(label))
            continue
        errors = capi_docs.problems(case_header, case_names, case_page)
        if message not in errors:
            failures.append("{} passed; got {}".format(label, errors[:3]))

    if failures:
        print("\n".join(failures), file=sys.stderr)
        return 1
    print("capi_docs catches {} kinds of omission".format(len(cases)))
    return 0


if __name__ == "__main__":
    raise SystemExit(main(sys.argv))
