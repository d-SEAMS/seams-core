#!/usr/bin/env python3
# Copyright (c) d-SEAMS developers
"""emit_api_symbols.splice must replace one symbol block and refuse any other layout.

Each case writes a small page to a temporary directory, runs splice on it,
and checks the result, or the error and that the page was left untouched.
"""

import sys
import tempfile
from pathlib import Path

sys.dont_write_bytecode = True
sys.path.insert(0, str(Path(__file__).resolve().parent))
import emit_api_symbols  # noqa: E402

BEGIN = emit_api_symbols.BEGIN
END = emit_api_symbols.END
BLOCK = "\n## Symbol reference\n\n### seams_new\n\nPurpose: New.\n"
WRAPPED = BEGIN + "\n" + BLOCK + END + "\n"
OLD = BEGIN + "\n\n## Symbol reference\n\n### seams_old\n\nPurpose: Old.\n" + END + "\n"


def refusal(begins, ends):
    return (
        "need exactly one API_SYMBOLS_BEGIN followed by one API_SYMBOLS_END, "
        "found API_SYMBOLS_BEGIN on line(s) {} and API_SYMBOLS_END on line(s) {}".format(begins, ends)
    )


CASES = [
    ("no markers", "# C\n\nIntro.\n", "# C\n\nIntro.\n\n" + WRAPPED, None),
    ("no markers and no final newline", "# C\n\nIntro.", "# C\n\nIntro.\n\n" + WRAPPED, None),
    ("one block", "# C\n\n" + OLD + "\n## After\n", "# C\n\n" + WRAPPED + "\n## After\n", None),
    ("one block at the end of the file", "# C\n\n" + OLD, "# C\n\n" + WRAPPED, None),
    ("BEGIN only", "# C\n\n" + BEGIN + "\n### seams_old\n", None, refusal("3", "none")),
    ("END only", "# C\n\n### seams_old\n" + END + "\n", None, refusal("none", "4")),
    ("two blocks", "# C\n\n" + OLD + "\nKept between blocks.\n\n" + OLD, None, refusal("3, 14", "10, 21")),
    ("END before BEGIN", "# C\n\n" + END + "\nText.\n" + BEGIN + "\n", None, refusal("5", "3")),
    ("two BEGIN and one END", "# C\n" + BEGIN + "\nA\n" + BEGIN + "\nB\n" + END + "\n", None, refusal("2, 4", "6")),
]


def run_case(directory, label, before, after, error):
    page = Path(directory) / "page.md"
    page.write_text(before)
    try:
        emit_api_symbols.splice(page, BLOCK)
    except ValueError as exc:
        if error is None:
            return ["{}: splice refused a valid page: {}".format(label, exc)]
        problems = []
        if str(exc) != "{}: {}".format(page, error):
            problems.append("{}: wrong error\n  got:  {}\n  want: {}: {}".format(label, exc, page, error))
        if page.read_text() != before:
            problems.append("{}: the refused page was rewritten".format(label))
        return problems
    if error is not None:
        return ["{}: splice accepted a page it must refuse".format(label)]
    got = page.read_text()
    if got != after:
        return ["{}: wrong page\n  got:  {!r}\n  want: {!r}".format(label, got, after)]
    emit_api_symbols.splice(page, BLOCK)
    if page.read_text() != got:
        return ["{}: a second splice changed the page".format(label)]
    return []


def main():
    failures = []
    with tempfile.TemporaryDirectory() as directory:
        for label, before, after, error in CASES:
            failures.extend(run_case(directory, label, before, after, error))
    if failures:
        print("\n".join(failures), file=sys.stderr)
        return 1
    print("emit_api_symbols.splice: {} cases".format(len(CASES)))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
