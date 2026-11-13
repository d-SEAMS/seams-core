#!/usr/bin/env python3
"""The science gallery and the paper use one census."""

import re
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "docs"))

from write_science_gallery import CSV, load  # noqa: E402


def fail(message):
    print(message, file=sys.stderr)
    return 1


def main():
    rows = dict(((figure, quantity), count) for figure, quantity, count in load())
    tutorials = (ROOT / "docs" / "examples" / "science_tutorials.py").read_text()
    expected = {
        ("interface", "hexagonal"): 352,
        ("interface", "water"): 434,
        ("interface", "interfacial"): 49,
        ("interface", "cubic"): 1,
        ("brine", "hexagonal"): 741,
        ("brine", "water"): 621,
        ("brine", "interfacial"): 89,
        ("brine", "cubic"): 0,
        ("sii", "512"): 16,
        ("sii", "51264"): 8,
        ("sh", "512"): 12,
        ("sh", "sh"): 8,
        ("sh", "51268"): 4,
    }
    for key, count in expected.items():
        if rows.get(key) != count:
            return fail(f"{CSV.name} has {rows.get(key)} for {key}, wanted {count}")
        if key[0] == "interface" and key[1] in {"hexagonal", "water", "interfacial", "cubic"}:
            token = f"{key[1]} {count}"
            if token not in tutorials:
                return fail(f"science_tutorials.py is missing {token}")
    if "brine tend hexagonal 741 water 621 interfacial 89 cubic 0" not in tutorials:
        return fail("the brine tutorial census does not match the CSV")
    macros = (ROOT / "paper" / "tex" / "figures" / "data" / "macros.tex").read_text()
    for cage, count in (("512", 16), ("51264", 8)):
        match = re.search(
            rf"dsv@w@sII@cages@{cage}\\endcsname\{{\\num\{{(\d+)}}\}}",
            macros,
        )
        if match is None or int(match.group(1)) != count:
            return fail(f"paper sII {cage} is not {count}")
    if (ROOT / "docs" / "write_science_gallery.py").read_text() and __name__:
        from write_science_gallery import main as write_main

        if write_main(["--check"]) != 0:
            return 1
    print(f"figure census matches the tutorials and the sII paper macros ({len(rows)} rows)")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
