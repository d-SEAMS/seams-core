#!/usr/bin/env python3
# Copyright (c) d-SEAMS developers
"""Fill the cutoff-and-padding how-to from the padding CSV.

Every integer and decimal in the written page is a cell of
docs/data/recall-padding-validation.csv. That file is
data/recall-padding-validation.csv from science/interface-nucleation
at 6ec9c283, the table scripts/recall_padding.py writes. The Sphinx
configuration calls this before it reads the book.
"""

import csv
import pathlib
import re
import sys

ROOT = pathlib.Path(__file__).resolve().parents[1]
CSV = ROOT / "docs" / "data" / "recall-padding-validation.csv"
OUT = ROOT / "docs" / "source" / "howto" / "cutoff-padding-tables.rst"

SUMMARY = (
    "lane",
    "temperature_K",
    "frame",
    "molecules",
    "chill_ice",
    "dseams_off",
    "dseams_on",
    "fp_off",
    "fp_on",
    "kind_changed",
    "agree_off",
    "agree_on",
    "source",
)
LAYER = (
    "lane",
    "frame",
    "layer",
    "z_lo_A",
    "molecules",
    "dseams_off",
    "dseams_on",
    "fp_off",
    "fp_on",
    "kind_changed",
)
TIME = (
    "molecules",
    "padding",
    "rep",
    "seconds",
    "cores_online",
    "cores_affinity",
    "cores_used",
    "source",
)


def load(path):
    with path.open(newline="") as handle:
        return list(csv.DictReader(handle))


def cell(text):
    return " ".join(text.split())


def table(headers, records):
    lines = [".. list-table::", "   :header-rows: 1", ""]
    lines.append("   * - " + headers[0])
    for name in headers[1:]:
        lines.append("     - " + name)
    for row in records:
        lines.append("   * - " + cell(row[headers[0]]))
        for name in headers[1:]:
            lines.append("     - " + cell(row[name]))
    return lines


def render(rows):
    gate = [row for row in rows if row["record"] == "gate"]
    if len(gate) != 1:
        raise SystemExit("the CSV needs one gate row")
    summary = [row for row in rows if row["record"] == "summary"]
    layers = [
        row
        for row in rows
        if row["record"] == "layer" and row["kind_changed"] not in ("", "0")
    ]
    times = [row for row in rows if row["record"] == "time"]
    if not summary or not layers or not times:
        raise SystemExit("the CSV is missing summary, layer, or time rows")
    lines = [
        "The bulk-liquid gate",
        "--------------------",
        "",
        cell(gate[0]["source"]),
        "",
        "False positives off are ``{}`` and false positives on are ``{}``.".format(
            cell(gate[0]["fp_off"]), cell(gate[0]["fp_on"])
        ),
        "",
        "Summary rows",
        "------------",
        "",
        "Each row is one ``summary`` record. ``frame`` is the number of",
        "frames that row adds.",
        "",
    ]
    lines.extend(table(SUMMARY, summary))
    lines.extend(
        [
            "",
            "Layers whose label changed",
            "--------------------------",
            "",
            "These are the ``layer`` records whose ``kind_changed`` cell is",
            "not ``0``. ``z_lo_A`` is the lower edge of the slab.",
            "",
        ]
    )
    lines.extend(table(LAYER, layers))
    lines.extend(
        [
            "",
            "Classify timings",
            "----------------",
            "",
            "``padding`` is ``0`` when the open-face completion is off and",
            "``1`` when it is on. ``seconds`` is one replicate.",
            "",
        ]
    )
    lines.extend(table(TIME, times))
    lines.append("")
    return "\n".join(lines)


def quoted_numbers(page, raw):
    missing = []
    for token in re.findall(r"\d+(?:\.\d+)?", page):
        if token not in raw:
            missing.append(token)
    return missing


def main(argv):
    del argv
    if not CSV.is_file():
        print(f"missing {CSV.relative_to(ROOT)}", file=sys.stderr)
        return 1
    raw = CSV.read_text(encoding="utf-8")
    banned = tuple("/" + name for name in ("tmp", "home", "opt", "cvmfs"))
    for prefix in banned:
        if prefix in raw:
            print("CSV contains a machine path", file=sys.stderr)
            return 1
    page = render(load(CSV))
    missing = quoted_numbers(page, raw)
    if missing:
        print("page quotes numbers that are not CSV cells:", ", ".join(missing), file=sys.stderr)
        return 1
    for prefix in banned:
        if prefix in page:
            print("page contains a machine path", file=sys.stderr)
            return 1
    OUT.parent.mkdir(parents=True, exist_ok=True)
    OUT.write_text(page, encoding="utf-8")
    return 0


if __name__ == "__main__":
    sys.exit(main(sys.argv[1:]))
