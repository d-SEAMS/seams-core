#!/usr/bin/env python3
# Copyright (c) d-SEAMS developers
"""Build the rc3 sitting tables from benchmarks/rc3-sitting.csv.

The median of an even sample is the mean of the two central values.
The interquartile range is the type-7 quartile gap, the same linear
position q*(n-1) that scripts/h2h_protocol.py uses. --check recomputes
both from samples_ms and refuses a row that does not match the stored
summary. The Sphinx configuration calls write_tables before it reads
the book, so the numbers are not copied by hand.
"""

import csv
import pathlib
import sys

ROOT = pathlib.Path(__file__).resolve().parents[1]
CSV = ROOT / "benchmarks" / "rc3-sitting.csv"
OUT = ROOT / "docs" / "source" / "explanation" / "performance-tables.rst"

ORDER = ("neighbours", "q4", "q6", "CHILL+", "rings", "cages", "PTM", "king")


def quantile_linear(sorted_vals, q):
    if len(sorted_vals) == 1:
        return sorted_vals[0]
    pos = q * (len(sorted_vals) - 1)
    lo = int(pos)
    frac = pos - lo
    if lo + 1 >= len(sorted_vals):
        return sorted_vals[-1]
    return sorted_vals[lo] * (1.0 - frac) + sorted_vals[lo + 1] * frac


def summarize(samples):
    ordered = sorted(float(v) for v in samples)
    count = len(ordered)
    mid = count // 2
    if count % 2 == 0:
        median = 0.5 * (ordered[mid - 1] + ordered[mid])
    else:
        median = ordered[mid]
    return count, median, quantile_linear(ordered, 0.75) - quantile_linear(ordered, 0.25)


def load():
    with CSV.open(newline="") as handle:
        return list(csv.DictReader(handle))


def check(rows):
    errors = []
    for index, row in enumerate(rows, start=2):
        samples = [float(part) for part in row["samples_ms"].split()]
        count, median, iqr = summarize(samples)
        if int(row["n"]) != count:
            errors.append("line {}: n {} != {} samples".format(index, row["n"], count))
        # samples_ms is already rounded to 0.001, so the recomputed
        # summary can move by one stored unit from the original.
        if abs(median - float(row["median_ms"])) > 1.5e-3:
            errors.append("line {}: median {} != {:.6f}".format(index, row["median_ms"], median))
        if abs(iqr - float(row["iqr_ms"])) > 1.5e-3:
            errors.append("line {}: iqr {} != {:.6f}".format(index, row["iqr_ms"], iqr))
    return errors


def cell(text):
    return text.replace("\n", " ")


def write_tables(rows, path):
    grouped = {}
    for row in rows:
        grouped.setdefault(row["operation"], []).append(row)
    names = [name for name in ORDER if name in grouped]
    names.extend(sorted(set(grouped) - set(names)))
    lines = [
        ".. list-table::",
        "   :header-rows: 0",
        "",
        "   * - These tables are written by ``docs/render_sitting.py`` from",
        "       ``benchmarks/rc3-sitting.csv``. Edit the CSV, not this file.",
        "",
    ]
    header = ("code", "version", "molecules", "threads", "cores", "n", "median_ms", "iqr_ms")
    for name in names:
        block = sorted(
            grouped[name],
            key=lambda row: (int(row["molecules"]), int(row["threads"]), row["code"], row["version"]),
        )
        lines.append(name)
        lines.append("^" * len(name))
        lines.append("")
        lines.append(".. csv-table::")
        lines.append("   :header: " + ", ".join(header))
        lines.append("")
        for row in block:
            lines.append("   " + ", ".join(cell(row[key]) for key in header))
        lines.append("")
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text("\n".join(lines) + "\n")


def main(argv):
    rows = load()
    errors = check(rows)
    if errors:
        print("\n".join(errors), file=sys.stderr)
        return 1
    if "--check" in argv:
        print("sitting csv ok", len(rows), "rows")
        return 0
    write_tables(rows, OUT)
    print("wrote", OUT.relative_to(ROOT), len(rows), "rows")
    return 0


if __name__ == "__main__":
    raise SystemExit(main(sys.argv[1:]))
