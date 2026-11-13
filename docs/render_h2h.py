#!/usr/bin/env python3
# Copyright (c) d-SEAMS developers
"""Build the performance head-to-head from the benchmark CSVs.

benchmarks/one-sitting-v3.csv is included when that file is in the
tree, so a later one-sitting measurement shows up on the next render.
--check writes the tables and refuses a cell the CSV does not hold.
"""

import csv
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
OUT = ROOT / "docs" / "source" / "explanation" / "performance-h2h-tables.rst"
SITTING_V3 = ROOT / "benchmarks" / "one-sitting-v3.csv"
NUMERIC = {"median_ms", "iqr_ms", "n"}
HEADER = {
    "code": "code",
    "version": "version",
    "operation": "operation",
    "molecules": "molecules",
    "threads": "threads",
    "phase": "phase",
    "cores": "cores",
    "median_ms": "median ms",
    "iqr_ms": "iqr ms",
    "n": "n",
}


def read_csv(path):
    with path.open(newline="") as handle:
        return list(csv.DictReader(handle))


def keep(row, operation, codes, molecules=None, threads=None):
    if row.get("operation") != operation:
        return False
    if codes and row.get("code") not in codes:
        return False
    if molecules is not None and row.get("molecules") not in molecules:
        return False
    if threads is not None and row.get("threads") not in threads:
        return False
    return True


def dedupe(rows, keys):
    """Drop a later row only when every kept field matches, including the median."""
    seen = set()
    kept = []
    for row in rows:
        token = tuple(row.get(key, "") for key in keys)
        if token in seen:
            continue
        seen.add(token)
        kept.append(row)
    return kept


def cell(row, key):
    raw = row.get(key, "")
    if key in {"median_ms", "iqr_ms"} and raw:
        return "{:.3f}".format(float(raw))
    if key == "n" and raw:
        return str(int(float(raw)))
    return raw


def table(title, csv_path, rows, columns):
    present = [key for key in columns if any(row.get(key, "") for row in rows)]
    lines = [
        title,
        "-" * len(title),
        "",
        "The rows are {}.".format(csv_path.relative_to(ROOT).as_posix()),
        "",
        ".. csv-table::",
        "   :header: {}".format(", ".join(HEADER[key] for key in present)),
        "",
    ]
    for row in rows:
        lines.append("   {}".format(", ".join(cell(row, key) for key in present)))
    lines.append("")
    return "\n".join(lines), present


def sections():
    found = []
    rc3 = ROOT / "benchmarks" / "rc3-sitting.csv"
    large = ROOT / "benchmarks" / "paper-cells" / "sitting-large.csv"
    watch = ROOT / "benchmarks" / "h2h-watch-v3.csv"
    sqq = ROOT / "benchmarks" / "h2h-watch-v6-sqq.csv"
    rings_csv = ROOT / "benchmarks" / "cutoff-rings.csv"
    rc3_rows = read_csv(rc3)
    large_rows = read_csv(large)
    jobs = (
        (
            "Neighbours",
            large,
            [row for row in large_rows if keep(
                row, "neighbours", {"d-SEAMS", "freud", "mdapy", "OVITO"},
                {"884736"}, {"1", "8"})],
            ("code", "version", "phase", "operation", "molecules", "threads", "median_ms", "n"),
        ),
        (
            "Rings",
            rc3,
            [row for row in rc3_rows if keep(
                row, "rings", {"d-SEAMS", "OVITO", "mdapy"},
                {"110592", "884736"}, {"1", "8"})],
            ("code", "version", "operation", "molecules", "threads", "median_ms", "iqr_ms", "n"),
        ),
        (
            "Rings against cycless",
            watch,
            dedupe(
                [row for row in read_csv(watch) if keep(
                    row, "rings", {"d-SEAMS", "cycless"}, {"1000", "110592"})],
                ("code", "operation", "molecules", "median_ms"),
            ),
            ("code", "operation", "molecules", "median_ms", "n"),
        ),
        (
            "Rings against SQQ",
            sqq,
            read_csv(sqq),
            ("code", "operation", "molecules", "threads", "median_ms", "n"),
        ),
        (
            "q6",
            rc3,
            [row for row in rc3_rows if keep(
                row, "q6", {"d-SEAMS", "freud", "mdapy"}, {"884736"}, {"1", "8"})],
            ("code", "version", "operation", "molecules", "threads", "median_ms", "n"),
        ),
        (
            "q6 against OVITO",
            large,
            [row for row in large_rows if keep(
                row, "q6", {"d-SEAMS", "OVITO"}, {"884736"}, {"1", "8"})],
            ("code", "version", "phase", "operation", "molecules", "threads", "median_ms", "n"),
        ),
        (
            "CHILL+",
            rc3,
            [row for row in rc3_rows if keep(row, "CHILL+", {"d-SEAMS", "OVITO", "mdapy"})],
            ("code", "version", "operation", "molecules", "threads", "median_ms", "n"),
        ),
        (
            "Cages",
            rc3,
            [row for row in rc3_rows if keep(row, "cages", {"d-SEAMS", "ClavisCage"})],
            ("code", "version", "operation", "molecules", "threads", "median_ms", "n"),
        ),
        (
            "PTM",
            rc3,
            [row for row in rc3_rows if keep(row, "PTM", {"d-SEAMS", "OVITO"})],
            ("code", "version", "operation", "molecules", "threads", "median_ms", "n"),
        ),
        (
            "Linkcell scaling",
            rings_csv,
            [row for row in read_csv(rings_csv) if row.get("operation") == "total"],
            ("code", "operation", "threads", "median_ms", "iqr_ms", "n"),
        ),
    )
    for title, path, rows, columns in jobs:
        if not rows:
            raise SystemExit("no rows for {}".format(title))
        found.append((path, *table(title, path, rows, columns)))
    if SITTING_V3.is_file():
        rows = [
            row for row in read_csv(SITTING_V3)
            if row.get("operation") in {
                "neighbours", "rings", "q6", "CHILL+", "cages", "PTM", "ptm",
            }
        ]
        columns = tuple(
            key for key in (
                "code", "version", "phase", "operation", "molecules", "threads",
                "cores", "median_ms", "iqr_ms", "n",
            )
            if rows and key in rows[0]
        )
        if rows and columns:
            found.append((SITTING_V3, *table("One sitting v3", SITTING_V3, rows, columns)))
    return found


def render():
    blocks = [
        ".. Generated by docs/render_h2h.py. Edit the CSVs, then rerender.",
        "",
    ]
    built = sections()
    for _path, text, _columns in built:
        blocks.append(text)
    if not SITTING_V3.is_file():
        blocks.append(
            "benchmarks/one-sitting-v3.csv is the next one-sitting file. "
            "This page adds its table when the file is in the tree."
        )
        blocks.append("")
    return "\n".join(blocks).rstrip() + "\n", built


def close_enough(shown, raw):
    if shown == raw:
        return True
    try:
        left = float(shown)
        right = float(raw)
    except ValueError:
        return False
    decimals = len(shown.split(".", 1)[1]) if "." in shown else 0
    return abs(left - right) <= 1.5 * 10 ** (-decimals)


def check(text, built):
    errors = []
    for path, block, columns in built:
        present = [key for key in columns if any(row.get(key, "") for row in [])]
        # Recover the columns that table() kept by reading the header line.
        header = None
        body = []
        for line in block.splitlines():
            if line.startswith("   :header:"):
                header = [part.strip() for part in line.split(":", 2)[-1].split(",")]
            elif line.startswith("   ") and header and not line.startswith("   :"):
                body.append([part.strip() for part in line.strip().split(",")])
        inverse = {label: key for key, label in HEADER.items()}
        keys = [inverse[label] for label in header]
        csv_rows = read_csv(path)
        for cells in body:
            record = dict(zip(keys, cells))
            matches = [
                row for row in csv_rows
                if all(row.get(key, "") == record[key] or (
                    key in NUMERIC and close_enough(record[key], row.get(key, ""))
                ) for key in keys if record.get(key, "") != "")
            ]
            # A numeric cell must match the same CSV row as the labels.
            labels = [key for key in keys if key not in NUMERIC]
            pool = [
                row for row in csv_rows
                if all(row.get(key, "") == record[key] for key in labels)
            ]
            if not pool:
                errors.append("{} has no row {}".format(path.name, record))
                continue
            for key in keys:
                if key not in NUMERIC:
                    continue
                if not any(close_enough(record[key], row.get(key, "")) for row in pool):
                    errors.append("{} {} is {}, CSV differs".format(path.name, key, record[key]))
        del matches
        del present
    if "one-sitting-v3.csv" not in text:
        errors.append("the page does not name one-sitting-v3.csv")
    return errors


def main(argv):
    text, built = render()
    OUT.parent.mkdir(parents=True, exist_ok=True)
    OUT.write_text(text)
    if "--check" in argv:
        errors = check(text, built)
        if errors:
            print("\n".join(errors), file=sys.stderr)
            return 1
        print("performance tables match {} csv sections".format(len(built)))
    return 0


if __name__ == "__main__":
    raise SystemExit(main(sys.argv[1:]))
