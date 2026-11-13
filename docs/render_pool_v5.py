#!/usr/bin/env python3
# Copyright (c) d-SEAMS developers
"""Render benchmarks/pool-audit-v5.txt into the benchmarks page.

The median of an even sample is the mean of the two central values.
--check recomputes that median from the sample lines and refuses a
stored summary that does not match.
"""

import pathlib
import re
import sys

ROOT = pathlib.Path(__file__).resolve().parents[1]
SRC = ROOT / "benchmarks" / "pool-audit-v5.txt"
CSV = ROOT / "benchmarks" / "pool-audit-v5.csv"
OUT = ROOT / "docs" / "source" / "howto" / "pool-audit-v5-tables.rst"
sys.path.insert(0, str(ROOT / "docs"))
import render_sitting  # noqa: E402

NEIGH = re.compile(
    r"neighbours threads (\d+) seams_ms ([0-9.]+) ovito_ms ([0-9.]+)"
)
RING = re.compile(
    r"rings threads (\d+) base_total ([0-9.]+) new_total ([0-9.]+) "
    r"base_link ([0-9.]+) new_link ([0-9.]+) base_search ([0-9.]+) new_search ([0-9.]+)"
)


def _median(samples):
    _count, median, _iqr = render_sitting.summarize(samples)
    return median


def parse(text):
    lines = text.splitlines()
    rows = []
    index = 0
    while index < len(lines):
        line = lines[index]
        neigh = NEIGH.match(line)
        ring = RING.match(line)
        if neigh:
            threads, seams_ms, ovito_ms = neigh.groups()
            seams_samples = [float(v) for v in lines[index + 1].split()[1:]]
            ovito_samples = [float(v) for v in lines[index + 2].split()[1:]]
            rows.append(
                {
                    "code": "d-SEAMS",
                    "operation": "neighbours",
                    "threads": threads,
                    "n": "8",
                    "median_ms": seams_ms,
                    "samples": seams_samples,
                }
            )
            rows.append(
                {
                    "code": "OVITO",
                    "operation": "neighbours",
                    "threads": threads,
                    "n": "8",
                    "median_ms": ovito_ms,
                    "samples": ovito_samples,
                }
            )
            index += 3
            continue
        if ring:
            threads, base_total, new_total, base_link, new_link, base_search, new_search = ring.groups()
            base_samples = [float(v) for v in lines[index + 1].split()[1:]]
            new_samples = [float(v) for v in lines[index + 2].split()[1:]]
            phases = (
                ("cad887b2", "total", base_total, base_samples),
                ("9fcaab36", "total", new_total, new_samples),
                ("cad887b2", "link", base_link, None),
                ("9fcaab36", "link", new_link, None),
                ("cad887b2", "search", base_search, None),
                ("9fcaab36", "search", new_search, None),
            )
            for code, operation, median, samples in phases:
                rows.append(
                    {
                        "code": code,
                        "operation": operation,
                        "threads": threads,
                        "n": "8",
                        "median_ms": median,
                        "samples": samples,
                    }
                )
            index += 3
            continue
        index += 1
    return rows


def check(rows):
    errors = []
    if len(rows) != 28:
        errors.append("pool-audit v5 has {} rows, expected 28".format(len(rows)))
    for row in rows:
        if not row["samples"]:
            continue
        got = _median(row["samples"])
        if abs(got - float(row["median_ms"])) > 1.5e-3:
            errors.append(
                "{} {} threads {} median {} != {:.3f}".format(
                    row["code"], row["operation"], row["threads"], row["median_ms"], got
                )
            )
    return errors


def write(rows, path):
    lines = [
        ".. list-table::",
        "   :header-rows: 0",
        "",
        "   * - Written by ``docs/render_pool_v5.py`` from ``benchmarks/pool-audit-v5.csv``.",
        "       The CSV is the columnar copy of ``benchmarks/pool-audit-v5.txt``.",
        "       Ice Ic, 884736 molecules, cutoff 3.5 angstrom. ``n`` is 8.",
        "       Neighbour medians are recomputed from the sample lines.",
        "       Ring ``total`` medians are recomputed the same way.",
        "       ``cad887b2`` is the previous full ring call. ``9fcaab36`` is this sitting.",
        "",
        "Neighbours against OVITO",
        "^^^^^^^^^^^^^^^^^^^^^^^^",
        "",
        ".. csv-table::",
        "   :header: code, operation, threads, n, median_ms",
        "",
    ]
    header = ("code", "operation", "threads", "n", "median_ms")
    for row in rows:
        if row["operation"] != "neighbours":
            continue
        lines.append("   " + ", ".join(row[key] for key in header))
    lines.extend(
        [
            "",
            "Full cutoff ring call",
            "^^^^^^^^^^^^^^^^^^^^^^",
            "",
            ".. csv-table::",
            "   :header: code, operation, threads, n, median_ms",
            "",
        ]
    )
    order = {"total": 0, "link": 1, "search": 2}
    rings = [row for row in rows if row["operation"] != "neighbours"]
    rings.sort(key=lambda row: (order[row["operation"]], row["code"], int(row["threads"])))
    for row in rings:
        lines.append("   " + ", ".join(row[key] for key in header))
    lines.append("")
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text("\n".join(lines) + "\n")


def write_csv(rows, path):
    header = "code,operation,threads,n,median_ms"
    body = [header]
    for row in rows:
        body.append(
            "{},{},{},{},{}".format(
                row["code"], row["operation"], row["threads"], row["n"], row["median_ms"]
            )
        )
    path.write_text("\n".join(body) + "\n")


def main(argv):
    rows = parse(SRC.read_text())
    errors = check(rows)
    if errors:
        print("\n".join(errors), file=sys.stderr)
        return 1
    if "--check" in argv:
        if CSV.is_file():
            committed = CSV.read_text().splitlines()[1:]
            fresh = [
                "{},{},{},{},{}".format(
                    row["code"], row["operation"], row["threads"], row["n"], row["median_ms"]
                )
                for row in rows
            ]
            if committed != fresh:
                print("pool-audit-v5.csv does not match the log", file=sys.stderr)
                return 1
        print("pool-audit v5 ok", len(rows), "rows")
        return 0
    write(rows, OUT)
    write_csv(rows, CSV)
    print("wrote", OUT.relative_to(ROOT))
    return 0


if __name__ == "__main__":
    raise SystemExit(main(sys.argv[1:]))
