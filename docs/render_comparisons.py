#!/usr/bin/env python3
# Copyright (c) d-SEAMS developers
"""Render the Comparisons tables from the vendored benchmark CSVs.

The CSVs are copies of the files named in docs/data/comparisons/SOURCE.txt.
--check writes the tables and refuses a cell that is not in its CSV.
"""

import csv
import io
import re
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
DATA = ROOT / "docs" / "data" / "comparisons"
OUT = ROOT / "docs" / "source" / "explanation"
SOURCE = (DATA / "SOURCE.txt").read_text(encoding="utf-8").strip()

CENSUS = re.compile(
    r"\b(51264|51262|512) both (\d+) d-SEAMS-only (\d+) TRACE-only (\d+)"
)
WRAP = re.compile(r"512-0 has .*")
TIP4P = "tip4p2005-4096"
EXACT = ("O-O edges", "4-rings", "5-rings", "6-rings")
LIKE = "ice label, SQQ rule on d-SEAMS chordless 6-rings"
WHOLE = "ice"


def read_csv(name):
    with (DATA / name).open(newline="", encoding="utf-8") as handle:
        rows = list(csv.DictReader(handle))
    if not rows:
        raise SystemExit("{} is empty".format(name))
    return rows


def emit_table(headers, records):
    buffer = io.StringIO()
    writer = csv.writer(buffer, lineterminator="\n")
    lines = [
        ".. csv-table::",
        "   :header: {}".format(",".join(headers)),
        "",
    ]
    for record in records:
        writer.writerow(["" if value is None else str(value) for value in record])
        lines.append("   {}".format(buffer.getvalue().strip()))
        buffer.seek(0)
        buffer.truncate(0)
    lines.append("")
    return "\n".join(lines)


def require(text, needle, where):
    if needle not in text:
        raise SystemExit("{} is missing from {}".format(needle, where))


def census(text):
    found = []
    for match in CENSUS.finditer(text):
        found.append(
            (
                match.group(1),
                int(match.group(2)),
                int(match.group(3)),
                int(match.group(4)),
            )
        )
    if not found:
        raise SystemExit("cage census missing: {}".format(text[:80]))
    return found


def summary_tables(rows):
    headers = [
        "comparator",
        "version",
        "operation",
        "molecules",
        "threads",
        "cores",
        "n",
        "d-SEAMS ms",
        "other ms",
        "ratio",
        "verdict",
        "agreement",
    ]
    records = []
    for row in rows:
        for key in ("version", "threads", "cores", "n"):
            if row.get(key, "") == "" and row.get("verdict") != "absent":
                raise SystemExit("summary row lacks {}: {}".format(key, row.get("competitor")))
        records.append(
            [
                row["competitor"],
                row["version"],
                row["operation"],
                row["molecules"],
                row["threads"],
                row["cores"],
                row["n"],
                row["dseams_median_ms"],
                row["other_median_ms"],
                row["ratio"],
                row["verdict"],
                row["agreement"],
            ]
        )
    body = [
        "Head-to-head summary",
        "--------------------",
        "",
        SOURCE,
        "",
        "Each row is one record of ``h2h-summary.csv``. "
        "The version, thread count, repeat count ``n`` and core count are the CSV columns.",
        "",
        emit_table(headers, records),
    ]
    return "\n".join(body)


def clavis_tables(rows):
    by_key = {}
    for row in rows:
        by_key[(row["molecules"], row["threads"])] = row
        for key in ("version", "threads", "cores", "n"):
            if not row.get(key):
                raise SystemExit("Clavis row lacks {}".format(key))
    time_rows = []
    census_rows = []
    notes = []
    for row in rows:
        time_rows.append(
            [
                row["version"],
                row["molecules"],
                row["threads"],
                row["cores"],
                row["n"],
                row["dseams_median_ms"],
                row["other_median_ms"],
                row["ratio"],
            ]
        )
        for signature, both, only_us, only_trace in census(row["agreement"]):
            census_rows.append(
                [
                    row["version"],
                    row["molecules"],
                    row["threads"],
                    row["cores"],
                    row["n"],
                    signature,
                    both,
                    only_us,
                    only_trace,
                ]
            )
    small = by_key[("46", "1")]
    parsed = census(small["agreement"])
    dropped = sum(item[2] for item in parsed)
    trace_only = sum(item[3] for item in parsed)
    if dropped != 7 or trace_only != 0:
        raise SystemExit("sI wrap census is {} drops, {} TRACE-only".format(dropped, trace_only))
    if ("46", "8") not in by_key or census(by_key[("46", "8")]["agreement"]) != parsed:
        raise SystemExit("the 8-thread sI census does not match the 1-thread census")
    wrap = WRAP.search(small["agreement"])
    if wrap is None:
        raise SystemExit("the periodic-wrap clause is not in the sI agreement")
    notes.append(
        "On the 46-molecule sI cell the census drops {} cages that only d-SEAMS counts, and TRACE counts none that d-SEAMS does not.".format(
            dropped
        )
    )
    notes.append(wrap.group(0))
    large = by_key.get(("110160", "1"))
    if large is None or "serial" not in large["agreement"]:
        raise SystemExit("the serial TRACE note is not in the 110160-molecule row")
    one = by_key[("110160", "1")]
    eight = by_key[("110160", "8")]
    notes.append(
        "On {} molecules, one thread, d-SEAMS takes {} ms and TRACE takes {} ms. "
        "At eight threads d-SEAMS takes {} ms and TRACE takes {} ms. {}".format(
            one["molecules"],
            one["dseams_median_ms"],
            one["other_median_ms"],
            eight["dseams_median_ms"],
            eight["other_median_ms"],
            large["agreement"].split(". ")[-1],
        )
    )
    body = [
        "ClavisCage TRACE",
        "----------------",
        "",
        SOURCE,
        "",
        "Times are the median milliseconds in ``h2h-clavis-cages.csv``. "
        "The cage counts are parsed from that file's agreement column.",
        "",
    ]
    body.extend(notes)
    body.append("")
    body.append(
        emit_table(
            [
                "version",
                "molecules",
                "threads",
                "cores",
                "n",
                "d-SEAMS ms",
                "TRACE ms",
                "ratio",
            ],
            time_rows,
        )
    )
    body.append(
        emit_table(
            [
                "version",
                "molecules",
                "threads",
                "cores",
                "n",
                "signature",
                "both",
                "d-SEAMS only",
                "TRACE only",
            ],
            census_rows,
        )
    )
    return "\n".join(body)


def sqq_tables(timing, agree):
    versions = sorted({(row["dseams_version"], row["sqq_version"]) for row in timing})
    if versions != [("2.14.0", "0.5.8")]:
        raise SystemExit("unexpected SQQ versions {}".format(versions))
    time_rows = []
    featured = {}
    for row in timing:
        for key in ("threads", "cores", "n", "dseams_version", "sqq_version"):
            if not row.get(key):
                raise SystemExit("SQQ timing row lacks {}".format(key))
        time_rows.append(
            [
                "{}/{}".format(row["dseams_version"], row["sqq_version"]),
                row["frame"],
                row["molecules"],
                row["operation"],
                row["threads"],
                row["cores"],
                row["n"],
                row["dseams_mean_ms"],
                row["sqq_mean_ms"],
                row["ratio_sqq_over_dseams"],
            ]
        )
        if row["frame"] == TIP4P:
            featured[(row["operation"], row["threads"])] = row
    tip = [row for row in agree if row["frame"] == TIP4P]
    if not tip:
        raise SystemExit("TIP4P/2005 frame is missing from the agreement CSV")
    exact = []
    for quantity in EXACT:
        match = [row for row in tip if row["quantity"] == quantity]
        if len(match) != 1:
            raise SystemExit("expected one {} row".format(quantity))
        row = match[0]
        if row["d_seams"] != row["sqq"] or row["agreement"] != "1.000000":
            raise SystemExit("{} does not match on {}".format(quantity, TIP4P))
        exact.append("{} {} = {}".format(quantity, row["d_seams"], row["sqq"]))
    like = [row for row in tip if row["quantity"] == LIKE]
    if len(like) != 1 or like[0]["agreement"] != "1.000000":
        raise SystemExit("chordless like-for-like labels do not match on {}".format(TIP4P))
    for operation, threads in (
        ("ice-chordless", "1"),
        ("ice-chordless", "8"),
        (WHOLE, "1"),
        (WHOLE, "8"),
    ):
        if (operation, threads) not in featured:
            raise SystemExit("missing {} threads {}".format(operation, threads))
    chord_1 = featured[("ice-chordless", "1")]
    chord_8 = featured[("ice-chordless", "8")]
    whole_1 = featured[(WHOLE, "1")]
    whole_8 = featured[(WHOLE, "8")]
    agree_rows = []
    for row in agree:
        agree_rows.append(
            [
                "2.14.0/0.5.8",
                row["frame"],
                row["molecules"],
                row["quantity"],
                row["d_seams"],
                row["sqq"],
                row["both"],
                row["agreement"],
            ]
        )
    body = [
        "SQQ",
        "---",
        "",
        SOURCE,
        "",
        "On the {}-molecule TIP4P/2005 frame the edges and rings match: {}.".format(
            tip[0]["molecules"], "; ".join(exact)
        ),
        "The chordless census is the like-for-like ice row. "
        "At one thread the SQQ/d-SEAMS ratio is {} (version {}, threads {}, n={}, cores {}) "
        "and at eight threads it is {} (version {}, threads {}, n={}, cores {}).".format(
            chord_1["ratio_sqq_over_dseams"],
            "2.14.0/0.5.8",
            chord_1["threads"],
            chord_1["n"],
            chord_1["cores"],
            chord_8["ratio_sqq_over_dseams"],
            "2.14.0/0.5.8",
            chord_8["threads"],
            chord_8["n"],
            chord_8["cores"],
        ),
        "The whole-tool ratio is CHILL+ against SQQ's ice stage. "
        "At one thread it is {} (version {}, threads {}, n={}, cores {}) "
        "and at eight threads it is {} (version {}, threads {}, n={}, cores {}).".format(
            whole_1["ratio_sqq_over_dseams"],
            "2.14.0/0.5.8",
            whole_1["threads"],
            whole_1["n"],
            whole_1["cores"],
            whole_8["ratio_sqq_over_dseams"],
            "2.14.0/0.5.8",
            whole_8["threads"],
            whole_8["n"],
            whole_8["cores"],
        ),
        "",
        "Agreement. The version is the build named on every timing row. "
        "This census is not itself a timed run.",
        "",
        emit_table(
            [
                "version",
                "frame",
                "molecules",
                "quantity",
                "d-SEAMS",
                "SQQ",
                "both",
                "agreement",
            ],
            agree_rows,
        ),
        "Timing. Every row names the version, the thread count, the repeat count and the core count.",
        "",
        emit_table(
            [
                "version",
                "frame",
                "molecules",
                "operation",
                "threads",
                "cores",
                "n",
                "d-SEAMS ms",
                "SQQ ms",
                "SQQ/d-SEAMS",
            ],
            time_rows,
        ),
    ]
    return "\n".join(body)


def reject_hosts(text):
    for banned in ("Intel", "Xeon", "cpu_model", "localhost"):
        if banned in text:
            raise SystemExit("rendered comparisons name a machine: {}".format(banned))


def main(argv):
    summary = read_csv("h2h-summary.csv")
    clavis = read_csv("h2h-clavis-cages.csv")
    timing = read_csv("h2h-sqq.csv")
    agree = read_csv("h2h-sqq-agree.csv")
    pages = {
        "comparisons-summary-tables.rst": summary_tables(summary),
        "clavis-trace-tables.rst": clavis_tables(clavis),
        "sqq-tables.rst": sqq_tables(timing, agree),
    }
    for name, text in pages.items():
        reject_hosts(text)
        for needle in ("323bfd3872e726b3fa384a67c6b1bedf627ce6a0", "c7f71310f27c5ce62d11e97fd82f9757fa4e8628"):
            require(text, needle, name)
        (OUT / name).write_text(text, encoding="utf-8")
    if "--check" in argv:
        written = "\n".join(pages.values())
        for row in summary:
            require(written, row["version"], "summary")
            if row["dseams_median_ms"]:
                require(written, row["dseams_median_ms"], "summary")
        for row in clavis:
            require(written, row["dseams_median_ms"], "clavis")
            require(written, row["other_median_ms"], "clavis")
        for row in timing:
            if row["frame"] == TIP4P:
                require(written, row["ratio_sqq_over_dseams"], "sqq")
        for row in agree:
            if row["frame"] == TIP4P and row["quantity"] in EXACT:
                require(written, row["d_seams"], "sqq agree")
    return 0


if __name__ == "__main__":
    sys.exit(main(sys.argv[1:]))
