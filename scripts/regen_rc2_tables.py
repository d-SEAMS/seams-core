#!/usr/bin/env python3
"""Fill the neighbour, CHILL+, ring and literature tables from a sitting log.

The timers stay in scripts/bench_neighbours_v2.py and scripts/h2h_literature.py.
This script only reads the text those runs wrote.
"""

from __future__ import annotations

import argparse
import csv
import re
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]

LINE = re.compile(
    r"(?P<tag>mode|kind)=(?P<name>\S+)\s+molecules=(?P<mol>\d+)\s+"
    r"threads=(?P<thr>\d+).*?median_ms=(?P<med>[-+0-9.nan]+)"
)


def load_lines(path: Path):
    rows = {}
    for line in path.read_text(encoding="utf-8").splitlines():
        match = LINE.search(line)
        if not match:
            continue
        key = (match.group("name"), int(match.group("mol")), int(match.group("thr")))
        med = match.group("med")
        rows[key] = None if med == "nan" else float(med)
    return rows


def cell(rows, name, molecules, threads):
    value = rows.get((name, molecules, threads))
    if value is None:
        return ""
    return f"{value:.3f}"


def markdown_tables(rows):
    tilt = [
        "| molecules | threads | d-SEAMS | vesin | freud | orthorhombic |",
        "| ---: | ---: | ---: | ---: | ---: | ---: |",
    ]
    for molecules, threads in (
        (10648, 1), (10648, 8), (110592, 1), (110592, 8), (884736, 1), (884736, 8),
    ):
        tilt.append(
            f"| {molecules} | {threads} | {cell(rows, 'tilt', molecules, threads)} | "
            f"{cell(rows, 'tilt-vesin', molecules, threads)} | "
            f"{cell(rows, 'tilt-freud', molecules, threads)} | "
            f"{cell(rows, 'ortho', molecules, threads)} |"
        )
    chill = [
        "| molecules | threads | d-SEAMS | mdapy | OVITO |",
        "| ---: | ---: | ---: | ---: | ---: |",
    ]
    rings = list(chill)
    diamond = [
        "| molecules | threads | d-SEAMS | mdapy diamond | OVITO identify diamond | mdapy PTM | OVITO PTM |",
        "| ---: | ---: | ---: | ---: | ---: | ---: | ---: |",
    ]
    for molecules, threads in ((110592, 1), (110592, 8), (884736, 1), (884736, 8)):
        chill.append(
            f"| {molecules} | {threads} | {cell(rows, 'chill', molecules, threads)} | "
            f"{cell(rows, 'mdapy-chill', molecules, threads)} | "
            f"{cell(rows, 'ovito-chill', molecules, threads)} |"
        )
        rings.append(
            f"| {molecules} | {threads} | {cell(rows, 'rings', molecules, threads)} | "
            f"{cell(rows, 'mdapy-rings', molecules, threads)} | "
            f"{cell(rows, 'ovito-rings', molecules, threads)} |"
        )
        diamond.append(
            f"| {molecules} | {threads} | {cell(rows, 'symmetry', molecules, threads)} | "
            f"{cell(rows, 'mdapy-diamond', molecules, threads)} | "
            f"{cell(rows, 'ovito-diamond', molecules, threads)} | "
            f"{cell(rows, 'mdapy-symmetry', molecules, threads)} | "
            f"{cell(rows, 'ovito-symmetry', molecules, threads)} |"
        )
    return "\n".join(tilt), "\n".join(chill), "\n".join(diamond), "\n".join(rings)


def replace_table(text, heading, table):
    pattern = re.compile(
        rf"(## {re.escape(heading)}\n\n(?:.*\n)*?)(\| molecules \|.*\n(?:\|.*\n)+)",
        re.M,
    )
    match = pattern.search(text)
    if not match:
        raise SystemExit(f"missing table under {heading}")
    return text[: match.start(2)] + table + "\n" + text[match.end(2) :]


def raw_block(path: Path):
    lines = []
    for line in path.read_text(encoding="utf-8").splitlines():
        if line.startswith("mode=") or line.startswith("kind="):
            lines.append(line)
    return "\n".join(lines)


def rewrite_markdown(md_path: Path, log_path: Path):
    rows = load_lines(log_path)
    tilt, chill, diamond, rings = markdown_tables(rows)
    text = md_path.read_text(encoding="utf-8")
    text = replace_table(text, "Triclinic neighbour list", tilt)
    text = replace_table(text, "CHILL+", chill)
    text = replace_table(text, "Diamond common-neighbour test", diamond)
    text = replace_table(text, "Rings", rings)
    text = re.sub(
        r"## Raw log\n\n```\n.*?\n```",
        "## Raw log\n\n```\n" + raw_block(log_path) + "\n```",
        text,
        count=1,
        flags=re.S,
    )
    text = text.replace(
        "NEIGH_BENCH=build-release/tests/bench_v2",
        "NEIGH_BENCH=bbdir-release/tests/bench_v2",
    )
    md_path.write_text(text, encoding="utf-8")


def literature_focus(csv_path: Path):
    """110592 ice Ic rows, one line per code."""
    out = []
    with csv_path.open(newline="") as handle:
        for row in csv.DictReader(handle):
            if row.get("molecules") != "110592" or row.get("phase") != "Ic":
                continue
            out.append(row)
    return out


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--v3-log", type=Path, required=True)
    parser.add_argument("--markdown", type=Path, default=ROOT / "markdown" / "neighbours-v3.md")
    parser.add_argument("--literature", type=Path)
    args = parser.parse_args()
    rewrite_markdown(args.markdown, args.v3_log)
    if args.literature and args.literature.is_file():
        for row in literature_focus(args.literature):
            print(
                "{threads} {operation} {code} median={median_ms} iqr={iqr_ms} detail={detail}".format(
                    **row
                )
            )


if __name__ == "__main__":
    main()
