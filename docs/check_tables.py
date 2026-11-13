#!/usr/bin/env python3
# Copyright (c) d-SEAMS developers
"""Fail when a hand-copied docs table disagrees with the CSV it names."""

import csv
import re
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
DOCS = ROOT / "docs"
SKIP = {"build", "doctrees", "__pycache__"}
CSV_RE = re.compile(r"(?<![\w./-])((?:[\w.-]+/)*[\w.-]+\.csv)")
NUMERIC = {"median_ms", "iqr_ms", "q1_ms", "q3_ms", "n", "warmup"}
ALIASES = {
    "median ms": "median_ms",
    "median_ms": "median_ms",
    "iqr ms": "iqr_ms",
    "iqr_ms": "iqr_ms",
    "q1 ms": "q1_ms",
    "q1_ms": "q1_ms",
    "q3 ms": "q3_ms",
    "q3_ms": "q3_ms",
    "code": "code",
    "version": "version",
    "molecules": "molecules",
    "operation": "operation",
    "threads": "threads",
    "cores": "cores",
    "n": "n",
    "warmup": "warmup",
    "kernel": "kernel",
    "crystal": "crystal",
    "cell": "cell",
    "phase": "phase",
    "tile": "tile",
}


def pages():
    found = []
    for path in DOCS.rglob("*"):
        if not path.is_file() or path.suffix not in {".org", ".md", ".rst"}:
            continue
        if SKIP.intersection(path.relative_to(DOCS).parts):
            continue
        found.append(path)
    return found


def clean_cell(text):
    text = text.strip()
    text = text.replace("~", "").replace("=", "")
    text = re.sub(r"\\[()]", "", text)
    return " ".join(text.split())


def header_name(cell):
    key = clean_cell(cell).lower()
    return ALIASES.get(key)


def parse_org_tables(lines):
    """Org and markdown pipe tables. Returns (start_line, headers, rows)."""
    tables = []
    index = 0
    while index < len(lines):
        line = lines[index]
        if not line.startswith("|") or index + 1 >= len(lines):
            index += 1
            continue
        sep = lines[index + 1]
        if not sep.startswith("|") or not re.fullmatch(r"\|[-+| ]+\|?", sep.strip()):
            index += 1
            continue
        headers = [header_name(cell) for cell in line.strip().strip("|").split("|")]
        rows = []
        cursor = index + 2
        while cursor < len(lines) and lines[cursor].startswith("|") and not re.fullmatch(
            r"\|[-+| ]+\|?", lines[cursor].strip()
        ):
            cells = [clean_cell(cell) for cell in lines[cursor].strip().strip("|").split("|")]
            rows.append((cursor + 1, cells))
            cursor += 1
        tables.append((index + 1, headers, rows, index))
        index = cursor
    return tables


def parse_csv_tables(lines):
    """Sphinx csv-table directives."""
    tables = []
    index = 0
    while index < len(lines):
        if lines[index].strip() != ".. csv-table::":
            index += 1
            continue
        headers = None
        cursor = index + 1
        while cursor < len(lines) and (lines[cursor].startswith("   ") or lines[cursor].strip() == ""):
            stripped = lines[cursor].strip()
            if stripped.startswith(":header:"):
                headers = [header_name(cell) for cell in stripped.split(":", 1)[1].split(",")]
            elif stripped and not stripped.startswith(":") and headers is not None:
                break
            cursor += 1
        rows = []
        while cursor < len(lines) and lines[cursor].startswith("   ") and lines[cursor].strip():
            cells = [clean_cell(cell) for cell in lines[cursor].strip().split(",")]
            rows.append((cursor + 1, cells))
            cursor += 1
        if headers:
            tables.append((index + 1, headers, rows, index))
        index = max(cursor, index + 1)
    return tables


def parse_grid_tables(lines):
    """Sphinx grid tables exported from Org."""
    tables = []
    index = 0
    while index < len(lines):
        if lines[index].strip() != ".. table::":
            index += 1
            continue
        cursor = index + 1
        while cursor < len(lines) and lines[cursor].strip() == "":
            cursor += 1
        if cursor >= len(lines) or not lines[cursor].strip().startswith("+"):
            index += 1
            continue
        border = lines[cursor]
        cursor += 1
        if cursor >= len(lines) or not lines[cursor].strip().startswith("|"):
            index += 1
            continue
        headers = [header_name(cell) for cell in lines[cursor].strip().strip("|").split("|")]
        cursor += 1
        rows = []
        while cursor < len(lines) and (
            lines[cursor].strip().startswith("+") or lines[cursor].strip().startswith("|")
        ):
            if lines[cursor].strip().startswith("|"):
                cells = [clean_cell(cell) for cell in lines[cursor].strip().strip("|").split("|")]
                rows.append((cursor + 1, cells))
            cursor += 1
        tables.append((index + 1, headers, rows, index))
        index = cursor
        del border
    return tables


def named_csv(lines, start):
    """The CSV named in the lines immediately above a table."""
    chosen = None
    for line in reversed(lines[max(0, start - 12) : start]):
        if line.startswith("*") or line.startswith("#"):
            break
        if line.startswith("|"):
            break
        for match in CSV_RE.finditer(line):
            path = ROOT / match.group(1)
            if path.is_file():
                chosen = path
                return chosen
    return chosen


def load_csv(path, cache):
    if path not in cache:
        with path.open(newline="") as handle:
            cache[path] = list(csv.DictReader(handle))
    return cache[path]


def as_number(text):
    try:
        return float(text)
    except ValueError:
        return None


def close_enough(cell, raw):
    if re.fullmatch(r"-?\d+", cell):
        return int(cell) == int(float(raw))
    got = as_number(cell)
    if got is None:
        return cell == raw
    decimals = len(cell.split(".", 1)[1]) if "." in cell else 0
    return abs(got - float(raw)) <= 1.5 * 10 ** (-decimals)


def check_table(path, start, headers, rows, csv_rows, errors):
    columns = [name for name in headers if name]
    if not any(name in NUMERIC for name in columns):
        return
    keys = [name for name in columns if name not in NUMERIC]
    if not keys:
        return
    index = {}
    for record in csv_rows:
        index.setdefault(tuple(record.get(key, "") for key in keys), []).append(record)
    for line_no, cells in rows:
        if len(cells) < len(headers):
            cells = cells + [""] * (len(headers) - len(cells))
        record_cells = {}
        for name, cell in zip(headers, cells):
            if name:
                record_cells[name] = cell
        key = tuple(record_cells.get(name, "") for name in keys)
        found = index.get(key)
        if not found:
            errors.append("{}:{} no CSV row for {}".format(path, line_no, key))
            continue
        for name in columns:
            if name not in NUMERIC and name not in keys:
                continue
            if name not in found[0]:
                continue
            cell = record_cells.get(name, "")
            if cell == "":
                continue
            if not any(close_enough(cell, record[name]) for record in found):
                errors.append(
                    "{}:{} {} is {}, CSV has {}".format(
                        path, line_no, name, cell, found[0][name]
                    )
                )


def main():
    errors = []
    cache = {}
    checked = 0
    for page in pages():
        lines = page.read_text().splitlines()
        file_csv = None
        for line in lines:
            if "from" in line or "table is" in line or "rows are" in line or "numbers below" in line:
                for match in CSV_RE.finditer(line):
                    path = ROOT / match.group(1)
                    if path.is_file():
                        file_csv = path
        org_tables = [(item, False) for item in parse_org_tables(lines)]
        grid_tables = [(item, False) for item in parse_grid_tables(lines)]
        csv_tables = [(item, True) for item in parse_csv_tables(lines)]
        for (start, headers, rows, origin), generated in org_tables + grid_tables + csv_tables:
            csv_path = named_csv(lines, origin)
            if csv_path is None and generated:
                csv_path = file_csv
            if csv_path is None:
                continue
            if not any(name in NUMERIC for name in headers if name):
                continue
            check_table(
                page.relative_to(ROOT),
                start,
                headers,
                rows,
                load_csv(csv_path, cache),
                errors,
            )
            checked += 1
    if errors:
        print("\n".join(errors), file=sys.stderr)
        print("{} tables checked, {} mismatches".format(checked, len(errors)), file=sys.stderr)
        return 1
    print("docs tables match their CSVs ({} tables)".format(checked))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
