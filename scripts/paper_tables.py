#!/usr/bin/env python3
"""Write one CSV per paper table from benchmarks/paper-cells/sitting.csv.

The sitting is the only timing source. These files are selections of
its rows. A table that is not a 110592-molecule or hydrate-cage timing
is left to its own rule.
"""

from __future__ import annotations

import csv
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]


def load(path: Path) -> list[dict]:
    with path.open(newline="") as handle:
        return list(csv.DictReader(handle))


def write(path: Path, rows: list[dict]) -> None:
    if not rows:
        raise SystemExit(f"no rows for {path.name}")
    with path.open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(rows[0]))
        writer.writeheader()
        writer.writerows(rows)
    print(f"wrote {path} ({len(rows)} rows)", flush=True)


def select(rows: list[dict], **want) -> list[dict]:
    found = []
    for row in rows:
        if all(row.get(key) == value for key, value in want.items()):
            found.append(row)
    return found


def main() -> int:
    out = Path(sys.argv[1] if len(sys.argv) > 1 else ROOT / "benchmarks/paper-cells")
    rows = load(out / "sitting.csv")
    ice = [row for row in rows if row["phase"] == "ice"]
    write(out / "tab-rings-214.csv", [
        row for row in ice
        if row["operation"] == "rings"
        and row["code"] in ("d-SEAMS", "rings-statistics-matter")
    ])
    write(out / "tab-ql-eq.csv", [
        row for row in ice
        if row["operation"] in ("q4", "q6", "q4-bar", "q6-bar", "q4+q6", "q4-core", "q6-core")
        and row["code"] == "d-SEAMS"
    ])
    write(out / "tab-neighbours.csv", [
        row for row in ice if row["operation"] == "neighbours"
    ])
    write(out / "tab-chill.csv", [
        row for row in ice if row["operation"] == "chill+"
    ])
    write(out / "frag-cli-q.csv", [
        row for row in ice if row["operation"] == "q4+q6" and row["code"] == "d-SEAMS"
    ])
    write(out / "frag-cli-rsm.csv", [
        row for row in ice
        if row["operation"] == "rings" and row["code"] == "rings-statistics-matter"
    ])
    write(out / "frag-lit-q6.csv", [
        row for row in ice if row["operation"] in ("q4", "q6", "q4-core", "q6-core")
    ])
    cages = [
        row for row in rows
        if row["operation"] == "cages" and row["phase"] != "ice"
    ]
    write(out / "frag-lit-cages.csv", cages)
    write(out / "tab-cages.csv", cages)
    sym = [
        row for row in ice
        if row["operation"] in ("classifier", "ptm-default", "ptm-ice")
    ]
    write(out / "symmetry-vs-ptm.csv", symmetry_rows(sym))
    write(out / "tab-sym-eq.csv", sym_eq_rows(sym))
    ice_ptm = [
        row for row in sym
        if row["operation"] == "ptm-ice" and row["threads"] == "8"
    ]
    write(out / "frag-cli-ptm.csv", ice_ptm)
    return 0


def symmetry_rows(rows: list[dict]) -> list[dict]:
    classifier = {
        (row["threads"]): float(row["median_ms"])
        for row in rows if row["operation"] == "classifier"
    }
    names = {
        "classifier": ("classifier", "classifier"),
        "ptm-default": ("PTM-default", "default"),
        "ptm-ice": ("PTM-ice", "ice"),
    }
    out = []
    for row in rows:
        code, kind = names[row["operation"]]
        median = float(row["median_ms"])
        base = classifier[row["threads"]]
        sink = crystal = ""
        if row["operation"] == "classifier":
            detail = row["detail"]
            sink = detail.split("sink=")[1].split()[0]
            crystal = detail.split("crystal=")[1].split()[0]
        out.append({
            "source_commit": row["source_commit"],
            "cpu_model": row["cpu_model"],
            "cpu_family": row["cpu_family"],
            "cpu_model_id": row["cpu_model_id"],
            "cpu_stepping": row["cpu_stepping"],
            "cpu_flags": row["cpu_flags"],
            "cpu_cache": row["cpu_cache"],
            "cpu_cores": row["cpu_cores"],
            "code": code,
            "template_set": kind,
            "molecules": row["molecules"],
            "threads": row["threads"],
            "affinity": row["affinity"],
            "warmup": row["warmup"],
            "n": row["n"],
            "median_ms": row["median_ms"],
            "iqr_ms": row["iqr_ms"],
            "q1_ms": row["q1_ms"],
            "q3_ms": row["q3_ms"],
            "sink": sink,
            "crystal": crystal,
            "ratio_this_over_classifier": f"{median / base:.6f}",
            "ratio_classifier_over_this": f"{base / median:.6f}",
            "samples_ms": row["samples_ms"],
            "detail": row["detail"],
        })
    return out


def sym_eq_rows(rows: list[dict]) -> list[dict]:
    cells = {
        "classifier": "classifier against PTM",
        "ptm-default": "PTM default templates",
        "ptm-ice": "PTM diamond templates",
    }
    codes = {
        "classifier": "symmetry",
        "ptm-default": "ovito-symmetry-default",
        "ptm-ice": "ovito-symmetry",
    }
    out = []
    for row in rows:
        molecules = int(row["molecules"])
        median = float(row["median_ms"])
        sink = ""
        if row["operation"] == "classifier":
            sink = row["detail"].split("sink=")[1].split()[0]
        out.append({
            "source_commit": row["source_commit"],
            "cpu_model": row["cpu_model"],
            "cpu_family": row["cpu_family"],
            "cpu_model_id": row["cpu_model_id"],
            "cpu_stepping": row["cpu_stepping"],
            "cpu_flags": row["cpu_flags"],
            "cpu_cache": row["cpu_cache"],
            "cpu_cores": row["cpu_cores"],
            "cell": cells[row["operation"]],
            "code": codes[row["operation"]],
            "molecules": row["molecules"],
            "threads": row["threads"],
            "affinity": row["affinity"],
            "warmup": row["warmup"],
            "n": row["n"],
            "median_ms": row["median_ms"],
            "iqr_ms": row["iqr_ms"],
            "q1_ms": row["q1_ms"],
            "q3_ms": row["q3_ms"],
            "us_per_atom": f"{median * 1000.0 / molecules:.6f}",
            "sink": sink,
            "samples_ms": row["samples_ms"],
            "detail": row["detail"],
        })
    return out


if __name__ == "__main__":
    sys.exit(main())
