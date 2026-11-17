#!/usr/bin/env python3
# Copyright (c) d-SEAMS developers
"""Rebuild h2h-cpc/summary.tsv from the raw per-repeat logs.

Quartiles are type 7. For an even count the median is the mean of the
two central samples. Warm-up calls are kept in the logs and left out of
the median.
"""

from __future__ import annotations

import math
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
RAW = ROOT / "h2h-cpc" / "raw"
OUT = ROOT / "h2h-cpc" / "summary.tsv"
MIN_REPEATS = 6

COLUMNS = (
    "comparison",
    "system",
    "molecules",
    "tool",
    "repeats",
    "median_ms",
    "q1_ms",
    "q3_ms",
    "iqr_ms",
    "histogram",
    "csv_sha256",
    "identity",
    "git_sha",
    "cutoff",
    "min_ring",
    "max_ring",
    "threads",
    "source",
    "freud_version",
    "ovito_version",
    "numpy_version",
    "scipy_version",
    "cpu_model",
)


def type7(values: list[float], probability: float) -> float:
    ordered = sorted(values)
    count = len(ordered)
    if count == 0:
        raise ValueError("empty sample")
    if count == 1:
        return ordered[0]
    position = (count - 1) * probability + 1.0
    index = int(math.floor(position))
    if index >= count:
        return ordered[-1]
    if index < 1:
        return ordered[0]
    weight = position - index
    return ordered[index - 1] + weight * (ordered[index] - ordered[index - 1])


def self_check() -> None:
    sample = [1.0, 2.0, 3.0, 4.0, 5.0, 6.0]
    median = type7(sample, 0.5)
    q1 = type7(sample, 0.25)
    q3 = type7(sample, 0.75)
    if abs(median - 3.5) > 1e-12 or abs(q1 - 2.25) > 1e-12 or abs(q3 - 4.75) > 1e-12:
        raise SystemExit(f"quartile self-check failed: {median} {q1} {q3}")


def parse_log(path: Path) -> dict[str, str]:
    record: dict[str, str] = {}
    for line in path.read_text(encoding="utf-8").splitlines():
        if not line.strip():
            continue
        key, _, value = line.partition(" ")
        record[key] = value
    return record


def collect() -> list[dict[str, str]]:
    grouped: dict[tuple[str, str, str], list[dict[str, str]]] = {}
    unsupported: list[dict[str, str]] = []
    for path in sorted(RAW.rglob("*.txt")):
        if path.name in {"environment.txt", "done.txt", "order.txt"}:
            continue
        record = parse_log(path)
        if record.get("status") == "unsupported":
            unsupported.append(record)
            continue
        if "time_ms" not in record or "phase" not in record:
            continue
        key = (record["comparison"], record["system"], record["tool"])
        grouped.setdefault(key, []).append(record)

    rows: list[dict[str, str]] = []
    identities: dict[tuple[str, str], dict[str, str]] = {}
    for key, records in sorted(grouped.items()):
        timed = [r for r in records if r.get("phase") == "timed"]
        if len(timed) < MIN_REPEATS:
            raise SystemExit(f"{key} has {len(timed)} timed repeats, need {MIN_REPEATS}")
        samples = [float(r["time_ms"]) for r in timed]
        shas = {r.get("csv_sha256", "") for r in records if r.get("csv_sha256")}
        if len(shas) != 1:
            raise SystemExit(f"{key} identity changed across repeats: {sorted(shas)}")
        sha = next(iter(shas))
        histograms = {r.get("histogram", r.get("tally", "")) for r in records}
        if len(histograms) != 1:
            raise SystemExit(f"{key} histogram changed across repeats: {sorted(histograms)}")
        comparison, system, tool = key
        identities.setdefault((comparison, system), {})[tool] = sha
        head = timed[0]
        median = type7(samples, 0.5)
        q1 = type7(samples, 0.25)
        q3 = type7(samples, 0.75)
        rows.append({
            "comparison": comparison,
            "system": system,
            "molecules": head["molecules"],
            "tool": tool,
            "repeats": str(len(samples)),
            "median_ms": f"{median:.6f}",
            "q1_ms": f"{q1:.6f}",
            "q3_ms": f"{q3:.6f}",
            "iqr_ms": f"{(q3 - q1):.6f}",
            "histogram": head.get("histogram", head.get("tally", "")),
            "csv_sha256": sha,
            "identity": "",
            "git_sha": head.get("git_sha", ""),
            "cutoff": head.get("cutoff", ""),
            "min_ring": head.get("min_ring", ""),
            "max_ring": head.get("max_ring", ""),
            "threads": head.get("threads", ""),
            "source": head.get("source", ""),
            "freud_version": head.get("freud_version", ""),
            "ovito_version": head.get("ovito_version", ""),
            "numpy_version": head.get("numpy_version", ""),
            "scipy_version": head.get("scipy_version", ""),
            "cpu_model": head.get("cpu_model", ""),
        })

    for record in unsupported:
        comparison = record["comparison"]
        system = record["system"]
        tool = record["tool"]
        identities.setdefault((comparison, system), {})[tool] = "unsupported"
        rows.append({
            "comparison": comparison,
            "system": system,
            "molecules": record.get("molecules", ""),
            "tool": tool,
            "repeats": "0",
            "median_ms": "",
            "q1_ms": "",
            "q3_ms": "",
            "iqr_ms": "",
            "histogram": "",
            "csv_sha256": "",
            "identity": "unsupported",
            "git_sha": record.get("git_sha", ""),
            "cutoff": record.get("cutoff", ""),
            "min_ring": record.get("min_ring", ""),
            "max_ring": record.get("max_ring", ""),
            "threads": record.get("threads", ""),
            "source": record.get("source", ""),
            "freud_version": record.get("freud_version", ""),
            "ovito_version": record.get("ovito_version", ""),
            "numpy_version": record.get("numpy_version", ""),
            "scipy_version": record.get("scipy_version", ""),
            "cpu_model": record.get("cpu_model", ""),
        })

    for row in rows:
        if row["identity"] == "unsupported":
            continue
        shas = identities[(row["comparison"], row["system"])]
        measured = {tool: sha for tool, sha in shas.items() if sha != "unsupported"}
        if measured and len(set(measured.values())) == 1:
            row["identity"] = "match"
        else:
            row["identity"] = "differ"
    return rows


def main() -> None:
    self_check()
    if not RAW.is_dir():
        raise SystemExit(f"no raw logs under {RAW.relative_to(ROOT)}")
    rows = collect()
    lines = ["\t".join(COLUMNS)]
    for row in rows:
        lines.append("\t".join(row[column] for column in COLUMNS))
    OUT.parent.mkdir(parents=True, exist_ok=True)
    OUT.write_text("\n".join(lines) + "\n", encoding="utf-8")
    print(f"wrote {OUT.relative_to(ROOT)} rows {len(rows)}")


if __name__ == "__main__":
    main()
