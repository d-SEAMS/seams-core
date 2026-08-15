#!/usr/bin/env python3
"""Assemble the Elja paper-bench outputs into one JSON object.

Usage: elja_bench_json.py <out_dir>

Reads the text outputs written by elja_paper_benches.sh and emits the JSON
shape the campaign spec asks for on stdout. Missing files become nulls rather
than errors, so a partial run still yields an inspectable object.
"""
import json
import pathlib
import re
import sys


def read(path):
    try:
        return path.read_text()
    except OSError:
        return ""


def parse_scaling(text):
    rows = {}
    for line in text.splitlines():
        m = re.match(
            r"\s*(\d+)\s+([\d.]+)\s+([\d.]+)\s+([\d.]+)\s+([\d.]+|skipped)",
            line,
        )
        if m:
            n = int(m.group(1))
            rows[n] = {
                "neighListO_ms": float(m.group(2)),
                "byIndex_ms": float(m.group(3)),
                "getCorrelPlus_ms": float(m.group(4)),
                "ringNetwork_ms": (
                    None if m.group(5) == "skipped" else float(m.group(5))
                ),
            }
    return rows


def parse_kv(text):
    out = {}
    for line in text.splitlines():
        m = re.match(r"([A-Za-z/_-]+)\s+([-\d.]+)\s*$", line.strip())
        if m:
            out[m.group(1)] = float(m.group(2))
    return out


def parse_overhead(text):
    rows = {}
    for line in text.splitlines():
        m = re.match(r"\s*(\d+)\s+([\d.]+)\s+([\d.]+)\s+([\d.]+)", line)
        if m:
            rows[int(m.group(1))] = {
                "vesin_ms": float(m.group(2)),
                "neighListO_ms": float(m.group(3)),
            }
    return rows


def parse_strong(text):
    for line in text.splitlines():
        m = re.match(r"\s*(\d+)\s+(\d+)\s+([\d.]+)", line)
        if m:
            return float(m.group(3))
    return None


def main():
    out_dir = pathlib.Path(sys.argv[1])
    conditions = {}
    for line in read(out_dir / "conditions.txt").splitlines():
        if ":" in line:
            key, _, value = line.partition(":")
            conditions[key.strip()] = value.strip()

    result = {
        "conditions": conditions,
        "tip": {
            "scaling": parse_scaling(read(out_dir / "tip-scaling.txt")),
            "cages": parse_kv(read(out_dir / "tip-cages.txt")),
            "overhead": parse_overhead(read(out_dir / "tip-overhead.txt")),
            "strong": {
                t: parse_strong(read(out_dir / f"tip-strong-t{t}.txt"))
                for t in (1, 2, 4)
            },
        },
        "base": {
            "scaling": parse_scaling(read(out_dir / "base-scaling.txt")),
            "cages": parse_kv(read(out_dir / "base-cages.txt")),
        },
    }
    json.dump(result, sys.stdout, indent=2, sort_keys=True)
    print()


if __name__ == "__main__":
    main()
