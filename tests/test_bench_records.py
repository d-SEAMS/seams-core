#!/usr/bin/env python3
# Copyright (c) d-SEAMS developers
"""Checks for committed benchmark records."""

from __future__ import annotations

import subprocess
import csv
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]


def _tracked_bench_files() -> list[Path]:
    done = subprocess.run(
        ["git", "-C", str(ROOT), "ls-files", "benchmarks"],
        check=True,
        text=True,
        capture_output=True,
    )
    rows = []
    for line in done.stdout.splitlines():
        if not line:
            continue
        path = ROOT / line
        if path.suffix in {".csv", ".txt"}:
            rows.append(path)
    return rows


def test_committed_bench_records_are_lf_only():
    for path in _tracked_bench_files():
        body = path.read_bytes()
        assert b"\r\n" not in body, f"{path} has CRLF line endings"
        assert b"\r" not in body, f"{path} has CR bytes"


def _object_exists(spec: str) -> bool:
    done = subprocess.run(
        ["git", "-C", str(ROOT), "cat-file", "-e", spec],
        check=False,
        text=True,
        capture_output=True,
    )
    return done.returncode == 0


def test_sitting_rows_have_resolvable_provenance():
    samples = ROOT / "benchmarks" / "fas230-samples.csv"
    rows = list(csv.DictReader(samples.open(newline="", encoding="utf-8")))
    assert rows, "samples csv is empty"
    for row in rows:
        base_ref = row["base_ref"]
        h2h_blob = row["h2h_fas230_blob"]
        fas_blob = row["h2h_fastatomstruct_blob"]
        assert _object_exists(f"{base_ref}^{{commit}}"), base_ref
        assert _object_exists(f"{h2h_blob}^{{blob}}"), h2h_blob
        assert _object_exists(f"{fas_blob}^{{blob}}"), fas_blob


def test_thread_probe_has_at_least_three_samples_per_arm():
    probe = ROOT / "benchmarks" / "fas230-thread-probe.csv"
    rows = list(csv.DictReader(probe.open(newline="", encoding="utf-8")))
    assert rows, "thread probe csv is empty"
    by_key = {}
    for row in rows:
        key = (row["arm"], row["size"], row["operation"])
        by_key[key] = by_key.get(key, 0) + 1
        assert int(row["os_threads"]) >= 1, row
    for key, count in by_key.items():
        assert count >= 3, (key, count)
