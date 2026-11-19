#!/usr/bin/env python3
# Copyright (c) d-SEAMS developers
"""Checks for committed benchmark records."""

from __future__ import annotations

import subprocess
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
