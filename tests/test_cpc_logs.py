# Copyright (c) d-SEAMS developers
"""The CPC records are text logs, not built programs."""

import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
RECORDS = ROOT / "benchmarks" / "cpc" / "records"


def main() -> int:
    if not RECORDS.is_dir():
        print("missing", RECORDS, file=sys.stderr)
        return 1
    files = [p for p in RECORDS.rglob("*") if p.is_file()]
    if len(files) < 10:
        print("too few CPC records", len(files), file=sys.stderr)
        return 1
    logs = [p for p in files if p.suffix == ".log"]
    if not logs:
        print("no raw .log files", file=sys.stderr)
        return 1
    claims = RECORDS / "CLAIMS.tsv"
    text = claims.read_text(encoding="utf-8")
    rows = [line for line in text.splitlines() if line.strip()]
    if len(rows) < 2 or "\t" not in rows[0]:
        print("CLAIMS.tsv is not a table", file=sys.stderr)
        return 1
    for path in files:
        blob = path.read_bytes()
        if b"\0" in blob:
            print("binary record", path, file=sys.stderr)
            return 1
    print(f"cpc records {len(files)} logs {len(logs)} claims {len(rows) - 1}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
