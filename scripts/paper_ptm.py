#!/usr/bin/env python3
"""End-to-end diamond classifier against two-template PTM.

One sitting. Ideal ice Ic, 110592 molecules, eight threads, one warm-up
and eight calls. PTM is PolyhedralTemplateMatchingModifier with RMSD
cutoff 0.1 and only the cubic and hexagonal diamond templates. The
classifier is bench_v2 mode symmetry: the cutoff neighbour list plus
diamondFromCsr. NEIGH_ONLY on scripts/bench_neighbours_v2.py selects
those two jobs and no others.
"""

from __future__ import annotations

import argparse
import csv
import os
import subprocess
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "scripts"))
import paper_stats as stats

MOLECULES = 110592


def parse_line(line: str) -> dict:
    fields = {}
    for token in line.split():
        if "=" not in token:
            continue
        key, value = token.split("=", 1)
        fields[key] = value
    fields["kind"] = fields.get("kind") or fields.get("mode")
    samples = [float(piece) for piece in fields["samples_ms"].split(",") if piece]
    return {
        "kind": fields["kind"],
        "molecules": int(fields["molecules"]),
        "threads": int(fields["threads"]),
        "sink": int(float(fields["sink"])),
        "samples": samples,
    }


def row_for(cell: str, parsed: dict, summary: dict, cores: int, affinity: str) -> dict:
    median = summary["median_ms"]
    return {
        "cell": cell,
        "code": parsed["kind"],
        "molecules": parsed["molecules"],
        "threads": parsed["threads"],
        "cores": cores,
        "affinity": affinity,
        "warmup": 1,
        "n": summary["n"],
        "median_ms": f"{median:.6f}",
        "iqr_ms": f"{summary['iqr_ms']:.6f}",
        "q1_ms": f"{summary['q1_ms']:.6f}",
        "q3_ms": f"{summary['q3_ms']:.6f}",
        "us_per_atom": f"{median * 1000.0 / parsed['molecules']:.6f}",
        "sink": parsed["sink"],
        "samples_ms": stats.fmt_samples(parsed["samples"]),
    }


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--out", type=Path, default=ROOT / "benchmarks/paper-cells")
    parser.add_argument("--bench", default=str(ROOT / "bbdir-release/tests/bench_v2"))
    parser.add_argument("--threads", type=int, default=8)
    parser.add_argument("--repeats", type=int, default=8)
    args = parser.parse_args()
    if args.repeats < 8:
        raise SystemExit("n must be at least 8")
    if not Path(args.bench).is_file():
        raise SystemExit(f"missing {args.bench}")
    args.out.mkdir(parents=True, exist_ok=True)
    raw = args.out / "raw-ptm.txt"
    env = os.environ.copy()
    env["NEIGH_ONLY"] = f"symmetry:{MOLECULES}:{args.threads}"
    env["NEIGH_BENCH"] = args.bench
    done = subprocess.run(
        [sys.executable, str(ROOT / "scripts/bench_neighbours_v2.py"), str(raw)],
        cwd=ROOT, env=env, text=True, capture_output=True,
    )
    if done.returncode != 0:
        sys.stderr.write(done.stdout)
        sys.stderr.write(done.stderr)
        raise SystemExit(f"bench_neighbours_v2 exited {done.returncode}")
    lines = [
        line for line in raw.read_text(encoding="utf-8").splitlines()
        if "samples_ms=" in line
    ]
    if len(lines) != 2:
        raise SystemExit(f"expected two jobs, got {len(lines)}")
    parsed = [parse_line(line) for line in lines]
    by_kind = {item["kind"]: item for item in parsed}
    if set(by_kind) != {"symmetry", "ovito-symmetry"}:
        raise SystemExit(f"unexpected jobs {sorted(by_kind)}")
    cores, affinity = stats.cores_now()
    rows = []
    for kind, cell in (
        ("symmetry", "classifier against PTM"),
        ("ovito-symmetry", "PTM diamond templates"),
    ):
        item = by_kind[kind]
        if item["molecules"] != MOLECULES or item["threads"] != args.threads:
            raise SystemExit(f"{kind} was not the 110592 / {args.threads} job")
        if len(item["samples"]) != args.repeats:
            raise SystemExit(f"{kind} wrote {len(item['samples'])} samples")
        if item["sink"] != MOLECULES:
            raise SystemExit(f"{kind} labelled {item['sink']} of {MOLECULES}")
        summary = stats.summarize(item["samples"])
        rows.append(row_for(cell, item, summary, cores, affinity))
        print(f"{cell} {summary['median_ms']:.4f} ms iqr {summary['iqr_ms']:.4f}", flush=True)
    header = [
        f"cores {cores}",
        f"affinity {affinity}",
        "ptm rmsd_cutoff 0.1 templates cubic-diamond hexagonal-diamond",
        "classifier bench_v2 symmetry cutoff-csr diamondFromCsr",
        "sitting one",
    ]
    raw.write_text("\n".join(header + lines) + "\n", encoding="utf-8")
    sym = args.out / "tab-sym-eq.csv"
    with sym.open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(rows[0]))
        writer.writeheader()
        writer.writerows(rows)
    ptm = rows[1]
    frag = args.out / "frag-cli-ptm.csv"
    with frag.open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=[
            "cell", "code", "molecules", "threads", "cores", "affinity",
            "warmup", "n", "median_ms", "iqr_ms", "samples_ms", "detail",
        ])
        writer.writeheader()
        writer.writerow({
            "cell": "PTM diamond templates",
            "code": "OVITO PTM",
            "molecules": ptm["molecules"],
            "threads": ptm["threads"],
            "cores": ptm["cores"],
            "affinity": ptm["affinity"],
            "warmup": ptm["warmup"],
            "n": ptm["n"],
            "median_ms": ptm["median_ms"],
            "iqr_ms": ptm["iqr_ms"],
            "samples_ms": ptm["samples_ms"],
            "detail": "rmsd 0.1, templates 6 and 7; IQR is q3-q1",
        })
    print(f"wrote {sym}", flush=True)
    return 0


if __name__ == "__main__":
    sys.exit(main())
