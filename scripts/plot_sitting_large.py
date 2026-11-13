#!/usr/bin/env python3
"""Strong scaling of the 884736-molecule sitting.

Speedup is the one-thread median divided by the median at that width.
Efficiency is the speedup divided by the width. A cell that is not a
time (a rival that ran out of memory, or a call over 120 s) is left
off the curve. The host name is not read.
"""

from __future__ import annotations

import argparse
import csv
from collections import defaultdict
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]

OPERATIONS = (
    "rings",
    "q4",
    "q6",
    "q4+q6",
    "neighbours",
    "chill+",
    "symmetry",
    "ptm-default",
    "ptm-ice",
)


def load(path: Path) -> list[dict]:
    with path.open(newline="") as handle:
        return list(csv.DictReader(handle))


def numeric(text: str):
    try:
        return float(text)
    except (TypeError, ValueError):
        return None


def series(rows: list[dict]) -> dict[tuple[str, str], dict[int, float]]:
    found: dict[tuple[str, str], dict[int, float]] = defaultdict(dict)
    for row in rows:
        value = numeric(row.get("median_ms", ""))
        if value is None or value <= 0:
            continue
        if row.get("operation") not in OPERATIONS:
            continue
        found[(row["operation"], row["code"])][int(row["threads"])] = value
    return found


def draw(rows: list[dict], dest: Path) -> None:
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    curves = series(rows)
    operations = [op for op in OPERATIONS if any(key[0] == op for key in curves)]
    if not operations:
        raise SystemExit("no timed rows to plot")
    colors = {
        "d-SEAMS": "#0072B2",
        "rings-statistics-matter": "#E69F00",
        "mdapy": "#009E73",
        "vesin": "#CC79A7",
        "freud": "#D55E00",
        "OVITO": "#56B4E9",
    }
    fig, axes = plt.subplots(
        len(operations), 2, figsize=(9.2, 2.15 * len(operations)), sharex=True, squeeze=False)
    machine = rows[0].get("cpu_model", "")
    commit = rows[0].get("source_commit", "")[:12]
    fig.suptitle(
        f"884736 molecules, ideal ice Ic\n{machine}, commit {commit}",
        fontsize=11)
    widths = [1, 2, 4, 8]
    for row_index, operation in enumerate(operations):
        speed = axes[row_index][0]
        eff = axes[row_index][1]
        speed.plot(widths, widths, color="#999999", linewidth=0.8, label="ideal")
        eff.axhline(1.0, color="#999999", linewidth=0.8)
        codes = sorted({code for (op, code) in curves if op == operation})
        for code in codes:
            times = curves[(operation, code)]
            base = times.get(1)
            if base is None:
                continue
            xs = [width for width in widths if width in times]
            ys = [base / times[width] for width in xs]
            speed.plot(xs, ys, marker="o", color=colors.get(code, "#000000"), label=code)
            eff.plot(
                xs, [value / width for value, width in zip(ys, xs)],
                marker="o", color=colors.get(code, "#000000"), label=code)
        speed.set_ylabel("speedup")
        eff.set_ylabel("efficiency")
        speed.set_title(operation, loc="left", fontsize=10)
        if row_index == 0:
            speed.legend(fontsize=7, loc="upper left", frameon=False)
    for axis in axes[-1]:
        axis.set_xlabel("threads")
        axis.set_xticks(widths)
    fig.tight_layout()
    dest.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(dest, dpi=140)
    print(f"wrote {dest}")


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--csv", type=Path,
        default=ROOT / "benchmarks/paper-cells/sitting-large.csv")
    parser.add_argument(
        "--out", type=Path,
        default=ROOT / "benchmarks/paper-cells/sitting-large-scaling.png")
    args = parser.parse_args()
    draw(load(args.csv), args.out)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
