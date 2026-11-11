#!/usr/bin/env python3
"""Scaling figure for the spatial-domain split.

Reads benchmarks/rank-scale.csv and, when present,
benchmarks/h2h-parallel.csv. Writes PDF and SVG next to the paper figures.
"""

from __future__ import annotations

import csv
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt

ROOT = Path(__file__).resolve().parents[1]
RANK = ROOT / "benchmarks" / "rank-scale.csv"
H2H = ROOT / "benchmarks" / "h2h-parallel.csv"
OUT = ROOT / "paper" / "tex" / "figures"

# Okabe-Ito
COLOURS = {
    "neigh": "#0072B2",
    "chill": "#D55E00",
    "q6": "#009E73",
    "rings": "#CC79A7",
}
LABELS = {
    "neigh": "neighbours",
    "chill": "CHILL+",
    "q6": "q6",
    "rings": "rings",
}


def best_rank_rows(molecules: int):
    rows = []
    with RANK.open() as handle:
        for row in csv.DictReader(handle):
            if int(row["molecules"]) != molecules:
                continue
            if int(row["threads_per_rank"]) != 1:
                continue
            rows.append(row)
    return rows


def main():
    fig, axes = plt.subplots(1, 2, figsize=(7.2, 3.4), constrained_layout=True)
    for ax, molecules, title in (
        (axes[0], 110592, "110592 molecules"),
        (axes[1], 884736, "884736 molecules"),
    ):
        by_op = {}
        for row in best_rank_rows(molecules):
            by_op.setdefault(row["operation"], []).append(row)
        ax.plot([1, 8], [1, 8], color="#000000", lw=0.8, ls="--", label="ideal")
        for op in ("neigh", "chill", "q6", "rings"):
            series = sorted(by_op.get(op, []), key=lambda r: int(r["ranks"]))
            if not series:
                continue
            xs = [int(r["ranks"]) for r in series]
            ys = [float(r["speedup"]) for r in series]
            ax.plot(xs, ys, marker="o", color=COLOURS[op], label=LABELS[op], lw=1.4)
        ax.set_xlabel("Ranks (1 thread each)")
        ax.set_ylabel("Speedup vs 1 rank")
        ax.set_title(title)
        ax.set_xticks([1, 2, 4, 8])
        ax.set_ylim(0, 9)
        ax.spines["top"].set_visible(False)
        ax.spines["right"].set_visible(False)
    axes[0].legend(frameon=False, fontsize=8)
    pdf = OUT / "strong-scale.pdf"
    svg = OUT / "strong-scale.svg"
    fig.savefig(pdf)
    fig.savefig(svg)
    print(f"wrote {pdf}")
    print(f"wrote {svg}")


if __name__ == "__main__":
    main()
