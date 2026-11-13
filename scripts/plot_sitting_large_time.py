#!/usr/bin/env python3
"""Time against threads for the 884736-molecule sitting.

The vertical axis is logarithmic. Each operation is its own panel and
each code is a line. The PDF is written where --out says and is not
part of the source tree. The host name is not read.
"""

from __future__ import annotations

import argparse
import csv
from collections import defaultdict
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]


def load(path: Path) -> list[dict]:
    with path.open(newline="") as handle:
        return list(csv.DictReader(handle))


def numeric(text: str):
    try:
        return float(text)
    except (TypeError, ValueError):
        return None


def draw(rows: list[dict], dest: Path) -> None:
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    curves: dict[str, dict[tuple[str, str], dict[int, float]]] = defaultdict(lambda: defaultdict(dict))
    for row in rows:
        if str(row.get("molecules")) != "884736":
            continue
        value = numeric(row.get("median_ms", ""))
        if value is None or value <= 0:
            continue
        phase = row.get("phase") or ""
        curves[row["operation"]][(phase, row["code"])][int(row["threads"])] = value
    operations = [op for op in curves if curves[op]]
    if not operations:
        raise SystemExit("no 884736 timings to plot")
    colors = {
        "d-SEAMS": "#0072B2",
        "rings-statistics-matter": "#E69F00",
        "mdapy": "#009E73",
        "vesin": "#CC79A7",
        "freud": "#D55E00",
        "OVITO": "#56B4E9",
        "pyscal3": "#000000",
        "PLUMED": "#F0E442",
    }
    fig, axes = plt.subplots(
        len(operations), 1, figsize=(7.2, 1.7 * len(operations)), sharex=True, squeeze=False)
    commit = rows[0].get("source_commit", "")[:12]
    fig.suptitle(f"884736 molecules, time against threads, commit {commit}", fontsize=11)
    widths = [1, 2, 4, 8]
    for index, operation in enumerate(operations):
        axis = axes[index][0]
        axis.set_yscale("log")
        for (phase, code), times in sorted(curves[operation].items()):
            xs = [width for width in widths if width in times]
            if not xs:
                continue
            label = "{} {}".format(phase, code).strip()
            axis.plot(
                xs, [times[width] for width in xs],
                marker="o", color=colors.get(code, "#333333"),
                linestyle="--" if phase == "liquid" else "-",
                label=label,
            )
        axis.set_ylabel("ms")
        axis.set_title(operation, loc="left", fontsize=10)
        axis.legend(fontsize=7, loc="upper right", frameon=False, ncol=2)
    axes[-1][0].set_xlabel("threads")
    axes[-1][0].set_xticks(widths)
    fig.tight_layout()
    dest.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(dest)
    print(f"wrote {dest}")


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--csv", type=Path, default=ROOT / "benchmarks/paper-cells/sitting-large.csv")
    parser.add_argument(
        "--out", type=Path,
        default=ROOT / "benchmarks/paper-cells/sitting-large-time.pdf")
    args = parser.parse_args()
    draw(load(args.csv), args.out)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
