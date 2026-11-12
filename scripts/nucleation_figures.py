#!/usr/bin/env python3
# Copyright (c) d-SEAMS developers

"""Plot a nucleation-pathway TSV.

The table is the one written by nucleation-pathway. Time is
(frame - first) times --ps-per-frame. The default, 10, is one
1000-step dump interval at the 10 fs mW step.
"""

from __future__ import annotations

import argparse
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt


def load(path: Path):
    lines = path.read_text().splitlines()
    header = lines[0].split("\t")
    rows = []
    for line in lines[1:]:
        if not line.strip():
            continue
        parts = line.split("\t")
        row = {}
        for key, value in zip(header, parts):
            row[key] = float(value)
        rows.append(row)
    return rows


def series(rows, key):
    return [row[key] for row in rows]


def fraction(rows, keys):
    out = []
    for row in rows:
        den = row["nop"]
        num = sum(row[key] for key in keys)
        out.append(num / den if den else 0.0)
    return out


def stack(rows, cubic, hexagon):
    out = []
    for row in rows:
        den = row[cubic] + row[hexagon]
        out.append(row[cubic] / den if den else float("nan"))
    return out


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("tsv", type=Path)
    parser.add_argument("--out", type=Path, required=True)
    parser.add_argument("--ps-per-frame", type=float, default=10.0)
    parser.add_argument("--title", default="")
    args = parser.parse_args()
    rows = load(args.tsv)
    if not rows:
        raise SystemExit(f"{args.tsv} has no frames")
    t0 = rows[0]["frame"]
    time = [(row["frame"] - t0) * args.ps_per_frame for row in rows]
    ring = []
    for row in rows:
        den = row["ring5"] + row["ring6"]
        ring.append(row["ring6"] / den if den else 0.0)

    fig, axes = plt.subplots(3, 1, figsize=(7.2, 8.4), sharex=True)
    axes[0].plot(time, ring, color="black", label="n6 / (n5 + n6)")
    axes[0].plot(time, fraction(rows, ["chill_cubic", "chill_hex", "chill_inter"]),
                 label="CHILL+ extended")
    axes[0].plot(time, fraction(rows, ["ovito_cubic", "ovito_hex", "ovito_inter"]),
                 label="ovito-chill+ extended")
    axes[0].plot(time, fraction(rows, ["q6bar_ice"]), label="q6-bar")
    axes[0].plot(time, fraction(rows, ["sym_cubic", "sym_hex"]), label="diamond")
    axes[0].set_ylabel("fraction")
    axes[0].legend(fontsize=8, loc="best")

    axes[1].plot(time, series(rows, "n_chill_core"), label="CHILL+ core")
    axes[1].plot(time, series(rows, "n_chill_ext"), label="CHILL+ extended")
    axes[1].plot(time, series(rows, "n_ovito_core"), label="ovito core")
    axes[1].plot(time, series(rows, "n_ovito_ext"), label="ovito extended")
    axes[1].plot(time, series(rows, "n_q6"), label="q6-bar")
    axes[1].plot(time, series(rows, "n_sym"), label="diamond")
    axes[1].set_ylabel("largest cluster")
    axes[1].legend(fontsize=8, loc="best")

    axes[2].plot(time, stack(rows, "c_chill_ext", "h_chill_ext"), label="CHILL+ extended")
    axes[2].plot(time, stack(rows, "c_ovito_ext", "h_ovito_ext"), label="ovito extended")
    axes[2].plot(time, stack(rows, "c_q6", "h_q6"), label="q6-bar, by diamond")
    axes[2].plot(time, stack(rows, "c_sym", "h_sym"), label="diamond")
    axes[2].set_ylabel("Ic / (Ic + Ih)")
    axes[2].set_xlabel("time / ps")
    axes[2].set_ylim(-0.05, 1.05)
    axes[2].legend(fontsize=8, loc="upper right")
    if args.title:
        axes[0].set_title(args.title)
    fig.tight_layout()
    args.out.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(args.out, dpi=140)
    print(args.out)


if __name__ == "__main__":
    main()
