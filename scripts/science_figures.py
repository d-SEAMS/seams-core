#!/usr/bin/env python3
# Copyright (c) d-SEAMS developers

"""Publication figures for the basal interface and the nucleation census.

Reads the CSV files written by science-front. Writes PDF only.
"""

from __future__ import annotations

import argparse
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np

ROOT = Path(__file__).resolve().parents[1]
DATA = ROOT / "paper" / "tex" / "figures" / "data" / "science"

CHILL = {
    0: ("cubic", "#0072B2"),
    1: ("hexagonal", "#E69F00"),
    3: ("interfacial", "#009E73"),
    2: ("liquid", "#C8C8C8"),
    4: ("clathrate", "#56B4E9"),
    5: ("interfacial clathrate", "#CC79A7"),
}


def read_table(path: Path):
    header = None
    rows = []
    for line in path.read_text().splitlines():
        if not line or line.startswith("#"):
            continue
        parts = line.split()
        if header is None:
            header = parts
            continue
        rows.append(parts)
    return header, rows


def column(header, rows, name: str) -> np.ndarray:
    index = header.index(name)
    return np.array([float(row[index]) for row in rows])


def style() -> None:
    plt.rcParams.update(
        {
            "font.size": 9,
            "axes.labelsize": 10,
            "axes.linewidth": 0.6,
            "xtick.direction": "in",
            "ytick.direction": "in",
            "pdf.fonttype": 42,
            "ps.fonttype": 42,
            "font.family": "DejaVu Sans",
        }
    )


def save(fig, path: Path) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    fig.tight_layout()
    fig.savefig(path, format="pdf")
    plt.close(fig)


def snapshot(path: Path) -> None:
    header, rows = read_table(DATA / "interface-atoms.csv")
    layer_header, layer_rows = read_table(DATA / "interface-layers.csv")
    axis = int(column(layer_header, layer_rows, "axis")[0])
    coords = [column(header, rows, name) for name in ("x", "y", "z")]
    normal = coords[axis]
    lateral = coords[(axis + 1) % 3]
    labels = column(header, rows, "chill").astype(int)
    fig, ax = plt.subplots(figsize=(7.0, 4.4))
    for code, (name, colour) in CHILL.items():
        mask = labels == code
        if not np.any(mask):
            continue
        ax.scatter(
            lateral[mask],
            normal[mask],
            s=8,
            c=colour,
            linewidths=0,
            label=name,
            rasterized=True,
        )
    ax.set_xlabel("lateral coordinate (A)")
    ax.set_ylabel("basal normal (A)")
    ax.set_aspect("equal", adjustable="box")
    ax.legend(frameon=False, loc="upper right", markerscale=2)
    save(fig, path)


def profiles(path: Path) -> None:
    header, rows = read_table(DATA / "interface-layers.csv")
    frame = column(header, rows, "frame")
    first = frame == frame.min()
    last = frame == frame.max()
    front_header, front_rows = read_table(DATA / "interface-fronts.csv")
    vel_header, vel_rows = read_table(DATA / "interface-velocity.csv")
    fig, axes = plt.subplots(1, 2, figsize=(7.2, 3.6))
    ax = axes[0]
    for mask, style_name in ((first, "-"), (last, "--")):
        mid = column(header, rows, "mid")[mask]
        n = column(header, rows, "n")[mask]
        solid = (
            column(header, rows, "cubic")[mask]
            + column(header, rows, "hexagonal")[mask]
            + column(header, rows, "interfacial")[mask]
        )
        frac = np.divide(solid, n, out=np.zeros_like(solid), where=n > 0)
        q6 = column(header, rows, "q6bar")[mask]
        label = "first" if style_name == "-" else "last"
        ax.plot(mid, frac, style_name, color="#0072B2", label=f"CHILL+ solid, {label}")
        ax.plot(mid, q6, style_name, color="#D55E00", label=f"q6-bar, {label}")
    ax.set_xlabel("basal normal (A)")
    ax.set_ylabel("solid fraction and q6-bar")
    ax.set_ylim(-0.05, 1.05)
    ax.legend(frameon=False, fontsize=7)
    ax = axes[1]
    time = column(front_header, front_rows, "time_ps")
    ax.plot(time, column(front_header, front_rows, "pos0"), color="#0072B2", label="front 0")
    ax.plot(time, column(front_header, front_rows, "pos1"), color="#E69F00", label="front 1")
    ax.plot(time, column(front_header, front_rows, "width0"), color="#0072B2", ls="--", label="width 0")
    ax.plot(time, column(front_header, front_rows, "width1"), color="#E69F00", ls="--", label="width 1")
    vel = {row[0]: row for row in vel_rows}
    if "mean_front" in vel:
        slope = float(vel["mean_front"][vel_header.index("slope_A_per_ps")])
        lo = float(vel["mean_front"][vel_header.index("lo")])
        hi = float(vel["mean_front"][vel_header.index("hi")])
        ax.text(
            0.04,
            0.04,
            f"mean front {slope:.3f} A/ps\n[{lo:.3f}, {hi:.3f}]",
            transform=ax.transAxes,
            fontsize=8,
            va="bottom",
        )
    ax.set_xlabel("time (ps)")
    ax.set_ylabel("position and width (A)")
    ax.legend(frameon=False, fontsize=7)
    save(fig, path)


def clusters(path: Path) -> None:
    header, rows = read_table(DATA / "nucleation-clusters.csv")
    summary_header, summary_rows = read_table(DATA / "nucleation-summary.csv")
    frame = column(header, rows, "frame")
    fig, ax = plt.subplots(figsize=(7.0, 4.0))
    series = (
        ("n_ext", "CHILL+ extended", "#0072B2"),
        ("n_core", "CHILL+ core", "#E69F00"),
        ("n_q6", "q6-bar", "#009E73"),
        ("n_shell", "shell label", "#CC79A7"),
    )
    for name, label, colour in series:
        ax.plot(frame, column(header, rows, name), color=colour, label=label)
    for row in summary_rows:
        if row[0] != "chill_extended":
            continue
        onset = int(float(row[summary_header.index("frame")]))
        n_star = float(row[summary_header.index("n_star")])
        lo = float(row[summary_header.index("lo")])
        hi = float(row[summary_header.index("hi")])
        ax.axhline(n_star, color="#0072B2", lw=0.7, ls=":")
        ax.axhspan(lo, hi, color="#0072B2", alpha=0.12)
        ax.axvline(onset, color="#333333", lw=0.6, ls=":")
        bound = int(float(row[summary_header.index("upper_bound")])) == 1
        title = "upper bound" if bound else "critical size"
        ax.text(
            0.98,
            0.04,
            f"{title} {n_star:.1f} [{lo:.1f}, {hi:.1f}]",
            transform=ax.transAxes,
            ha="right",
            va="bottom",
            fontsize=8,
        )
    ax.set_xlabel("frame")
    ax.set_ylabel("largest cluster (molecules)")
    ax.legend(frameon=False)
    save(fig, path)


def read_literature():
    path = DATA / "growth-literature.csv"
    rows = []
    if not path.is_file():
        return rows
    header = None
    for line in path.read_text().splitlines():
        if not line or line.startswith("#"):
            continue
        parts = line.split("|")
        if header is None:
            header = parts
            continue
        rows.append(dict(zip(header, parts)))
    return rows


def speeds(path: Path) -> None:
    header, rows = read_table(DATA / "interface-speed.csv")
    temp = column(header, rows, "temperature_K")
    rate = column(header, rows, "mean_front_A_per_ns")
    lo = column(header, rows, "lo_A_per_ns")
    hi = column(header, rows, "hi_A_per_ns")
    fig, ax = plt.subplots(figsize=(4.6, 3.6))
    ax.errorbar(
        temp,
        rate,
        yerr=np.vstack([rate - lo, hi - rate]),
        fmt="o",
        color="#0072B2",
        ms=5,
        lw=1,
        capsize=2,
        label="this work, basal, restrained ice",
    )
    for row in read_literature():
        if row.get("kind") != "measured" or not row.get("temperature_K") or not row.get("rate_A_per_ns"):
            continue
        sigma = float(row["rate_sigma"]) if row.get("rate_sigma") else 0.0
        ax.errorbar(
            [float(row["temperature_K"])],
            [float(row["rate_A_per_ns"])],
            yerr=[[sigma], [sigma]] if sigma else None,
            fmt="s",
            color="#D55E00",
            ms=5,
            capsize=2,
            label="Shi 2024 basal maximum",
        )
    ax.set_xlabel("temperature (K)")
    ax.set_ylabel("front speed (A/ns)")
    ax.legend(frameon=False, fontsize=7)
    save(fig, path)


def nucleus(path: Path) -> None:
    header, rows = read_table(DATA / "nucleus-atoms.csv")
    keep = column(header, rows, "in_extended") > 0.5
    if not np.any(keep):
        keep = np.ones(len(rows), dtype=bool)
    labels = column(header, rows, "chill").astype(int)[keep]
    shell = column(header, rows, "shell").astype(int)[keep]
    x = column(header, rows, "x")[keep]
    y = column(header, rows, "y")[keep]
    z = column(header, rows, "z")[keep]
    fig, axes = plt.subplots(1, 2, figsize=(7.2, 3.8))
    for ax, lateral, normal, xlabel, ylabel in (
        (axes[0], x, y, "x (A)", "y (A)"),
        (axes[1], x, z, "x (A)", "z (A)"),
    ):
        for code, (name, colour) in CHILL.items():
            mask = labels == code
            if not np.any(mask):
                continue
            edge = np.where(shell[mask] == 2, "#222222", "none")
            ax.scatter(
                lateral[mask],
                normal[mask],
                s=10,
                c=colour,
                linewidths=0.4,
                edgecolors=edge,
                label=name,
                rasterized=True,
            )
        ax.set_xlabel(xlabel)
        ax.set_ylabel(ylabel)
        ax.set_aspect("equal", adjustable="datalim")
    handles, labels_legend = axes[1].get_legend_handles_labels()
    axes[1].legend(handles, labels_legend, frameon=False, fontsize=7, loc="best")
    save(fig, path)


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--out", type=Path, required=True)
    args = parser.parse_args()
    style()
    out = args.out
    profiles(out / "interface-profiles.pdf")
    clusters(out / "cluster-size.pdf")
    if (DATA / "interface-speed.csv").is_file():
        speeds(out / "growth-rate.pdf")
    print(out)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
