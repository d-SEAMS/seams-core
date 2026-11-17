#!/usr/bin/env python3
"""How the three cluster series of each climb move together.

For every climb directory (one seed, with frames.tsv from
nucleation-pathway) this writes the Pearson and Spearman correlation
of n_chill_core, n_chill_ext and n_q6solid over the frames from 1 ns,
and the Pearson correlation of their 100 ps increments. A slow common
drift makes the levels correlate even when the series do not track
each other; the increments do not carry that drift. The maxima and the
last values are written next to them.

    python3 scripts/committor_order.py \\
        data/science-frames/committor-climb-0-of-20/mw-4096/seed-* \\
        --out data/committor-shards/order-0-of-20.csv
"""

from __future__ import annotations

import argparse
import csv
import math
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "scripts"))
import committor_shoot as cs  # noqa: E402

EQUIL_PS = cs.EQUIL_PS
STEP_PS = 100.0
PAIRS = (("n_chill_core", "n_chill_ext"), ("n_chill_core", "n_q6solid"),
         ("n_chill_ext", "n_q6solid"))
SHORT = {"n_chill_core": "core", "n_chill_ext": "ext", "n_q6solid": "q6solid"}


def pearson(x: list[float], y: list[float]) -> float:
    n = len(x)
    if n < 3:
        return float("nan")
    mx = sum(x) / n
    my = sum(y) / n
    sxy = sum((a - mx) * (b - my) for a, b in zip(x, y))
    sxx = sum((a - mx) ** 2 for a in x)
    syy = sum((b - my) ** 2 for b in y)
    if sxx <= 0 or syy <= 0:
        return float("nan")
    return sxy / math.sqrt(sxx * syy)


def ranks(values: list[float]) -> list[float]:
    order = sorted(range(len(values)), key=lambda i: values[i])
    out = [0.0] * len(values)
    i = 0
    while i < len(order):
        j = i
        while j + 1 < len(order) and values[order[j + 1]] == values[order[i]]:
            j += 1
        mean = (i + j) / 2.0
        for k in range(i, j + 1):
            out[order[k]] = mean
        i = j + 1
    return out


def spearman(x: list[float], y: list[float]) -> float:
    return pearson(ranks(x), ranks(y))


def increments(values: list[float], times: list[float], step: float) -> list[float]:
    out = []
    j = 0
    for i, t in enumerate(times):
        while j < len(times) and times[j] < t + step:
            j += 1
        if j >= len(times):
            break
        if abs(times[j] - t - step) < 1e-6:
            out.append(values[j] - values[i])
    return out


def climb_order(directory: Path) -> dict[str, str]:
    rows = cs.read_frames(directory / "frames.tsv")
    times = cs.frame_times(rows)
    keep = [i for i, t in enumerate(times) if t >= EQUIL_PS]
    series = {name: cs.column(rows, name) for name in cs.SERIES}
    out = {
        "seed": directory.name.removeprefix("seed-"),
        "frames": str(len(keep)),
        "t_lo_ps": cs.fmt(times[keep[0]]) if keep else "",
        "t_hi_ps": cs.fmt(times[keep[-1]]) if keep else "",
    }
    for name in cs.SERIES:
        values = series[name]
        whole = all(v == v for v in values)
        out[f"{SHORT[name]}_max"] = cs.fmt(max(values)) if whole else ""
        out[f"{SHORT[name]}_end"] = cs.fmt(values[-1])
    kept_times = [times[i] for i in keep]
    for a, b in PAIRS:
        x = [series[a][i] for i in keep]
        y = [series[b][i] for i in keep]
        key = f"{SHORT[a]}_{SHORT[b]}"
        if not keep or any(v != v for v in x + y):
            out[f"pearson_{key}"] = out[f"spearman_{key}"] = out[f"dpearson_{key}"] = ""
            continue
        out[f"pearson_{key}"] = f"{pearson(x, y):.4f}"
        out[f"spearman_{key}"] = f"{spearman(x, y):.4f}"
        dx = increments(x, kept_times, STEP_PS)
        dy = increments(y, kept_times, STEP_PS)
        out[f"dpearson_{key}"] = f"{pearson(dx, dy):.4f}"
    return out


def fields() -> list[str]:
    names = ["seed", "frames", "t_lo_ps", "t_hi_ps"]
    for name in cs.SERIES:
        names += [f"{SHORT[name]}_max", f"{SHORT[name]}_end"]
    for kind in ("pearson", "spearman", "dpearson"):
        names += [f"{kind}_{SHORT[a]}_{SHORT[b]}" for a, b in PAIRS]
    return names


def self_test() -> None:
    if abs(pearson([1, 2, 3, 4], [2, 4, 6, 8]) - 1.0) > 1e-12:
        raise SystemExit("pearson of a line is not 1")
    if abs(spearman([1, 2, 3, 4], [1, 8, 27, 64]) - 1.0) > 1e-12:
        raise SystemExit("spearman of a monotone map is not 1")
    if ranks([3, 1, 3, 2]) != [2.5, 0.0, 2.5, 1.0]:
        raise SystemExit(f"ties {ranks([3, 1, 3, 2])}")
    times = [10.0 * i for i in range(30)]
    if increments([float(i) for i in range(30)], times, 100.0) != [10.0] * 20:
        raise SystemExit("100 ps increments of a ramp are not 10")
    print("self-test ok")


def main() -> None:
    if "--self-test" in sys.argv:
        self_test()
        return
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("climbs", nargs="+", type=Path)
    parser.add_argument("--out", type=Path, required=True)
    parser.add_argument("--commit", default="")
    args = parser.parse_args()
    commit = args.commit or cs.head_commit()
    rows = [climb_order(d) for d in sorted(args.climbs) if (d / "frames.tsv").is_file()]
    args.out.parent.mkdir(parents=True, exist_ok=True)
    with args.out.open("w", newline="") as handle:
        handle.write(f"# commit {commit}\n")
        handle.write(f"# frames from {cs.fmt(EQUIL_PS)} ps; increments over "
                     f"{cs.fmt(STEP_PS)} ps; order parameter {cs.ORDER}\n")
        writer = csv.DictWriter(handle, fieldnames=fields(), lineterminator="\n")
        writer.writeheader()
        writer.writerows(rows)
    print(args.out.read_text())


if __name__ == "__main__":
    main()
