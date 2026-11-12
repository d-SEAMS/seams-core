#!/usr/bin/env python3
"""Front velocity against temperature from coexistence CSVs.

Each CSV is one trajectory written by front_csv.py. The ordinary
least-squares velocity is half the slope of ice thickness against
time on every frame after --skip-ns. That standard error treats
frames as independent.

The gate uses a coarser error. The trajectory is split into
half-nanosecond windows, each window is fitted on its own, and the
reported uncertainty is the standard error of those window
velocities. The gate passes only when the cold side grows and the
warm side melts, each by more than twice that standard error.
Two seeds at the same temperature are averaged before the test.
"""

from __future__ import annotations

import argparse
import math
import re
import sys
from pathlib import Path


def read_series(path: Path):
    rows = []
    for line in path.read_text().splitlines():
        if not line or line.startswith("#") or line.startswith("t_ns"):
            continue
        parts = line.split(",")
        rows.append((float(parts[0]), float(parts[1]), float(parts[2])))
    return rows


def slope(xs, ys):
    n = len(xs)
    if n < 3:
        return float("nan"), float("nan")
    mx = sum(xs) / n
    my = sum(ys) / n
    ssxx = sum((x - mx) ** 2 for x in xs)
    ssxy = sum((x - mx) * (y - my) for x, y in zip(xs, ys))
    if ssxx <= 0:
        return float("nan"), float("nan")
    b = ssxy / ssxx
    sse = sum((y - (my + b * (x - mx))) ** 2 for x, y in zip(xs, ys))
    err = math.sqrt((sse / (n - 2)) / ssxx)
    return b, err


def temperature(path: Path) -> float:
    match = re.search(r"(\d+(?:\.\d+)?)", path.stem)
    if not match:
        raise SystemExit(f"no temperature in {path.name}")
    return float(match.group(1))


def window_velocities(rows, width):
    if not rows:
        return []
    start = rows[0][0]
    end = rows[-1][0]
    velocities = []
    left = start
    while left < end - 1e-9:
        right = min(left + width, end)
        last = abs(right - end) < 1e-9
        block = [
            row
            for row in rows
            if (left <= row[0] < right) or (last and abs(row[0] - end) < 1e-9)
        ]
        if len(block) >= 3:
            thick, _err = slope([row[0] for row in block], [row[2] for row in block])
            velocities.append(thick / 2.0)
        left = right
    return velocities


def mean_se(values):
    n = len(values)
    if n == 0:
        return float("nan"), float("nan")
    mean = sum(values) / n
    if n == 1:
        return mean, float("nan")
    var = sum((value - mean) ** 2 for value in values) / (n * (n - 1))
    return mean, math.sqrt(var)


def one_run(path: Path, skip: float, width: float):
    rows = read_series(path)
    kept = [row for row in rows if row[0] >= skip]
    thick, thick_err = slope([row[0] for row in kept], [row[2] for row in kept])
    dn, dn_err = slope([row[0] for row in kept], [row[1] for row in kept])
    block_v, block_se = mean_se(window_velocities(rows, width))
    return {
        "path": path,
        "T": temperature(path),
        "n": len(kept),
        "t_end": rows[-1][0] if rows else float("nan"),
        "v": thick / 2.0,
        "ve": thick_err / 2.0,
        "dn": dn,
        "dne": dn_err,
        "block_v": block_v,
        "block_se": block_se,
        "n_block": len(window_velocities(rows, width)),
    }


def mean_group(runs):
    n = len(runs)
    v = sum(run["v"] for run in runs) / n
    dn = sum(run["dn"] for run in runs) / n
    block_v = sum(run["block_v"] for run in runs) / n
    if n == 1:
        ve, dne = runs[0]["ve"], runs[0]["dne"]
        block_se = runs[0]["block_se"]
        kind = "slope"
    else:
        ve = math.sqrt(sum((run["v"] - v) ** 2 for run in runs) / (n * (n - 1)))
        dne = math.sqrt(sum((run["dn"] - dn) ** 2 for run in runs) / (n * (n - 1)))
        block_se = math.sqrt(
            sum((run["block_v"] - block_v) ** 2 for run in runs) / (n * (n - 1))
        )
        kind = "replicate"
    return {
        "T": runs[0]["T"],
        "nseed": n,
        "v": v,
        "ve": ve,
        "dn": dn,
        "dne": dne,
        "block_v": block_v,
        "block_se": block_se,
        "kind": kind,
        "t_end": min(run["t_end"] for run in runs),
    }


def crossing(cold, warm, key, err_key):
    dv = warm[key] - cold[key]
    if dv == 0 or math.isnan(dv):
        return float("nan"), float("nan")
    span = warm["T"] - cold["T"]
    t0 = cold["T"] - cold[key] * span / dv
    d_dv1 = -span * warm[key] / dv**2
    d_dv2 = span * cold[key] / dv**2
    ce = 0.0 if math.isnan(cold[err_key]) else cold[err_key]
    we = 0.0 if math.isnan(warm[err_key]) else warm[err_key]
    err = math.sqrt((d_dv1 * ce) ** 2 + (d_dv2 * we) ** 2)
    return t0, err


def grows(point) -> bool:
    return point["block_v"] > 2.0 * point["block_se"]


def melts(point) -> bool:
    return point["block_v"] < -2.0 * point["block_se"]


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("csvs", nargs="+", type=Path)
    parser.add_argument("--skip-ns", type=float, default=0.06)
    parser.add_argument("--block-ns", type=float, default=0.5)
    args = parser.parse_args()
    runs = [one_run(path, args.skip_ns, args.block_ns) for path in args.csvs]
    by_t = {}
    for run in runs:
        by_t.setdefault(run["T"], []).append(run)
    points = [mean_group(by_t[t]) for t in sorted(by_t)]
    print(
        f"# skip_ns {args.skip_ns} block_ns {args.block_ns}  "
        "v is half the thickness slope, A/ns"
    )
    print(
        f"{'T':>8} {'seeds':>5} {'v':>10} {'se':>10} "
        f"{'block':>10} {'bse':>10} {'dNdt':>10} {'se':>10} {'t_end':>8}"
    )
    for point in points:
        print(
            f"{point['T']:8.1f} {point['nseed']:5d} "
            f"{point['v']:10.4f} {point['ve']:10.4f} "
            f"{point['block_v']:10.4f} {point['block_se']:10.4f} "
            f"{point['dn']:10.2f} {point['dne']:10.2f} "
            f"{point['t_end']:8.2f}"
        )
    if len(points) < 2:
        return
    cold, warm = points[0], points[-1]
    t0, t0e = crossing(cold, warm, "v", "ve")
    tb, tbe = crossing(cold, warm, "block_v", "block_se")
    cold_grows = grows(cold)
    warm_melts = melts(warm)
    print(f"# ols_zero_K {t0:.2f} {t0e:.2f}")
    print(f"# block_zero_K {tb:.2f} {tbe:.2f}")
    print(f"# cold_grows {str(cold_grows).lower()} warm_melts {str(warm_melts).lower()}")
    passed = cold_grows and warm_melts
    print(f"# gate {'pass' if passed else 'fail'}")
    if not passed:
        sys.exit(2)


if __name__ == "__main__":
    main()
