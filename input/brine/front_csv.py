#!/usr/bin/env python3
"""Write a small per-frame CSV from a coexistence dump.

Frames stay where they are. The CSV is the series of ice thickness
and crystal count, plus one summary line with the slope and its
standard error along that single trajectory.
"""

from __future__ import annotations

import argparse
import math
import subprocess
from pathlib import Path


def parse_frames(text: str):
    rows = []
    for line in text.splitlines():
        if not line.startswith("frame "):
            continue
        parts = line.split()
        kv = dict(zip(parts[0::2], parts[1::2]))
        rows.append(
            (
                float(kv["t_ns"]),
                float(kv["n_crystal"]),
                float(kv["thickness"]),
                float(kv["c_liquid"]),
                float(kv["k"]),
                int(kv["n_ice_ion"]),
                int(kv["n_front_ion"]),
                int(kv["n_liquid_ion"]),
            )
        )
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


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("dump", type=Path)
    parser.add_argument("--brine-front", required=True)
    parser.add_argument("--out", type=Path, required=True)
    parser.add_argument("--skip-ns", type=float, default=0.06)
    parser.add_argument("--dt-ps", default="0.002")
    parser.add_argument("--every", default="10000")
    args = parser.parse_args()
    text = subprocess.check_output(
        [
            args.brine_front,
            str(args.dump),
            "--dt-ps",
            args.dt_ps,
            "--every",
            args.every,
            "--bins",
            "40",
        ],
        text=True,
    )
    rows = parse_frames(text)
    kept = [row for row in rows if row[0] >= args.skip_ns]
    thick_slope, thick_err = slope([r[0] for r in kept], [r[2] for r in kept])
    v, ve = thick_slope / 2.0, thick_err / 2.0
    dn, dne = slope([r[0] for r in kept], [r[1] for r in kept])
    args.out.parent.mkdir(parents=True, exist_ok=True)
    lines = [
        "t_ns,n_crystal,thickness_A,c_liquid,k,n_ice_ion,n_front_ion,n_liquid_ion"
    ]
    for row in rows:
        lines.append(
            f"{row[0]:.4f},{row[1]:.0f},{row[2]:.6f},{row[3]:.6f},{row[4]:.6f},{row[5]},{row[6]},{row[7]}"
        )
    lines.append(
        f"# skip_ns {args.skip_ns} n {len(kept)} "
        f"front_A_per_ns {v:.6f} {ve:.6f} dN_per_ns {dn:.4f} {dne:.4f}"
    )
    args.out.write_text("\n".join(lines) + "\n")
    print(args.out)
    print(lines[-1])


if __name__ == "__main__":
    main()
