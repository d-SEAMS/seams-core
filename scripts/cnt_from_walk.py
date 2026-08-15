#!/usr/bin/env python3
"""CNT-style read of walk_seeded output.

Niu et al. PRL 122, 245501 (2019) report Nc = 314 +- 20.
This script takes n_max(t) and reports barrier recrossings and
whether clusters above Nc tend to grow.
"""
from __future__ import annotations

import argparse
import sys


def main() -> int:
    p = argparse.ArgumentParser()
    p.add_argument("table")
    p.add_argument("--nc", type=int, default=314)
    p.add_argument("--window", type=int, default=20)
    args = p.parse_args()

    frames: list[int] = []
    nmax: list[int] = []
    nice: list[int] = []
    nclus: list[int] = []
    with open(args.table) as fh:
        for line in fh:
            if not line.strip() or line.startswith("#"):
                continue
            w = line.split()
            frames.append(int(w[0]))
            nice.append(int(w[6]))
            nmax.append(int(w[7]))
            nclus.append(int(w[8]))

    if len(frames) < 3:
        print("too few samples", file=sys.stderr)
        return 1

    nc = args.nc
    above = [n >= nc for n in nmax]
    first_up = next((frames[i] for i, a in enumerate(above) if a), None)
    last_down = None
    for i in range(1, len(above)):
        if above[i - 1] and not above[i]:
            last_down = frames[i]
    last_up = next(
        (frames[i] for i in range(len(above) - 1, -1, -1) if above[i]), None
    )
    recross = 0
    for i in range(1, len(above)):
        if above[i] != above[i - 1]:
            recross += 1

    grow = 0
    shrink = 0
    lo, hi = nc - args.window, nc + args.window
    for i in range(len(nmax) - 1):
        if lo <= nmax[i] <= hi:
            if nmax[i + 1] > nmax[i]:
                grow += 1
            elif nmax[i + 1] < nmax[i]:
                shrink += 1

    print(f"samples {len(frames)} frames {frames[0]}..{frames[-1]}")
    print(f"n_max range {min(nmax)}..{max(nmax)}  n_ice range {min(nice)}..{max(nice)}")
    print(f"n_clus range {min(nclus)}..{max(nclus)}")
    print(f"Nc {nc}")
    print(f"first n_max>=Nc frame {first_up}")
    print(f"last drop below Nc frame {last_down}")
    print(f"last n_max>=Nc frame {last_up}")
    print(f"crossings of Nc {recross}")
    print(f"near Nc [+/-{args.window}]: grow {grow} shrink {shrink}")
    if grow + shrink:
        print(f"P(grow|near Nc) {grow / (grow + shrink):.3f}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
