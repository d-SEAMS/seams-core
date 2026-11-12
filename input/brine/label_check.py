#!/usr/bin/env python3
"""Compare published CHILL+ crystal labels with freud's Steinhardt q6.

freud  does not implement CHILL+. On one frame the oxygen q6 values are
split by the d-SEAMS crystal label (hexagonal or cubic). Agreement of a
fixed threshold q6 > 0.45 with that label is reported as a count, not
as a second implementation of CHILL+.
"""

from __future__ import annotations

import argparse
import subprocess
from pathlib import Path

import freud
import numpy as np


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("dump")
    parser.add_argument("--brine-front", required=True)
    parser.add_argument("--q6", type=float, default=0.45)
    parser.add_argument("--cutoff", type=float, default=3.5)
    args = parser.parse_args()
    label_path = Path(args.dump).with_suffix(".chill.txt")
    subprocess.run(
        [
            args.brine_front,
            args.dump,
            "--max-frames",
            "1",
            "--labels",
            str(label_path),
            "--cutoff",
            str(args.cutoff),
        ],
        check=True,
        stdout=subprocess.DEVNULL,
    )
    rows = []
    for line in label_path.read_text().splitlines():
        if not line or line.startswith("#"):
            continue
        parts = line.split()
        rows.append((float(parts[2]), float(parts[3]), float(parts[4]), int(parts[6])))
    positions = np.array([[row[0], row[1], row[2]] for row in rows], dtype=float)
    crystal = np.array([row[3] == 1 for row in rows])
    # The box is not in the label file. Read it from the first dump frame.
    bounds = []
    capture = False
    for line in Path(args.dump).read_text().splitlines():
        if line.startswith("ITEM: BOX BOUNDS"):
            capture = True
            bounds = []
            continue
        if capture:
            if line.startswith("ITEM:"):
                break
            lo, hi = line.split()[:2]
            bounds.append(float(hi) - float(lo))
    box = np.diag(np.array(bounds, dtype=float))
    steinhardt = freud.order.Steinhardt(l=6)
    steinhardt.compute((box, positions), {"r_max": args.cutoff, "exclude_ii": True})
    order = np.asarray(steinhardt.particle_order, dtype=float)
    freud_crystal = order > args.q6
    both = int(np.count_nonzero(crystal & freud_crystal))
    chill_only = int(np.count_nonzero(crystal & ~freud_crystal))
    freud_only = int(np.count_nonzero(~crystal & freud_crystal))
    neither = int(np.count_nonzero(~crystal & ~freud_crystal))
    n = len(order)
    agree = both + neither

    def median(mask):
        values = order[mask]
        if values.size == 0:
            return float("nan")
        return float(np.median(values))

    print(f"oxygens {n}")
    print(f"chill_crystal {int(crystal.sum())} freud_crystal {int(freud_crystal.sum())}")
    print(f"q6_median_crystal {median(crystal):.4f} q6_median_other {median(~crystal):.4f}")
    print(f"threshold {args.q6} agree {agree} of {n} fraction {agree / n:.4f}")
    print(f"both {both} chill_only {chill_only} freud_only {freud_only} neither {neither}")


if __name__ == "__main__":
    main()
