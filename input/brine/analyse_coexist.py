#!/usr/bin/env python3
"""Block means of a direct-coexistence log and, if a dump exists, of CHILL+.

The LAMMPS timestep is 2 fs. Production starts at the last reset, which
this script takes to be the thermo table whose first step is 0 and whose
temperature is near the target, after the melt at 300 K.
"""

from __future__ import annotations

import argparse
import math
import statistics
import subprocess
from pathlib import Path


def thermo_tables(text):
    tables = []
    current = None
    for line in text.splitlines():
        parts = line.split()
        if len(parts) >= 5 and parts[0] == "Step" and parts[1] == "Temp":
            if current:
                tables.append(current)
            current = []
            continue
        if current is None or len(parts) < 4 or not parts[0].isdigit():
            continue
        try:
            current.append(
                {
                    "step": int(parts[0]),
                    "temp": float(parts[1]),
                    "pe": float(parts[2]),
                    "density": float(parts[3]),
                }
            )
        except ValueError:
            continue
    if current:
        tables.append(current)
    return tables


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


def blocks(values, n_blocks):
    if len(values) < n_blocks:
        n_blocks = len(values)
    if n_blocks < 2:
        return float("nan"), float("nan")
    width = len(values) // n_blocks
    means = [statistics.fmean(values[i * width : (i + 1) * width]) for i in range(n_blocks)]
    mean = statistics.fmean(means)
    var = sum((m - mean) ** 2 for m in means) / (len(means) - 1)
    return mean, math.sqrt(var / len(means))


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("work", type=Path)
    parser.add_argument("--skip-ps", type=float, default=50.0)
    parser.add_argument("--brine-front", default="")
    args = parser.parse_args()
    log = (args.work / "lammps.log").read_text(errors="replace")
    tables = thermo_tables(log)
    # Melt is near 300 K. Production is the last table.
    prod = tables[-1] if tables else []
    rows = [row for row in prod if row["step"] * 0.002 >= args.skip_ps]
    print(f"production_rows {len(prod)} used {len(rows)}")
    if not rows:
        return
    t_ns = [row["step"] * 2e-6 for row in rows]
    pe = [row["pe"] for row in rows]
    b, err = slope(t_ns, pe)
    pe_mean, pe_se = blocks(pe, 6)
    print(f"pe_kcal {pe_mean:.2f} {pe_se:.2f}")
    print(f"pe_slope_kcal_per_ns {b:.2f} {err:.2f}")
    print(f"t_ns {t_ns[0]:.3f} {t_ns[-1]:.3f}")
    dump = args.work / "dump.lammpstrj"
    if args.brine_front and dump.is_file() and dump.stat().st_size > 0:
        report = args.work / "front.txt"
        subprocess.run(
            [
                args.brine_front,
                str(dump),
                "--dt-ps",
                "0.002",
                "--every",
                "10000",
                "--bins",
                "40",
                "--time",
                "0",
            ],
            check=True,
            stdout=report.open("w"),
        )
        summary = next(line for line in report.read_text().splitlines() if line.startswith("summary "))
        print(summary)


if __name__ == "__main__":
    main()
