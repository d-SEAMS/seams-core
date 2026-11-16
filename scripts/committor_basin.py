#!/usr/bin/env python3
"""Committor table from the foot of each replica's final strict climb.

A replica that ends on 2N, and never passes it, was cut by the old wall.
It is marked truncated and left out of the committed fraction. The
barrier is the size at the foot of the climb, not 2N.

The series live in data/committor-series.tsv. The table is
data/committor-table.csv.
"""

from __future__ import annotations

import argparse
import csv
import subprocess
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]


def final_ascent(series: list[float]) -> tuple[bool, float]:
    if not series:
        return False, float("nan")
    imax = max(range(len(series)), key=lambda i: series[i])
    foot = imax
    while foot > 0 and series[foot - 1] < series[foot]:
        foot -= 1
    if not series[imax] > series[foot]:
        return False, float("nan")
    return True, series[foot]


def committor(series: list[float], molecules: int) -> dict[str, str]:
    climbed, barrier = final_ascent(series)
    maximum = max(series) if series else float("nan")
    final = series[-1] if series else float("nan")
    cap = 2.0 * molecules
    cut = bool(series) and final == cap and maximum == cap
    if cut:
        fate = "truncated"
    elif climbed and final >= barrier:
        fate = "committed"
    else:
        fate = "uncommitted"
    return {
        "barrier": "" if barrier != barrier else f"{barrier:.6f}",
        "final_nmax": f"{final:.6f}",
        "fate": fate,
        "cut_at_2n": "1" if cut else "0",
    }


def load_series(path: Path) -> list[tuple[str, int, list[float]]]:
    groups: dict[str, tuple[int, list[float]]] = {}
    with path.open() as handle:
        header = handle.readline().rstrip("\n").split("\t")
        idx = {name: i for i, name in enumerate(header)}
        for line in handle:
            if not line.strip():
                continue
            parts = line.rstrip("\n").split("\t")
            replica = parts[idx["replica"]]
            molecules = int(parts[idx["molecules"]])
            nmax = float(parts[idx["nmax"]])
            if replica not in groups:
                groups[replica] = (molecules, [])
            groups[replica][1].append(nmax)
    return [(name, groups[name][0], groups[name][1]) for name in groups]


def head_commit() -> str:
    done = subprocess.run(
        ["git", "rev-parse", "HEAD"], cwd=ROOT, text=True, capture_output=True, check=True
    )
    return done.stdout.strip()


def self_test() -> None:
    grown = committor([2, 2, 8, 40], 100)
    if grown["fate"] != "committed" or grown["barrier"] != "2.000000" or grown["cut_at_2n"] != "0":
        raise SystemExit(f"grown {grown}")
    past = committor([10, 20, 100, 250, 400], 100)
    if past["fate"] != "committed" or float(past["barrier"]) >= 200:
        raise SystemExit(f"past {past}")
    cut = committor([10, 40, 80, 200], 100)
    if cut["fate"] != "truncated" or cut["cut_at_2n"] != "1":
        raise SystemExit(f"cut {cut}")
    flat = committor([4, 4, 4], 100)
    if flat["fate"] != "uncommitted":
        raise SystemExit(f"flat {flat}")
    print("self-test ok")


def main() -> None:
    if "--self-test" in sys.argv:
        self_test()
        return
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--series", type=Path, required=True)
    parser.add_argument("--out", type=Path, required=True)
    parser.add_argument("--commit", default="")
    args = parser.parse_args()
    commit = args.commit or head_commit()
    rows = []
    used = 0
    committed = 0
    for replica, molecules, series in load_series(args.series):
        call = committor(series, molecules)
        if call["fate"] != "truncated":
            used += 1
            committed += call["fate"] == "committed"
        rows.append({
            "replica": replica,
            "molecules": str(molecules),
            **call,
            "commit": commit,
        })
    fraction = committed / used if used else float("nan")
    fields = ["replica", "molecules", "barrier", "final_nmax", "fate", "cut_at_2n", "commit"]
    args.out.parent.mkdir(parents=True, exist_ok=True)
    with args.out.open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields, lineterminator="\n")
        writer.writeheader()
        writer.writerows(rows)
    summary = args.out.with_name("committor-summary.csv")
    with summary.open("w", newline="") as handle:
        writer = csv.DictWriter(
            handle,
            fieldnames=["n_replicas", "n_used", "n_committed", "p_commit", "commit"],
            lineterminator="\n",
        )
        writer.writeheader()
        writer.writerow({
            "n_replicas": str(len(rows)),
            "n_used": str(used),
            "n_committed": str(committed),
            "p_commit": f"{fraction:.6f}",
            "commit": commit,
        })
    print(f"p_commit {fraction:.6f} used {used} truncated {len(rows) - used}")
    print(f"wrote {args.out}")
    print(f"wrote {summary}")


if __name__ == "__main__":
    main()
