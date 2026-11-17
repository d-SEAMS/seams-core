#!/usr/bin/env python3
"""SI density-window table from the Ih and ice VI constants in phase.hpp.

Ice is below the midpoint of those two number densities. The old 0.028
ceiling is listed so the table shows it is not the split the code uses.
"""

from __future__ import annotations

import argparse
import csv
import re
import subprocess
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
HEADER = ROOT / "src" / "include" / "internal" / "phase.hpp"
NA = 6.02214076e23
MOLAR = 18.015


def constants() -> tuple[float, float]:
    text = HEADER.read_text()
    ih = float(re.search(r"kIceIhMassDensity = ([0-9.]+);", text).group(1))
    vi = float(re.search(r"kIceVIMassDensity = ([0-9.]+);", text).group(1))
    return ih, vi


def number_density(mass: float) -> float:
    return mass * NA / (MOLAR * 1.0e24)


def kind(rho: float, ice_max: float) -> str:
    if rho < ice_max:
        return "ice"
    if rho < 0.040:
        return "lda"
    if rho < 0.070:
        return "mda"
    return "hda"


def self_test() -> None:
    ih, vi = constants()
    ice_max = 0.5 * (number_density(ih) + number_density(vi))
    if not (number_density(ih) < ice_max < number_density(vi)):
        raise SystemExit("the midpoint does not separate Ih and ice VI")
    if kind(number_density(ih), ice_max) != "ice":
        raise SystemExit("ideal Ih is outside the ice window")
    if kind(number_density(vi), ice_max) == "ice":
        raise SystemExit("ice VI is inside the ice window")
    if abs(ice_max - 0.028) < 1e-6:
        raise SystemExit("the window is still the old 0.028 constant")
    print("self-test ok")


def head_commit() -> str:
    done = subprocess.run(
        ["git", "rev-parse", "HEAD"], cwd=ROOT, text=True, capture_output=True, check=True
    )
    return done.stdout.strip()


def main() -> None:
    if "--self-test" in sys.argv:
        self_test()
        return
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--out", type=Path, required=True)
    parser.add_argument("--commit", default="")
    args = parser.parse_args()
    ih, vi = constants()
    ih_rho = number_density(ih)
    vi_rho = number_density(vi)
    ice_max = 0.5 * (ih_rho + vi_rho)
    commit = args.commit or head_commit()
    rows = [
        ("ice Ih", f"{ih_rho:.8f}", kind(ih_rho, ice_max)),
        ("ice VI", f"{vi_rho:.8f}", kind(vi_rho, ice_max)),
        ("midpoint", f"{ice_max:.8f}", "boundary"),
        ("old ceiling 0.028", "0.02800000", kind(0.028, ice_max)),
        ("0.033", "0.03300000", kind(0.033, ice_max)),
        ("0.039", "0.03900000", kind(0.039, ice_max)),
    ]
    args.out.parent.mkdir(parents=True, exist_ok=True)
    with args.out.open("w", newline="") as handle:
        writer = csv.writer(handle, lineterminator="\n")
        writer.writerow(["phase", "rho_per_A3", "kind", "ice_max", "commit"])
        for phase, rho, label in rows:
            writer.writerow([phase, rho, label, f"{ice_max:.8f}", commit])
    print(f"ice_max {ice_max:.8f}")
    print(f"wrote {args.out}")


if __name__ == "__main__":
    main()
