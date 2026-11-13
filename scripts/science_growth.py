#!/usr/bin/env python3
# Copyright (c) d-SEAMS developers

"""Collect basal front speeds at 230, 240 and 250 K.

science-front writes slope in angstrom per picosecond. This table stores
angstrom per nanosecond so it can sit next to the published TIP4P/Ice rates
in growth-literature.csv. The literature file is not rewritten.
"""

from __future__ import annotations

from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
DATA = ROOT / "paper" / "tex" / "figures" / "data" / "science"
TEMPS = (230, 240, 250)


def read_pairs(path: Path) -> dict[str, str]:
    rows = {}
    for line in path.read_text().splitlines():
        if not line or line.startswith("key ") or line.startswith("quantity "):
            continue
        parts = line.split()
        rows[parts[0]] = " ".join(parts[1:])
    return rows


def main() -> int:
    lines = [
        "temperature_K prod_ps dump_ps mean_front_A_per_ns lo_A_per_ns hi_A_per_ns "
        "thickness_A_per_ns n_frames"
    ]
    for temp in TEMPS:
        folder = DATA / f"interface-{temp}K"
        velocity = read_pairs(folder / "interface-velocity.csv")
        meta_path = DATA / f"interface-meta-{temp}K.csv"
        meta = read_pairs(meta_path) if meta_path.is_file() else {}
        mean = velocity["mean_front"].split()
        thick = velocity["thickness"].split()
        # mean_front row is "slope lo hi n" after the name was used as the key.
        slope, lo, hi, nframes = (float(mean[0]), float(mean[1]), float(mean[2]), mean[3])
        lines.append(
            f"{temp} {meta.get('prod_ps', '')} {meta.get('dump_ps', '')} "
            f"{1000 * slope:.6f} {1000 * lo:.6f} {1000 * hi:.6f} "
            f"{1000 * float(thick[0]):.6f} {nframes}"
        )
    out = DATA / "interface-speed.csv"
    out.write_text("\n".join(lines) + "\n")
    print(out.relative_to(ROOT).as_posix())
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
