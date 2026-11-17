# Copyright (c) d-SEAMS developers

"""mW growth-speed curve on the three ice faces.

The cell, the 500 K melt, the 2 ps Nose-Hoover periods and the 10 fs
step are the ones in espinosa_faces.py. That deck is called as it is.
240 K and 260 K keep the production lengths written there. The other
temperatures are added here, long enough for the slower basal face.

Basal, primary prism and secondary prism, at 220, 230, 240, 250, 260
and 265 K, four velocity seeds each. Frames every 5 ps go to si4-data/
(or $SI4_DATA) and are not committed.

    python3 input/interface/growth_curve.py
"""

from __future__ import annotations

import os
import sys
import time
import traceback
from concurrent.futures import ThreadPoolExecutor
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))

import espinosa_faces as faces  # noqa: E402

# Step counts at 10 fs. 60000 is 600 ps. 240 K and 260 K are left at
# the values in espinosa_faces.NSTEPS. 220 K and 230 K sit near the
# Table II maximum, so they use the 240 K length plus a margin. 250 K
# and 265 K are slower on the basal face, so the production is longer.
EXTRA_STEPS = {220: 80000, 230: 80000, 250: 120000, 265: 200000}
FACES = ("prism2", "prism1", "basal")
TEMPERATURES = (240, 260, 250, 230, 220, 265)
SEEDS = (1, 2, 3, 4)


def planned():
    for temperature in TEMPERATURES:
        for face in FACES:
            for seed in SEEDS:
                yield (face, temperature, seed, 10.0)


def main() -> None:
    for temperature, steps in EXTRA_STEPS.items():
        faces.NSTEPS[temperature] = steps
    missing = [temperature for temperature in TEMPERATURES if temperature not in faces.NSTEPS]
    if missing:
        raise SystemExit(f"no production length for {missing}")
    out = faces.OUT
    out.mkdir(parents=True, exist_ok=True)
    binary = faces.proto.lammps()
    cores = os.cpu_count() or 1
    work = list(planned())
    print(f"cores {cores} jobs {len(work)} data {out}", flush=True)
    for temperature in TEMPERATURES:
        print(f"steps T={temperature} {faces.NSTEPS[temperature]}", flush=True)
    failed = []
    t0 = time.perf_counter()

    def run(job):
        face, temperature, seed, step = job
        try:
            faces.one(job, binary)
            return job, None
        except Exception:
            return job, traceback.format_exc(limit=8)

    with ThreadPoolExecutor(max_workers=cores) as pool:
        for job, err in pool.map(run, work):
            if err:
                failed.append((job, err))
                print(f"fail {job[0]} T={job[1]} s={job[2]}\n{err}", flush=True)
    wall = time.perf_counter() - t0
    present = []
    absent = []
    for job in work:
        name = faces.tag(*job)
        dest = out / f"{name}.lammpstrj"
        if dest.exists() and dest.stat().st_size > 10000:
            present.append(name)
        else:
            absent.append(name)
    lines = [
        f"cores {cores}",
        f"wall_s {wall:.1f}",
        f"present {len(present)}",
        f"absent {len(absent)}",
        f"failed {len(failed)}",
    ]
    lines.extend(f"have {name}" for name in present)
    lines.extend(f"missing {name}" for name in absent)
    (out / "growth-curve-status.txt").write_text("\n".join(lines) + "\n")
    print(
        f"campaign {wall:.1f}s present {len(present)} absent {len(absent)} failed {len(failed)}",
        flush=True,
    )
    if failed:
        raise SystemExit(1)


if __name__ == "__main__":
    main()
