# Copyright (c) d-SEAMS developers

"""mW slab growth at 240, 250, 260 and 270 K, four velocity seeds each.

The starting coordinates are input/traj/interface_f0.lammpstrj.
The ensemble is NVT at that volume, the mW Stillinger-Weber potential,
a 10 fs step, and one OpenMP thread per run. Frames are written every
10 ps for 800 ps. Trajectories stay under data/interface/campaign and
are not committed.

    python3 input/interface/campaign.py
"""

from __future__ import annotations

import os
import shutil
import subprocess
import time
from pathlib import Path

ROOT = Path(__file__).resolve().parents[2]
DECK = Path(__file__).resolve().parent
OUT = ROOT / "data" / "interface" / "campaign"
TEMPERATURES = (240, 250, 260, 270)
SEEDS = (1, 2, 3, 4)
NSTEPS = 80000  # 800 ps at 10 fs
DUMP_EVERY = 1000  # 10 ps


def lammps() -> str:
    for name in ("lmp", "lmp_serial", "lammps"):
        found = shutil.which(name)
        if found:
            return found
    raise SystemExit("LAMMPS was not found")


def deck(temperature: int, seed: int) -> str:
    # The velocity seed is temperature * 10 + replica, so the sixteen
    # runs do not share a random stream.
    velocity = temperature * 10 + seed
    return f"""# mW slab, NVT, {temperature} K, replica {seed}.
units real
atom_style atomic
boundary p p p
atom_modify map array
region box block 0 1 0 1 0 1
create_box 1 box
read_dump interface_f0.lammpstrj 0 x y z box yes add yes
mass 1 18.015
pair_style sw
pair_coeff * * mW.sw mW
neighbor 2.0 bin
neigh_modify every 1 delay 0 check yes
timestep 10.0
thermo 10000
velocity all create {temperature}.0 {velocity} mom yes rot yes dist gaussian
fix 1 all nvt temp {temperature}.0 {temperature}.0 100.0
dump 1 all custom {DUMP_EVERY} growth.lammpstrj id type x y z
dump_modify 1 sort id first yes
run {NSTEPS}
"""


def one(temperature: int, seed: int, binary: str) -> None:
    dest = OUT / f"t{temperature}-s{seed}.lammpstrj"
    if dest.exists() and dest.stat().st_size > 1000:
        print(f"keep {dest.name}", flush=True)
        return
    work = OUT / f"work-t{temperature}-s{seed}"
    if work.exists():
        shutil.rmtree(work)
    work.mkdir(parents=True)
    shutil.copy(ROOT / "input" / "mw215" / "mW.sw", work / "mW.sw")
    shutil.copy(ROOT / "input" / "traj" / "interface_f0.lammpstrj",
                work / "interface_f0.lammpstrj")
    (work / "in.growth").write_text(deck(temperature, seed))
    env = os.environ.copy()
    env["OMP_NUM_THREADS"] = "1"
    t0 = time.perf_counter()
    subprocess.run(
        [binary, "-in", "in.growth", "-log", "log.lammps"],
        cwd=work,
        env=env,
        check=True,
        stdout=subprocess.DEVNULL,
    )
    shutil.copy(work / "growth.lammpstrj", dest)
    shutil.rmtree(work)
    print(f"wrote {dest.name} {time.perf_counter() - t0:.1f}s", flush=True)


def main() -> None:
    OUT.mkdir(parents=True, exist_ok=True)
    binary = lammps()
    jobs = [(t, s) for t in TEMPERATURES for s in SEEDS]
    # One process per core. Each LAMMPS run is pinned to one thread.
    cores = os.cpu_count() or 1
    print(f"cores {cores} jobs {len(jobs)} lammps {binary}", flush=True)
    t0 = time.perf_counter()
    import concurrent.futures
    with concurrent.futures.ThreadPoolExecutor(max_workers=cores) as pool:
        list(pool.map(lambda job: one(job[0], job[1], binary), jobs))
    print(f"campaign {time.perf_counter() - t0:.1f}s", flush=True)


if __name__ == "__main__":
    main()
