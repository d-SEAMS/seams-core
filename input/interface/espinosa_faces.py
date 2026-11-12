# Copyright (c) d-SEAMS developers

"""mW direct-coexistence runs for the Table I candidates.

Same 4000-molecule cell, melt, and 2 ps Nose-Hoover periods as
espinosa_protocol.py. The secondary prism is repeated at 10 fs.
The basal and primary prism faces use the same molecule count with
the long axis along that normal. One secondary-prism seed is also
run at a 5 fs step, with the dump interval kept at 5 ps.

Trajectories go to si4-data/ (or $SI4_DATA) and are not committed.

    python3 input/interface/espinosa_faces.py
"""

from __future__ import annotations

import os
import shutil
import subprocess
import time
from concurrent.futures import ThreadPoolExecutor
from pathlib import Path

import espinosa_protocol as proto

ROOT = proto.ROOT
SW = proto.SW
OUT = Path(os.environ.get("SI4_DATA", "si4-data"))
PRESSURE_ATM = proto.PRESSURE_ATM
# 10 fs steps. 240 K freezes in a few hundred picoseconds. 260 K
# and the basal face are given longer.
NSTEPS = {240: 60000, 260: 100000}
DUMP_EVERY = 500  # 5 ps at 10 fs


def slab(face: str):
    sites, cell = proto.motif()
    if face == "prism2":
        counts = (5, 5, 10)
        axis = 2
    elif face == "prism1":
        counts = (10, 5, 5)
        axis = 0
    elif face == "basal":
        counts = (5, 10, 5)
        axis = 1
    else:
        raise SystemExit(face)
    atoms = []
    for ix in range(counts[0]):
        for iy in range(counts[1]):
            for iz in range(counts[2]):
                shift = (ix * cell[0], iy * cell[1], iz * cell[2])
                for site in sites:
                    atoms.append(tuple(site[k] + shift[k] for k in range(3)))
    lengths = tuple(counts[k] * cell[k] for k in range(3))
    if len(atoms) != 4000:
        raise SystemExit(f"{face}: expected 4000 molecules, found {len(atoms)}")
    if lengths[axis] <= max(lengths[k] for k in range(3) if k != axis):
        raise SystemExit(f"{face}: growth axis is not the long axis")
    return atoms, lengths, axis


def deck(face: str, temperature: int, seed: int, axis: int, lengths, step_fs: float) -> str:
    lo = 0.25 * lengths[axis]
    hi = 0.75 * lengths[axis]
    if axis == 0:
        block = f"{lo:.6f} {hi:.6f} INF INF INF INF"
    elif axis == 1:
        block = f"INF INF {lo:.6f} {hi:.6f} INF INF"
    else:
        block = f"INF INF INF INF {lo:.6f} {hi:.6f}"
    velocity = temperature * 10 + seed + 100 * axis
    # Keep frames 5 ps apart whatever the step is.
    dump_every = max(1, int(round(5000.0 / step_fs)))
    steps = int(round(NSTEPS[temperature] * 10.0 / step_fs))
    return f"""# {face}, {temperature} K, replica {seed}, step {step_fs:g} fs
# Nose-Hoover periods 2 ps, anisotropic, 1 bar.
units real
atom_style atomic
boundary p p p
atom_modify map array
read_data slab.data
pair_style sw
pair_coeff * * mW.sw mW
neighbor 2.0 bin
neigh_modify every 1 delay 0 check yes
timestep {step_fs:.1f}
region crystal block {block}
group ice region crystal
group liquid subtract all ice
velocity liquid create 500.0 {velocity} mom yes rot yes dist gaussian
velocity ice set 0.0 0.0 0.0
fix hold ice setforce 0.0 0.0 0.0
fix melt liquid nvt temp 500.0 500.0 100.0
run 10000
unfix melt
fix cool liquid nvt temp {temperature}.0 {temperature}.0 2000.0
run 1000
unfix cool
unfix hold
reset_timestep 0
velocity all create {temperature}.0 {velocity} mom yes rot yes dist gaussian
fix 1 all npt temp {temperature}.0 {temperature}.0 2000.0 &
    x {PRESSURE_ATM:.8f} {PRESSURE_ATM:.8f} 2000.0 &
    y {PRESSURE_ATM:.8f} {PRESSURE_ATM:.8f} 2000.0 &
    z {PRESSURE_ATM:.8f} {PRESSURE_ATM:.8f} 2000.0
thermo_style custom step temp pe press lx ly lz
thermo {dump_every}
dump 1 all custom {dump_every} growth.lammpstrj id type x y z
dump_modify 1 sort id first yes
run {steps}
"""


def tag(face: str, temperature: int, seed: int, step_fs: float) -> str:
    if abs(step_fs - 10.0) < 1.0e-6:
        return f"{face}-t{temperature}-s{seed}"
    return f"{face}-dt{step_fs:.0f}-t{temperature}-s{seed}"


def one(job, binary: str) -> None:
    face, temperature, seed, step_fs = job
    name = tag(face, temperature, seed, step_fs)
    dest = OUT / f"{name}.lammpstrj"
    if dest.exists() and dest.stat().st_size > 10000:
        print(f"keep {dest.name}", flush=True)
        return
    atoms, lengths, axis = slab(face)
    work = OUT / f"work-{name}"
    if work.exists():
        shutil.rmtree(work)
    work.mkdir(parents=True)
    proto.write_data(work / "slab.data", atoms, lengths)
    shutil.copy(SW, work / "mW.sw")
    (work / "in.growth").write_text(deck(face, temperature, seed, axis, lengths, step_fs))
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
    shutil.copy(work / "log.lammps", OUT / f"{name}.log")
    shutil.rmtree(work)
    print(f"wrote {dest.name} {time.perf_counter() - t0:.1f}s", flush=True)


def jobs():
    out = []
    for temperature in (240, 260):
        for seed in (1, 2, 3, 4):
            out.append(("prism2", temperature, seed, 10.0))
        for seed in (1, 2):
            out.append(("prism1", temperature, seed, 10.0))
            out.append(("basal", temperature, seed, 10.0))
    out.append(("prism2", 240, 1, 5.0))
    return out


def main() -> None:
    OUT.mkdir(parents=True, exist_ok=True)
    for face in ("prism2", "prism1", "basal"):
        atoms, lengths, axis = slab(face)
        print(
            f"{face} n {len(atoms)} axis {axis} box {lengths[0]:.2f} {lengths[1]:.2f} {lengths[2]:.2f}",
            flush=True,
        )
    binary = proto.lammps()
    cores = os.cpu_count() or 1
    work = jobs()
    print(f"cores {cores} jobs {len(work)}", flush=True)
    t0 = time.perf_counter()
    with ThreadPoolExecutor(max_workers=cores) as pool:
        list(pool.map(lambda job: one(job, binary), work))
    print(f"campaign {time.perf_counter() - t0:.1f}s", flush=True)


if __name__ == "__main__":
    main()
