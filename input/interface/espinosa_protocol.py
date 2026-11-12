# Copyright (c) d-SEAMS developers

"""Secondary-prism mW growth at the size in Espinosa et al.

J. Chem. Phys. 145, 211922 (2016), page 211922-2: 4000 molecules,
about half of them ice, box about 3 x 3 x 9 nm, anisotropic
barostat. Page 211922-4: mW was run in LAMMPS with a setup similar
to the TIP4P/ICE runs, whose thermostat and barostat each had a
relaxation time of about 2 ps. This deck uses Nose-Hoover NpT with
those periods. LAMMPS real units have no velocity-rescaling fix of
the Gromacs kind. The time step is 10 fs. The 2 fs step in that
paper is stated for the Gromacs TIP4P/ICE runs.

The oxygen cell is the GenIce ice Ih cell. Growth is along z, the
secondary-prism normal (planes 2.25 A apart). Five by five by ten
cells are 4000 molecules. The box is about 39 x 37 x 90 A.

Trajectories go to si4-data/ (or $SI4_DATA) and are not committed.

    python3 input/interface/espinosa_protocol.py
"""

from __future__ import annotations

import os
import shutil
import subprocess
import time
from concurrent.futures import ThreadPoolExecutor
from pathlib import Path

ROOT = Path(__file__).resolve().parents[2]
ICE = ROOT / "input" / "traj" / "ice1h.lammpstrj"
SW = ROOT / "input" / "mw215" / "mW.sw"
OUT = Path(os.environ.get("SI4_DATA", "si4-data"))
TEMPERATURES = (240, 260)
SEEDS = (1, 2, 3, 4)
# 1 bar, in atmospheres (LAMMPS real).
PRESSURE_ATM = 1.0 / 1.01325
# 10 fs steps. 240 K freezes in a few hundred picoseconds on this
# face. 260 K is given twice that.
NSTEPS = {240: 80000, 260: 160000}
DUMP_EVERY = 500  # 5 ps


def lammps() -> str:
    for name in ("lmp", "lmp_serial", "lammps"):
        found = shutil.which(name)
        if found:
            return found
    raise SystemExit("LAMMPS was not found")


def motif() -> tuple[list[tuple[float, float, float]], tuple[float, float, float]]:
    lines = ICE.read_text().splitlines()
    box: list[float] = []
    atoms: list[tuple[float, float, float]] = []
    mode = None
    for line in lines:
        if line.startswith("ITEM: BOX"):
            mode = "box"
            continue
        if line.startswith("ITEM: ATOMS"):
            mode = "atoms"
            continue
        if line.startswith("ITEM:"):
            mode = None
            continue
        if mode == "box" and len(box) < 3:
            lo, hi = line.split()[:2]
            box.append(float(hi) - float(lo))
        elif mode == "atoms":
            parts = line.split()
            if len(parts) >= 5:
                atoms.append((float(parts[2]), float(parts[3]), float(parts[4])))
    cell = tuple(length / 2.0 for length in box)
    kept = [
        (x, y, z)
        for x, y, z in atoms
        if x < cell[0] - 1.0e-6 and y < cell[1] - 1.0e-6 and z < cell[2] - 1.0e-6
    ]
    if len(kept) != 16:
        raise SystemExit(f"expected 16 ice sites, found {len(kept)}")
    return kept, cell  # type: ignore[return-value]


def slab() -> tuple[list[tuple[float, float, float]], tuple[float, float, float]]:
    sites, cell = motif()
    counts = (5, 5, 10)
    atoms = []
    for ix in range(counts[0]):
        for iy in range(counts[1]):
            for iz in range(counts[2]):
                shift = (ix * cell[0], iy * cell[1], iz * cell[2])
                for site in sites:
                    atoms.append(tuple(site[k] + shift[k] for k in range(3)))
    lengths = tuple(counts[k] * cell[k] for k in range(3))
    if len(atoms) != 4000:
        raise SystemExit(f"expected 4000 molecules, found {len(atoms)}")
    return atoms, lengths  # type: ignore[return-value]


def write_data(path: Path, atoms, lengths) -> None:
    lines = [
        "mW secondary prism, 4000 molecules",
        "",
        f"{len(atoms)} atoms",
        "1 atom types",
        "",
        f"0.0 {lengths[0]:.8f} xlo xhi",
        f"0.0 {lengths[1]:.8f} ylo yhi",
        f"0.0 {lengths[2]:.8f} zlo zhi",
        "",
        "Masses",
        "",
        "1 18.015",
        "",
        "Atoms",
        "",
    ]
    for i, (x, y, z) in enumerate(atoms, start=1):
        lines.append(f"{i} 1 {x:.6f} {y:.6f} {z:.6f}")
    path.write_text("\n".join(lines) + "\n")


def deck(temperature: int, seed: int, lengths) -> str:
    lo = 0.25 * lengths[2]
    hi = 0.75 * lengths[2]
    velocity = temperature * 10 + seed
    return f"""# secondary prism, {temperature} K, replica {seed}
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
timestep 10.0
region crystal block INF INF INF INF {lo:.6f} {hi:.6f}
group ice region crystal
group liquid subtract all ice
velocity liquid create 500.0 {velocity} mom yes rot yes dist gaussian
velocity ice set 0.0 0.0 0.0
fix hold ice setforce 0.0 0.0 0.0
# Melt with a short thermostat period so the liquid actually leaves
# the crystal wells. The 2 ps period is used for the growth run.
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
dump 1 all custom {DUMP_EVERY} growth.lammpstrj id type x y z
dump_modify 1 sort id first yes
thermo 2000
run {NSTEPS[temperature]}
"""


def one(temperature: int, seed: int, binary: str, atoms, lengths) -> None:
    dest = OUT / f"prism2-t{temperature}-s{seed}.lammpstrj"
    if dest.exists() and dest.stat().st_size > 10000:
        print(f"keep {dest.name}", flush=True)
        return
    work = OUT / f"work-t{temperature}-s{seed}"
    if work.exists():
        shutil.rmtree(work)
    work.mkdir(parents=True)
    write_data(work / "slab.data", atoms, lengths)
    shutil.copy(SW, work / "mW.sw")
    (work / "in.growth").write_text(deck(temperature, seed, lengths))
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
    # keep the log, it is small and records the thermostat
    shutil.copy(work / "log.lammps", OUT / f"prism2-t{temperature}-s{seed}.log")
    shutil.rmtree(work)
    print(f"wrote {dest.name} {time.perf_counter() - t0:.1f}s", flush=True)


def main() -> None:
    OUT.mkdir(parents=True, exist_ok=True)
    atoms, lengths = slab()
    print(
        f"n {len(atoms)} box {lengths[0]:.2f} {lengths[1]:.2f} {lengths[2]:.2f} A",
        flush=True,
    )
    binary = lammps()
    cores = os.cpu_count() or 1
    jobs = [(t, s) for t in TEMPERATURES for s in SEEDS]
    print(f"cores {cores} jobs {len(jobs)}", flush=True)
    t0 = time.perf_counter()
    with ThreadPoolExecutor(max_workers=cores) as pool:
        list(pool.map(lambda job: one(job[0], job[1], binary, atoms, lengths), jobs))
    print(f"campaign {time.perf_counter() - t0:.1f}s", flush=True)


if __name__ == "__main__":
    main()
