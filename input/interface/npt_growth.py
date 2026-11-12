# Copyright (c) d-SEAMS developers

"""Direct-coexistence mW growth on the basal and secondary-prism faces.

The oxygen sites are the GenIce ice Ih cell in input/traj/ice1h.lammpstrj.
That cell has the basal normal along y (bilayers 0.92 A apart, then a
2.76 A gap) and the secondary-prism normal along z (planes 2.25 A
apart). Each face is an orthorhombic slab of 8064 molecules, ice in
the middle half and liquid in the two caps, so the periodic cell has
two interfaces. The run is anisotropic NpT at 1 bar. Trajectories are
written under data/interface/npt and are not committed.

    python3 input/interface/npt_growth.py
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
OUT = ROOT / "data" / "interface" / "npt"
TEMPERATURES = (240, 250, 260, 270)
SEEDS = (1, 2, 3, 4)
FACES = ("basal", "prism2")
# 1 bar in atmospheres, which is the pressure unit of LAMMPS real.
PRESSURE_ATM = 1.0 / 1.01325
# 10 fs steps. 240 and 250 K freeze inside 600 ps. 260 and 270 K
# are slower, so those runs are longer.
NSTEPS = {240: 60000, 250: 80000, 260: 140000, 270: 200000}
DUMP_EVERY = 1000  # 10 ps


def lammps() -> str:
    for name in ("lmp", "lmp_serial", "lammps"):
        found = shutil.which(name)
        if found:
            return found
    raise SystemExit("LAMMPS was not found")


def read_ice() -> tuple[list[tuple[float, float, float]], tuple[float, float, float]]:
    lines = ICE.read_text().splitlines()
    box = []
    atoms = []
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
            if len(parts) < 5:
                continue
            atoms.append((float(parts[2]), float(parts[3]), float(parts[4])))
    return atoms, (box[0], box[1], box[2])


def motif() -> tuple[list[tuple[float, float, float]], tuple[float, float, float]]:
    atoms, box = read_ice()
    cell = tuple(length / 2.0 for length in box)
    kept = []
    for x, y, z in atoms:
        if x < cell[0] - 1.0e-6 and y < cell[1] - 1.0e-6 and z < cell[2] - 1.0e-6:
            kept.append((x, y, z))
    if len(kept) != 16:
        raise SystemExit(f"expected 16 sites in the ice cell, found {len(kept)}")
    return kept, cell  # type: ignore[return-value]


def slab(face: str) -> tuple[list[tuple[float, float, float]], tuple[float, float, float], int]:
    sites, cell = motif()
    # Six cells in each lateral direction and fourteen along the normal
    # gives 16 * 6 * 6 * 14 = 8064 molecules.
    if face == "basal":
        counts = (6, 14, 6)
        axis = 1
    elif face == "prism2":
        counts = (6, 6, 14)
        axis = 2
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
    return atoms, lengths, axis


def write_data(path: Path, atoms: list[tuple[float, float, float]], lengths: tuple[float, float, float]) -> None:
    lines = [
        "mW ice slab",
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


def deck(face: str, temperature: int, seed: int, axis: int, lengths: tuple[float, float, float]) -> str:
    lo = 0.25 * lengths[axis]
    hi = 0.75 * lengths[axis]
    if axis == 1:
        block = f"INF INF {lo:.6f} {hi:.6f} INF INF"
    else:
        block = f"INF INF INF INF {lo:.6f} {hi:.6f}"
    velocity = temperature * 10 + seed + (0 if face == "basal" else 100)
    return f"""# {face} face, {temperature} K, replica {seed}, NpT at 1 bar.
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
region crystal block {block}
group ice region crystal
group liquid subtract all ice
velocity liquid create 400.0 {velocity} mom yes rot yes dist gaussian
velocity ice set 0.0 0.0 0.0
fix hold ice setforce 0.0 0.0 0.0
fix melt liquid nvt temp 400.0 320.0 100.0
run 5000
unfix melt
fix cool liquid nvt temp {temperature}.0 {temperature}.0 200.0
run 2000
unfix cool
unfix hold
reset_timestep 0
velocity all create {temperature}.0 {velocity} mom yes rot yes dist gaussian
fix 1 all npt temp {temperature}.0 {temperature}.0 200.0 &
    x {PRESSURE_ATM:.8f} {PRESSURE_ATM:.8f} 2000.0 &
    y {PRESSURE_ATM:.8f} {PRESSURE_ATM:.8f} 2000.0 &
    z {PRESSURE_ATM:.8f} {PRESSURE_ATM:.8f} 2000.0
dump 1 all custom {DUMP_EVERY} growth.lammpstrj id type x y z
dump_modify 1 sort id first yes
thermo 5000
run {NSTEPS[temperature]}
"""


def one(face: str, temperature: int, seed: int, binary: str) -> None:
    dest = OUT / f"{face}-t{temperature}-s{seed}.lammpstrj"
    if dest.exists() and dest.stat().st_size > 10000:
        print(f"keep {dest.name}", flush=True)
        return
    atoms, lengths, axis = slab(face)
    work = OUT / f"work-{face}-t{temperature}-s{seed}"
    if work.exists():
        shutil.rmtree(work)
    work.mkdir(parents=True)
    write_data(work / "slab.data", atoms, lengths)
    shutil.copy(SW, work / "mW.sw")
    (work / "in.growth").write_text(deck(face, temperature, seed, axis, lengths))
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
    print(f"wrote {dest.name} {time.perf_counter() - t0:.1f}s n {len(atoms)}", flush=True)


def main() -> None:
    OUT.mkdir(parents=True, exist_ok=True)
    atoms, lengths, _axis = slab("basal")
    print(f"basal n {len(atoms)} box {lengths[0]:.2f} {lengths[1]:.2f} {lengths[2]:.2f}")
    atoms, lengths, _axis = slab("prism2")
    print(f"prism2 n {len(atoms)} box {lengths[0]:.2f} {lengths[1]:.2f} {lengths[2]:.2f}")
    binary = lammps()
    cores = os.cpu_count() or 1
    jobs = [(face, t, s) for face in FACES for t in TEMPERATURES for s in SEEDS]
    print(f"cores {cores} jobs {len(jobs)} pressure_atm {PRESSURE_ATM:.6f}", flush=True)
    t0 = time.perf_counter()
    with ThreadPoolExecutor(max_workers=cores) as pool:
        list(pool.map(lambda job: one(job[0], job[1], job[2], binary), jobs))
    print(f"campaign {time.perf_counter() - t0:.1f}s", flush=True)


if __name__ == "__main__":
    main()
