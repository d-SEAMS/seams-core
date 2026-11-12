#!/usr/bin/env python3
"""Rebuild the 200 ps mW nucleation seeds.

A 12x12x12 ice Ic cell (13824 molecules) is melted for 40 ps at 500 K.
The 300 molecules nearest the centre are put back, liquid within 2.4 A of
that core is dropped, and the core is held fixed while the rest runs at
220 K for 200 ps. Velocity seeds are 91001, 92002 and 93003. LAMMPS is
the 22 Jul 2025 binary, pair_style sw, timestep 5 fs. Frames are not
committed.
"""

from __future__ import annotations

import math
import os
import subprocess
import sys

LMP = "/tmp/h2h-watch-venv/lib/python3.13/site-packages/lammps/lmp"
D = 2.75
MW_SW = """\
# DATE: 2009-03-01 UNITS: real CONTRIBUTOR: Molinero and Moore, J. Phys. Chem. B 113, 4008 (2009)
mW mW mW 6.189 2.3925 1.80 23.15 1.20 -0.333333333333
         7.049556277 0.6022245584 4.0 0.0 0.0
"""


def ice_ic(nrep):
    lattice = 4.0 * D / math.sqrt(3.0)
    basis = (
        (0, 0, 0), (0.5, 0.5, 0), (0.5, 0, 0.5), (0, 0.5, 0.5),
        (0.25, 0.25, 0.25), (0.75, 0.75, 0.25), (0.75, 0.25, 0.75), (0.25, 0.75, 0.75),
    )
    side = nrep * lattice
    atoms = []
    for ix in range(nrep):
        for iy in range(nrep):
            for iz in range(nrep):
                for bx, by, bz in basis:
                    atoms.append(((ix + bx) * lattice, (iy + by) * lattice, (iz + bz) * lattice))
    return side, atoms


def core_indices(atoms, side, target=300):
    center = side * 0.5
    scored = []
    for i, pos in enumerate(atoms):
        d2 = sum((pos[k] - center) ** 2 for k in range(3))
        scored.append((d2, i))
    scored.sort()
    return [i for _, i in scored[:target]]


def write_data(path, side, atoms):
    with open(path, "w", encoding="utf-8") as out:
        out.write("mW nucleation\n\n")
        out.write(f"{len(atoms)} atoms\n1 atom types\n\n")
        out.write(f"0.0 {side} xlo xhi\n0.0 {side} ylo yhi\n0.0 {side} zlo zhi\n\n")
        out.write("Masses\n\n1 18.01528\n\nAtoms # atomic\n\n")
        for i, pos in enumerate(atoms, start=1):
            out.write(f"{i} 1 {pos[0]} {pos[1]} {pos[2]}\n")


def run(workdir, text):
    os.makedirs(workdir, exist_ok=True)
    path = os.path.join(workdir, "in.mw")
    with open(path, "w", encoding="utf-8") as out:
        out.write(text)
    env = os.environ.copy()
    env["LD_LIBRARY_PATH"] = "/tmp/h2h-watch-venv/lib"
    env["OMP_NUM_THREADS"] = "1"
    log = os.path.join(workdir, "log.lammps")
    cmd = [
        "/tmp/h2h-watch-venv/bin/mpirun", "-np", "8",
        LMP, "-in", "in.mw",
    ]
    with open(log, "w", encoding="utf-8") as err:
        subprocess.run(cmd, cwd=workdir, check=True, stdout=err, stderr=subprocess.STDOUT, env=env)


def read_dump(path):
    chunks = open(path, encoding="utf-8").read().split("ITEM: TIMESTEP")
    body = chunks[-1].splitlines()
    atoms = []
    reading = False
    for line in body:
        if line.startswith("ITEM: ATOMS"):
            reading = True
            continue
        if reading and line.strip():
            parts = line.split()
            atoms.append((int(parts[0]), float(parts[2]), float(parts[3]), float(parts[4])))
    atoms.sort()
    return [(xyz[1], xyz[2], xyz[3]) for xyz in atoms]


def main():
    outdir = sys.argv[1]
    steps = int(sys.argv[2]) if len(sys.argv) > 2 else 40000
    os.makedirs(outdir, exist_ok=True)
    side, crystal = ice_ic(12)
    core = core_indices(crystal, side, 300)
    print(f"side={side:.4f} n={len(crystal)} core={len(core)}", flush=True)
    work = os.path.join(outdir, "melt")
    os.makedirs(work, exist_ok=True)
    with open(os.path.join(work, "mW.sw"), "w", encoding="utf-8") as out:
        out.write(MW_SW)
    write_data(os.path.join(work, "data.mw"), side, crystal)
    run(work, "\n".join([
        "units real", "atom_style atomic", "boundary p p p",
        "read_data data.mw",
        "pair_style sw", "pair_coeff * * mW.sw mW",
        "neighbor 2.0 bin", "neigh_modify every 1 delay 0 check yes",
        "timestep 5.0", "thermo 2000", "thermo_style custom step temp pe",
        "velocity all create 500.0 50011 dist gaussian",
        "fix bath all nvt temp 500.0 500.0 200.0",
        f"run {8000 if steps > 1000 else steps}",
        "write_dump all custom melt.lammpstrj id type xu yu zu",
        "",
    ]))
    if steps <= 1000:
        print("probe melt done", flush=True)
        return
    melted = read_dump(os.path.join(work, "melt.lammpstrj"))
    core_pos = [crystal[i] for i in core]

    def overlaps(pos):
        for core_xyz in core_pos:
            d2 = 0.0
            for a, b in zip(pos, core_xyz):
                d = a - b
                d -= side * math.floor(d / side + 0.5)
                d2 += d * d
            if d2 < 2.4 * 2.4:
                return True
        return False

    mixed = []
    core_ids = []
    for i, pos in enumerate(crystal):
        if i in set(core):
            core_ids.append(len(mixed))
            mixed.append(pos)
        else:
            wrapped = tuple(x - side * math.floor(x / side) for x in melted[i])
            if overlaps(wrapped):
                continue
            mixed.append(wrapped)
    print(f"mixed={len(mixed)} core={len(core_ids)}", flush=True)
    for seed in (91001, 92002, 93003):
        seed_dir = os.path.join(outdir, f"seed-{seed}")
        os.makedirs(seed_dir, exist_ok=True)
        with open(os.path.join(seed_dir, "mW.sw"), "w", encoding="utf-8") as out:
            out.write(MW_SW)
        write_data(os.path.join(seed_dir, "data-core.mw"), side, mixed)
        ids = " ".join(str(i + 1) for i in core_ids)
        run(seed_dir, "\n".join([
            "units real", "atom_style atomic", "boundary p p p",
            "read_data data-core.mw",
            "pair_style sw", "pair_coeff * * mW.sw mW",
            "neighbor 2.0 bin", "neigh_modify every 1 delay 0 check yes",
            "timestep 5.0", "thermo 4000", "thermo_style custom step temp pe",
            f"group core id {ids}",
            "group outer subtract all core",
            "fix hold core setforce 0.0 0.0 0.0",
            "min_style cg",
            "minimize 1.0e-4 1.0e-6 40 400",
            f"velocity outer create 220.0 {seed} dist gaussian",
            "fix bath outer nvt temp 220.0 220.0 250.0",
            f"run {steps}",
            "write_dump all custom last.lammpstrj id type xu yu zu",
            "",
        ]))
        coords = read_dump(os.path.join(seed_dir, "last.lammpstrj"))
        frame = os.path.join(outdir, f"seed-{seed}-200ps.xyz")
        with open(frame, "w", encoding="utf-8") as out:
            out.write(f"{len(coords)}\nside {side}\n")
            for x, y, z in coords:
                out.write(f"O {x} {y} {z}\n")
        print(f"wrote {frame} n={len(coords)}", flush=True)


if __name__ == "__main__":
    main()
