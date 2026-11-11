#!/usr/bin/env python3
# Copyright (c) d-SEAMS developers
# SPDX-License-Identifier: MIT
"""Small mW frames for the symmetry kernel.

Stacking-disordered ice is eight bilayers on the hexagonal lattice, site
sequence ABABCACB, with A=(1/3,2/3), B=(2/3,1/3), C=(0,0). A step between
two different sites keeps the oxygen-oxygen bond at 2.75 A. The nucleus is
a 4x4x4 ice Ic cell whose molecules outside an 8 A sphere are heated, then
the whole cell is held at 250 K. Neither trajectory is written into the tree.

Usage: mw_frames.py outdir [--ideal-only]
"""

from __future__ import annotations

import math
import os
import subprocess
import sys

D = 2.75
C_HEX = D / 0.375
A_HEX = C_HEX / math.sqrt(8.0 / 3.0)
SITES = {
    "A": (1.0 / 3.0, 2.0 / 3.0),
    "B": (2.0 / 3.0, 1.0 / 3.0),
    "C": (0.0, 0.0),
}
SEQUENCE = "ABABCACB"
MW_SW = """\
# DATE: 2009-03-01 UNITS: real CONTRIBUTOR: Molinero and Moore, J. Phys. Chem. B 113, 4008 (2009)
# epsilon kcal/mol, sigma Angstrom
mW mW mW 6.189 2.3925 1.80 23.15 1.20 -0.333333333333
         7.049556277 0.6022245584 4.0 0.0 0.0
"""


def wrap(value):
    return value - math.floor(value)


def cart(frac, vecs):
    out = [0.0, 0.0, 0.0]
    for i in range(3):
        for k in range(3):
            out[k] += frac[i] * vecs[i][k]
    return out


def stacking(nx=4, ny=4):
    n_bi = len(SEQUENCE)
    height = n_bi * 0.5
    vecs = (
        (nx * A_HEX, 0.0, 0.0),
        (-0.5 * ny * A_HEX, ny * A_HEX * math.sqrt(3.0) / 2.0, 0.0),
        (0.0, 0.0, height * C_HEX),
    )
    atoms = []
    for ix in range(nx):
        for iy in range(ny):
            for layer, site in enumerate(SEQUENCE):
                sx, sy = SITES[site]
                for z_unit, which in ((1.0 / 16.0, 0), (7.0 / 16.0, 1)):
                    frac = (
                        wrap((sx + ix) / nx),
                        wrap((sy + iy) / ny),
                        wrap((layer * 0.5 + z_unit) / height),
                    )
                    atoms.append((cart(frac, vecs), layer, 0))
    return vecs, atoms


def ice_ic(nrep=4):
    lattice = 4.0 * D / math.sqrt(3.0)
    basis = (
        (0, 0, 0),
        (0.5, 0.5, 0),
        (0.5, 0, 0.5),
        (0, 0.5, 0.5),
        (0.25, 0.25, 0.25),
        (0.75, 0.75, 0.25),
        (0.75, 0.25, 0.75),
        (0.25, 0.75, 0.75),
    )
    side = nrep * lattice
    vecs = ((side, 0.0, 0.0), (0.0, side, 0.0), (0.0, 0.0, side))
    center = side * 0.5
    atoms = []
    for ix in range(nrep):
        for iy in range(nrep):
            for iz in range(nrep):
                for bx, by, bz in basis:
                    pos = ((ix + bx) * lattice, (iy + by) * lattice, (iz + bz) * lattice)
                    dx = pos[0] - center
                    dy = pos[1] - center
                    dz = pos[2] - center
                    core = 1 if dx * dx + dy * dy + dz * dz <= 64.0 else 0
                    atoms.append((pos, -1, core))
    return vecs, atoms


def write_frame(path, name, temperature, role, vecs, atoms):
    with open(path, "w", encoding="utf-8") as out:
        out.write(f"name {name}\n")
        out.write(f"temperature {temperature}\n")
        out.write(f"role {role}\n")
        out.write(f"n {len(atoms)}\nbox\n")
        for row in vecs:
            out.write(f"{row[0]} {row[1]} {row[2]}\n")
        out.write("atoms\n")
        for pos, layer, core in atoms:
            out.write(f"{pos[0]} {pos[1]} {pos[2]} {layer} {core}\n")


def write_data(path, vecs, atoms):
    a, b, c = vecs
    with open(path, "w", encoding="utf-8") as out:
        out.write("mW frame\n\n")
        out.write(f"{len(atoms)} atoms\n")
        out.write("1 atom types\n\n")
        out.write(f"0.0 {a[0]} xlo xhi\n")
        out.write(f"0.0 {b[1]} ylo yhi\n")
        out.write(f"0.0 {c[2]} zlo zhi\n")
        out.write(f"{b[0]} {c[0]} {c[1]} xy xz yz\n\n")
        out.write("Masses\n\n1 18.01528\n\n")
        out.write("Atoms # atomic\n\n")
        for i, (pos, _layer, _core) in enumerate(atoms, start=1):
            out.write(f"{i} 1 {pos[0]} {pos[1]} {pos[2]}\n")


def write_sw(path):
    with open(path, "w", encoding="utf-8") as out:
        out.write(MW_SW)


def lammps_input(data, dump, temperature, steps, seed, nucleus=False):
    lines = [
        "units real",
        "atom_style atomic",
        "boundary p p p",
        f"read_data {data}",
        "pair_style sw",
        "pair_coeff * * mW.sw mW",
        "neighbor 2.0 bin",
        "neigh_modify every 1 delay 0 check yes",
        "timestep 5.0",
        "thermo 1000",
        "thermo_style custom step temp pe",
        "min_style cg",
        "minimize 1.0e-6 1.0e-8 400 4000",
        "reset_timestep 0",
    ]
    if nucleus:
        lines += [
            "variable cx equal (xhi+xlo)/2",
            "variable cy equal (yhi+ylo)/2",
            "variable cz equal (zhi+zlo)/2",
            "region core sphere ${cx} ${cy} ${cz} 8",
            "group core region core",
            "group outer subtract all core",
            f"velocity outer create 360.0 {seed} dist gaussian",
            f"velocity core create 230.0 {seed + 1} dist gaussian",
            "fix hot outer nvt temp 360.0 360.0 250.0",
            "fix cold core nvt temp 230.0 230.0 250.0",
            "run 2000",
            "unfix hot",
            "unfix cold",
            "fix bath all nvt temp 250.0 250.0 250.0",
            "run 2000",
        ]
    else:
        lines += [
            f"velocity all create {temperature} {seed} dist gaussian",
            f"fix bath all nvt temp {temperature} {temperature} 250.0",
            f"run {steps}",
        ]
    lines += [f"write_dump all custom {dump} id type xu yu zu"]
    return "\n".join(lines) + "\n"


def run_lammps(workdir, script_text):
    path = os.path.join(workdir, "in.mw")
    with open(path, "w", encoding="utf-8") as out:
        out.write(script_text)
    log = os.path.join(workdir, "log.lammps")
    with open(log, "w", encoding="utf-8") as err:
        subprocess.run(["lmp", "-in", "in.mw"], cwd=workdir, check=True, stdout=err, stderr=subprocess.STDOUT)


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


def thermal(workdir, name, temperature, role, vecs, atoms, steps, seed, nucleus=False):
    os.makedirs(workdir, exist_ok=True)
    write_sw(os.path.join(workdir, "mW.sw"))
    write_data(os.path.join(workdir, "data.mw"), vecs, atoms)
    dump = "last.lammpstrj"
    run_lammps(
        workdir,
        lammps_input("data.mw", dump, temperature, steps, seed, nucleus=nucleus),
    )
    coords = read_dump(os.path.join(workdir, dump))
    if len(coords) != len(atoms):
        raise SystemExit(f"{name}: dump has {len(coords)} atoms, data has {len(atoms)}")
    moved = [(coords[i], atoms[i][1], atoms[i][2]) for i in range(len(atoms))]
    return moved


def main():
    outdir = sys.argv[1] if len(sys.argv) > 1 else "/tmp/sym-v3"
    ideal_only = "--ideal-only" in sys.argv
    os.makedirs(outdir, exist_ok=True)
    vecs, atoms = stacking()
    write_frame(os.path.join(outdir, "stack-ideal.frame"), "stack-ideal", 0, "crystal", vecs, atoms)
    print(f"ideal stacking atoms={len(atoms)}", flush=True)
    if ideal_only:
        return
    for temperature, seed in ((220.0, 22011), (270.0, 27011)):
        moved = thermal(
            os.path.join(outdir, f"run-{int(temperature)}"),
            f"stack-{int(temperature)}",
            temperature,
            "crystal",
            vecs,
            atoms,
            4000,
            seed,
        )
        write_frame(
            os.path.join(outdir, f"stack-{int(temperature)}.frame"),
            f"stack-{int(temperature)}",
            temperature,
            "crystal",
            vecs,
            moved,
        )
        print(f"wrote stack-{int(temperature)}", flush=True)
    nvecs, natoms = ice_ic()
    moved = thermal(
        os.path.join(outdir, "run-nucleus"),
        "nucleus-250",
        250.0,
        "nucleus",
        nvecs,
        natoms,
        0,
        25011,
        nucleus=True,
    )
    write_frame(
        os.path.join(outdir, "nucleus-250.frame"),
        "nucleus-250",
        250.0,
        "nucleus",
        nvecs,
        moved,
    )
    print(f"wrote nucleus atoms={len(natoms)}", flush=True)


if __name__ == "__main__":
    main()
