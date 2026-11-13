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


def stacking(nx=4, ny=4, sequence=SEQUENCE):
    n_bi = len(sequence)
    height = n_bi * 0.5
    vecs = (
        (nx * A_HEX, 0.0, 0.0),
        (-0.5 * ny * A_HEX, ny * A_HEX * math.sqrt(3.0) / 2.0, 0.0),
        (0.0, 0.0, height * C_HEX),
    )
    atoms = []
    for ix in range(nx):
        for iy in range(ny):
            for layer, site in enumerate(sequence):
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


def nucleus_frame(outdir):
    """Melt the cell, put an ice Ic core back, and let the interface sit for 1 ps."""
    vecs, atoms = ice_ic()
    work = os.path.join(outdir, "run-melt")
    os.makedirs(work, exist_ok=True)
    write_sw(os.path.join(work, "mW.sw"))
    write_data(os.path.join(work, "data.mw"), vecs, atoms)
    melt = "\n".join(
        [
            "units real",
            "atom_style atomic",
            "boundary p p p",
            "read_data data.mw",
            "pair_style sw",
            "pair_coeff * * mW.sw mW",
            "neighbor 2.0 bin",
            "neigh_modify every 1 delay 0 check yes",
            "timestep 5.0",
            "thermo 2000",
            "thermo_style custom step temp pe",
            "velocity all create 500.0 50011 dist gaussian",
            "fix bath all nvt temp 500.0 500.0 200.0",
            "run 8000",
            "write_dump all custom melt.lammpstrj id type xu yu zu",
            "",
        ]
    )
    run_lammps(work, melt)
    melted = read_dump(os.path.join(work, "melt.lammpstrj"))
    side = vecs[0][0]
    core_pos = []
    for atom in atoms:
        if atom[2] == 1:
            core_pos.append(atom[0])

    def wrap(pos):
        return tuple(x - side * math.floor(x / side) for x in pos)

    def overlaps(pos):
        for core in core_pos:
            d2 = 0.0
            for a, b in zip(pos, core):
                d = a - b
                d -= side * math.floor(d / side + 0.5)
                d2 += d * d
            if d2 < 2.4 * 2.4:
                return True
        return False

    mixed = []
    for i, atom in enumerate(atoms):
        if atom[2] == 1:
            mixed.append((atom[0], -1, 1))
        else:
            pos = wrap(melted[i])
            if overlaps(pos):
                continue
            mixed.append((pos, -1, 0))
    core_ids = [i for i, atom in enumerate(mixed) if atom[2] == 1]
    write_data(os.path.join(work, "data-core.mw"), vecs, mixed)
    relax = "\n".join(
        [
            "units real",
            "atom_style atomic",
            "boundary p p p",
            "read_data data-core.mw",
            "pair_style sw",
            "pair_coeff * * mW.sw mW",
            "neighbor 2.0 bin",
            "neigh_modify every 1 delay 0 check yes",
            "timestep 5.0",
            "thermo 200",
            "thermo_style custom step temp pe",
            "group core id " + " ".join(str(i + 1) for i in core_ids),
            "group outer subtract all core",
            "fix hold core setforce 0.0 0.0 0.0",
            "min_style cg",
            "minimize 1.0e-4 1.0e-6 40 400",
            "velocity outer create 250.0 25011 dist gaussian",
            "fix bath outer nvt temp 250.0 250.0 200.0",
            "run 400",
            "unfix hold",
            "unfix bath",
            "write_dump all custom last.lammpstrj id type xu yu zu",
            "",
        ]
    )
    run_lammps(work, relax)
    coords = read_dump(os.path.join(work, "last.lammpstrj"))
    moved = [(coords[i], mixed[i][1], mixed[i][2]) for i in range(len(mixed))]
    write_frame(os.path.join(outdir, "nucleus-250.frame"), "nucleus-250", 250.0, "nucleus", vecs, moved)
    print(f"wrote nucleus atoms={len(mixed)} core={len(core_ids)}", flush=True)


def hydrate_sI(nrep=5):
    """Guest-free sI oxygen lattice. Sites are the standard cell in tests/hydrate_sites.inc."""
    here = os.path.dirname(os.path.abspath(__file__))
    inc = os.path.join(here, "..", "tests", "hydrate_sites.inc")
    text = open(inc, encoding="utf-8").read()
    body = text.split("kSIxyz")[1].split("kSII")[0]
    sites = []
    for line in body.splitlines():
        stripped = line.strip()
        if not stripped.startswith("{"):
            continue
        parts = [float(x) for x in stripped.strip("{},").split(",") if x.strip()]
        if len(parts) == 3:
            sites.append(tuple(parts))
    if len(sites) != 46:
        raise SystemExit(f"sI site table has {len(sites)} rows")
    cell = 12.238818320447
    side = nrep * cell
    vecs = ((side, 0.0, 0.0), (0.0, side, 0.0), (0.0, 0.0, side))
    atoms = []
    for ix in range(nrep):
        for iy in range(nrep):
            for iz in range(nrep):
                for sx, sy, sz in sites:
                    atoms.append(((ix * cell + sx, iy * cell + sy, iz * cell + sz), -1, 0))
    return vecs, atoms


def quench(workdir, vecs, atoms, steps=2000):
    """Steepest-descent inherent structure of a thermal frame."""
    os.makedirs(workdir, exist_ok=True)
    write_sw(os.path.join(workdir, "mW.sw"))
    write_data(os.path.join(workdir, "data.mw"), vecs, atoms)
    script = "\n".join(
        [
            "units real",
            "atom_style atomic",
            "boundary p p p",
            "read_data data.mw",
            "pair_style sw",
            "pair_coeff * * mW.sw mW",
            "neighbor 2.0 bin",
            "neigh_modify every 1 delay 0 check yes",
            "min_style sd",
            f"minimize 0.0 1.0e-4 {steps} {steps * 10}",
            "write_dump all custom quench.lammpstrj id type xu yu zu",
            "",
        ]
    )
    run_lammps(workdir, script)
    coords = read_dump(os.path.join(workdir, "quench.lammpstrj"))
    if len(coords) != len(atoms):
        raise SystemExit(f"quench atom count {len(coords)} != {len(atoms)}")
    return [(coords[i], atoms[i][1], atoms[i][2]) for i in range(len(atoms))]


def roc_frames(outdir):
    """Thermal ensemble from 200 K to 270 K and its steepest-descent quench."""
    systems = []
    systems.append(("ic", ice_ic(4)))
    systems.append(("ih", stacking(sequence="ABABABAB")))
    systems.append(("stack", stacking()))
    temps = (200.0, 220.0, 240.0, 260.0, 270.0)
    for name, (vecs, atoms) in systems:
        for temperature in temps:
            tag = f"{name}-{int(temperature)}"
            moved = thermal(
                os.path.join(outdir, f"run-{tag}"),
                tag,
                temperature,
                "crystal",
                vecs,
                atoms,
                2000,
                int(temperature) * 100 + len(name),
            )
            write_frame(
                os.path.join(outdir, f"{tag}.frame"),
                tag,
                temperature,
                "crystal",
                vecs,
                moved,
            )
            quenched = quench(os.path.join(outdir, f"run-{tag}-q"), vecs, moved)
            write_frame(
                os.path.join(outdir, f"{tag}-q.frame"),
                f"{tag}-q",
                0.0,
                "quench",
                vecs,
                quenched,
            )
            print(f"wrote {tag} and quench", flush=True)
    # The large cells already on disk, quenched the same way.
    for path in (
        "ice-4096-220.frame",
        "stack-220.frame",
        "stack-270.frame",
        "nucleus-250.frame",
    ):
        full = os.path.join(outdir, path)
        if not os.path.exists(full):
            continue
        lines = open(full, encoding="utf-8").read().splitlines()
        n = int(lines[3].split()[1])
        vecs = tuple(tuple(float(x) for x in lines[5 + i].split()) for i in range(3))
        atoms = []
        for line in lines[9 : 9 + n]:
            parts = line.split()
            atoms.append(
                (
                    (float(parts[0]), float(parts[1]), float(parts[2])),
                    int(parts[3]),
                    int(parts[4]),
                )
            )
        tag = path.replace(".frame", "")
        quenched = quench(os.path.join(outdir, f"run-{tag}-q"), vecs, atoms)
        write_frame(
            os.path.join(outdir, f"{tag}-q.frame"),
            f"{tag}-q",
            0.0,
            "quench",
            vecs,
            quenched,
        )
        print(f"quenched {tag}", flush=True)


def hydrate_sII(nrep=2):
    here = os.path.dirname(os.path.abspath(__file__))
    inc = os.path.join(here, "..", "tests", "hydrate_sites.inc")
    text = open(inc, encoding="utf-8").read()
    body = text.split("kSIIxyz")[1]
    sites = []
    for line in body.splitlines():
        stripped = line.strip()
        if not stripped.startswith("{"):
            continue
        parts = [float(x) for x in stripped.strip("{},").split(",") if x.strip()]
        if len(parts) == 3:
            sites.append(tuple(parts))
    if len(sites) != 136:
        raise SystemExit(f"sII site table has {len(sites)} rows")
    cell = 6.20577 * 2.80
    side = nrep * cell
    vecs = ((side, 0.0, 0.0), (0.0, side, 0.0), (0.0, 0.0, side))
    atoms = []
    for ix in range(nrep):
        for iy in range(nrep):
            for iz in range(nrep):
                for sx, sy, sz in sites:
                    atoms.append(((ix * cell + sx, iy * cell + sy, iz * cell + sz), -1, 0))
    return vecs, atoms


def cage_frames(outdir):
    """Thermal sI and sII, plus the steepest-descent inherent structure."""
    for name, build, nrep, seed in (
        ("sI", hydrate_sI, 2, 22071),
        ("sII", hydrate_sII, 2, 22081),
    ):
        vecs, atoms = build(nrep)
        for temperature in (220.0, 250.0):
            tag = f"{name}-{int(temperature)}"
            moved = thermal(
                os.path.join(outdir, f"run-{tag}"),
                tag,
                temperature,
                "hydrate",
                vecs,
                atoms,
                2000,
                seed + int(temperature),
            )
            write_frame(
                os.path.join(outdir, f"{tag}.frame"),
                tag,
                temperature,
                "hydrate",
                vecs,
                moved,
            )
            quenched = quench(os.path.join(outdir, f"run-{tag}-q"), vecs, moved)
            write_frame(
                os.path.join(outdir, f"{tag}-q.frame"),
                f"{tag}-q",
                0.0,
                "quench",
                vecs,
                quenched,
            )
            print(f"wrote {tag} atoms={len(moved)} and quench", flush=True)


def large_frames(outdir):
    """Thermal ice of at least 4096 molecules, and a thermal sI cell."""
    vecs, atoms = ice_ic(8)
    moved = thermal(
        os.path.join(outdir, "run-ice4096"),
        "ice-4096-220",
        220.0,
        "crystal",
        vecs,
        atoms,
        2000,
        409611,
    )
    write_frame(
        os.path.join(outdir, "ice-4096-220.frame"),
        "ice-4096-220",
        220.0,
        "crystal",
        vecs,
        moved,
    )
    print(f"wrote ice-4096-220 atoms={len(moved)}", flush=True)
    vecs, atoms = hydrate_sI(5)
    moved = thermal(
        os.path.join(outdir, "run-si5750"),
        "sI-5750-220",
        220.0,
        "hydrate",
        vecs,
        atoms,
        2000,
        575011,
    )
    write_frame(
        os.path.join(outdir, "sI-5750-220.frame"),
        "sI-5750-220",
        220.0,
        "hydrate",
        vecs,
        moved,
    )
    print(f"wrote sI-5750-220 atoms={len(moved)}", flush=True)


def main():
    outdir = sys.argv[1] if len(sys.argv) > 1 and not sys.argv[1].startswith("--") else "/tmp/sym-v3"
    if "--large" in sys.argv:
        os.makedirs(outdir, exist_ok=True)
        large_frames(outdir)
        return
    if "--roc" in sys.argv:
        os.makedirs(outdir, exist_ok=True)
        roc_frames(outdir)
        return
    if "--cages" in sys.argv:
        os.makedirs(outdir, exist_ok=True)
        cage_frames(outdir)
        return
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
    nucleus_frame(outdir)


if __name__ == "__main__":
    main()
