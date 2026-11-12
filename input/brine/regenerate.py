#!/usr/bin/env python3
"""Build NaCl brine frames with TIP4P/2005 and Madrid-2019 scaled charges.

Writes data files and LAMMPS dumps under input/brine/work/. Those frames
are not committed. `--check` only tests the ice Ih cell.

    python3 input/brine/regenerate.py --check
    python3 input/brine/regenerate.py --lmp /path/to/lmp
"""

from __future__ import annotations

import argparse
import math
import os
import random
import shutil
import subprocess
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parent
DECK = ROOT / "in.tip4p2005-madrid2019.lmp"
WORK = ROOT / "work"

# Orthogonal ice Ih, 16 oxygens, from GenIce's 1h cell (absolute angstrom).
# Every oxygen has four neighbours between 2.76 and 2.77 A.
ICE_BOX = (7.84813412606925, 7.37735062301457, 9.06573834219084)
ICE_OXYGEN = (
    (1.328, 1.802, 3.38),
    (5.267, 4.524, 1.109),
    (6.58, 5.442, 3.365),
    (5.267, 4.542, 5.629),
    (2.623, 0.877, 5.644),
    (2.667, 5.488, 5.625),
    (5.241, 1.756, 1.12),
    (5.241, 1.774, 5.64),
    (1.354, 4.588, 7.888),
    (1.354, 4.57, 3.369),
    (2.667, 5.47, 1.105),
    (2.623, 0.858, 1.124),
    (6.537, 0.831, 3.384),
    (6.537, 0.849, 7.903),
    (6.581, 5.461, 7.884),
    (1.328, 1.82, 7.899),
)

OH = 0.9572
HOH = math.radians(104.52)
Q_O = -1.1128
Q_H = 0.5564
Q_NA = 0.85
Q_CL = -0.85
MOLAL_FACTOR = 18.01528 / 1000.0
RHO = 0.997 * 0.602214


def mic(p, q, box):
    acc = 0.0
    for k in range(3):
        d = q[k] - p[k]
        d -= box[k] * round(d / box[k])
        acc += d * d
    return math.sqrt(acc)


def check_ice():
    bad = 0
    for i, p in enumerate(ICE_OXYGEN):
        n = sum(1 for j, q in enumerate(ICE_OXYGEN) if i != j and mic(p, q, ICE_BOX) < 3.0)
        if n != 4:
            bad += 1
    if bad or len(ICE_OXYGEN) != 16:
        raise SystemExit(f"ice Ih cell is not 4-coordinated ({bad} sites)")


def replicate(nx, ny, nz):
    pts = []
    for ix in range(nx):
        for iy in range(ny):
            for iz in range(nz):
                shift = (ix * ICE_BOX[0], iy * ICE_BOX[1], iz * ICE_BOX[2])
                for p in ICE_OXYGEN:
                    pts.append((p[0] + shift[0], p[1] + shift[1], p[2] + shift[2]))
    box = (nx * ICE_BOX[0], ny * ICE_BOX[1], nz * ICE_BOX[2])
    return pts, box


def rotation(rng):
    u1, u2, u3 = rng.random(), rng.random(), rng.random()
    q = (
        math.sqrt(1 - u1) * math.sin(2 * math.pi * u2),
        math.sqrt(1 - u1) * math.cos(2 * math.pi * u2),
        math.sqrt(u1) * math.sin(2 * math.pi * u3),
        math.sqrt(u1) * math.cos(2 * math.pi * u3),
    )
    x, y, z, w = q
    return (
        (1 - 2 * (y * y + z * z), 2 * (x * y - z * w), 2 * (x * z + y * w)),
        (2 * (x * y + z * w), 1 - 2 * (x * x + z * z), 2 * (y * z - x * w)),
        (2 * (x * z - y * w), 2 * (y * z + x * w), 1 - 2 * (x * x + y * y)),
    )


def apply_rot(rot, v):
    return tuple(sum(rot[i][k] * v[k] for k in range(3)) for i in range(3))


def add(a, b):
    return (a[0] + b[0], a[1] + b[1], a[2] + b[2])


def tip4p_at(origin, rng):
    h1 = (OH, 0.0, 0.0)
    h2 = (OH * math.cos(HOH), OH * math.sin(HOH), 0.0)
    rot = rotation(rng)
    return origin, add(origin, apply_rot(rot, h1)), add(origin, apply_rot(rot, h2))


def ion_count(n_water, molality):
    return max(1, int(round(molality * n_water * MOLAL_FACTOR)))


def liquid_origins(n_water, rng):
    volume = n_water * 18.01528 / RHO
    side = volume ** (1.0 / 3.0)
    nside = math.ceil(n_water ** (1.0 / 3.0))
    spacing = side / nside
    origins = []
    for ix in range(nside):
        for iy in range(nside):
            for iz in range(nside):
                origins.append(
                    ((ix + 0.5) * spacing, (iy + 0.5) * spacing, (iz + 0.5) * spacing)
                )
                if len(origins) == n_water:
                    return origins, (side, side, side)
    raise RuntimeError("grid shorter than the water count")


def place_ions(n_pairs, box, oxygens, rng, zmin=None, gap=2.2):
    """Na and Cl at least `gap` angstrom from every oxygen. A cell list
    keeps the search on the neighbouring bins."""
    cell = gap
    dims = [max(1, math.ceil(box[k] / cell)) for k in range(3)]

    def key(pos):
        return tuple(min(dims[k] - 1, max(0, int(pos[k] / cell))) for k in range(3))

    buckets = {}
    for pos in oxygens:
        buckets.setdefault(key(pos), []).append(pos)

    def free(pos):
        center = key(pos)
        for ix in range(center[0] - 1, center[0] + 2):
            for iy in range(center[1] - 1, center[1] + 2):
                for iz in range(center[2] - 1, center[2] + 2):
                    for other in buckets.get((ix % dims[0], iy % dims[1], iz % dims[2]), ()):
                        if mic(pos, other, box) <= gap:
                            return False
        return True

    ions = []
    zlo = 0.0 if zmin is None else min(zmin, box[2] * 0.95)
    span = max(box[2] - zlo, box[2] * 0.05)
    for _kind in range(2):
        made = 0
        tries = 0
        limit = max(2000, n_pairs * 800)
        while made < n_pairs and tries < limit:
            tries += 1
            pos = (
                rng.random() * box[0],
                rng.random() * box[1],
                zlo + rng.random() * span,
            )
            if free(pos):
                ions.append(pos)
                buckets.setdefault(key(pos), []).append(pos)
                made += 1
        if made < n_pairs:
            raise RuntimeError(f"placed {made} of {n_pairs} ions")
    return ions[:n_pairs], ions[n_pairs:]


def build_bulk(n_water, molality, seed):
    rng = random.Random(seed)
    origins, box = liquid_origins(n_water, rng)
    waters = [tip4p_at(o, rng) for o in origins]
    n_pairs = ion_count(n_water, molality)
    # Ions need a gap. Stretch the cell a little so 2 mol/kg still fits.
    scale = (1.0 + 2.0 * n_pairs / n_water) ** (1.0 / 3.0)
    box = tuple(length * scale for length in box)
    waters = [tuple(tuple(coord * scale for coord in atom) for atom in trip) for trip in waters]
    na, cl = place_ions(n_pairs, box, [trip[0] for trip in waters], rng)
    return waters, na, cl, box, 0


def build_slab(n_target, molality, seed, reject=False):
    rng = random.Random(seed)
    best = None
    for nx in range(2, 30):
        for ny in range(2, 30):
            for nz in (2, 4, 6, 8, 10):
                n = 16 * nx * ny * nz
                if best is None or abs(n - n_target) < abs(best[0] - n_target):
                    best = (n, nx, ny, nz)
                if n == n_target:
                    break
    _n, nx, ny, nz = best
    origins, box = replicate(nx, ny, nz)
    mid = box[2] * 0.5
    ice_o = [p for p in origins if p[2] < mid]
    liq_o = [p for p in origins if p[2] >= mid]
    waters = [tip4p_at(o, rng) for o in ice_o + liq_o]
    n_pairs = ion_count(len(waters), molality)
    zmin = (mid + box[2]) * 0.5 if reject else None
    na, cl = place_ions(n_pairs, box, [trip[0] for trip in waters], rng, zmin=zmin)
    return waters, na, cl, box, len(ice_o)


def write_data(path, waters, na, cl, box):
    n_water = len(waters)
    n_atom = n_water * 3 + len(na) + len(cl)
    n_bond = n_water * 2
    n_angle = n_water
    lines = [
        "LAMMPS data file for TIP4P/2005 NaCl brine\n\n",
        f"{n_atom} atoms\n",
        f"{n_bond} bonds\n",
        f"{n_angle} angles\n\n",
        "4 atom types\n1 bond types\n1 angle types\n\n",
        f"0.0 {box[0]:.6f} xlo xhi\n",
        f"0.0 {box[1]:.6f} ylo yhi\n",
        f"0.0 {box[2]:.6f} zlo zhi\n\n",
        "Masses\n\n1 15.9994\n2 1.008\n3 22.9898\n4 35.453\n\n",
        "Atoms # full\n\n",
    ]
    atom_id = 1
    mol = 1
    bonds = []
    angles = []
    for ox, h1, h2 in waters:
        o_id, h1_id, h2_id = atom_id, atom_id + 1, atom_id + 2
        lines.append(f"{o_id} {mol} 1 {Q_O:.4f} {ox[0]:.5f} {ox[1]:.5f} {ox[2]:.5f}\n")
        lines.append(f"{h1_id} {mol} 2 {Q_H:.4f} {h1[0]:.5f} {h1[1]:.5f} {h1[2]:.5f}\n")
        lines.append(f"{h2_id} {mol} 2 {Q_H:.4f} {h2[0]:.5f} {h2[1]:.5f} {h2[2]:.5f}\n")
        bonds.append((o_id, h1_id))
        bonds.append((o_id, h2_id))
        angles.append((h1_id, o_id, h2_id))
        atom_id += 3
        mol += 1
    for pos in na:
        lines.append(f"{atom_id} {mol} 3 {Q_NA:.4f} {pos[0]:.5f} {pos[1]:.5f} {pos[2]:.5f}\n")
        atom_id += 1
        mol += 1
    for pos in cl:
        lines.append(f"{atom_id} {mol} 4 {Q_CL:.4f} {pos[0]:.5f} {pos[1]:.5f} {pos[2]:.5f}\n")
        atom_id += 1
        mol += 1
    lines.append("\nBonds\n\n")
    for i, (a, b) in enumerate(bonds, start=1):
        lines.append(f"{i} 1 {a} {b}\n")
    lines.append("\nAngles\n\n")
    for i, (a, b, c) in enumerate(angles, start=1):
        lines.append(f"{i} 1 {a} {b} {c}\n")
    path.write_text("".join(lines))


def write_input(path, freeze_n):
    text = DECK.read_text()
    if freeze_n > 0:
        text = text.replace(
            "velocity all create ${temp} 12345 dist gaussian\n",
            "group ice molecule 1:{n}\n"
            "group mobile subtract all ice\n"
            "fix freeze ice setforce 0.0 0.0 0.0\n"
            "velocity mobile create ${{temp}} 12345 dist gaussian\n".format(n=freeze_n),
        )
    path.write_text(text)


def find_lmp(explicit):
    if explicit:
        return explicit
    env = os.environ.get("LAMMPS_BIN")
    if env:
        return env
    found = shutil.which("lmp")
    if found:
        return found
    eessi = (
        "/cvmfs/software.eessi.io/versions/2023.06/software/linux/x86_64/"
        "intel/sapphirerapids/software/LAMMPS/2Aug2023_update2-foss-2023a-kokkos/bin/lmp"
    )
    if Path(eessi).is_file():
        return eessi
    raise SystemExit("lmp not found; pass --lmp")


def run_one(lmp, name, waters, na, cl, box, freeze_n, nsteps):
    system = WORK / name
    system.mkdir(parents=True, exist_ok=True)
    data = system / "data.brine"
    write_data(data, waters, na, cl, box)
    script = system / "in.brine"
    write_input(script, freeze_n)
    env = os.environ.copy()
    env["OMP_NUM_THREADS"] = "1"
    cmd = [
        lmp,
        "-in",
        str(script),
        "-var",
        "datafile",
        "data.brine",
        "-var",
        "nsteps",
        str(nsteps),
    ]
    log = system / "lammps.log"
    with log.open("w") as handle:
        proc = subprocess.run(cmd, cwd=system, env=env, stdout=handle, stderr=subprocess.STDOUT)
    if proc.returncode != 0:
        tail = log.read_text(errors="replace")[-2000:]
        raise SystemExit(f"LAMMPS failed for {name}\n{tail}")
    dump = system / "dump.lammpstrj"
    if not dump.is_file() or dump.stat().st_size == 0:
        raise SystemExit(f"{name} did not write dump.lammpstrj")
    print(f"{name} waters={len(waters)} ions={len(na) + len(cl)} freeze={freeze_n}")


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--check", action="store_true")
    parser.add_argument("--lmp", default="")
    parser.add_argument("--nsteps", type=int, default=0)
    parser.add_argument("--sizes", default="1000,8000,50000")
    args = parser.parse_args()
    check_ice()
    if args.check:
        print("ice Ih cell ok: 16 oxygens, four neighbours each")
        return
    lmp = find_lmp(args.lmp)
    sizes = [int(part) for part in args.sizes.split(",") if part]
    WORK.mkdir(parents=True, exist_ok=True)
    for n_water in sizes:
        for molality in (0.5, 1.0, 2.0):
            tag = f"nacl-n{n_water}-m{molality:.1f}"
            waters, na, cl, box, freeze_n = build_bulk(n_water, molality, seed=n_water * 10 + int(molality * 10))
            run_one(lmp, tag, waters, na, cl, box, freeze_n, args.nsteps)
        waters, na, cl, box, freeze_n = build_slab(n_water, 1.0, seed=1000 + n_water, reject=False)
        run_one(lmp, f"slab-n{n_water}", waters, na, cl, box, freeze_n, args.nsteps)
        waters, na, cl, box, freeze_n = build_slab(n_water, 1.0, seed=2000 + n_water, reject=True)
        run_one(lmp, f"front-n{n_water}", waters, na, cl, box, freeze_n, args.nsteps)


if __name__ == "__main__":
    main()
