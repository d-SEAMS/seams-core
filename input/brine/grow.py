#!/usr/bin/env python3
"""Madrid-2019 ice/brine trajectories. Frames stay in the run directory.

The deck is TIP4P/2005 plus the Madrid-2019 NaCl parameters
(Zeron, Abascal and Vega, J. Chem. Phys. 151, 134504 (2019),
doi:10.1063/1.5121392). Ion charges are +/- 0.85.

    python3 input/brine/grow.py --ice-check
    python3 input/brine/grow.py --crystal --lmp /path/to/lmp
    python3 input/brine/grow.py --grow --lmp /path/to/lmp --prod-ps 400

--crystal is the Table V check: rock-salt NaCl with the scaled charges
at 298.15 K and 1 bar, published density 2.050 g/cm^3.
"""

from __future__ import annotations

import argparse
import os
import math
import random
import shutil
import subprocess
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))
import regenerate as regen

OUT = Path(os.environ.get("BRINE_DIR", "brine-runs"))
KCAL_PER_KJ = 1.0 / 4.184
OH = regen.OH
HOH = regen.HOH


def mic_vec(p, q, box):
    return tuple(q[k] - p[k] - box[k] * round((q[k] - p[k]) / box[k]) for k in range(3))


def norm(v):
    length = math.sqrt(sum(c * c for c in v))
    if length < 1e-10:
        return (0.0, 0.0, 1.0)
    return tuple(c / length for c in v)


def cross(a, b):
    return (
        a[1] * b[2] - a[2] * b[1],
        a[2] * b[0] - a[0] * b[2],
        a[0] * b[1] - a[1] * b[0],
    )


def add(a, b):
    return (a[0] + b[0], a[1] + b[1], a[2] + b[2])


def wrap(p, box):
    return tuple(p[k] - box[k] * math.floor(p[k] / box[k]) for k in range(3))


def bond_graph(oxygens, box, cutoff=3.0):
    neigh = [[] for _ in oxygens]
    for i, origin in enumerate(oxygens):
        for j in range(i + 1, len(oxygens)):
            if regen.mic(origin, oxygens[j], box) < cutoff:
                neigh[i].append(j)
                neigh[j].append(i)
    return neigh


def key_of(i, j):
    return (i, j) if i < j else (j, i)


def assign_ice_rules(neigh, rng):
    n = len(neigh)
    bonds = []
    seen = set()
    for i, partners in enumerate(neigh):
        for j in partners:
            key = key_of(i, j)
            if key not in seen:
                seen.add(key)
                bonds.append(key)

    def hydrogen_count(owner, i):
        return sum(1 for j in neigh[i] if owner[key_of(i, j)] == i)

    for _attempt in range(30):
        owner = {bond: bond[0] if rng.random() < 0.5 else bond[1] for bond in bonds}

        def repair():
            from collections import deque

            for _step in range(n * 4):
                excess = [i for i in range(n) if hydrogen_count(owner, i) > 2]
                if not excess:
                    return True
                start = excess[0]
                previous = {start: None}
                queue = deque([start])
                target = None
                while queue:
                    current = queue.popleft()
                    if current != start and hydrogen_count(owner, current) < 2:
                        target = current
                        break
                    for other in neigh[current]:
                        if owner[key_of(current, other)] == current and other not in previous:
                            previous[other] = current
                            queue.append(other)
                if target is None:
                    return False
                node = target
                while previous[node] is not None:
                    parent = previous[node]
                    owner[key_of(parent, node)] = node
                    node = parent
            return all(hydrogen_count(owner, i) == 2 for i in range(n))

        if repair() and all(hydrogen_count(owner, i) == 2 for i in range(n)):
            return owner
    raise RuntimeError("ice rules did not converge")


def orient_water(origin, toward_a, toward_b):
    bisector = norm(add(norm(toward_a), norm(toward_b)))
    normal = norm(cross(toward_a, toward_b))
    side = cross(normal, bisector)
    h1 = (OH, 0.0, 0.0)
    h2 = (OH * math.cos(HOH), OH * math.sin(HOH), 0.0)
    local_bis = norm((h1[0] + h2[0], h1[1] + h2[1], 0.0))
    local_normal = (0.0, 0.0, 1.0)
    local_side = cross(local_normal, local_bis)

    def rotate(vector):
        c_bis = sum(vector[i] * local_bis[i] for i in range(3))
        c_side = sum(vector[i] * local_side[i] for i in range(3))
        c_normal = sum(vector[i] * local_normal[i] for i in range(3))
        return tuple(
            c_bis * bisector[i] + c_side * side[i] + c_normal * normal[i] for i in range(3)
        )

    return origin, add(origin, rotate(h1)), add(origin, rotate(h2))


def periodic_ice(nx, ny, nz, seed):
    rng = random.Random(seed)
    oxygens, box = regen.replicate(nx, ny, nz)
    neigh = bond_graph(oxygens, box)
    bad = [i for i, partners in enumerate(neigh) if len(partners) != 4]
    if bad:
        raise RuntimeError(f"{len(bad)} ice oxygens are not 4-coordinated")
    owner = assign_ice_rules(neigh, rng)
    waters = []
    for i, origin in enumerate(oxygens):
        acceptors = [j for j in neigh[i] if owner[key_of(i, j)] == i]
        if len(acceptors) != 2:
            raise RuntimeError("ice oxygen does not own two hydrogens")
        v1 = mic_vec(origin, oxygens[acceptors[0]], box)
        v2 = mic_vec(origin, oxygens[acceptors[1]], box)
        oxygen, h1, h2 = orient_water(origin, v1, v2)
        waters.append((wrap(oxygen, box), wrap(h1, box), wrap(h2, box)))
    return waters, box


def liquid_waters(n_water, box_xy, z0, z1, rng):
    height = z1 - z0
    if n_water <= 0 or height <= 0:
        return []
    spacing = (18.01528 / regen.RHO) ** (1.0 / 3.0)
    n_x = max(1, round(box_xy[0] / spacing))
    n_y = max(1, round(box_xy[1] / spacing))
    per_layer = n_x * n_y
    n_z = max(1, math.ceil(n_water / per_layer))
    dx = box_xy[0] / n_x
    dy = box_xy[1] / n_y
    dz = height / n_z
    origins = []
    for iz in range(n_z):
        for iy in range(n_y):
            for ix in range(n_x):
                if len(origins) == n_water:
                    break
                origins.append(
                    ((ix + 0.5) * dx, (iy + 0.5) * dy, z0 + (iz + 0.5) * dz)
                )
    return [regen.tip4p_at(origin, rng) for origin in origins]


def pairs_for(n_water, molality):
    if molality <= 0 or n_water <= 0:
        return 0
    return max(1, int(round(molality * n_water * regen.MOLAL_FACTOR)))


def build_growth(nx, ny, nz, molality, seed, gap=2.6):
    """Ice Ih slab, liquid of the same water count, ions only in the liquid.

    The ice is built as a periodic crystal, so every oxygen owns two
    hydrogens. The box is then opened along z and the liquid is inserted.
    Hydrogens that used to bond across the old z boundary point into the gap.
    """
    rng = random.Random(seed)
    ice, ice_box = periodic_ice(nx, ny, nz, seed)
    n_liquid = len(ice)
    area = ice_box[0] * ice_box[1]
    liquid_height = n_liquid * 18.01528 / (regen.RHO * area)
    z0 = ice_box[2] + gap
    z1 = z0 + liquid_height
    box = (ice_box[0], ice_box[1], z1 + gap)
    liquid = liquid_waters(n_liquid, box, z0, z1, rng)
    waters = ice + liquid
    n_pairs = pairs_for(n_liquid, molality)
    oxygens = [trip[0] for trip in waters]
    # Ions see the full box for the minimum image, but they are drawn
    # only inside the liquid slab.
    ions = []
    if n_pairs:
        cell = 2.2
        dims = [max(1, math.ceil(box[k] / cell)) for k in range(3)]

        def bucket(pos):
            return tuple(min(dims[k] - 1, max(0, int(pos[k] / cell))) for k in range(3))

        buckets = {}
        for pos in oxygens:
            buckets.setdefault(bucket(pos), []).append(pos)

        def free(pos):
            center = bucket(pos)
            for ix in range(center[0] - 1, center[0] + 2):
                for iy in range(center[1] - 1, center[1] + 2):
                    for iz in range(center[2] - 1, center[2] + 2):
                        for other in buckets.get((ix % dims[0], iy % dims[1], iz % dims[2]), ()):
                            if regen.mic(pos, other, box) <= 2.2:
                                return False
            return True

        for _kind in range(2):
            made = 0
            tries = 0
            while made < n_pairs and tries < n_pairs * 2000:
                tries += 1
                pos = (
                    rng.random() * box[0],
                    rng.random() * box[1],
                    z0 + rng.random() * (z1 - z0),
                )
                if free(pos):
                    ions.append(pos)
                    buckets.setdefault(bucket(pos), []).append(pos)
                    made += 1
            if made < n_pairs:
                raise RuntimeError(f"placed {made} of {n_pairs} ions at {molality} mol/kg")
    na, cl = ions[:n_pairs], ions[n_pairs:]
    oxygens = [trip[0] for trip in waters]
    for i, origin in enumerate(oxygens):
        for other in oxygens[i + 1 :]:
            if regen.mic(origin, other, box) < 2.3:
                raise RuntimeError("two oxygens are closer than 2.3 A")
    return waters, na, cl, box, len(ice)


def write_growth_input(path, n_ice, temperature, equil_steps, prod_steps, dump_every, dump_name):
    text = f"""# TIP4P/2005 + Madrid-2019. Ice molecules 1 through {n_ice} are
# pinned for the equilibration, then released.
units real
atom_style full
boundary p p p
read_data data.brine

pair_style lj/cut/tip4p/long 1 2 1 1 0.1546 10.0
pair_modify mix arithmetic
kspace_style pppm/tip4p 1.0e-4
pair_coeff 1 1 0.1852 3.1589
pair_coeff 2 2 0.0 1.0
pair_coeff 3 3 0.351853250 2.21737
pair_coeff 4 4 0.018384799 4.69906
pair_coeff 1 3 0.189624283 2.60838
pair_coeff 1 4 0.014814293 4.23867
pair_coeff 3 4 0.343856119 3.00512
pair_coeff 1 2 0.0 1.0
pair_coeff 2 3 0.0 1.0
pair_coeff 2 4 0.0 1.0

bond_style harmonic
bond_coeff 1 450.0 0.9572
angle_style harmonic
angle_coeff 1 55.0 104.52
special_bonds lj/coul 0.0 0.0 0.0

neighbor 2.0 bin
neigh_modify delay 0 every 1 check yes
timestep 2.0
thermo_style custom step temp pe density press
thermo 500

group ice molecule 1:{n_ice}
group mobile subtract all ice
fix shake all shake 0.0001 20 0 b 1 a 1
fix freeze ice setforce 0.0 0.0 0.0
minimize 1.0e-4 1.0e-6 40 400
velocity mobile create {temperature:.2f} 12345 dist gaussian
velocity ice set 0.0 0.0 0.0
fix eq mobile nvt temp {temperature:.2f} {temperature:.2f} 200.0
run {equil_steps}
unfix eq
unfix freeze
velocity ice create {temperature:.2f} 54321 dist gaussian
reset_timestep 0
fix prod all npt temp {temperature:.2f} {temperature:.2f} 200.0 z 1.0 1.0 1000.0
dump dmp all custom {dump_every} {dump_name} id mol type x y z
dump_modify dmp sort id
run {prod_steps}
"""
    path.write_text(text)


def run_lammps(lmp, work, script):
    log = work / "lammps.log"
    env = dict(**{k: v for k, v in __import__("os").environ.items()})
    env["OMP_NUM_THREADS"] = "1"
    with log.open("w") as handle:
        proc = subprocess.run(
            [lmp, "-in", script.name],
            cwd=work,
            env=env,
            stdout=handle,
            stderr=subprocess.STDOUT,
        )
    if proc.returncode != 0:
        tail = log.read_text(errors="replace")[-2500:]
        raise SystemExit(f"LAMMPS failed in {work}\n{tail}")
    return log


def thermo_rows(log):
    rows = []
    for line in log.read_text(errors="replace").splitlines():
        parts = line.split()
        if len(parts) < 5:
            continue
        try:
            step = int(parts[0])
            temp = float(parts[1])
            density = float(parts[3])
        except ValueError:
            continue
        rows.append((step, temp, density))
    return rows


def block_mean(values, blocks):
    if not values:
        return float("nan"), float("nan")
    if len(values) < blocks:
        blocks = len(values)
    width = len(values) // blocks
    means = []
    for i in range(blocks):
        chunk = values[i * width : (i + 1) * width]
        if chunk:
            means.append(sum(chunk) / len(chunk))
    mean = sum(means) / len(means)
    if len(means) < 2:
        return mean, float("nan")
    var = sum((item - mean) ** 2 for item in means) / (len(means) - 1)
    return mean, math.sqrt(var / len(means))


def ice_check():
    waters, box = periodic_ice(2, 2, 2, seed=1)
    # Every oxygen has two hydrogens within 1.15 A, and no two oxygens
    # are closer than 2.4 A.
    close_h = 0
    for oxygen, h1, h2 in waters:
        if regen.mic(oxygen, h1, box) > 1.15 or regen.mic(oxygen, h2, box) > 1.15:
            close_h += 1
    near = 0
    oxygens = [trip[0] for trip in waters]
    for i, origin in enumerate(oxygens):
        for other in oxygens[i + 1 :]:
            if regen.mic(origin, other, box) < 2.4:
                near += 1
    if close_h or near:
        raise SystemExit(f"ice geometry failed: bad hydrogens {close_h}, close oxygens {near}")
    print(f"ice rules ok: {len(waters)} molecules, box {box[0]:.2f} {box[1]:.2f} {box[2]:.2f}")


def crystal(lmp, steps):
    n = 4
    # Table V density 2.050 g/cm^3 is looser than the experimental
    # 5.64 A cell. Start outside the steep Na-Cl core and let NPT settle.
    a = 6.4
    na_basis = ((0, 0, 0), (0.5, 0.5, 0), (0.5, 0, 0.5), (0, 0.5, 0.5))
    cl_basis = ((0.5, 0, 0), (0, 0.5, 0), (0, 0, 0.5), (0.5, 0.5, 0.5))
    sodium, chloride = [], []
    for ix in range(n):
        for iy in range(n):
            for iz in range(n):
                origin = (ix, iy, iz)
                for basis in na_basis:
                    sodium.append(tuple((origin[k] + basis[k]) * a for k in range(3)))
                for basis in cl_basis:
                    chloride.append(tuple((origin[k] + basis[k]) * a for k in range(3)))
    length = n * a
    work = OUT / "crystal-madrid2019"
    work.mkdir(parents=True, exist_ok=True)
    lines = [
        "Madrid-2019 NaCl crystal, scaled charges\n\n",
        f"{len(sodium) + len(chloride)} atoms\n\n",
        "2 atom types\n\n",
        f"0.0 {length:.6f} xlo xhi\n",
        f"0.0 {length:.6f} ylo yhi\n",
        f"0.0 {length:.6f} zlo zhi\n\n",
        "Masses\n\n1 22.989769\n2 35.453\n\n",
        "Atoms # charge\n\n",
    ]
    atom = 1
    for pos in sodium:
        lines.append(f"{atom} 1 0.85 {pos[0]:.5f} {pos[1]:.5f} {pos[2]:.5f}\n")
        atom += 1
    for pos in chloride:
        lines.append(f"{atom} 2 -0.85 {pos[0]:.5f} {pos[1]:.5f} {pos[2]:.5f}\n")
        atom += 1
    (work / "data.crystal").write_text("".join(lines))
    script = f"""units real
atom_style charge
boundary p p p
read_data data.crystal
pair_style lj/cut/coul/long 10.0
pair_modify mix arithmetic
kspace_style pppm 1.0e-4
pair_coeff 1 1 0.351853250 2.21737
pair_coeff 2 2 0.018384799 4.69906
pair_coeff 1 2 0.343856119 3.00512
neighbor 2.0 bin
neigh_modify delay 0 every 1 check yes
timestep 1.0
thermo_style custom step temp pe density press
thermo 200
minimize 1.0e-4 1.0e-6 200 2000
velocity all create 298.15 12345 dist gaussian
fix eq all npt temp 298.15 298.15 100.0 iso 1.0 1.0 1000.0
run {steps}
"""
    (work / "in.crystal").write_text(script)
    log = run_lammps(lmp, work, work / "in.crystal")
    rows = thermo_rows(log)
    # Drop the first half. Density is g/cm^3. Table V quotes 2050 kg/m^3.
    sample = [row[2] for row in rows[len(rows) // 2 :]]
    mean, err = block_mean(sample, 6)
    print(
        f"crystal density {mean:.4f} +/- {err:.4f} g/cm^3 "
        f"({mean * 1000:.1f} +/- {err * 1000:.1f} kg/m^3), "
        f"published 2.050 g/cm^3, n_thermo {len(sample)}"
    )
    (work / "density.txt").write_text(f"{mean:.6f} {err:.6f} {len(sample)}\n")
    return mean, err


def grow(lmp, brine_front, nx, ny, nz, prod_ps, equil_ps, molalities, temperatures):
    OUT.mkdir(parents=True, exist_ok=True)
    dump_every = 5000  # 10 ps at 2 fs
    equil_steps = int(round(equil_ps * 500))
    prod_steps = int(round(prod_ps * 500))
    summary = OUT / "summary.txt"
    lines = []
    for temperature in temperatures:
        for molality in molalities:
            name = f"m{molality:.1f}-T{temperature:.0f}"
            work = OUT / name
            if work.exists():
                shutil.rmtree(work)
            work.mkdir(parents=True)
            waters, na, cl, box, n_ice = build_growth(nx, ny, nz, molality, seed=1000 + int(temperature) + int(molality * 100))
            regen.write_data(work / "data.brine", waters, na, cl, box)
            write_growth_input(
                work / "in.growth",
                n_ice,
                temperature,
                equil_steps,
                prod_steps,
                dump_every,
                "dump.lammpstrj",
            )
            print(f"running {name} waters={len(waters)} ice={n_ice} ions={len(na)+len(cl)}", flush=True)
            run_lammps(lmp, work, work / "in.growth")
            report = work / "front.txt"
            cmd = [
                brine_front,
                str(work / "dump.lammpstrj"),
                "--dt-ps",
                "0.002",
                "--every",
                str(dump_every),
                "--time",
                "8",
                "--bins",
                "20",
            ]
            with report.open("w") as handle:
                proc = subprocess.run(cmd, stdout=handle, stderr=subprocess.STDOUT)
            if proc.returncode != 0:
                tail = report.read_text(errors="replace")[-1500:]
                raise SystemExit(f"brine_front failed for {name}\n{tail}")
            text = report.read_text()
            summary_line = next(line for line in text.splitlines() if line.startswith("summary "))
            timing = next((line for line in text.splitlines() if line.startswith("timing_ms ")), "")
            record = f"{name} {summary_line} {timing}"
            print(record, flush=True)
            lines.append(record)
            summary.write_text("\n".join(lines) + "\n")
    print(f"wrote {summary}")


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--ice-check", action="store_true")
    parser.add_argument("--crystal", action="store_true")
    parser.add_argument("--grow", action="store_true")
    parser.add_argument("--lmp", default=os.environ.get("LAMMPS_BIN", "lmp"))
    parser.add_argument("--brine-front", default="")
    parser.add_argument("--crystal-steps", type=int, default=40000)
    parser.add_argument("--prod-ps", type=float, default=400)
    parser.add_argument("--equil-ps", type=float, default=20)
    parser.add_argument("--nx", type=int, default=3)
    parser.add_argument("--ny", type=int, default=3)
    parser.add_argument("--nz", type=int, default=2)
    parser.add_argument("--molalities", default="0,0.5,1,2")
    parser.add_argument("--temperatures", default="236,220")
    args = parser.parse_args()
    if args.ice_check:
        ice_check()
    if args.crystal:
        crystal(args.lmp, args.crystal_steps)
    if args.grow:
        binary = args.brine_front
        if not binary:
            raise SystemExit("pass --brine-front")
        molalities = [float(part) for part in args.molalities.split(",") if part]
        temperatures = [float(part) for part in args.temperatures.split(",") if part]
        grow(
            args.lmp,
            binary,
            args.nx,
            args.ny,
            args.nz,
            args.prod_ps,
            args.equil_ps,
            molalities,
            temperatures,
        )
    if not (args.ice_check or args.crystal or args.grow):
        parser.print_help()


if __name__ == "__main__":
    main()
