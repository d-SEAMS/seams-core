#!/usr/bin/env python3
"""Direct-coexistence slabs for Madrid-2019 + TIP4P/2005.

The geometry follows Lamas, Vega and Noya, J. Chem. Phys. 156, 134503
(2022): an ice Ih block in contact with a liquid, cross section at
least 30 A, at least 2000 waters, cutoff 10 A, 2 fs, 1 bar. Frames
are written only under /tmp.

    python3 input/brine/coexist.py --bench --threads 8
    python3 input/brine/coexist.py --run --tag pure-245 --temp 245 --molality 0 \
        --threads 8 --melt-ps 20 --prod-ps 2
"""

from __future__ import annotations

import argparse
import random
import shutil
import subprocess
import sys
import time
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))
import grow
import regenerate as regen

OUT = Path("/tmp/brine-science-v2")
LMP = "/tmp/lammps/build-omp/lmp"


def build_slab(nx, ny, nz, n_liquid, molality, seed, gap=3.0):
    rng = random.Random(seed)
    ice, ice_box = grow.periodic_ice(nx, ny, nz, seed)
    area = ice_box[0] * ice_box[1]
    height = n_liquid * 18.01528 / (regen.RHO * area)
    z0 = ice_box[2] + gap
    z1 = z0 + height
    box = (ice_box[0], ice_box[1], z1 + gap)
    liquid = grow.liquid_waters(n_liquid, box, z0, z1, rng)
    # A perfect liquid grid crystallises too readily. Jitter the oxygens
    # by a few tenths of an angstrom and keep the molecule rigid.
    jittered = []
    for oxygen, h1, h2 in liquid:
        shift = tuple((rng.random() - 0.5) * 0.6 for _ in range(3))
        jittered.append(
            (
                tuple(oxygen[k] + shift[k] for k in range(3)),
                tuple(h1[k] + shift[k] for k in range(3)),
                tuple(h2[k] + shift[k] for k in range(3)),
            )
        )
    waters = ice + jittered
    n_pairs = 0 if molality <= 0 else int(round(molality * n_liquid * regen.MOLAL_FACTOR))
    ions = []
    if n_pairs:
        oxygens = [trip[0] for trip in waters]

        def far(pos):
            return all(regen.mic(pos, other, box) > 2.3 for other in oxygens[::3])

        # Spot-check against every oxygen is quadratic and slow at this
        # size. The cell list in grow.place is not exported, so sample
        # against a stride and then against the ice surface band.
        surface = [p for p in oxygens if p[2] > ice_box[2] - 4.0 or p[2] < 4.0]
        made = 0
        guard = 0
        while made < 2 * n_pairs and guard < n_pairs * 5000:
            guard += 1
            pos = (
                rng.random() * box[0],
                rng.random() * box[1],
                z0 + 2.0 + rng.random() * max(z1 - z0 - 4.0, 1.0),
            )
            if all(regen.mic(pos, other, box) > 2.4 for other in surface):
                ions.append(pos)
                surface.append(pos)
                made += 1
        if made < 2 * n_pairs:
            raise RuntimeError(f"placed {made} of {2 * n_pairs} ions")
    na, cl = ions[:n_pairs], ions[n_pairs:]
    return waters, na, cl, box, len(ice)


def write_input(path, n_ice, temperature, melt_steps, cool_steps, prod_steps, dump_every, seed, threads):
    path.write_text(
        f"""# Madrid-2019 + TIP4P/2005 direct coexistence.
# Ice molecules 1..{n_ice} are pinned while the liquid is melted at 300 K,
# then the whole slab runs at {temperature:.2f} K.
units real
atom_style full
boundary p p p
processors * * *
package omp {threads}
suffix omp
read_data data.brine

pair_style lj/cut/tip4p/long 1 2 1 1 0.1546 10.0
pair_modify mix arithmetic tail yes
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
# Shake and the minimiser fight each other on a fresh interface.
# The harmonic bonds take the overlaps out, then Shake is turned on.
fix freeze ice setforce 0.0 0.0 0.0
minimize 1.0e-3 1.0e-4 40 200
fix shake all shake 0.0001 20 0 b 1 a 1
velocity mobile create 300.0 {seed} dist gaussian
velocity ice set 0.0 0.0 0.0
fix melt mobile nvt temp 300.0 300.0 200.0
run {melt_steps}
unfix melt
fix cool mobile nvt temp {temperature:.2f} {temperature:.2f} 200.0
run {cool_steps}
unfix cool
unfix freeze
velocity ice create {temperature:.2f} {seed + 17} dist gaussian
reset_timestep 0
fix prod all npt temp {temperature:.2f} {temperature:.2f} 200.0 &
    x 1.0 1.0 2000.0 y 1.0 1.0 2000.0 z 1.0 1.0 2000.0
dump dmp all custom {dump_every} dump.lammpstrj id mol type x y z vx vy vz
dump_modify dmp sort id
# 0.25 ns at a 2 fs step. Two files so one is always complete.
restart 125000 restart.a restart.b
run {prod_steps}
"""
    )


def run_lammps(work, threads):
    log = work / "lammps.log"
    env = dict(**{k: v for k, v in __import__("os").environ.items()})
    env["OMP_NUM_THREADS"] = str(threads)
    started = time.time()
    with log.open("w") as handle:
        proc = subprocess.run(
            [LMP, "-sf", "omp", "-pk", f"omp {threads}", "-in", "in.coexist"],
            cwd=work,
            env=env,
            stdout=handle,
            stderr=subprocess.STDOUT,
        )
    wall = time.time() - started
    (work / "wall.txt").write_text(f"{wall:.3f}\n")
    if proc.returncode != 0:
        tail = log.read_text(errors="replace")[-2000:]
        raise SystemExit(f"LAMMPS failed in {work} after {wall:.1f} s\n{tail}")
    return wall


def loop_time(log):
    """Return (steps, loop seconds) from the last 'Loop time' block."""
    text = log.read_text(errors="replace").splitlines()
    steps = None
    loop = None
    for line in text:
        if line.startswith("Loop time of"):
            # Loop time of 1.23 on 1 procs for 500 steps with ...
            parts = line.split()
            loop = float(parts[3])
            steps = int(parts[8])
    return steps, loop


def prepare(tag, nx, ny, nz, n_liquid, molality, temperature, seed, threads, melt_ps, prod_ps, dump_ps):
    work = OUT / tag
    if work.exists():
        shutil.rmtree(work)
    work.mkdir(parents=True)
    waters, na, cl, box, n_ice = build_slab(nx, ny, nz, n_liquid, molality, seed)
    regen.write_data(work / "data.brine", waters, na, cl, box)
    melt_steps = int(round(melt_ps * 500))
    prod_steps = int(round(prod_ps * 500))
    cool_steps = 200 if tag == "bench" else 5000
    dump_every = max(500, int(round(dump_ps * 500)))
    write_input(
        work / "in.coexist",
        n_ice,
        temperature,
        melt_steps,
        cool_steps,
        prod_steps,
        dump_every,
        seed,
        threads,
    )
    meta = (
        f"waters {len(waters)} ice {n_ice} liquid {len(waters) - n_ice} "
        f"ions {len(na) + len(cl)} box {box[0]:.3f} {box[1]:.3f} {box[2]:.3f} "
        f"m {molality} T {temperature} seed {seed}\n"
    )
    (work / "meta.txt").write_text(meta)
    print(meta.strip(), flush=True)
    return work


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--bench", action="store_true")
    parser.add_argument("--run", action="store_true")
    parser.add_argument("--tag", default="bench")
    parser.add_argument("--temp", type=float, default=245.0)
    parser.add_argument("--molality", type=float, default=0.0)
    parser.add_argument("--seed", type=int, default=247)
    parser.add_argument("--threads", type=int, default=8)
    parser.add_argument("--nx", type=int, default=4)
    parser.add_argument("--ny", type=int, default=5)
    parser.add_argument("--nz", type=int, default=4)
    parser.add_argument("--n-liquid", type=int, default=2048)
    parser.add_argument("--melt-ps", type=float, default=20.0)
    parser.add_argument("--prod-ps", type=float, default=0.2)
    parser.add_argument("--dump-ps", type=float, default=10.0)
    parser.add_argument("--steps", type=int, default=400)
    args = parser.parse_args()
    if args.bench:
        # Time the production integrator only. Melt is skipped by a short
        # prod after a token melt so the neighbour list is realistic.
        work = prepare(
            "bench",
            args.nx,
            args.ny,
            args.nz,
            args.n_liquid,
            0.0,
            250.0,
            1,
            args.threads,
            melt_ps=0.2,
            prod_ps=args.steps * 0.002,
            dump_ps=1000.0,
        )
        wall = run_lammps(work, args.threads)
        steps, loop = loop_time(work / "lammps.log")
        rate = steps / loop if loop else float("nan")
        print(
            f"bench threads {args.threads} steps {steps} loop_s {loop:.2f} "
            f"steps_per_s {rate:.2f} wall_s {wall:.1f} "
            f"ns_per_day {rate * 2e-6 * 86400:.2f}",
            flush=True,
        )
        return
    if args.run:
        work = prepare(
            args.tag,
            args.nx,
            args.ny,
            args.nz,
            args.n_liquid,
            args.molality,
            args.temp,
            args.seed,
            args.threads,
            args.melt_ps,
            args.prod_ps,
            args.dump_ps,
        )
        wall = run_lammps(work, args.threads)
        steps, loop = loop_time(work / "lammps.log")
        print(f"{args.tag} wall_s {wall:.1f} last_loop_steps {steps} last_loop_s {loop}", flush=True)


if __name__ == "__main__":
    main()
