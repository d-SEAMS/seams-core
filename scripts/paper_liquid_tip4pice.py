#!/usr/bin/env python3
"""Write oxygen coordinates for a TIP4P/ice liquid.

The nonbonded parameters are the TIP4P/ice force field shipped with
the in-situ adapter. A 4096-molecule box is melted at 300 K. The
sitting tiles that frame 6 by 6 by 6, which is 884736 molecules.
GenIce is not required. The LAMMPS run stays outside the source tree.
The oxygen list is plain text.
"""

from __future__ import annotations

import argparse
import math
import os
import random
import subprocess
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "input" / "brine"))
import regenerate as brine  # noqa: E402

# adapters/python/dseams_insitu/tip4pice.xml
EPSILON_KCAL = 0.882164484 / 4.184
SIGMA = 3.1668
Q_H = 0.5897
Q_M = -1.1794
QDIST = 0.1577
OH = 0.9572
HOH = math.radians(104.52)
N_WATER = 4096
STEPS = 5000
TEMPERATURE = 300.0
SEED = 884736


def tip4p_ice_at(origin, rng):
    h1 = (OH, 0.0, 0.0)
    h2 = (OH * math.cos(HOH), OH * math.sin(HOH), 0.0)
    rot = brine.rotation(rng)
    return origin, brine.add(origin, brine.apply_rot(rot, h1)), brine.add(origin, brine.apply_rot(rot, h2))


def write_data(path: Path, waters, box) -> None:
    n_water = len(waters)
    lines = [
        "TIP4P/ice liquid\n\n",
        f"{n_water * 3} atoms\n",
        f"{n_water * 2} bonds\n",
        f"{n_water} angles\n\n",
        "2 atom types\n1 bond types\n1 angle types\n\n",
        f"0.0 {box[0]:.6f} xlo xhi\n",
        f"0.0 {box[1]:.6f} ylo yhi\n",
        f"0.0 {box[2]:.6f} zlo zhi\n\n",
        "Masses\n\n1 15.9994\n2 1.008\n\n",
        "Atoms # full\n\n",
    ]
    atom_id = 1
    bonds = []
    angles = []
    for mol, (ox, h1, h2) in enumerate(waters, start=1):
        o_id, h1_id, h2_id = atom_id, atom_id + 1, atom_id + 2
        lines.append(f"{o_id} {mol} 1 {Q_M:.4f} {ox[0]:.5f} {ox[1]:.5f} {ox[2]:.5f}\n")
        lines.append(f"{h1_id} {mol} 2 {Q_H:.4f} {h1[0]:.5f} {h1[1]:.5f} {h1[2]:.5f}\n")
        lines.append(f"{h2_id} {mol} 2 {Q_H:.4f} {h2[0]:.5f} {h2[1]:.5f} {h2[2]:.5f}\n")
        bonds.append((o_id, h1_id))
        bonds.append((o_id, h2_id))
        angles.append((h1_id, o_id, h2_id))
        atom_id += 3
    lines.append("\nBonds\n\n")
    for i, (a, b) in enumerate(bonds, start=1):
        lines.append(f"{i} 1 {a} {b}\n")
    lines.append("\nAngles\n\n")
    for i, (a, b, c) in enumerate(angles, start=1):
        lines.append(f"{i} 1 {a} {b} {c}\n")
    path.write_text("".join(lines))


def deck(steps: int) -> str:
    return f"""units real
atom_style full
boundary p p p
atom_modify map array
read_data data.tip4pice
pair_style lj/cut/tip4p/long 1 2 1 1 {QDIST:.4f} 10.0
pair_modify mix arithmetic
kspace_style pppm/tip4p 1.0e-4
pair_coeff 1 1 {EPSILON_KCAL:.8f} {SIGMA:.4f}
pair_coeff 2 2 0.0 1.0
bond_style harmonic
bond_coeff 1 450.0 {OH:.4f}
angle_style harmonic
angle_coeff 1 55.0 104.52
special_bonds lj/coul 0.0 0.0 0.0
reset_atoms image all
neighbor 2.0 bin
timestep 2.0
thermo 500
fix shake all shake 0.0001 20 0 b 1 a 1
minimize 1.0e-4 1.0e-6 100 400
velocity all create {TEMPERATURE:.1f} {SEED} mom yes rot yes dist gaussian
fix bath all nvt temp {TEMPERATURE:.1f} {TEMPERATURE:.1f} 200.0
run {steps}
group oxygen type 1
write_dump oxygen custom oxy.lammpstrj id type x y z modify sort id
"""


def parse_dump(path: Path):
    lines = path.read_text().splitlines()
    box = []
    atoms = []
    mode = ""
    for line in lines:
        if line.startswith("ITEM: BOX BOUNDS"):
            mode = "box"
            box = []
            continue
        if line.startswith("ITEM: ATOMS"):
            mode = "atoms"
            atoms = []
            continue
        if line.startswith("ITEM:"):
            mode = ""
            continue
        if mode == "box":
            lo, hi = line.split()[:2]
            box.append(float(hi) - float(lo))
        elif mode == "atoms":
            parts = line.split()
            atoms.append((int(parts[0]), float(parts[2]), float(parts[3]), float(parts[4])))
    atoms.sort()
    return box, atoms


def write_oxygens(path: Path, box, atoms) -> None:
    side = box[0]
    if max(abs(box[i] - side) for i in range(3)) > 1e-4:
        raise RuntimeError(f"box is not cubic: {box}")
    rows = [f"{len(atoms)} {side:.10f} {side:.10f} {side:.10f}"]
    for _ident, x, y, z in atoms:
        x -= side * math.floor(x / side)
        y -= side * math.floor(y / side)
        z -= side * math.floor(z / side)
        rows.append(f"{x:.6f} {y:.6f} {z:.6f}")
    path.write_text("\n".join(rows) + "\n")


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--lmp", default=os.environ.get("LMP", "lmp"))
    parser.add_argument("--out", type=Path, default=ROOT / "benchmarks/paper-cells/tip4pice-liquid-4096.txt")
    parser.add_argument(
        "--work", type=Path, default=ROOT / "benchmarks/paper-cells/.liquid-work")
    parser.add_argument("--steps", type=int, default=STEPS)
    parser.add_argument("--threads", type=int, default=4)
    args = parser.parse_args()
    args.work.mkdir(parents=True, exist_ok=True)
    rng = random.Random(SEED)
    origins, box = brine.liquid_origins(N_WATER, rng)
    waters = [tip4p_ice_at(origin, rng) for origin in origins]
    write_data(args.work / "data.tip4pice", waters, box)
    (args.work / "in.liquid").write_text(deck(args.steps))
    env = dict(os.environ)
    env["UCX_TLS"] = "tcp,self,sm"
    env["OMPI_MCA_btl"] = "^openib"
    env["OMPI_MCA_pml"] = "ob1"
    cmd = [
        args.lmp, "-sf", "omp", "-pk", "omp", str(args.threads), "-nocite",
        "-in", "in.liquid", "-log", "log.lammps", "-screen", "none",
    ]
    subprocess.run(cmd, cwd=args.work, env=env, check=True)
    box, atoms = parse_dump(args.work / "oxy.lammpstrj")
    if len(atoms) != N_WATER:
        raise RuntimeError(f"dump has {len(atoms)} oxygens, need {N_WATER}")
    args.out.parent.mkdir(parents=True, exist_ok=True)
    write_oxygens(args.out, box, atoms)
    print(f"wrote {args.out} molecules {len(atoms)} side {box[0]:.6f}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
