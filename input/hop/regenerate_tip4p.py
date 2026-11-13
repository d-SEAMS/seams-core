#!/usr/bin/env python3
# Copyright (c) d-SEAMS developers

"""Regenerate production-stride TIP4P/2005 ice/water trajectories.

The slab is 4096 molecules: 2048 of proton-ordered ice Ih and 2048 of
liquid water. The ice cell is the 16-molecule GenIce Ih cell used by
input/brine. The lateral box is 31.4 by 59.0 A, and the ice is about
36 A thick. A 3 A gap separates the crystal from the liquid and the
liquid from its periodic image.

TIP4P/2005 melts at 252 K (Abascal and Vega, J. Chem. Phys. 123,
234505, 2005), which is the interface temperature used here. The
liquid is melted at 300 K with the ice pinned, cooled to 252 K, and
the whole slab is then equilibrated in NPT at 1 atm with the barostat
on z only. The step is 2 fs.

mW production strides were 0.1, 1 and 10 ps on a 10 fs step. TIP4P
water is integrated ten times finer, and an oxygen moves farther in
one picosecond than an mW site does. The dumps here are every 0.01,
0.1 and 1 ps (5, 50 and 500 steps). 0.01 ps is a few steps, 0.1 ps is
a fraction of a hydrogen-bond lifetime, and 1 ps is a production dump
for an ice/water interface. Each file has the equilibrated frame and
200 transitions. Only oxygens are written. The dumps stay outside git.

Usage:
  python3 input/hop/regenerate_tip4p.py
  python3 input/hop/regenerate_tip4p.py --lmp lmp --strides 5,50,500 --transitions 200
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

ROOT = Path(__file__).resolve().parents[2]
DECK = Path(__file__).resolve().parent
OUT = ROOT / "data" / "hop"
BRINE = ROOT / "input" / "brine"

sys.path.insert(0, str(BRINE))
import grow  # noqa: E402
import regenerate as brine  # noqa: E402

# 4 x 8 x 4 of the 16-molecule cell is 2048 ice molecules.
NX, NY, NZ = 4, 8, 4
N_ICE = 16 * NX * NY * NZ
GAP = 3.0
# Abascal and Vega, J. Chem. Phys. 123, 234505 (2005).
TEMPERATURE = 252.0
TIMESTEP_FS = 2.0
# 0.01 ps, 0.1 ps, 1 ps at a 2 fs step.
DEFAULT_STRIDES = "5,50,500"


def stride_ps(steps: int) -> float:
    return steps * TIMESTEP_FS * 1.0e-3


def min_image(p, q, box):
    acc = 0.0
    for k in range(3):
        d = q[k] - p[k]
        d -= box[k] * round(d / box[k])
        acc += d * d
    return math.sqrt(acc)


def cell_pairs(origins, box, cutoff):
    """True when two origins lie closer than cutoff."""
    cell = cutoff
    dims = [max(1, math.ceil(box[k] / cell)) for k in range(3)]

    def key(pos):
        return tuple(min(dims[k] - 1, max(0, int(pos[k] / cell))) for k in range(3))

    buckets = {}
    for i, pos in enumerate(origins):
        buckets.setdefault(key(pos), []).append(i)
    for i, pos in enumerate(origins):
        center = key(pos)
        for ix in range(center[0] - 1, center[0] + 2):
            for iy in range(center[1] - 1, center[1] + 2):
                for iz in range(center[2] - 1, center[2] + 2):
                    for j in buckets.get((ix % dims[0], iy % dims[1], iz % dims[2]), ()):
                        if j <= i:
                            continue
                        if min_image(pos, origins[j], box) < cutoff:
                            return True
    return False


def nearest(origin, other, box):
    """Place other on the periodic image closest to origin."""
    return tuple(
        other[k] - box[k] * round((other[k] - origin[k]) / box[k]) for k in range(3)
    )


def whole_molecule(ox, h1, h2, box):
    """Oxygen inside the box, hydrogens on the nearby image.

    A hydrogen may sit just outside a face. The bond is then the
    chemical one, and LAMMPS stores that crossing as an image flag.
    """
    shift = tuple(box[k] * math.floor(ox[k] / box[k]) for k in range(3))
    ox = tuple(ox[k] - shift[k] for k in range(3))
    h1 = nearest(ox, tuple(h1[k] - shift[k] for k in range(3)), box)
    h2 = nearest(ox, tuple(h2[k] - shift[k] for k in range(3)), box)
    return ox, h1, h2


def molecules_intact(waters) -> None:
    bad = 0
    for ox, h1, h2 in waters:
        d1 = math.dist(ox, h1)
        d2 = math.dist(ox, h2)
        if not (0.90 < d1 < 1.05 and 0.90 < d2 < 1.05):
            bad += 1
    if bad:
        raise SystemExit(f"{bad} molecules are not rigid TIP4P water")


def build_slab(seed: int):
    rng = random.Random(seed)
    ice, ice_box = grow.periodic_ice(NX, NY, NZ, seed)
    if len(ice) != N_ICE:
        raise SystemExit(f"ice count {len(ice)} != {N_ICE}")
    # periodic_ice wraps each atom into the crystal cell, which sends a
    # hydrogen that crossed a face to the far side of its oxygen.
    ice = [whole_molecule(ox, h1, h2, ice_box) for ox, h1, h2 in ice]
    n_liquid = len(ice)
    area = ice_box[0] * ice_box[1]
    volume = n_liquid * 18.01528 / brine.RHO
    height = volume / area
    z0 = ice_box[2] + GAP
    z1 = z0 + height
    box = (ice_box[0], ice_box[1], z1 + GAP)
    liquid = grow.liquid_waters(n_liquid, box, z0, z1, rng)
    jittered = []
    for oxygen, h1, h2 in liquid:
        placed = None
        for _try in range(40):
            shift = tuple((rng.random() - 0.5) * 0.5 for _ in range(3))
            moved = tuple(oxygen[k] + shift[k] for k in range(3))
            if z0 + 0.4 <= moved[2] <= z1 - 0.4 and 0.0 <= moved[0] < box[0] and 0.0 <= moved[1] < box[1]:
                placed = shift
                break
        if placed is None:
            placed = (0.0, 0.0, 0.0)
        jittered.append(
            tuple(tuple(atom[k] + placed[k] for k in range(3)) for atom in (oxygen, h1, h2))
        )
    waters = [whole_molecule(*mol, box) for mol in ice + jittered]
    molecules_intact(waters)
    origins = [mol[0] for mol in waters]
    if cell_pairs(origins, box, 2.2):
        raise SystemExit("two oxygens are closer than 2.2 A")
    return waters, box, len(ice), ice_box, height


def write_data(path: Path, waters, box) -> None:
    n_water = len(waters)
    lines = [
        "LAMMPS data file for a TIP4P/2005 ice/water slab\n\n",
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
        lines.append(
            f"{o_id} {mol} 1 {brine.Q_O:.4f} {ox[0]:.5f} {ox[1]:.5f} {ox[2]:.5f}\n"
        )
        lines.append(
            f"{h1_id} {mol} 2 {brine.Q_H:.4f} {h1[0]:.5f} {h1[1]:.5f} {h1[2]:.5f}\n"
        )
        lines.append(
            f"{h2_id} {mol} 2 {brine.Q_H:.4f} {h2[0]:.5f} {h2[1]:.5f} {h2[2]:.5f}\n"
        )
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


def run_lammps(lmp: str, deck: str, variables: dict[str, str], log_name: str, threads: int) -> None:
    cmd = [lmp, "-sf", "omp", "-pk", "omp", str(threads), "-nocite",
           "-in", deck, "-log", log_name, "-screen", "none"]
    for key, value in variables.items():
        cmd.extend(["-var", key, value])
    env = os.environ.copy()
    env["OMP_NUM_THREADS"] = str(threads)
    # A single rank. The conda-forge build is MPI, and an unset UCX
    # transport hangs on a host without RDMA. OpenMP is the pair style.
    env.setdefault("UCX_TLS", "tcp,self,sm")
    env.setdefault("OMPI_MCA_btl", "^openib")
    env.setdefault("OMPI_MCA_pml", "ob1")
    env.setdefault("OMPI_MCA_hwloc_base_binding_policy", "none")
    subprocess.run(cmd, cwd=OUT, env=env, check=True)


def count_frames(path: Path) -> tuple[int, int]:
    frames = 0
    atoms = None
    pending = False
    with path.open() as handle:
        for line in handle:
            if line.startswith("ITEM: TIMESTEP"):
                frames += 1
                pending = False
            elif line.startswith("ITEM: NUMBER OF ATOMS"):
                pending = True
            elif pending:
                atoms = int(line)
                pending = False
    return frames, atoms if atoms is not None else 0


def thermo_tail(log_path: Path, n: int = 5) -> list[tuple[int, float, float]]:
    """Last thermo samples as (step, temperature, density)."""
    rows = []
    for line in log_path.read_text(errors="replace").splitlines():
        parts = line.split()
        if len(parts) < 5:
            continue
        try:
            step = int(parts[0])
            temp = float(parts[1])
            density = float(parts[4])
        except ValueError:
            continue
        if step < 0 or temp < 1.0 or temp > 5000.0:
            continue
        rows.append((step, temp, density))
    return rows[-n:]


def lammps_banner(log_path: Path) -> str:
    for line in log_path.read_text(errors="replace").splitlines():
        text = line.strip()
        if text.startswith("LAMMPS"):
            return text
    return "LAMMPS"


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--lmp", default=os.environ.get("LMP", "lmp"))
    parser.add_argument("--strides", default=DEFAULT_STRIDES)
    parser.add_argument("--transitions", type=int, default=200)
    parser.add_argument("--nmelt", type=int, default=10000)
    parser.add_argument("--ncool", type=int, default=5000)
    parser.add_argument("--nequil", type=int, default=10000)
    parser.add_argument("--seed", type=int, default=2522005)
    parser.add_argument("--temp", type=float, default=TEMPERATURE)
    parser.add_argument("--threads", type=int, default=os.cpu_count() or 1)
    args = parser.parse_args()
    if args.transitions < 1:
        raise SystemExit("need at least one transition")
    strides = [int(piece) for piece in args.strides.split(",") if piece.strip()]
    if not strides or any(stride < 1 for stride in strides):
        raise SystemExit(f"bad strides {args.strides}")
    lmp = args.lmp
    if shutil.which(lmp) is None and not Path(lmp).is_file():
        raise SystemExit(f"lmp not found: {lmp}")

    OUT.mkdir(parents=True, exist_ok=True)
    waters, box, n_ice, ice_box, height = build_slab(args.seed)
    write_data(OUT / "data.tip4p", waters, box)
    shutil.copy(DECK / "in.tip4p.equilibrate", OUT / "in.tip4p.equilibrate")
    shutil.copy(DECK / "in.tip4p.stride", OUT / "in.tip4p.stride")
    print(
        f"slab molecules {len(waters)} ice {n_ice} "
        f"box {box[0]:.3f} {box[1]:.3f} {box[2]:.3f} "
        f"ice_z {ice_box[2]:.3f} liquid_z {height:.3f}",
        flush=True,
    )

    run_lammps(
        lmp,
        "in.tip4p.equilibrate",
        {
            "nice": str(n_ice),
            "nmelt": str(args.nmelt),
            "ncool": str(args.ncool),
            "nequil": str(args.nequil),
            "seed": str(args.seed),
            "iseed": str(args.seed + 17),
            "temp": f"{args.temp:.2f}",
        },
        "log.tip4p-equil",
        args.threads,
    )
    banner = lammps_banner(OUT / "log.tip4p-equil")
    tail = thermo_tail(OUT / "log.tip4p-equil")
    if tail:
        _step, temp, density = tail[-1]
        print(
            f"equilibrated T {temp:.2f} K density {density:.4f} g/cm^3 "
            f"threads {args.threads}",
            flush=True,
        )
    rows = []
    for stride in strides:
        run_lammps(
            lmp,
            "in.tip4p.stride",
            {
                "stride": str(stride),
                "ntrans": str(args.transitions),
                "temp": f"{args.temp:.2f}",
            },
            f"log.tip4p-stride-{stride}",
            args.threads,
        )
        produced = OUT / "prod.lammpstrj"
        if not produced.is_file():
            raise SystemExit("LAMMPS did not write prod.lammpstrj")
        dest = OUT / f"tip4p-stride-{stride}.lammpstrj"
        produced.replace(dest)
        frames, atoms = count_frames(dest)
        size = dest.stat().st_size
        print(
            f"stride {stride} ps {stride_ps(stride):.4f} frames {frames} "
            f"oxygens {atoms} bytes {size} "
            f"path data/hop/tip4p-stride-{stride}.lammpstrj",
            flush=True,
        )
        if frames < args.transitions + 1:
            raise SystemExit(
                f"stride {stride} wrote {frames} frames, need {args.transitions + 1}"
            )
        if atoms != len(waters):
            raise SystemExit(f"stride {stride} has {atoms} atoms, expected {len(waters)} oxygens")
        rows.append((stride, frames, atoms, size))

    command = (
        "python3 input/hop/regenerate_tip4p.py "
        f"--strides {args.strides} --transitions {args.transitions} "
        f"--nmelt {args.nmelt} --ncool {args.ncool} --nequil {args.nequil} "
        f"--seed {args.seed} --temp {args.temp:.2f} --threads {args.threads}"
    )
    manifest = [
        f"command {command}",
        f"lammps {banner}",
        "model TIP4P/2005",
        "reference Abascal and Vega, J. Chem. Phys. 123, 234505 (2005)",
        f"molecules {len(waters)}",
        f"ice {n_ice}",
        f"liquid {len(waters) - n_ice}",
        f"replicate {NX} {NY} {NZ}",
        f"box_A {box[0]:.6f} {box[1]:.6f} {box[2]:.6f}",
        f"ice_z_A {ice_box[2]:.6f}",
        f"liquid_z_A {height:.6f}",
        f"gap_A {GAP:.1f}",
        f"seed {args.seed}",
        f"nmelt {args.nmelt}",
        f"ncool {args.ncool}",
        f"nequil {args.nequil}",
        f"temperature_K {args.temp:.2f}",
        "melt_temperature_K 300",
        "ensemble NPT",
        "pressure_atm 1",
        "thermostat Nose-Hoover",
        "barostat Nose-Hoover z",
        "tdamp_ps 0.2",
        "pdamp_ps 2",
        f"timestep_fs {TIMESTEP_FS:.1f}",
        "epsilon_kcal_mol 0.1852",
        "sigma_A 3.1589",
        "q_H 0.5564",
        "q_M -1.1128",
        "d_OM_A 0.1546",
        "r_OH_A 0.9572",
        "angle_deg 104.52",
        "lj_cutoff_A 10.0",
        "kspace pppm/tip4p 1.0e-4",
        "dump oxygen",
        "oxygen_type 1",
        f"threads {args.threads}",
    ]
    if tail:
        _step, temp, density = tail[-1]
        manifest.append(f"equilibrated_temperature_K {temp:.4f}")
        manifest.append(f"equilibrated_density_g_cm3 {density:.6f}")
    for stride, frames, atoms, size in rows:
        manifest.append(
            f"trajectory stride {stride} stride_ps {stride_ps(stride):.4f} "
            f"frames {frames} transitions {frames - 1} oxygens {atoms} "
            f"bytes {size} data/hop/tip4p-stride-{stride}.lammpstrj"
        )
    (OUT / "tip4p-manifest.txt").write_text("\n".join(manifest) + "\n")
    print(command)


if __name__ == "__main__":
    main()
