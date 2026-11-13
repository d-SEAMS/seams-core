#!/usr/bin/env python3
"""Build a direct-coexistence ice/brine box for LAMMPS.

TIP4P/2005 water with the Madrid 2019 scaled-charge ions (Zeron, Abascal
and Vega, J. Chem. Phys. 151, 134504 (2019); topology files by
Fernandez Sedano, Zenodo 10.5281/zenodo.18685416, CC-BY-4.0). GenIce
generates a hydrogen-disordered ice Ih block with TIP4P geometry; the upper
half of the block along z is marked liquid and NaCl pairs replace water
molecules there. The LAMMPS input melts the liquid half with the ice half
held, then lets the front move at the run temperature.

Atom types: 1 O (TIP4P/2005), 2 H, 3 Na, 4 Cl. atom_style full, one bond
type (O-H) and one angle type (H-O-H) for SHAKE; the M site is implicit in
pair_style lj/cut/tip4p/long. Molecule ids: one per water, one per ion.

Charges: O -1.1128, H 0.5564, Na 0.85, Cl -0.85.
"""

from __future__ import annotations

import argparse
import math
import subprocess

import numpy as np

Q_O, Q_H, Q_NA, Q_CL = -1.1128, 0.5564, 0.85, -0.85
M_O, M_H, M_NA, M_CL = 15.9994, 1.00794, 22.9898, 35.453
GENICE_KIND = "1h"


def genice_tip4p(exe, rep):
    cmd = [exe, GENICE_KIND, "--rep", *map(str, rep), "--water", "tip4p", "--format", "gromacs"]
    out = subprocess.run(cmd, check=True, capture_output=True, text=True).stdout.splitlines()
    n = int(out[1])
    atoms = []
    for line in out[2 : 2 + n]:
        name = line[10:15].strip()
        resid = int(line[0:5])
        xyz = np.array([float(line[20:28]), float(line[28:36]), float(line[36:44])]) * 10.0
        atoms.append((resid, name, xyz))
    cell = np.array([float(x) for x in out[2 + n].split()[:3]]) * 10.0
    # group by residue: O first, then the two H (skip the M site)
    mols = {}
    for resid, name, xyz in atoms:
        mols.setdefault(resid, {})
        if name.startswith("O"):
            mols[resid]["O"] = xyz
        elif name.startswith("H"):
            mols[resid].setdefault("H", []).append(xyz)
    water = [(m["O"], m["H"][0], m["H"][1]) for m in mols.values() if "O" in m and len(m.get("H", [])) == 2]
    return water, cell


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--genice", default="genice2")
    ap.add_argument("--rep", nargs=3, type=int, default=[4, 4, 6])
    ap.add_argument("--pairs", type=int, default=20, help="NaCl pairs in the liquid half")
    ap.add_argument("--seed", type=int, default=7)
    ap.add_argument("--data", required=True)
    ap.add_argument("--groups", required=True, help="LAMMPS include with ice/liquid groups")
    args = ap.parse_args()

    water, cell = genice_tip4p(args.genice, args.rep)
    rng = np.random.default_rng(args.seed)
    z = np.array([o[2] for o, _, _ in water])
    liquid = z >= 0.5 * cell[2]
    liquid_idx = np.nonzero(liquid)[0]
    if 2 * args.pairs > len(liquid_idx):
        raise SystemExit("more ions than liquid molecules")
    ion_sites = rng.choice(liquid_idx, size=2 * args.pairs, replace=False)
    na_sites = set(ion_sites[: args.pairs].tolist())
    cl_sites = set(ion_sites[args.pairs :].tolist())

    lines_atoms, lines_bonds, lines_angles = [], [], []
    aid = 0
    mol = 0
    ice_ids, liq_ids, ion_ids = [], [], []
    for i, (o, h1, h2) in enumerate(water):
        mol += 1
        if i in na_sites or i in cl_sites:
            aid += 1
            t, q = (3, Q_NA) if i in na_sites else (4, Q_CL)
            lines_atoms.append(f"{aid} {mol} {t} {q:.4f} {o[0]:.6f} {o[1]:.6f} {o[2]:.6f}")
            ion_ids.append(aid)
            continue
        oid = aid + 1
        for t, q, xyz in ((1, Q_O, o), (2, Q_H, h1), (2, Q_H, h2)):
            aid += 1
            lines_atoms.append(f"{aid} {mol} {t} {q:.4f} {xyz[0]:.6f} {xyz[1]:.6f} {xyz[2]:.6f}")
        lines_bonds.append(f"{len(lines_bonds) + 1} 1 {oid} {oid + 1}")
        lines_bonds.append(f"{len(lines_bonds) + 1} 1 {oid} {oid + 2}")
        lines_angles.append(f"{len(lines_angles) + 1} 1 {oid + 1} {oid} {oid + 2}")
        (liq_ids if liquid[i] else ice_ids).extend((oid, oid + 1, oid + 2))

    with open(args.data, "w") as fh:
        fh.write("ice Ih slab and brine, TIP4P/2005 with Madrid 2019 NaCl (charges scaled 0.85)\n\n")
        fh.write(f"{aid} atoms\n{len(lines_bonds)} bonds\n{len(lines_angles)} angles\n")
        fh.write("4 atom types\n1 bond types\n1 angle types\n\n")
        fh.write(f"0.0 {cell[0]:.6f} xlo xhi\n0.0 {cell[1]:.6f} ylo yhi\n0.0 {cell[2]:.6f} zlo zhi\n\n")
        fh.write(f"Masses\n\n1 {M_O}\n2 {M_H}\n3 {M_NA}\n4 {M_CL}\n\n")
        fh.write("Atoms # full\n\n" + "\n".join(lines_atoms) + "\n\n")
        fh.write("Bonds\n\n" + "\n".join(lines_bonds) + "\n\n")
        fh.write("Angles\n\n" + "\n".join(lines_angles) + "\n")

    def ranges(ids):
        ids = sorted(ids)
        out, start, prev = [], ids[0], ids[0]
        for x in ids[1:]:
            if x != prev + 1:
                out.append(f"{start}:{prev}" if start != prev else f"{start}")
                start = x
            prev = x
        out.append(f"{start}:{prev}" if start != prev else f"{start}")
        return out

    def group_lines(name, ids):
        chunks = ranges(ids)
        lines = []
        for k in range(0, len(chunks), 40):
            lines.append(f"group {name} id " + " ".join(chunks[k : k + 40]))
        return lines

    with open(args.groups, "w") as fh:
        fh.write("# groups written by brine_system.py\n")
        fh.write("\n".join(group_lines("ice", ice_ids)) + "\n")
        fh.write("\n".join(group_lines("liquid", liq_ids + ion_ids)) + "\n")
        fh.write("group ions type 3 4\ngroup oxygens type 1\ngroup water type 1 2\n")
        fh.write(f"variable nwater equal {len(water) - 2 * args.pairs}\n")
        fh.write(f"variable nions equal {2 * args.pairs}\n")
    n_w = len(water) - 2 * args.pairs
    molality = args.pairs / (n_w * M_O / 1000.0 + n_w * 2 * M_H / 1000.0)
    print(f"water={n_w} pairs={args.pairs} molality={molality:.3f} mol/kg cell={cell.round(3).tolist()} "
          f"ice_atoms={len(ice_ids)} liquid_atoms={len(liq_ids) + len(ion_ids)}")


if __name__ == "__main__":
    main()
