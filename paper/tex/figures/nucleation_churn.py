"""Bond churn between consecutive frames of a LAMMPS dump on the 3.5 A graph.

For each pair of consecutive frames, counts the bonds of the cutoff graph
that appear or disappear and the atoms whose neighbour row changes, and
tests them against the two full-pass rules of RingUpdater (changed atoms
above n/5, changed bonds above 0.15 per molecule). Quartiles interpolate
linearly between order statistics. numpy only.

Usage: python3 nucleation_churn.py DUMP > data/nucleation-churn.txt
"""

import sys

import numpy as np

CUTOFF_A = 3.5
ATOM_FRACTION = 0.2
BOND_CROSSOVER = 0.15
# Lines of a LAMMPS dump frame after its "ITEM: TIMESTEP" line: the atom count,
# the first of three box bounds, the atom column names and the first atom.
DUMP_COUNT_LINE = 3
DUMP_BOX_LINE = 5
DUMP_COLUMNS_LINE = 8
DUMP_FIRST_ATOM_LINE = 9
# A frame gives one box bound line and one position column per axis.
AXES = ("x", "y", "z")
# Atoms whose distances to every atom are computed in one array, which bounds
# the memory of that array.
CHUNK_ATOMS = 256
# The churn rows give the first and third quartiles of the transitions.
QUARTILES = (0.25, 0.75)


def frames(path):
    lines = open(path).read().splitlines()
    index = 0
    while index < len(lines):
        count = int(lines[index + DUMP_COUNT_LINE])
        bounds = [lines[index + DUMP_BOX_LINE + k].split() for k in range(len(AXES))]
        low = np.array([float(row[0]) for row in bounds])
        high = np.array([float(row[1]) for row in bounds])
        columns = lines[index + DUMP_COLUMNS_LINE].split()[2:]
        xyz = [columns.index(axis) for axis in AXES]
        ids = columns.index("id")
        rows = [lines[index + DUMP_FIRST_ATOM_LINE + k].split() for k in range(count)]
        rows.sort(key=lambda row: int(row[ids]))
        positions = np.array([[float(row[c]) for c in xyz] for row in rows])
        yield positions - low, high - low
        index += DUMP_FIRST_ATOM_LINE + count


def bonds(positions, box):
    found = set()
    for start in range(0, len(positions), CHUNK_ATOMS):
        delta = positions[start:start + CHUNK_ATOMS, None, :] - positions[None, :, :]
        delta -= box * np.round(delta / box)
        close = (delta**2).sum(axis=2) <= CUTOFF_A**2
        for i, j in zip(*np.nonzero(close)):
            if start + i < j:
                found.add((start + i, j))
    return found


def main():
    previous = None
    rows = []
    for positions, box in frames(sys.argv[1]):
        current = bonds(positions, box)
        if previous is not None:
            changed = previous ^ current
            atoms = {a for pair in changed for a in pair}
            n = len(positions)
            rows.append((len(changed), len(atoms), len(changed) / n,
                         len(atoms) > int(ATOM_FRACTION * n),
                         len(changed) / n > BOND_CROSSOVER))
        previous = current
    per_bond = sorted(r[2] for r in rows)
    changed_atoms = sorted(r[1] for r in rows)
    bond_q1, bond_q3 = np.quantile(per_bond, QUARTILES)
    atom_q1, atom_q3 = np.quantile(changed_atoms, QUARTILES)
    print("transitions", len(rows))
    print("molecules", len(positions))
    print("median_bonds_per_molecule", f"{per_bond[len(per_bond) // 2]:.3f}")
    print("q1_bonds_per_molecule", f"{bond_q1:.3f}")
    print("q3_bonds_per_molecule", f"{bond_q3:.3f}")
    print("min_bonds_per_molecule", f"{per_bond[0]:.3f}")
    print("max_bonds_per_molecule", f"{per_bond[-1]:.3f}")
    print("median_changed_atoms", changed_atoms[len(rows) // 2])
    print("q1_changed_atoms", f"{atom_q1:g}")
    print("q3_changed_atoms", f"{atom_q3:g}")
    print("over_atom_rule", sum(r[3] for r in rows))
    print("over_bond_rule", sum(r[4] for r in rows))
    print("over_either_rule", sum(r[3] or r[4] for r in rows))


if __name__ == "__main__":
    main()
