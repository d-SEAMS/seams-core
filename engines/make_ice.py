# Copyright (c) d-SEAMS developers
"""Write the small ice decks the engine tests read. No trajectory archives."""

import math
from pathlib import Path

ROOT = Path(__file__).resolve().parent


def diamond(nrep):
    a = 2.75 * 4.0 / math.sqrt(3.0)
    fcc = ((0, 0, 0), (0.5, 0.5, 0), (0.5, 0, 0.5), (0, 0.5, 0.5))
    basis = ((0, 0, 0), (0.25, 0.25, 0.25))
    points = []
    for i in range(nrep):
        for j in range(nrep):
            for k in range(nrep):
                for shift in fcc:
                    for site in basis:
                        points.append(
                            (
                                (i + shift[0] + site[0]) * a,
                                (j + shift[1] + site[1]) * a,
                                (k + shift[2] + site[2]) * a,
                            )
                        )
    return points, nrep * a


def write_lammps(path, points, length):
    lines = [
        "cubic ice",
        "",
        f"{len(points)} atoms",
        "1 atom types",
        "",
        f"0.0 {length:.8f} xlo xhi",
        f"0.0 {length:.8f} ylo yhi",
        f"0.0 {length:.8f} zlo zhi",
        "",
        "Masses",
        "",
        "1 16.0",
        "",
        "Atoms # atomic",
        "",
    ]
    for i, (x, y, z) in enumerate(points, start=1):
        lines.append(f"{i} 1 {x:.8f} {y:.8f} {z:.8f}")
    lines.append("")
    path.write_text("\n".join(lines))


def write_gro(path, points, length):
    chunks = []
    for frame, shift in enumerate((0.0, 0.01)):
        rows = [f"cubic ice frame {frame}", f"{len(points):5d}"]
        for i, (x, y, z) in enumerate(points, start=1):
            rows.append(
                f"{1:5d}ICE     O{i:5d}{(x + shift) / 10:8.4f}{(y) / 10:8.4f}{(z) / 10:8.4f}"
            )
        rows.append(f"{length / 10:10.5f}{length / 10:10.5f}{length / 10:10.5f}")
        chunks.append("\n".join(rows))
    path.write_text("\n".join(chunks) + "\n")


def write_xyz(path, points, length):
    # plumed driver reads xyz as nm and hands the CV Angstrom under UNITS LENGTH=A.
    rows = [str(len(points)),
            f"{length / 10:.8f} 0 0 0 {length / 10:.8f} 0 0 0 {length / 10:.8f}"]
    for x, y, z in points:
        rows.append(f"O {x / 10:.8f} {y / 10:.8f} {z / 10:.8f}")
    path.write_text("\n".join(rows) + "\n")


def main():
    big, length = diamond(2)
    small, _ = diamond(1)
    write_lammps(ROOT / "lammps" / "ice.data", big, length)
    write_gro(ROOT / "gromacs" / "ice.gro", big, length)
    write_xyz(ROOT / "plumed" / "ice.xyz", big, length)


if __name__ == "__main__":
    main()
