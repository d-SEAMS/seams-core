# Copyright (c) d-SEAMS developers
"""Proton-ordered TIP4P/Ice ice Ic slab for the engine end-to-end runs.

Oxygens sit on a diamond lattice. Sublattice A donates along (1,1,1) and
(1,-1,-1), sublattice B along (1,-1,1) and (1,1,-1), so every O-O bond
carries exactly one proton and the Bernal-Fowler rules hold without a
proton-disorder search. The slab fills the lower half of the box; the
engine's own solvent tool fills the upper half with liquid.
"""
import argparse
import math

A_IC = 6.416  # Angstrom, ice Ic at the TIP4P/Ice density (0.906 g/cm3)
D_OH = 0.9572
D_OM = 0.1577
HALF_HOH = math.radians(104.52 / 2)

_BASIS_A = [(0, 0, 0), (0, 0.5, 0.5), (0.5, 0, 0.5), (0.5, 0.5, 0)]
_DONATE = {
    "A": [(1, 1, 1), (1, -1, -1)],
    "B": [(1, -1, 1), (1, 1, -1)],
}


def _unit(v):
    n = math.sqrt(sum(c * c for c in v))
    return tuple(c / n for c in v)


def _water(o, donors):
    u1, u2 = (_unit(d) for d in donors)
    b = _unit(tuple(x + y for x, y in zip(u1, u2)))
    p = _unit(tuple(x - y for x, y in zip(u1, u2)))
    c, s = math.cos(HALF_HOH), math.sin(HALF_HOH)
    h1 = tuple(oc + D_OH * (c * bc + s * pc) for oc, bc, pc in zip(o, b, p))
    h2 = tuple(oc + D_OH * (c * bc - s * pc) for oc, bc, pc in zip(o, b, p))
    m = tuple(oc + D_OM * bc for oc, bc in zip(o, b))
    return o, h1, h2, m


def slab(nx, ny, nz, a=A_IC):
    """Return (molecules, box) in Angstrom; the box is twice the slab in z."""
    shift = a / 8
    mols = []
    for i in range(nx):
        for j in range(ny):
            for k in range(nz):
                for f in _BASIS_A:
                    for lat, off in (("A", 0.0), ("B", 0.25)):
                        o = tuple((c + idx + off) * a + shift
                                  for c, idx in zip(f, (i, j, k)))
                        mols.append(_water(o, _DONATE[lat]))
    return mols, (nx * a, ny * a, 2 * nz * a)


def write_gro(path, mols, box):
    names = ("OW", "HW1", "HW2", "MW")
    with open(path, "w") as fh:
        fh.write("TIP4P/Ice ice Ic slab\n%d\n" % (4 * len(mols)))
        serial = 0
        for r, mol in enumerate(mols, 1):
            for name, x in zip(names, mol):
                serial += 1
                fh.write("%5d%-5s%5s%5d%8.3f%8.3f%8.3f\n" % (
                    r % 100000, "SOL", name, serial % 100000,
                    x[0] / 10, x[1] / 10, x[2] / 10))
        fh.write("%10.5f%10.5f%10.5f\n" % tuple(b / 10 for b in box))


def write_pdb(path, mols, box):
    names = ("O", "H1", "H2", "M")
    elements = ("O", "H", "H", "")
    with open(path, "w") as fh:
        fh.write("CRYST1%9.3f%9.3f%9.3f  90.00  90.00  90.00 P 1           1\n" % box)
        serial = 0
        for r, mol in enumerate(mols, 1):
            for name, el, x in zip(names, elements, mol):
                serial += 1
                fh.write("HETATM%5d %-4s HOH A%4d    %8.3f%8.3f%8.3f  1.00  0.00          %2s\n" % (
                    serial % 100000, name, r % 10000, x[0], x[1], x[2], el))
        fh.write("TER\nEND\n")


def main():
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("out")
    ap.add_argument("--cells", type=int, nargs=3, default=(6, 6, 3))
    ap.add_argument("--format", choices=("gro", "pdb"), default="gro")
    args = ap.parse_args()
    mols, box = slab(*args.cells)
    (write_gro if args.format == "gro" else write_pdb)(args.out, mols, box)


if __name__ == "__main__":
    main()
