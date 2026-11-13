"""Rings, cages, signature cages, and symmetry agree on one ideal lattice."""

import math
import os
import sys

sys.path.insert(0, os.path.join(os.path.dirname(__file__), "..", "bindings", "python"))

try:
    import numpy as np
except ImportError:
    sys.stderr.write("numpy is not installed\n")
    raise SystemExit(77)

from frames import cages, rings, signature_cages, symmetry


def diamond(nrep=2):
    oo = 2.75
    a = oo * 4.0 / math.sqrt(3.0)
    basis = ((0, 0, 0), (0.25, 0.25, 0.25))
    fcc = ((0, 0, 0), (0.5, 0.5, 0), (0.5, 0, 0.5), (0, 0.5, 0.5))
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
    return np.asarray(points, dtype=np.float64), (nrep * a, nrep * a, nrep * a)


def main():
    pos, box = diamond()
    try:
        hist = rings(pos, box)
    except OSError as exc:
        # The interpreter and libseams_lang were not built for the same
        # runtime (a conda libstdc++ or a different libpython). Skip.
        sys.stderr.write(f"python/library mismatch: {exc}\n")
        raise SystemExit(77) from exc
    if hist.shape != (7,) or int(hist[6]) <= 0:
        raise SystemExit("ideal lattice has no six-rings")
    per, counts = cages(pos, box)
    if per.shape != (pos.shape[0],) or sum(counts) != pos.shape[0]:
        raise SystemExit("cage populations do not cover the lattice")
    n_cages, n_atoms = signature_cages(pos, box, spec="512")
    if n_cages < 0 or n_atoms < 0:
        raise SystemExit("signature census is negative")
    labels = symmetry(pos, box)
    if labels.shape != (pos.shape[0],):
        raise SystemExit("symmetry labels do not match the lattice")


if __name__ == "__main__":
    try:
        main()
    except OSError as exc:
        sys.stderr.write(f"python/library mismatch: {exc}\n")
        raise SystemExit(77) from exc
