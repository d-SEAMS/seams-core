"""The ctypes wrapper in bindings/python/chill_plus.py labels an ideal lattice."""

import math
import os
import sys

sys.path.insert(0, os.path.join(os.path.dirname(__file__), "..", "bindings", "python"))

try:
    import numpy as np
except ImportError:
    sys.stderr.write("numpy is not installed\n")
    raise SystemExit(77)

from chill_plus import chill_plus


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
    cubic = 0
    for rule in ("chill+", "ovito-chill+"):
        labels, counts = chill_plus(pos, box, rule=rule)
        if labels.shape != (pos.shape[0],) or int(counts.sum()) != pos.shape[0]:
            raise SystemExit("bad shape for " + rule)
        if not np.all(labels == cubic):
            raise SystemExit(rule + " is not cubic on the ideal lattice")
    text = 'rule = "ovito-chill+"\n'
    labels, _ = chill_plus(pos, box, rule=text)
    if not np.all(labels == cubic):
        raise SystemExit("document text did not select ovito-chill+")
    try:
        chill_plus(pos, box, rule="not-a-rule")
    except RuntimeError:
        return
    raise SystemExit("an unknown rule was accepted")


if __name__ == "__main__":
    main()
