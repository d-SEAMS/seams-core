# Copyright (c) d-SEAMS developers
"""Label one cubic-ice cell and print the ice count."""

import math

import numpy as np

from dseams import PUBLISHED, OVITO, chill_plus, version


def diamond(nrep=2):
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
    return np.asarray(points), (nrep * a, nrep * a, nrep * a)


def main():
    pos, box = diamond()
    print("dseams", version())
    for rule in (PUBLISHED, OVITO):
        labels, counts = chill_plus(pos, box, rule=rule)
        print(rule, int(counts[0]), labels.shape[0])


if __name__ == "__main__":
    main()
