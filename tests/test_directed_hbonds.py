# Copyright (c) d-SEAMS developers
"""Hand-checked directed hydrogen bonds through the Python constructor."""

import math
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "bindings" / "python"))

import seams_hbonds


def fail(message):
    print(message, file=sys.stderr)
    raise SystemExit(1)


def close(got, want, tol, label):
    if abs(got - want) > tol:
        fail("{} got {} want {}".format(label, got, want))


def angle(a, b):
    dot = a[0] * b[0] + a[1] * b[1] + a[2] * b[2]
    na = math.sqrt(a[0] * a[0] + a[1] * a[1] + a[2] * a[2])
    nb = math.sqrt(b[0] * b[0] + b[1] * b[1] + b[2] * b[2])
    cosine = max(-1.0, min(1.0, dot / (na * nb)))
    return math.degrees(math.acos(cosine))


def triclinic(xi, yi, zi, xj, yj, zj):
    lx = 10.0
    ly = 8.660254037844386
    lz = 10.0
    xy = 5.0
    sy_i = yi / ly
    sy_j = yj / ly
    dsx = (xi - xy * sy_i) / lx - (xj - xy * sy_j) / lx
    dsy = sy_i - sy_j
    dsz = zi / lz - zj / lz
    dsx -= round(dsx)
    dsy -= round(dsy)
    dsz -= round(dsz)
    return (lx * dsx + xy * dsy, ly * dsy, lz * dsz)


def linear():
    frame = seams_hbonds.Frame.from_clouds(
        [0, 0, 0, 2.8, 0, 0],
        [1, 0, 0],
        [0, 0, 0, 20, 20, 20, 0, 0, 0],
        [8, 8],
        [1, 2],
        [2],
        [1],
    )
    if frame.count() != 3:
        fail("linear count")
    oo = (2.8, 0.0, 0.0)
    oh = (1.8, 0.0, 0.0)
    close(angle(oo, oh), 0.0, 1e-12, "linear oracle")
    edges = frame.directed_hbonds([8], [8])
    if len(edges) != 1 or (edges[0].donor, edges[0].hydrogen, edges[0].acceptor) != (0, 2, 1):
        fail("linear edge {}".format(edges))
    close(edges[0].d_DA, 2.8, 1e-12, "linear d_DA")
    close(edges[0].d_AH, 1.8, 1e-12, "linear d_AH")
    close(edges[0].angle, 0.0, 1e-8, "linear angle")
    if frame.directed_hbonds([], [8]):
        fail("empty donor set")
    frame.close()


def triclinic_bond():
    ly = 8.660254037844386
    frame = seams_hbonds.Frame.from_clouds(
        [5.5, 6.360254037844386, 1.0, 0.5, 0.5, 1.0],
        [0.5, 0.1, 1.0],
        [0, 0, 0, 10, ly, 10, 5, 0, 0],
        [1, 1],
        [1, 2],
        [2],
        [1],
    )
    oo = triclinic(0.5, 0.5, 1.0, 5.5, 6.360254037844386, 1.0)
    oh = triclinic(0.5, 0.5, 1.0, 0.5, 0.1, 1.0)
    close(oo[1], 2.8, 1e-10, "triclinic oo")
    close(math.sqrt(oh[0] ** 2 + oh[1] ** 2 + oh[2] ** 2), 0.4, 1e-10, "triclinic oh")
    expect = angle(oo, oh)
    edges = frame.directed_hbonds([1], [1])
    if len(edges) != 1:
        fail("triclinic edges {}".format(len(edges)))
    close(edges[0].d_DA, 2.8, 1e-10, "triclinic d_DA")
    close(edges[0].d_AH, 0.4, 1e-10, "triclinic d_AH")
    close(edges[0].angle, expect, 1e-8, "triclinic angle")
    frame.close()


def framework():
    frame = seams_hbonds.Frame.from_clouds(
        [5, 5, 5, 5, 7.8, 5, 7.8, 5, 5, 5, 5, 7.8],
        [6, 5, 5, 5, 6, 5, 5, 6.8, 5, 5, 5, 6.8],
        [0, 0, 0, 30, 30, 30, 0, 0, 0],
        [8, 8, 10, 11],
        [1, 2, 50, 70],
        [2, 2, 2, 2],
        [1, 1, 2, 70],
    )
    water = frame.directed_hbonds([8], [8])
    if [(e.donor, e.hydrogen, e.acceptor) for e in water] != [(0, 5, 1), (1, 6, 0)]:
        fail("water edges {}".format(water))
    carbonyl = frame.directed_hbonds([8], [10])
    if len(carbonyl) != 1 or (carbonyl[0].donor, carbonyl[0].hydrogen, carbonyl[0].acceptor) != (0, 4, 2):
        fail("carbonyl edge {}".format(carbonyl))
    close(carbonyl[0].angle, 0.0, 1e-8, "carbonyl angle")
    amine = frame.directed_hbonds([11], [8])
    if len(amine) != 1 or (amine[0].donor, amine[0].acceptor) != (3, 0):
        fail("amine edge {}".format(amine))
    if frame.directed_hbonds([8], [99]):
        fail("wrong acceptor type")
    right = seams_hbonds.Frame.from_clouds(
        [0, 0, 0, 2.5, 0, 0],
        [2.5, 0.8, 0],
        [0, 0, 0, 20, 20, 20, 0, 0, 0],
        [1, 1],
        [1, 2],
        [2],
        [1],
    )
    if abs(angle((2.5, 0, 0), (0, -0.8, 0)) - 90.0) > 1e-8:
        fail("right-angle oracle")
    if right.directed_hbonds([1], [1]):
        fail("right angle was kept")
    right.close()
    same = seams_hbonds.Frame.from_clouds(
        [0, 0, 0, 2.8, 0, 0],
        [1, 0, 0],
        [0, 0, 0, 20, 20, 20, 0, 0, 0],
        [1, 1],
        [7, 7],
        [2],
        [7],
    )
    if same.directed_hbonds([1], [1], exclude_same_mol=True):
        fail("same molecule was kept")
    if len(same.directed_hbonds([1], [1], exclude_same_mol=False)) != 1:
        fail("same molecule was dropped")
    same.close()
    frame.close()


def short_metadata_is_rejected():
    try:
        seams_hbonds.Frame.from_clouds(
            [0, 0, 0, 2.8, 0, 0],
            [1, 0, 0],
            [0, 0, 0, 20, 20, 20, 0, 0, 0],
            [8],
            [1, 2],
            [2],
            [1],
        )
    except seams_hbonds.SeamsError as exc:
        if exc.status != 4:
            fail("short heavy_types status {}".format(exc.status))
    else:
        fail("short heavy_types was accepted")


def main():
    linear()
    triclinic_bond()
    framework()
    short_metadata_is_rejected()
    print("directed hbonds ok")


if __name__ == "__main__":
    main()
