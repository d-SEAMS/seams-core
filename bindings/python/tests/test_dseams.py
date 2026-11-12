# Copyright (c) d-SEAMS developers
"""The C API, reached through nanobind, on one cubic cell."""

import math
import unittest

import numpy as np

try:
    from ase import Atoms
except ImportError:
    Atoms = None

import dseams
from dseams import CUTOFF, OVITO, PUBLISHED


def diamond(nrep=2):
    a = 2.75 * 4.0 / math.sqrt(3.0)
    fcc = ((0, 0, 0), (0.5, 0.5, 0), (0.5, 0, 0.5), (0, 0.5, 0.5))
    basis = ((0, 0, 0), (0.25, 0.25, 0.25))
    points = [
        (
            (i + shift[0] + site[0]) * a,
            (j + shift[1] + site[1]) * a,
            (k + shift[2] + site[2]) * a,
        )
        for i in range(nrep)
        for j in range(nrep)
        for k in range(nrep)
        for shift in fcc
        for site in basis
    ]
    length = nrep * a
    return np.asarray(points, dtype=np.float64), length


class ApiTest(unittest.TestCase):
    def test_both_rules_and_order(self):
        pos, length = diamond()
        dseams.set_num_threads(1)
        self.assertEqual(dseams.num_threads(), 1)
        for rule in (PUBLISHED, OVITO):
            labels, counts = dseams.chill_plus(pos, (length, length, length), rule=rule)
            self.assertEqual(labels.shape, (pos.shape[0],))
            self.assertTrue(np.all(labels == 0))
            self.assertEqual(int(counts[0]), pos.shape[0])
        q4, q6 = dseams.q4q6(pos, (length, length, length))
        self.assertGreater(float(q6.mean()), 0.3)
        self.assertTrue(np.all(np.isfinite(q4)))
        offsets, cols = dseams.neighbours(pos, (length, length, length), graph=CUTOFF)
        degree = np.diff(offsets)
        self.assertTrue(np.all(degree == 4))
        self.assertEqual(cols.shape[0], int(degree.sum()))
        sizes, sym, hist, fraction = dseams.rings(pos, (length, length, length))
        self.assertGreater(int(hist[6]), 0)
        self.assertEqual(int(np.sum(sizes == 6)), int(hist[6]))
        self.assertGreater(fraction, 0.5)
        self.assertEqual(sym.shape, sizes.shape)

    def test_ase_and_mdanalysis(self):
        pos, length = diamond()
        if Atoms is None:
            self.skipTest("ASE is not installed")
        atoms = Atoms(symbols=["O"] * len(pos), positions=pos, cell=[length, length, length], pbc=True)
        apos, box = dseams.from_ase(atoms)
        labels, _ = dseams.chill_plus(apos, box)
        self.assertTrue(np.all(labels == 0))
        try:
            import MDAnalysis as mda
        except ImportError:
            self.skipTest("MDAnalysis is not installed")

        u = mda.Universe.empty(len(pos), trajectory=True)
        u.atoms.positions = pos
        u.dimensions = [length, length, length, 90.0, 90.0, 90.0]
        mpos, mbox = dseams.from_mdanalysis(u)
        labels, _ = dseams.chill_plus(mpos, mbox)
        self.assertTrue(np.all(labels == 0))

    def test_triclinic_box_is_accepted(self):
        pos, length = diamond()
        box = np.array([0, 0, 0, length, length, length, 0.3, 0.0, 0.0])
        labels, _ = dseams.chill_plus(pos, box, rule=PUBLISHED)
        self.assertEqual(labels.shape, (len(pos),))


if __name__ == "__main__":
    unittest.main()
