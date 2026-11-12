# Copyright (c) d-SEAMS developers
"""OpenMM reporter. Positions are read in Angstrom and passed to the C API."""

from __future__ import annotations

import time

import numpy as np
import openmm.unit as unit

import dseams


class DseamsReporter:
    def __init__(self, reportInterval, cutoff=3.5, rule=dseams.PUBLISHED):
        self.interval = int(reportInterval)
        self.cutoff = float(cutoff)
        self.rule = rule
        self.calls = 0
        self.ice = 0
        self.atoms = 0
        self.seconds = 0.0

    def describeNextReport(self, simulation):
        steps = self.interval - simulation.currentStep % self.interval
        return {"steps": steps, "periodic": True, "include": ["positions"]}

    def report(self, simulation, state):
        pos = state.getPositions(asNumpy=True).value_in_unit(unit.angstrom)
        box = state.getPeriodicBoxVectors(asNumpy=True).value_in_unit(unit.angstrom)
        cell = np.zeros(9, dtype=np.float64)
        cell[3] = box[0, 0]
        cell[4] = box[1, 1]
        cell[5] = box[2, 2]
        cell[6] = box[1, 0]
        cell[7] = box[2, 0]
        cell[8] = box[2, 1]
        t0 = time.perf_counter()
        labels, counts = dseams.chill_plus(pos, cell, rule=self.rule, cutoff=self.cutoff)
        self.seconds += time.perf_counter() - t0
        self.calls += 1
        self.ice += int(counts[0])
        self.atoms += int(labels.shape[0])

    def overhead(self):
        return self.seconds / self.calls if self.calls else 0.0
