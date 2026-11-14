# Copyright (c) d-SEAMS developers
"""Two steps of a frozen cubic-ice cell through the one OpenMM reporter."""

import math
import sys

import numpy as np
import openmm as mm
import openmm.app as app
import openmm.unit as unit

from dseams_insitu.engine import InSitu
from dseams_insitu.openmm import DseamsReporter


def diamond(nrep=1):
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
    return np.asarray(points, dtype=np.float64), nrep * a


def main():
    pos_a, length = diamond(1)
    n = len(pos_a)
    system = mm.System()
    for _ in range(n):
        system.addParticle(16.0 * unit.dalton)
    box = length * unit.angstrom
    system.setDefaultPeriodicBoxVectors((box, 0, 0), (0, box, 0), (0, 0, box))
    force = mm.CustomExternalForce("0")
    for i in range(n):
        force.addParticle(i, [])
    system.addForce(force)
    topology = app.Topology()
    chain = topology.addChain()
    residue = topology.addResidue("ICE", chain)
    element = app.element.oxygen
    for _ in range(n):
        topology.addAtom("O", element, residue)
    integrator = mm.VerletIntegrator(0.001 * unit.picoseconds)
    platform = mm.Platform.getPlatformByName("Reference")
    simulation = app.Simulation(topology, system, integrator, platform)
    simulation.context.setPositions(pos_a * unit.angstrom)
    engine = InSitu()
    reporter = DseamsReporter(engine, interval=1)
    simulation.reporters.append(reporter)
    simulation.step(2)
    if len(engine.frames) < 2:
        print("openmm frames", len(engine.frames), file=sys.stderr)
        return 1
    for frame in engine.frames:
        labels = frame["labels"]
        if labels.shape != (n,) or not np.all(labels == 0):
            print("openmm labels", labels, file=sys.stderr)
            return 1
    print(f"openmm calls {len(engine.frames)} ice {n * len(engine.frames)}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
