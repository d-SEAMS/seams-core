# Copyright (c) d-SEAMS developers
"""OpenMM reporter. Positions come from the context, on the host."""


class DseamsReporter:
    """OpenMM reporter that calls the in-situ engine.

    Positions are read from the context on the host. The reporter
    writes no trajectory of its own.
    """

    def __init__(self, engine, dump, interval=1, types=None):
        self.engine = engine
        self.dump = dump
        self.interval = int(interval)
        self.types = types

    def describeNextReport(self, simulation):
        step = int(simulation.currentStep)
        remain = self.interval - (step % self.interval)
        return {"steps": remain, "periodic": None, "include": ["positions"]}

    def report(self, simulation, state):
        del state
        import time

        from openmm.unit import angstrom

        t0 = time.perf_counter()
        pulled = simulation.context.getState(getPositions=True, enforcePeriodicBox=True)
        pos = pulled.getPositions(asNumpy=True).value_in_unit(angstrom)
        vecs = pulled.getPeriodicBoxVectors(asNumpy=True).value_in_unit(angstrom)
        box = self.engine.box_from_rows(vecs)
        copy_s = time.perf_counter() - t0
        self.engine.observe(
            int(simulation.currentStep), pos, box, None, types=self.types, copy_s=copy_s
        )
