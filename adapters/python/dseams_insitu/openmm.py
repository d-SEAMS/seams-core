# Copyright (c) d-SEAMS developers
"""OpenMM reporter. Positions come from the context, on the host."""

import time

import numpy as np


class DseamsReporter:
    """OpenMM reporter that calls the in-situ engine.

    Positions are read from the context on the host. With dump set, each
    report appends one LAMMPS dump frame holding the coordinates handed to
    d-SEAMS, in angstrom, and the CHILL+ label of each site. Atom ids are
    1-based OpenMM particle indices. seconds and calls time the reports.
    """

    def __init__(self, engine, dump=None, interval=1, types=None):
        self.engine = engine
        self.dump = dump
        self.interval = int(interval)
        self.types = types
        self.seconds = 0.0
        self.calls = 0
        if types is None:
            self._ids = None
        else:
            self._ids = np.flatnonzero(np.asarray(types) == 1) + 1
        if dump:
            open(dump, "w").close()

    def describeNextReport(self, simulation):
        step = int(simulation.currentStep)
        remain = self.interval - (step % self.interval)
        return {"steps": remain, "periodic": None, "include": ["positions"]}

    def report(self, simulation, state):
        del state
        from openmm.unit import angstrom

        t0 = time.perf_counter()
        pulled = simulation.context.getState(getPositions=True, enforcePeriodicBox=True)
        pos = pulled.getPositions(asNumpy=True).value_in_unit(angstrom)
        vecs = pulled.getPeriodicBoxVectors(asNumpy=True).value_in_unit(angstrom)
        box = self.engine.box_from_rows(vecs)
        copy_s = time.perf_counter() - t0
        step = int(simulation.currentStep)
        labels, _, _ = self.engine.observe(step, pos, box, None, types=self.types, copy_s=copy_s)
        self.seconds += time.perf_counter() - t0
        self.calls += 1
        if self.dump:
            self._write(step, np.asarray(pos, dtype=np.float64), box, labels)

    def _write(self, step, pos, box, labels):
        ids = self._ids if self._ids is not None else np.arange(1, pos.shape[0] + 1)
        # tolist() gives Python floats; repr() of those round-trips exactly.
        sites = pos[ids - 1].tolist()
        box = [float(v) for v in box]
        lx, ly, lz, xy, xz, yz = box[3:9]
        xlo = box[0] + min(0.0, xy, xz, xy + xz)
        xhi = box[0] + lx + max(0.0, xy, xz, xy + xz)
        ylo = box[1] + min(0.0, yz)
        yhi = box[1] + ly + max(0.0, yz)
        with open(self.dump, "a") as fh:
            fh.write(f"ITEM: TIMESTEP\n{step}\nITEM: NUMBER OF ATOMS\n{len(ids)}\n")
            fh.write("ITEM: BOX BOUNDS xy xz yz pp pp pp\n")
            fh.write(f"{xlo!r} {xhi!r} {xy!r}\n{ylo!r} {yhi!r} {xz!r}\n")
            fh.write(f"{box[2]!r} {box[2] + lz!r} {yz!r}\n")
            fh.write("ITEM: ATOMS id type x y z label\n")
            for atom, (x, y, z), label in zip(ids.tolist(), sites, labels.tolist()):
                fh.write(f"{atom} 1 {x!r} {y!r} {z!r} {label}\n")
