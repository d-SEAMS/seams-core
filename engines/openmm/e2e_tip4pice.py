# Copyright (c) d-SEAMS developers
"""TIP4P/Ice ice Ic slab and liquid in OpenMM with the d-SEAMS reporter.

usage: e2e_tip4pice.py BUILD [--steps S] [--record FILE]

Needs the openmm Python package. DseamsReporter pulls positions from the
context every --every steps, classifies the oxygens through
libseams_lang, and dumps the coordinates it classified with their labels.
The CLI then classifies that dump and every label must agree.
"""

import argparse
import importlib.util
import sys
import tempfile
from pathlib import Path

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[1]
# Appended: engines/openmm would otherwise shadow the openmm package.
sys.path.append(str(HERE.parent))
import e2e_common as e2e  # noqa: E402
import tip4pice_slab  # noqa: E402


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("build")
    ap.add_argument("--cells", type=int, nargs=3, default=(6, 6, 3))
    ap.add_argument("--steps", type=int, default=200)
    ap.add_argument("--every", type=int, default=10)
    ap.add_argument("--threads", type=int, default=2)
    ap.add_argument("--record", default="")
    args = ap.parse_args()

    spec = importlib.util.find_spec("openmm")
    if spec is None or spec.origin is None:
        e2e.skip("the openmm Python package was not found")
    lib = e2e.find_build_file(args.build, "libseams_lang.so")
    if lib is None:
        sys.exit(f"libseams_lang.so was not built under {args.build}")
    cli = e2e.seams_cli(args.build)

    import numpy as np
    import openmm
    from openmm import app, unit

    sys.path.insert(0, str(ROOT / "adapters" / "python"))
    from dseams_insitu import DseamsReporter, InSitu

    with tempfile.TemporaryDirectory(prefix="dseams-openmm-") as tmp:
        work = Path(tmp)
        mols, box = tip4pice_slab.slab(*args.cells)
        tip4pice_slab.write_pdb(work / "ice.pdb", mols, box)
        pdb = app.PDBFile(str(work / "ice.pdb"))
        ff = app.ForceField(str(ROOT / "adapters" / "python" / "dseams_insitu" / "tip4pice.xml"))
        lx, ly, lz = (b / 10 for b in box)
        half = lz / 2
        # The liquid is a separate box stacked on the slab, so no water
        # lands in the lattice interstitials.
        height = half - 0.2
        liquid = app.Modeller(app.Topology(), [])
        liquid.addSolvent(ff, model="tip4pew", boxSize=openmm.Vec3(lx, ly, height) * unit.nanometer)
        # addSolvent centres its box on the origin.
        lifted = [openmm.Vec3(p.x + lx / 2, p.y + ly / 2, p.z + height / 2 + half + 0.1)
                  for p in liquid.positions.value_in_unit(unit.nanometer)]
        model = app.Modeller(pdb.topology, pdb.positions)
        model.add(liquid.topology, lifted * unit.nanometer)
        model.topology.setPeriodicBoxVectors(
            [openmm.Vec3(lx, 0, 0), openmm.Vec3(0, ly, 0), openmm.Vec3(0, 0, lz)] * unit.nanometer)
        system = ff.createSystem(model.topology, nonbondedMethod=app.PME,
                                 nonbondedCutoff=0.9 * unit.nanometer, rigidWater=True)
        integrator = openmm.LangevinMiddleIntegrator(250 * unit.kelvin, 1 / unit.picosecond,
                                                     0.002 * unit.picoseconds)
        integrator.setRandomNumberSeed(2014)
        platform = openmm.Platform.getPlatformByName("CPU")
        sim = app.Simulation(model.topology, system, integrator, platform,
                             {"Threads": str(args.threads)})
        sim.context.setPositions(model.positions)
        sim.minimizeEnergy(maxIterations=500)
        sim.context.setVelocitiesToTemperature(250 * unit.kelvin, 2014)

        kinds = np.array([1 if atom.name == "O" else 0 for atom in model.topology.atoms()],
                         dtype=np.int32)
        engine = InSitu(str(lib), cutoff=3.5)
        engine.keep_frames = False
        reporter = DseamsReporter(engine, dump=str(work / "insitu.lammpstrj"),
                                  interval=args.every, types=kinds)
        sim.reporters.append(reporter)
        sim.step(args.steps)

        insitu = e2e.read_dump(work / "insitu.lammpstrj", "label")
        offline = e2e.offline_labels(cli, work / "insitu.lammpstrj", len(insitu), 3.5, work,
                                     extra=("--type", "1"))
        frames, atom_frames, bad = e2e.compare(insitu, offline)
        if [s for s, _ in insitu] != list(range(args.every, args.steps + 1, args.every)):
            sys.exit("the reporter did not fire on every interval")
        record = {
            "engine": "openmm",
            "version": openmm.Platform.getOpenMMVersion(),
            "system": f"TIP4P/Ice ice Ic slab and liquid, {int(kinds.sum())} molecules",
            "steps": args.steps,
            "frames": frames,
            "atom_frames": atom_frames,
            "labels": e2e.census(insitu),
            "insitu_s_per_frame": reporter.seconds / max(reporter.calls, 1),
        }
        e2e.finish(args.record, record, bad)


if __name__ == "__main__":
    main()
