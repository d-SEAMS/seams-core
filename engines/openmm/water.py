# Copyright (c) d-SEAMS developers
"""TIP4P/Ew boxes for the OpenMM reporter.

``box_nm=2`` is a few hundred waters. ``box_nm=4.93`` is about four
thousand. Both are built when the caller runs. The coordinates are not
stored in the repository.
"""

from openmm import LangevinMiddleIntegrator, Platform, Vec3
from openmm.app import ForceField, HBonds, Modeller, PME, Simulation, Topology
from openmm.unit import kelvin, nanometer, picosecond, picoseconds


def oxygen_indices(topology):
    from openmm.app import element

    return [
        i
        for i, atom in enumerate(topology.atoms())
        if atom.element == element.oxygen
    ]


def _platform():
    for name in ("CPU", "Reference"):
        try:
            return Platform.getPlatformByName(name), name
        except Exception:
            continue
    raise RuntimeError("OpenMM has neither the CPU nor the Reference platform")


def build_simulation(box_nm=2.0, seed=1, threads=1):
    """TIP4P/Ew waters in a cubic box of edge ``box_nm`` nanometres."""
    forcefield = ForceField("tip4pew.xml")
    modeller = Modeller(Topology(), [])
    modeller.addSolvent(
        forcefield,
        model="tip4pew",
        boxSize=Vec3(box_nm, box_nm, box_nm) * nanometer,
    )
    system = forcefield.createSystem(
        modeller.topology,
        nonbondedMethod=PME,
        nonbondedCutoff=0.9 * nanometer,
        constraints=HBonds,
        rigidWater=True,
    )
    integrator = LangevinMiddleIntegrator(
        300 * kelvin, 1 / picosecond, 0.002 * picoseconds
    )
    integrator.setRandomNumberSeed(seed)
    platform, name = _platform()
    properties = {"Threads": str(threads)} if name == "CPU" else {}
    simulation = Simulation(
        modeller.topology,
        system,
        integrator,
        platform,
        properties,
    )
    simulation.context.setPositions(modeller.positions)
    simulation.minimizeEnergy()
    simulation.context.setVelocitiesToTemperature(300 * kelvin, seed)
    return simulation, oxygen_indices(modeller.topology)


def save_state(simulation):
    return simulation.context.getState(getPositions=True, getVelocities=True)


def restore(simulation, state, seed=1):
    simulation.context.setPositions(state.getPositions())
    simulation.context.setVelocities(state.getVelocities())
    simulation.context.setTime(0.0)
    simulation.currentStep = 0
    simulation.integrator.setRandomNumberSeed(seed)
