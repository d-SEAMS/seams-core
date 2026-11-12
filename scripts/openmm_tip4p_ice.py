# Copyright (c) d-SEAMS developers
"""TIP4P/Ice in OpenMM, with the in-situ reporter, checked against the file reader.

Usage:
    python3 scripts/openmm_tip4p_ice.py LIBSEAMS_LANG SEAMS LOG

4096 molecules on the CPU, then 32768. The OpenMM thread count and the
classifier's OpenMP team are the same number. One match trajectory is
written and read back. The 4096-molecule timing is 5000 steps at report
intervals 100, 500 and 1000. The 32768-molecule timing is shorter and
reports every 100 steps.
"""

from __future__ import annotations

import ctypes
import os
import re
import subprocess
import sys
import time
from pathlib import Path

import numpy as np

HERE = Path(__file__).resolve().parent
ROOT = HERE.parent
sys.path.insert(0, str(ROOT / "adapters" / "python"))

from dseams_insitu.engine import InSitu, _NAMES
from dseams_insitu.openmm import DseamsReporter

MOLECULES = int(os.environ.get("OPENMM_TIP4P_MOLECULES", "4096"))
STEPS = int(os.environ.get("OPENMM_TIP4P_STEPS", "5000"))
REPEATS = int(os.environ.get("OPENMM_TIP4P_REPEATS", "6"))
INTERVALS = tuple(
    int(item)
    for item in os.environ.get("OPENMM_TIP4P_INTERVALS", "100,500,1000").split(",")
    if item
)
THREADS = tuple(
    int(item)
    for item in os.environ.get("OPENMM_TIP4P_THREADS", "1,8").split(",")
    if item
)
MATCH_STEPS = int(os.environ.get("OPENMM_TIP4P_MATCH_STEPS", "400"))
MATCH_INTERVAL = int(os.environ.get("OPENMM_TIP4P_MATCH_INTERVAL", "100"))
LARGE_MOLECULES = int(os.environ.get("OPENMM_TIP4P_LARGE", "32768"))
LARGE_STEPS = int(os.environ.get("OPENMM_TIP4P_LARGE_STEPS", "1000"))
LARGE_REPEATS = int(os.environ.get("OPENMM_TIP4P_LARGE_REPEATS", "6"))
LARGE_INTERVAL = int(os.environ.get("OPENMM_TIP4P_LARGE_INTERVAL", "100"))
DT_NS = 2.0e-6
CUTOFF = 3.5
OUT = Path(os.environ.get("OPENMM_TIP4P_DIR", "openmm-tip4p"))


def water_box(n_mol, threads):
    from openmm.app import ForceField, Modeller, PME, Topology
    from openmm import LangevinMiddleIntegrator, Platform
    from openmm.unit import kelvin, nanometer, picosecond

    ff = ForceField(str(ROOT / "adapters" / "python" / "dseams_insitu" / "tip4pice.xml"))
    modeller = Modeller(Topology(), [])
    modeller.addSolvent(ff, model="tip4pew", numAdded=n_mol)
    span = modeller.topology.getUnitCellDimensions().value_in_unit(nanometer)
    # PME needs the cutoff inside half the shortest edge.
    cutoff = min(0.9, 0.49 * min(span))
    system = ff.createSystem(
        modeller.topology,
        nonbondedMethod=PME,
        nonbondedCutoff=cutoff * nanometer,
        constraints=None,
        rigidWater=True,
    )
    integrator = LangevinMiddleIntegrator(260 * kelvin, 1.0 / picosecond, 0.002 * picosecond)
    platform = Platform.getPlatformByName("CPU")
    from openmm.app import Simulation

    simulation = Simulation(
        modeller.topology, system, integrator, platform, {"Threads": str(int(threads))}
    )
    simulation.context.setPositions(modeller.positions)
    simulation.integrator.setRandomNumberSeed(1)
    return simulation


def oxygen_index(topology):
    index = []
    for atom in topology.atoms():
        index.append(1 if atom.element is not None and atom.element.symbol == "O" else 0)
    return np.asarray(index, dtype=np.int32)


def num(value):
    return format(float(value), ".17g")


def write_dump(path, frames):
    with path.open("w") as out:
        for number, frame in enumerate(frames, start=1):
            pos = frame["positions"]
            box = frame["box"]
            n = pos.shape[0]
            out.write(f"ITEM: TIMESTEP\n{number}\n")
            out.write(f"ITEM: NUMBER OF ATOMS\n{n}\n")
            lx, ly, lz = box[3], box[4], box[5]
            xy, xz, yz = box[6], box[7], box[8]
            tilted = abs(xy) > 0.0 or abs(xz) > 0.0 or abs(yz) > 0.0
            if tilted:
                out.write("ITEM: BOX BOUNDS xy xz yz pp pp pp\n")
                out.write(f"0 {num(lx)} {num(xy)}\n")
                out.write(f"0 {num(ly)} {num(xz)}\n")
                out.write(f"0 {num(lz)} {num(yz)}\n")
            else:
                out.write("ITEM: BOX BOUNDS pp pp pp\n")
                out.write(f"0 {num(lx)}\n")
                out.write(f"0 {num(ly)}\n")
                out.write(f"0 {num(lz)}\n")
            out.write("ITEM: ATOMS id type x y z\n")
            for i in range(n):
                x, y, z = pos[i]
                out.write(f"{i + 1} 1 {num(x)} {num(y)} {num(z)}\n")


def snapshot(simulation):
    state = simulation.context.getState(getPositions=True, getVelocities=True)
    return state


def set_classifier_threads(n):
    """The OpenMP team inside CHILL+ and the ring search."""
    os.environ["OMP_NUM_THREADS"] = str(int(n))
    gomp = ctypes.CDLL("libgomp.so.1")
    gomp.omp_set_num_threads(int(n))


def ns_per_day(steps, seconds):
    if seconds <= 0.0:
        return 0.0
    return steps * DT_NS * 86400.0 / seconds


def restore(simulation, state):
    simulation.context.setPositions(state.getPositions())
    simulation.context.setVelocities(state.getVelocities())
    simulation.currentStep = 0
    simulation.integrator.setRandomNumberSeed(1)


def run_steps(simulation, steps, reporter):
    simulation.reporters.clear()
    if reporter is not None:
        simulation.reporters.append(reporter)
    t0 = time.perf_counter()
    simulation.step(steps)
    return time.perf_counter() - t0


def histogram(labels):
    counts = {name: 0 for name in _NAMES}
    for value in labels:
        counts[_NAMES[int(value)]] += 1
    return counts


def cli_counts(seams, dump, frame):
    env = os.environ.copy()
    env["OMP_NUM_THREADS"] = "1"
    env["NO_COLOR"] = "1"
    proc = subprocess.run(
        [
            seams,
            "chill-plus",
            str(dump),
            "--type",
            "1",
            "--cutoff",
            "3.5",
            "--graph",
            "cutoff",
            "--frame",
            str(frame),
            "--last",
            str(frame),
        ],
        check=True,
        capture_output=True,
        text=True,
        env=env,
    )
    text = re.sub(r"\x1b\[[0-9;]*m", "", proc.stdout)
    parts = text.split()
    found = {}
    i = 0
    while i < len(parts):
        if parts[i] in _NAMES and i + 1 < len(parts) and parts[i + 1].isdigit():
            found[parts[i]] = int(parts[i + 1])
            i += 2
            continue
        i += 1
    return found, text


def check_skin(library, positions, box):
    """A small drift must keep the cutoff pairs a fresh list would build."""
    fresh = InSitu(library, skin=0.0)
    fresh.keep_frames = False
    base_labels, _base_counts, base_rings = fresh.analyse(positions, box)
    held = InSitu(library, skin=0.5)
    held.keep_frames = False
    held.analyse(positions, box)
    labels, _counts, rings = held.analyse(positions, box)
    if held._phase.reused != 1:
        raise SystemExit("identical coordinates rebuilt the neighbour list")
    if not np.array_equal(labels, base_labels) or not np.array_equal(rings, base_rings):
        raise SystemExit("reused list changed an unmoved frame")
    rng = np.random.default_rng(215)
    shift = rng.normal(size=positions.shape)
    peak = np.max(np.linalg.norm(shift, axis=1))
    moved = np.array(positions, copy=True)
    moved += shift * (0.02 / peak)
    expect_labels, _expect_counts, expect_rings = fresh.analyse(moved, box)
    labels, _counts, rings = held.analyse(moved, box)
    if held._phase.reused != 1:
        raise SystemExit("a 0.02 Å drift rebuilt the neighbour list")
    if not np.array_equal(labels, expect_labels) or not np.array_equal(rings, expect_rings):
        raise SystemExit("skin reuse does not match a fresh cutoff list")
    far = np.array(positions, copy=True)
    far[0, 0] += 2.0
    expect_labels, _expect_counts, expect_rings = fresh.analyse(far, box)
    labels, _counts, rings = held.analyse(far, box)
    if held._phase.reused != 0:
        raise SystemExit("a 2 Å drift reused the neighbour list")
    if not np.array_equal(labels, expect_labels) or not np.array_equal(rings, expect_rings):
        raise SystemExit("rebuilt list does not match a fresh cutoff list")


def mean(values):
    values = list(values)
    return sum(values) / len(values)


def main(argv):
    if len(argv) < 4:
        print(__doc__, file=sys.stderr)
        return 2
    library, seams, log_path = argv[1], argv[2], Path(argv[3])
    os.environ["SEAMS_LANG_LIB"] = library
    os.environ["OMP_NUM_THREADS"] = "1"
    OUT.mkdir(parents=True, exist_ok=True)
    log_path.parent.mkdir(parents=True, exist_ok=True)
    lines = []

    def say(text):
        print(text, flush=True)
        lines.append(text)
        log_path.write_text("\n".join(lines) + "\n")

    say(f"molecules {MOLECULES}")
    say(f"steps {STEPS} intervals {','.join(str(v) for v in INTERVALS)} repeats {REPEATS}")
    say(f"threads {','.join(str(v) for v in THREADS)}")
    say("platform CPU")
    say(f"library {library}")
    say(f"seams {seams}")
    say("classifier-threads follow the OpenMM thread count")
    set_classifier_threads(1)

    from openmm import NonbondedForce
    from openmm.unit import nanometer

    prepared = {}
    for threads in THREADS:
        simulation = water_box(MOLECULES, threads)
        if threads == THREADS[0]:
            say(f"atoms {simulation.system.getNumParticles()}")
            edges = simulation.topology.getUnitCellDimensions().value_in_unit(nanometer)
            say("box-nm " + " ".join(f"{float(v):.6f}" for v in edges))
            for force in simulation.system.getForces():
                if isinstance(force, NonbondedForce):
                    say(
                        "pme-cutoff-nm "
                        f"{force.getCutoffDistance().value_in_unit(nanometer):.6f}"
                    )
            simulation.minimizeEnergy(maxIterations=40)
            kinds = oxygen_index(simulation.topology)
            say(f"oxygens {int(kinds.sum())}")
            state = snapshot(simulation)
            prepared["kinds"] = kinds
            prepared["positions"] = None
        else:
            simulation.context.setPositions(prepared["positions"])
            simulation.context.setVelocities(prepared["velocities"])
        if prepared["positions"] is None:
            held = snapshot(simulation)
            prepared["positions"] = held.getPositions()
            prepared["velocities"] = held.getVelocities()
        prepared[threads] = simulation

    kinds = prepared["kinds"]
    match_sim = prepared[THREADS[0]]

    engine = InSitu(library, skin=0.0)
    pulled = match_sim.context.getState(getPositions=True, enforcePeriodicBox=True)
    from openmm.unit import angstrom

    pos = pulled.getPositions(asNumpy=True).value_in_unit(angstrom)
    oxy = pos[kinds == 1]
    box = engine.box_from_rows(
        pulled.getPeriodicBoxVectors(asNumpy=True).value_in_unit(angstrom)
    )
    say("skin-check")
    check_skin(library, np.array(oxy, copy=True), np.array(box, copy=True))
    say("skin-check exact")

    reporter = DseamsReporter(
        engine, str(OUT / "reporter.dump"), interval=MATCH_INTERVAL, types=kinds
    )
    say("match-run")
    engine.frames.clear()
    engine.reset_profile()
    elapsed = run_steps(match_sim, MATCH_STEPS, reporter)
    say(f"match-run seconds {elapsed:.6f} reports {len(engine.frames)}")
    dump = OUT / "oxygens.lammpstrj"
    write_dump(dump, engine.frames)
    for index, frame in enumerate(engine.frames, start=1):
        labels, counts, rings = engine.analyse_file(dump, index)
        if not np.array_equal(labels, frame["labels"]):
            raise SystemExit(f"frame {index} CHILL+ labels differ")
        if not np.array_equal(counts, frame["counts"]):
            raise SystemExit(f"frame {index} CHILL+ counts differ")
        if not np.array_equal(rings, frame["rings"]):
            raise SystemExit(f"frame {index} ring counts differ")
        cli, raw = cli_counts(seams, dump, index)
        live = histogram(frame["labels"])
        for name, number in live.items():
            if number == 0:
                continue
            if cli.get(name) != number:
                raise SystemExit(f"frame {index} CLI {name} {cli.get(name)} != {number}\n{raw}")
        say(
            f"frame {index} step {frame['step']} labels-match rings "
            + " ".join(f"{size}:{int(frame['rings'][size])}" for size in range(3, 7))
        )
    say("match exact")

    def rewind(simulation):
        simulation.context.setPositions(prepared["positions"])
        simulation.context.setVelocities(prepared["velocities"])
        simulation.currentStep = 0
        simulation.integrator.setRandomNumberSeed(1)

    for threads in THREADS:
        set_classifier_threads(threads)
        simulation = prepared[threads]
        say(f"threads {threads}")
        bare = []
        for trial in range(1, REPEATS + 1):
            rewind(simulation)
            seconds = run_steps(simulation, STEPS, None)
            bare.append(seconds)
            say(f"bare {threads} {trial} {seconds:.6f} ns/day {ns_per_day(STEPS, seconds):.6f}")
        bare_mean = mean(bare)
        say(f"bare-mean {threads} {bare_mean:.6f} ns/day {ns_per_day(STEPS, bare_mean):.6f}")
        for interval in INTERVALS:
            engine.keep_frames = False
            engine.frames.clear()
            reporter.interval = interval
            reported = []
            profiles = []
            for trial in range(1, REPEATS + 1):
                rewind(simulation)
                engine.reset_profile()
                seconds = run_steps(simulation, STEPS, reporter)
                reported.append(seconds)
                profiles.append(dict(engine.profile))
                say(
                    f"reporter {threads} {interval} {trial} {seconds:.6f} "
                    f"ns/day {ns_per_day(STEPS, seconds):.6f} "
                    f"reports {engine.profile['reports']} reused {engine.profile['reused']}"
                )
            rep_mean = mean(reported)
            reports = profiles[-1]["reports"]
            if reports <= 0:
                raise SystemExit(f"interval {interval} produced no reports in {STEPS} steps")
            per_report = (rep_mean - bare_mean) / reports * 1.0e3
            bare_rate = ns_per_day(STEPS, bare_mean)
            rep_rate = ns_per_day(STEPS, rep_mean)
            slowdown = (bare_rate - rep_rate) / bare_rate if bare_rate else 0.0

            def phase(name):
                return mean(item[name] / item["reports"] for item in profiles) * 1.0e3

            say(f"reporter-mean {threads} {interval} {rep_mean:.6f} ns/day {rep_rate:.6f}")
            say(f"per-report-ms {threads} {interval} {per_report:.6f}")
            say(f"slowdown {threads} {interval} {slowdown:.6f}")
            say(
                f"split-ms {threads} {interval} "
                f"copy {phase('copy'):.6f} neighbours {phase('neighbours'):.6f} "
                f"chill {phase('chill'):.6f} rings {phase('rings'):.6f} "
                f"reused {mean(item['reused'] / item['reports'] for item in profiles):.6f}"
            )

    say(f"large-molecules {LARGE_MOLECULES}")
    say(f"large-steps {LARGE_STEPS} interval {LARGE_INTERVAL} repeats {LARGE_REPEATS}")
    large = {}
    for threads in THREADS:
        set_classifier_threads(1)
        simulation = water_box(LARGE_MOLECULES, threads)
        if "positions" not in large:
            simulation.minimizeEnergy(maxIterations=25)
            held = snapshot(simulation)
            large["positions"] = held.getPositions()
            large["velocities"] = held.getVelocities()
            large["kinds"] = oxygen_index(simulation.topology)
            say(f"large-atoms {simulation.system.getNumParticles()}")
            say(f"large-oxygens {int(large['kinds'].sum())}")
        else:
            simulation.context.setPositions(large["positions"])
            simulation.context.setVelocities(large["velocities"])
        large[threads] = simulation
    large_engine = InSitu(library, skin=0.0)
    large_reporter = DseamsReporter(
        large_engine, str(OUT / "large.dump"), interval=LARGE_INTERVAL, types=large["kinds"]
    )

    def rewind_large(simulation):
        simulation.context.setPositions(large["positions"])
        simulation.context.setVelocities(large["velocities"])
        simulation.currentStep = 0
        simulation.integrator.setRandomNumberSeed(1)

    set_classifier_threads(1)
    large_engine.keep_frames = True
    large_engine.frames.clear()
    say("large-match")
    run_steps(large[THREADS[0]], LARGE_INTERVAL, large_reporter)
    write_dump(OUT / "large-oxygens.lammpstrj", large_engine.frames)
    if large_engine.frames:
        labels, counts, rings = large_engine.analyse_file(OUT / "large-oxygens.lammpstrj", 1)
        frame = large_engine.frames[0]
        if not np.array_equal(labels, frame["labels"]) or not np.array_equal(rings, frame["rings"]):
            raise SystemExit("32768-molecule frame does not match the file reader")
        say(
            f"large-frame step {frame['step']} labels-match rings "
            + " ".join(f"{size}:{int(frame['rings'][size])}" for size in range(3, 7))
        )
    say("large-match exact")

    for threads in THREADS:
        set_classifier_threads(threads)
        simulation = large[threads]
        say(f"large-threads {threads}")
        bare = []
        for trial in range(1, LARGE_REPEATS + 1):
            rewind_large(simulation)
            seconds = run_steps(simulation, LARGE_STEPS, None)
            bare.append(seconds)
            say(
                f"large-bare {threads} {trial} {seconds:.6f} "
                f"ns/day {ns_per_day(LARGE_STEPS, seconds):.6f}"
            )
        bare_mean = mean(bare)
        say(
            f"large-bare-mean {threads} {bare_mean:.6f} "
            f"ns/day {ns_per_day(LARGE_STEPS, bare_mean):.6f}"
        )
        large_engine.keep_frames = False
        large_engine.frames.clear()
        large_reporter.interval = LARGE_INTERVAL
        reported = []
        profiles = []
        for trial in range(1, LARGE_REPEATS + 1):
            rewind_large(simulation)
            large_engine.reset_profile()
            seconds = run_steps(simulation, LARGE_STEPS, large_reporter)
            reported.append(seconds)
            profiles.append(dict(large_engine.profile))
            say(
                f"large-reporter {threads} {trial} {seconds:.6f} "
                f"ns/day {ns_per_day(LARGE_STEPS, seconds):.6f} "
                f"reports {large_engine.profile['reports']}"
            )
        rep_mean = mean(reported)
        reports = profiles[-1]["reports"]
        if reports <= 0:
            raise SystemExit("32768-molecule run produced no reports")
        per_report = (rep_mean - bare_mean) / reports * 1.0e3
        bare_rate = ns_per_day(LARGE_STEPS, bare_mean)
        rep_rate = ns_per_day(LARGE_STEPS, rep_mean)
        slowdown = (bare_rate - rep_rate) / bare_rate if bare_rate else 0.0

        def large_phase(name):
            return mean(item[name] / item["reports"] for item in profiles) * 1.0e3

        say(f"large-reporter-mean {threads} {rep_mean:.6f} ns/day {rep_rate:.6f}")
        say(f"large-per-report-ms {threads} {per_report:.6f}")
        say(f"large-slowdown {threads} {slowdown:.6f}")
        say(
            f"large-split-ms {threads} "
            f"copy {large_phase('copy'):.6f} neighbours {large_phase('neighbours'):.6f} "
            f"chill {large_phase('chill'):.6f} rings {large_phase('rings'):.6f}"
        )
    log_path.write_text("\n".join(lines) + "\n")
    return 0


if __name__ == "__main__":
    raise SystemExit(main(sys.argv))
