# Copyright (c) d-SEAMS developers
"""TIP4P/Ice seeded growth in OpenMM, classified in place.

A cubic-ice slab sits in a taller orthorhombic cell. OpenMM fills the
empty half with TIP4P/Ice and integrates at 255 K, below the model's
melting point. The in-situ reporter logs the CHILL+ ice count and the
primitive-ring census on the oxygen sites. Saved frames stay outside
the tree. An offline read of that dump must reproduce every label and
every ring count.

Usage:
    python3 scripts/openmm_tip4p_seed.py LIBSEAMS_LANG SEAMS LOG

The science trajectory is one nanosecond (500000 steps of 2 fs) on the
requested thread count. A separate pair of bare and reporter runs, on
1 and 8 threads, records ns/day and the reporter overhead. Set
OPENMM_SEED_DIR to choose the frame directory. The default is
`openmm-tip4p-seed`.
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

D_OO = 2.75
R_OH = 0.9572
ANGLE = 104.52
W_O = 0.7308332982689929
W_H = 0.13458335086550355
CELLS = int(os.environ.get("OPENMM_SEED_CELLS", "3"))
TEMP_K = float(os.environ.get("OPENMM_SEED_TEMP", "255"))
STEPS = int(os.environ.get("OPENMM_SEED_STEPS", "500000"))
INTERVAL = int(os.environ.get("OPENMM_SEED_INTERVAL", "5000"))
SCIENCE_THREADS = int(os.environ.get("OPENMM_SEED_THREADS", "8"))
TIME_STEPS = int(os.environ.get("OPENMM_SEED_TIME_STEPS", "50000"))
TIME_INTERVAL = int(os.environ.get("OPENMM_SEED_TIME_INTERVAL", "100"))
TIME_THREADS = tuple(
    int(item)
    for item in os.environ.get("OPENMM_SEED_TIME_THREADS", "1,8").split(",")
    if item
)
DT_NS = 2.0e-6
OUT = Path(os.environ.get("OPENMM_SEED_DIR", "openmm-tip4p-seed"))


def set_classifier_threads(n):
    os.environ["OMP_NUM_THREADS"] = str(int(n))
    gomp = ctypes.CDLL("libgomp.so.1")
    gomp.omp_set_num_threads(int(n))


def ns_per_day(steps, seconds):
    if seconds <= 0.0:
        return 0.0
    return steps * DT_NS * 86400.0 / seconds


def diamond_oxygens(n_cell, d_oo=D_OO):
    """Oxygen sites of cubic ice, one conventional cell repeated n_cell times."""
    a = 4.0 * d_oo / np.sqrt(3.0)
    basis = np.array(
        [
            [0.0, 0.0, 0.0],
            [0.5, 0.5, 0.0],
            [0.5, 0.0, 0.5],
            [0.0, 0.5, 0.5],
            [0.25, 0.25, 0.25],
            [0.75, 0.75, 0.25],
            [0.75, 0.25, 0.75],
            [0.25, 0.75, 0.75],
        ],
        dtype=np.float64,
    )
    pts = []
    for ix in range(n_cell):
        for iy in range(n_cell):
            for iz in range(n_cell):
                origin = np.array([ix, iy, iz], dtype=np.float64)
                pts.append((basis + origin) * a)
    pos = np.concatenate(pts, axis=0)
    return pos, np.array([a * n_cell, a * n_cell, a * n_cell], dtype=np.float64)


def pair_list(pos, box, wrap_z):
    n = pos.shape[0]
    inv = 1.0 / box
    cutoff2 = (D_OO * 1.15) ** 2
    nbrs = [[] for _ in range(n)]
    for i in range(n):
        delta = pos[i + 1 :] - pos[i]
        delta[:, 0] -= box[0] * np.round(delta[:, 0] * inv[0])
        delta[:, 1] -= box[1] * np.round(delta[:, 1] * inv[1])
        if wrap_z:
            delta[:, 2] -= box[2] * np.round(delta[:, 2] * inv[2])
        dist2 = np.einsum("ij,ij->i", delta, delta)
        hits = np.flatnonzero((dist2 > 0.04) & (dist2 < cutoff2))
        for h in hits:
            j = i + 1 + int(h)
            nbrs[i].append(j)
            nbrs[j].append(i)
    return nbrs


def orient_donors(nbrs, seed=1):
    """Point two hydrogens out of each bulk oxygen, without two donors on one bond."""
    n = len(nbrs)
    edges = []
    for i, row in enumerate(nbrs):
        for j in row:
            if i < j:
                edges.append((i, j))
    rng = np.random.default_rng(seed)
    direction = {edge: (1 if rng.random() < 0.5 else -1) for edge in edges}

    def out_degree(vertex):
        count = 0
        for other in nbrs[vertex]:
            a, b = (vertex, other) if vertex < other else (other, vertex)
            arrow = direction[(a, b)]
            if (vertex == a and arrow == 1) or (vertex == b and arrow == -1):
                count += 1
        return count

    def outgoing(vertex):
        found = []
        for other in nbrs[vertex]:
            a, b = (vertex, other) if vertex < other else (other, vertex)
            arrow = direction[(a, b)]
            if (vertex == a and arrow == 1) or (vertex == b and arrow == -1):
                found.append(other)
        return found

    def reverse_edge(a, b):
        key = (a, b) if a < b else (b, a)
        direction[key] = -direction[key]

    guard = 0
    while True:
        high = [v for v in range(n) if out_degree(v) > 2]
        if not high:
            break
        guard += 1
        if guard > n * 8:
            raise RuntimeError("ice-rule orientation did not converge")
        parent = {high[0]: None}
        queue = [high[0]]
        sink = None
        while queue and sink is None:
            cur = queue.pop(0)
            for nxt in outgoing(cur):
                if nxt in parent:
                    continue
                parent[nxt] = cur
                if out_degree(nxt) < 2:
                    sink = nxt
                    break
                queue.append(nxt)
        if sink is None:
            raise RuntimeError("an over-coordinated oxygen has no path to a deficit")
        while parent[sink] is not None:
            prev = parent[sink]
            reverse_edge(prev, sink)
            sink = prev
    donors = []
    for v in range(n):
        donors.append(outgoing(v)[:2])
    return donors


def mic(delta, box, wrap_z):
    out = np.array(delta, dtype=np.float64, copy=True)
    out[0] -= box[0] * np.round(out[0] / box[0])
    out[1] -= box[1] * np.round(out[1] / box[1])
    if wrap_z:
        out[2] -= box[2] * np.round(out[2] / box[2])
    return out


def tip4p_sites(oxygen, toward):
    """Rigid TIP4P geometry whose bisector follows the two donor directions."""
    u = toward[0] / np.linalg.norm(toward[0])
    v = toward[1] / np.linalg.norm(toward[1])
    bis = u + v
    bis /= np.linalg.norm(bis)
    side = u - np.dot(u, bis) * bis
    side /= np.linalg.norm(side)
    half = np.deg2rad(ANGLE / 2.0)
    d1 = np.cos(half) * bis + np.sin(half) * side
    d2 = np.cos(half) * bis - np.sin(half) * side
    h1 = oxygen + R_OH * d1
    h2 = oxygen + R_OH * d2
    msite = W_O * oxygen + W_H * h1 + W_H * h2
    return h1, h2, msite


def seed_molecules(n_cell):
    pos, box = diamond_oxygens(n_cell)
    # Fully periodic crystal: every oxygen has four neighbours.
    crystal_nbrs = pair_list(pos, box, wrap_z=True)
    degrees = sorted(len(row) for row in crystal_nbrs)
    if degrees[0] != 4 or degrees[-1] != 4:
        raise RuntimeError(f"cubic ice coordination is {degrees[0]}..{degrees[-1]}, expected 4")
    # Open the cell along z. The two faces become the crystal/liquid interfaces.
    donors = orient_donors(pair_list(pos, box, wrap_z=False))
    tall = np.array([box[0], box[1], box[2] * 2.0], dtype=np.float64)
    molecules = []
    for i, oxygen in enumerate(pos):
        vectors = [mic(pos[j] - oxygen, box, wrap_z=False) for j in donors[i][:2]]
        if len(vectors) < 2:
            # A surface oxygen can be left with one donor. The second
            # hydrogen stays in the molecular plane, off the first bond,
            # so the rigid TIP4P angle is still defined.
            outward = np.array(
                [0.0, 0.0, 1.0 if oxygen[2] > box[2] / 2.0 else -1.0],
                dtype=np.float64,
            )
            if not vectors:
                vectors = [outward, np.array([1.0, 0.0, 0.0])]
            else:
                first = vectors[0] / np.linalg.norm(vectors[0])
                aux = (
                    np.array([1.0, 0.0, 0.0])
                    if abs(first[0]) < 0.9
                    else np.array([0.0, 1.0, 0.0])
                )
                vectors.append(np.cross(first, aux))
        h1, h2, msite = tip4p_sites(oxygen, vectors)
        molecules.append((oxygen, h1, h2, msite))
    return molecules, tall, pos, box


def prepare_box(n_cell):
    """One topology. The crystal is the lower half; solvent fills the rest."""
    from openmm.app import Element, ForceField, Modeller, Topology
    from openmm import Vec3
    from openmm.unit import nanometer

    molecules, tall, crystal, crystal_box = seed_molecules(n_cell)
    ff = ForceField(str(ROOT / "adapters" / "python" / "dseams_insitu" / "tip4pice.xml"))
    topology = Topology()
    chain = topology.addChain()
    oxygen = Element.getBySymbol("O")
    hydrogen = Element.getBySymbol("H")
    coords = np.empty((len(molecules) * 4, 3), dtype=np.float64)
    at = 0
    for o, h1, h2, msite in molecules:
        residue = topology.addResidue("HOH", chain)
        ao = topology.addAtom("O", oxygen, residue)
        ah1 = topology.addAtom("H1", hydrogen, residue)
        ah2 = topology.addAtom("H2", hydrogen, residue)
        topology.addAtom("M", None, residue)
        topology.addBond(ao, ah1)
        topology.addBond(ao, ah2)
        coords[at : at + 4] = (o, h1, h2, msite)
        at += 4
    modeller = Modeller(topology, (coords / 10.0) * nanometer)
    # boxSize fills the existing orthorhombic cell. The crystal already
    # occupies the lower half, so the added water is the liquid reservoir.
    box_nm = Vec3(float(tall[0]) / 10.0, float(tall[1]) / 10.0, float(tall[2]) / 10.0)
    modeller.addSolvent(ff, model="tip4pew", boxSize=box_nm * nanometer)
    span = modeller.topology.getUnitCellDimensions().value_in_unit(nanometer)
    cutoff = min(0.9, 0.49 * float(min(span)))
    return ff, modeller, cutoff, crystal, crystal_box


def make_simulation(ff, modeller, cutoff, threads):
    from openmm.app import PME, Simulation
    from openmm import LangevinMiddleIntegrator, Platform
    from openmm.unit import kelvin, nanometer, picosecond

    system = ff.createSystem(
        modeller.topology,
        nonbondedMethod=PME,
        nonbondedCutoff=cutoff * nanometer,
        constraints=None,
        rigidWater=True,
    )
    integrator = LangevinMiddleIntegrator(TEMP_K * kelvin, 1.0 / picosecond, 0.002 * picosecond)
    platform = Platform.getPlatformByName("CPU")
    simulation = Simulation(
        modeller.topology,
        system,
        integrator,
        platform,
        {"Threads": str(int(threads))},
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


def append_frame(handle, frame, number):
    pos = frame["positions"]
    box = frame["box"]
    n = pos.shape[0]
    handle.write(f"ITEM: TIMESTEP\n{number}\n")
    handle.write(f"ITEM: NUMBER OF ATOMS\n{n}\n")
    lx, ly, lz = box[3], box[4], box[5]
    handle.write("ITEM: BOX BOUNDS pp pp pp\n")
    handle.write(f"0 {num(lx)}\n")
    handle.write(f"0 {num(ly)}\n")
    handle.write(f"0 {num(lz)}\n")
    handle.write("ITEM: ATOMS id type x y z\n")
    for i in range(n):
        handle.write(
            f"{i + 1} 1 {num(pos[i, 0])} {num(pos[i, 1])} {num(pos[i, 2])}\n"
        )


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


def snapshot(simulation):
    return simulation.context.getState(getPositions=True, getVelocities=True)


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
    simulation.step(int(steps))
    return time.perf_counter() - t0


class CensusReporter(DseamsReporter):
    """Logs the ice count and the ring census, and appends the oxygen frame."""

    def __init__(self, engine, dump, interval, types, log, handle):
        super().__init__(engine, dump, interval=interval, types=types)
        self.log = log
        self.handle = handle
        self.saved = 0
        self.live = []

    def report(self, simulation, state):
        super().report(simulation, state)
        frame = self.engine.frames[-1]
        self.saved += 1
        append_frame(self.handle, frame, self.saved)
        self.handle.flush()
        counts = frame["counts"]
        rings = frame["rings"]
        ice = int(counts[0] + counts[1])
        line = (
            f"census step {int(frame['step'])} ice {ice} "
            + " ".join(f"{name} {int(counts[i])}" for i, name in enumerate(_NAMES))
            + " rings "
            + " ".join(f"{size}:{int(rings[size])}" for size in range(3, 7))
        )
        self.live.append((frame["labels"].copy(), counts.copy(), rings.copy()))
        self.engine.frames.clear()
        self.log(line)


def crystal_census(library, pos, box):
    engine = InSitu(library, skin=0.0)
    engine.keep_frames = False
    raw = np.zeros(9, dtype=np.float64)
    raw[3], raw[4], raw[5] = box
    _labels, counts, rings = engine.analyse(pos, raw)
    return counts, rings


def main(argv):
    if len(argv) < 4:
        print(__doc__, file=sys.stderr)
        return 2
    library, seams, log_path = argv[1], argv[2], Path(argv[3])
    os.environ["SEAMS_LANG_LIB"] = library
    OUT.mkdir(parents=True, exist_ok=True)
    log_path.parent.mkdir(parents=True, exist_ok=True)
    lines = []

    def say(text):
        print(text, flush=True)
        lines.append(text)
        log_path.write_text("\n".join(lines) + "\n")

    say(f"cores {os.cpu_count()}")
    say(f"affinity {len(os.sched_getaffinity(0))}")
    say(f"cells {CELLS} temperature-K {TEMP_K}")
    say(f"science-steps {STEPS} interval {INTERVAL} threads {SCIENCE_THREADS}")
    say(f"time-steps {TIME_STEPS} time-interval {TIME_INTERVAL}")
    say("platform CPU")
    say(f"library {library}")
    say(f"seams {seams}")
    say(f"frames {OUT}")
    set_classifier_threads(1)

    from openmm.unit import nanometer

    say("build")
    ff, modeller, cutoff, crystal, crystal_box = prepare_box(CELLS)
    counts, rings = crystal_census(library, crystal, crystal_box)
    say(
        "periodic-crystal "
        + " ".join(f"{name} {int(counts[i])}" for i, name in enumerate(_NAMES))
        + " rings "
        + " ".join(f"{size}:{int(rings[size])}" for size in range(3, 7))
    )
    say("minimize")
    simulation = make_simulation(ff, modeller, cutoff, 1)
    simulation.minimizeEnergy(maxIterations=80)
    kinds = oxygen_index(simulation.topology)
    edges = simulation.topology.getUnitCellDimensions().value_in_unit(nanometer)
    say(f"molecules {int(kinds.sum())} atoms {simulation.system.getNumParticles()}")
    say("box-nm " + " ".join(f"{float(v):.6f}" for v in edges))
    say(f"pme-cutoff-nm {cutoff:.6f}")
    state = snapshot(simulation)

    def sim_for(threads):
        if threads == 1:
            return simulation
        built = make_simulation(ff, modeller, cutoff, threads)
        built.context.setPositions(state.getPositions())
        built.context.setVelocities(state.getVelocities())
        return built

    prepared = {}

    def cached(threads):
        if threads not in prepared:
            prepared[threads] = sim_for(threads)
        return prepared[threads]

    say("timing")
    for threads in TIME_THREADS:
        set_classifier_threads(threads)
        sim = cached(threads)
        restore(sim, state)
        bare = run_steps(sim, TIME_STEPS, None)
        restore(sim, state)
        engine = InSitu(library, skin=0.0)
        engine.keep_frames = False
        reporter = DseamsReporter(engine, interval=TIME_INTERVAL, types=kinds)
        spent = run_steps(sim, TIME_STEPS, reporter)
        reports = max(engine.profile["reports"], 1)
        extra = spent - bare
        say(
            f"time threads {threads} bare-s {bare:.6f} bare-ns-day {ns_per_day(TIME_STEPS, bare):.6f} "
            f"reporter-s {spent:.6f} reporter-ns-day {ns_per_day(TIME_STEPS, spent):.6f} "
            f"reports {engine.profile['reports']} ms-per-report {extra / reports * 1000.0:.6f} "
            f"slowdown {extra / bare:.6f}"
        )
        prof = engine.profile
        nrep = max(prof["reports"], 1)
        say(
            "time-split-ms "
            f"{threads} copy {prof['copy'] / nrep * 1000.0:.6f} "
            f"neighbours {prof['neighbours'] / nrep * 1000.0:.6f} "
            f"chill {prof['chill'] / nrep * 1000.0:.6f} "
            f"rings {prof['rings'] / nrep * 1000.0:.6f}"
        )

    say("science")
    set_classifier_threads(SCIENCE_THREADS)
    sim = cached(SCIENCE_THREADS)
    restore(sim, state)
    engine = InSitu(library, skin=0.0)
    dump_path = OUT / "oxygens.lammpstrj"
    with dump_path.open("w") as handle:
        reporter = CensusReporter(
            engine, None, INTERVAL, kinds, say, handle
        )
        t0 = time.perf_counter()
        sim.reporters.clear()
        sim.reporters.append(reporter)
        # Step in chunks so the census hits the log as the trajectory advances.
        done = 0
        chunk = max(INTERVAL, 1)
        while done < STEPS:
            n = min(chunk, STEPS - done)
            sim.step(n)
            done += n
        elapsed = time.perf_counter() - t0
    say(
        f"science-s {elapsed:.6f} science-ns-day {ns_per_day(STEPS, elapsed):.6f} "
        f"frames {reporter.saved}"
    )
    say("offline")
    set_classifier_threads(1)
    checker = InSitu(library, skin=0.0)
    checker.keep_frames = False
    for index, (labels, counts, rings) in enumerate(reporter.live, start=1):
        file_labels, file_counts, file_rings = checker.analyse_file(dump_path, index)
        if not np.array_equal(file_labels, labels):
            raise SystemExit(f"frame {index} CHILL+ labels differ")
        if not np.array_equal(file_counts, counts):
            raise SystemExit(f"frame {index} CHILL+ counts differ")
        if not np.array_equal(file_rings, rings):
            raise SystemExit(f"frame {index} ring counts differ")
        cli, raw = cli_counts(seams, dump_path, index)
        live = histogram(labels)
        for name, number in live.items():
            if number == 0:
                continue
            if cli.get(name) != number:
                raise SystemExit(f"frame {index} CLI {name} {cli.get(name)} != {number}\n{raw}")
        say(f"offline {index} exact")
    say("match exact")
    say("EXIT:0")
    return 0


if __name__ == "__main__":
    sys.exit(main(sys.argv))
