# Copyright (c) d-SEAMS developers
"""Serial wall time of plain OpenMM against the reporter.

Two TIP4P/Ew boxes, report intervals 10, 100 and 1000, n paired trials.
One warmup of each leg is discarded per box. Trials run one after another
on one thread. The overhead of a trial is the in-situ time minus the plain
time.
"""

import os
import statistics
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT / "adapters" / "python"))
sys.path.insert(0, str(ROOT / "engines" / "openmm"))

from dseams_insitu.engine import InSitu  # noqa: E402
from dseams_insitu.openmm import DseamsReporter  # noqa: E402

from water import build_simulation, restore, save_state  # noqa: E402

LIB = Path(os.environ.get("SEAMS_LANG_LIB", "libseams_lang"))
BOXES_NM = (2.0, 4.93)
INTERVALS = (10, 100, 1000)
STEPS = 1000
N = 6


def hinges(values):
    """Median, and the Tukey hinge (upper-half median minus lower-half median)."""
    ordered = sorted(values)
    half = len(ordered) // 2
    lower = ordered[:half]
    upper = ordered[half:] if len(ordered) % 2 == 0 else ordered[half + 1 :]
    return statistics.median(ordered), statistics.median(upper) - statistics.median(lower)


def time_plain(simulation, state, steps):
    simulation.reporters.clear()
    restore(simulation, state)
    import time

    start = time.perf_counter()
    simulation.step(steps)
    return time.perf_counter() - start


def time_insitu(simulation, state, lib, interval, steps):
    import time

    simulation.reporters.clear()
    restore(simulation, state)
    engine = InSitu(library=str(lib), cutoff=3.5)
    simulation.reporters.append(DseamsReporter(engine, "unused", interval=interval))
    start = time.perf_counter()
    simulation.step(steps)
    elapsed = time.perf_counter() - start
    simulation.reporters.clear()
    return elapsed, len(engine.frames)


def run_box(box_nm, lib, steps, intervals, trials):
    simulation, _oxygens = build_simulation(box_nm, threads=1)
    waters = simulation.topology.getNumResidues()
    state = save_state(simulation)
    time_plain(simulation, state, steps)
    time_insitu(simulation, state, lib, intervals[-1], steps)
    rows = []
    for interval in intervals:
        plain = []
        insitu = []
        for _ in range(trials):
            plain.append(time_plain(simulation, state, steps))
            elapsed, _reports = time_insitu(simulation, state, lib, interval, steps)
            insitu.append(elapsed)
        overhead = [b - a for a, b in zip(plain, insitu)]
        median, hinge = hinges(overhead)
        rows.append((interval, median, hinge))
    return waters, rows


def main():
    lib = LIB
    if not lib.is_file() and str(lib) == "libseams_lang":
        built = ROOT / "build" / "bindings" / "libseams_lang.so"
        if built.is_file():
            lib = built
    for box in BOXES_NM:
        waters, rows = run_box(box, lib, STEPS, INTERVALS, N)
        print(f"waters {waters} box_nm {box}")
        for interval, median, hinge in rows:
            print(f"interval {interval} overhead_s {median:.6f} hinge_s {hinge:.6f}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
