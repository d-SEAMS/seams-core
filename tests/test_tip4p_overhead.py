# Copyright (c) d-SEAMS developers
"""The TIP4P overhead script's hinge, and one reporter step on a small box."""

import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "adapters" / "python"))
sys.path.insert(0, str(ROOT / "engines" / "openmm"))

from bench_overhead import hinges, time_insitu, time_plain  # noqa: E402
from dseams_insitu.engine import InSitu  # noqa: E402
from water import build_simulation, save_state  # noqa: E402


def main() -> int:
    median, hinge = hinges([1.0, 2.0, 3.0, 4.0])
    if median != 2.5 or abs(hinge - 2.0) > 1e-12:
        print(f"hinges {median} {hinge}", file=sys.stderr)
        return 1
    lib = os_lib()
    simulation, _oxygens = build_simulation(2.0, threads=1)
    if simulation.topology.getNumResidues() < 8:
        print("box has too few waters", file=sys.stderr)
        return 1
    state = save_state(simulation)
    plain = time_plain(simulation, state, 1)
    elapsed, reports = time_insitu(simulation, state, lib, 1, 2)
    if plain < 0.0 or elapsed < 0.0 or reports < 1:
        print(f"plain {plain} insitu {elapsed} reports {reports}", file=sys.stderr)
        return 1
    print(f"tip4p waters {simulation.topology.getNumResidues()} reports {reports}")
    return 0


def os_lib() -> str:
    import os

    env = os.environ.get("SEAMS_LANG_LIB")
    if env:
        return env
    built = ROOT / "build" / "bindings" / "libseams_lang.so"
    if built.is_file():
        return str(built)
    return "libseams_lang"


if __name__ == "__main__":
    raise SystemExit(main())
