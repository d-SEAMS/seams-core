# Copyright (c) d-SEAMS developers
"""The TIP4P overhead script's hinge, and one reporter step on a small box."""

import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "adapters" / "python"))
sys.path.insert(0, str(ROOT / "engines" / "openmm"))


def main() -> int:
    try:
        from bench_overhead import hinges, time_insitu, time_plain  # noqa: E402
        from dseams_insitu.engine import InSitu  # noqa: F401, E402
        from water import build_simulation, save_state  # noqa: E402
    except ModuleNotFoundError as exc:
        if getattr(exc, "name", "") == "openmm" or "openmm" in str(exc):
            print("openmm is not installed", file=sys.stderr)
            return 77
        raise
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

    sys.path.insert(0, str(ROOT / "scripts"))
    from shared_library import built_library

    env = os.environ.get("SEAMS_LANG_LIB")
    if env:
        return env
    for directory in (ROOT / "build" / "bindings", ROOT / "bbdir" / "bindings"):
        found = built_library(directory, "seams_lang")
        if found is not None:
            return str(found)
    return "libseams_lang"


if __name__ == "__main__":
    raise SystemExit(main())
