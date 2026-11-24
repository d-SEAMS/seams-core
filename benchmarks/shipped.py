"""Airspeed Velocity arms for the release comparison.

Each timed method is one shipped command or one bench_stages process on one
core. Peak memory is the resident set Airspeed Velocity samples around that
command. The stage clocks, the type 7 quartiles and the 1e6 abort are the
Snakemake rule ``shipped_campaign``.
"""

import os
import subprocess
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
ENV = os.environ.copy()
ENV.update({
    "OMP_NUM_THREADS": "1",
    "OMP_DYNAMIC": "FALSE",
    "OVITO_THREAD_COUNT": "1",
})


def _run(cmd):
    subprocess.run(cmd, cwd=ROOT, env=ENV, check=True)


def _need(path, what):
    """A missing frame or binary skips the arm. seams cages exits zero without a frame."""
    target = Path(path)
    if not target.is_absolute():
        target = ROOT / target
    if not target.is_file():
        raise NotImplementedError(f"{what} is not available: {path}")
    return str(target)


class ShippedOrth:
    """110592-molecule ice. ``graph`` is cutoff or the default seeded graph."""

    params = ["cutoff", "seeded"]
    param_names = ["graph"]

    def setup(self, graph):
        seams = _need(Path("bbdir-release") / "seams", "the release seams command")
        frame = _need("input/traj/orth-1e5-ring-scale.lammpstrj", "the 110592-molecule frame")
        args = ["cages"] if graph == "seeded" else ["cages", "--graph", "cutoff"]
        self.cmd = ["taskset", "-c", "0", seams, *args, frame]

    def time_cages(self, graph):
        _run(self.cmd)

    def peakmem_cages(self, graph):
        _run(self.cmd)


class ShippedStages:
    """bench_stages process. Its stdout splits rings from neighbours+rings."""

    def setup(self):
        self.cmd = [
            "taskset", "-c", "0",
            _need(Path("bbdir-release") / "tests" / "bench_stages", "release bench_stages"),
            _need("input/traj/orth-1e5-ring-scale.lammpstrj", "the 110592-molecule frame"),
            "1", "1", "1",
        ]

    def time_process(self):
        _run(self.cmd)

    def peakmem_process(self):
        _run(self.cmd)


class ShippedMillion:
    """995328-molecule ice. The base build is not this checkout."""

    params = ["cutoff", "seeded"]
    param_names = ["graph"]

    def setup(self, graph):
        frame = _need(os.environ.get("SEAMS_1E6", "v37-frames/ice-1e6.lammpstrj"),
                      "the 995328-molecule frame")
        seams = _need(Path("bbdir-release") / "seams", "the release seams command")
        args = ["cages"] if graph == "seeded" else ["cages", "--graph", "cutoff"]
        self.cmd = ["taskset", "-c", "0", seams, *args, frame]

    def time_cages(self, graph):
        _run(self.cmd)

    def peakmem_cages(self, graph):
        _run(self.cmd)
