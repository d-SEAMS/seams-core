"""Airspeed Velocity arms for the literature head-to-head.

Each timed method is one call. The median, the interquartile range and
the per-atom agreement are the Snakemake rule ``literature_h2h``, which
runs ``scripts/h2h_literature.py`` for eight samples in one session.
Serial is cpu 0. Eight threads are cpus 0-7.
"""

import importlib.util
import os
import subprocess
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
PY = os.environ.get("H2H_PYTHON", sys.executable)
LIB = os.environ.get(
    "SEAMS_CORE_LIB", str(ROOT / "bbdir-release" / "src" / "libyodaLib.so"))
ENV = os.environ.copy()
ENV.update({
    "OMP_DYNAMIC": "FALSE",
    "OPENBLAS_NUM_THREADS": "1",
    "MKL_NUM_THREADS": "1",
    "NUMEXPR_NUM_THREADS": "1",
})


def _run(cmd):
    subprocess.run(
        cmd, cwd=ROOT, env=ENV, check=True,
        stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL,
    )


def _need(module=None):
    """A missing library or tool skips the arm instead of failing the sample."""
    if not Path(LIB).is_file():
        raise NotImplementedError(f"{Path(LIB).name} is not built")
    if module and importlib.util.find_spec(module) is None:
        raise NotImplementedError(f"{module} is not installed")


def _command(code, operation, threads, molecules):
    return [
        PY, str(ROOT / "scripts" / "h2h_literature.py"),
        "--once", code, operation,
        "--threads", str(threads),
        "--molecules", str(molecules),
        "--phase", "Ic",
        "--lib", LIB,
    ]


class LiteratureMdapy:
    params = [["neighbours", "q4", "q6", "diamond", "cna", "rings"], [1, 8], [128, 4096, 110592]]
    param_names = ["operation", "threads", "molecules"]

    def setup(self, operation, threads, molecules):
        _need("mdapy")
        self.cmd = _command("mdapy", operation, threads, molecules)

    def time_call(self, operation, threads, molecules):
        _run(self.cmd)


class LiteraturePyscal3:
    params = [["neighbours", "q4", "q6", "diamond", "cna"], [1, 8], [128, 4096, 110592]]
    param_names = ["operation", "threads", "molecules"]

    def setup(self, operation, threads, molecules):
        _need("pyscal3")
        self.cmd = _command("pyscal3", operation, threads, molecules)

    def time_call(self, operation, threads, molecules):
        _run(self.cmd)


class LiteratureVesin:
    params = [[1, 8], [128, 4096, 110592]]
    param_names = ["threads", "molecules"]

    def setup(self, threads, molecules):
        _need("vesin")
        self.cmd = _command("vesin", "neighbours", threads, molecules)

    def time_call(self, threads, molecules):
        _run(self.cmd)


class LiteratureIra:
    params = [[1, 8], [128, 4096, 110592]]
    param_names = ["threads", "molecules"]

    def setup(self, threads, molecules):
        _need()
        self.cmd = _command("IRA", "symmetry", threads, molecules)

    def time_call(self, threads, molecules):
        _run(self.cmd)


class LiteratureSqq:
    params = [[1, 8], [128, 4096, 110592]]
    param_names = ["threads", "molecules"]

    def setup(self, threads, molecules):
        _need("sqq")
        self.cmd = _command("SQQ", "rings", threads, molecules)

    def time_call(self, threads, molecules):
        _run(self.cmd)


class LiteratureClavis:
    """TRACE cages on the GenIce sI (46) and sII (136) cells, not on ice."""

    params = [[1, 8], [46, 136]]
    param_names = ["threads", "molecules"]

    def setup(self, threads, molecules):
        _need("clavis_engine")
        phase = "sI" if molecules == 46 else "sII"
        self.cmd = [
            PY, str(ROOT / "scripts" / "h2h_literature.py"),
            "--once", "ClavisCage", "cages",
            "--threads", str(threads),
            "--molecules", str(molecules),
            "--phase", phase,
            "--lib", LIB,
        ]

    def time_call(self, threads, molecules):
        _run(self.cmd)


class LiteraturePlumed:
    params = [["q4", "q6", "coordination"], [1, 8], [128, 4096, 110592]]
    param_names = ["operation", "threads", "molecules"]

    def setup(self, operation, threads, molecules):
        _need("plumed")
        self.cmd = _command("PLUMED", operation, threads, molecules)

    def time_call(self, operation, threads, molecules):
        _run(self.cmd)
