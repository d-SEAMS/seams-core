"""Airspeed Velocity arms for the related-work watch.

Each timed method is one serial call on core 0. The median, the
interquartile range and the label agreement are the Snakemake rule
``related_watch``, which runs ``scripts/h2h_watch.py``. Set H2H_PYTHON
to a interpreter that has mdapy and fastatomstruct, and H2H_LIB to
libseams_lang.so.
"""

import os
import subprocess
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
PY = os.environ.get("H2H_PYTHON", "python3")
LIB = os.environ.get(
    "H2H_LIB", str(ROOT / "bbdir-release/bindings/libseams_lang.so"))
ENV = os.environ.copy()
ENV.update({
    "OMP_NUM_THREADS": "1",
    "OMP_DYNAMIC": "FALSE",
    "OPENBLAS_NUM_THREADS": "1",
    "MKL_NUM_THREADS": "1",
    "RAYON_NUM_THREADS": "1",
    "FASTATOMSTRUCT_NO_TELEMETRY": "1",
})

ARMS = [
    "dseams-chill",
    "mdapy-chill",
    "dseams-q6",
    "fastatomstruct-q6",
    "dseams-rings",
    "fastatomstruct-rings",
]
MOLECULES = [1000, 10648, 110592]


def _run(cmd):
    subprocess.run(
        cmd, cwd=ROOT, env=ENV, check=True,
        stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL,
    )


class RelatedWatch:
    params = [ARMS, MOLECULES]
    param_names = ["arm", "molecules"]

    def setup(self, arm, molecules):
        self.cmd = [
            "taskset", "-c", "0", PY, str(ROOT / "scripts/h2h_watch.py"),
            "--once", arm, "--molecules", str(molecules), "--lib", LIB,
        ]

    def time_call(self, arm, molecules):
        _run(self.cmd)
