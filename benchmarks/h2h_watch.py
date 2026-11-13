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
    "mdapy-rings",
    "dseams-neighbours",
    "vesin-neighbours",
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
        if arm == "fastatomstruct-rings" and molecules > 50000:
            raise NotImplementedError(
                "King rings on this cell were still running after 35 min")
        self.cmd = [
            "taskset", "-c", "0", PY, str(ROOT / "scripts/h2h_watch.py"),
            "--once", arm, "--molecules", str(molecules), "--lib", LIB,
        ]

    def time_call(self, arm, molecules):
        _run(self.cmd)


# One call per arm. The median of six calls is the Snakemake rule
# related_watch, which runs scripts/h2h_watch.py.
RelatedWatch.time_call.repeat = 1
RelatedWatch.time_call.number = 1
RelatedWatch.time_call.warmup_time = 0


# The published medians are benchmarks/rc3-sitting.csv, from
# snakemake rc3_sitting. This arm is one call of the same worker.
RC3_KINDS = [
    "dseams-neighbours",
    "vesin-neighbours",
    "dseams-q4",
    "mdapy-q4",
    "freud-q4",
    "plumed-q4",
    "dseams-q6",
    "mdapy-q6",
    "freud-q6",
    "plumed-q6",
    "dseams-chill",
    "mdapy-chill",
    "ovito-chill",
    "dseams-rings",
    "ovito-rings",
    "mdapy-rings",
    "ringsstat-rings",
    "dseams-ptm",
    "ovito-ptm",
]
RC3_MOLECULES = [110592, 884736]
RC3_THREADS = [1, 8]


class Rc3Sitting:
    params = [RC3_KINDS, RC3_MOLECULES, RC3_THREADS]
    param_names = ["kind", "molecules", "threads"]

    def setup(self, kind, molecules, threads):
        import sys
        import numpy as np
        sys_path = str(ROOT / "scripts")
        if sys_path not in sys.path:
            sys.path.insert(0, sys_path)
        import h2h_protocol as proto
        repetitions = 24 if molecules == 110592 else 48
        pos, side = proto.ice_ic(repetitions)
        frame = ROOT / "bbdir" / "rc3-asv" / f"ic-{molecules}.npz"
        frame.parent.mkdir(parents=True, exist_ok=True)
        if not frame.is_file():
            np.savez(frame, positions=np.ascontiguousarray(pos), side=np.float64(side))
        cores = "0" if threads == 1 else "0-7"
        self.cmd = [
            "taskset", "-c", cores, PY, str(ROOT / "scripts/h2h_rc3_sitting.py"),
            "--worker", kind, "--molecules", str(molecules),
            "--threads", str(threads), "--frame", str(frame), "--calls", "1",
        ]

    def time_call(self, kind, molecules, threads):
        _run(self.cmd)


Rc3Sitting.time_call.repeat = 1
Rc3Sitting.time_call.number = 1
Rc3Sitting.time_call.warmup_time = 0
