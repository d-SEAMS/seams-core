"""Airspeed Velocity arms for the campaign head-to-heads.

Two head-to-heads, neighbour lists and ring statistics, on the frames and
with the tools of scripts/h2h_campaign.py. The arms are the same adapter
classes the campaign worker uses, so a sample here and a sample there are
the same call.

One sample is one process and one call. Setup imports the tool, builds the
frame and the tool's input object. The timed call is the first call of that
process, with no warm-up before it. Run it with the sampling attributes
number=1, repeat=1, processes=1 and warmup_time=0, one core, and the thread
variables below. ``scripts/h2h_campaign.py lanes`` does that, one invocation
per sample, in interleaved rounds. The Snakemake rules ``h2h_neighbours`` and
``h2h_rings`` run the campaign worker over the same cells.

peakmem_call is the peak resident set of the whole process, imports
included, as Airspeed Velocity reports it.
"""

import os
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]

for _name, _value in (
    ("OMP_NUM_THREADS", "1"),
    ("OMP_DYNAMIC", "FALSE"),
    ("MKL_NUM_THREADS", "1"),
    ("OPENBLAS_NUM_THREADS", "1"),
    ("NUMEXPR_NUM_THREADS", "1"),
    ("RAYON_NUM_THREADS", "1"),
    ("OVITO_THREAD_COUNT", "1"),
    ("QT_QPA_PLATFORM", "offscreen"),
    ("FASTATOMSTRUCT_NO_TELEMETRY", "1"),
    ("FASTATOMSTRUCT_SUPPRESS_TELEMETRY_NOTICE", "1"),
):
    os.environ.setdefault(_name, _value)

sys.path.insert(0, str(ROOT / "scripts"))

import h2h_campaign as campaign  # noqa: E402

SYSTEMS = list(campaign.SYSTEMS)


class _Arm:
    """One cell of a head-to-head: a frame, a workload and an arm."""

    workload = None

    def setup(self, *params):
        workload, system, arm = self._cell(*params)
        if arm not in campaign.ARMS[workload]:
            raise NotImplementedError(f"{arm} has no {workload} arm")
        pos, lengths = campaign.build_frame(system)
        if not campaign.applicable(arm, workload, int(pos.shape[0])):
            raise NotImplementedError(f"{arm} is not run on {system} for {workload}")
        self.workload = workload
        self.adapter = campaign.adapter_for(arm)
        self.adapter.prepare(pos, lengths, workload)
        self.out = None

    def time_call(self, *params):
        self.out = self.adapter.call(self.workload)

    def peakmem_call(self, *params):
        self.out = self.adapter.call(self.workload)

    def teardown(self, *params):
        self.out = None


class NeighbourLists(_Arm):
    """Cutoff pairs and four nearest neighbours of every atom."""

    params = [["cutoff", "knn"], SYSTEMS, list(campaign.ARMS["cutoff"])]
    param_names = ["workload", "system", "arm"]

    def _cell(self, workload, system, arm):
        return workload, system, arm


class RingStatistics(_Arm):
    """Franzblau primitive rings of size 3 to 6 on the 3.5 A graph."""

    params = [SYSTEMS, list(campaign.ARMS["rings"])]
    param_names = ["system", "arm"]

    def _cell(self, system, arm):
        return "rings", system, arm
