"""seams_rings_find on 4096-molecule ideal ice Ih.

One call per method. The Snakemake rule ``rings_ih`` records the median
of eight calls, at one thread and at eight, in one process. The line
it writes names the core count and does not name the host.
"""

import os
import subprocess
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
ENV = os.environ.copy()
ENV["OMP_DYNAMIC"] = "FALSE"
ENV["OPENBLAS_NUM_THREADS"] = "1"


def _bench():
    built = ROOT / "bbdir-release" / "tests" / "bench_rings_ih"
    if built.is_file():
        return built
    return ROOT / "build" / "tests" / "bench_rings_ih"


def _run(threads):
    env = ENV.copy()
    env["OMP_NUM_THREADS"] = str(threads)
    subprocess.run(
        [
            str(_bench()),
            "--molecules",
            "4096",
            "--threads",
            str(threads),
            "--reps",
            "1",
            "--warmup",
            "0",
        ],
        cwd=ROOT,
        env=env,
        check=True,
        stdout=subprocess.DEVNULL,
        stderr=subprocess.DEVNULL,
    )


class RingsIh:
    def setup(self):
        if not _bench().is_file():
            raise RuntimeError(f"missing {_bench()}")

    def time_one_thread(self):
        _run(1)

    def time_eight_threads(self):
        _run(8)
