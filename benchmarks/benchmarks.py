"""Airspeed Velocity benchmarks for the paper's timed arms.

Imports and frame setup run outside the timed region. Each timed method
is one serial call of the arm on one core.
"""

import os
import subprocess
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
CORE = ["taskset", "-c", "0"]
ENV = os.environ.copy()
ENV.update({
    "OMP_NUM_THREADS": "1",
    "OMP_DYNAMIC": "FALSE",
    "OVITO_THREAD_COUNT": "1",
    "OPENBLAS_NUM_THREADS": "1",
})


def _run(cmd):
    subprocess.run(
        cmd,
        cwd=ROOT,
        env=ENV,
        check=True,
        stdout=subprocess.DEVNULL,
        stderr=subprocess.DEVNULL,
    )


class NeighbourArm:
    def setup(self):
        self.cmd = CORE + [
            str(ROOT / "bbdir/tests/bench_neigh_only"),
            "--traj", "input/traj/genice_sH.lammpstrj",
            "--frame", "1",
            "--type", "1",
            "--reps", "1",
            "--warmup", "0",
            "--cutoff", "3.5",
            "--k", "4",
            "--threads", "1",
        ]

    def time_mutual_neighbours(self):
        _run(self.cmd)


class RingArm:
    def setup(self):
        self.cmd = CORE + [
            str(ROOT / "bbdir/tests/bench_stages"),
            "input/traj/genice_sH.lammpstrj",
            "1",
            "1",
            "1",
        ]

    def time_primitive_rings(self):
        _run(self.cmd)


class CampaignArms:
    def setup(self):
        bench = str(ROOT / "bbdir/tests/bench_paper")
        frame = "input/traj/genice_sI.lammpstrj"
        modes = {
            "chill": "chill",
            "q6": "q6-full",
            "seeded": "seeded-adj",
            "rings": "rings",
        }
        self.commands = {
            name: CORE + [bench, frame, "1", "0", "1", mode]
            for name, mode in modes.items()
        }

    def time_chill(self):
        _run(self.commands["chill"])

    def time_q6(self):
        _run(self.commands["q6"])

    def time_seeded_cages(self):
        _run(self.commands["seeded"])

    def time_rings(self):
        _run(self.commands["rings"])


def _scale_plan():
    import importlib.util
    path = ROOT / "scripts" / "scale_campaign.py"
    spec = importlib.util.spec_from_file_location("scale_campaign", path)
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


_SCALE = _scale_plan()


class ScaleProtocol:
    """Twelve-call strong and weak points. Airspeed Velocity runs each once."""

    # One axis. A list of points is not a Cartesian product of those points.
    params = [[
        (point["kind"], point["model"], point["analysis"], point["threads"],
         point["reps"][0], point["reps"][1], point["reps"][2])
        for point in _SCALE.points()
    ]]
    param_names = ["point"]

    def setup(self, point):
        kind, model, analysis, threads, nx, ny, nz = point
        spec = _SCALE.MODELS[model]
        frame = _SCALE.dump_path(model, (nx, ny, nz))
        if not frame.is_file():
            _SCALE.prepare_dump(spec, model, (nx, ny, nz))
        cpus = "0" if threads == 1 else f"0-{threads - 1}"
        self.cmd = [
            "taskset", "-c", cpus, str(_SCALE.BENCH), "run",
            "--traj", str(frame), "--type", spec["type"],
            "--analysis", analysis, "--threads", str(threads),
            "--reps", str(_SCALE.REPS), "--out-dir", "bbdir/scale/asv-tmp",
        ]
        self.log_name = (
            f"{kind}-{model}-{analysis}-t{threads}-{nx}x{ny}x{nz}.log"
        )
        self.env = os.environ.copy()
        self.env.update({
            "OMP_NUM_THREADS": str(threads),
            "OMP_DYNAMIC": "FALSE",
            "OMP_PROC_BIND": "close",
            "RAYON_NUM_THREADS": str(threads),
            "OPENBLAS_NUM_THREADS": "1",
            "MKL_NUM_THREADS": "1",
        })

    def time_twelve_calls(self, point):
        head = subprocess.check_output(
            ["git", "rev-parse", "HEAD"], cwd=ROOT, text=True).strip()
        tree = subprocess.check_output(
            ["git", "rev-parse", f"{head}^{{tree}}"], cwd=ROOT, text=True).strip()
        shown = " ".join(self.cmd)
        done = subprocess.run(
            self.cmd, cwd=ROOT, env=self.env, check=False,
            stdout=subprocess.PIPE, stderr=subprocess.STDOUT, text=True)
        log_dir = ROOT / "bbdir" / "scale" / "asv-logs"
        log_dir.mkdir(parents=True, exist_ok=True)
        text = "\n".join([
            f"HEAD {head}",
            f"TREE {tree}",
            f"COMMAND {shown}",
            done.stdout.rstrip("\n"),
            f"EXIT {done.returncode}",
            "",
        ])
        (log_dir / self.log_name).write_text(text, encoding="utf-8")
        if done.returncode != 0:
            raise subprocess.CalledProcessError(done.returncode, self.cmd, done.stdout)


ScaleProtocol.time_twelve_calls.repeat = 1
ScaleProtocol.time_twelve_calls.number = 1
ScaleProtocol.time_twelve_calls.warmup_time = 0
