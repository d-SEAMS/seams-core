#!/usr/bin/env python3
"""One-thread and eight-thread ring search against OVITO FindRings.

The cutoff graph is built before the timer on both sides. d-SEAMS times
primitive::ringNetwork. OVITO times FindRingsModifier on bonds from
CreateBondsModifier(cutoff). A fresh OVITO pipeline is used for every
sample so the modifier does not return a cached frame.

Eight timed calls after two warmups, interleaved, in one driver run.
The median of an even sample is the mean of the two central values.
The record has the core count and the CPU model, not a host name.
"""

from __future__ import annotations

import argparse
import importlib.util
import json
import os
import subprocess
import sys
import tempfile
from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]
WARMUP = 2
REPEATS = 8
CUTOFF = 3.5
MAX_RING = 6


def load_h2h():
    path = ROOT / "tests" / "bench_h2h_rings.py"
    spec = importlib.util.spec_from_file_location("bench_h2h_rings", path)
    module = importlib.util.module_from_spec(spec)
    sys.modules["bench_h2h_rings"] = module
    spec.loader.exec_module(module)
    return module


def cpu_model() -> str:
    try:
        text = Path("/proc/cpuinfo").read_text(encoding="utf-8")
    except OSError:
        return "unknown"
    for line in text.splitlines():
        if line.startswith("model name"):
            name = line.split(":", 1)[1].strip()
            return " ".join(name.replace("(R)", "").replace("(r)", "").split())
    return "unknown"


def median(samples: list[float]) -> float:
    ordered = sorted(samples)
    mid = len(ordered) // 2
    if len(ordered) % 2:
        return ordered[mid]
    return 0.5 * (ordered[mid - 1] + ordered[mid])


def parse_hist(text: str) -> dict[int, int]:
    hist: dict[int, int] = {}
    if text.strip() == "none":
        return hist
    for item in text.split():
        size, count = item.split(":")
        hist[int(size)] = int(count)
    return hist


class DseamsRepl:
    def __init__(self, bench: str, traj: str, threads: int):
        cores = "0" if threads == 1 else "0-7"
        env = os.environ.copy()
        env.pop("LD_LIBRARY_PATH", None)
        env["OMP_NUM_THREADS"] = str(threads)
        env["OMP_DYNAMIC"] = "FALSE"
        env["OMP_PROC_BIND"] = "close"
        self.proc = subprocess.Popen(
            ["taskset", "-c", cores, bench, "--repl", traj, "1", "1",
             f"{CUTOFF:.8f}", str(MAX_RING)],
            stdin=subprocess.PIPE,
            stdout=subprocess.PIPE,
            text=True,
            env=env,
        )
        assert self.proc.stdout is not None
        line = self.proc.stdout.readline()
        if not line.startswith("ready "):
            raise RuntimeError(f"bench did not start: {line!r}")
        parts = line.split()
        self.nop = int(parts[parts.index("nop") + 1])
        hist_at = parts.index("hist")
        self.hist = parse_hist(" ".join(parts[hist_at + 1:]))

    def sample(self) -> float:
        assert self.proc.stdin is not None and self.proc.stdout is not None
        self.proc.stdin.write("go\n")
        self.proc.stdin.flush()
        line = self.proc.stdout.readline().strip()
        if not line.startswith("ms "):
            raise RuntimeError(f"bad sample line {line!r}")
        return float(line.split()[1])

    def close(self) -> None:
        if self.proc.stdin is not None:
            self.proc.stdin.write("quit\n")
            self.proc.stdin.flush()
        self.proc.wait(timeout=30)


class OvitoRepl:
    """Child interpreter so OVITO_THREAD_COUNT is read at import."""

    def __init__(self, python: str, traj: str, threads: int):
        cores = "0" if threads == 1 else "0-7"
        env = os.environ.copy()
        env["OVITO_THREAD_COUNT"] = str(threads)
        env["OMP_NUM_THREADS"] = str(threads)
        self.proc = subprocess.Popen(
            ["taskset", "-c", cores, python, str(Path(__file__).resolve()),
             "--ovito-worker", traj],
            stdin=subprocess.PIPE,
            stdout=subprocess.PIPE,
            text=True,
            env=env,
        )
        assert self.proc.stdout is not None
        line = self.proc.stdout.readline()
        if not line.startswith("ready "):
            err = ""
            if self.proc.stderr is not None:
                err = self.proc.stderr.read()
            raise RuntimeError(f"ovito worker failed: {line!r} {err}")
        payload = json.loads(line[len("ready "):])
        self.nop = int(payload["nop"])
        self.hist = {int(k): int(v) for k, v in payload["hist"].items()}
        self.version = str(payload["version"])

    def sample(self) -> float:
        assert self.proc.stdin is not None and self.proc.stdout is not None
        self.proc.stdin.write("go\n")
        self.proc.stdin.flush()
        line = self.proc.stdout.readline().strip()
        if not line.startswith("ms "):
            raise RuntimeError(f"bad ovito sample {line!r}")
        return float(line.split()[1])

    def close(self) -> None:
        if self.proc.stdin is not None:
            self.proc.stdin.write("quit\n")
            self.proc.stdin.flush()
        self.proc.wait(timeout=60)


def ovito_worker_main(traj: str) -> None:
    from ovito.data import DataCollection
    from ovito.modifiers import CreateBondsModifier, FindRingsModifier
    from ovito.pipeline import Pipeline, StaticSource
    import ovito

    h2h = load_h2h()
    frame = h2h.read_lammps_frame(traj, 1)
    data = DataCollection()
    data.create_cell(frame.H, pbc=(True, True, True))
    parts = data.create_particles(count=frame.positions.shape[0])
    parts.create_property("Position", data=frame.positions - frame.origin)
    bond_pipe = Pipeline(source=StaticSource(data=data))
    bond_pipe.modifiers.append(CreateBondsModifier(cutoff=CUTOFF))
    bonded = bond_pipe.compute()
    probe = Pipeline(source=StaticSource(data=bonded))
    probe.modifiers.append(FindRingsModifier(maximum_ring_size=MAX_RING))
    out = probe.compute()
    table = out.tables["ring-size-histogram"]
    hist = {
        str(int(size)): int(count)
        for size, count in zip(table["Ring Size"][...], table["Count"][...])
        if int(count) != 0
    }
    sys.stdout.write("ready " + json.dumps({
        "nop": int(frame.positions.shape[0]),
        "hist": hist,
        "version": ".".join(str(part) for part in ovito.version),
    }) + "\n")
    sys.stdout.flush()
    for line in sys.stdin:
        if line.startswith("q"):
            break
        if not line.startswith("g"):
            continue
        ring_pipe = Pipeline(source=StaticSource(data=bonded))
        ring_pipe.modifiers.append(FindRingsModifier(maximum_ring_size=MAX_RING))
        import time
        t0 = time.perf_counter()
        ring_pipe.compute()
        t1 = time.perf_counter()
        sys.stdout.write(f"ms {(t1 - t0) * 1000.0:.6f}\n")
        sys.stdout.flush()


def prepare_frames(h2h, work: Path) -> list[tuple[str, str, int]]:
    small = ROOT / "input" / "traj" / "genice_sH.lammpstrj"
    cubic_path = ROOT / "input" / "traj" / "mW_cubic.lammpstrj"
    cubic = h2h.read_lammps_frame(str(cubic_path), 1)
    replica = h2h.replicate_frame(cubic, 3, 3, 3)
    replica_path = work / "mW-cubic-3x3x3.lammpstrj"
    h2h.write_lammps_dump(str(replica_path), replica, 1)
    return [
        ("tri-136", str(small), 136),
        ("orth-4k", str(cubic_path), 4096),
        ("orth-1e5", str(replica_path), 110592),
    ]


def run_dataset(bench: str, python: str, name: str, traj: str, molecules: int,
                threads: int) -> dict:
    dseams = DseamsRepl(bench, traj, threads)
    ovito = OvitoRepl(python, traj, threads)
    try:
        if dseams.nop != molecules or ovito.nop != molecules:
            raise SystemExit(
                f"{name} molecule count d-SEAMS {dseams.nop} OVITO {ovito.nop}"
            )
        if dseams.hist != ovito.hist:
            raise SystemExit(
                f"{name} census differs d-SEAMS {dseams.hist} OVITO {ovito.hist}"
            )
        samples = {"dseams": [], "ovito": []}
        for i in range(WARMUP + REPEATS):
            order = ("dseams", "ovito") if i % 2 == 0 else ("ovito", "dseams")
            for tool in order:
                ms = dseams.sample() if tool == "dseams" else ovito.sample()
                if i >= WARMUP:
                    samples[tool].append(ms)
        return {
            "molecules": molecules,
            "census": {str(k): v for k, v in sorted(dseams.hist.items())},
            "dseams_ms": samples["dseams"],
            "ovito_ms": samples["ovito"],
            "dseams_median_ms": median(samples["dseams"]),
            "ovito_median_ms": median(samples["ovito"]),
            "ovito_version": ovito.version,
        }
    finally:
        dseams.close()
        ovito.close()


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--bench", required=False)
    ap.add_argument("--python", default=sys.executable)
    ap.add_argument("--json", dest="json_path")
    ap.add_argument("--ovito-worker", metavar="TRAJ")
    args = ap.parse_args()
    if args.ovito_worker:
        ovito_worker_main(args.ovito_worker)
        return
    if not args.bench or not args.json_path:
        raise SystemExit("--bench and --json are required")

    h2h = load_h2h()
    record: dict = {
        "nproc": os.cpu_count() or 1,
        "cpu_model": cpu_model(),
        "cutoff_A": CUTOFF,
        "max_ring_size": MAX_RING,
        "repeats": REPEATS,
        "warmup": WARMUP,
        "median": "mean of the two central samples",
        "timing": (
            "cutoff graph built before the timer; "
            "d-SEAMS ringNetwork; OVITO FindRings on a fresh pipeline"
        ),
        "threads": {},
    }
    with tempfile.TemporaryDirectory(prefix="rings-h2h-") as work:
        frames = prepare_frames(h2h, Path(work))
        for threads in (1, 8):
            arm = {}
            for name, traj, molecules in frames:
                print(f"[rings] {name} threads {threads}", flush=True)
                arm[name] = run_dataset(
                    args.bench, args.python, name, traj, molecules, threads)
                row = arm[name]
                print(
                    f"[rings] {name} t{threads} "
                    f"dseams {row['dseams_median_ms']:.4f} "
                    f"ovito {row['ovito_median_ms']:.4f} "
                    f"census {row['census']}",
                    flush=True,
                )
            record["threads"][str(threads)] = arm
    text = json.dumps(record, indent=2) + "\n"
    Path(args.json_path).write_text(text, encoding="utf-8")
    print(f"Wrote {args.json_path}")


if __name__ == "__main__":
    main()
