# Copyright (c) d-SEAMS developers

"""Time per-frame interface analysis against a pure-Python q6bar.

One session. The same frames are scored by interface_order growth
(the 2015 and OVITO CHILL+ rules, thresholds unchanged), by
interface_order q6 (Lechner-Dellago q6bar), and by
scripts/q6bar_baseline.py. n replicates, then the median and the
interquartile range. taskset pins the run to 1 core and to 8 cores.
The core count is the size of the affinity mask. No host name is read.

    python3 scripts/q6bar_timing.py --bin bbdir/tests/interface_order
"""

from __future__ import annotations

import argparse
import os
import subprocess
import time
from pathlib import Path

import espinosa_candidates as cand
import q6bar_baseline as baseline

ROOT = Path(__file__).resolve().parents[1]
N_REPEAT = 8
N_FRAMES = 8


def quantile_iqr(samples: list[float]) -> tuple[float, float, float, float]:
    ordered = sorted(samples)
    median = cand.quantile7(ordered, 0.5)
    q1 = cand.quantile7(ordered, 0.25)
    q3 = cand.quantile7(ordered, 0.75)
    return median, q1, q3, q3 - q1


def extract_frames(src: Path, dest: Path, n_frames: int) -> int:
    text = src.read_text().split("ITEM: TIMESTEP\n")
    # split keeps a leading header fragment before the first timestep.
    chunks = [block for block in text[1:] if block.strip()]
    if len(chunks) < n_frames:
        raise SystemExit(f"{src.name} has {len(chunks)} frames, need {n_frames}")
    if len(chunks) == n_frames:
        chosen = list(range(n_frames))
    else:
        # Spread the frames across the trajectory so the sample is not one instant.
        step = (len(chunks) - 1) / (n_frames - 1)
        chosen = [int(round(i * step)) for i in range(n_frames)]
    parts = ["ITEM: TIMESTEP\n" + chunks[i] for i in chosen]
    dest.parent.mkdir(parents=True, exist_ok=True)
    dest.write_text("".join(parts))
    if not dest.read_text().endswith("\n"):
        dest.write_text(dest.read_text() + "\n")
    return len(chosen)


def affinity_cores(mask: str) -> int:
    script = "import os; print(len(os.sched_getaffinity(0)))"
    env = os.environ.copy()
    env["OMP_NUM_THREADS"] = "1"
    done = subprocess.run(
        ["taskset", "-c", mask, "python3", "-c", script],
        check=True,
        capture_output=True,
        text=True,
        env=env,
    )
    return int(done.stdout.strip().splitlines()[-1])


def run_timed(mask: str, cores: int, cmd: list[str]) -> float:
    env = os.environ.copy()
    env["OMP_NUM_THREADS"] = str(cores)
    env["OMP_PROC_BIND"] = "close"
    env["OMP_PLACES"] = "cores"
    # The pure-Python baseline splits frames across this many processes.
    # d-SEAMS reads the same variable and ignores it.
    env["Q6_WORKERS"] = str(cores)
    t0 = time.perf_counter()
    subprocess.run(
        ["taskset", "-c", mask, *cmd],
        check=True,
        capture_output=True,
        env=env,
    )
    return time.perf_counter() - t0


def choose_dump(data: Path) -> Path:
    found = cand.discover(data)
    # Prefer a secondary-prism run at 240 K, the temperature already scored.
    found.sort(key=lambda item: (item[0] != "prism2", item[1] != 240, item[1], item[2]))
    for _face, _temperature, _seed, step, path, _log in found:
        if abs(step - 10.0) < 1.0e-6 and path.stat().st_size > 10000:
            return path
    raise SystemExit(f"no trajectory in {data}")


def summarise(tool: str, cores: int, n_frames: int, samples: list[float]) -> str:
    median, q1, q3, iqr = quantile_iqr(samples)
    per = [sample / n_frames for sample in samples]
    p_med, p_q1, p_q3, p_iqr = quantile_iqr(per)
    return ",".join(
        [
            tool,
            str(cores),
            str(len(samples)),
            str(n_frames),
            f"{median:.6f}",
            f"{q1:.6f}",
            f"{q3:.6f}",
            f"{iqr:.6f}",
            f"{p_med:.6f}",
            f"{p_q1:.6f}",
            f"{p_q3:.6f}",
            f"{p_iqr:.6f}",
        ]
    )


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--bin", type=Path, default=ROOT / "bbdir" / "tests" / "interface_order")
    parser.add_argument("--data", type=Path, default=Path(os.environ.get("SI4_DATA", "si4-data")))
    parser.add_argument(
        "--out",
        type=Path,
        default=ROOT / "docs" / "orgmode" / "explanation" / "growth-curve-timing.csv",
    )
    parser.add_argument("--repeats", type=int, default=N_REPEAT)
    parser.add_argument("--frames", type=int, default=N_FRAMES)
    args = parser.parse_args()
    if args.repeats < 8 or args.frames < 8:
        raise SystemExit("n and the frame count are at least 8")
    src = choose_dump(args.data)
    sample = args.data / "timing" / "sample.lammpstrj"
    n_frames = extract_frames(src, sample, args.frames)
    masks = (("0", 1), ("0-7", 8))
    header = (
        "tool,cores,n,n_frames,median_s,q1_s,q3_s,iqr_s,"
        "per_frame_median_s,per_frame_q1_s,per_frame_q3_s,per_frame_iqr_s"
    )
    rows = [header]
    py = "python3"
    commands = {
        "dseams-growth": [str(args.bin), "growth", str(sample), "5.0"],
        "dseams-q6": [str(args.bin), "q6", str(sample), "5.0"],
        "python-q6bar": [py, str(ROOT / "scripts" / "q6bar_baseline.py"), str(sample)],
    }
    notes = [
        f"# source_frames {src.name}",
        f"# repeats {args.repeats}",
        f"# frames {n_frames}",
    ]
    for mask, cores in masks:
        seen = affinity_cores(mask)
        notes.append(f"# affinity cores {cores} mask {mask} sched {seen}")
        if seen != cores:
            raise SystemExit(f"taskset {mask} gave {seen} cores, expected {cores}")
        for tool, cmd in commands.items():
            samples = []
            for _ in range(args.repeats):
                samples.append(run_timed(mask, cores, cmd))
            rows.append(summarise(tool, cores, n_frames, samples))
            print(rows[-1], flush=True)
    notes.append("# median and IQR are over the repeats. IQR is the type-7 third quartile minus the first.")
    notes.append("# per_frame divides each repeat by the frame count. The same frames are used for every tool.")
    notes.append("# dseams-growth is interface_order growth. dseams-q6 is interface_order q6.")
    notes.append("# python-q6bar is scripts/q6bar_baseline.py, one process, the Lechner-Dellago average.")
    notes.append("# OMP_NUM_THREADS equals the core count. No host name is recorded.")
    args.out.parent.mkdir(parents=True, exist_ok=True)
    args.out.write_text("\n".join(rows + notes) + "\n")
    if args.out.stat().st_size > 80000:
        raise SystemExit(f"{args.out} exceeds 80 KB")
    print(f"wrote {args.out.name} ({args.out.stat().st_size} bytes)")


if __name__ == "__main__":
    main()
