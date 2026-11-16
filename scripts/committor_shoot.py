#!/usr/bin/env python3
"""Committor thresholds for the 4096-molecule mW lane.

The commit and uncommit cuts are fixed here, before any trajectory is
read. A shot commits when the largest cluster reaches at least twice
the shooting size and at least 50 molecules. It uncommits when the
largest cluster falls to the liquid baseline, the 95th percentile of
the largest cluster in pure liquid at the same temperature. A shot
that has done neither after 1 ns is unresolved.

p_B(n) is the committed fraction among resolved shots that started at
size n. n* is the size where a logistic fit of p_B crosses one half.
The 64-molecule deck is a fixture for the driver, not this table.

One VM runs one shard of the climbs. Shard k of K is the seeds whose
index modulo K equals k, with k starting at 0:

    export LMP=/path/to/lmp
    export SEAMS_LIB=/path/to/lib
    snakemake -j1 --config seeds=20261111-20261130 shard=0/20 \\
        data/committor-shards/climb-0-of-20.csv

The matching shot shard is data/committor-shards/shot-0-of-20.csv with
the same --config. science_committor_merge reads every shot shard and
writes p_B(n) and n*. Until those shards exist, the paper tables stay
the single line 'pending: 4096-molecule p_B(n)'.
"""

from __future__ import annotations

import argparse
import csv
import math
import os
import random
import subprocess
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]

# Cuts chosen before any shot outcome is known.
COMMIT_MULTIPLE = 2.0
COMMIT_FLOOR = 50
LIQUID_PERCENTILE = 95.0
SHOT_NS = 1.0
N_SHOOT_SIZES = 5
N_SHOTS = 20
BOOT_DRAWS = 10000
BOOT_SEED = 20261111

# Published homogeneous lane, input/nucleation/in.homogeneous.
# 10 fs steps. 10 ns is 1e6 steps. One shot of 1 ns is 1e5 steps.
TIMESTEP_FS = 10.0
CLIMB_STEPS = 1_000_000
SHOT_STEPS = 100_000
N_SEEDS = 20


def commit_size(n_shoot: float) -> float:
    return max(COMMIT_MULTIPLE * n_shoot, COMMIT_FLOOR)


def type7(values: list[float], p: float) -> float:
    ordered = sorted(values)
    n = len(ordered)
    if n == 0:
        raise ValueError("empty sample")
    if n == 1:
        return ordered[0]
    h = (n - 1) * (p / 100.0 if p > 1.0 else p)
    # p may be 95 or 0.95. Percentile helpers below pass 0..1.
    lo = int(math.floor(h))
    hi = int(math.ceil(h))
    if hi >= n:
        hi = n - 1
    weight = h - lo
    return ordered[lo] * (1.0 - weight) + ordered[hi] * weight


def liquid_baseline(nmax: list[float]) -> float:
    """95th percentile of the largest cluster in pure liquid."""
    if not nmax:
        raise ValueError("liquid sample is empty")
    ordered = sorted(nmax)
    h = (len(ordered) - 1) * (LIQUID_PERCENTILE / 100.0)
    lo = int(math.floor(h))
    hi = int(math.ceil(h))
    weight = h - lo
    return ordered[lo] * (1.0 - weight) + ordered[hi] * weight


def shot_fate(series: list[float], n_shoot: float, liquid_p95: float) -> str:
    """Walk one shot in time. The cuts do not depend on later shots."""
    goal = commit_size(n_shoot)
    for size in series:
        if size >= goal:
            return "committed"
        if size <= liquid_p95:
            return "uncommitted"
    return "unresolved"


def check_nmax(nmax: float, molecules: int) -> None:
    if nmax > molecules:
        raise ValueError(f"nmax {nmax} is larger than N {molecules}")


def shooting_sizes(series: list[float], n_sizes: int = N_SHOOT_SIZES) -> list[int]:
    """Sizes around the plateau that precedes the final climb.

    The plateau is the run of samples immediately before the final
    strict ascent whose values stay within 20 percent of their median.
    The returned sizes are evenly spaced across that run. They are not
    fit to a committor.
    """
    if len(series) < 3 or n_sizes < 1:
        return []
    imax = max(range(len(series)), key=lambda i: series[i])
    foot = imax
    while foot > 0 and series[foot - 1] < series[foot]:
        foot -= 1
    window = series[:foot] if foot > 0 else series[:1]
    if len(window) >= 3:
        med = sorted(window)[len(window) // 2]
        band = [value for value in window if abs(value - med) <= 0.2 * max(med, 1.0)]
        if len(band) >= 3:
            window = band
    if len(window) == 1:
        return [int(round(window[0]))]
    lo = min(window)
    hi = max(window)
    if hi <= lo:
        return [int(round(lo))]
    raw = [lo + (hi - lo) * i / (n_sizes - 1) for i in range(n_sizes)]
    sizes = []
    for value in raw:
        size = int(round(value))
        if size not in sizes:
            sizes.append(size)
    return sizes


def logistic(n: float, n_star: float, width: float) -> float:
    if width <= 0.0:
        return 1.0 if n >= n_star else 0.0
    return 1.0 / (1.0 + math.exp(-(n - n_star) / width))


def fit_nstar(points: list[tuple[float, float]]) -> tuple[float, float]:
    """Least squares of p_B(n) to a logistic. points are (n, p)."""
    if len(points) < 2:
        raise ValueError("need two sizes to fit n*")
    best = (points[0][0], 1.0, 1.0e99)
    ns = [n for n, _p in points]
    lo = min(ns) - 10.0
    hi = max(ns) + 10.0
    for i in range(81):
        n_star = lo + (hi - lo) * i / 80.0
        for j in range(1, 41):
            width = j * (hi - lo) / 40.0
            err = 0.0
            for n, p in points:
                diff = logistic(n, n_star, width) - p
                err += diff * diff
            if err < best[2]:
                best = (n_star, width, err)
    return best[0], best[1]


def p_of_n(shots: list[tuple[float, str]]) -> list[tuple[float, float, int, int]]:
    """(n, p_B, committed, resolved) with unresolved shots left out of p_B."""
    groups: dict[float, list[str]] = {}
    for n_shoot, fate in shots:
        groups.setdefault(n_shoot, []).append(fate)
    rows = []
    for n_shoot in sorted(groups):
        fates = groups[n_shoot]
        resolved = [fate for fate in fates if fate != "unresolved"]
        committed = sum(fate == "committed" for fate in resolved)
        p = committed / len(resolved) if resolved else float("nan")
        rows.append((n_shoot, p, committed, len(resolved)))
    return rows


def bootstrap_nstar(shots: list[tuple[float, str]], n_draw: int = BOOT_DRAWS,
                    seed: int = BOOT_SEED) -> tuple[float, float, float]:
    points = [(n, p) for n, p, _c, resolved in p_of_n(shots) if resolved > 0 and p == p]
    n_star, _width = fit_nstar(points)
    rng = random.Random(seed)
    draws = []
    n = len(shots)
    for _ in range(n_draw):
        sample = [shots[rng.randrange(n)] for _k in range(n)]
        sampled = [(size, p) for size, p, _c, resolved in p_of_n(sample) if resolved > 0 and p == p]
        if len(sampled) < 2:
            continue
        try:
            star, _w = fit_nstar(sampled)
        except ValueError:
            continue
        draws.append(star)
    if not draws:
        return n_star, float("nan"), float("nan")
    draws.sort()
    lo = type7_unit(draws, 0.025)
    hi = type7_unit(draws, 0.975)
    return n_star, lo, hi


def type7_unit(values: list[float], p: float) -> float:
    ordered = sorted(values)
    n = len(ordered)
    if n == 1:
        return ordered[0]
    h = (n - 1) * p
    lo = int(math.floor(h))
    hi = min(int(math.ceil(h)), n - 1)
    weight = h - lo
    return ordered[lo] * (1.0 - weight) + ordered[hi] * weight


def parse_seed_list(text: str) -> list[int]:
    seeds: list[int] = []
    for piece in text.split(","):
        piece = piece.strip()
        if not piece:
            continue
        if "-" in piece:
            lo, hi = piece.split("-", 1)
            seeds.extend(range(int(lo), int(hi) + 1))
        else:
            seeds.append(int(piece))
    if len(seeds) != len(set(seeds)):
        raise ValueError(f"repeated seed in {text}")
    return seeds


def parse_shard(text: str) -> tuple[int, int]:
    if text.count("/") != 1:
        raise ValueError(f"shard must be k/K, not {text}")
    k_text, n_text = text.split("/", 1)
    k, n = int(k_text), int(n_text)
    if n < 1 or k < 0 or k >= n:
        raise ValueError(f"shard {text} is outside 0..{n - 1}")
    return k, n


def take_shard(items: list, shard: str) -> list:
    k, n = parse_shard(shard)
    return [item for index, item in enumerate(items) if index % n == k]


def shard_tag(shard: str) -> str:
    k, n = parse_shard(shard)
    return f"{k}-of-{n}"


def eta_hours(steps_per_s: float) -> dict[str, float]:
    if steps_per_s <= 0.0:
        raise ValueError("steps/s must be positive")
    one = CLIMB_STEPS / steps_per_s / 3600.0
    climbs = N_SEEDS * one
    shots = N_SEEDS * N_SHOOT_SIZES * N_SHOTS * SHOT_STEPS / steps_per_s / 3600.0
    return {"one_climb_h": one, "twenty_climbs_h": climbs, "shooting_h": shots}


def write_rate(path: Path, steps: int, seconds: float, commit: str) -> None:
    rate = steps / seconds
    hours = eta_hours(rate)
    command = (
        "$LMP -in input/nucleation/in.homogeneous "
        "-var n 16 -var nmelt 500 -var nquench 2000 -var ndump 2000 "
        "-var nskip 0 -var Tq 180 -var seed 20261111"
    )
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", newline="") as handle:
        handle.write(f"# commit {commit}\n")
        handle.write(f"# command {command}\n")
        writer = csv.writer(handle, lineterminator="\n")
        writer.writerow([
            "molecules", "quench_steps", "loop_s", "steps_per_s",
            "eta_10ns_h", "eta_20climbs_h", "eta_shooting_h", "commit",
        ])
        writer.writerow([
            4096, steps, f"{seconds:.4f}", f"{rate:.3f}",
            f"{hours['one_climb_h']:.3f}", f"{hours['twenty_climbs_h']:.3f}",
            f"{hours['shooting_h']:.3f}", commit,
        ])


PENDING = "pending: 4096-molecule p_B(n)\n"


def write_pending(table: Path, summary: Path) -> None:
    table.write_text(PENDING)
    summary.write_text(PENDING)


def read_shot_shard(path: Path) -> list[tuple[float, str]]:
    shots = []
    for line in path.read_text().splitlines():
        if not line or line.startswith("#") or line.startswith("seed,"):
            continue
        parts = line.split(",")
        shots.append((float(parts[1]), parts[3]))
    return shots


def merge_shot_dir(directory: Path) -> list[tuple[float, str]]:
    shots: list[tuple[float, str]] = []
    for path in sorted(directory.glob("shot-*-of-*.csv")):
        shots.extend(read_shot_shard(path))
    return shots


def coverage_note(directory: Path | None) -> str:
    if directory is None or not directory.is_dir():
        return ""
    present: dict[int, Path] = {}
    quiet: list[int] = []
    width = None
    for path in sorted(directory.glob("shot-*-of-*.csv")):
        parts = path.stem.split("-")
        if len(parts) != 4 or parts[0] != "shot" or parts[2] != "of":
            continue
        k = int(parts[1])
        width = int(parts[3])
        present[k] = path
        head = "\n".join(path.read_text().splitlines()[:8])
        if "# nucleation no" in head:
            quiet.append(k)
    if width is None:
        return ""
    missing = [i for i in range(width) if i not in present]
    return (
        f"# shards present {sorted(present)} of {width}; "
        f"missing {missing}; nucleation_none {sorted(quiet)}\n"
    )


def write_merged(table: Path, summary: Path, shots: list[tuple[float, str]], commit: str,
                 command: str, directory: Path | None = None) -> None:
    if not shots:
        write_pending(table, summary)
        return
    note = f"# commit {commit}\n# command {command}\n" + coverage_note(directory)
    points = p_of_n(shots)
    usable = [(n, p) for n, p, _c, resolved in points if resolved > 0 and p == p]
    unresolved = sum(fate == "unresolved" for _n, fate in shots)
    resolved = sum(fate != "unresolved" for _n, fate in shots)
    if len(usable) < 2:
        star_s = lo_s = hi_s = ""
    else:
        n_star, lo, hi = bootstrap_nstar(shots)
        star_s, lo_s, hi_s = f"{n_star:.3f}", f"{lo:.3f}", f"{hi:.3f}"
    table.parent.mkdir(parents=True, exist_ok=True)
    with table.open("w", newline="") as handle:
        handle.write(note)
        writer = csv.writer(handle, lineterminator="\n")
        writer.writerow(["shoot_n", "p_B", "committed", "resolved", "commit"])
        for n_shoot, p, committed, n_resolved in points:
            writer.writerow([
                f"{n_shoot:.0f}",
                "" if p != p else f"{p:.6f}",
                committed,
                n_resolved,
                commit,
            ])
    with summary.open("w", newline="") as handle:
        handle.write(note)
        writer = csv.writer(handle, lineterminator="\n")
        writer.writerow(["n_star", "n_star_lo", "n_star_hi", "n_resolved", "n_unresolved", "commit"])
        writer.writerow([star_s, lo_s, hi_s, resolved, unresolved, commit])


def self_test() -> None:
    if commit_size(10) != 50:
        raise SystemExit("a small seed still has to reach 50")
    if commit_size(40) != 80:
        raise SystemExit("twice 40 is 80")
    if not (1.0 < liquid_baseline([0, 0, 0, 1, 2]) < 2.0):
        raise SystemExit(f"baseline {liquid_baseline([0, 0, 0, 1, 2])}")
    base = liquid_baseline([0, 0, 1, 1, 2])
    if shot_fate([20, 30, 60], 20, base) != "committed":
        raise SystemExit("60 should commit a shot from 20")
    if shot_fate([20, 12, base], 20, base) != "uncommitted":
        raise SystemExit("the liquid baseline should uncommit")
    if shot_fate([20, 30, 40], 20, base) != "unresolved":
        raise SystemExit("a shot inside the window is unresolved")
    try:
        check_nmax(160, 64)
    except ValueError:
        pass
    else:
        raise SystemExit("nmax above N was accepted")
    sizes = shooting_sizes([5, 6, 40, 42, 41, 39, 40, 80, 150, 200])
    if len(sizes) < 3 or min(sizes) < 30 or max(sizes) > 50:
        raise SystemExit(f"shooting sizes left the plateau: {sizes}")
    shots = []
    for n, p in ((20, 0.0), (40, 0.25), (60, 0.75), (80, 1.0)):
        committed = int(round(20 * p))
        shots.extend((n, "committed") for _ in range(committed))
        shots.extend((n, "uncommitted") for _ in range(20 - committed))
    star, lo, hi = bootstrap_nstar(shots, n_draw=200, seed=1)
    if not (30.0 < star < 70.0 and lo <= star <= hi):
        raise SystemExit(f"n* {star} interval {lo} {hi}")
    hours = eta_hours(160.0)
    if not (1.5 < hours["one_climb_h"] < 2.0 and hours["shooting_h"] > hours["twenty_climbs_h"]):
        raise SystemExit(f"eta {hours}")
    seeds = parse_seed_list("20261111-20261130")
    if len(seeds) != 20:
        raise SystemExit("the published seed list is not 20 long")
    parts = [take_shard(seeds, f"{k}/20") for k in range(20)]
    if any(len(part) != 1 for part in parts):
        raise SystemExit("a 20-way shard is not one seed")
    flat = [seed for part in parts for seed in part]
    if sorted(flat) != seeds or len(set(flat)) != 20:
        raise SystemExit("shards overlap or drop a seed")
    if set(take_shard(seeds, "0/4")) & set(take_shard(seeds, "1/4")):
        raise SystemExit("shard 0/4 meets shard 1/4")
    import tempfile
    with tempfile.TemporaryDirectory() as tmp:
        folder = Path(tmp)
        (folder / "shot-0-of-2.csv").write_text(
            "seed,shoot_n,shot,fate\n20261111,20,0,uncommitted\n20261111,80,0,committed\n"
        )
        (folder / "shot-1-of-2.csv").write_text(
            "seed,shoot_n,shot,fate\n20261112,20,0,uncommitted\n20261112,80,0,committed\n"
        )
        table = folder / "table.csv"
        summary = folder / "summary.csv"
        write_merged(table, summary, merge_shot_dir(folder), "abc",
                     "$LMP -in input/nucleation/in.homogeneous", folder)
        text = summary.read_text()
        if "$LMP" not in text or "n_star" not in text:
            raise SystemExit("merge did not record n*")
        if "/cvmfs" in text:
            raise SystemExit("merge wrote a library path")
    try:
        parse_shard("5/4")
    except ValueError:
        pass
    else:
        raise SystemExit("shard 5/4 was accepted")
    paper = ROOT / "data"
    for name in ("committor-table.csv", "committor-summary.csv"):
        if (paper / name).read_text() != PENDING:
            raise SystemExit(f"{name} is not the pending header")
    fixture = ROOT / "data" / "committor-merge-fixture"
    banner = "# synthetic test fixture: hand-constructed, not simulation output"
    for name in ("shot-0-of-3.csv", "shot-1-of-3.csv"):
        if not (fixture / name).read_text().startswith(banner):
            raise SystemExit(f"{name} is not labelled as a synthetic fixture")
    if "# nucleation no" not in (fixture / "shot-1-of-3.csv").read_text():
        raise SystemExit("the empty climb is not marked")
    merged = merge_shot_dir(fixture)
    if not any(fate == "unresolved" for _n, fate in merged):
        raise SystemExit("the fixture's unresolved shot was dropped")
    if any(size == 0 for size, _fate in merged):
        raise SystemExit("the no-nucleation shard contributed a shot")
    note = coverage_note(fixture)
    if "missing [2]" not in note or "nucleation_none [1]" not in note:
        raise SystemExit(f"partial shards were not reported: {note}")
    with tempfile.TemporaryDirectory() as merged_tmp:
        out_table = Path(merged_tmp) / "table.csv"
        out_summary = Path(merged_tmp) / "summary.csv"
        write_merged(out_table, out_summary, merged, "fixture", "$LMP", fixture)
        text = out_summary.read_text()
        if "missing [2]" not in text or "nucleation_none [1]" not in text:
            raise SystemExit(f"partial merge dropped the shard note: {text}")
        row = text.strip().splitlines()[-1].split(",")
        if row[0] == "" or row[4] == "0":
            raise SystemExit(f"partial merge did not fit n*: {text}")
    only_quiet = ROOT / "data" / "committor-merge-fixture"
    quiet_dir_shots = read_shot_shard(only_quiet / "shot-1-of-3.csv")
    if quiet_dir_shots:
        raise SystemExit("a no-nucleation shard yielded shots")
    print("self-test ok")


def head_commit() -> str:
    done = subprocess.run(
        ["git", "rev-parse", "HEAD"], cwd=ROOT, text=True, capture_output=True, check=True
    )
    return done.stdout.strip()


def write_climb_manifest(path: Path, seeds: str, shard: str, commit: str) -> None:
    chosen = take_shard(parse_seed_list(seeds), shard)
    path.parent.mkdir(parents=True, exist_ok=True)
    command = (
        "$LMP -in input/nucleation/in.homogeneous -var n 16 "
        f"--seeds {seeds} --shard {shard}"
    )
    lines = [f"# commit {commit}", f"# command {command}", f"# shard {shard}", "seed"]
    lines.extend(str(seed) for seed in chosen)
    path.write_text("\n".join(lines) + "\n")


def main() -> None:
    if "--self-test" in sys.argv:
        self_test()
        return
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--out", type=Path)
    parser.add_argument("--quench-steps", type=int)
    parser.add_argument("--quench-seconds", type=float)
    parser.add_argument("--merge", type=Path)
    parser.add_argument("--table", type=Path)
    parser.add_argument("--summary", type=Path)
    parser.add_argument("--print-shard", action="store_true")
    parser.add_argument("--write-climb", type=Path)
    parser.add_argument("--write-shot", type=Path)
    parser.add_argument("--seeds", default="20261111-20261130")
    parser.add_argument("--shard", default="0/1")
    parser.add_argument("--commit", default="")
    args = parser.parse_args()
    commit = args.commit or head_commit()
    if args.print_shard:
        chosen = take_shard(parse_seed_list(args.seeds), args.shard)
        print(",".join(str(seed) for seed in chosen))
        return
    if args.write_climb is not None:
        write_climb_manifest(args.write_climb, args.seeds, args.shard, commit)
        return
    if args.write_shot is not None:
        lmp = os.environ.get("LMP", "")
        lib = os.environ.get("SEAMS_LIB", "")
        if not lmp:
            raise SystemExit("LMP is unset")
        if not lib:
            raise SystemExit("SEAMS_LIB is unset")
        raise SystemExit(
            "shot shard waits on finished 4096 climbs; "
            f"no fates are invented for shard {args.shard}"
        )
    if args.merge is not None:
        if args.table is None or args.summary is None:
            raise SystemExit("--merge needs --table and --summary")
        shots = merge_shot_dir(args.merge) if args.merge.is_dir() else []
        command = "$LMP -in input/nucleation/in.homogeneous -var n 16 --seeds 20261111-20261130"
        write_merged(args.table, args.summary, shots, commit, command, args.merge)
        return
    if args.out is None or args.quench_steps is None or args.quench_seconds is None:
        raise SystemExit("need --out, --quench-steps and --quench-seconds")
    write_rate(args.out, args.quench_steps, args.quench_seconds, commit)
    print(args.out.read_text())


if __name__ == "__main__":
    main()
