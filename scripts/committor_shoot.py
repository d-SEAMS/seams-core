#!/usr/bin/env python3
"""Committor driver for the 4096-molecule mW lane at 180 K.

The order parameter is n_chill_core, the largest cluster of published
CHILL+ cubic and hexagonal ice from nucleation-pathway. Every climb
and shot table also carries n_chill_ext (cubic, hexagonal and
interfacial) and n_q6solid (the ten Wolde q6 solid cluster, parameters
in src/include/internal/nucleation_pathway.hpp). Only n_chill_core
enters the rules below.

All cuts are fixed here, before any shot outcome is read.

Nucleation. A climb has nucleated when n_chill_core passes 300 and
then stays at or above 300 for 1 ns. The onset is the first frame of
that passage. A passage that the end of the run cuts short is not an
onset: the climb is 'not nucleated, extend' and is continued. Any
other climb that has not nucleated is 'no-nucleation'. The one rule,
nucleation_onset(), decides the climb outcome and which climbs are
shot; only nucleated climbs are shot.

Uncommit cut. The baseline is always the 95th percentile of
n_chill_core in a separate 2 ns liquid at 180 K, over its frames from
1 ns, with a velocity seed derived from the climb seed. The climb
frames before onset are never used: they already hold the growth.

Shooting points. For n in 50, 75, 100, 150, 200 and 300, the last
frame at or before the onset where n_chill_core crosses up to n.
Sizes below the liquid fluctuation ceiling are kept: p_B there should
be near 0.

Shots. 20 per point, fresh Maxwell-Boltzmann velocities at 180 K from
a seed derived from (climb seed, n, shot index). n_chill_core is read
every 10 ps. A shot commits at max(2n, 300), uncommits at or below
the baseline, and is unresolved after 1 ns.

p_B(n) is the committed fraction among resolved shots that started at
size n. n* is the size where a logistic fit of p_B crosses one half.
The 64-molecule deck is a fixture for the driver, not this table.

One VM runs one shard. Shard k of K is the seeds whose index modulo K
equals k, with k starting at 0. LMP is a launcher command line and
LMP_NP adds mpirun (scripts/lmp_launch.py). Build the census in
release; a debug build writes the same frames about twice as
slowly. Snakemake 9 reads every word after --config as a value, so
the target comes first:

    meson setup bbdir --buildtype=release && ninja -C bbdir nucleation-pathway
    export LMP=lmp LMP_NP=$(nproc --all)
    export SEAMS_LIB=$PWD/bbdir/src
    snakemake -j1 data/committor-shards/climb-0-of-20.csv \\
        --config seeds=20261111-20261130 shard=0/20
    snakemake -j1 data/committor-shards/shot-0-of-20.csv \\
        --config seeds=20261111-20261130 shard=0/20

A climb that has not nucleated is continued 10 ns at a time, up to
40 ns, by the extension rule, which rewrites the same climb table.
Segments run in order, 1 to 3:

    snakemake -j1 data/committor-shards/extend-0-of-20.3.log \\
        --config seeds=20261111-20261130 shard=0/20

science_committor_merge reads every shot shard and writes p_B(n) and
n*. Until those shards exist, the paper tables stay the single line
'pending: 4096-molecule p_B(n)'.
"""

from __future__ import annotations

import argparse
import csv
import hashlib
import math
import os
import random
import shutil
import signal
import subprocess
import sys
import time
from dataclasses import dataclass, replace
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "scripts"))
import lmp_launch  # noqa: E402

# Cuts chosen before any shot outcome is known.
COMMIT_MULTIPLE = 2.0
COMMIT_FLOOR = 300
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

ORDER = "n_chill_core"
SERIES = ("n_chill_core", "n_chill_ext", "n_q6solid")
NUCLEUS_CUT = 300
HOLD_PS = 1000.0
SHOOT_SIZES = (50, 75, 100, 150, 200, 300)
EQUIL_PS = 1000.0
TEMPERATURE = 180.0
N_CELL = 16
MELT_STEPS = 10_000
DUMP_STEPS = 1000
LIQUID_STEPS = 200_000
EXTEND_STEPS = 1_000_000
CLIMB_LIMIT_STEPS = 4_000_000
FRAME_PS = DUMP_STEPS * TIMESTEP_FS / 1000.0


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


def nucleation_onset(series: list[float], cut: float = NUCLEUS_CUT,
                     hold_ps: float = HOLD_PS,
                     times_ps: list[float] | None = None) -> int | None:
    """First frame of a passage above cut after which the series stays at
    or above cut for hold_ps. A passage cut short by the end of the run
    is not an onset. times_ps defaults to frames FRAME_PS apart."""
    times = times_ps if times_ps is not None else [FRAME_PS * i for i in range(len(series))]
    i = 0
    n = len(series)
    while i < n:
        if series[i] > cut and (i == 0 or series[i - 1] <= cut):
            j = i
            while j < n and series[j] >= cut and times[j] - times[i] < hold_ps:
                j += 1
            if j < n and series[j] >= cut:
                return i
            i = j
        i += 1
    return None


def open_passage(series: list[float], cut: float = NUCLEUS_CUT) -> bool:
    """The run ends on a passage above cut that has not yet fallen below."""
    j = len(series)
    while j > 0 and series[j - 1] >= cut:
        j -= 1
    return any(value > cut for value in series[j:])


EXTEND = "not nucleated, extend"


def climb_outcome(series: list[float], times_ps: list[float], cut: float = NUCLEUS_CUT,
                  hold_ps: float = HOLD_PS, limit_ps: float = math.inf) -> str:
    """nucleated; 'not nucleated, extend' when the run ends inside a passage
    not yet held for hold_ps and is short of the limit; else no-nucleation."""
    if nucleation_onset(series, cut, hold_ps, times_ps) is not None:
        return "nucleated"
    if open_passage(series, cut) and times_ps[-1] < limit_ps:
        return EXTEND
    return "no-nucleation"


def shooting_points(series: list[float], onset: int,
                    sizes: tuple[int, ...] = SHOOT_SIZES) -> list[tuple[int, int]]:
    """(n, frame index) at the last up-crossing of n at or before onset."""
    points = []
    for n in sizes:
        for j in range(min(onset, len(series) - 1), 0, -1):
            if series[j] >= n and series[j - 1] < n:
                points.append((n, j))
                break
    return points


def derived_seed(*parts: object) -> int:
    """LAMMPS velocity seed, 1..9e8, fixed by its parts."""
    text = ":".join(str(part) for part in parts)
    digest = hashlib.sha256(text.encode()).digest()
    return int.from_bytes(digest[:4], "big") % 900_000_000 + 1


def shooting_sizes(series: list[float], n_sizes: int = N_SHOOT_SIZES) -> list[int]:
    """Sizes around the plateau that precedes the final climb.

    The plateau is the run of samples immediately before the final
    strict ascent whose values stay within 20 percent of their median.
    The returned sizes are evenly spaced across that run. They are not
    fit to a committor. A climb that has not nucleated has no plateau
    before a climb and returns nothing.
    """
    if len(series) < 3 or n_sizes < 1:
        return []
    if nucleation_onset(series) is None:
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


CLIMB_FIELDS = [
    "seed", "outcome", "onset_ps",
    "n_chill_core_max", "t_max_ps", "n_chill_ext_max", "n_q6solid_max",
    "n_chill_core_end", "n_chill_ext_end", "n_q6solid_end",
    "t_end_ps", "frames", "segments", "extension",
]


def read_climb_outcomes(path: Path) -> dict[int, dict[str, str]]:
    """Climb rows by seed. Reads the outcome columns, or the older table
    that lists seeds and marks the outcome with '# nucleation no|yes'."""
    lines = path.read_text().splitlines()
    comments = [line for line in lines if line.startswith("#")]
    body = [line for line in lines if line and not line.startswith("#")]
    if not body:
        return {}
    header = body[0].split(",")
    if "outcome" in header:
        return {int(row["seed"]): row for row in csv.DictReader(body)}
    marked = ""
    if any(line.strip() == "# nucleation no" for line in comments):
        marked = "no-nucleation"
    elif any(line.strip() == "# nucleation yes" for line in comments):
        marked = "nucleated"
    rows = {}
    for line in body[1:]:
        seed = int(line.split(",")[0])
        rows[seed] = {"seed": str(seed), "outcome": marked}
    return rows


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
        climb = directory / f"climb-{k}-of-{width}.csv"
        outcomes = read_climb_outcomes(climb) if climb.is_file() else {}
        if "# nucleation no" in head or (
            outcomes and all(row["outcome"] != "nucleated" for row in outcomes.values())
        ):
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


def self_test_rules() -> None:
    # Frames are 10 ps apart, so 1 ns of hold is 100 frames past the onset.
    held = [1, 100] + [301] * 101
    if nucleation_onset(held) != 2:
        raise SystemExit("301 held for 1 ns is nucleation at frame 2")
    if nucleation_onset([1, 100] + [301] * 100) is not None:
        raise SystemExit("990 ps above 300 is short of the hold")
    if nucleation_onset([1, 100, 301, 450, 610, 900]) is not None:
        raise SystemExit("a run that ends above 300 before the hold is not nucleated")
    if nucleation_onset([1, 320, 290] + [310] * 101) != 3:
        raise SystemExit("a dip below 300 restarts the hold")
    if nucleation_onset([1, 320] + [300] * 100) != 1:
        raise SystemExit("300 itself is not below 300")
    if nucleation_onset([1, 200] + [300] * 200) is not None:
        raise SystemExit("300 is not past 300")
    if nucleation_onset([1, 301, 301, 301], times_ps=[0.0, 10.0, 500.0, 1010.0]) != 1:
        raise SystemExit("the hold is read from the frame times")
    # Shard 0, seed 20261111 at 30 ns: 274 then 311 on the last frame.
    end_run = [150, 525, 274, 311]
    times = [29840.0, 29870.0, 29990.0, 30000.0]
    if climb_outcome(end_run, times, limit_ps=40000.0) != EXTEND:
        raise SystemExit("an onset on the last frame was not sent to extend")
    if climb_outcome(end_run, times, limit_ps=30000.0) != "no-nucleation":
        raise SystemExit("a climb at the limit was sent to extend")
    if climb_outcome([1, 90, 160, 120], times, limit_ps=40000.0) != "no-nucleation":
        raise SystemExit("a climb below 300 at the end is not an open passage")
    if climb_outcome(held, [10.0 * i for i in range(len(held))]) != "nucleated":
        raise SystemExit("a held passage is nucleated")
    climb = [1, 30, 60, 40, 55, 80, 70, 49, 51, 90, 120, 160, 210, 320] + [480] * 100
    onset = nucleation_onset(climb)
    points = dict(shooting_points(climb, onset))
    want = {50: 8, 75: 9, 100: 10, 150: 11, 200: 12, 300: 13}
    if onset != 13 or points != want:
        raise SystemExit(f"last crossings {points} onset {onset}")
    if shooting_points([1, 30, 301, 650], 2) != [(50, 2), (75, 2), (100, 2), (150, 2),
                                                  (200, 2), (300, 2)]:
        raise SystemExit("a jump crosses every size at once")
    if derived_seed(20261125, 50, 0) != derived_seed(20261125, 50, 0):
        raise SystemExit("a derived seed moved")
    seeds = {derived_seed(20261125, n, k) for n in SHOOT_SIZES for k in range(N_SHOTS)}
    if len(seeds) != len(SHOOT_SIZES) * N_SHOTS or min(seeds) < 1:
        raise SystemExit("derived seeds collide")
    import tempfile
    with tempfile.TemporaryDirectory() as tmp:
        folder = Path(tmp)
        old = folder / "climb-1-of-3.csv"
        old.write_text("# commit x\n# nucleation no\nseed\n20261112\n")
        if read_climb_outcomes(old)[20261112]["outcome"] != "no-nucleation":
            raise SystemExit("the old comment form was not read")
        new = folder / "climb-0-of-3.csv"
        new.write_text("# commit x\n" + ",".join(CLIMB_FIELDS) + "\n"
                       "20261111,no-nucleation,,160,9000,1069,349,"
                       "120,1000,300,10000,1001,1,\n")
        row = read_climb_outcomes(new)[20261111]
        if row["outcome"] != "no-nucleation" or row["n_chill_core_max"] != "160":
            raise SystemExit(f"climb columns {row}")
        for k in (0, 1):
            (folder / f"shot-{k}-of-3.csv").write_text("seed,shoot_n,shot,fate\n")
        note = coverage_note(folder)
        (folder / "climb-1-of-3.csv").write_text(
            "# commit x\n" + ",".join(CLIMB_FIELDS) + "\n"
            '20261112,"not nucleated, extend",,311,30000,1779,845,'
            "311,1779,845,30000,3001,3,restart\n")
        if read_climb_outcomes(folder / "climb-1-of-3.csv")[20261112]["outcome"] != EXTEND:
            raise SystemExit("the extend outcome did not survive the csv")
        note = coverage_note(folder)
        if "nucleation_none [0, 1]" not in note:
            raise SystemExit(f"climb outcomes did not reach the merge: {note}")
        frames = folder / "frames.tsv"
        frames.write_text("frame\ttimestep\tn_chill_core\tn_chill_ext\n"
                          "1\t10000\t3\t9\n2\t11000\t5\t12\n")
        rows = read_frames(frames)
        if has_series(rows):
            raise SystemExit("a census without n_q6solid was taken as complete")
        rows.append({"frame": "9", "timestep": "12000", "n_chill_core": "7",
                     "n_chill_ext": "15", "n_q6solid": "4"})
        write_frames(frames, rows)
        again = read_frames(frames)
        if [r["frame"] for r in again] != ["1", "2", "3"] or again[0]["n_q6solid"] != "":
            raise SystemExit(f"frames were not merged: {again}")
        if frame_times(again) != [0.0, 10.0, 20.0]:
            raise SystemExit("frame times are not 10 ps apart")
        row = climb_row(20261111, again, folder)
        if row["outcome"] != "no-nucleation" or row["n_q6solid_max"] != "" or row["t_end_ps"] != "20":
            raise SystemExit(f"climb row {row}")
        dump = folder / "a.lammpstrj"
        frame_text = ("ITEM: TIMESTEP\n{t}\nITEM: NUMBER OF ATOMS\n2\n"
                      "ITEM: BOX BOUNDS pp pp pp\n-1 9\n0 10\n0 10\n"
                      "ITEM: ATOMS id type x y z\n2 1 {x} 1 1\n1 1 0.5 0.5 0.5\n")
        dump.write_text(frame_text.format(t=0, x=1) + frame_text.format(t=1000, x=2))
        step, bounds, atoms = read_dump_frame(dump)
        if step != 1000 or bounds[0] != (-1.0, 9.0) or sorted(atoms)[1][2] != 2.0:
            raise SystemExit("the last dump frame was not read")
        if read_dump_frame(dump, 0)[2][0][2] != 1.0 or read_dump_frame(dump, 5) is not None:
            raise SystemExit("a dump frame by step was not read")
        write_data(folder / "a.data", bounds, atoms)
        if "-1.0000000000 9.0000000000 xlo xhi" not in (folder / "a.data").read_text():
            raise SystemExit("the data file lost the box origin")


def self_test() -> None:
    if commit_size(50) != 300:
        raise SystemExit("a small seed still has to reach 300")
    if commit_size(200) != 400:
        raise SystemExit("twice 200 is 400")
    if not (1.0 < liquid_baseline([0, 0, 0, 1, 2]) < 2.0):
        raise SystemExit(f"baseline {liquid_baseline([0, 0, 0, 1, 2])}")
    base = liquid_baseline([0, 0, 1, 1, 2])
    if shot_fate([20, 30, 60, 300], 20, base) != "committed":
        raise SystemExit("300 should commit a shot from 20")
    if shot_fate([150, 299], 150, base) != "unresolved":
        raise SystemExit("299 is short of twice 150")
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
    sizes = shooting_sizes([5, 6, 40, 42, 41, 39, 40, 80, 150, 350] + [700] * 101)
    if len(sizes) < 3 or min(sizes) < 30 or max(sizes) > 50:
        raise SystemExit(f"shooting sizes left the plateau: {sizes}")
    if shooting_sizes([5, 6, 40, 42, 41, 39, 40, 80, 150, 200]) != []:
        raise SystemExit("a climb that stops at 200 has shooting sizes")
    # Shard 17: fluctuations from 30 to 92 over 10 ns, no climb.
    wander = [1, 5, 20] + [30 + (37 * i) % 62 for i in range(200)] + [92, 61]
    if shooting_sizes(wander) != [] or nucleation_onset(wander) is not None:
        raise SystemExit("a liquid that wanders to 92 was taken as nucleated")
    self_test_rules()
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


@dataclass(frozen=True)
class Plan:
    """Run lengths and the nucleation cuts. The defaults are the paper
    lane. The commit and uncommit rules are not part of a plan."""

    n_cell: int = N_CELL
    temperature: float = TEMPERATURE
    cut: float = NUCLEUS_CUT
    hold_ps: float = HOLD_PS
    sizes: tuple[int, ...] = SHOOT_SIZES
    shots: int = N_SHOTS
    shot_steps: int = SHOT_STEPS
    dump_steps: int = DUMP_STEPS
    melt_steps: int = MELT_STEPS
    liquid_steps: int = LIQUID_STEPS
    extend_steps: int = EXTEND_STEPS
    limit_steps: int = CLIMB_LIMIT_STEPS
    equil_ps: float = EQUIL_PS

    @property
    def molecules(self) -> int:
        return self.n_cell ** 3

    @property
    def limit_ps(self) -> float:
        return self.limit_steps * TIMESTEP_FS / 1000.0


PAPER = Plan()
SHOT_FIELDS = ["seed", "shoot_n", "shot", "fate", "t_end_ps", "nmax_end",
               "n_chill_ext_end", "n_q6solid_end"]
BASELINE_FIELDS = ["seed", "baseline", "baseline_source", "window_lo_ps", "window_hi_ps",
                   "frames_used", "liquid_seed", "onset_ps"]
POINT_FIELDS = ["seed", "shoot_n", "t_ps", "timestep", "n_chill_core", "n_chill_ext",
                "n_q6solid"]
SEGMENT_FIELDS = ["segment", "source", "velocity_seed", "step0", "steps"]


def rel(path: Path) -> str:
    try:
        return str(Path(path).resolve().relative_to(ROOT))
    except ValueError:
        return Path(path).name


def fmt(value: float | int | str | None) -> str:
    if value is None or value == "":
        return ""
    if isinstance(value, str):
        return value
    if value != value:
        return ""
    if float(value).is_integer():
        return str(int(value))
    return f"{value:.3f}"


def seed_dir(frames_root: Path, seed: int, plan: Plan = PAPER) -> Path:
    return frames_root / f"mw-{plan.molecules}" / f"seed-{seed}"


def read_frames(path: Path) -> list[dict[str, str]]:
    with path.open() as handle:
        return list(csv.DictReader(handle, delimiter="\t"))


def write_frames(path: Path, rows: list[dict[str, str]]) -> None:
    fields: list[str] = []
    for row in rows:
        for key in row:
            if key not in fields:
                fields.append(key)
    first = int(rows[0]["frame"]) if rows and rows[0].get("frame") else 1
    tmp = path.with_name(path.name + ".tmp")
    with tmp.open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields, delimiter="\t",
                                lineterminator="\n", restval="")
        writer.writeheader()
        for index, row in enumerate(rows):
            writer.writerow({**row, "frame": str(first + index)})
    tmp.replace(path)


def frame_times(rows: list[dict[str, str]]) -> list[float]:
    step0 = int(rows[0]["timestep"])
    return [(int(row["timestep"]) - step0) * TIMESTEP_FS / 1000.0 for row in rows]


def column(rows: list[dict[str, str]], name: str) -> list[float]:
    out = []
    for row in rows:
        text = row.get(name) or ""
        out.append(float(text) if text != "" else float("nan"))
    return out


def has_series(rows: list[dict[str, str]]) -> bool:
    return bool(rows) and all(all((row.get(name) or "") != "" for row in rows)
                              for name in SERIES)


def check_release(build_dir: Path) -> None:
    """The census is the release build. Debug gives the same series slower."""
    import json
    info = build_dir / "meson-info" / "intro-buildoptions.json"
    if not info.is_file():
        raise SystemExit(f"{build_dir} is not a meson build; "
                         f"run meson setup {build_dir} --buildtype=release")
    options = {item["name"]: item["value"] for item in json.loads(info.read_text())}
    if options.get("buildtype") != "release":
        raise SystemExit(f"{build_dir} is a {options.get('buildtype')} build; "
                         f"run meson configure {build_dir} --buildtype=release")


def pathway_path(explicit: str = "") -> str:
    chosen = explicit or os.environ.get("NUCLEATION_PATHWAY", "") or str(
        ROOT / "bbdir" / "nucleation-pathway")
    if not Path(chosen).is_file():
        raise SystemExit(f"nucleation-pathway not found: {chosen}")
    return chosen


def census_env(threads: int) -> dict[str, str]:
    env = os.environ.copy()
    lib = os.environ.get("SEAMS_LIB", "")
    if lib:
        env["LD_LIBRARY_PATH"] = lib
    env["OMP_NUM_THREADS"] = str(threads)
    return env


def census(pathway: str, dump: Path, out: Path, threads: int = 1) -> list[dict[str, str]]:
    """n_chill_core and the other columns of nucleation-pathway, one row
    per frame of dump."""
    done = subprocess.run(
        [pathway, str(dump), "--out", str(out), "--type", "1", "--threads", str(threads)],
        cwd=ROOT, env=census_env(threads), capture_output=True, text=True,
    )
    if done.returncode != 0:
        raise SystemExit(f"census failed on {dump}\n{done.stderr[-1500:]}")
    return read_frames(out)


def lmp_env() -> dict[str, str]:
    env = os.environ.copy()
    env.pop("LD_LIBRARY_PATH", None)
    env["OMP_NUM_THREADS"] = "1"
    return env


def lmp_args(lmp: list[str], deck: str, variables: dict[str, object], log: str) -> list[str]:
    args = [*lmp, "-in", deck]
    for key, value in variables.items():
        args += ["-var", key, str(value)]
    return args + ["-log", log, "-screen", "none"]


def run_lammps(lmp: list[str], work: Path, deck: str, variables: dict[str, object],
               log: str = "log.lammps") -> None:
    with (work / (Path(log).stem + ".out")).open("w") as handle:
        proc = subprocess.run(lmp_args(lmp, deck, variables, log), cwd=work, env=lmp_env(),
                              stdout=handle, stderr=subprocess.STDOUT)
    if proc.returncode != 0:
        path = work / log
        tail = path.read_text(errors="replace")[-1500:] if path.is_file() else ""
        raise SystemExit(f"LAMMPS failed in {work}\n{tail}")


def stage(work: Path, deck: str) -> None:
    work.mkdir(parents=True, exist_ok=True)
    shutil.copy(ROOT / "input" / "mw215" / "mW.sw", work / "mW.sw")
    shutil.copy(ROOT / "input" / "nucleation" / deck, work / deck)


Frame = tuple[int, list[tuple[float, float]], list[tuple[int, int, float, float, float]]]


def read_dump_frame(path: Path, timestep: int | None = None) -> Frame | None:
    """The frame at timestep, or the last frame. None when timestep is absent."""
    offset = None
    with path.open() as handle:
        while True:
            where = handle.tell()
            line = handle.readline()
            if not line:
                break
            if not line.startswith("ITEM: TIMESTEP"):
                continue
            step = int(handle.readline())
            handle.readline()
            count = int(handle.readline())
            for _ in range(5 + count):
                handle.readline()
            if timestep is None or step == timestep:
                offset = where
                if timestep is not None:
                    break
        if offset is None:
            return None
        handle.seek(offset)
        handle.readline()
        step = int(handle.readline())
        handle.readline()
        count = int(handle.readline())
        handle.readline()
        bounds = []
        for _ in range(3):
            lo, hi = handle.readline().split()[:2]
            bounds.append((float(lo), float(hi)))
        header = handle.readline().split()[2:]
        idx = {name: i for i, name in enumerate(header)}
        atoms = []
        for _ in range(count):
            parts = handle.readline().split()
            atoms.append((int(parts[idx["id"]]), int(parts[idx["type"]]),
                          float(parts[idx["x"]]), float(parts[idx["y"]]),
                          float(parts[idx["z"]])))
    return step, bounds, atoms


def write_data(path: Path, bounds: list[tuple[float, float]],
               atoms: list[tuple[int, int, float, float, float]]) -> None:
    lines = ["mW frame for a committor run", "", f"{len(atoms)} atoms", "1 atom types", ""]
    for (lo, hi), axis in zip(bounds, "xyz"):
        lines.append(f"{lo:.10f} {hi:.10f} {axis}lo {axis}hi")
    lines += ["", "Masses", "", "1 18.015", "", "Atoms", ""]
    for atom in sorted(atoms):
        lines.append(f"{atom[0]} {atom[1]} {atom[2]:.8f} {atom[3]:.8f} {atom[4]:.8f}")
    path.write_text("\n".join(lines) + "\n")


def climb_dumps(directory: Path) -> list[Path]:
    extensions = sorted(directory.glob("extend-*.lammpstrj"),
                        key=lambda p: int(p.stem.split("-")[1]))
    return [directory / "quench.lammpstrj", *extensions]


def read_segments(directory: Path) -> list[dict[str, str]]:
    path = directory / "segments.tsv"
    if not path.is_file():
        return []
    return read_frames(path)


def append_segment(directory: Path, row: dict[str, object]) -> None:
    path = directory / "segments.tsv"
    fresh = not path.is_file()
    with path.open("a", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=SEGMENT_FIELDS, delimiter="\t",
                                lineterminator="\n")
        if fresh:
            writer.writeheader()
        writer.writerow({key: row.get(key, "") for key in SEGMENT_FIELDS})


def recensus(directory: Path, pathway: str, threads: int) -> list[dict[str, str]]:
    """Census every dump of a climb again and rebuild frames.tsv."""
    merged: list[dict[str, str]] = []
    for dump in climb_dumps(directory):
        if not dump.is_file():
            raise SystemExit(f"missing {dump}; the census cannot be rebuilt")
        rows = census(pathway, dump, directory / f"{dump.stem}.census.tsv", threads)
        last = int(merged[-1]["timestep"]) if merged else None
        merged.extend(row for row in rows if last is None or int(row["timestep"]) > last)
    write_frames(directory / "frames.tsv", merged)
    return merged


def climb_frames(directory: Path, pathway: str | None, threads: int) -> list[dict[str, str]]:
    path = directory / "frames.tsv"
    if not path.is_file():
        raise SystemExit(f"no census for this climb: {path}")
    rows = read_frames(path)
    if not has_series(rows) and pathway:
        rows = recensus(directory, pathway, threads)
    return rows


def climb_row(seed: int, rows: list[dict[str, str]], directory: Path,
              plan: Plan = PAPER) -> dict[str, str]:
    times = frame_times(rows)
    core = column(rows, ORDER)
    for size in core:
        check_nmax(size, plan.molecules)
    onset = nucleation_onset(core, plan.cut, plan.hold_ps, times)
    top = max(range(len(core)), key=lambda i: core[i])
    ext = column(rows, "n_chill_ext")
    q6 = column(rows, "n_q6solid")
    segments = read_segments(directory)
    extension = ";".join(
        "restart" if seg["source"] == "restart" else f"dump-mb:{seg['velocity_seed']}"
        for seg in segments
    )
    return {
        "seed": str(seed),
        "outcome": climb_outcome(core, times, plan.cut, plan.hold_ps, plan.limit_ps),
        "onset_ps": fmt(times[onset]) if onset is not None else "",
        "n_chill_core_max": fmt(core[top]),
        "t_max_ps": fmt(times[top]),
        "n_chill_ext_max": fmt(max(ext)) if all(v == v for v in ext) else "",
        "n_q6solid_max": fmt(max(q6)) if all(v == v for v in q6) else "",
        "n_chill_core_end": fmt(core[-1]),
        "n_chill_ext_end": fmt(ext[-1]),
        "n_q6solid_end": fmt(q6[-1]),
        "t_end_ps": fmt(times[-1]),
        "frames": str(len(rows)),
        "segments": str(1 + len(segments)),
        "extension": extension,
    }


def rule_note(plan: Plan) -> list[str]:
    return [
        f"# order {ORDER}",
        f"# rule nucleated when {ORDER} passes {fmt(plan.cut)} and stays at or above "
        f"{fmt(plan.cut)} for {fmt(plan.hold_ps)} ps; a passage cut short by the end of the "
        f"run is '{EXTEND}' until {fmt(plan.limit_ps)} ps",
    ]


def write_climb(path: Path, seeds: str, shard: str, frames_root: Path, commit: str,
                pathway: str | None = None, threads: int = 1, plan: Plan = PAPER) -> None:
    chosen = take_shard(parse_seed_list(seeds), shard)
    rows = []
    for seed in chosen:
        directory = seed_dir(frames_root, seed, plan)
        rows.append(climb_row(seed, climb_frames(directory, pathway, threads), directory, plan))
    command = (
        f"python3 scripts/committor_shoot.py --write-climb {rel(path)} "
        f"--frames {rel(frames_root)} --seeds {seeds} --shard {shard}"
    )
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", newline="") as handle:
        handle.write(f"# commit {commit}\n# command {command}\n# shard {shard}\n")
        handle.write("\n".join(rule_note(plan)) + "\n")
        writer = csv.DictWriter(handle, fieldnames=CLIMB_FIELDS, lineterminator="\n")
        writer.writeheader()
        writer.writerows(rows)


def extend_seed(seed: int, directory: Path, lmp: list[str], pathway: str, threads: int,
                plan: Plan = PAPER) -> str:
    """Continue one climb that has not nucleated by plan.extend_steps."""
    rows = climb_frames(directory, pathway, threads)
    core = column(rows, ORDER)
    if nucleation_onset(core, plan.cut, plan.hold_ps, frame_times(rows)) is not None:
        return "nucleated"
    step0 = int(rows[0]["timestep"])
    step_end = int(rows[-1]["timestep"])
    done = step_end - step0
    if done >= plan.limit_steps:
        return "at the limit"
    steps = min(plan.extend_steps, plan.limit_steps - done)
    segment = len(read_segments(directory)) + 1
    previous = "quench" if segment == 1 else f"extend-{segment - 1}"
    stage(directory, "in.extend")
    out = f"extend-{segment}"
    variables: dict[str, object] = {
        "nrun": steps, "ndump": plan.dump_steps, "Tq": plan.temperature, "out": out,
    }
    restart = directory / f"{previous}.restart"
    if restart.is_file():
        variables.update(source="restart", restart=restart.name)
        source, velocity_seed = "restart", ""
    else:
        frame = read_dump_frame(directory / f"{previous}.lammpstrj")
        if frame is None or frame[0] != step_end:
            raise SystemExit(f"{previous}.lammpstrj does not end at step {step_end}")
        write_data(directory / "extend.data", frame[1], frame[2])
        velocity_seed = derived_seed(seed, "extend", segment)
        variables.update(source="data", datafile="extend.data", step0=step_end,
                         seed=velocity_seed)
        source = "dump-mb"
    run_lammps(lmp, directory, "in.extend", variables, log=f"{out}.log")
    fresh = census(pathway, directory / f"{out}.lammpstrj",
                   directory / f"{out}.census.tsv", threads)
    rows.extend(row for row in fresh if int(row["timestep"]) > step_end)
    write_frames(directory / "frames.tsv", rows)
    append_segment(directory, {"segment": segment, "source": source,
                               "velocity_seed": velocity_seed, "step0": step_end,
                               "steps": steps})
    return f"{source} {steps} steps"


def wait_frame(path: Path, molecules: int, proc: subprocess.Popen) -> None:
    need = 9 + molecules
    while True:
        if path.is_file():
            data = path.read_bytes()
            if data.endswith(b"\n") and data.count(b"\n") >= need:
                return
        if proc.poll() is not None:
            if path.is_file():
                data = path.read_bytes()
                if data.endswith(b"\n") and data.count(b"\n") >= need:
                    return
            log = path.parent / "log.lammps"
            tail = log.read_text(errors="replace")[-1500:] if log.is_file() else ""
            raise SystemExit(f"LAMMPS stopped before {path.name} in {path.parent}\n{tail}")
        time.sleep(0.2)


def stop(proc: subprocess.Popen) -> None:
    if proc.poll() is not None:
        return
    for sig, wait in ((signal.SIGTERM, 15), (signal.SIGKILL, 15)):
        try:
            os.killpg(proc.pid, sig)
        except ProcessLookupError:
            return
        try:
            proc.wait(timeout=wait)
            return
        except subprocess.TimeoutExpired:
            continue


def run_shot(lmp: list[str], work: Path, data: Path, velocity_seed: int, n_shoot: int,
             baseline: float, pathway: str, threads: int = 1,
             plan: Plan = PAPER) -> dict[str, str]:
    """One shot. LAMMPS runs the whole window; it is stopped as soon as a
    10 ps frame resolves the shot."""
    if work.exists():
        shutil.rmtree(work)
    stage(work, "in.shot")
    shutil.copy(data, work / "shot.data")
    variables = {"datafile": "shot.data", "nrun": plan.shot_steps, "ndump": plan.dump_steps,
                 "Tq": plan.temperature, "seed": velocity_seed}
    trace: list[tuple[float, dict[str, float]]] = []
    fate = "unresolved"
    with (work / "lammps.out").open("w") as handle:
        proc = subprocess.Popen(lmp_args(lmp, "in.shot", variables, "log.lammps"), cwd=work,
                                env=lmp_env(), stdout=handle, stderr=subprocess.STDOUT,
                                start_new_session=True)
        try:
            for k in range(1, plan.shot_steps // plan.dump_steps + 1):
                step = k * plan.dump_steps
                frame = work / f"frame.{step}.lammpstrj"
                wait_frame(frame, plan.molecules, proc)
                row = census(pathway, frame, work / "check.tsv", threads)[-1]
                values = {name: float(row[name]) for name in SERIES}
                check_nmax(values[ORDER], plan.molecules)
                trace.append((step * TIMESTEP_FS / 1000.0, values))
                frame.unlink()
                fate = shot_fate([values[ORDER]], n_shoot, baseline)
                if fate != "unresolved":
                    break
        finally:
            stop(proc)
    for leftover in work.glob("frame.*.lammpstrj"):
        leftover.unlink()
    with (work / "trace.tsv").open("w") as handle:
        handle.write("t_ps\t" + "\t".join(SERIES) + "\n")
        for t_ps, values in trace:
            handle.write(fmt(t_ps) + "\t" + "\t".join(fmt(values[n]) for n in SERIES) + "\n")
    t_end, last = trace[-1]
    return {
        "fate": fate,
        "t_end_ps": fmt(t_end),
        "nmax_end": fmt(last[ORDER]),
        "n_chill_ext_end": fmt(last["n_chill_ext"]),
        "n_q6solid_end": fmt(last["n_q6solid"]),
    }


def liquid_run(seed: int, directory: Path, lmp: list[str], pathway: str, threads: int,
               plan: Plan = PAPER) -> tuple[float, list[float], int]:
    """Baseline from a separate 2 ns liquid at the climb temperature."""
    work = directory / "liquid-2ns"
    liquid_seed = derived_seed(seed, "liquid")
    if not (work / "frames.tsv").is_file():
        stage(work, "in.homogeneous")
        run_lammps(lmp, work, "in.homogeneous", {
            "n": plan.n_cell, "nmelt": plan.melt_steps, "nquench": plan.liquid_steps,
            "ndump": plan.dump_steps, "nskip": 0, "Tq": plan.temperature, "seed": liquid_seed,
        })
        census(pathway, work / "quench.lammpstrj", work / "frames.tsv", threads)
    rows = read_frames(work / "frames.tsv")
    times = frame_times(rows)
    core = column(rows, ORDER)
    onset = nucleation_onset(core, plan.cut, plan.hold_ps, times)
    end = len(core) if onset is None else onset
    frames = [i for i in range(end) if times[i] >= plan.equil_ps]
    if not frames:
        raise SystemExit(f"the liquid for seed {seed} has no frames after {plan.equil_ps} ps")
    return liquid_baseline([core[i] for i in frames]), [times[i] for i in frames], liquid_seed


def seed_baseline(seed: int, onset_ps: float, directory: Path, lmp: list[str],
                  pathway: str, threads: int, plan: Plan = PAPER) -> dict[str, str]:
    """The uncommit cut always comes from the separate liquid, never from
    the climb, whose frames before onset already hold growth."""
    value, used, liquid_seed = liquid_run(seed, directory, lmp, pathway, threads, plan)
    return {
        "seed": str(seed), "baseline": fmt(value), "baseline_source": "liquid-2ns",
        "window_lo_ps": fmt(used[0]), "window_hi_ps": fmt(used[-1]),
        "frames_used": str(len(used)), "liquid_seed": str(liquid_seed),
        "onset_ps": fmt(onset_ps),
    }


def shooting_frame(directory: Path, timestep: int, out: Path) -> None:
    if out.is_file():
        return
    for dump in climb_dumps(directory):
        if not dump.is_file():
            continue
        frame = read_dump_frame(dump, timestep)
        if frame is not None:
            out.parent.mkdir(parents=True, exist_ok=True)
            write_data(out, frame[1], frame[2])
            return
    raise SystemExit(f"step {timestep} is in no dump of {directory}")


def sibling(path: Path, stem: str) -> Path:
    tag = path.name.split("-", 1)[1] if "-" in path.name else path.name
    return path.with_name(f"{stem}-{tag}")


def write_rows(path: Path, note: list[str], fields: list[str], rows: list[dict]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", newline="") as handle:
        handle.write("\n".join(note) + "\n")
        writer = csv.DictWriter(handle, fieldnames=fields, lineterminator="\n")
        writer.writeheader()
        writer.writerows(rows)


def write_shots(out: Path, seeds: str, shard: str, frames_root: Path, commit: str,
                lmp: list[str], pathway: str, threads: int = 1, plan: Plan = PAPER,
                baseline_out: Path | None = None, points_out: Path | None = None) -> None:
    """Shoot every nucleated climb of the shard. Rows go to a .partial file
    as each shot ends, so a stopped shard resumes where it was."""
    chosen = take_shard(parse_seed_list(seeds), shard)
    baseline_out = baseline_out or sibling(out, "baseline")
    points_out = points_out or sibling(out, "points")
    note = [
        f"# commit {commit}",
        f"# command python3 scripts/committor_shoot.py --write-shot {rel(out)} "
        f"--frames {rel(frames_root)} --seeds {seeds} --shard {shard}",
        f"# shard {shard}",
        *rule_note(plan),
        f"# commit at max({fmt(COMMIT_MULTIPLE)}n, {COMMIT_FLOOR}); uncommit at the baseline; "
        f"check every {fmt(plan.dump_steps * TIMESTEP_FS / 1000.0)} ps; "
        f"unresolved after {fmt(plan.shot_steps * TIMESTEP_FS / 1000.0)} ps",
        f"# launcher {' '.join(Path(w).name if '/' in w else w for w in lmp)}",
    ]
    partial = out.with_name(out.name + ".partial")
    done: set[tuple[int, int, int]] = set()
    if partial.is_file():
        for line in partial.read_text().splitlines():
            if line and not line.startswith("#") and not line.startswith("seed,"):
                parts = line.split(",")
                done.add((int(parts[0]), int(parts[1]), int(parts[2])))
    else:
        write_rows(partial, note, SHOT_FIELDS, [])
    baselines = []
    points = []
    for seed in chosen:
        directory = seed_dir(frames_root, seed, plan)
        rows = climb_frames(directory, pathway, threads)
        core = column(rows, ORDER)
        times = frame_times(rows)
        onset = nucleation_onset(core, plan.cut, plan.hold_ps, times)
        if onset is None:
            outcome = climb_outcome(core, times, plan.cut, plan.hold_ps, plan.limit_ps)
            print(f"seed {seed} {outcome}, max {fmt(max(core))}; no shots", flush=True)
            continue
        base = seed_baseline(seed, times[onset], directory, lmp, pathway, threads, plan)
        baselines.append(base)
        baseline = float(base["baseline"])
        for n_shoot, index in shooting_points(core, onset, plan.sizes):
            row = rows[index]
            points.append({
                "seed": str(seed), "shoot_n": str(n_shoot), "t_ps": fmt(times[index]),
                "timestep": row["timestep"], "n_chill_core": row[ORDER],
                "n_chill_ext": row.get("n_chill_ext", ""), "n_q6solid": row.get("n_q6solid", ""),
            })
            data = directory / "shots" / f"n{n_shoot}-step{row['timestep']}" / "shot.data"
            shooting_frame(directory, int(row["timestep"]), data)
            for shot in range(plan.shots):
                if (seed, n_shoot, shot) in done:
                    continue
                result = run_shot(lmp, data.parent / f"shot-{shot}", data,
                                  derived_seed(seed, n_shoot, shot), n_shoot, baseline,
                                  pathway, threads, plan)
                with partial.open("a", newline="") as handle:
                    csv.DictWriter(handle, fieldnames=SHOT_FIELDS, lineterminator="\n").writerow(
                        {"seed": seed, "shoot_n": n_shoot, "shot": shot, **result})
                print(f"seed {seed} n {n_shoot} shot {shot} {result['fate']} "
                      f"t {result['t_end_ps']} ps nmax {result['nmax_end']}", flush=True)
    write_rows(baseline_out, note[:3], BASELINE_FIELDS, baselines)
    write_rows(points_out, note[:4], POINT_FIELDS, points)
    partial.replace(out)


def main() -> None:
    if "--self-test" in sys.argv:
        self_test()
        return
    if len(sys.argv) == 3 and sys.argv[1] == "--check-release":
        check_release(Path(sys.argv[2]))
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
    parser.add_argument("--extend-climb", type=Path,
                        help="continue every climb of the shard that has not nucleated, "
                             "rewrite the climb table and write this log")
    parser.add_argument("--climb", type=Path, help="climb table of the shard")
    parser.add_argument("--frames", type=Path, help="climb frames, the regenerate.py --dest")
    parser.add_argument("--pathway", default="", help="nucleation-pathway; NUCLEATION_PATHWAY")
    parser.add_argument("--threads", type=int, default=1, help="census threads")
    parser.add_argument("--seeds", default="20261111-20261130")
    parser.add_argument("--shard", default="0/1")
    parser.add_argument("--commit", default="")
    args = parser.parse_args()
    commit = args.commit or head_commit()
    if args.print_shard:
        chosen = take_shard(parse_seed_list(args.seeds), args.shard)
        print(",".join(str(seed) for seed in chosen))
        return
    if args.write_climb is not None or args.write_shot is not None or args.extend_climb:
        if args.frames is None:
            raise SystemExit("--frames is the climb directory, the regenerate.py --dest")
    if args.write_climb is not None:
        write_climb(args.write_climb, args.seeds, args.shard, args.frames, commit,
                    pathway_path(args.pathway), args.threads)
        print(args.write_climb.read_text())
        return
    if args.extend_climb is not None:
        if args.climb is None:
            raise SystemExit("--extend-climb needs --climb")
        lmp = lmp_launch.launcher()
        pathway = pathway_path(args.pathway)
        lines = []
        for seed in take_shard(parse_seed_list(args.seeds), args.shard):
            result = extend_seed(seed, seed_dir(args.frames, seed), lmp, pathway, args.threads)
            lines.append(f"seed {seed} {result}")
            print(lines[-1], flush=True)
        write_climb(args.climb, args.seeds, args.shard, args.frames, commit, pathway,
                    args.threads)
        args.extend_climb.parent.mkdir(parents=True, exist_ok=True)
        args.extend_climb.write_text(f"# commit {commit}\n" + "\n".join(lines) + "\n")
        return
    if args.write_shot is not None:
        lmp = lmp_launch.launcher()
        write_shots(args.write_shot, args.seeds, args.shard, args.frames, commit, lmp,
                    pathway_path(args.pathway), args.threads)
        return
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
