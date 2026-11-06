#!/usr/bin/env python3
"""Head-to-head ring-statistics timings for d-SEAMS vs OVITO.

Protocol mirrors neighbour timing:
- fresh process per sample
- interleaved order with per-repetition rotation
- n >= 6 from confidence target
- median as mean of the middle pair
"""

from __future__ import annotations

import argparse
import json
import math
import os
import platform
import re
import subprocess
import sys
import time
import warnings
from dataclasses import dataclass
from pathlib import Path
from typing import Callable

import numpy as np

MIN_REPEATS = 6
WARMUP_CALLS = 2
HEXAGONAL_RING_SIZE = 6


@dataclass(frozen=True)
class Frame:
    positions: np.ndarray
    H: np.ndarray
    origin: np.ndarray
    is_triclinic: bool
    source: str


@dataclass(frozen=True)
class DatasetSpec:
    name: str
    path: str
    frame: int
    atom_type: int
    rep_xyz: tuple[int, int, int]


def quantile_linear(sorted_vals: list[float], q: float) -> float:
    if not sorted_vals:
        return 0.0
    if len(sorted_vals) == 1:
        return sorted_vals[0]
    pos = q * (len(sorted_vals) - 1)
    lo = int(pos)
    frac = pos - lo
    if lo + 1 >= len(sorted_vals):
        return sorted_vals[-1]
    return sorted_vals[lo] * (1.0 - frac) + sorted_vals[lo + 1] * frac


def type7(values: list[float], probability: float) -> float:
    ordered = sorted(float(value) for value in values)
    count = len(ordered)
    if count < MIN_REPEATS:
        raise ValueError(f"need n >= {MIN_REPEATS}, got {count}")
    position = (count - 1) * probability + 1.0
    index = int(math.floor(position))
    if index >= count:
        return ordered[-1]
    if index < 1:
        return ordered[0]
    weight = position - index
    return ordered[index - 1] + weight * (ordered[index] - ordered[index - 1])


def summarize_samples(samples_ms: list[float]) -> dict[str, float]:
    if len(samples_ms) < 2:
        raise ValueError("need at least two samples")
    s = sorted(samples_ms)
    mid = len(s) // 2
    if len(s) % 2 == 0:
        median = 0.5 * (s[mid - 1] + s[mid])
    else:
        median = s[mid]
    q1 = quantile_linear(s, 0.25)
    q3 = quantile_linear(s, 0.75)
    return {
        "median_ms": float(median),
        "q1_ms": float(q1),
        "q3_ms": float(q3),
        "iqr_ms": float(q3 - q1),
    }


def derive_repeats(confidence_target: float) -> int:
    base = int(math.ceil(1.0 - math.log2(1.0 - confidence_target)))
    n = max(MIN_REPEATS, base)
    if n % 2 != 0:
        n += 1
    return n


def derive_warmup(repeats: int) -> int:
    if repeats < WARMUP_CALLS:
        raise ValueError(f"repeats={repeats} must be >= warmup calls {WARMUP_CALLS}")
    # One cold-start sink plus one stabilization call.
    return WARMUP_CALLS


def self_check_summary_rules() -> None:
    even = summarize_samples([1.0, 2.0, 3.0, 4.0, 5.0, 6.0])
    if abs(even["median_ms"] - 3.5) > 1e-12:
        raise AssertionError(f"even-n median regression: {even['median_ms']}")
    odd = summarize_samples([1.0, 2.0, 3.0, 4.0, 9.0])
    if abs(odd["median_ms"] - 3.0) > 1e-12:
        raise AssertionError(f"odd-n median regression: {odd['median_ms']}")


def rotate(values: list[str], shift: int) -> list[str]:
    if not values:
        return []
    k = shift % len(values)
    return values[k:] + values[:k]


def display_cmd(cmd: list[str]) -> str:
    cwd = os.getcwd()
    shown = []
    for part in cmd:
        if part.startswith(cwd + os.sep):
            shown.append(os.path.relpath(part, cwd))
        else:
            shown.append(part)
    return " ".join(shown)


def pinned(cmd: list[str]) -> list[str]:
    return ["taskset", "-c", "0", *cmd]


def git_sha() -> str:
    return subprocess.check_output(["git", "rev-parse", "HEAD"], text=True).strip()


def git_dirty() -> str:
    unstaged = subprocess.call(
        ["git", "diff", "--quiet", "--ignore-submodules", "HEAD"]
    )
    staged = subprocess.call(
        ["git", "diff", "--cached", "--quiet", "--ignore-submodules", "HEAD"]
    )
    return "dirty" if unstaged != 0 or staged != 0 else "clean"


def published_cpu_model(name: str) -> str:
    """Drop the register mark so the record reads 'Intel Xeon Processor'."""
    return " ".join(name.replace("(R)", "").split())


def cpu_model_name() -> str:
    try:
        with Path("/proc/cpuinfo").open("r", encoding="utf-8") as fh:
            for line in fh:
                if line.startswith("model name"):
                    parts = line.split(":", 1)
                    if len(parts) == 2:
                        return published_cpu_model(parts[1].strip())
    except OSError:
        pass
    return "unknown"


def runtime_machine_header() -> dict[str, object]:
    return {
        "uname_m": platform.machine(),
        "nproc": int(os.cpu_count() or 1),
        "cpu_model": cpu_model_name(),
    }


def parse_box(bounds_header: str, bounds_lines: list[str]) -> tuple[np.ndarray, np.ndarray, bool]:
    tric = "xy" in bounds_header.lower()
    if not tric:
        xlo, xhi = map(float, bounds_lines[0].split()[:2])
        ylo, yhi = map(float, bounds_lines[1].split()[:2])
        zlo, zhi = map(float, bounds_lines[2].split()[:2])
        H = np.array(
            [[xhi - xlo, 0.0, 0.0], [0.0, yhi - ylo, 0.0], [0.0, 0.0, zhi - zlo]],
            dtype=np.float64,
        )
        return H, np.array([xlo, ylo, zlo], dtype=np.float64), False
    xlo_b, xhi_b, xy = map(float, bounds_lines[0].split()[:3])
    ylo_b, yhi_b, xz = map(float, bounds_lines[1].split()[:3])
    zlo_b, zhi_b, yz = map(float, bounds_lines[2].split()[:3])
    xlo = xlo_b - min(0.0, xy, xz, xy + xz)
    xhi = xhi_b - max(0.0, xy, xz, xy + xz)
    ylo = ylo_b - min(0.0, yz)
    yhi = yhi_b - max(0.0, yz)
    H = np.array(
        [[xhi - xlo, xy, xz], [0.0, yhi - ylo, yz], [0.0, 0.0, zhi_b - zlo_b]],
        dtype=np.float64,
    )
    return H, np.array([xlo, ylo, zlo_b], dtype=np.float64), True


def read_lammps_frame(path: str, target_frame: int) -> Frame:
    p = Path(path)
    if not p.is_file():
        raise FileNotFoundError(path)
    with p.open("r", encoding="utf-8") as fh:
        while True:
            line = fh.readline()
            if not line:
                raise ValueError(f"frame {target_frame} not found in {path}")
            if not line.startswith("ITEM: TIMESTEP"):
                continue
            timestep = int(fh.readline().strip())
            if not fh.readline().startswith("ITEM: NUMBER OF ATOMS"):
                raise ValueError(f"malformed dump: {path}")
            n_atoms = int(fh.readline().strip())
            header = fh.readline().strip()
            bounds = [fh.readline().strip(), fh.readline().strip(), fh.readline().strip()]
            atoms_header = fh.readline().strip()
            cols = atoms_header.split()[2:]
            colmap = {name: idx for idx, name in enumerate(cols)}
            if timestep != target_frame and target_frame != 1:
                for _ in range(n_atoms):
                    fh.readline()
                continue
            pts = np.empty((n_atoms, 3), dtype=np.float64)
            ids = np.empty(n_atoms, dtype=np.int64)
            for i in range(n_atoms):
                row = fh.readline().split()
                ids[i] = int(row[colmap["id"]])
                pts[i, 0] = float(row[colmap["x"]])
                pts[i, 1] = float(row[colmap["y"]])
                pts[i, 2] = float(row[colmap["z"]])
            order = np.argsort(ids, kind="mergesort")
            pts = pts[order]
            H, origin, tric = parse_box(header, bounds)
            return Frame(
                positions=pts,
                H=H,
                origin=origin,
                is_triclinic=tric,
                source=f"{path}@{target_frame}",
            )


def replicate_frame(frame: Frame, nx: int, ny: int, nz: int) -> Frame:
    if nx == 1 and ny == 1 and nz == 1:
        return frame
    a = frame.H[:, 0]
    b = frame.H[:, 1]
    c = frame.H[:, 2]
    offsets = []
    for ix in range(nx):
        for iy in range(ny):
            for iz in range(nz):
                offsets.append(ix * a + iy * b + iz * c)
    pts = np.concatenate([frame.positions + off for off in offsets], axis=0)
    H = frame.H.copy()
    H[:, 0] *= nx
    H[:, 1] *= ny
    H[:, 2] *= nz
    return Frame(
        positions=pts,
        H=H,
        origin=frame.origin.copy(),
        is_triclinic=frame.is_triclinic,
        source=f"{frame.source} rep={nx}x{ny}x{nz}",
    )


def write_lammps_dump(path: str, frame: Frame, atom_type: int) -> None:
    p = Path(path)
    p.parent.mkdir(parents=True, exist_ok=True)
    H = frame.H
    lx = float(H[0, 0])
    xy = float(H[0, 1])
    xz = float(H[0, 2])
    ly = float(H[1, 1])
    yz = float(H[1, 2])
    lz = float(H[2, 2])
    xlo = float(frame.origin[0])
    ylo = float(frame.origin[1])
    zlo = float(frame.origin[2])
    xmin = min(0.0, xy, xz, xy + xz)
    xmax = max(0.0, xy, xz, xy + xz)
    ymin = min(0.0, yz)
    ymax = max(0.0, yz)
    xlo_b = xlo + xmin
    xhi_b = xlo + lx + xmax
    ylo_b = ylo + ymin
    yhi_b = ylo + ly + ymax
    zlo_b = zlo
    zhi_b = zlo + lz
    with p.open("w", encoding="utf-8") as fh:
        fh.write("ITEM: TIMESTEP\n1\n")
        fh.write("ITEM: NUMBER OF ATOMS\n")
        fh.write(f"{frame.positions.shape[0]}\n")
        fh.write("ITEM: BOX BOUNDS xy xz yz pp pp pp\n")
        fh.write(f"{xlo_b:.12f} {xhi_b:.12f} {xy:.12f}\n")
        fh.write(f"{ylo_b:.12f} {yhi_b:.12f} {xz:.12f}\n")
        fh.write(f"{zlo_b:.12f} {zhi_b:.12f} {yz:.12f}\n")
        fh.write("ITEM: ATOMS id type x y z\n")
        for i, pos in enumerate(frame.positions, start=1):
            fh.write(f"{i:d} {atom_type:d} {pos[0]:.12f} {pos[1]:.12f} {pos[2]:.12f}\n")


def save_frame_npz(frame: Frame, path: str) -> None:
    Path(path).parent.mkdir(parents=True, exist_ok=True)
    np.savez(
        path,
        positions=frame.positions,
        H=frame.H,
        origin=frame.origin,
        triclinic=np.array([1 if frame.is_triclinic else 0], dtype=np.int8),
        source=np.array([frame.source], dtype=object),
    )


def load_frame_npz(path: str) -> Frame:
    with np.load(path, allow_pickle=True) as data:
        return Frame(
            positions=data["positions"],
            H=data["H"],
            origin=data["origin"],
            is_triclinic=bool(int(data["triclinic"][0])),
            source=str(data["source"][0]),
        )


def parse_bench_stages(out: str) -> tuple[int, int, float, float, float]:
    rings_pat = re.compile(r"^# .* rings ([0-9]+) six ([0-9]+) ddc ([0-9]+).*$", re.MULTILINE)
    neigh_pat = re.compile(r"^neigh_ms ([0-9.eE+-]+)$", re.MULTILINE)
    ring_pat = re.compile(r"^rings_ms ([0-9.eE+-]+)$", re.MULTILINE)
    io_pat = re.compile(r"^io_ms ([0-9.eE+-]+)$", re.MULTILINE)
    m_rings = rings_pat.search(out)
    m_neigh = neigh_pat.search(out)
    m_ring = ring_pat.search(out)
    m_io = io_pat.search(out)
    if not (m_rings and m_neigh and m_ring and m_io):
        raise RuntimeError(f"cannot parse bench_stages output\n{out}")
    return (
        int(m_rings.group(1)),
        int(m_rings.group(2)),
        float(m_neigh.group(1)),
        float(m_ring.group(1)),
        float(m_io.group(1)),
    )


def run_dseams_once(
    bench_stages_bin: str,
    traj: str,
    frame: int,
    atom_type: int,
    threads: int,
    git_sha_value: str,
    dirty_value: str,
) -> dict[str, object]:
    cmd = pinned([bench_stages_bin, traj, str(frame), str(atom_type), "1"])
    env = os.environ.copy()
    # The bench is linked against the toolchain loader. A system library
    # directory on the loader path pulls the wrong libc.
    env.pop("LD_LIBRARY_PATH", None)
    env["OMP_NUM_THREADS"] = str(threads)
    env["OMP_DYNAMIC"] = "FALSE"
    env["SEAMS_GIT_SHA"] = git_sha_value
    env["SEAMS_GIT_DIRTY"] = dirty_value
    out = subprocess.check_output(cmd, text=True, env=env)
    rings_total, rings_six, neigh_ms, ring_ms, io_ms = parse_bench_stages(out)
    cmd_text = (
        f"OMP_NUM_THREADS={threads} OMP_DYNAMIC=FALSE "
        f"SEAMS_GIT_SHA={git_sha_value} "
        f"SEAMS_GIT_DIRTY={dirty_value} "
        + display_cmd(cmd)
    )
    return {
        "command": cmd_text,
        "rings_total": rings_total,
        "rings_size6": rings_six,
        "ring_ms": ring_ms,
        "neigh_ms": neigh_ms,
        "io_ms": io_ms,
        "total_ms": neigh_ms + ring_ms,
    }


def drop_file_cache(path: str) -> None:
    fd = os.open(path, os.O_RDONLY)
    try:
        os.posix_fadvise(fd, 0, 0, os.POSIX_FADV_DONTNEED)
    finally:
        os.close(fd)


def read_dseams_file(bench_stages_bin: str, traj: str, frame: int, atom_type: int) -> float:
    cmd = pinned([bench_stages_bin, traj, str(frame), str(atom_type), "1"])
    env = os.environ.copy()
    env.pop("LD_LIBRARY_PATH", None)
    env["OMP_NUM_THREADS"] = "1"
    env["OMP_DYNAMIC"] = "FALSE"
    env["SEAMS_TIME_READ"] = "1"
    out = subprocess.check_output(cmd, text=True, env=env)
    _total, _six, _neigh, _ring, io_ms = parse_bench_stages(out)
    return float(io_ms)


def read_ovito_file(script_path: str, traj: str) -> float:
    cmd = pinned([
        sys.executable,
        script_path,
        "--worker-ovito",
        "--frame-dump",
        traj,
        "--ovito-stage",
        "setup",
        "--frame-npz",
        traj,
    ])
    out = json.loads(subprocess.check_output(cmd, text=True))
    return float(out["ms"])


def setup_reads(bench_stages_bin: str, script_path: str, traj: str, frame: int, atom_type: int) -> dict[str, object]:
    # Twelve fresh processes on each side. The first read follows an advice
    # to drop the file from the page cache. The other eleven are the
    # amortized reads, summarized by the type 7 median.
    drop_file_cache(traj)
    dseams_cold = read_dseams_file(bench_stages_bin, traj, frame, atom_type)
    drop_file_cache(traj)
    ovito_cold = read_ovito_file(script_path, traj)
    dseams_samples = []
    ovito_samples = []
    for _ in range(11):
        dseams_samples.append(read_dseams_file(bench_stages_bin, traj, frame, atom_type))
        ovito_samples.append(read_ovito_file(script_path, traj))
    return {
        "dseams_cold": dseams_cold,
        "ovito_cold": ovito_cold,
        "dseams_samples": dseams_samples,
        "ovito_samples": ovito_samples,
        "dseams": type7(dseams_samples, 0.5),
        "ovito": type7(ovito_samples, 0.5),
    }


def run_ovito_worker(
    script_path: str,
    frame_npz: str,
    cutoff: float,
    max_ring_size: int,
    ovito_threads: int,
    stage: str,
) -> dict[str, object]:
    cmd = [
        sys.executable,
        script_path,
        "--worker-ovito",
        "--frame-npz",
        frame_npz,
        "--cutoff",
        f"{cutoff:.8f}",
        "--max-ring-size",
        str(max_ring_size),
        "--ovito-threads",
        str(ovito_threads),
        "--ovito-stage",
        stage,
    ]
    pinned_cmd = pinned(cmd)
    out = json.loads(subprocess.check_output(pinned_cmd, text=True))
    out["command"] = display_cmd(pinned_cmd)
    return out


def run_interleaved(
    tools: list[str],
    warmup_runs: int,
    repeats: int,
    run_tool: Callable[[str], dict[str, object]],
    label: str,
) -> dict[str, object]:
    records: dict[str, dict[str, object]] = {}
    active: list[str] = []

    for tool in tools:
        res = run_tool(tool)
        if "unsupported" in res:
            continue
        active.append(tool)
        records[tool] = {
            "first_call_ms": float(res["first_ms"]),
            "warmup_runs": warmup_runs,
            "repeats": repeats,
            "samples_ms": [],
            "count": int(res["count"]),
            "command": res.get("command"),
        }

    if not active:
        return {"unsupported": "no active tools"}

    for warm in range(warmup_runs):
        for tool in rotate(active, warm):
            res = run_tool(tool)
            if "unsupported" in res:
                raise RuntimeError(f"{label} warm-up became unsupported for {tool}")
            if int(res["count"]) != int(records[tool]["count"]):
                raise RuntimeError(f"{label} warm-up count drift for {tool}")

    for rep in range(repeats):
        for tool in rotate(active, rep):
            res = run_tool(tool)
            if "unsupported" in res:
                raise RuntimeError(f"{label} sample became unsupported for {tool}")
            if int(res["count"]) != int(records[tool]["count"]):
                raise RuntimeError(f"{label} sample count drift for {tool}")
            records[tool]["samples_ms"].append(float(res["sample_ms"]))

    for tool in active:
        records[tool].update(summarize_samples(records[tool]["samples_ms"]))

    return records


def _ring_histogram(out) -> tuple[int, int]:
    hist = out.tables["ring-size-histogram"]
    sizes = hist["Ring Size"][...]
    counts = hist["Count"][...]
    total = int(np.sum(counts))
    six = (
        int(counts[np.where(sizes == HEXAGONAL_RING_SIZE)[0][0]])
        if np.any(sizes == HEXAGONAL_RING_SIZE)
        else 0
    )
    return total, six


def worker_ovito(args: argparse.Namespace) -> None:
    warnings.filterwarnings("ignore", message=".*OVITO.*PyPI")
    os.environ["OVITO_THREAD_COUNT"] = str(args.ovito_threads)
    from ovito.data import DataCollection
    from ovito.modifiers import CreateBondsModifier, FindRingsModifier
    from ovito.pipeline import Pipeline, StaticSource

    stage = args.ovito_stage
    if stage == "setup" and args.frame_dump:
        from ovito.io import import_file

        t0 = time.perf_counter()
        data = import_file(args.frame_dump).compute()
        t1 = time.perf_counter()
        count = int(data.particles.count) if data.particles is not None else 0
        print(json.dumps({"ms": (t1 - t0) * 1000.0, "count": count, "rings_size6": 0}))
        return
    frame = load_frame_npz(args.frame_npz)
    t_setup0 = time.perf_counter()
    data = DataCollection()
    data.create_cell(frame.H, pbc=(True, True, True))
    parts = data.create_particles(count=frame.positions.shape[0])
    parts.create_property("Position", data=frame.positions - frame.origin)
    t_setup1 = time.perf_counter()
    if stage == "setup":
        print(json.dumps({
            "ms": (t_setup1 - t_setup0) * 1000.0,
            "count": int(frame.positions.shape[0]),
            "rings_size6": 0,
        }))
        return

    if stage == "rings":
        # Bonds are built before the timer. The timed call is FindRings only.
        bond_pipe = Pipeline(source=StaticSource(data=data))
        bond_pipe.modifiers.append(CreateBondsModifier(cutoff=args.cutoff))
        bonded = bond_pipe.compute()
        ring_pipe = Pipeline(source=StaticSource(data=bonded))
        ring_pipe.modifiers.append(FindRingsModifier(maximum_ring_size=args.max_ring_size))
        t0 = time.perf_counter()
        out = ring_pipe.compute()
        t1 = time.perf_counter()
    elif stage == "bonds":
        pipe = Pipeline(source=StaticSource(data=data))
        pipe.modifiers.append(CreateBondsModifier(cutoff=args.cutoff))
        pipe.modifiers.append(FindRingsModifier(maximum_ring_size=args.max_ring_size))
        t0 = time.perf_counter()
        out = pipe.compute()
        t1 = time.perf_counter()
    else:
        raise SystemExit(f"unknown ovito stage {stage}")
    total, six = _ring_histogram(out)
    print(json.dumps({"ms": (t1 - t0) * 1000.0, "count": total, "rings_size6": six}))


def main() -> None:
    self_check_summary_rules()

    ap = argparse.ArgumentParser()
    ap.add_argument("--worker-ovito", action="store_true")
    ap.add_argument("--frame-npz")
    ap.add_argument("--bench-stages-bin")
    ap.add_argument("--output-json")
    ap.add_argument("--cutoff", type=float, default=3.5)
    ap.add_argument("--max-ring-size", type=int, default=HEXAGONAL_RING_SIZE)
    ap.add_argument("--confidence", type=float, default=0.95)
    ap.add_argument("--repeats", type=int, default=0)
    ap.add_argument("--warmup", type=int, default=-1)
    ap.add_argument("--max-molecules", type=int, default=0)
    ap.add_argument("--seed", type=int, default=7)
    ap.add_argument("--dseams-threads", type=int, default=1)
    ap.add_argument("--ovito-threads", type=int, default=1)
    ap.add_argument("--frame-dump")
    ap.add_argument("--ovito-stage", choices=("setup", "rings", "bonds"), default="bonds")
    args = ap.parse_args()

    if args.worker_ovito:
        worker_ovito(args)
        return

    if not args.bench_stages_bin or not args.output_json:
        raise SystemExit("--bench-stages-bin and --output-json are required")

    sha = git_sha()
    dirty = git_dirty()
    repeats = args.repeats if args.repeats > 0 else derive_repeats(args.confidence)
    warmup = args.warmup if args.warmup >= 0 else derive_warmup(repeats)
    machine = runtime_machine_header()

    print(
        f"# header git_sha {sha} dirty {dirty} seed {args.seed} repeats {repeats} warmup {warmup}",
        flush=True,
    )
    print(
        "# machine "
        f"uname_m {machine['uname_m']} nproc {machine['nproc']} cpu_model {machine['cpu_model']}",
        flush=True,
    )

    generated_specs: list[DatasetSpec] = []
    for name, path, rep in [
        ("orth-1e5", "input/traj/mW_cubic.lammpstrj", (3, 3, 3)),
        ("tri-1e5", "input/traj/genice_sH.lammpstrj", (10, 10, 10)),
    ]:
        base = read_lammps_frame(path, 1)
        rep_frame = replicate_frame(base, *rep)
        if args.max_molecules > 0 and int(rep_frame.positions.shape[0]) > args.max_molecules:
            continue
        out_path = f"input/traj/{name}-ring-scale.lammpstrj"
        write_lammps_dump(out_path, rep_frame, atom_type=1)
        generated_specs.append(DatasetSpec(name, out_path, 1, 1, (1, 1, 1)))

    datasets = [
        DatasetSpec("orth-4k", "input/traj/mW_cubic.lammpstrj", 1, 1, (1, 1, 1)),
        DatasetSpec("tri-136", "input/traj/genice_sH.lammpstrj", 1, 1, (1, 1, 1)),
        *generated_specs,
    ]

    out: dict[str, object] = {
        "git_sha": sha,
        "git_dirty": dirty,
        "seed": args.seed,
        "cutoff": args.cutoff,
        "max_ring_size": args.max_ring_size,
        "machine": machine,
        "protocol": {
            "confidence_target": args.confidence,
            "repeat_lower_bound": "n >= ceil(1 - log2(1 - confidence_target))",
            "median_rule": "type 7 median and quartiles",
            "repeats": repeats,
            "warmup_runs": warmup,
            "warmup_rule": (
                "fixed 2 calls (one cold-start sink plus one stabilization call)"
            ),
            "fresh_process_per_sample": True,
            "interleaving": "per repetition, rotate tool order",
            "first_call_reported_separately": True,
        },
        "thread_policy": {
            "dseams_threads_reported": [1, args.dseams_threads],
            "ovito_threads": args.ovito_threads,
        },
        "datasets": {},
    }

    frame_store = Path("input/.bench_frames")
    frame_store.mkdir(parents=True, exist_ok=True)

    for spec in datasets:
        frame = read_lammps_frame(spec.path, spec.frame)
        n_atoms = int(frame.positions.shape[0])
        if args.max_molecules > 0 and n_atoms > args.max_molecules:
            print(f"[ring] {spec.name} n={n_atoms} skipped", flush=True)
            continue
        frame_npz = str(frame_store / f"ring-{spec.name}.npz")
        save_frame_npz(frame, frame_npz)

        print(f"[ring] {spec.name} n={n_atoms}", flush=True)

        setup_reads_row = setup_reads(
            args.bench_stages_bin,
            str(Path(__file__)),
            spec.path,
            spec.frame,
            spec.atom_type,
        )
        setup_dseams = run_dseams_once(
            args.bench_stages_bin,
            spec.path,
            spec.frame,
            spec.atom_type,
            1,
            sha,
            dirty,
        )
        print(
            f"[ring] {spec.name} setup cold dseams_ms {setup_reads_row['dseams_cold']:.6f} "
            f"ovito_ms {setup_reads_row['ovito_cold']:.6f} "
            f"amortized dseams_ms {setup_reads_row['dseams']:.6f} "
            f"ovito_ms {setup_reads_row['ovito']:.6f}",
            flush=True,
        )

        def run_arm(tool: str, work: str) -> dict[str, object]:
            if tool.startswith("dseams"):
                threads = 1 if tool.endswith("omp1") else args.dseams_threads
                res = run_dseams_once(
                    args.bench_stages_bin,
                    spec.path,
                    spec.frame,
                    spec.atom_type,
                    threads,
                    sha,
                    dirty,
                )
                sample = float(res["ring_ms"] if work == "rings" else res["total_ms"])
                return {
                    "first_ms": sample,
                    "sample_ms": sample,
                    "count": int(res["rings_total"]),
                    "command": res["command"],
                }
            ov = run_ovito_worker(
                script_path=str(Path(__file__)),
                frame_npz=frame_npz,
                cutoff=args.cutoff,
                max_ring_size=args.max_ring_size,
                ovito_threads=args.ovito_threads,
                stage="rings" if work == "rings" else "bonds",
            )
            return {
                "first_ms": float(ov["ms"]),
                "sample_ms": float(ov["ms"]),
                "count": int(ov["count"]),
                "command": ov["command"],
            }

        def timed_pair(work: str) -> dict[str, object]:
            stage = run_interleaved(
                ["dseams_omp1", "ovito"],
                warmup,
                repeats,
                lambda tool: run_arm(tool, work),
                f"{spec.name}/{work}",
            )
            dseams_count = int(stage["dseams_omp1"]["count"])
            ovito_count = int(stage["ovito"]["count"])
            if dseams_count != ovito_count:
                raise SystemExit(
                    f"{spec.name} {work} counts differ: dseams {dseams_count} ovito {ovito_count}"
                )
            print(f"[ring] {spec.name} {work} counts equal {dseams_count}", flush=True)
            for tool in ("dseams_omp1", "ovito"):
                samples = [float(value) for value in stage[tool]["samples_ms"]]
                median = type7(samples, 0.5)
                shown = " ".join(f"{value:.6f}" for value in samples)
                print(
                    f"[ring] {spec.name} {work} {tool} n={len(samples)} "
                    f"median_ms {median:.6f} samples_ms {shown}",
                    flush=True,
                )
            return stage

        # Rings alone: d-SEAMS ring search against OVITO FindRings on bonds
        # built outside the timer. Bonds plus rings times both steps.
        ring_stage = timed_pair("rings")
        bond_stage = timed_pair("bonds")

        out["datasets"][spec.name] = {
            "source": frame.source,
            "n_atoms": n_atoms,
            "triclinic": bool(frame.is_triclinic),
            "rings_total": int(setup_dseams["rings_total"]),
            "rings_size6": int(setup_dseams["rings_size6"]),
            "setup_ms": setup_reads_row,
            "ring_stage_ms": ring_stage,
            "bonds_and_rings_ms": bond_stage,
        }

    Path(args.output_json).write_text(json.dumps(out, indent=2), encoding="utf-8")
    print(f"Wrote {args.output_json}")


if __name__ == "__main__":
    main()
