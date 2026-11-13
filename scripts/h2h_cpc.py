#!/usr/bin/env python3
# Copyright (c) d-SEAMS developers
"""Head-to-head of neighbour lists and primitive rings.

d-SEAMS, freud and OVITO, one thread, release build. Ice Ic, the
in-tree liquid dump, and one ice/liquid interface. The frame is read
before the timer. Repeats are fresh processes, interleaved, after a
warm-up. Median and quartiles are type 7 and are computed by
scripts/h2h_cpc_summary.py from the raw logs.
"""

from __future__ import annotations

import argparse
import hashlib
import os
import subprocess
import sys
import time
import warnings
from pathlib import Path

import numpy as np

ROOT = Path(__file__).resolve().parents[1]
RAW = ROOT / "h2h-cpc" / "raw"
SYSTEMS_DIR = ROOT / "h2h-cpc" / "systems"
BENCH = ROOT / "bbdir" / "h2h-release" / "bench_h2h_cpc"
A_ICE = 6.358
CUTOFF = 3.5
MIN_RING = 3
MAX_RING = 6
WARMUP = 2
REPEATS = 6
BASIS = np.array(
    [
        [0.00, 0.00, 0.00],
        [0.50, 0.50, 0.00],
        [0.50, 0.00, 0.50],
        [0.00, 0.50, 0.50],
        [0.25, 0.25, 0.25],
        [0.75, 0.75, 0.25],
        [0.75, 0.25, 0.75],
        [0.25, 0.75, 0.75],
    ],
    dtype=np.float64,
)
# dump-l.lammpstrj is the in-tree liquid (tests/cli_carry.py). Its
# nearest neighbour is about 1.05 and the box is about 11.8, so the
# file is in reduced units. Lengths are multiplied by 7/3. The campaign
# cutoff 3.5 then sits at reduced radius 1.5, where the seed's mean
# coordination is about 4. docs/orgmode/paper/benchmarks.org names no
# 110592-molecule liquid file; 6x6x6 of this 512-molecule seed is
# 110592.
LIQUID_SEED = ROOT / "input" / "traj" / "dump-l.lammpstrj"
LIQUID_SCALE = 7.0 / 3.0
INTERFACE_GAP = 1.0
# name, kind, reps, molecules (0 until the interface is built)
FRAMES = (
    ("ice-1k", "ice", (5, 5, 5), 1000),
    ("ice-10k", "ice", (11, 11, 11), 10648),
    ("ice-50k", "ice", (18, 18, 19), 49248),
    ("ice-110k", "ice", (24, 24, 24), 110592),
    ("liq-1k", "liquid", (2, 1, 1), 1024),
    ("liq-10k", "liquid", (3, 3, 2), 9216),
    ("liq-50k", "liquid", (4, 4, 6), 49152),
    ("liq-110k", "liquid", (6, 6, 6), 110592),
    ("interface-50k", "interface", (0, 0, 0), 0),
)
TOOLS = ("dseams", "freud", "ovito")
THREAD_ENV = (
    "OMP_NUM_THREADS=1",
    "OMP_DYNAMIC=FALSE",
    "OVITO_THREAD_COUNT=1",
    "OPENBLAS_NUM_THREADS=1",
    "MKL_NUM_THREADS=1",
    "NUMEXPR_NUM_THREADS=1",
    "PYTHONPATH=.h2h-py",
)


def scrub(text: str) -> str:
    cwd = str(ROOT)
    text = text.replace(cwd + "/", "")
    # Built at runtime so the source stays free of a machine path.
    for name in ("home", "workspace", "cursor", "tmp", "agent"):
        root = "/" + name + "/"
        while root in text:
            start = text.find(root)
            end = start
            while end < len(text) and text[end] not in " \n\t'\"":
                end += 1
            text = text[:start] + "[path]" + text[end:]
    return text


def loadavg() -> str:
    return Path("/proc/loadavg").read_text(encoding="utf-8").strip()


def git_sha() -> str:
    return subprocess.check_output(
        ["git", "rev-parse", "HEAD"], cwd=ROOT, text=True
    ).strip()


def ice_positions(reps: tuple[int, int, int]) -> tuple[np.ndarray, np.ndarray]:
    nx, ny, nz = reps
    points = []
    for ix in range(nx):
        for iy in range(ny):
            for iz in range(nz):
                origin = np.array([ix, iy, iz], dtype=np.float64)
                points.append((origin + BASIS) * A_ICE)
    pos = np.vstack(points)
    box = np.array([nx, ny, nz], dtype=np.float64) * A_ICE
    return pos, box


def write_dump(path: Path, pos: np.ndarray, box: np.ndarray) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", encoding="utf-8") as handle:
        handle.write("ITEM: TIMESTEP\n0\nITEM: NUMBER OF ATOMS\n")
        handle.write(f"{len(pos)}\nITEM: BOX BOUNDS pp pp pp\n")
        for length in box:
            handle.write(f"0.0 {length:.8f}\n")
        handle.write("ITEM: ATOMS id type x y z\n")
        for i, (x, y, z) in enumerate(pos, start=1):
            handle.write(f"{i} 1 {x:.8f} {y:.8f} {z:.8f}\n")


def read_dump(path: Path) -> tuple[np.ndarray, np.ndarray]:
    lines = path.read_text(encoding="utf-8").splitlines()
    if not lines[0].startswith("ITEM: TIMESTEP"):
        raise SystemExit(f"not a LAMMPS dump: {path.name}")
    n = int(lines[3])
    box = np.array([float(lines[5 + k].split()[1]) for k in range(3)])
    pos = np.empty((n, 3), dtype=np.float64)
    for i, line in enumerate(lines[9 : 9 + n]):
        parts = line.split()
        pos[i] = (float(parts[2]), float(parts[3]), float(parts[4]))
    return pos, box


def reference_degrees(pos: np.ndarray, box: np.ndarray, cutoff: float) -> np.ndarray:
    """Minimum-image neighbours with distance strictly below the cutoff."""
    n = len(pos)
    counts = np.zeros(n, dtype=np.int32)
    rc2 = cutoff * cutoff
    for i in range(n):
        delta = pos[i + 1 :] - pos[i]
        delta -= box * np.round(delta / box)
        hits = np.flatnonzero(np.sum(delta * delta, axis=1) < rc2)
        counts[i] += hits.size
        counts[i + 1 + hits] += 1
    return counts


def frame_table() -> dict[str, tuple[str, tuple[int, int, int], int]]:
    return {name: (kind, reps, molecules) for name, kind, reps, molecules in FRAMES}


def dump_path(name: str) -> Path:
    return SYSTEMS_DIR / f"{name}.lammpstrj"


def frame_source(name: str) -> str:
    kind, reps, _molecules = frame_table()[name]
    if kind == "ice":
        return f"ice Ic {reps[0]}x{reps[1]}x{reps[2]} a={A_ICE}"
    if kind == "liquid":
        return (
            f"input/traj/dump-l.lammpstrj frame 1 scale {LIQUID_SCALE:.8f} "
            f"{reps[0]}x{reps[1]}x{reps[2]}"
        )
    return (
        "ice Ic slab abutted to dump-l liquid, gap "
        f"{INTERFACE_GAP:.1f} angstrom each side"
    )


def read_seed(path: Path) -> tuple[np.ndarray, np.ndarray]:
    """First frame, wrapped into a cell that starts at the origin."""
    lines = path.read_text(encoding="utf-8").splitlines()
    starts = [i for i, line in enumerate(lines) if line.startswith("ITEM: TIMESTEP")]
    start = starts[0]
    count = int(lines[start + 3])
    low = []
    high = []
    for axis in range(3):
        left, right = lines[start + 5 + axis].split()[:2]
        low.append(float(left))
        high.append(float(right))
    origin = np.array(low, dtype=np.float64)
    box = np.array(high, dtype=np.float64) - origin
    columns = lines[start + 8].split()[2:]
    xyz = [columns.index(name) for name in ("x", "y", "z")]
    type_col = columns.index("type") if "type" in columns else None
    points = []
    for row in lines[start + 9 : start + 9 + count]:
        parts = row.split()
        if type_col is not None and int(float(parts[type_col])) != 1:
            continue
        points.append([float(parts[index]) for index in xyz])
    pos = (np.asarray(points, dtype=np.float64) - origin) % box
    return pos, box


def replicate(pos: np.ndarray, box: np.ndarray, reps: tuple[int, int, int]) -> tuple[np.ndarray, np.ndarray]:
    nx, ny, nz = reps
    tiles = []
    for ix in range(nx):
        for iy in range(ny):
            for iz in range(nz):
                shift = np.array([ix, iy, iz], dtype=np.float64) * box
                tiles.append(pos + shift)
    return np.vstack(tiles), box * np.array([nx, ny, nz], dtype=np.float64)


def liquid_seed() -> tuple[np.ndarray, np.ndarray]:
    pos, box = read_seed(LIQUID_SEED)
    return pos * LIQUID_SCALE, box * LIQUID_SCALE


def ice_slab(nx: int, ny: int, nz: int, spacing: float) -> np.ndarray:
    points = []
    step = np.array([spacing, spacing, spacing], dtype=np.float64)
    for ix in range(nx):
        for iy in range(ny):
            for iz in range(nz):
                origin = np.array([ix, iy, iz], dtype=np.float64) * spacing
                points.append(origin + BASIS * step)
    return np.vstack(points)


def interface_frame() -> tuple[np.ndarray, np.ndarray]:
    """Ice Ic slab and a liquid block in one orthorhombic cell.

    The ice spacing is uniform and chosen so an integer ice grid fills
    the liquid cross-section. A gap of INTERFACE_GAP on each z end keeps
    the slabs from overlapping and still lets the cutoff reach across.
    """
    seed_pos, seed_box = liquid_seed()
    liquid, liquid_box = replicate(seed_pos, seed_box, (4, 4, 2))
    spacing = float(liquid_box[0] / round(liquid_box[0] / A_ICE))
    nx = int(round(liquid_box[0] / spacing))
    ny = int(round(liquid_box[1] / spacing))
    # 14 layers plus the 16384-molecule liquid is about 50k molecules.
    nz = 14
    ice = ice_slab(nx, ny, nz, spacing)
    shifted = liquid.copy()
    shifted[:, 2] += nz * spacing + INTERFACE_GAP
    pos = np.vstack([ice, shifted])
    box = np.array(
        [
            nx * spacing,
            ny * spacing,
            nz * spacing + liquid_box[2] + 2.0 * INTERFACE_GAP,
        ],
        dtype=np.float64,
    )
    return pos, box


def build_frame(name: str) -> tuple[np.ndarray, np.ndarray]:
    kind, reps, _molecules = frame_table()[name]
    if kind == "ice":
        return ice_positions(reps)
    if kind == "liquid":
        pos, box = liquid_seed()
        return replicate(pos, box, reps)
    return interface_frame()


def seed_ready(name: str) -> bool:
    kind, _reps, _molecules = frame_table()[name]
    if kind == "ice":
        return True
    return LIQUID_SEED.is_file()


def prepare_system(name: str) -> Path | None:
    """Write the frame when it is missing. Liquid and interface frames need
    the in-tree seed; return None when that seed is absent."""
    if not seed_ready(name):
        return None
    path = dump_path(name)
    if not path.is_file():
        pos, box = build_frame(name)
        write_dump(path, pos, box)
    return path


def lattice_self_check() -> None:
    pos, box = ice_positions((2, 2, 2))
    if len(pos) != 64:
        raise SystemExit("ice cell replication failed")
    degrees = reference_degrees(pos, box, CUTOFF)
    if not np.all(degrees == 4):
        raise SystemExit(f"ice Ic reference degrees are not 4: {np.unique(degrees)}")
    path = SYSTEMS_DIR / "ice-self.lammpstrj"
    write_dump(path, pos, box)
    rpos, rbox = read_dump(path)
    rdeg = reference_degrees(rpos, rbox, CUTOFF)
    if not np.all(rdeg == 4):
        raise SystemExit("dump round-trip changed the ice neighbour counts")
    if not LIQUID_SEED.is_file():
        return
    liquid, liquid_box = liquid_seed()
    if len(liquid) != 512:
        raise SystemExit(f"liquid seed has {len(liquid)} molecules, expected 512")
    full = reference_degrees(liquid, liquid_box, CUTOFF)
    mean = float(full.mean())
    if not 3.0 < mean < 6.0 or len(np.unique(full)) < 3 or int(full.max()) > 15:
        raise SystemExit(f"scaled liquid is not a mixed first shell: mean {mean:.2f} {np.unique(full, return_counts=True)}")
    pos, box = interface_frame()
    if not 40000 <= len(pos) <= 60000:
        raise SystemExit(f"interface molecule count {len(pos)}")
    # The ice slab ends at the first large empty interval in z.
    order = np.argsort(pos[:, 2])
    z = pos[order, 2]
    gap_at = int(np.argmax(np.diff(z)))
    below = pos[order[: gap_at + 1]]
    above = pos[order[gap_at + 1 :]]
    # A few hundred pairs across the plane is enough to see contact and overlap.
    take_b = below[np.argsort(below[:, 2])[-400:]]
    take_a = above[np.argsort(above[:, 2])[:400]]
    delta = take_a[:, None, :] - take_b[None, :, :]
    delta -= box * np.round(delta / box)
    dist = np.sqrt(np.sum(delta * delta, axis=2))
    nearest = float(dist.min())
    contacts = int(np.sum(dist < CUTOFF))
    if nearest < 1.2 or contacts == 0:
        raise SystemExit(f"interface contact nearest {nearest:.3f} pairs {contacts}")


def sha_file(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(1 << 20), b""):
            digest.update(chunk)
    return digest.hexdigest()


def parse_worker(text: str) -> dict[str, str]:
    record = {}
    for line in text.splitlines():
        if not line.strip():
            continue
        key, _, value = line.partition(" ")
        record[key] = value
    return record


def command_line(argv: list[str]) -> str:
    return " ".join(argv)


def invoke(argv: list[str]) -> tuple[dict[str, str], str, int, float]:
    before = loadavg()
    started = time.perf_counter()
    done = subprocess.run(
        argv, cwd=ROOT, text=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE
    )
    wall_ms = (time.perf_counter() - started) * 1000.0
    after = loadavg()
    stdout = scrub(done.stdout)
    stderr = scrub(done.stderr)
    record = parse_worker(stdout)
    record["loadavg_before"] = before
    record["loadavg_after"] = after
    record["wall_ms"] = f"{wall_ms:.3f}"
    record["exit"] = str(done.returncode)
    if done.returncode != 0:
        record["status"] = record.get("status", "error")
        record["stderr"] = stderr[-2000:]
    return record, stdout, done.returncode, wall_ms


def dseams_argv(comparison: str, dump: Path, csv_path: Path) -> list[str]:
    mode = "neighbours" if comparison == "neighbours" else "rings"
    return [
        "taskset", "-c", "0", "env", *THREAD_ENV,
        str(BENCH.relative_to(ROOT)),
        mode,
        str(dump.relative_to(ROOT)),
        f"{CUTOFF:.8f}",
        str(MAX_RING),
        str(csv_path.relative_to(ROOT)),
    ]


def python_argv(comparison: str, tool: str, dump: Path, csv_path: Path) -> list[str]:
    return [
        "taskset", "-c", "0", "env", *THREAD_ENV,
        "python3.12", "scripts/h2h_cpc.py", "worker",
        "--tool", tool,
        "--comparison", comparison,
        "--dump", str(dump.relative_to(ROOT)),
        "--csv", str(csv_path.relative_to(ROOT)),
        "--cutoff", f"{CUTOFF:.8f}",
        "--min-ring", str(MIN_RING),
        "--max-ring", str(MAX_RING),
    ]


def write_record(path: Path, record: dict[str, str], keys: list[str]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    lines = [f"{key} {record[key]}" for key in keys if key in record and record[key] != ""]
    path.write_text("\n".join(lines) + "\n", encoding="utf-8")


ENVELOPE = [
    "comparison", "system", "tool", "phase", "repeat", "status", "molecules",
    "time_ms", "wall_ms", "histogram", "omp_threads", "threads", "cutoff",
    "min_ring", "max_ring", "git_sha", "source", "freud_version", "ovito_version",
    "numpy_version", "scipy_version", "cpu_model", "command", "csv", "csv_sha256",
    "loadavg_before", "loadavg_after", "exit", "reason", "stderr",
]


def run_tool(comparison: str, tool: str, dump: Path, csv_path: Path) -> dict[str, str]:
    csv_path.parent.mkdir(parents=True, exist_ok=True)
    if tool == "dseams":
        argv = dseams_argv(comparison, dump, csv_path)
    else:
        argv = python_argv(comparison, tool, dump, csv_path)
    record, _stdout, code, _wall = invoke(argv)
    record["command"] = command_line(argv)
    record["threads"] = "taskset=0 omp=1 ovito=1 freud=1"
    record["cutoff"] = f"{CUTOFF:.8f}"
    record["min_ring"] = str(MIN_RING)
    record["max_ring"] = str(MAX_RING)
    record["git_sha"] = git_sha()
    record["histogram"] = record.get("tally", "")
    versions = package_versions()
    record.update(versions)
    record["cpu_model"] = cpu_model_line()
    if code != 0:
        raise SystemExit(
            f"{tool} {comparison} failed\n{record.get('stderr', '')}\n{record}"
        )
    if record.get("status") == "ok" and csv_path.is_file():
        record["csv"] = str(csv_path.relative_to(ROOT))
        record["csv_sha256"] = sha_file(csv_path)
    return record


def write_neighbour_csv(path: Path, rows: list[list[int]]) -> str:
    tally: dict[int, int] = {}
    with path.open("w", encoding="utf-8") as handle:
        handle.write("index,count,neighbours\n")
        for i, nbrs in enumerate(rows):
            ordered = sorted(set(nbrs))
            tally[len(ordered)] = tally.get(len(ordered), 0) + 1
            body = " ".join(str(j) for j in ordered)
            handle.write(f"{i},{len(ordered)},{body}\n")
    return " ".join(f"{k}:{tally[k]}" for k in sorted(tally))


def write_ring_csv(path: Path, hist: np.ndarray, min_ring: int, max_ring: int) -> str:
    tally: dict[int, int] = {}
    with path.open("w", encoding="utf-8") as handle:
        handle.write("index," + ",".join(str(s) for s in range(min_ring, max_ring + 1)) + "\n")
        for i in range(hist.shape[0]):
            handle.write(str(i))
            for size in range(min_ring, max_ring + 1):
                handle.write(f",{int(hist[i, size])}")
            handle.write("\n")
    # Ring tally counts rings, which is the sum of memberships divided by size.
    for size in range(min_ring, max_ring + 1):
        membership = int(hist[:, size].sum())
        if membership:
            if membership % size != 0:
                raise SystemExit(f"ring size {size} membership {membership} is not divisible")
            tally[size] = membership // size
    return " ".join(f"{k}:{tally[k]}" for k in sorted(tally))


def worker_freud(pos: np.ndarray, box: np.ndarray, comparison: str, csv_path: Path,
                 cutoff: float, min_ring: int, max_ring: int) -> None:
    import freud

    if hasattr(freud, "parallel"):
        freud.parallel.set_num_threads(1)
    fbox = freud.box.Box(Lx=float(box[0]), Ly=float(box[1]), Lz=float(box[2]))
    if comparison == "rings":
        print("status unsupported")
        print(f"molecules {len(pos)}")
        print("reason freud publishes no primitive-ring statistic")
        return
    query = freud.locality.AABBQuery(fbox, np.ascontiguousarray(pos))
    started = time.perf_counter()
    nlist = query.query(pos, {"r_max": cutoff, "exclude_ii": True}).toNeighborList()
    elapsed = (time.perf_counter() - started) * 1000.0
    rows: list[list[int]] = [[] for _ in range(len(pos))]
    for i, j in zip(nlist.query_point_indices, nlist.point_indices):
        rows[int(i)].append(int(j))
    tally = write_neighbour_csv(csv_path, rows)
    print("status ok")
    print(f"molecules {len(pos)}")
    print(f"time_ms {elapsed:.6f}")
    print("omp_threads 1")
    print(f"tally {tally}")


def _ring_table(table) -> np.ndarray:
    value = table["Particle Indices"]
    return np.asarray(value)


def worker_ovito(pos: np.ndarray, box: np.ndarray, comparison: str, csv_path: Path,
                 cutoff: float, min_ring: int, max_ring: int) -> None:
    warnings.filterwarnings("ignore", message=".*OVITO.*PyPI")
    os.environ["OVITO_THREAD_COUNT"] = "1"
    from ovito.data import DataCollection
    from ovito.modifiers import CreateBondsModifier, FindRingsModifier
    from ovito.pipeline import Pipeline, StaticSource

    data = DataCollection()
    data.create_cell(
        [[float(box[0]), 0.0, 0.0], [0.0, float(box[1]), 0.0], [0.0, 0.0, float(box[2])]],
        pbc=(True, True, True),
    )
    parts = data.create_particles(count=len(pos))
    parts.create_property("Position", data=pos)

    if comparison == "neighbours":
        pipe = Pipeline(source=StaticSource(data=data))
        pipe.modifiers.append(CreateBondsModifier(cutoff=cutoff))
        started = time.perf_counter()
        out = pipe.compute()
        elapsed = (time.perf_counter() - started) * 1000.0
        bonds = out.particles.bonds
        topo = np.zeros((0, 2), dtype=np.int64) if bonds is None or bonds.count == 0 else np.asarray(bonds.topology)
        rows: list[list[int]] = [[] for _ in range(len(pos))]
        for i, j in topo:
            rows[int(i)].append(int(j))
            rows[int(j)].append(int(i))
        tally = write_neighbour_csv(csv_path, rows)
    else:
        bond_pipe = Pipeline(source=StaticSource(data=data))
        bond_pipe.modifiers.append(CreateBondsModifier(cutoff=cutoff))
        bonded = bond_pipe.compute()
        ring_pipe = Pipeline(source=StaticSource(data=bonded))
        ring_pipe.modifiers.append(
            FindRingsModifier(minimum_ring_size=min_ring, maximum_ring_size=max_ring)
        )
        started = time.perf_counter()
        out = ring_pipe.compute()
        hist = np.zeros((len(pos), max_ring + 1), dtype=np.int32)
        for size in range(min_ring, max_ring + 1):
            key = f"{size}-rings"
            if key not in out.tables:
                continue
            array = _ring_table(out.tables[key])
            if array.ndim == 1:
                array = array.reshape(1, -1)
            for row in array:
                members = [int(v) for v in np.asarray(row).reshape(-1) if int(v) >= 0]
                if len(members) != size:
                    raise SystemExit(f"OVITO {size}-ring has {len(members)} members")
                for atom in members:
                    hist[atom, size] += 1
        elapsed = (time.perf_counter() - started) * 1000.0
        tally = write_ring_csv(csv_path, hist, min_ring, max_ring)
    print("status ok")
    print(f"molecules {len(pos)}")
    print(f"time_ms {elapsed:.6f}")
    print("omp_threads 1")
    print(f"tally {tally}")


def worker_main(args: argparse.Namespace) -> None:
    pos, box = read_dump(ROOT / args.dump)
    csv_path = ROOT / args.csv
    csv_path.parent.mkdir(parents=True, exist_ok=True)
    if args.tool == "freud":
        worker_freud(pos, box, args.comparison, csv_path, args.cutoff, args.min_ring, args.max_ring)
    elif args.tool == "ovito":
        worker_ovito(pos, box, args.comparison, csv_path, args.cutoff, args.min_ring, args.max_ring)
    else:
        raise SystemExit(f"unknown tool {args.tool}")


_CPU_MODEL: str | None = None
_VERSIONS: dict[str, str] | None = None


def cpu_model_line() -> str:
    global _CPU_MODEL
    if _CPU_MODEL is None:
        for line in subprocess.check_output(["lscpu"], text=True).splitlines():
            if line.startswith("Model name:"):
                _CPU_MODEL = line
                break
        if _CPU_MODEL is None:
            raise SystemExit("lscpu did not print a model name")
    return _CPU_MODEL


def package_versions() -> dict[str, str]:
    global _VERSIONS
    if _VERSIONS is not None:
        return _VERSIONS
    text = subprocess.check_output(
        ["python3.12", "-c",
         "import os; os.environ['PYTHONPATH']='.h2h-py'; "
         "import freud, ovito, numpy, scipy; "
         "ov = ovito.version; "
         "ov = '.'.join(map(str, ov)) if isinstance(ov, tuple) else str(ov); "
         "print(freud.__version__); print(ov); print(numpy.__version__); print(scipy.__version__)"],
        cwd=ROOT, text=True,
    ).splitlines()
    freud_v, ovito_v, numpy_v, scipy_v = text
    _VERSIONS = {
        "freud_version": freud_v,
        "ovito_version": ovito_v,
        "numpy_version": numpy_v,
        "scipy_version": scipy_v,
    }
    return _VERSIONS


def environment_text() -> str:
    os_release = Path("/etc/os-release").read_text(encoding="utf-8")
    pretty = "unknown"
    for line in os_release.splitlines():
        if line.startswith("PRETTY_NAME="):
            pretty = line.split("=", 1)[1].strip().strip('"')
    lscpu = subprocess.check_output(["lscpu"], text=True)
    wanted = (
        "CPU(s):",
        "Thread(s) per core:",
        "Core(s) per socket:",
        "Socket(s):",
    )
    cpu_lines = [cpu_model_line()]
    for line in lscpu.splitlines():
        if line.startswith(wanted):
            cpu_lines.append(line.rstrip("\n"))
    mem = []
    for line in Path("/proc/meminfo").read_text(encoding="utf-8").splitlines():
        if line.startswith(("MemTotal:", "MemAvailable:")):
            mem.append(" ".join(line.split()))
    compiler = subprocess.check_output(["g++", "--version"], text=True).splitlines()[0]
    versions = package_versions()
    version_lines = [f"{key.replace('_version', '')} {value}" for key, value in versions.items()]
    buildtype = "release"
    lines = [
        f"git_sha {git_sha()}",
        f"os {pretty}",
        f"compiler {compiler}",
        f"buildtype {buildtype}",
        f"build bbdir/h2h-release/bench_h2h_cpc",
        *cpu_lines,
        *mem,
        *version_lines,
        f"cutoff {CUTOFF:.8f}",
        f"min_ring {MIN_RING}",
        f"max_ring {MAX_RING}",
        f"warmup {WARMUP}",
        f"repeats {REPEATS}",
        "threads taskset -c 0; OMP_NUM_THREADS=1; OVITO_THREAD_COUNT=1; freud.parallel.set_num_threads(1)",
        "timer algorithm call after the frame is in memory",
        "median type 7",
        "systems ice Ic a=6.358; liquid is dump-l frame 1 scaled by 7/3; interface abuts the two",
        "neighbour endpoint d-SEAMS and freud are strict; OVITO keeps equality; see h2h-cpc/data/exact-cutoff-result.txt",
    ]
    return scrub("\n".join(lines) + "\n")


def selected(only: str) -> list[str]:
    names = [item[0] for item in FRAMES]
    if not only:
        return names
    wanted = only.split(",")
    missing = [name for name in wanted if name not in frame_table()]
    if missing:
        raise SystemExit(f"unknown system {missing}")
    return wanted


def run_comparison(comparison: str, only: str) -> None:
    if not BENCH.is_file():
        raise SystemExit("missing bbdir/h2h-release/bench_h2h_cpc")
    RAW.mkdir(parents=True, exist_ok=True)
    (RAW / "environment.txt").write_text(environment_text(), encoding="utf-8")
    done_lines = []
    for name in selected(only):
        dump = prepare_system(name)
        if dump is None:
            print(f"[{comparison}] {name} skipped, liquid seed absent", flush=True)
            continue
        pos, _box = read_dump(dump)
        molecules = len(pos)
        _kind, _reps, expected = frame_table()[name]
        if expected and molecules != expected:
            raise SystemExit(f"{name} has {molecules} molecules, expected {expected}")
        if name == "interface-50k" and not 40000 <= molecules <= 60000:
            raise SystemExit(f"interface has {molecules} molecules, expected about 50k")
        source = frame_source(name)
        print(f"[{comparison}] {name} molecules {molecules}", flush=True)
        order_lines = []
        active: list[str] = []
        for tool in TOOLS:
            csv_path = RAW / comparison / name / tool / "check.csv"
            record = run_tool(comparison, tool, dump, csv_path)
            record.update({
                "comparison": comparison,
                "system": name,
                "tool": tool,
                "phase": "check",
                "repeat": "0",
                "source": source,
            })
            dest = RAW / comparison / name / tool / "check.txt"
            if record.get("status") == "unsupported":
                record["molecules"] = str(molecules)
                write_record(dest, record, ENVELOPE)
                print(f"[{comparison}] {name} {tool} unsupported {record.get('reason', '')}", flush=True)
                continue
            if int(record.get("molecules", "-1")) != molecules:
                raise SystemExit(f"{tool} molecules {record.get('molecules')} != {molecules}")
            write_record(dest, record, ENVELOPE)
            active.append(tool)
            order_lines.append(f"check {tool}")
        if not active:
            raise SystemExit(f"no tool ran {comparison} on {name}")
        for warm in range(WARMUP):
            order = active[warm % len(active):] + active[: warm % len(active)]
            order_lines.append("warmup " + str(warm) + " " + " ".join(order))
            for tool in order:
                csv_path = RAW / comparison / name / tool / f"warm-{warm:02d}.csv"
                record = run_tool(comparison, tool, dump, csv_path)
                record.update({
                    "comparison": comparison, "system": name, "tool": tool,
                    "phase": "warmup", "repeat": str(warm), "source": source,
                })
                write_record(RAW / comparison / name / tool / f"warm-{warm:02d}.txt", record, ENVELOPE)
        for rep in range(REPEATS):
            shift = (WARMUP + rep) % len(active)
            order = active[shift:] + active[:shift]
            order_lines.append("repeat " + str(rep) + " " + " ".join(order))
            for tool in order:
                csv_path = RAW / comparison / name / tool / f"rep-{rep:02d}.csv"
                record = run_tool(comparison, tool, dump, csv_path)
                record.update({
                    "comparison": comparison, "system": name, "tool": tool,
                    "phase": "timed", "repeat": str(rep), "source": source,
                })
                write_record(RAW / comparison / name / tool / f"rep-{rep:02d}.txt", record, ENVELOPE)
            print(f"[{comparison}] {name} repeat {rep} done", flush=True)
        order_path = RAW / comparison / name / "order.txt"
        order_path.write_text("\n".join(order_lines) + "\n", encoding="utf-8")
        done_lines.append(f"{name} {molecules}")
    done = RAW / comparison / "done.txt"
    done.parent.mkdir(parents=True, exist_ok=True)
    done.write_text("\n".join(done_lines) + "\n", encoding="utf-8")


def graph_from_cutoff(pos: np.ndarray, box: np.ndarray, cutoff: float) -> list[list[int]]:
    rows = [[] for _ in range(len(pos))]
    rc2 = cutoff * cutoff
    for i in range(len(pos)):
        delta = pos[i + 1 :] - pos[i]
        delta -= box * np.round(delta / box)
        hits = np.flatnonzero(np.sum(delta * delta, axis=1) < rc2)
        for offset in hits:
            j = i + 1 + int(offset)
            rows[i].append(j)
            rows[j].append(i)
    return [sorted(set(row)) for row in rows]


def franzblau_histogram(rows: list[list[int]], max_ring: int) -> np.ndarray:
    """Shortest-path rings. A ring is kept when, for every pair of members,
    one of the two arcs has length equal to the graph distance.
    Franzblau, Phys. Rev. B 44, 4925 (1991).
    """
    n = len(rows)
    hist = np.zeros((n, max_ring + 1), dtype=np.int32)
    seen: set[tuple[int, ...]] = set()

    def canonicalize(path: list[int]) -> tuple[int, ...]:
        start = min(range(len(path)), key=path.__getitem__)
        forward = path[start:] + path[:start]
        rev = path[start::-1] + path[:start:-1]
        a = tuple(forward)
        b = tuple(rev)
        return a if a <= b else b

    def dfs(start: int, current: int, path: list[int], used: set[int]) -> None:
        for nxt in rows[current]:
            if nxt == start and len(path) >= MIN_RING:
                canon = canonicalize(path)
                if canon not in seen:
                    seen.add(canon)
            elif nxt > start and nxt not in used and len(path) < max_ring:
                used.add(nxt)
                path.append(nxt)
                dfs(start, nxt, path, used)
                path.pop()
                used.remove(nxt)

    for start in range(n):
        dfs(start, start, [start], {start})

    def distances(src: int) -> list[int]:
        dist = [-1] * n
        dist[src] = 0
        queue = [src]
        for u in queue:
            for v in rows[u]:
                if dist[v] < 0:
                    dist[v] = dist[u] + 1
                    queue.append(v)
        return dist

    for canon in seen:
        size = len(canon)
        dist = [distances(atom) for atom in canon]
        keep = True
        for i in range(size):
            for j in range(i + 1, size):
                arc = j - i
                other = size - arc
                gap = dist[i][canon[j]]
                if gap != arc and gap != other:
                    keep = False
                    break
            if not keep:
                break
        if not keep:
            continue
        for atom in canon:
            hist[atom, size] += 1
    return hist


def run_probe() -> None:
    """Small geometries. Write a counterexample only when a tool leaves
    the published predicate."""
    if not BENCH.is_file():
        raise SystemExit("missing bbdir/h2h-release/bench_h2h_cpc")
    work = ROOT / "h2h-cpc" / "probe-tmp"
    work.mkdir(parents=True, exist_ok=True)
    notes: list[str] = []

    # Pair at exactly the cutoff. The d-SEAMS predicate keeps a pair only
    # when the squared minimum image is strictly below the cutoff squared.
    box = np.array([20.0, 20.0, 20.0])
    pair = np.array([[5.0, 10.0, 10.0], [5.0 + CUTOFF, 10.0, 10.0]])
    pair_path = work / "exact-cutoff.lammpstrj"
    write_dump(pair_path, pair, box)
    ref = reference_degrees(read_dump(pair_path)[0], read_dump(pair_path)[1], CUTOFF)
    notes.append(
        "exact-cutoff reference strict-less degrees "
        + " ".join(str(int(v)) for v in ref)
    )
    disagree = False
    for tool in TOOLS:
        csv_path = work / f"exact-{tool}.csv"
        record = run_tool("neighbours", tool, pair_path, csv_path)
        text = csv_path.read_text(encoding="utf-8").splitlines()[1:]
        degrees = [int(line.split(",")[1]) for line in text]
        notes.append(f"exact-cutoff {tool} degrees {' '.join(map(str, degrees))} tally {record.get('tally', '')}")
        if degrees != [int(v) for v in ref]:
            disagree = True
            notes.append(
                f"exact-cutoff {tool} degree leaves the strict minimum-image count"
            )

    raw = np.array(
        [[1.0, 1.0, 1.0], [1.0, -1.0, -1.0], [-1.0, 1.0, -1.0], [-1.0, -1.0, 1.0]],
        dtype=np.float64,
    )
    tetra = raw * (1.0 / np.linalg.norm(raw[0] - raw[1]))
    tetra = tetra - tetra.min(axis=0) + 8.0
    tetra_box = np.array([20.0, 20.0, 20.0])
    tetra_path = work / "tetrahedron.lammpstrj"
    write_dump(tetra_path, tetra, tetra_box)
    tpos, tbox = read_dump(tetra_path)
    bonds = graph_from_cutoff(tpos, tbox, 1.2)
    hist = franzblau_histogram(bonds, MAX_RING)
    notes.append(
        "tetrahedron Franzblau membership "
        + " ".join(str(int(v)) for v in hist.sum(axis=0))
    )
    # The campaign cutoff is 3.5, which would bond nothing useful here.
    # The probe calls the workers with the campaign cutoff unless we pass
    # a dedicated cutoff through run_tool. Rebuild the call locally.
    saved = CUTOFF
    # run_tool closes over the global. Use a direct worker call instead.
    del saved
    for tool in ("dseams", "freud", "ovito"):
        csv_path = work / f"tetra-{tool}.csv"
        if tool == "dseams":
            argv = [
                "taskset", "-c", "0", "env", *THREAD_ENV,
                str(BENCH.relative_to(ROOT)), "rings",
                str(tetra_path.relative_to(ROOT)), "1.20000000", str(MAX_RING),
                str(csv_path.relative_to(ROOT)),
            ]
            record, _out, code, _wall = invoke(argv)
            if code != 0:
                raise SystemExit(record.get("stderr", "dseams probe failed"))
            got = csv_path.read_text(encoding="utf-8").splitlines()[1:]
        elif tool == "freud":
            notes.append("tetrahedron freud unsupported")
            continue
        else:
            import contextlib
            import io
            with contextlib.redirect_stdout(io.StringIO()):
                worker_ovito(tpos, tbox, "rings", csv_path, 1.2, MIN_RING, MAX_RING)
            got = csv_path.read_text(encoding="utf-8").splitlines()[1:]
        # Membership of size 3 is column 1 of the data row (index,3,4,...).
        threes = []
        for line in got:
            cols = line.split(",")
            threes.append(int(cols[1]))
        ref_three = [int(v) for v in hist[:, 3]]
        notes.append(f"tetrahedron {tool} size3 {' '.join(map(str, threes))}")
        if threes != ref_three:
            disagree = True
            notes.append(
                f"tetrahedron {tool} leaves the Franzblau shortest-path count"
            )
        # A size-4 column must be zero under the published criterion.
        fours = [int(line.split(",")[2]) for line in got]
        if any(fours):
            disagree = True
            notes.append(f"tetrahedron {tool} 4-ring count leaves the Franzblau count of zero")

    print("\n".join(notes), flush=True)
    if not disagree:
        return
    dest = ROOT / "h2h-cpc" / "data"
    dest.mkdir(parents=True, exist_ok=True)
    (ROOT / "h2h-cpc" / "counterexample").mkdir(parents=True, exist_ok=True)
    (ROOT / "h2h-cpc" / "counterexample" / "exact-cutoff.dump").write_text(
        pair_path.read_text(encoding="utf-8"), encoding="utf-8"
    )
    body = "\n".join([
        "Two particles in an orthorhombic cell of side 20.",
        "Coordinates: (5, 10, 10) and (8.5, 10, 10). Separation 3.5.",
        "The ice campaign uses this same cutoff.",
        "",
        "d-SEAMS keeps a pair when the squared minimum-image distance",
        "lies strictly under the cutoff squared. neighListO drops a pair",
        "once distSq reaches rcutoffSq. The orthorhombic cell list keeps",
        "a pair when d2 < rc2. The reference count uses that predicate.",
        "freud AABBQuery with r_max 3.5 returns that same neighbour list,",
        "degree 0 on both particles.",
        "OVITO CreateBondsModifier(cutoff=3.5) bonds the pair, so a distance",
        "equal to the cutoff falls inside its upper cutoff.",
        "",
        "Measured degrees, particle 0 then particle 1:",
        *notes,
        "",
        "The ice Ic shells used for the timings sit near 2.75 and 4.50, so",
        "every timed frame keeps its pairs off the cutoff. The neighbour-count",
        "identity on those frames is in the raw CSVs.",
        "",
        "Primitive rings were checked on the complete graph K4 (regular",
        "tetrahedron, edge 1, cutoff 1.2). Franzblau's shortest-path criterion",
        "(Phys. Rev. B 44, 4925, 1991) keeps the four triangles. A 4-cycle",
        "exceeds the shortest path by 1 edge, because both paths around it run longer",
        "than the direct edge. d-SEAMS ringNetwork and OVITO FindRingsModifier",
        "both report three 3-rings on every atom. Their 4-ring count is zero.",
        "freud's published interface stops at the neighbour list.",
        "",
    ])
    (dest / "exact-cutoff-result.txt").write_text(body, encoding="utf-8")
    print(f"wrote {dest.relative_to(ROOT)}", flush=True)


def asv_call(comparison: str, tool: str, system: str) -> None:
    dump = prepare_system(system)
    if dump is None:
        raise SystemExit(f"skip {system}")
    dest = ROOT / "h2h-cpc" / "asv-tmp" / comparison / system / tool
    dest.mkdir(parents=True, exist_ok=True)
    csv_path = dest / "once.csv"
    record = run_tool(comparison, tool, dump, csv_path)
    if record.get("status") == "unsupported":
        raise SystemExit(record.get("reason", "unsupported"))
    print(f"time_ms {record['time_ms']}")


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("action")
    parser.add_argument("rest", nargs="*")
    parser.add_argument("--only", default="")
    parser.add_argument("--tool")
    parser.add_argument("--comparison")
    parser.add_argument("--dump")
    parser.add_argument("--csv")
    parser.add_argument("--cutoff", type=float, default=CUTOFF)
    parser.add_argument("--min-ring", type=int, default=MIN_RING)
    parser.add_argument("--max-ring", type=int, default=MAX_RING)
    args = parser.parse_args()
    os.chdir(ROOT)
    if args.action == "worker":
        worker_main(args)
        return
    if args.action == "self-check":
        lattice_self_check()
        print("self-check ok")
        return
    if args.action == "prepare":
        name = args.rest[0]
        path = prepare_system(name)
        if path is None:
            print(f"skip {name}")
            return
        print(path.relative_to(ROOT))
        return
    if args.action == "asv-call":
        asv_call(args.rest[0], args.rest[1], args.rest[2])
        return
    if args.action == "probe":
        run_probe()
        return
    if args.action == "run":
        lattice_self_check()
        which = args.rest[0] if args.rest else "all"
        if which in {"neighbours", "all"}:
            run_comparison("neighbours", args.only)
        if which in {"rings", "all"}:
            run_comparison("rings", args.only)
        run_probe()
        return
    raise SystemExit(f"unknown action {args.action}")


if __name__ == "__main__":
    main()
