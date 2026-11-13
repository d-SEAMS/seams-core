#!/usr/bin/env python3
"""Cage head-to-head on GenIce sI, sII and sH.

Frames are written under benchmarks/genice-frames and are not tracked.
GenIce builds each cell. d-SEAMS reads that dump. ClavisCage TRACE and
cycless 0.7 read the same coordinates from memory.

d-SEAMS signatures are 512, 51262 and 51264, one ring search in
cmdCages. The library call is seams_frame_cages on the seeded graph.
cycless follows GenIce's cage assessment: cycles up to 8 with fractional
positions, so a spanning cycle is dropped, then polyhedra up to 22
faces. quick is used when it returns the same cages as the full search.
ClavisCage is TRACE at 0.35 nm and 35 degrees, with hydrogens on an
ice-rule orientation of the oxygen graph. A tilted cell is passed as
six box values so the minimum image keeps the tilt. ClavisCage and
cycless have no thread pool; the process is still pinned and
RAYON_NUM_THREADS is set.
"""

from __future__ import annotations

import argparse
import ctypes
import json
import logging
import os
import re
import subprocess
import sys
import time
from collections import Counter
from pathlib import Path

import numpy as np

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "scripts"))

import paper_stats  # noqa: E402

CUTOFF = 3.5
SIGNATURES = ("512", "51262", "51264")
# Smallest cubic replica of each GenIce cell with at least 110592 molecules,
# plus the unit cell and the 2x2x2 cell.
CELLS = (
    ("sI", (1, 1, 1)),
    ("sI", (2, 2, 2)),
    ("sI", (14, 14, 14)),
    ("sII", (1, 1, 1)),
    ("sII", (2, 2, 2)),
    ("sII", (10, 10, 10)),
    ("sH", (1, 1, 1)),
    ("sH", (2, 2, 2)),
    ("sH", (12, 12, 12)),
)
SIG_LINE = re.compile(r"signature\s+(\S+)\s+cages\s+(\d+)")
SIG_TOKEN = {
    "5:12": "512",
    "5:12,6:2": "51262",
    "5:12,6:4": "51264",
}


def cpu_model() -> str:
    with open("/proc/cpuinfo", encoding="utf-8") as handle:
        for line in handle:
            if line.startswith("model name"):
                return line.split(":", 1)[1].strip()
    return ""


def lammps_box(cell: np.ndarray) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    """LAMMPS bounds, lattice lengths and tilt for a row-vector cell."""
    length = np.array([cell[0, 0], cell[1, 1], cell[2, 2]], dtype=np.float64)
    tilt = np.array([cell[1, 0], cell[2, 0], cell[2, 1]], dtype=np.float64)
    xy, xz, yz = (float(v) for v in tilt)
    xmin = min(0.0, xy, xz, xy + xz)
    xmax = max(0.0, xy, xz, xy + xz)
    ymin = min(0.0, yz)
    ymax = max(0.0, yz)
    bound = np.array(
        [length[0] + xmax - xmin, length[1] + ymax - ymin, length[2]],
        dtype=np.float64,
    )
    low = np.array([xmin, ymin, 0.0], dtype=np.float64)
    return low, bound, tilt


def generate(name: str, rep: tuple[int, int, int]):
    logging.disable(logging.CRITICAL)
    from genice2.genice import GenIce
    from genice2.plugin import safe_import

    lattice = safe_import("lattice", name).Lattice()
    ice = GenIce(lattice, reshape=np.diag(list(rep)))
    ice.Stage1()
    cell = np.asarray(ice.repcell.mat, dtype=np.float64) * 10.0
    pos = np.ascontiguousarray(np.asarray(ice.reppositions, dtype=np.float64) @ cell)
    return pos, cell


def write_dump(path: Path, pos: np.ndarray, cell: np.ndarray) -> None:
    low, bound, tilt = lammps_box(cell)
    lines = [
        "ITEM: TIMESTEP",
        "0",
        "ITEM: NUMBER OF ATOMS",
        str(len(pos)),
        "ITEM: BOX BOUNDS xy xz yz pp pp pp",
        f"{low[0]:.10f} {low[0] + bound[0]:.10f} {tilt[0]:.10f}",
        f"{low[1]:.10f} {low[1] + bound[1]:.10f} {tilt[1]:.10f}",
        f"{low[2]:.10f} {low[2] + bound[2]:.10f} {tilt[2]:.10f}",
        "ITEM: ATOMS id type x y z",
    ]
    body = "\n".join(
        f"{i + 1} 1 {p[0]:.8f} {p[1]:.8f} {p[2]:.8f}" for i, p in enumerate(pos)
    )
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text("\n".join(lines) + "\n" + body + "\n", encoding="utf-8")


def read_dump(path: Path):
    lines = path.read_text(encoding="utf-8").splitlines()
    count = int(lines[3])
    bounds = []
    tilts = []
    for row in lines[5:8]:
        parts = row.split()
        bounds.append((float(parts[0]), float(parts[1])))
        tilts.append(float(parts[2]))
    low = np.array([bounds[0][0], bounds[1][0], bounds[2][0]])
    span = np.array([bounds[d][1] - bounds[d][0] for d in range(3)])
    xy, xz, yz = tilts
    xmin = min(0.0, xy, xz, xy + xz)
    xmax = max(0.0, xy, xz, xy + xz)
    ymin = min(0.0, yz)
    ymax = max(0.0, yz)
    lx = span[0] - (xmax - xmin)
    ly = span[1] - (ymax - ymin)
    lz = span[2]
    cell = np.array([[lx, 0.0, 0.0], [xy, ly, 0.0], [xz, yz, lz]], dtype=np.float64)
    pos = np.empty((count, 3), dtype=np.float64)
    for i, line in enumerate(lines[9 : 9 + count]):
        parts = line.split()
        pos[i] = (float(parts[2]), float(parts[3]), float(parts[4]))
    return pos, cell


def min_image_edges(pos: np.ndarray, cell: np.ndarray, cutoff: float = CUTOFF):
    """Undirected pairs and the i->j minimum-image vector."""
    n = len(pos)
    inv = np.linalg.inv(cell)
    frac = np.mod(pos @ inv, 1.0)
    lengths = np.linalg.norm(cell, axis=1)
    ncell = np.maximum(1, np.floor(lengths / cutoff).astype(np.int32))
    reach = np.empty(3, dtype=np.int32)
    for dim in range(3):
        width = lengths[dim] / float(ncell[dim])
        reach[dim] = max(1, int(np.ceil(cutoff / width - 1e-12)))
    # A skewed cell can put a contact two bins away along a long diagonal.
    reach = np.minimum(reach + 1, ncell)
    bins = np.floor(frac * ncell).astype(np.int32) % ncell
    buckets: dict[tuple[int, int, int], list[int]] = {}
    for i, b in enumerate(bins):
        buckets.setdefault((int(b[0]), int(b[1]), int(b[2])), []).append(i)
    cut2 = cutoff * cutoff
    edges = []
    vecs: dict[tuple[int, int], np.ndarray] = {}
    for i in range(n):
        bi = bins[i]
        seen = set()
        for ox in range(-int(reach[0]), int(reach[0]) + 1):
            for oy in range(-int(reach[1]), int(reach[1]) + 1):
                for oz in range(-int(reach[2]), int(reach[2]) + 1):
                    key = (
                        int((bi[0] + ox) % ncell[0]),
                        int((bi[1] + oy) % ncell[1]),
                        int((bi[2] + oz) % ncell[2]),
                    )
                    if key in seen:
                        continue
                    seen.add(key)
                    for j in buckets.get(key, ()):
                        if j <= i:
                            continue
                        delta = frac[j] - frac[i]
                        delta -= np.round(delta)
                        cart = delta @ cell
                        if float(cart @ cart) < cut2:
                            edges.append((i, j))
                            vecs[(i, j)] = cart
                            vecs[(j, i)] = -cart
    return edges, vecs


def orient_donors(n: int, edges) -> list[list[int]]:
    remaining = [set() for _ in range(n)]
    for i, j in edges:
        remaining[i].add(j)
        remaining[j].add(i)
    donors = [[] for _ in range(n)]
    for start in range(n):
        if not remaining[start]:
            continue
        stack = [start]
        trail = []
        while stack:
            node = stack[-1]
            if remaining[node]:
                nxt = remaining[node].pop()
                remaining[nxt].discard(node)
                stack.append(nxt)
            else:
                trail.append(stack.pop())
        trail.reverse()
        for a, b in zip(trail, trail[1:]):
            donors[a].append(b)
    return donors


def clavis_waters(pos: np.ndarray, cell: np.ndarray, edges, vecs) -> str:
    donors = orient_donors(len(pos), edges)
    waters = []
    for i, point in enumerate(pos):
        placed = []
        for other in donors[i][:2]:
            bond = vecs[(i, other)]
            placed.append((point + bond / np.linalg.norm(bond) * 0.96) / 10.0)
        while len(placed) < 2:
            placed.append(point / 10.0 + np.array([0.096, 0.0, 0.0] if not placed else [0.0, 0.096, 0.0]))
        oxygen = (point / 10.0).tolist()
        waters.append(
            {
                "oxygen": oxygen,
                "hydrogen1": placed[0].tolist(),
                "hydrogen2": placed[1].tolist(),
            }
        )
    return json.dumps(waters)


def clavis_box(cell: np.ndarray) -> list[float]:
    """Nanometres. Six values when the cell is tilted, else three."""
    length = np.array([cell[0, 0], cell[1, 1], cell[2, 2]]) / 10.0
    tilt = np.array([cell[1, 0], cell[2, 0], cell[2, 1]]) / 10.0
    if np.any(np.abs(tilt) > 1e-8):
        return [float(v) for v in (*length, *tilt)]
    return [float(v) for v in length]


def face_expr(cage, rings) -> str:
    counts = Counter(len(rings[rid]) for rid in cage)
    return " ".join(f"{size}^{counts[size]}" for size in sorted(counts))


def cycless_census(graph, frac, quick: bool):
    from cycless.cycles import cycles_iter
    from cycless.polyhed import polyhedra_iter

    rings = [[int(node) for node in ring] for ring in cycles_iter(graph, 8, pos=frac)]
    cages = list(polyhedra_iter(rings, 22, quick=quick))
    return rings, Counter(face_expr(cage, rings) for cage in cages)


class Dseams:
    def __init__(self, pos: np.ndarray, cell: np.ndarray):
        lib_path = ROOT / "bbdir-release" / "bindings" / "libseams_lang.so"
        self.lib = ctypes.CDLL(str(lib_path))
        self.lib.seams_frame_from_arrays.restype = ctypes.c_void_p
        self.lib.seams_frame_from_arrays.argtypes = [
            ctypes.POINTER(ctypes.c_double),
            ctypes.c_int,
            ctypes.POINTER(ctypes.c_double),
            ctypes.POINTER(ctypes.c_int),
            ctypes.POINTER(ctypes.c_int),
        ]
        self.lib.seams_frame_free.argtypes = [ctypes.c_void_p]
        self.lib.seams_frame_cages.argtypes = [
            ctypes.c_void_p,
            ctypes.c_double,
            ctypes.c_int,
            ctypes.c_int,
            ctypes.c_char_p,
            ctypes.c_int,
            ctypes.c_void_p,
            ctypes.POINTER(ctypes.c_int),
            ctypes.POINTER(ctypes.c_int),
            ctypes.POINTER(ctypes.c_int),
            ctypes.POINTER(ctypes.c_int),
        ]
        self.lib.seams_frame_signature.restype = ctypes.c_int
        self.lib.seams_frame_signature.argtypes = [
            ctypes.c_void_p,
            ctypes.c_double,
            ctypes.c_int,
            ctypes.c_char_p,
            ctypes.c_char_p,
            ctypes.POINTER(ctypes.c_int),
            ctypes.POINTER(ctypes.c_int),
        ]
        n = len(pos)
        low, bound, tilt = lammps_box(cell)
        box = np.zeros(9, dtype=np.float64)
        box[0:3] = low
        box[3:6] = bound
        box[6:9] = tilt
        types = np.ones(n, dtype=np.int32)
        ids = np.arange(1, n + 1, dtype=np.int32)
        self._keep = (pos, box, types, ids)
        self.frame = self.lib.seams_frame_from_arrays(
            pos.ctypes.data_as(ctypes.POINTER(ctypes.c_double)),
            n,
            box.ctypes.data_as(ctypes.POINTER(ctypes.c_double)),
            types.ctypes.data_as(ctypes.POINTER(ctypes.c_int)),
            ids.ctypes.data_as(ctypes.POINTER(ctypes.c_int)),
        )
        if not self.frame:
            raise RuntimeError("the language library rejected the frame")

    def close(self) -> None:
        if self.frame:
            self.lib.seams_frame_free(self.frame)
            self.frame = None

    def seeded(self) -> tuple[int, int, int, int]:
        hexagonal = ctypes.c_int()
        cubic = ctypes.c_int()
        water = ctypes.c_int()
        both = ctypes.c_int()
        rc = self.lib.seams_frame_cages(
            self.frame,
            CUTOFF,
            1,
            4,
            b"seeded",
            0,
            None,
            hexagonal,
            cubic,
            water,
            both,
        )
        if rc != 0:
            raise RuntimeError("seams_frame_cages failed")
        return hexagonal.value, cubic.value, water.value, both.value

    def signature(self, spec: str) -> int:
        cages = ctypes.c_int()
        atoms = ctypes.c_int()
        rc = self.lib.seams_frame_signature(
            self.frame,
            CUTOFF,
            1,
            b"cutoff",
            spec.encode(),
            ctypes.byref(cages),
            ctypes.byref(atoms),
        )
        if rc != 0:
            raise RuntimeError(f"seams_frame_signature {spec} failed")
        return int(cages.value)


def cmd_signature(dump: Path) -> dict[str, int]:
    proc = subprocess.run(
        [
            str(ROOT / "bbdir-release" / "seams"),
            "cages",
            str(dump),
            "--type",
            "1",
            "--graph",
            "cutoff",
            "--signature",
            "512+51262+51264",
        ],
        cwd=ROOT,
        check=True,
        capture_output=True,
        text=True,
    )
    found = {}
    for token, count in SIG_LINE.findall(proc.stdout):
        key = SIG_TOKEN.get(token)
        if key:
            found[key] = int(count)
    if set(found) != set(SIGNATURES):
        raise RuntimeError("cmdCages did not report 512, 51262 and 51264")
    return found


def prepare(dump: Path):
    import networkx as nx

    pos, cell = read_dump(dump)
    edges, vecs = min_image_edges(pos, cell)
    degree = Counter()
    for i, j in edges:
        degree[i] += 1
        degree[j] += 1
    if any(degree[i] != 4 for i in range(len(pos))):
        raise RuntimeError("the oxygen graph is not 4-regular")
    frac = np.mod(pos @ np.linalg.inv(cell), 1.0)
    graph = nx.Graph()
    graph.add_nodes_from(range(len(pos)))
    graph.add_edges_from(edges)
    waters = clavis_waters(pos, cell, edges, vecs)
    box = clavis_box(cell)
    dseams = Dseams(pos, cell)
    return {
        "pos": pos,
        "cell": cell,
        "frac": frac,
        "graph": graph,
        "waters": waters,
        "box": box,
        "dseams": dseams,
        "dump": dump,
        "edges": len(edges),
    }


def time_call(fn) -> float:
    start = time.perf_counter()
    fn()
    return (time.perf_counter() - start) * 1e3


def worker(dump: Path, samples: int) -> None:
    import clavis_engine
    from cycless.cycles import cycles_iter
    from cycless.polyhed import polyhedra_iter

    bundle = prepare(dump)
    dseams = bundle["dseams"]

    def seeded():
        dseams.seeded()

    def command():
        cmd_signature(dump)

    def clavis():
        clavis_engine.analyze_waters(
            bundle["waters"], "[]", bundle["box"], "trace", 0.35, 35.0
        )

    def cycless():
        rings = [
            [int(node) for node in ring]
            for ring in cycles_iter(bundle["graph"], 8, pos=bundle["frac"])
        ]
        list(polyhedra_iter(rings, 22, quick=True))

    arms = (
        ("dseams_cages", seeded),
        ("dseams_cmd", command),
        ("clavis", clavis),
        ("cycless", cycless),
    )
    for _name, fn in arms:
        fn()
    taken = {name: [] for name, _fn in arms}
    for sample in range(samples):
        for name, fn in arms:
            taken[name].append(time_call(fn))
            print(f"sample {sample} {name} {taken[name][-1]:.3f}", file=sys.stderr, flush=True)
    dseams.close()
    json.dump(taken, sys.stdout)
    sys.stdout.write("\n")


def census_of(dump: Path, winding: bool) -> dict:
    import clavis_engine
    import networkx as nx
    from cycless.cycles import cycles_iter

    pos, cell = read_dump(dump)
    edges, vecs = min_image_edges(pos, cell)
    frac = np.mod(pos @ np.linalg.inv(cell), 1.0)
    graph = nx.Graph()
    graph.add_nodes_from(range(len(pos)))
    graph.add_edges_from(edges)
    dseams = Dseams(pos, cell)
    signatures = {spec: dseams.signature(spec) for spec in (*SIGNATURES, "51268", "sh")}
    seeded = dseams.seeded()
    command = cmd_signature(dump)
    dseams.close()
    # The timed arm uses quick=True. On a unit cell that result is checked
    # against the full polyhedra search and against the ring list that keeps
    # spanning cycles. A replica keeps the quick census.
    rings_pos, cages_pos = cycless_census(graph, frac, quick=True)
    if winding:
        _rings_full, cages_full = cycless_census(graph, frac, quick=False)
        quick_same = cages_full == cages_pos
        rings_all = list(cycles_iter(graph, 8))
    else:
        quick_same = True
        rings_all = []
    waters = clavis_waters(pos, cell, edges, vecs)
    raw = json.loads(
        clavis_engine.analyze_waters(waters, "[]", clavis_box(cell), "trace", 0.35, 35.0)
    )
    clavis = {
        "512": len(raw["cages_512"]),
        "51262": len(raw["cages_51262"]),
        "51264": len(raw["cages_51264"]),
    }
    return {
        "molecules": len(pos),
        "edges": len(edges),
        "dseams": signatures,
        "dseams_cmd": command,
        "seeded": {
            "hex": seeded[0],
            "cubic": seeded[1],
            "water": seeded[2],
            "both": seeded[3],
        },
        "cycless": dict(cages_pos),
        "cycless_quick_same": quick_same,
        "rings_pos": len(rings_pos),
        "rings_pos_sizes": dict(Counter(len(ring) for ring in rings_pos)),
        "rings_all": len(rings_all),
        "rings_all_sizes": dict(Counter(len(ring) for ring in rings_all)),
        "clavis": clavis,
    }


def fmt_counts(mapping: dict) -> str:
    return " ".join(f"{key}={mapping[key]}" for key in mapping)


def package_version(name: str) -> str:
    try:
        from importlib.metadata import version
    except ImportError:
        return ""
    try:
        return version(name)
    except Exception:
        return ""


def write_report(path: Path, rows: list[dict]) -> None:
    cores, affinity = paper_stats.cores_now()
    lines = [
        f"cores {cores}",
        f"model {cpu_model()}",
        f"affinity {affinity}",
        "n 8",
        "genice2 " + package_version("genice2"),
        "cycless " + package_version("cycless"),
        "clavis_engine " + package_version("clavis-engine"),
        "cutoff_A 3.5",
        "clavis trace 0.35 nm 35 deg ice-rule hydrogens",
        "cycless cycles_iter maxsize 8 fractional positions, polyhedra_iter maxnfaces 22 quick",
        "dseams franzblau rings include winding shortest paths; cage faces are the contractible ones",
        "winding rings are compared on the unit cell; a replica says cycless_all unit_cell",
        "clavis and cycless have no thread pool; the process is pinned and RAYON_NUM_THREADS is set",
        "dseams keeps one thread below 2048 molecules",
        "each arm is interleaved inside one width",
    ]
    for row in rows:
        cell = row["cell"]
        rep = " ".join(str(v) for v in row["rep"])
        counted = row["census"]
        lines.append(
            f"cell {cell} rep {rep} molecules {counted['molecules']} edges {counted['edges']}"
        )
        lines.append(
            "counts dseams "
            + " ".join(f"{spec}={counted['dseams'][spec]}" for spec in (*SIGNATURES, "51268", "sh"))
        )
        lines.append("counts dseams_cmd " + fmt_counts(counted["dseams_cmd"]))
        lines.append(
            "counts seeded hex={hex} cubic={cubic} water={water} both={both}".format(
                **counted["seeded"]
            )
        )
        lines.append("counts clavis " + fmt_counts(counted["clavis"]))
        lines.append("counts cycless " + fmt_counts(counted["cycless"]))
        if counted["rings_all"]:
            spanned = f"cycless_all {counted['rings_all']} {counted['rings_all_sizes']}"
        else:
            spanned = "cycless_all unit_cell"
        lines.append(
            f"rings cycless_pos {counted['rings_pos']} {counted['rings_pos_sizes']} "
            f"{spanned} quick_same {int(counted['cycless_quick_same'])}"
        )
        for width, taken in row["times"]:
            lines.append(f"threads {width}")
            for arm, samples in taken.items():
                summary = paper_stats.summarize(samples)
                lines.append(
                    f"arm {arm} median_ms {summary['median_ms']:.3f} "
                    f"iqr_ms {summary['iqr_ms']:.3f} n {summary['n']}"
                )
                lines.append(f"samples {arm} {paper_stats.fmt_samples(samples)}")
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text("\n".join(lines) + "\n", encoding="utf-8")


def sitting(out: Path, samples: int) -> None:
    frames = ROOT / "benchmarks" / "genice-frames"
    frames.mkdir(parents=True, exist_ok=True)
    rows = []
    for name, rep in CELLS:
        dump = frames / f"{name}-{rep[0]}x{rep[1]}x{rep[2]}.lammpstrj"
        print(f"generate {name} {rep}", file=sys.stderr, flush=True)
        pos, cell = generate(name, rep)
        write_dump(dump, pos, cell)
        print(f"census {name} {rep} n {len(pos)}", file=sys.stderr, flush=True)
        counted = census_of(dump, winding=rep == (1, 1, 1))
        if not counted["cycless_quick_same"]:
            raise RuntimeError(f"{name} {rep}: quick polyhedra differ")
        times = []
        for width in (1, 8):
            print(f"time {name} {rep} threads {width}", file=sys.stderr, flush=True)
            proc = subprocess.run(
                [
                    "taskset",
                    "-c",
                    paper_stats.taskset_cores(width),
                    sys.executable,
                    str(Path(__file__).resolve()),
                    "--worker",
                    str(dump.relative_to(ROOT)),
                    str(samples),
                ],
                cwd=ROOT,
                env=paper_stats.pin_env(width),
                check=True,
                stdout=subprocess.PIPE,
                text=True,
            )
            taken = json.loads(proc.stdout)
            times.append((width, taken))
        rows.append({"cell": name, "rep": rep, "census": counted, "times": times})
        write_report(out, rows)


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--out", type=Path, default=Path("benchmarks/cage-h2h-v8.txt"))
    parser.add_argument("--samples", type=int, default=8)
    parser.add_argument("--worker", nargs=2)
    args = parser.parse_args()
    if args.worker:
        dump = ROOT / args.worker[0]
        worker(dump, int(args.worker[1]))
        return
    if args.samples < 8:
        raise SystemExit("the sitting uses at least 8 samples")
    sitting(args.out, args.samples)


if __name__ == "__main__":
    main()
