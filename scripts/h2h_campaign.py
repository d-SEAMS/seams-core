#!/usr/bin/env python3
# Copyright (c) d-SEAMS developers
"""Head-to-head campaign: neighbour lists and ring statistics, cold calls.

Three workloads run on oxygen-only periodic frames of ice and water. Every
tool gets the same frame, the same cutoff and the same question.

- cutoff: every pair closer than 3.5 A.
- knn: the four nearest neighbours of every atom, candidates inside 5.0 A.
- rings: Franzblau primitive rings of size 3 to 6 on the 3.5 A graph.

One sample is one process. The process is pinned to one core with every
thread-count variable set to one. It imports the tool, loads the frame,
builds the tool's own input object, and makes a single timed call. There
is no warm-up call and no repeated call, so the first call of a fresh
process is what is timed. The timed call builds the neighbour structure
from scratch and ends with the output as arrays or a ring histogram. The
loaded frame and the tool's input object are built before the timer.

Rounds interleave the arms of each cell. Odd rounds run the cells and the
arms in reverse order. Every cell carries a second d-SEAMS arm, seams_aa,
that is the same binary run through the same code: the A/A band is the
spread between the two.

The record is JSON Lines and append-only. Every row carries the host id,
CPU model, lscpu summary, kernel and code revision of the run.

    h2h_campaign.py run --head-to-head neighbours --record R.jsonl
    h2h_campaign.py summary --record R.jsonl
"""

from __future__ import annotations

import argparse
import ast
import ctypes
import dataclasses
import datetime
import gc
import hashlib
import importlib.metadata
import importlib.util
import itertools
import json
import math
import os
import platform
import re
import resource
import signal
import subprocess
import sys
import time
import uuid
import warnings
from pathlib import Path

import numpy as np

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "scripts"))

import h2h_protocol as proto  # noqa: E402

SCHEMA = "seams-h2h-campaign/1"
CUTOFF_A = proto.CUTOFF_A
KNN = 4
KNN_REACH_A = CUTOFF_A + 1.5
MAX_RING = proto.MAX_RING
MIN_ROUNDS = proto.REPEATS
ROUNDS = 8
TIMEOUT_S = 900

HEAD_TO_HEADS = {
    "neighbours": ("cutoff", "knn"),
    "rings": ("rings",),
}

# Frames. Ice is the orthorhombic Ih cell of h2h_freud_ovito.ice_positions
# with ncell^3 * 8 molecules. The mixed frames are scripts/mixed_phase_frames.
# The mW frames are tracked mW water dumps. sII is the GenIce hydrate cell.
ICE_CELLS = {"ice512": 4, "ice4096": 8, "ice32768": 16, "ice110592": 24}
MIXED = ("interface", "stacking", "melted")
DUMPS = {
    "mw_interface": "input/traj/interface_f0.lammpstrj",
    "mw_nucleation": "input/traj/nucleation_f1.lammpstrj",
    "sII": "input/traj/genice_sII.lammpstrj",
}
SYSTEMS = (*ICE_CELLS, *MIXED, *DUMPS)

TOOL_OF = {
    "seams": "seams",
    "seams_aa": "seams",
    "freud": "freud",
    "ovito": "ovito",
    "vesin": "vesin",
    "mdapy": "mdapy",
    "pyscal3": "pyscal3",
    "rsm": "rings-statistics-matter",
    "sqq": "sqq",
    "cycless": "cycless",
    "fastatomstruct": "fastatomstruct",
}
DIST_OF = {
    "freud": "freud-analysis",
    "ovito": "ovito",
    "vesin": "vesin",
    "mdapy": "mdapy",
    "pyscal3": "pyscal3",
    "rings-statistics-matter": "rings-statistics-matter",
    "sqq": "sqq",
    "cycless": "cycless",
    "fastatomstruct": "fastatomstruct",
}
ARMS = {
    "cutoff": ("seams", "seams_aa", "freud", "ovito", "vesin", "mdapy", "pyscal3"),
    "knn": ("seams", "seams_aa", "freud", "ovito", "mdapy", "pyscal3"),
    "rings": ("seams", "seams_aa", "ovito", "rsm", "mdapy", "sqq", "fastatomstruct", "cycless"),
}
# Largest frame an arm is attempted on. A pure-Python census that needs
# minutes on a few hundred molecules is not repeated on thousands.
LIMITS = {("cycless", "rings"): 1000}
# Rounds an arm runs on one cell when a single sample costs about ten
# minutes. Six is the floor.
ROUND_CAPS = {("fastatomstruct", "rings", "ice110592"): MIN_ROUNDS}
# What the timed call covers when it is not the same for every arm.
ARM_NOTES = {
    "cycless/rings": "takes a graph, not positions: the 3.5 A graph is built before the timer",
    "ovito/rings": "CreateBonds then FindRings in one pipeline; the table has no 3-rings",
    "sqq/rings": "reports no 3-rings",
    "seams/knn": "undirected union of the 4-nearest-neighbour relations, rows may exceed 4",
    "freud/knn": "strict per-atom 4-nearest lists, symmetrised outside the timer to compare",
    "ovito/knn": "strict per-atom 4-nearest lists, symmetrised outside the timer to compare",
    "mdapy/knn": "strict per-atom 4-nearest lists, symmetrised outside the timer to compare",
    "pyscal3/knn": "strict per-atom 4-nearest lists, symmetrised outside the timer to compare",
}
# Ring sizes a tool reports. OVITO and SQQ do not report 3-rings.
RING_SIZES = {"ovito": (4, 5, 6), "sqq": (4, 5, 6)}
REFERENCE = "seams"
TWIN = "seams_aa"
ORACLE_LIMIT = 4096

PINNED_ENV = {
    "OMP_NUM_THREADS": "1",
    "OMP_DYNAMIC": "FALSE",
    "OMP_WAIT_POLICY": "PASSIVE",
    "MKL_NUM_THREADS": "1",
    "OPENBLAS_NUM_THREADS": "1",
    "NUMEXPR_NUM_THREADS": "1",
    "RAYON_NUM_THREADS": "1",
    "OVITO_THREAD_COUNT": "1",
    "QT_QPA_PLATFORM": "offscreen",
    "FASTATOMSTRUCT_NO_TELEMETRY": "1",
    "FASTATOMSTRUCT_SUPPRESS_TELEMETRY_NOTICE": "1",
    "PYTHONHASHSEED": "0",
    "PYTHONWARNINGS": "ignore",
}
CARRIED_ENV = ("PATH", "HOME", "LANG", "LC_ALL", "VIRTUAL_ENV", "SEAMS_LIB", "TMPDIR")
OVITO_LIBS_ENV = "H2H_OVITO_LD_LIBRARY_PATH"


# ---------------------------------------------------------------- frames


def read_dump(path: Path):
    """First frame of an orthorhombic LAMMPS dump, moved to a zero origin."""
    lines = Path(path).read_text(encoding="utf-8").splitlines()
    count = int(lines[3])
    if "xy" in lines[4]:
        raise ValueError(f"{path.name} is triclinic")
    lows, highs = [], []
    for line in lines[5:8]:
        lo, hi = line.split()[:2]
        lows.append(float(lo))
        highs.append(float(hi))
    columns = lines[8].split()[2:]
    where = [columns.index(name) for name in ("x", "y", "z")]
    types = columns.index("type") if "type" in columns else None
    rows = [line.split() for line in lines[9 : 9 + count]]
    if types is not None:
        rows = [row for row in rows if int(row[types]) == 1]
    pos = np.array([[float(row[w]) for w in where] for row in rows], dtype=np.float64)
    lengths = np.array(highs) - np.array(lows)
    pos = np.mod(pos - np.array(lows), lengths)
    return np.ascontiguousarray(pos), tuple(float(v) for v in lengths)


def build_frame(name: str):
    import h2h_freud_ovito as base

    if name in ICE_CELLS:
        pos, lengths = base.ice_positions(ICE_CELLS[name])
    elif name in MIXED:
        import mixed_phase_frames

        pos, lengths = mixed_phase_frames.FRAMES[name]()
    elif name in DUMPS:
        pos, lengths = read_dump(ROOT / DUMPS[name])
    else:
        raise KeyError(name)
    pos = np.ascontiguousarray(pos, dtype=np.float64)
    return pos, tuple(float(v) for v in lengths)


def frame_digest(pos, lengths) -> str:
    held = hashlib.blake2b(digest_size=16)
    held.update(np.ascontiguousarray(pos, dtype=np.float64).tobytes())
    held.update(np.asarray(lengths, dtype=np.float64).tobytes())
    return held.hexdigest()


def write_frames(directory: Path, names) -> dict:
    directory.mkdir(parents=True, exist_ok=True)
    digests = {}
    for name in names:
        pos, lengths = build_frame(name)
        np.savez(directory / f"{name}.npz", positions=pos, lengths=np.asarray(lengths))
        digests[name] = {"molecules": int(pos.shape[0]), "digest": frame_digest(pos, lengths)}
    return digests


def load_frame(directory: Path, name: str):
    with np.load(Path(directory) / f"{name}.npz") as data:
        return np.ascontiguousarray(data["positions"]), tuple(float(v) for v in data["lengths"])


# --------------------------------------------------------------- outputs


def digest_of(array: np.ndarray) -> str:
    return hashlib.blake2b(np.ascontiguousarray(array).tobytes(), digest_size=16).hexdigest()


def pair_keys(left, right, n: int) -> np.ndarray:
    """Sorted unique i*n+j keys of the undirected pairs, i < j."""
    a = np.minimum(left, right).astype(np.int64, copy=False)
    b = np.maximum(left, right).astype(np.int64, copy=False)
    keep = a != b
    return np.unique(a[keep] * np.int64(n) + b[keep])


def degrees(keys: np.ndarray, n: int) -> np.ndarray:
    return np.bincount(np.concatenate((keys // n, keys % n)), minlength=n)


def canonical(workload: str, raw, n: int, tool: str):
    """Reduce a tool's raw output to a common form. Not timed."""
    kind = raw[0]
    if kind == "hist":
        hist = {int(k): int(v) for k, v in raw[1].items() if 3 <= int(k) <= MAX_RING}
        for size in RING_SIZES.get(tool, range(3, MAX_RING + 1)):
            hist.setdefault(size, 0)
        sizes = sorted(s for s in hist if s in RING_SIZES.get(tool, range(3, MAX_RING + 1)))
        return {
            "form": "histogram",
            "hist": {str(s): hist[s] for s in sizes},
            "total": int(sum(hist[s] for s in sizes)),
        }, None
    if kind == "csr":
        off = np.asarray(raw[1], dtype=np.int64)
        col = np.asarray(raw[2], dtype=np.int64)
        left = np.repeat(np.arange(off.shape[0] - 1, dtype=np.int64), np.diff(off))
        right = col[: left.shape[0]]
        row_len = np.diff(off)
        form = "csr"
    elif kind == "rows":
        idx = np.asarray(raw[1], dtype=np.int64)
        left = np.repeat(np.arange(idx.shape[0], dtype=np.int64), idx.shape[1])
        right = idx.reshape(-1)
        row_len = np.full(idx.shape[0], idx.shape[1], dtype=np.int64)
        form = "rows"
    else:
        left = np.asarray(raw[1], dtype=np.int64)
        right = np.asarray(raw[2], dtype=np.int64)
        row_len = np.bincount(left, minlength=n)
        form = "pairs"
    keys = pair_keys(left, right, n)
    deg = degrees(keys, n)
    out = {
        "form": form,
        "directed_entries": int(left.shape[0]),
        "n_pairs": int(keys.shape[0]),
        "pairs_digest": digest_of(keys),
        "degree_digest": digest_of(deg),
        "degree_min": int(deg.min()) if n else 0,
        "degree_max": int(deg.max()) if n else 0,
        "row_len_hist": {
            str(k): int(v) for k, v in enumerate(np.bincount(row_len)) if v
        },
    }
    return out, keys


# --------------------------------------------------------------- adapters


class Adapter:
    """prepare() builds the tool's input. call() is the one timed call."""

    def prepare(self, pos, lengths, workload):
        raise NotImplementedError

    def call(self, workload):
        raise NotImplementedError


class SeamsAdapter(Adapter):
    def prepare(self, pos, lengths, workload):
        import h2h_freud_ovito as base
        import h2h_rings
        import release_build

        self.base = base
        path = base.library_path()
        self.build = release_build.ensure_release(str(path))[1]
        lib = base.load_seams(path)
        h2h_rings.bind_rings(lib)
        base.require(lib, lib.seams_set_num_threads(1), "threads")
        self.lib = lib
        self.system = base.frame_system(lib, pos, lengths)

    def call(self, workload):
        lib, base, system = self.lib, self.base, self.system
        if workload == "rings":
            rings = ctypes.c_void_p()
            base.require(
                lib,
                lib.seams_rings_find(system, CUTOFF_A, 1, 0, MAX_RING, ctypes.byref(rings)),
                "rings",
            )
            hist = (ctypes.c_int * (MAX_RING + 1))()
            base.require(
                lib, lib.seams_rings_histogram(rings, hist, MAX_RING + 1), "histogram"
            )
            lib.seams_rings_destroy(rings)
            return ("hist", {k: int(hist[k]) for k in range(3, MAX_RING + 1)})
        found = ctypes.c_void_p()
        if workload == "cutoff":
            status = lib.seams_neighbours_build(system, CUTOFF_A, 1, 0, ctypes.byref(found))
        else:
            status = lib.seams_neighbours_knn(system, KNN, CUTOFF_A, 1, 0, ctypes.byref(found))
        base.require(lib, status, workload)
        off, col = base.copy_csr(lib, found)
        return ("csr", off, col)


class FreudAdapter(Adapter):
    def prepare(self, pos, lengths, workload):
        import freud

        freud.parallel.set_num_threads(1)
        self.freud = freud
        self.box = freud.box.Box(*lengths)
        self.pos = pos

    def call(self, workload):
        query = self.freud.locality.AABBQuery(self.box, self.pos)
        if workload == "cutoff":
            spec = {"r_max": CUTOFF_A, "exclude_ii": True}
        else:
            spec = {
                "mode": "nearest",
                "num_neighbors": KNN,
                "r_max": KNN_REACH_A,
                "exclude_ii": True,
            }
        found = query.query(self.pos, spec).toNeighborList()
        return ("pairs", np.asarray(found.query_point_indices), np.asarray(found.point_indices))


class OvitoAdapter(Adapter):
    def prepare(self, pos, lengths, workload):
        import ovito.modifiers as modifiers
        from ovito.data import DataCollection
        from ovito.pipeline import Pipeline, StaticSource

        data = DataCollection()
        data.create_cell(np.diag(lengths), pbc=(True, True, True))
        particles = data.create_particles(count=int(pos.shape[0]))
        particles.create_property("Position", data=pos)
        self.data = data
        if workload == "knn":
            return
        self.pipe = Pipeline(source=StaticSource(data=data))
        self.pipe.modifiers.append(modifiers.CreateBondsModifier(cutoff=CUTOFF_A))
        if workload == "rings":
            self.pipe.modifiers.append(modifiers.FindRingsModifier(maximum_ring_size=MAX_RING))

    def call(self, workload):
        if workload == "knn":
            from ovito.data import NearestNeighborFinder

            idx, _vectors = NearestNeighborFinder(KNN, self.data).find_all()
            return ("rows", np.asarray(idx))
        out = self.pipe.compute()
        if workload == "cutoff":
            topology = np.asarray(out.particles.bonds.topology)
            return ("pairs", topology[:, 0], topology[:, 1])
        table = out.tables["ring-size-histogram"]
        sizes = np.asarray(table["Ring Size"][...])
        counts = np.asarray(table["Count"][...])
        return ("hist", {int(s): int(c) for s, c in zip(sizes, counts)})


class VesinAdapter(Adapter):
    def prepare(self, pos, lengths, workload):
        import vesin

        self.calc = vesin.NeighborList(cutoff=CUTOFF_A, full_list=True, n_threads=1)
        self.box = np.diag(np.asarray(lengths, dtype=np.float64))
        self.pos = pos

    def call(self, workload):
        left, right = self.calc.compute(
            points=self.pos, box=self.box, periodic=True, quantities="ij", copy=False
        )
        return ("pairs", left, right)


class MdapyAdapter(Adapter):
    def prepare(self, pos, lengths, workload):
        import mdapy

        mdapy.set_num_threads(1)
        if mdapy.get_num_threads() != 1:
            raise RuntimeError("mdapy did not accept one thread")
        self.mdapy = mdapy
        self.system = mdapy.System(pos, box=lengths)

    def call(self, workload):
        if workload == "rings":
            from mdapy.analysis import RingStatistics

            found = RingStatistics(max_size=MAX_RING, min_size=3, cutoff=CUTOFF_A).compute(
                self.system
            )
            return ("hist", {k: int(found.counts[k]) for k in range(3, MAX_RING + 1)})
        if workload == "cutoff":
            neigh = self.mdapy.NeighborList.build(self.system, cutoff=CUTOFF_A)
        else:
            neigh = self.mdapy.NeighborList.build(self.system, num_neighbors=KNN)
        return ("csr", np.asarray(neigh.offsets), np.asarray(neigh.j))


class Pyscal3Adapter(Adapter):
    def prepare(self, pos, lengths, workload):
        import pyscal3
        from ase import Atoms
        from pyscal3._bridge import neighbor_arrays

        pyscal3.set_num_threads(1)
        self.pyscal3 = pyscal3
        self.arrays = neighbor_arrays
        self.atoms = Atoms(positions=pos, cell=list(lengths), pbc=True)

    def call(self, workload):
        if workload == "cutoff":
            self.pyscal3.find_neighbors(
                self.atoms, method="cutoff", cutoff=CUTOFF_A, store_rows=False
            )
        else:
            self.pyscal3.find_neighbors(
                self.atoms, method="number", nmax=KNN, store_rows=False
            )
        offsets, values = self.arrays(self.atoms, "neighbors")
        return ("csr", np.asarray(offsets), np.asarray(values["neighbors"]))


class RingsStatisticsAdapter(Adapter):
    def prepare(self, pos, lengths, workload):
        import rings_statistics_matter
        from ase import Atoms

        self.module = rings_statistics_matter
        self.atoms = Atoms(
            symbols=["O"] * int(pos.shape[0]), positions=pos, cell=list(lengths), pbc=True
        )

    def call(self, workload):
        _counts, blocks = self.module.ring_statistics(
            self.atoms, cutoff=CUTOFF_A, maxlvl=MAX_RING + 2, supercell=False,
            deduplicate=True, backend="rust", num_threads=1,
        )
        hist = {}
        for block in blocks:
            arr = np.asarray(block)
            if arr.size and arr.ndim == 2:
                hist[int(arr.shape[1])] = int(arr.shape[0])
        return ("hist", hist)


class SqqAdapter(Adapter):
    def prepare(self, pos, lengths, workload):
        from sqq.api import analyze_frame, load_config
        from sqq.models import Atom, Frame

        nm = pos / 10.0
        atoms = [
            Atom(index=i, resid=i + 1, resname="SOL", atomname="O", atomid=i + 1,
                 xyz=nm[i], molecule_id=i + 1)
            for i in range(nm.shape[0])
        ]
        self.frame = Frame(name="h2h", atoms=atoms, box=np.asarray(lengths) / 10.0)
        self.analyze = analyze_frame
        self.config = load_config({
            "engine": "py",
            "graph": {"mode": "oo", "oo_cutoff_nm": CUTOFF_A / 10.0},
            "ring": {"sizes": [4, 5, 6], "definition": "shortest_path", "chordless": True},
            "half_cage": {"enabled": False},
            "quasi_cage": {"enabled": False},
            "cage": {"enabled": False},
            "order_parameter": {"enabled": []},
            "parallel": {"workers": 1, "math_threads": 1},
        }, engine="py")

    def call(self, workload):
        result = self.analyze(self.frame, self.config, engine="py")
        return ("hist", {size: len(result.rings.get(size, [])) for size in (4, 5, 6)})


class FastatomstructAdapter(Adapter):
    def prepare(self, pos, lengths, workload):
        import fastatomstruct as fs

        fs.disable_telemetry()
        fs.suppress_telemetry_notice()
        fs.set_num_threads(1)
        if fs.get_num_threads() != 1:
            raise RuntimeError("fastatomstruct did not accept one thread")
        self.fs = fs
        cell = np.diag(np.asarray(lengths, dtype=np.float64))
        self.atoms = fs.Atoms(pos, cell, [fs.Element.O] * int(pos.shape[0]))
        self.cutoff = fs.analysis.distance.CutoffMode.Fixed(CUTOFF_A)

    def call(self, workload):
        sizes, counts = self.fs.analysis.rings.ring_size_distribution(
            self.atoms, self.cutoff, MAX_RING, self.fs.analysis.rings.RingMethod.King
        )
        return ("hist", {int(s): int(round(float(c))) for s, c in zip(sizes, counts)})


class CyclessAdapter(Adapter):
    """The bond graph is built before the timer, as in the one-pair sessions."""

    def prepare(self, pos, lengths, workload):
        import freud
        import networkx as nx
        from cycless.cycles import cycles_iter

        freud.parallel.set_num_threads(1)
        found = freud.locality.AABBQuery(freud.box.Box(*lengths), pos).query(
            pos, {"r_max": CUTOFF_A, "exclude_ii": True}
        ).toNeighborList()
        graph = nx.Graph()
        graph.add_nodes_from(range(int(pos.shape[0])))
        graph.add_edges_from(
            (int(i), int(j))
            for i, j in zip(found.query_point_indices, found.point_indices)
            if int(i) < int(j)
        )
        self.graph = graph
        self.cycles = cycles_iter

    def call(self, workload):
        hist = {k: 0 for k in range(3, MAX_RING + 1)}
        for cycle in self.cycles(self.graph, MAX_RING):
            if len(cycle) <= MAX_RING:
                hist[len(cycle)] += 1
        return ("hist", hist)


ADAPTERS = {
    "seams": SeamsAdapter,
    "freud": FreudAdapter,
    "ovito": OvitoAdapter,
    "vesin": VesinAdapter,
    "mdapy": MdapyAdapter,
    "pyscal3": Pyscal3Adapter,
    "rings-statistics-matter": RingsStatisticsAdapter,
    "sqq": SqqAdapter,
    "fastatomstruct": FastatomstructAdapter,
    "cycless": CyclessAdapter,
}


def adapter_for(arm: str) -> Adapter:
    return ADAPTERS[TOOL_OF[arm]]()


def applicable(arm: str, workload: str, molecules: int | None = None) -> bool:
    if arm not in ARMS[workload]:
        return False
    limit = LIMITS.get((arm, workload))
    return not (limit is not None and molecules is not None and molecules > limit)


def tool_version(tool: str, build: str | None = None) -> str:
    if tool == "seams":
        return f"{head_revision()[:8]} {build or ''}".strip()
    try:
        return importlib.metadata.version(DIST_OF[tool])
    except importlib.metadata.PackageNotFoundError:
        return "unknown"


# ------------------------------------------------------------ measurement


def proc_status(key: str) -> int:
    with open("/proc/self/status", encoding="ascii") as handle:
        for line in handle:
            if line.startswith(key + ":"):
                return int(line.split()[1])
    raise KeyError(key)


def reset_peak() -> bool:
    try:
        with open("/proc/self/clear_refs", "w", encoding="ascii") as handle:
            handle.write("5")
    except OSError:
        return False
    return True


def cpus_allowed() -> str:
    return proto.cpus_allowed()


def run_worker(args) -> int:
    warnings.simplefilter("ignore")
    frames = Path(args.frames)
    pos, lengths = load_frame(frames, args.system)
    n = int(pos.shape[0])
    adapter = adapter_for(args.arm)
    adapter.prepare(pos, lengths, args.workload)
    tool = TOOL_OF[args.arm]

    gc.collect()
    gc.disable()
    peak_reset = reset_peak()
    rss_before = proc_status("VmRSS")
    t0 = time.perf_counter_ns()
    c0 = time.process_time_ns()
    raw = adapter.call(args.workload)
    c1 = time.process_time_ns()
    t1 = time.perf_counter_ns()
    peak = proc_status("VmHWM")
    rss_after = proc_status("VmRSS")
    threads = proc_status("Threads")
    gc.enable()

    summary, keys = canonical(args.workload, raw, n, tool)
    if args.canon and keys is not None:
        out = Path(args.canon)
        out.mkdir(parents=True, exist_ok=True)
        np.save(out / f"{args.arm}.{args.workload}.{args.system}.npy", keys)
    row = {
        "arm": args.arm,
        "tool": tool,
        "tool_version": tool_version(tool, getattr(adapter, "build", None)),
        "workload": args.workload,
        "system": args.system,
        "molecules": n,
        "frame_digest": frame_digest(pos, lengths),
        "wall_ms": (t1 - t0) / 1e6,
        "cpu_ms": (c1 - c0) / 1e6,
        "rss_before_kb": rss_before,
        "rss_peak_kb": peak,
        "rss_after_kb": rss_after,
        "peak_resettable": peak_reset,
        "maxrss_kb": int(resource.getrusage(resource.RUSAGE_SELF).ru_maxrss),
        "os_threads": threads,
        "cpus_allowed": cpus_allowed(),
        "output": summary,
    }
    sys.stdout.write(json.dumps(row, sort_keys=True, allow_nan=False) + "\n")
    sys.stdout.flush()
    return 0


# ------------------------------------------------------------- provenance

PATH_TOKEN = re.compile(r"(?<![\w<>])/(?:[\w.+@=-]+/)*[\w.+@=-]*")


def scrub(text: str) -> str:
    """Replace absolute paths so that a record carries no machine paths."""
    text = text.replace(str(ROOT), "<repo>")
    return PATH_TOKEN.sub("<path>", text)


def git(*args: str) -> str:
    done = subprocess.run(
        ["git", *args], cwd=ROOT, capture_output=True, text=True, check=False,
    )
    return done.stdout.strip()


def head_revision() -> str:
    return os.environ.get("H2H_CODE_REVISION") or git("rev-parse", "HEAD")


ENGINE_PATHS = ("src", "include", "bindings", "meson.build", "meson_options.txt",
                "subprojects/packagefiles", "subprojects/linkcell.wrap",
                "subprojects/minimage.wrap", "subprojects/readcon-core.wrap")


def engine_dirty() -> bool:
    return bool(git("status", "--porcelain", "--untracked-files=no", "--", *ENGINE_PATHS))


def harness_blob() -> str:
    return git("hash-object", "scripts/h2h_campaign.py")


def lscpu_lines() -> list[str]:
    done = subprocess.run(
        ["lscpu"], capture_output=True, text=True, check=False, env={**os.environ, "LC_ALL": "C"},
    )
    return [line.rstrip() for line in done.stdout.splitlines() if line.strip()]


def lscpu_fields(lines: list[str]) -> dict:
    fields = {}
    for line in lines:
        key, _, value = line.partition(":")
        fields[key.strip()] = value.strip()
    return fields


def machine_summary(lines: list[str] | None = None) -> dict:
    lines = lines if lines is not None else lscpu_lines()
    f = lscpu_fields(lines)
    model = f.get("Model name", "unknown")
    short = (
        f"{model}; family {f.get('CPU family', '?')} model {f.get('Model', '?')} "
        f"stepping {f.get('Stepping', '?')}; {f.get('CPU(s)', '?')} cpus; "
        f"{f.get('Thread(s) per core', '?')} thread per core; "
        f"L2 {f.get('L2 cache', '?')}; L3 {f.get('L3 cache', '?')}; "
        f"hypervisor {f.get('Hypervisor vendor', 'none')}"
    )
    digest = hashlib.sha256("\n".join(lines).encode("utf-8")).hexdigest()
    return {"cpu_model": model, "lscpu": short, "lscpu_sha256": digest}


def host_id() -> str:
    """An opaque id for this machine: the host name is not recorded."""
    for source in ("/etc/machine-id", "/var/lib/dbus/machine-id"):
        try:
            seed = Path(source).read_text(encoding="ascii").strip()
            break
        except OSError:
            continue
    else:
        seed = platform.node()
    return "bench-" + hashlib.sha256(seed.encode("utf-8")).hexdigest()[:12]


def provenance(run_id: str, stamp: str) -> dict:
    summary = machine_summary()
    return {
        "schema": SCHEMA,
        "run_id": run_id,
        "date": stamp,
        "host": host_id(),
        "cpu_model": summary["cpu_model"],
        "lscpu": summary["lscpu"],
        "lscpu_sha256": summary["lscpu_sha256"],
        "kernel": platform.release(),
        "code_revision": head_revision(),
        "code_dirty": engine_dirty(),
        "harness_blob": harness_blob(),
    }


def read_text(path: str) -> str:
    try:
        return Path(path).read_text(encoding="ascii").strip()
    except OSError:
        return "unavailable"


STACK = ("linkcell", "minimage", "readcon-core")


def stack_revisions() -> dict:
    """The pinned and the checked-out revision of each canonical subproject."""
    out = {}
    for name in STACK:
        pinned = ""
        wrap = ROOT / "subprojects" / f"{name}.wrap"
        if wrap.is_file():
            for line in wrap.read_text(encoding="utf-8").splitlines():
                key, _, value = line.partition("=")
                if key.strip() == "revision":
                    pinned = value.strip()
        checkout = ROOT / "subprojects" / name
        head = ""
        if (checkout / ".git").exists():
            head = subprocess.run(
                ["git", "-C", str(checkout), "rev-parse", "HEAD"],
                capture_output=True, text=True, check=False).stdout.strip()
        out[name] = {"wrap_revision": pinned, "checkout": head}
    return out


def machine_row() -> dict:
    lines = lscpu_lines()
    return {
        "type": "machine",
        "lscpu_lines": lines,
        "kernel_version": platform.version(),
        "machine": platform.machine(),
        "python": platform.python_version(),
        "numpy": np.__version__,
        "memory_total_kb": int(
            next((l.split()[1] for l in read_text("/proc/meminfo").splitlines()
                  if l.startswith("MemTotal")), 0)
        ),
        "transparent_hugepage": read_text("/sys/kernel/mm/transparent_hugepage/enabled"),
        "cpufreq_governor": read_text("/sys/devices/system/cpu/cpu0/cpufreq/scaling_governor"),
        "stack": stack_revisions(),
        "tools": {
            tool: tool_version(tool) for tool in sorted(DIST_OF)
        },
    }


class Record:
    """Append-only JSON Lines. Rows are never rewritten or removed."""

    def __init__(self, path, prov: dict):
        self.path = Path(path)
        self.prov = prov
        self.path.parent.mkdir(parents=True, exist_ok=True)

    def append(self, row: dict) -> dict:
        full = {**self.prov, **row}
        line = json.dumps(full, sort_keys=True, separators=(",", ":"), allow_nan=False)
        with open(self.path, "a", encoding="utf-8") as handle:
            handle.write(line + "\n")
            handle.flush()
            os.fsync(handle.fileno())
        return full


def read_rows(path):
    rows = []
    with open(path, encoding="utf-8") as handle:
        for line in handle:
            line = line.strip()
            if line:
                rows.append(json.loads(line))
    return rows


def has_plan(path, run_id: str) -> bool:
    """Whether a plan row of this run is already in the record."""
    return Path(path).exists() and any(
        r.get("type") == "plan" and r.get("run_id") == run_id for r in read_rows(path))


# ----------------------------------------------------------------- driver


def pick_core(spec: str) -> str:
    if spec != "auto":
        return spec
    allowed = sorted(os.sched_getaffinity(0))
    return str(allowed[1] if len(allowed) > 1 else allowed[0])


def worker_env(arm: str) -> dict:
    env = {k: os.environ[k] for k in CARRIED_ENV if k in os.environ}
    env.update(PINNED_ENV)
    if TOOL_OF[arm] == "ovito" and os.environ.get(OVITO_LIBS_ENV):
        env["LD_LIBRARY_PATH"] = os.environ[OVITO_LIBS_ENV]
    return env


def cells_for(workloads, systems, frames: dict):
    out = []
    for workload in workloads:
        for system in systems:
            n = frames[system]["molecules"]
            arms = tuple(a for a in ARMS[workload] if applicable(a, workload, n))
            out.append((workload, system, arms))
    return out


def schedule(cells, rounds: int):
    """Yield (round, position, order, workload, system, arm)."""
    for r in range(rounds):
        flip = r % 2 == 1
        sequence = list(reversed(cells)) if flip else list(cells)
        for workload, system, arms in sequence:
            order = list(reversed(arms)) if flip else list(arms)
            for position, arm in enumerate(order):
                yield r, position, "rev" if flip else "fwd", workload, system, arm


def read_cpu_times() -> list[tuple[int, int]]:
    rows = []
    with open("/proc/stat", encoding="ascii") as handle:
        for line in handle:
            if line.startswith("cpu") and line[3].isdigit():
                parts = [int(v) for v in line.split()[1:]]
                idle = parts[3] + parts[4]
                rows.append((idle, sum(parts)))
    return rows


def idle_snapshot(core: int, window: float = 1.0) -> dict:
    first = read_cpu_times()
    time.sleep(window)
    second = read_cpu_times()
    idle = [
        100.0 * (b[0] - a[0]) / max(b[1] - a[1], 1) for a, b in zip(first, second)
    ]
    others = [v for i, v in enumerate(idle) if i != core]
    with open("/proc/loadavg", encoding="ascii") as handle:
        load = [float(v) for v in handle.read().split()[:3]]
    return {
        "loadavg": load,
        "pinned_core_idle_pct": round(idle[core], 2) if core < len(idle) else None,
        "other_cores_idle_pct_min": round(min(others), 2) if others else None,
    }


def launch(argv, env, timeout):
    start = time.perf_counter()
    proc = subprocess.Popen(
        argv, cwd=ROOT, env=env, stdin=subprocess.DEVNULL,
        stdout=subprocess.PIPE, stderr=subprocess.PIPE, text=True,
        start_new_session=True,
    )
    try:
        out, err = proc.communicate(timeout=timeout)
        status = "ok" if proc.returncode == 0 else "error"
    except subprocess.TimeoutExpired:
        os.killpg(proc.pid, signal.SIGKILL)
        out, err = proc.communicate()
        status = "timeout"
    return status, proc.returncode, out, err, (time.perf_counter() - start) * 1e3


def parse_worker(out: str):
    for line in reversed(out.splitlines()):
        line = line.strip()
        if line.startswith("{"):
            return json.loads(line)
    return None


def run_cells(args, workloads, record: Record, raw: Path, work: Path):
    systems = args.systems or list(SYSTEMS)
    frames = write_frames(work / "frames", systems)
    core = pick_core(args.core)
    python = args.python or sys.executable
    cells = cells_for(workloads, systems, frames)
    dead: set = set()
    done: set = set()
    if args.resume and record.path.exists():
        for row in read_rows(record.path):
            if row.get("type") == "sample" and row.get("run_id") == record.prov["run_id"]:
                done.add((row["round"], row["workload"], row["system"], row["arm"]))
                if row["status"] == "timeout":
                    dead.add((row["workload"], row["system"], row["arm"]))
    arm_filter = set(args.arms or [])
    raw.mkdir(parents=True, exist_ok=True)
    (raw / "stderr").mkdir(exist_ok=True)
    log = open(raw / "driver.log", "a", encoding="utf-8")
    seq = 0
    current_round = -1
    for r, position, order, workload, system, arm in schedule(cells, args.rounds):
        if arm_filter and arm not in arm_filter and arm != REFERENCE and arm != TWIN:
            continue
        if r != current_round:
            current_round = r
            record.append({"type": "round", "round": r, **idle_snapshot(int(core))})
        key = (workload, system, arm)
        if (r, workload, system, arm) in done:
            continue
        base_row = {
            "type": "sample", "round": r, "position": position, "order": order,
            "workload": workload, "system": system, "arm": arm,
            "molecules": frames[system]["molecules"], "cell": f"{workload}/{system}",
            "core": core,
        }
        if r >= ROUND_CAPS.get((arm, workload, system), args.rounds):
            record.append({**base_row, "status": "skipped",
                           "detail": f"capped at {ROUND_CAPS[(arm, workload, system)]} rounds"})
            continue
        if key in dead:
            record.append({**base_row, "status": "skipped",
                           "detail": "an earlier round of this cell timed out"})
            continue
        argv = ["taskset", "-c", core, python, str(Path(__file__).resolve()), "worker",
                "--arm", arm, "--workload", workload, "--system", system,
                "--frames", str(work / "frames")]
        if r == 0:
            argv += ["--canon", str(work / "canon")]
        stamp = datetime.datetime.now(datetime.timezone.utc).strftime("%Y-%m-%dT%H:%M:%SZ")
        status, code, out, err, process_ms = launch(argv, worker_env(arm), args.timeout)
        seq += 1
        parsed = parse_worker(out) if status == "ok" else None
        if status == "ok" and parsed is None:
            status = "error"
        label = f"r{r} {workload} {system} {arm}"
        if err.strip():
            (raw / "stderr" / f"{label.replace(' ', '_').replace('/', '-')}.txt").write_text(
                scrub(err), encoding="utf-8")
        row = {**base_row, "utc": stamp, "status": status, "process_wall_ms": round(process_ms, 1)}
        if status == "ok":
            row.update(parsed)
            row["position"] = position
            log.write(
                f"{label} order {order} pos {position} wall_ms {parsed['wall_ms']:.4f} "
                f"cpu_ms {parsed['cpu_ms']:.4f} rss_peak_kb {parsed['rss_peak_kb']} "
                f"rss_delta_kb {parsed['rss_peak_kb'] - parsed['rss_before_kb']}\n"
            )
        else:
            detail = scrub("\n".join(err.strip().splitlines()[-3:]))[:400]
            row.update({"detail": detail, "returncode": code, "timeout_s": args.timeout})
            if status == "timeout":
                dead.add(key)
            log.write(f"{label} order {order} pos {position} {status} {detail!r}\n")
        record.append(row)
        log.flush()
        if args.verbose:
            print(log_line(row), flush=True)
    log.close()
    return frames, cells


def log_line(row: dict) -> str:
    if row["status"] != "ok":
        return f"r{row['round']} {row['cell']} {row['arm']} {row['status']}"
    return (
        f"r{row['round']} {row['cell']} {row['arm']} wall_ms {row['wall_ms']:.3f} "
        f"peak_mib {row['rss_peak_kb'] / 1024:.1f}"
    )


def oracle_keys(workload: str, pos, lengths):
    """Exact answer by minimum-image brute force. Small frames only."""
    n = pos.shape[0]
    box = np.asarray(lengths)
    left, right = [], []
    for i in range(n):
        d = pos - pos[i]
        d -= box * np.round(d / box)
        r2 = np.einsum("ij,ij->i", d, d)
        r2[i] = np.inf
        if workload == "cutoff":
            near = np.nonzero(r2 < CUTOFF_A * CUTOFF_A)[0]
        else:
            reach = np.nonzero(r2 < KNN_REACH_A * KNN_REACH_A)[0]
            near = reach[np.argsort(r2[reach], kind="stable")[:KNN]]
        left.append(np.full(near.shape[0], i, dtype=np.int64))
        right.append(near.astype(np.int64))
    return pair_keys(np.concatenate(left), np.concatenate(right), n)


def agreement_rows(record: Record, work: Path, frames: dict, workloads) -> list[dict]:
    run_id = record.prov["run_id"]
    every = read_rows(record.path)
    rows = [r for r in every
            if r.get("type") == "sample" and r.get("run_id") == run_id and r["status"] == "ok"]
    written = {(r["workload"], r["system"], r["arm"]) for r in every
               if r.get("type") == "agreement" and r.get("run_id") == run_id}
    out = []
    cells = sorted({(r["workload"], r["system"]) for r in rows})
    for workload, system in cells:
        if workload not in workloads:
            continue
        group = [r for r in rows if r["workload"] == workload and r["system"] == system]
        by_arm: dict = {}
        for r in group:
            by_arm.setdefault(r["arm"], []).append(r)
        n = frames[system]["molecules"]
        oracle = None
        if workload != "rings" and n <= ORACLE_LIMIT:
            pos, lengths = load_frame(work / "frames", system)
            oracle = oracle_keys(workload, pos, lengths)
        ref_keys = None
        ref_path = work / "canon" / f"{REFERENCE}.{workload}.{system}.npy"
        if ref_path.exists():
            ref_keys = np.load(ref_path)
        ref_out = by_arm.get(REFERENCE, [{}])[0].get("output")
        for arm, samples in sorted(by_arm.items()):
            if (workload, system, arm) in written:
                continue
            digests = {json.dumps(s["output"], sort_keys=True) for s in samples}
            row = {
                "type": "agreement", "workload": workload, "system": system,
                "molecules": n, "arm": arm, "reference": REFERENCE,
                "samples": len(samples), "samples_identical": len(digests) == 1,
                "output_form": samples[0]["output"]["form"],
            }
            mine = samples[0]["output"]
            if workload == "rings":
                sizes = [s for s in (3, 4, 5, 6)
                         if str(s) in mine["hist"] and ref_out and str(s) in ref_out["hist"]]
                diffs = {
                    str(s): {"reference": ref_out["hist"][str(s)], "arm": mine["hist"][str(s)]}
                    for s in sizes if mine["hist"][str(s)] != ref_out["hist"][str(s)]
                }
                row.update({
                    "sizes_compared": sizes,
                    "equal": not diffs,
                    "differences": diffs,
                    "counts": mine["hist"],
                })
            else:
                keys_path = work / "canon" / f"{arm}.{workload}.{system}.npy"
                keys = np.load(keys_path) if keys_path.exists() else None
                row.update({
                    "n_pairs": mine["n_pairs"],
                    "directed_entries": mine["directed_entries"],
                    "row_len_hist": mine["row_len_hist"],
                    "equal": bool(ref_out and mine["pairs_digest"] == ref_out["pairs_digest"]),
                })
                if keys is not None and ref_keys is not None:
                    ref_set = np.isin(ref_keys, keys, assume_unique=True)
                    arm_set = np.isin(keys, ref_keys, assume_unique=True)
                    only_ref = np.nonzero(~ref_set)[0]
                    only_arm = np.nonzero(~arm_set)[0]
                    ref_deg = degrees(ref_keys, n)
                    arm_deg = degrees(keys, n)
                    row.update({
                        "pairs_only_reference": int(only_ref.shape[0]),
                        "pairs_only_arm": int(only_arm.shape[0]),
                        "atoms_degree_differ": int(np.count_nonzero(ref_deg != arm_deg)),
                    })
                if oracle is not None and keys is not None:
                    row["matches_brute_force"] = bool(np.array_equal(oracle, keys))
                    row["pairs_brute_force"] = int(oracle.shape[0])
            out.append(row)
    for row in out:
        record.append(row)
    return out


def run_head_to_heads(args) -> int:
    names = list(HEAD_TO_HEADS) if args.head_to_head == "all" else [args.head_to_head]
    workloads = [w for name in names for w in HEAD_TO_HEADS[name]]
    if args.workloads:
        workloads = [w for w in workloads if w in args.workloads]
    stamp = args.stamp or datetime.date.today().isoformat()
    run_id = args.run_id or uuid.uuid4().hex[:12]
    if args.rounds < args.min_rounds:
        raise SystemExit(f"need at least {args.min_rounds} rounds, got {args.rounds}")
    if args.dry_run:
        systems = args.systems or list(SYSTEMS)
        total = 0
        for name in systems:
            pos, lengths = build_frame(name)
            for workload in workloads:
                arms = [a for a in ARMS[workload] if applicable(a, workload, pos.shape[0])]
                total += len(arms) * args.rounds
                print(f"{workload:7s} {name:14s} {pos.shape[0]:7d} {' '.join(arms)}")
        print(f"samples {total}")
        return 0
    import release_build
    import h2h_freud_ovito as base

    lib = base.library_path()
    label = release_build.ensure_release(str(lib))[1]
    prov = provenance(run_id, stamp)
    prov["build"] = label
    record = Record(args.record, prov)
    work = Path(args.work or ROOT / "bbdir-h2h" / stamp / "work")
    raw = Path(args.raw or ROOT / "bbdir-h2h" / stamp / "raw")
    if not (args.resume and has_plan(record.path, run_id)):
        record.append({
            "type": "plan", "workloads": workloads, "rounds": args.rounds,
            "systems": args.systems or list(SYSTEMS), "arms": {w: list(ARMS[w]) for w in workloads},
            "limits": {f"{a}/{w}": v for (a, w), v in LIMITS.items()},
            "round_caps": {"/".join(k): v for k, v in ROUND_CAPS.items()},
            "arm_notes": ARM_NOTES,
            "timeout_s": args.timeout, "cutoff_A": CUTOFF_A, "knn": KNN,
            "knn_reach_A": KNN_REACH_A, "max_ring": MAX_RING,
            "env": PINNED_ENV, "protocol": (
                "one cold call per fresh process, pinned to one core, no warm-up, "
                "arms interleaved per round, odd rounds reversed, seams_aa is the A/A twin"
            ),
        })
        record.append(machine_row())
    frames, cells = run_cells(args, workloads, record, raw, work)
    agreement_rows(record, work, frames, workloads)
    if args.marker:
        marker = Path(args.marker)
        marker.parent.mkdir(parents=True, exist_ok=True)
        marker.write_text(f"{run_id} {args.head_to_head}\n", encoding="utf-8")
    return 0


# ------------------------------------------------------------------ lanes


@dataclasses.dataclass(frozen=True)
class Lane:
    """One Airspeed Velocity benchmark method, run cell by cell."""

    name: str
    group: str  # kernel, rings, scaling or h2h
    benchmark: str  # name as asv lists it, without the parameters
    cpus: str = "one"  # one: the lane core. threads: 0 .. t-1. all: every allowed cpu.
    ovito: bool = False
    twin: bool = True  # an A/A twin beside every sample
    only: tuple = ()  # ((parameter name, (repr, ...)), ...) keeps those values only
    note: str = ""


H2H_ASV_SYSTEMS = ("'ice512'", "'ice110592'")
H2H_ASV_RINGS_ARMS = ("'seams'", "'seams_aa'", "'ovito'", "'rsm'", "'mdapy'")

LANES = (
    Lane("cutoff-ice-neighbours", "kernel", "benchmarks.CutoffIce.time_neighbours"),
    Lane("cutoff-ice-chill-plus", "kernel", "benchmarks.CutoffIce.time_published_chill_plus"),
    Lane("cutoff-ice-q6", "kernel", "benchmarks.CutoffIce.time_q6"),
    Lane("chill-plus-1000", "kernel", "benchmarks.ChillPlusIce.time_published_chill_plus"),
    Lane("ddc-1000", "kernel", "benchmarks.DdcArm.time_find_ddc",
         note="wall time of one bench_ddc_ladder process"),
    Lane("cutoff-ice-rings", "rings", "benchmarks.CutoffIce.time_rings"),
    Lane("rings-stages-sH", "rings", "benchmarks.RingArm.time_primitive_rings",
         note="wall time of one bench_stages process: read, lists, rings, cages"),
    Lane("rings-ih-1t", "rings", "rings_ih.RingsIh.time_one_thread",
         note="wall time of one bench_rings_ih process"),
    Lane("rings-ih-8t", "rings", "rings_ih.RingsIh.time_eight_threads", cpus="all",
         note="wall time of one bench_rings_ih process"),
    Lane("rings-one-thread-dseams", "rings", "benchmarks.RingsOneThread.time_dseams",
         ovito=True, note="persistent worker, one census already done in setup"),
    Lane("rings-one-thread-ovito", "rings", "benchmarks.RingsOneThread.time_ovito_findrings",
         ovito=True, note="persistent worker, one census already done in setup"),
    Lane("omp-scaling", "scaling", "benchmarks.OmpScaling.time_call", cpus="threads"),
    Lane("omp-session", "scaling", "benchmarks.OmpScale.time_session", cpus="threads",
         note="one call is eight samples of knn, q6, CHILL and a dump read"),
    Lane("classify-ice", "scaling", "benchmarks.ClassifyIce.time_classify", cpus="threads",
         note="wall time of one bench_shell_scale process"),
    Lane("local-q4-8t", "scaling", "benchmarks.LocalQl.time_q4", cpus="all"),
    Lane("local-q6-8t", "scaling", "benchmarks.LocalQl.time_q6", cpus="all"),
    Lane("h2h-neighbours-time", "h2h", "campaign_h2h.NeighbourLists.time_call", twin=False,
         only=(("system", H2H_ASV_SYSTEMS),)),
    Lane("h2h-rings-time", "h2h", "campaign_h2h.RingStatistics.time_call", twin=False,
         only=(("system", H2H_ASV_SYSTEMS), ("arm", H2H_ASV_RINGS_ARMS))),
    Lane("h2h-neighbours-peakmem", "h2h", "campaign_h2h.NeighbourLists.peakmem_call",
         twin=False, only=(("system", ("'ice110592'",)),)),
    Lane("h2h-rings-peakmem", "h2h", "campaign_h2h.RingStatistics.peakmem_call", twin=False,
         only=(("system", ("'ice110592'",)), ("arm", H2H_ASV_RINGS_ARMS))),
)
LANE_GROUPS = ("kernel", "rings", "scaling", "h2h")
INVENTORY = "inventory"
THREAD_VARS = (
    "OMP_NUM_THREADS", "OMP_DYNAMIC", "OMP_PROC_BIND", "OMP_PLACES", "OMP_WAIT_POLICY",
    "RAYON_NUM_THREADS", "OVITO_THREAD_COUNT", "OPENBLAS_NUM_THREADS", "MKL_NUM_THREADS",
    "NUMEXPR_NUM_THREADS",
)
ASV_ATTRIBUTES = ("number=1", "repeat=1", "processes=1", "warmup_time=0")


def select_lanes(spec: str) -> list[Lane]:
    wanted = [s for s in spec.split(",") if s] or ["kernel", "rings", "scaling"]
    chosen = []
    for lane in LANES:
        if lane.name in wanted or lane.group in wanted or "all" in wanted:
            chosen.append(lane)
    unknown = [w for w in wanted if w not in ("all", INVENTORY) and w not in LANE_GROUPS
               and w not in {lane.name for lane in LANES}]
    if unknown:
        raise SystemExit(f"unknown lane or group: {', '.join(unknown)}")
    return chosen


def inventory_lanes(spec: str, found: dict) -> list[Lane]:
    """Every benchmark asv lists that no lane covers, for one run each."""
    if INVENTORY not in spec.split(","):
        return []
    covered = {lane.benchmark for lane in LANES}
    return [Lane(name, INVENTORY, name, cpus="all", ovito=True, twin=False)
            for name in sorted(found) if name not in covered]


def asv_executable_dir() -> Path:
    spec = importlib.util.find_spec("asv")
    if spec is None or not spec.submodule_search_locations:
        raise SystemExit("asv is not importable")
    return Path(next(iter(spec.submodule_search_locations)))


def discover_benchmarks(python: str, benchmark_dir: Path, scratch: Path) -> dict:
    """The benchmarks asv finds, by name, with their parameter axes."""
    out = scratch / "discovered.json"
    done = subprocess.run(
        [python, str(asv_executable_dir() / "benchmark.py"), "discover",
         str(benchmark_dir), str(out)],
        cwd=ROOT, capture_output=True, text=True, check=False,
        env={**{k: os.environ[k] for k in CARRIED_ENV if k in os.environ}, **PINNED_ENV},
    )
    if done.returncode != 0 or not out.is_file():
        raise SystemExit(f"asv discovery failed: {scrub(done.stderr[-400:])}")
    return {b["name"]: b for b in json.loads(out.read_text(encoding="utf-8"))}


def h2h_cell_runs(params: dict) -> bool:
    """Whether the campaign has an arm for this workload, system and tool."""
    value = {k: ast.literal_eval(v) for k, v in params.items()}
    workload = value.get("workload", "rings")
    system, arm = value["system"], value["arm"]
    if arm not in ARMS[workload]:
        return False
    molecules = ICE_CELLS[system] ** 3 * 8 if system in ICE_CELLS else None
    return applicable(arm, workload, molecules)


PROBE_CELLS = 3


def probe(cells: list[dict]) -> list[dict]:
    """The first, middle and last cell of a long parameter grid."""
    if len(cells) <= PROBE_CELLS:
        return cells
    return [cells[i] for i in sorted({0, len(cells) // 2, len(cells) - 1})]


def lane_cells(lane: Lane, found: dict) -> list[dict]:
    bench = found.get(lane.benchmark)
    if bench is None:
        raise SystemExit(f"asv does not list {lane.benchmark}")
    names = list(bench.get("param_names") or [])
    axes = [list(axis) for axis in (bench.get("params") or [])]
    keep = dict(lane.only)
    cells = []
    for index, combo in enumerate(itertools.product(*axes) if axes else [()]):
        if any(name in keep and value not in keep[name]
               for name, value in zip(names, combo)):
            continue
        params = dict(zip(names, combo))
        if lane.group == "h2h" and not h2h_cell_runs(params):
            continue
        label = f"{lane.benchmark}({', '.join(combo)})" if combo else lane.benchmark
        cells.append({
            "lane": lane.name, "group": lane.group, "benchmark": lane.benchmark,
            "index": index, "label": label,
            "type": "peakmem" if bench.get("type") == "peakmemory" else bench.get("type", "time"),
            "params": params,
            "version": bench.get("version"),
        })
    return cells


def cpu_list(cpus: list[int]) -> str:
    if len(cpus) > 1 and cpus == list(range(cpus[0], cpus[-1] + 1)):
        return f"{cpus[0]}-{cpus[-1]}"
    return ",".join(str(c) for c in cpus)


def cell_cpus(lane: Lane, cell: dict, core: str, allowed: list[int]) -> str | None:
    """The --cpu-affinity value for a cell, or None when the machine is too small."""
    if lane.cpus == "one":
        return core
    if lane.cpus == "all":
        return cpu_list(allowed)
    threads = int(cell["params"].get("threads", "1"))
    if threads > len(allowed):
        return None
    return cpu_list(allowed[:threads])


def needs_ovito(lane: Lane, cell: dict) -> bool:
    """The OVITO library path is added for the cells that import OVITO and no others."""
    return lane.ovito or any("ovito" in value for value in cell["params"].values())


def lane_env(lane: Lane, cell: dict, home: Path, build: Path) -> dict:
    env = {k: os.environ[k] for k in CARRIED_ENV if k in os.environ}
    env["HOME"] = str(home)
    for name in THREAD_VARS:
        env.pop(name, None)
    env.update({
        "OPENBLAS_NUM_THREADS": "1", "MKL_NUM_THREADS": "1", "NUMEXPR_NUM_THREADS": "1",
        "QT_QPA_PLATFORM": "offscreen", "PYTHONHASHSEED": "0", "PYTHONWARNINGS": "ignore",
        "FASTATOMSTRUCT_NO_TELEMETRY": "1", "FASTATOMSTRUCT_SUPPRESS_TELEMETRY_NOTICE": "1",
        "SEAMS_BENCH_PREWARM": "0",
        "SEAMS_LANG_LIB": str(build / "bindings" / "libseams_lang.so"),
        "SEAMS_CORE_LIB": str(build / "src" / "libyodaLib.so"),
        "SEAMS_LIB": str(build / "src" / "libseams.so"),
        "SEAMS_RING_BENCH": str(build / "tests" / "bench_rings_h2h"),
    })
    if lane.cpus == "one":
        env.update({
            "OMP_NUM_THREADS": "1", "OMP_DYNAMIC": "FALSE", "RAYON_NUM_THREADS": "1",
            "OVITO_THREAD_COUNT": "1",
        })
    if lane.group == "h2h":
        env.update(PINNED_ENV)
    if needs_ovito(lane, cell) and os.environ.get(OVITO_LIBS_ENV):
        env["LD_LIBRARY_PATH"] = os.environ[OVITO_LIBS_ENV]
    return env


def write_asv_files(scratch: Path, host: str) -> tuple[Path, Path, Path]:
    """A config, an isolated home with a machine file, and a results directory."""
    home = scratch / "home"
    results = scratch / "results"
    home.mkdir(parents=True, exist_ok=True)
    results.mkdir(parents=True, exist_ok=True)
    fields = lscpu_fields(lscpu_lines())
    (home / ".asv-machine.json").write_text(json.dumps({
        host: {
            "arch": platform.machine(), "cpu": fields.get("Model name", "unknown"),
            "machine": host, "num_cpu": fields.get("CPU(s)", "unknown"),
            "os": f"{platform.system()} {platform.release()}",
            "ram": str(int(next((l.split()[1] for l in read_text("/proc/meminfo").splitlines()
                                 if l.startswith("MemTotal")), 0))),
        },
        "version": 1,
    }, indent=2), encoding="utf-8")
    conf = scratch / "asv.conf.json"
    conf.write_text(json.dumps({
        "version": 1, "project": "d-seams", "repo": str(ROOT), "branches": ["main"],
        "environment_type": "existing", "env_dir": str(scratch / "env"),
        "results_dir": str(results), "html_dir": str(scratch / "html"),
        "benchmark_dir": str(ROOT / "benchmarks"),
        "build_command": [], "install_command": [], "uninstall_command": [],
    }, indent=2), encoding="utf-8")
    return conf, home, results


def cell_regex(cell: dict) -> str:
    return "^" + re.escape(cell["label"]) + "$"


def clear_results(results: Path, host: str) -> None:
    """Asv merges into the file of the same commit. Each sample starts empty."""
    for stale in (results / host).glob("*-existing-*.json"):
        stale.unlink()


def harvest(results: Path, host: str, cell: dict):
    """The one value asv stored for a cell, or None."""
    files = sorted((results / host).glob("*-existing-*.json"))
    if not files:
        return None
    data = json.loads(files[0].read_text(encoding="utf-8"))
    entry = data.get("results", {}).get(cell["benchmark"])
    if entry is None:
        return None
    rec = dict(zip(data["result_columns"], entry))
    index = cell["index"]
    samples, result = rec.get("samples"), rec.get("result")
    samples = samples[index] if samples and index < len(samples) else None
    result = result[index] if isinstance(result, list) and index < len(result) else None
    if samples:
        return {"value": float(samples[0]), "samples": len(samples)}
    if isinstance(result, (int, float)) and not math.isnan(result):
        return {"value": float(result), "samples": 1}
    return None


def run_lanes(args) -> int:
    import release_build

    stamp = args.stamp or datetime.date.today().isoformat()
    run_id = args.run_id or uuid.uuid4().hex[:12]
    if args.rounds < args.min_rounds:
        raise SystemExit(f"need at least {args.min_rounds} rounds, got {args.rounds}")
    python = args.python or sys.executable
    allowed = sorted(os.sched_getaffinity(0))
    core = args.core if args.core != "auto" else str(allowed[0])
    work = Path(args.work or ROOT / "bbdir-h2h" / stamp / "work")
    raw = Path(args.raw or ROOT / "bbdir-h2h" / stamp / "raw")
    scratch = work / "asv"
    scratch.mkdir(parents=True, exist_ok=True)
    build = release_build.BUILD
    label = release_build.assert_release(build)
    found = discover_benchmarks(python, ROOT / "benchmarks", scratch)
    lanes = select_lanes(args.lanes) + inventory_lanes(args.lanes, found)
    cells = [cell for lane in lanes
             for cell in (probe(lane_cells(lane, found)) if lane.group == INVENTORY
                          else lane_cells(lane, found))]
    if args.cells:
        cells = [cell for cell in cells if re.search(args.cells, cell["label"])]
    by_name = {lane.name: lane for lane in lanes}
    if args.dry_run:
        for cell in cells:
            lane = by_name[cell["lane"]]
            print(f"{cell['lane']:26s} cpus {cell_cpus(lane, cell, core, allowed)} {cell['label']}")
        print(f"cells {len(cells)} invocations {sum(2 if by_name[c['lane']].twin else 1 for c in cells) * args.rounds}")
        return 0
    host = host_id()
    revision = head_revision()
    conf, home, results = write_asv_files(scratch, host)
    prov = provenance(run_id, stamp)
    prov["build"] = label
    record = Record(args.record, prov)
    raw.mkdir(parents=True, exist_ok=True)
    (raw / "stderr").mkdir(exist_ok=True)
    done: set = set()
    dead: dict = {}
    if args.resume and record.path.exists():
        for row in read_rows(record.path):
            if row.get("type") == "asv" and row.get("run_id") == run_id:
                done.add((row["round"], row["cell"], row["arm"]))
                if row["status"] in ("failed", "skipped", "timeout"):
                    dead[row["cell"]] = row["status"]
    if not (args.resume and has_plan(record.path, run_id)):
        record.append({
            "type": "plan", "kind": "lanes", "lanes": [lane.name for lane in lanes],
            "cells": len(cells), "rounds": args.rounds, "core": core,
            "asv": importlib.metadata.version("asv"),
            "asv_runner": importlib.metadata.version("asv_runner"),
            "attributes": list(ASV_ATTRIBUTES), "timeout_s": args.timeout,
            "prewarm": "SEAMS_BENCH_PREWARM=0: setup warm-up calls dropped",
            "protocol": (
                "one asv invocation per sample, number=1 repeat=1 processes=1 warmup_time=0, "
                "cells in order, odd rounds reversed, an A/A twin of every cell where twin is true, "
                "no quick mode"
            ),
            "notes": {lane.name: lane.note for lane in lanes if lane.note},
        })
        record.append(machine_row())

    def attempt(cell, lane, r, arm, position, order):
        cpus = cell_cpus(lane, cell, core, allowed)
        row = {
            "type": "asv", "lane": lane.name, "group": lane.group, "benchmark": cell["benchmark"],
            "cell": cell["label"], "params": cell["params"], "kind": cell["type"],
            "round": r, "position": position, "order": order, "arm": arm,
            "asv_benchmark_version": cell["version"],
        }
        if cpus is None:
            record.append({**row, "status": "skipped",
                           "detail": "more threads than allowed cpus"})
            dead[cell["label"]] = "skipped"
            return
        if cell["label"] in dead:
            record.append({**row, "status": "skipped", "detail": f"earlier round: {dead[cell['label']]}"})
            return
        row["cpus"] = cpus
        argv = [python, "-m", "asv", "--config", str(conf), "run", "--machine", host,
                "--python=same", "--set-commit-hash", revision, "--record-samples",
                "--show-stderr", "--cpu-affinity", cpus, "-b", cell_regex(cell)]
        for attribute in ASV_ATTRIBUTES:
            argv += ["-a", attribute]
        argv += ["-a", f"timeout={int(args.timeout)}"]
        env = lane_env(by_name[cell["lane"]], cell, home, build)
        utc = datetime.datetime.now(datetime.timezone.utc).strftime("%Y-%m-%dT%H:%M:%SZ")
        clear_results(results, host)
        status, code, out, err, process_ms = launch(argv, env, args.timeout * 2 + 120)
        found_value = harvest(results, host, cell) if status == "ok" else None
        text = out + "\n" + err
        label = f"r{r} {cell['label']} {arm}"
        if err.strip():
            (raw / "stderr" / (re.sub(r"[^\w.-]+", "_", label) + ".txt")).write_text(
                scrub(err), encoding="utf-8")
        row.update({"utc": utc, "process_wall_ms": round(process_ms, 1)})
        if found_value is not None:
            scale = 1.0 if cell["type"] == "peakmem" else 1e3
            row.update({"status": "ok", "value": found_value["value"],
                        "unit": "bytes" if cell["type"] == "peakmem" else "s",
                        "wall_ms": found_value["value"] * scale if cell["type"] == "time" else None})
            log.write(f"{label} order {order} pos {position} value {found_value['value']:.9g}\n")
        else:
            reason = "skipped" if "asv: skipped" in text else (
                "timeout" if status == "timeout" else "failed")
            detail = scrub("\n".join(line for line in text.splitlines()
                                     if line.strip())[-500:])
            row.update({"status": reason, "detail": detail, "returncode": code})
            dead[cell["label"]] = reason
            log.write(f"{label} order {order} pos {position} {reason} {detail[-200:]!r}\n")
        record.append(row)
        log.flush()
        if args.verbose:
            print(f"r{r} {cell['label']} {arm} {row['status']} "
                  f"{row.get('value', '')}", flush=True)

    log = open(raw / "driver.log", "a", encoding="utf-8")
    current_round = -1
    for r in range(args.rounds):
        flip = r % 2 == 1
        sequence = list(reversed(cells)) if flip else list(cells)
        if r != current_round:
            current_round = r
            record.append({"type": "round", "kind": "lanes", "round": r,
                           **idle_snapshot(int(core.split("-")[0].split(",")[0]))})
        for cell in sequence:
            lane = by_name[cell["lane"]]
            arms = ["A", "A2"] if lane.twin else ["A"]
            order = list(reversed(arms)) if flip else arms
            for position, arm in enumerate(order):
                if (r, cell["label"], arm) in done:
                    continue
                attempt(cell, lane, r, arm, position, "rev" if flip else "fwd")
    log.close()
    if args.marker:
        marker = Path(args.marker)
        marker.parent.mkdir(parents=True, exist_ok=True)
        marker.write_text(f"{run_id} lanes {args.lanes}\n", encoding="utf-8")
    return 0


# ---------------------------------------------------------------- summary


def quantile(sorted_vals, q):
    return proto.quantile_linear(sorted_vals, q)


def stats(values) -> dict:
    ordered = sorted(float(v) for v in values)
    n = len(ordered)
    if not n:
        return {"n": 0}
    mid = n // 2
    median = 0.5 * (ordered[mid - 1] + ordered[mid]) if n % 2 == 0 else ordered[mid]
    q1, q3 = quantile(ordered, 0.25), quantile(ordered, 0.75)
    return {"n": n, "median": median, "q1": q1, "q3": q3, "iqr": q3 - q1,
            "min": ordered[0], "max": ordered[-1]}


def sample_rows(rows, run_id=None):
    return [r for r in rows if r.get("type") == "sample" and r.get("status") == "ok"
            and (run_id is None or r.get("run_id") == run_id)]


def per_round(rows) -> dict:
    out: dict = {}
    for r in rows:
        out.setdefault((r["workload"], r["system"], r["arm"]), {})[r["round"]] = r
    return out


def ratio_band(pairs) -> dict:
    ratios = sorted(pairs)
    if not ratios:
        return {"n": 0}
    return {"n": len(ratios), "median": stats(ratios)["median"],
            "min": ratios[0], "max": ratios[-1]}


def cell_summary(rows, run_id=None) -> list[dict]:
    good = sample_rows(rows, run_id)
    table = per_round(good)
    cells = sorted({(w, s) for (w, s, _a) in table})
    out = []
    for workload, system in cells:
        mine = {a: table[(w, s, a)] for (w, s, a) in table if (w, s) == (workload, system)}
        ref = mine.get(REFERENCE, {})
        twin = mine.get(TWIN, {})
        shared = sorted(set(ref) & set(twin))
        aa = ratio_band([twin[r]["wall_ms"] / ref[r]["wall_ms"] for r in shared])
        for arm, byround in sorted(mine.items()):
            walls = [x["wall_ms"] for x in byround.values()]
            cpus = [x["cpu_ms"] for x in byround.values()]
            peaks = [x["rss_peak_kb"] for x in byround.values()]
            deltas = [x["rss_peak_kb"] - x["rss_before_kb"] for x in byround.values()]
            shared_ref = sorted(set(ref) & set(byround))
            paired = ratio_band([byround[r]["wall_ms"] / ref[r]["wall_ms"] for r in shared_ref])
            entry = {
                "workload": workload, "system": system, "arm": arm,
                "molecules": next(iter(byround.values()))["molecules"],
                "wall_ms": stats(walls), "cpu_ms": stats(cpus),
                "rss_peak_kb": stats(peaks), "rss_delta_kb": stats(deltas),
                "os_threads": max(x["os_threads"] for x in byround.values()),
                "paired_vs_reference": paired, "aa": aa,
            }
            if arm not in (REFERENCE, TWIN) and paired["n"] and aa["n"]:
                if paired["min"] > aa["max"]:
                    entry["verdict"] = "seams faster"
                elif paired["max"] < aa["min"]:
                    entry["verdict"] = "competitor faster"
                else:
                    entry["verdict"] = "within A/A band"
            out.append(entry)
    return out


def fmt_ms(s: dict) -> str:
    if not s.get("n"):
        return "-"
    return f"{s['median']:.3f} ({s['q1']:.3f}-{s['q3']:.3f})"


def render_summary(rows, run_id=None) -> str:
    cells = cell_summary(rows, run_id)
    out = []
    head = next((r for r in rows if r.get("type") in ("sample", "asv")), None)
    if head:
        out.append(f"host {head.get('host', '?')}; {head.get('lscpu', '?')}; "
                   f"kernel {head.get('kernel', '?')}")
        out.append(f"code_revision {head.get('code_revision', '?')}; "
                   f"build {head.get('build', '')}")
    out.append("")
    for workload in ("cutoff", "knn", "rings"):
        subset = [c for c in cells if c["workload"] == workload]
        if not subset:
            continue
        out.append(f"## {workload}: wall ms, median (Q1-Q3), n per cell")
        out.append("")
        out.append("| system | molecules | arm | n | wall ms | x seams | A/A band | "
                   "peak MiB | delta MiB | cpu ms | threads | verdict |")
        out.append("|---|---|---|---|---|---|---|---|---|---|---|---|")
        for c in subset:
            ref = next((x for x in subset
                        if x["system"] == c["system"] and x["arm"] == REFERENCE), None)
            ratio = (c["wall_ms"]["median"] / ref["wall_ms"]["median"]
                     if ref and ref["wall_ms"].get("n") else float("nan"))
            band = c["aa"]
            band_text = (f"{band['min']:.3f}-{band['max']:.3f}" if band.get("n") else "-")
            out.append(
                f"| {c['system']} | {c['molecules']} | {c['arm']} | {c['wall_ms']['n']} | "
                f"{fmt_ms(c['wall_ms'])} | {ratio:.2f} | {band_text} | "
                f"{c['rss_peak_kb']['median'] / 1024:.1f} | "
                f"{c['rss_delta_kb']['median'] / 1024:.1f} | "
                f"{c['cpu_ms']['median']:.3f} | {c['os_threads']} | "
                f"{c.get('verdict', '')} |"
            )
        out.append("")
    agreement = [r for r in rows if r.get("type") == "agreement"
                 and (run_id is None or r.get("run_id") == run_id)]
    if agreement:
        out.append("## agreement with seams (same cell, round 0 output)")
        out.append("")
        out.append("| workload | system | arm | equal | detail |")
        out.append("|---|---|---|---|---|")
        for r in agreement:
            if r["arm"] in (REFERENCE, TWIN):
                continue
            if r["workload"] == "rings":
                detail = (f"sizes {r['sizes_compared']}; " +
                          (json.dumps(r["differences"], sort_keys=True) or "")
                          if not r["equal"] else f"sizes {r['sizes_compared']}")
            else:
                detail = (f"pairs {r['n_pairs']}; only reference {r.get('pairs_only_reference', '?')}; "
                          f"only arm {r.get('pairs_only_arm', '?')}; "
                          f"rows {json.dumps(r['row_len_hist'], sort_keys=True)}")
                if "matches_brute_force" in r:
                    detail += f"; brute force {'match' if r['matches_brute_force'] else 'differs'}"
            out.append(f"| {r['workload']} | {r['system']} | {r['arm']} | "
                       f"{'yes' if r['equal'] else 'NO'} | {detail} |")
        out.append("")
    gaps: dict = {}
    for r in rows:
        if (r.get("type") == "sample" and r.get("status") != "ok"
                and (run_id is None or r.get("run_id") == run_id)):
            gaps.setdefault((r["workload"], r["system"], r["arm"]), []).append(r)
    if gaps:
        out.append("## samples that did not complete")
        out.append("")
        out.append("| workload | system | arm | samples | status | detail |")
        out.append("|---|---|---|---|---|---|")
        for (workload, system, arm), items in sorted(gaps.items()):
            status = "/".join(sorted({i["status"] for i in items}))
            detail = next((i.get("detail", "") for i in items if i.get("detail")), "")
            detail = " ".join(detail.split())[:110].replace("|", "/")
            out.append(f"| {workload} | {system} | {arm} | {len(items)} | {status} | {detail} |")
        out.append("")
    present = {(c["arm"], c["workload"]) for c in cells}
    notes = [(key, text) for key, text in ARM_NOTES.items() if tuple(key.split("/")) in present]
    if notes:
        out.append("## what the timed call covers")
        out.append("")
        for key, text in notes:
            out.append(f"- {key}: {text}")
        out.append("")
    return "\n".join(out) + "\n"


def lane_stat_rows(rows, run_id=None) -> list[dict]:
    asv = [r for r in rows if r.get("type") == "asv"
           and (run_id is None or r.get("run_id") == run_id)]
    cells: dict = {}
    for r in asv:
        cells.setdefault((r["lane"], r["cell"]), []).append(r)
    out = []
    for (lane, cell), group in cells.items():
        head = group[0]
        entry = {"lane": lane, "group": head["group"], "cell": cell,
                 "params": head["params"], "kind": head["kind"]}
        ok = [r for r in group if r["status"] == "ok"]
        if not ok:
            entry["status"] = "/".join(sorted({r["status"] for r in group}))
            entry["detail"] = next((r.get("detail", "") for r in group if r.get("detail")), "")
            out.append(entry)
            continue
        scale = 1.0 / (1024 * 1024) if head["kind"] == "peakmem" else 1e3
        arms: dict = {}
        for r in ok:
            arms.setdefault(r["arm"], {})[r["round"]] = r["value"] * scale
        entry["unit"] = "MiB" if head["kind"] == "peakmem" else "ms"
        entry["A"] = stats(list(arms.get("A", {}).values()))
        if "A2" in arms:
            entry["A2"] = stats(list(arms["A2"].values()))
            shared = sorted(set(arms["A"]) & set(arms["A2"]))
            entry["aa"] = ratio_band([arms["A2"][r] / arms["A"][r] for r in shared if arms["A"][r]])
        out.append(entry)
    return out


def render_lanes(rows, run_id=None) -> str:
    cells = lane_stat_rows(rows, run_id)
    if not cells:
        return ""
    out = []
    for group, unit in itertools.product(LANE_GROUPS, ("ms", "MiB")):
        subset = [c for c in cells if c["group"] == group and "A" in c and c["unit"] == unit]
        if not subset:
            continue
        out.append(f"## asv lanes: {group}, {unit}, median (Q1-Q3); A and A2 are the same lane twice")
        out.append("")
        out.append("| lane | cell | n | A | n | A2 | A2/A min-max |")
        out.append("|---|---|---|---|---|---|---|")
        for c in subset:
            two = c.get("A2")
            band = c.get("aa", {})
            label = ", ".join(f"{k}={v}" for k, v in c["params"].items()) or "-"
            out.append(
                f"| {c['lane']} | {label} | {c['A']['n']} | {fmt_ms(c['A'])} | "
                f"{two['n'] if two else '-'} | {fmt_ms(two) if two else '-'} | "
                f"{(str(round(band['min'], 3)) + '-' + str(round(band['max'], 3))) if band.get('n') else '-'} |"
            )
        out.append("")
    speed: dict = {}
    for c in cells:
        if "A" not in c or "threads" not in c["params"]:
            continue
        rest = tuple(sorted((k, v) for k, v in c["params"].items() if k != "threads"))
        speed.setdefault((c["lane"], rest), {})[int(c["params"]["threads"])] = c["A"]["median"]
    if speed:
        out.append("## asv lanes: speed-up over one thread, ratio of A medians")
        out.append("")
        counts = sorted({t for v in speed.values() for t in v})
        out.append("| lane | cell | " + " | ".join(f"{t} thr" for t in counts) + " |")
        out.append("|---|---|" + "---|" * len(counts))
        for (lane, rest), by_threads in speed.items():
            base = by_threads.get(1)
            label = ", ".join(f"{k}={v}" for k, v in rest) or "-"
            cols = [f"{base / by_threads[t]:.2f}" if base and t in by_threads else "-"
                    for t in counts]
            out.append(f"| {lane} | {label} | " + " | ".join(cols) + " |")
        out.append("")
    failed = [c for c in cells if "A" not in c and c["group"] != INVENTORY]
    if failed:
        out.append("## asv lanes: cells with no sample")
        out.append("")
        for c in failed:
            out.append(f"- {c['lane']} {c['cell']}: {c['status']}: {c['detail'][-160:]}")
        out.append("")
    listed = sorted((c for c in cells if c["group"] == INVENTORY), key=lambda c: c["cell"])
    if listed:
        out.append("## asv inventory: every other benchmark in the suite, one run of the first, "
                   "middle and last cell of its grid")
        out.append("")
        out.append("| benchmark | status | value | detail |")
        out.append("|---|---|---|---|")
        for c in listed:
            if "A" in c:
                value = f"{c['A']['median']:.3f} {c['unit']}"
                out.append(f"| {c['cell']} | ran | {value} | |")
            else:
                detail = " ".join(c["detail"].split())[-120:].replace("|", "/")
                out.append(f"| {c['cell']} | {c['status']} | - | {detail} |")
        out.append("")
    return "\n".join(out) + "\n"


def summary_command(args) -> int:
    rows = read_rows(args.record)
    text = render_summary(rows, args.run_id) + render_lanes(rows, args.run_id)
    if args.out:
        Path(args.out).parent.mkdir(parents=True, exist_ok=True)
        Path(args.out).write_text(text, encoding="utf-8")
    else:
        sys.stdout.write(text)
    return 0


# -------------------------------------------------------------------- CLI


def main(argv=None) -> int:
    parser = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    sub = parser.add_subparsers(dest="command", required=True)

    w = sub.add_parser("worker", help="one cold sample; prints one JSON line")
    w.add_argument("--arm", required=True, choices=sorted(TOOL_OF))
    w.add_argument("--workload", required=True, choices=sorted(ARMS))
    w.add_argument("--system", required=True)
    w.add_argument("--frames", required=True)
    w.add_argument("--canon", default="")

    r = sub.add_parser("run", help="interleaved rounds over every cell")
    r.add_argument("--head-to-head", choices=("neighbours", "rings", "all"), default="all")
    r.add_argument("--workloads", nargs="*")
    r.add_argument("--systems", nargs="*", choices=SYSTEMS)
    r.add_argument("--arms", nargs="*", choices=sorted(TOOL_OF))
    r.add_argument("--rounds", type=int, default=ROUNDS)
    r.add_argument("--min-rounds", type=int, default=MIN_ROUNDS, help=argparse.SUPPRESS)
    r.add_argument("--core", default="auto")
    r.add_argument("--timeout", type=float, default=TIMEOUT_S)
    r.add_argument("--stamp", default="")
    r.add_argument("--run-id", default="")
    r.add_argument("--record", required=True)
    r.add_argument("--raw", default="")
    r.add_argument("--work", default="")
    r.add_argument("--marker", default="")
    r.add_argument("--python", default="")
    r.add_argument("--resume", action="store_true")
    r.add_argument("--dry-run", action="store_true")
    r.add_argument("--verbose", action="store_true")

    n = sub.add_parser("lanes", help="asv lanes, one invocation per sample, interleaved rounds")
    n.add_argument("--lanes", default="kernel,rings,scaling",
                   help="lane names or groups: kernel, rings, scaling, h2h, all")
    n.add_argument("--rounds", type=int, default=ROUNDS)
    n.add_argument("--min-rounds", type=int, default=MIN_ROUNDS, help=argparse.SUPPRESS)
    n.add_argument("--core", default="auto")
    n.add_argument("--cells", default="", help="keep the cells whose label matches this regex")
    n.add_argument("--timeout", type=float, default=TIMEOUT_S)
    n.add_argument("--stamp", default="")
    n.add_argument("--run-id", default="")
    n.add_argument("--record", required=True)
    n.add_argument("--raw", default="")
    n.add_argument("--work", default="")
    n.add_argument("--marker", default="")
    n.add_argument("--python", default="")
    n.add_argument("--resume", action="store_true")
    n.add_argument("--dry-run", action="store_true")
    n.add_argument("--verbose", action="store_true")

    s = sub.add_parser("summary", help="tables from a record")
    s.add_argument("--record", required=True)
    s.add_argument("--run-id", default=None)
    s.add_argument("--out", default="")

    args = parser.parse_args(argv)
    if args.command == "worker":
        return run_worker(args)
    if args.command == "run":
        return run_head_to_heads(args)
    if args.command == "lanes":
        return run_lanes(args)
    return summary_command(args)


if __name__ == "__main__":
    raise SystemExit(main())
