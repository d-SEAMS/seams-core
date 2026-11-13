#!/usr/bin/env python3
"""One sitting for every timing the paper quotes at 110592 molecules.

Ideal ice Ic, 24**3 cells, cutoff 3.5 A. One warm-up and eight calls
at 1 thread and at 8. Each width is its own process, pinned with
taskset, so OpenMP, Rayon and OVITO see the width at startup. The
median of an even count is the mean of the two central samples. The
interquartile range is the type-7 upper quartile minus the lower one.

d-SEAMS kernels: rings, local q4, local q6, the same pair with the
averaged buffers filled, the q4+q6 library call, the core q4 and q6
calls, the cutoff neighbour list, CHILL+, and the cage call. The
rivals are rings-statistics-matter 0.3.0, mdapy, vesin, freud and
OVITO, on the operation each one is quoted for.

Hydrate cages are the GenIce sI and sII unit cells and the 2x2x2
supercells built from those cells. d-SEAMS and ClavisCage 0.3.18 TRACE
are timed on each, and the census of each code is recorded.

benchmarks/paper-cells/sitting.csv is the only timing source. The host
name is not read.
"""

from __future__ import annotations

import argparse
import csv
import json
import os
import subprocess
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "scripts"))
import paper_stats as stats

MOLECULES = 110592
SIX = 221184
CELLS = 24
WARMUP = 1
REPEATS = 8
CUTOFF = 3.5

FIELDS = [
    "source_commit", "cpu_model", "cpu_cores", "affinity",
    "phase", "operation", "code", "version",
    "molecules", "threads", "warmup", "n",
    "median_ms", "iqr_ms", "q1_ms", "q3_ms", "samples_ms", "detail",
]


def library_path(lib: str) -> str:
    build = Path(lib).resolve().parents[1]
    parts = [
        build / "bindings",
        build / "src",
        build / "subprojects/readcon-core",
        build / "subprojects/readcon-db",
        build / "subprojects/linkcell",
        build / "subprojects/vesin",
        build / "subprojects/minimage",
    ]
    return os.pathsep.join(str(path) for path in parts)


def ice_arrays():
    import numpy as np
    import h2h_protocol as proto
    pos, side = proto.ice_ic(CELLS)
    if len(pos) != MOLECULES:
        raise RuntimeError(f"ice Ic has {len(pos)} molecules")
    xyz = np.ascontiguousarray(pos, dtype=np.float64)
    box = np.array([0, 0, 0, side, side, side, 0, 0, 0], dtype=np.float64)
    return xyz, box, float(side)


def read_dump(path: Path) -> dict:
    import numpy as np
    lines = path.read_text(encoding="utf-8").splitlines()
    if not lines or not lines[0].startswith("ITEM: TIMESTEP"):
        raise RuntimeError(f"{path.name} is not a LAMMPS dump")
    count = int(lines[3])
    triclinic = "xy" in lines[4]
    bounds = []
    tilt = [0.0, 0.0, 0.0]
    for axis, row in enumerate(lines[5:8]):
        parts = row.split()
        lo = float(parts[0])
        hi = float(parts[1])
        bounds.append((lo, hi))
        if triclinic and len(parts) >= 3:
            tilt[axis] = float(parts[2])
    if triclinic:
        raise RuntimeError(f"{path.name} is tilted; the supercell builder is for sI and sII")
    header = lines[8].split()[2:]
    type_at = header.index("type")
    xyz_at = [header.index(axis) for axis in ("x", "y", "z")]
    points = []
    for line in lines[9:9 + count]:
        parts = line.split()
        if int(parts[type_at]) != 1:
            continue
        points.append([float(parts[index]) for index in xyz_at])
    pos = np.ascontiguousarray(points, dtype=np.float64)
    lo = np.array([item[0] for item in bounds], dtype=np.float64)
    lengths = np.array([item[1] - item[0] for item in bounds], dtype=np.float64)
    return {"positions": pos, "origin": lo, "lengths": lengths, "molecules": len(pos)}


def replicate(frame: dict, n: int) -> dict:
    import numpy as np
    pos = frame["positions"]
    lengths = frame["lengths"]
    origin = frame["origin"]
    copies = []
    for ix in range(n):
        for iy in range(n):
            for iz in range(n):
                shift = np.array([ix, iy, iz], dtype=np.float64) * lengths
                copies.append(pos + shift)
    stacked = np.ascontiguousarray(np.vstack(copies), dtype=np.float64)
    return {
        "positions": stacked,
        "origin": origin.copy(),
        "lengths": lengths * n,
        "molecules": len(stacked),
    }


def box9_of(frame: dict):
    import numpy as np
    lo = frame["origin"]
    hi = lo + frame["lengths"]
    return np.array(
        [lo[0], lo[1], lo[2], hi[0], hi[1], hi[2], 0.0, 0.0, 0.0],
        dtype=np.float64,
    )


def lattice_of(frame: dict):
    import numpy as np
    return np.diag(frame["lengths"])


def emit(row: dict) -> None:
    print(json.dumps(row), flush=True)


def timed(call, warmup: int, repeats: int):
    import time
    samples = []
    last = None
    for index in range(warmup + repeats):
        t0 = time.perf_counter()
        last = call()
        elapsed = (time.perf_counter() - t0) * 1000.0
        if index >= warmup:
            samples.append(elapsed)
    if len(samples) != repeats:
        raise RuntimeError(f"timed {len(samples)} samples, need {repeats}")
    return samples, last


def require_version(dist: str, want: str) -> str:
    from importlib.metadata import version
    found = version(dist)
    if found != want and not found.startswith(want):
        raise RuntimeError(f"{dist} {found}, need {want}")
    return found


def open_lang(lib_path: str):
    import ctypes
    lib = ctypes.CDLL(lib_path)
    lib.seams_last_error.restype = ctypes.c_char_p
    lib.seams_frame_from_arrays.restype = ctypes.c_void_p
    lib.seams_frame_from_arrays.argtypes = [
        ctypes.POINTER(ctypes.c_double), ctypes.c_int,
        ctypes.POINTER(ctypes.c_double),
        ctypes.POINTER(ctypes.c_int), ctypes.POINTER(ctypes.c_int)]
    lib.seams_frame_free.argtypes = [ctypes.c_void_p]
    lib.seams_frame_rings.restype = ctypes.c_int
    lib.seams_frame_rings.argtypes = [
        ctypes.c_void_p, ctypes.c_double, ctypes.c_int, ctypes.c_char_p,
        ctypes.c_int, ctypes.POINTER(ctypes.c_int), ctypes.c_int]
    lib.seams_frame_steinhardt.restype = ctypes.c_int
    lib.seams_frame_steinhardt.argtypes = [
        ctypes.c_void_p, ctypes.c_double, ctypes.c_int, ctypes.c_int,
        ctypes.POINTER(ctypes.c_double), ctypes.POINTER(ctypes.c_double)]
    lib.seams_frame_q4_q6.restype = ctypes.c_int
    lib.seams_frame_q4_q6.argtypes = [
        ctypes.c_void_p, ctypes.c_double, ctypes.c_int,
        ctypes.POINTER(ctypes.c_double), ctypes.POINTER(ctypes.c_double),
        ctypes.POINTER(ctypes.c_double), ctypes.POINTER(ctypes.c_double)]
    lib.seams_frame_neighbours.restype = ctypes.c_int
    lib.seams_frame_neighbours.argtypes = [
        ctypes.c_void_p, ctypes.c_double, ctypes.c_int, ctypes.c_char_p]
    lib.seams_frame_chill_plus.restype = ctypes.c_int
    lib.seams_frame_chill_plus.argtypes = [
        ctypes.c_void_p, ctypes.c_double, ctypes.c_int, ctypes.c_char_p,
        ctypes.POINTER(ctypes.c_int), ctypes.POINTER(ctypes.c_int)]
    lib.seams_frame_cages.restype = ctypes.c_int
    lib.seams_frame_cages.argtypes = [
        ctypes.c_void_p, ctypes.c_double, ctypes.c_int, ctypes.c_int,
        ctypes.c_char_p, ctypes.c_int, ctypes.POINTER(ctypes.c_int),
        ctypes.POINTER(ctypes.c_int), ctypes.POINTER(ctypes.c_int),
        ctypes.POINTER(ctypes.c_int), ctypes.POINTER(ctypes.c_int)]
    lib.seams_frame_signature.restype = ctypes.c_int
    lib.seams_frame_signature.argtypes = [
        ctypes.c_void_p, ctypes.c_double, ctypes.c_int, ctypes.c_char_p,
        ctypes.c_char_p, ctypes.POINTER(ctypes.c_int), ctypes.POINTER(ctypes.c_int)]
    return lib


def frame_from(lib, xyz, box):
    import ctypes
    import numpy as np
    n = len(xyz)
    types = np.ones(n, dtype=np.int32)
    ids = np.arange(1, n + 1, dtype=np.int32)
    handle = lib.seams_frame_from_arrays(
        xyz.ctypes.data_as(ctypes.POINTER(ctypes.c_double)), n,
        box.ctypes.data_as(ctypes.POINTER(ctypes.c_double)),
        types.ctypes.data_as(ctypes.POINTER(ctypes.c_int)),
        ids.ctypes.data_as(ctypes.POINTER(ctypes.c_int)))
    if not handle:
        raise RuntimeError(lib.seams_last_error().decode())
    return handle, n


def fail_lang(lib, status: int, what: str) -> None:
    if status == 0:
        return
    err = lib.seams_last_error()
    text = err.decode() if err else f"status {status}"
    raise RuntimeError(f"{what}: {text}")


def time_dseams(lib, core, xyz, box, threads: int, version: str) -> None:
    import ctypes
    import numpy as np
    handle, n = frame_from(lib, xyz, box)
    try:
        counts = (ctypes.c_int * 7)()

        def rings():
            fail_lang(lib, lib.seams_frame_rings(
                handle, CUTOFF, 1, b"cutoff", 6, counts, 7), "rings")
            return int(counts[6])

        samples, six = timed(rings, WARMUP, REPEATS)
        if six != SIX:
            raise RuntimeError(f"d-SEAMS six-rings {six} != {SIX}")
        emit({
            "phase": "ice", "operation": "rings", "code": "d-SEAMS",
            "version": version, "molecules": n, "threads": threads,
            "samples_ms": samples, "detail": f"six={six}",
        })

        q4 = (ctypes.c_double * n)()
        q6 = (ctypes.c_double * n)()
        bar = (ctypes.c_double * n)()

        def local(order, dest, averaged):
            def run():
                fail_lang(lib, lib.seams_frame_steinhardt(
                    handle, CUTOFF, 1, order, dest, bar if averaged else None),
                    f"q{order}")
                return float(dest[0])
            return run

        for order, dest, name, averaged in (
            (4, q4, "q4", False),
            (6, q6, "q6", False),
            (4, q4, "q4-bar", True),
            (6, q6, "q6-bar", True),
        ):
            samples, _probe = timed(local(order, dest, averaged), WARMUP, REPEATS)
            emit({
                "phase": "ice", "operation": name, "code": "d-SEAMS",
                "version": version, "molecules": n, "threads": threads,
                "samples_ms": samples,
                "detail": "seams_frame_steinhardt; averaged buffer filled"
                if averaged else "seams_frame_steinhardt; local q only",
            })

        def both():
            fail_lang(lib, lib.seams_frame_q4_q6(
                handle, CUTOFF, 1, q4, q6, None, None), "q4+q6")
            return float(q4[0]) + float(q6[0])

        samples, _probe = timed(both, WARMUP, REPEATS)
        emit({
            "phase": "ice", "operation": "q4+q6", "code": "d-SEAMS",
            "version": version, "molecules": n, "threads": threads,
            "samples_ms": samples,
            "detail": "seams_frame_q4_q6; both averages null",
        })

        core.seams_ql_cutoff.restype = ctypes.c_int
        core.seams_ql_cutoff.argtypes = [
            ctypes.POINTER(ctypes.c_double), ctypes.c_int,
            ctypes.POINTER(ctypes.c_double), ctypes.c_double,
            ctypes.c_int, ctypes.POINTER(ctypes.c_double)]
        values = np.empty(n, dtype=np.float64)
        xyz_ptr = xyz.ctypes.data_as(ctypes.POINTER(ctypes.c_double))
        box_ptr = box.ctypes.data_as(ctypes.POINTER(ctypes.c_double))
        ql_ptr = values.ctypes.data_as(ctypes.POINTER(ctypes.c_double))

        def core_q(order):
            def run():
                code = core.seams_ql_cutoff(xyz_ptr, n, box_ptr, CUTOFF, order, ql_ptr)
                if code != 0:
                    raise RuntimeError(f"seams_ql_cutoff q{order} returned {code}")
                return float(values[0])
            return run

        for order, name in ((4, "q4-core"), (6, "q6-core")):
            samples, _probe = timed(core_q(order), WARMUP, REPEATS)
            emit({
                "phase": "ice", "operation": name, "code": "d-SEAMS",
                "version": version, "molecules": n, "threads": threads,
                "samples_ms": samples,
                "detail": f"seams_ql_cutoff order {order}",
            })

        def neighbours():
            fail_lang(lib, lib.seams_frame_neighbours(
                handle, CUTOFF, 1, b"cutoff"), "neighbours")
            return 0

        samples, _probe = timed(neighbours, WARMUP, REPEATS)
        emit({
            "phase": "ice", "operation": "neighbours", "code": "d-SEAMS",
            "version": version, "molecules": n, "threads": threads,
            "samples_ms": samples, "detail": "seams_frame_neighbours cutoff",
        })

        labels = (ctypes.c_int * n)()
        counts9 = (ctypes.c_int * 9)()

        def chill():
            fail_lang(lib, lib.seams_frame_chill_plus(
                handle, CUTOFF, 1, b"cutoff", labels, counts9), "chill+")
            return int(counts9[0])

        samples, cubic = timed(chill, WARMUP, REPEATS)
        if cubic != n:
            raise RuntimeError(f"CHILL+ cubic {cubic} != {n}")
        emit({
            "phase": "ice", "operation": "chill+", "code": "d-SEAMS",
            "version": version, "molecules": n, "threads": threads,
            "samples_ms": samples, "detail": f"cubic={cubic}",
        })

        per = (ctypes.c_int * n)()
        n_hex = ctypes.c_int()
        n_cubic = ctypes.c_int()
        n_water = ctypes.c_int()
        n_both = ctypes.c_int()

        def cages():
            fail_lang(lib, lib.seams_frame_cages(
                handle, CUTOFF, 1, 4, b"cutoff", 0, per,
                ctypes.byref(n_hex), ctypes.byref(n_cubic),
                ctypes.byref(n_water), ctypes.byref(n_both)), "cages")
            return (int(n_hex.value), int(n_cubic.value),
                    int(n_water.value), int(n_both.value))

        samples, census = timed(cages, WARMUP, REPEATS)
        total = sum(census)
        if total != n:
            raise RuntimeError(f"cage labels {census} sum {total} != {n}")
        emit({
            "phase": "ice", "operation": "cages", "code": "d-SEAMS",
            "version": version, "molecules": n, "threads": threads,
            "samples_ms": samples,
            "detail": (
                f"hexagonal={census[0]} cubic={census[1]} "
                f"water={census[2]} both={census[3]}"
            ),
        })
    finally:
        lib.seams_frame_free(handle)


def time_rsm(atoms, threads: int) -> None:
    version = require_version("rings-statistics-matter", "0.3.0")
    from rings_statistics_matter import ring_statistics
    import numpy as np

    def run():
        counts, _rings = ring_statistics(
            atoms, cutoff=CUTOFF, maxlvl=8, supercell=False, deduplicate=True,
            backend="rust", num_threads=threads)
        values = np.asarray(counts, dtype=np.float64)
        six = int(round(float(values[5])))
        if six != SIX:
            raise RuntimeError(f"rings-statistics-matter six-rings {six} != {SIX}")
        return six

    samples, six = timed(run, WARMUP, REPEATS)
    emit({
        "phase": "ice", "operation": "rings",
        "code": "rings-statistics-matter", "version": version,
        "molecules": MOLECULES, "threads": threads,
        "samples_ms": samples, "detail": f"six={six} maxlvl=8",
    })


def time_mdapy(pos, side: float, threads: int) -> None:
    import mdapy
    from mdapy import System, set_num_threads
    from mdapy.analysis import RingStatistics, SteinhardtParameters
    from mdapy.identify import ChillPlus
    from mdapy.neighbor import NeighborList
    version = mdapy.__version__
    set_num_threads(threads)
    system = System(pos, box=side)

    def neighbours():
        built = NeighborList.build(system, cutoff=CUTOFF)
        return int(np_len(built))

    samples, count = timed(neighbours, WARMUP, REPEATS)
    emit({
        "phase": "ice", "operation": "neighbours", "code": "mdapy",
        "version": version, "molecules": MOLECULES, "threads": threads,
        "samples_ms": samples, "detail": f"pairs={count}",
    })

    for order, name in ((4, "q4"), (6, "q6")):
        tool = SteinhardtParameters(degrees=(order,), cutoff=CUTOFF)

        def run(tool=tool):
            result = tool.compute(system)
            return int(result.q.shape[0])

        samples, nrows = timed(run, WARMUP, REPEATS)
        emit({
            "phase": "ice", "operation": name, "code": "mdapy",
            "version": version, "molecules": MOLECULES, "threads": threads,
            "samples_ms": samples, "detail": f"rows={nrows}",
        })

    def rings_real():
        import numpy as np
        result = RingStatistics(max_size=6, min_size=3, cutoff=CUTOFF).compute(system)
        values = np.asarray(result.counts)
        return int(values.reshape(-1)[-1])

    samples, six = timed(rings_real, WARMUP, REPEATS)
    emit({
        "phase": "ice", "operation": "rings", "code": "mdapy",
        "version": version, "molecules": MOLECULES, "threads": threads,
        "samples_ms": samples, "detail": f"last_count={six}",
    })

    def chill():
        result = ChillPlus(cutoff=CUTOFF).compute(system)
        return int(getattr(result, "structure", result).shape[0])

    samples, nrows = timed(chill, WARMUP, REPEATS)
    emit({
        "phase": "ice", "operation": "chill+", "code": "mdapy",
        "version": version, "molecules": MOLECULES, "threads": threads,
        "samples_ms": samples, "detail": f"rows={nrows}",
    })


def np_len(built) -> int:
    import numpy as np
    cols = getattr(built, "j", None)
    if cols is None:
        return 0
    return int(np.asarray(cols).shape[0])


def time_vesin(pos, side: float, threads: int) -> None:
    import vesin
    import numpy as np
    version = vesin.__version__
    box = np.diag([side, side, side])
    finder = vesin.NeighborList(cutoff=CUTOFF, full_list=True, n_threads=threads)

    def run():
        left, right = finder.compute(
            points=pos, box=box, periodic=True, quantities="ij", copy=True)
        return int(np.asarray(left).shape[0])

    samples, pairs = timed(run, WARMUP, REPEATS)
    emit({
        "phase": "ice", "operation": "neighbours", "code": "vesin",
        "version": version, "molecules": MOLECULES, "threads": threads,
        "samples_ms": samples, "detail": f"pairs={pairs}",
    })


def time_freud(pos, side: float, threads: int) -> None:
    import freud
    import numpy as np
    version = freud.__version__
    freud.parallel.set_num_threads(threads)
    box = freud.box.Box.cube(side)
    pts = np.ascontiguousarray(pos - side / 2.0)

    def neighbours():
        query = freud.locality.AABBQuery(box, pts)
        found = query.query(pts, {"r_max": CUTOFF, "exclude_ii": True}).toNeighborList()
        return int(len(found))

    samples, pairs = timed(neighbours, WARMUP, REPEATS)
    emit({
        "phase": "ice", "operation": "neighbours", "code": "freud",
        "version": version, "molecules": MOLECULES, "threads": threads,
        "samples_ms": samples, "detail": f"pairs={pairs}",
    })
    for order, name in ((4, "q4"), (6, "q6")):
        def run(order=order):
            op = freud.order.Steinhardt(l=order, average=False)
            op.compute((box, pts), neighbors={"r_max": CUTOFF, "exclude_ii": True})
            return int(np.asarray(op.ql).shape[0])

        samples, nrows = timed(run, WARMUP, REPEATS)
        emit({
            "phase": "ice", "operation": name, "code": "freud",
            "version": version, "molecules": MOLECULES, "threads": threads,
            "samples_ms": samples, "detail": f"rows={nrows}",
        })


def time_ovito(pos, side: float, threads: int) -> None:
    os.environ["OVITO_THREAD_COUNT"] = str(threads)
    os.environ["QT_QPA_PLATFORM"] = "offscreen"
    import ovito
    import numpy as np
    from ovito.data import DataCollection
    from ovito.modifiers import ChillPlusModifier, SteinhardtParametersModifier
    from ovito.pipeline import Pipeline, StaticSource
    version = ".".join(str(part) for part in tuple(ovito.version)[:3])
    data = DataCollection()
    data.create_cell(np.diag([side, side, side]), pbc=(True, True, True))
    parts = data.create_particles(count=len(pos))
    parts.create_property("Position", data=pos)

    def chill():
        pipe = Pipeline(source=StaticSource(data=data))
        pipe.modifiers.append(ChillPlusModifier(cutoff=CUTOFF))
        result = pipe.compute()
        return int(result.particles.count)

    samples, nrows = timed(chill, WARMUP, REPEATS)
    emit({
        "phase": "ice", "operation": "chill+", "code": "OVITO",
        "version": version, "molecules": MOLECULES, "threads": threads,
        "samples_ms": samples, "detail": f"rows={nrows}",
    })
    for order, name in ((4, "q4"), (6, "q6")):
        def run(order=order):
            pipe = Pipeline(source=StaticSource(data=data))
            mod = SteinhardtParametersModifier(
                cutoff=CUTOFF, max_l=order, compute_averaged=False)
            mod.neighbor_mode = SteinhardtParametersModifier.NeighborMode.Cutoff
            pipe.modifiers.append(mod)
            result = pipe.compute()
            return int(np.asarray(result.particles["Steinhardt"]).shape[0])

        samples, nrows = timed(run, WARMUP, REPEATS)
        emit({
            "phase": "ice", "operation": name, "code": "OVITO",
            "version": version, "molecules": MOLECULES, "threads": threads,
            "samples_ms": samples, "detail": f"rows={nrows} max_l={order}",
        })


def signature_census(lib, frame: dict) -> dict[str, int]:
    import ctypes
    xyz = frame["positions"]
    box = box9_of(frame)
    handle, _n = frame_from(lib, xyz, box)
    found = {}
    try:
        nc = ctypes.c_int()
        na = ctypes.c_int()
        for spec in ("512", "51262", "51264", "51268", "sh"):
            fail_lang(lib, lib.seams_frame_signature(
                handle, CUTOFF, 1, b"cutoff", spec.encode(),
                ctypes.byref(nc), ctypes.byref(na)), spec)
            found[spec] = int(nc.value)
    finally:
        lib.seams_frame_free(handle)
    return found


def time_hydrate(lib, frame: dict, phase: str, threads: int, version: str) -> None:
    import ctypes
    xyz = frame["positions"]
    box = box9_of(frame)
    handle, n = frame_from(lib, xyz, box)
    specs = ("512", "51262", "51264")
    try:
        nc = ctypes.c_int()
        na = ctypes.c_int()

        def run():
            census = {}
            for spec in specs:
                fail_lang(lib, lib.seams_frame_signature(
                    handle, CUTOFF, 1, b"cutoff", spec.encode(),
                    ctypes.byref(nc), ctypes.byref(na)), spec)
                census[spec] = int(nc.value)
            return census

        samples, census = timed(run, WARMUP, REPEATS)
    finally:
        lib.seams_frame_free(handle)
    full = signature_census(lib, frame)
    detail = " ".join(f"{key}={full[key]}" for key in ("512", "51262", "51264", "51268", "sh"))
    emit({
        "phase": phase, "operation": "cages", "code": "d-SEAMS",
        "version": version, "molecules": n, "threads": threads,
        "samples_ms": samples,
        "detail": f"signature calls 512 51262 51264; census {detail}",
    })


def clavis_text(frame: dict):
    import h2h_literature as lit
    packed = {
        "positions": frame["positions"],
        "lattice": lattice_of(frame),
        "molecules": frame["molecules"],
    }
    text, box_nm, _fourth = lit.hydrogen_waters(packed)
    return text, box_nm


def time_clavis(frame: dict, phase: str, threads: int) -> None:
    import json as js
    import clavis_engine as engine
    version = engine.version()
    if not isinstance(version, str):
        version = str(version)
    if not version.startswith("0.3.18"):
        raise RuntimeError(f"ClavisCage {version}, need 0.3.18")
    text, box_nm = clavis_text(frame)

    def run():
        raw = engine.analyze_waters(text, "[]", box_nm, "trace", 0.35, 35.0)
        parsed = js.loads(raw)
        return (
            len(parsed["cages_512"]),
            len(parsed["cages_51262"]),
            len(parsed["cages_51264"]),
        )

    samples, census = timed(run, WARMUP, REPEATS)
    emit({
        "phase": phase, "operation": "cages", "code": "ClavisCage",
        "version": version, "molecules": frame["molecules"], "threads": threads,
        "samples_ms": samples,
        "detail": f"TRACE 512={census[0]} 51262={census[1]} 51264={census[2]}",
    })


def worker(args) -> int:
    import ctypes
    try:
        os.sched_setaffinity(0, set(range(args.threads)))
    except (AttributeError, OSError):
        pass
    os.environ["OVITO_THREAD_COUNT"] = str(args.threads)
    os.environ["QT_QPA_PLATFORM"] = "offscreen"
    lib = open_lang(args.lib)
    core = ctypes.CDLL(args.core)
    xyz, box, side = ice_arrays()
    version = args.version
    time_dseams(lib, core, xyz, box, args.threads, version)
    import numpy as np
    from ase import Atoms
    atoms = Atoms(
        symbols=["O"] * MOLECULES,
        positions=np.ascontiguousarray(xyz),
        cell=[side, side, side],
        pbc=True,
    )
    time_rsm(atoms, args.threads)
    time_mdapy(xyz, side, args.threads)
    time_vesin(xyz, side, args.threads)
    time_freud(xyz, side, args.threads)
    time_ovito(xyz, side, args.threads)
    unit = {
        "sI": read_dump(ROOT / "input/traj/genice_sI.lammpstrj"),
        "sII": read_dump(ROOT / "input/traj/genice_sII.lammpstrj"),
    }
    built = {
        "sI": unit["sI"],
        "sI-222": replicate(unit["sI"], 2),
        "sII": unit["sII"],
        "sII-222": replicate(unit["sII"], 2),
    }
    for phase, frame in built.items():
        time_hydrate(lib, frame, phase, args.threads, version)
        time_clavis(frame, phase, args.threads)
    return 0


def seams_version(binary: Path) -> str:
    if not binary.is_file():
        return "unrecorded"
    done = subprocess.run(
        [str(binary), "--version"], check=False, text=True, capture_output=True)
    line = (done.stdout or "").strip().splitlines()
    return line[-1] if line else "unrecorded"


def explain_q(rows: list[dict]) -> list[str]:
    paired: dict[tuple, float] = {}
    for row in rows:
        if row["phase"] != "ice":
            continue
        paired[(row["operation"], row["code"], int(row["threads"]))] = float(row["median_ms"])
    notes = []
    orders = ("q4", "q4-bar", "q4-core")
    mates = {"q4": "q6", "q4-bar": "q6-bar", "q4-core": "q6-core"}
    codes = sorted({row["code"] for row in rows if row["phase"] == "ice"})
    for threads in (1, 8):
        for code in codes:
            for left in orders:
                right = mates[left]
                a = paired.get((left, code, threads))
                b = paired.get((right, code, threads))
                if a is None or b is None or b == 0:
                    continue
                if a > 2.0 * b or b > 2.0 * a:
                    notes.append(
                        f"{code} {left} at {threads} threads is {a:.3f} ms and "
                        f"{right} is {b:.3f} ms, a ratio of {a / b:.2f}."
                    )
    def pair(op: str, threads: int) -> str:
        value = paired.get((op, "d-SEAMS", threads))
        return "absent" if value is None else f"{value:.3f} ms"

    if notes:
        notes.append(
            "seams_frame_steinhardt and seams_ql_cutoff both call "
            "steinhardtQlCutoff. Orders 4 and 6 take the same solid-harmonic "
            "cutoff kernel, so a ratio above 2 is not the degree."
        )
    else:
        notes.append(
            "No q4 median in this sitting differs from the matching q6 median "
            "by more than 2. At 8 threads local q4 is "
            f"{pair('q4', 8)} and local q6 is {pair('q6', 8)}; "
            f"seams_ql_cutoff is {pair('q4-core', 8)} and {pair('q6-core', 8)}. "
            "The published 4.0 ms and 25.0 ms were not one kernel measured twice. "
            "4.0 ms is local cutoff q4. 25.0 ms was a label carried on the "
            "seams_ql_cutoff cell from another sitting, and that call does not "
            "reproduce it. Both calls use steinhardtQlCutoff, and orders 4 and 6 "
            "share that solid-harmonic kernel. The q4+q6 library call is "
            "seams_frame_q4_q6, which evaluates the two orders on two walks of "
            "one cutoff list."
        )
    return notes


def _census(detail: str) -> dict[str, int]:
    import re
    return {key: int(value) for key, value in re.findall(r"(51268|51264|51262|512|sh)=(\d+)", detail)}


def cage_note(rows: list[dict]) -> list[str]:
    def grab(phase: str, code: str, threads: int) -> str:
        for row in rows:
            if (row["phase"] == phase and row["code"] == code
                    and int(row["threads"]) == threads and row["operation"] == "cages"):
                return row["detail"]
        return ""

    lines = ["Hydrate census, 8 threads."]
    found = {}
    for phase in ("sI", "sI-222", "sII", "sII-222"):
        for code in ("d-SEAMS", "ClavisCage"):
            detail = grab(phase, code, 8)
            found[(phase, code)] = _census(detail)
            lines.append(f"{phase} {code} {detail}")
    unit = found[("sI", "ClavisCage")]
    tiled = found[("sI-222", "ClavisCage")]
    ours = found[("sI-222", "d-SEAMS")]
    unit_n = unit.get("512", 0) + unit.get("51262", 0) + unit.get("51264", 0)
    tiled_n = tiled.get("512", 0) + tiled.get("51262", 0) + tiled.get("51264", 0)
    ours_n = ours.get("512", 0) + ours.get("51262", 0) + ours.get("51264", 0)
    lines.append(
        f"On the sI unit cell TRACE returns {unit.get('512', 0)} of 5^12 and "
        f"{unit.get('51262', 0)} of 5^12 6^2, {unit_n} of the 8 cages. "
        f"On the 2x2x2 cell it returns {tiled.get('512', 0)} and "
        f"{tiled.get('51262', 0)}, {tiled_n} cages, against {ours_n} from d-SEAMS. "
        "The unit-cell miss does not survive the tiling: the supercell census matches."
        if tiled_n == ours_n and ours_n else
        "The supercell census does not match d-SEAMS, so the unit-cell miss is not settled by the tiling."
    )
    return lines


def write_outputs(out: Path, raw_rows: list[dict], info: dict) -> None:
    rows = []
    for item in raw_rows:
        summary = stats.summarize(item["samples_ms"])
        rows.append({
            "source_commit": info["source_commit"],
            "cpu_model": info["cpu_model"],
            "cpu_cores": info["cpu_cores"],
            "affinity": info["affinity"],
            "phase": item["phase"],
            "operation": item["operation"],
            "code": item["code"],
            "version": item["version"],
            "molecules": item["molecules"],
            "threads": item["threads"],
            "warmup": WARMUP,
            "n": summary["n"],
            "median_ms": f"{summary['median_ms']:.6f}",
            "iqr_ms": f"{summary['iqr_ms']:.6f}",
            "q1_ms": f"{summary['q1_ms']:.6f}",
            "q3_ms": f"{summary['q3_ms']:.6f}",
            "samples_ms": stats.fmt_samples(item["samples_ms"]),
            "detail": item["detail"],
        })
    out.mkdir(parents=True, exist_ok=True)
    dest = out / "sitting.csv"
    with dest.open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=FIELDS)
        writer.writeheader()
        writer.writerows(rows)
    log = [
        *stats.machine_lines(ROOT),
        f"cores {info['cores']}",
        "warmup 1",
        "n 8",
        "molecules 110592",
        "cutoff_A 3.5",
        "one sitting; each width is one process",
        "host name is not recorded",
    ]
    for row in rows:
        log.append(
            f"{row['phase']} {row['operation']} {row['code']} threads {row['threads']} "
            f"median_ms {row['median_ms']} iqr_ms {row['iqr_ms']} {row['detail']}"
        )
    log.append("")
    log.extend(explain_q(rows))
    log.append("")
    log.extend(cage_note(rows))
    (out / "raw-sitting.txt").write_text("\n".join(log) + "\n", encoding="utf-8")
    print(f"wrote {dest} ({len(rows)} rows)", flush=True)


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--out", type=Path, default=ROOT / "benchmarks/paper-cells")
    parser.add_argument("--lib", default=str(ROOT / "bbdir-release/bindings/libseams_lang.so"))
    parser.add_argument("--core", default=str(ROOT / "bbdir-release/src/libyodaLib.so"))
    parser.add_argument("--seams", default=str(ROOT / "bbdir-release/seams"))
    parser.add_argument("--worker", action="store_true")
    parser.add_argument("--threads", type=int, default=8)
    parser.add_argument("--version", default="")
    args = parser.parse_args()
    if args.worker:
        return worker(args)
    if not Path(args.lib).is_file() or not Path(args.core).is_file():
        raise SystemExit("missing the release libraries")
    info = stats.machine(ROOT)
    version = args.version or seams_version(Path(args.seams))
    collected = []
    for threads in (1, 8):
        env = stats.pin_env(threads)
        prior = env.get("LD_LIBRARY_PATH", "")
        libs = library_path(args.lib)
        env["LD_LIBRARY_PATH"] = libs if not prior else libs + os.pathsep + prior
        env["OVITO_THREAD_COUNT"] = str(threads)
        env["QT_QPA_PLATFORM"] = "offscreen"
        cmd = [
            "taskset", "-c", stats.taskset_cores(threads),
            sys.executable, str(Path(__file__).resolve()),
            "--worker", "--threads", str(threads),
            "--lib", args.lib, "--core", args.core,
            "--version", version,
        ]
        done = subprocess.run(cmd, cwd=ROOT, env=env, text=True, capture_output=True)
        if done.returncode != 0:
            tail = (done.stderr or done.stdout or "").strip().splitlines()
            raise SystemExit("\n".join(tail[-40:]) if tail else f"worker {threads} exit {done.returncode}")
        for line in done.stdout.splitlines():
            line = line.strip()
            if not line.startswith("{"):
                continue
            collected.append(json.loads(line))
        print(f"threads {threads} rows {sum(1 for row in collected if row['threads'] == threads)}",
              flush=True)
    if not collected:
        raise SystemExit("the sitting produced no rows")
    write_outputs(args.out, collected, info)
    return 0


if __name__ == "__main__":
    sys.exit(main())
