#!/usr/bin/env python3
"""One process: d-SEAMS, freud, and OVITO on the same ice Ih cells.

Eight samples of knn, q6, and read at 4096, 32768, and 110592 molecules.
d-SEAMS and freud run at one thread and at eight. OVITO reads
OVITO_THREAD_COUNT when it starts and has no runtime setter, so this
process sets that to 8 before importing it. The record is the core count.
The host name is not written.
"""

from __future__ import annotations

import ctypes
import math
import os
import statistics
import tempfile
import time
from pathlib import Path

import numpy as np

ROOT = Path(__file__).resolve().parents[1]
SAMPLES = 8
SIZES = (8, 16, 24)  # orthorhombic cells: 4096, 32768, 110592 molecules


def library_path() -> Path:
    env = os.environ.get("SEAMS_LIB", "")
    if env:
        return Path(env)
    for candidate in (
        ROOT / "bbdir-release" / "src" / "libseams.so",
        ROOT / "build" / "src" / "libseams.so",
    ):
        if candidate.is_file():
            return candidate
    raise SystemExit("libseams.so was not found; set SEAMS_LIB")


def load_seams(path: Path):
    lib = ctypes.CDLL(str(path))
    lib.seams_set_num_threads.argtypes = [ctypes.c_int]
    lib.seams_set_num_threads.restype = ctypes.c_int
    lib.seams_system_create.restype = ctypes.c_void_p
    lib.seams_system_destroy.argtypes = [ctypes.c_void_p]
    lib.seams_system_count.argtypes = [ctypes.c_void_p]
    lib.seams_system_count.restype = ctypes.c_int
    lib.seams_system_set_frame.argtypes = [
        ctypes.c_void_p,
        ctypes.POINTER(ctypes.c_double),
        ctypes.c_int,
        ctypes.POINTER(ctypes.c_double),
        ctypes.POINTER(ctypes.c_int),
        ctypes.POINTER(ctypes.c_int),
    ]
    lib.seams_system_set_frame.restype = ctypes.c_int
    lib.seams_system_read.argtypes = [
        ctypes.c_void_p, ctypes.c_char_p, ctypes.c_int, ctypes.c_int
    ]
    lib.seams_system_read.restype = ctypes.c_int
    lib.seams_neighbours_knn.argtypes = [
        ctypes.c_void_p, ctypes.c_int, ctypes.c_double, ctypes.c_int, ctypes.c_int,
        ctypes.POINTER(ctypes.c_void_p),
    ]
    lib.seams_neighbours_knn.restype = ctypes.c_int
    lib.seams_neighbours_destroy.argtypes = [ctypes.c_void_p]
    lib.seams_neighbours_count.argtypes = [ctypes.c_void_p]
    lib.seams_neighbours_count.restype = ctypes.c_int
    lib.seams_steinhardt.argtypes = [
        ctypes.c_void_p, ctypes.c_double, ctypes.c_int, ctypes.c_int,
        ctypes.POINTER(ctypes.c_double), ctypes.c_size_t,
        ctypes.POINTER(ctypes.c_double), ctypes.c_size_t,
    ]
    lib.seams_steinhardt.restype = ctypes.c_int
    lib.seams_last_error.restype = ctypes.c_char_p
    return lib


def ice_positions(ncell: int):
    bond = 2.76
    a = bond * math.sqrt(8.0 / 3.0)
    c_len = a * math.sqrt(8.0 / 3.0)
    b = a * math.sqrt(3.0)
    basis = (
        (1.0 / 3.0, 2.0 / 3.0, 0.0),
        (2.0 / 3.0, 1.0 / 3.0, 0.5),
        (1.0 / 3.0, 2.0 / 3.0, 3.0 / 8.0),
        (2.0 / 3.0, 1.0 / 3.0, 7.0 / 8.0),
    )
    lx, ly, lz = ncell * a, ncell * b, ncell * c_len
    coords = []
    for ix in range(ncell):
        for iy in range(ncell):
            for iz in range(ncell):
                for shift in (0.0, 0.5):
                    for site in basis:
                        coords.append((
                            (shift * a + site[0] * a - site[1] * a / 2.0 + ix * a) % lx,
                            (shift * b + site[1] * b / 2.0 + iy * b) % ly,
                            (site[2] + iz) * c_len,
                        ))
    pos = np.asarray(coords, dtype=np.float64)
    return pos, (lx, ly, lz)


def write_dump(path: Path, pos: np.ndarray, lengths) -> None:
    lx, ly, lz = lengths
    n = pos.shape[0]
    lines = [
        "ITEM: TIMESTEP",
        "0",
        "ITEM: NUMBER OF ATOMS",
        str(n),
        "ITEM: BOX BOUNDS pp pp pp",
        f"0 {lx}",
        f"0 {ly}",
        f"0 {lz}",
        "ITEM: ATOMS id type x y z",
    ]
    for i in range(n):
        lines.append(f"{i + 1} 1 {pos[i, 0]:.8f} {pos[i, 1]:.8f} {pos[i, 2]:.8f}")
    path.write_text("\n".join(lines) + "\n", encoding="ascii")


def samples_of(fn, n=SAMPLES) -> list[float]:
    fn()
    out = []
    for _ in range(n):
        t0 = time.perf_counter_ns()
        fn()
        out.append((time.perf_counter_ns() - t0) / 1e6)
    return out


def emit(lines: list[str], molecules: int, code: str, op: str, threads: int, got: list[float]) -> None:
    raw = " ".join(f"{v:.4f}" for v in got)
    lines.append(
        f"arm {code} op {op} threads {threads} molecules {molecules} "
        f"samples_ms {raw}"
    )
    lines.append(
        f"arm {code} op {op} threads {threads} molecules {molecules} "
        f"median_ms {statistics.median(got):.4f}"
    )


def run_session() -> list[str]:
    os.environ["OVITO_THREAD_COUNT"] = "8"
    import freud
    from ovito.data import (
        DataCollection,
        NearestNeighborFinder,
        Particles,
        SimulationCell,
    )
    from ovito.io import import_file
    from ovito.modifiers import SteinhardtParametersModifier
    from ovito.pipeline import Pipeline, StaticSource

    lib = load_seams(library_path())
    lines = [
        f"cores {os.cpu_count() or 1}",
        f"samples {SAMPLES}",
        "ovito_thread_count 8",
        "note freud has no lammpstrj reader; read is d-SEAMS and OVITO",
        "note OVITO_THREAD_COUNT is read when OVITO starts",
    ]
    for ncell in SIZES:
        pos, lengths = ice_positions(ncell)
        n = int(pos.shape[0])
        lx, ly, lz = lengths
        lines.append(f"molecules {n}")
        with tempfile.TemporaryDirectory(prefix="seams-h2h-") as folder:
            dump = Path(folder) / "ice.lammpstrj"
            write_dump(dump, pos, lengths)
            encoded = os.fsencode(dump)
            box9 = (ctypes.c_double * 9)(0, 0, 0, lx, ly, lz, 0, 0, 0)
            types = (ctypes.c_int * n)(*([1] * n))
            ids = (ctypes.c_int * n)(*range(1, n + 1))
            ptr = pos.ctypes.data_as(ctypes.POINTER(ctypes.c_double))
            sys = lib.seams_system_create()
            if lib.seams_system_set_frame(sys, ptr, n, box9, types, ids) != 0:
                raise SystemExit("d-SEAMS set_frame failed")
            q = (ctypes.c_double * n)()
            bar = (ctypes.c_double * n)()
            nb = ctypes.c_void_p()

            def dseams_threads(width: int) -> None:
                if lib.seams_set_num_threads(width) != 0:
                    raise SystemExit("seams_set_num_threads failed")

            def dseams_knn():
                if lib.seams_neighbours_knn(sys, 4, 3.5, 1, 1, ctypes.byref(nb)) != 0:
                    raise SystemExit("d-SEAMS knn failed")
                if lib.seams_neighbours_count(nb) != n:
                    raise SystemExit("d-SEAMS knn count")
                lib.seams_neighbours_destroy(nb)

            def dseams_q6():
                if lib.seams_steinhardt(sys, 3.5, 1, 6, q, n, bar, n) != 0:
                    raise SystemExit("d-SEAMS q6 failed")
                if not math.isfinite(q[0]):
                    raise SystemExit("d-SEAMS q6 value")

            def dseams_read():
                fresh = lib.seams_system_create()
                if lib.seams_system_read(fresh, encoded, 1, 1) != 0:
                    raise SystemExit("d-SEAMS read failed")
                if lib.seams_system_count(fresh) != n:
                    raise SystemExit("d-SEAMS read count")
                lib.seams_system_destroy(fresh)

            fbox = freud.box.Box(lx, ly, lz)

            def freud_knn():
                query = freud.locality.AABBQuery(fbox, pos)
                nlist = query.query(pos, {
                    "mode": "nearest",
                    "num_neighbors": 4,
                    "r_max": 5.0,
                    "exclude_ii": True,
                }).toNeighborList()
                if int(nlist.neighbor_counts[0]) < 1:
                    raise SystemExit("freud knn empty")

            def freud_q6():
                order = freud.order.Steinhardt(l=6, average=True)
                order.compute((fbox, pos), neighbors={"r_max": 3.5, "exclude_ii": True})
                if order.ql.shape[0] != n or order.particle_order.shape[0] != n:
                    raise SystemExit("freud q6 shape")

            cell = SimulationCell(pbc=(True, True, True))
            matrix = np.zeros((3, 4))
            matrix[0, 0], matrix[1, 1], matrix[2, 2] = lx, ly, lz
            cell[...] = matrix
            data = DataCollection()
            data.objects.append(cell)
            particles = Particles()
            particles.create_property("Position", data=pos)
            data.objects.append(particles)
            modifier = SteinhardtParametersModifier(
                neighbor_mode=SteinhardtParametersModifier.NeighborMode.Cutoff,
                cutoff=3.5,
                max_l=6,
                compute_averaged=True,
            )
            pipe = Pipeline(source=StaticSource(data=data))
            pipe.modifiers.append(modifier)

            def ovito_knn():
                idx, _vecs = NearestNeighborFinder(4, data).find_all()
                if idx.shape != (n, 4):
                    raise SystemExit("OVITO knn shape")

            def ovito_q6():
                modifier.notify_object_changed()
                out = pipe.compute()
                if out.particles["Steinhardt"].shape[0] != n:
                    raise SystemExit("OVITO q6 shape")
                if out.particles["Smooth Steinhardt"].shape[0] != n:
                    raise SystemExit("OVITO q6bar shape")

            def ovito_read():
                loaded = import_file(os.fspath(dump)).compute()
                if loaded.particles.count != n:
                    raise SystemExit("OVITO read count")

            for width in (1, 8):
                dseams_threads(width)
                emit(lines, n, "dseams", "knn", width, samples_of(dseams_knn))
                emit(lines, n, "dseams", "q6", width, samples_of(dseams_q6))
                emit(lines, n, "dseams", "read", width, samples_of(dseams_read))
                freud.parallel.set_num_threads(width)
                emit(lines, n, "freud", "knn", width, samples_of(freud_knn))
                emit(lines, n, "freud", "q6", width, samples_of(freud_q6))
            emit(lines, n, "ovito", "knn", 8, samples_of(ovito_knn))
            emit(lines, n, "ovito", "q6", 8, samples_of(ovito_q6))
            emit(lines, n, "ovito", "read", 8, samples_of(ovito_read))
            lib.seams_system_destroy(sys)
    return lines


def main() -> None:
    for line in run_session():
        print(line)


if __name__ == "__main__":
    main()
