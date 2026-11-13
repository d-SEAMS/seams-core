#!/usr/bin/env python3
"""Interleaved ABAB timings for cutoff, kNN, q6 and CHILL+.

Pin the process before launch (taskset -c 0 for one thread, taskset -c 0-7
for eight). Each d-SEAMS sample is followed at once by the same sample of
one other code. n is 20. A win is a strictly smaller median and a strictly
smaller minimum.
"""

from __future__ import annotations

import ctypes
import os
import sys
import time

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))

import numpy as np

from h2h_freud_ovito import (
    CUTOFF,
    KNN,
    SEAMS_CHILL_OVITO,
    SEAMS_GRAPH_CUTOFF,
    core_count,
    frame_system,
    ice_positions,
    library_path,
    load_seams,
    require,
    set_ovito_threads,
)


N = 20
SIZES = (8, 16, 24)


def stats(samples):
    arr = np.asarray(samples, dtype=np.float64)
    q1, med, q3 = np.percentile(arr, [25, 50, 75])
    return float(med), float(q3 - q1), float(arr.min())


def quiet() -> None:
    """Let the previous arm drop its workers before the next sample starts."""
    from PySide6.QtCore import QThreadPool

    pool = QThreadPool.globalInstance()
    pool.waitForDone(1000)


def time_once(fn) -> float:
    quiet()
    t0 = time.perf_counter_ns()
    fn()
    return (time.perf_counter_ns() - t0) / 1e6


def main():
    width = int(sys.argv[1]) if len(sys.argv) > 1 else 1
    os.environ.setdefault("QT_QPA_PLATFORM", "offscreen")
    os.environ["OVITO_THREAD_COUNT"] = str(width)
    os.environ["OMP_NUM_THREADS"] = str(width)
    os.environ.setdefault("OMP_PROC_BIND", "close")
    os.environ.setdefault("OMP_PLACES", "cores")
    os.environ.setdefault("OMP_WAIT_POLICY", "PASSIVE")

    import freud
    import mdapy
    import ovito
    from ovito.data import DataCollection, NearestNeighborFinder
    from ovito.modifiers import (
        ChillPlusModifier,
        CreateBondsModifier,
        SteinhardtParametersModifier,
    )
    from ovito.pipeline import Pipeline, StaticSource

    lib = load_seams(library_path())
    require(lib, lib.seams_set_num_threads(width), "threads")
    mdapy.set_num_threads(width)
    freud.parallel.set_num_threads(width)
    set_ovito_threads(width)
    print(f"cores {core_count()}", flush=True)
    print(f"pinned_threads {width}", flush=True)
    print(f"samples {N}", flush=True)
    print("order ABAB", flush=True)
    print(
        f"freud {freud.__version__} ovito {'.'.join(str(v) for v in ovito.version)} "
        f"mdapy {mdapy.__version__}",
        flush=True,
    )

    losses = []
    for ncell in SIZES:
        pos, lengths = ice_positions(ncell)
        n = int(pos.shape[0])
        cloud = frame_system(lib, pos, lengths)
        q6b = (ctypes.c_double * n)()
        b6b = (ctypes.c_double * n)()
        lab = (ctypes.c_int * n)()
        cnt = (ctypes.c_int * 9)()
        fbox = freud.box.Box(lengths[0], lengths[1], lengths[2])
        mdapy_system = mdapy.System(pos, box=lengths)

        data = DataCollection()
        data.create_cell(np.diag(lengths), pbc=(True, True, True))
        parts = data.create_particles(count=n)
        parts.create_property("Position", data=np.ascontiguousarray(pos))
        bonds_mod = CreateBondsModifier(cutoff=CUTOFF)
        bond_pipe = Pipeline(source=StaticSource(data=data))
        bond_pipe.modifiers.append(bonds_mod)
        stein = SteinhardtParametersModifier(
            neighbor_mode=SteinhardtParametersModifier.NeighborMode.Cutoff,
            cutoff=CUTOFF,
            max_l=6,
            compute_averaged=True,
        )
        stein_pipe = Pipeline(source=StaticSource(data=data))
        stein_pipe.modifiers.append(stein)
        chill = ChillPlusModifier(cutoff=CUTOFF)
        chill_pipe = Pipeline(source=StaticSource(data=data))
        chill_pipe.modifiers.append(chill)

        def d_cutoff():
            handle = ctypes.c_void_p()
            require(
                lib,
                lib.seams_neighbours_build(cloud, CUTOFF, 1, SEAMS_GRAPH_CUTOFF, ctypes.byref(handle)),
                "cutoff",
            )
            lib.seams_neighbours_destroy(handle)

        def d_knn():
            handle = ctypes.c_void_p()
            require(lib, lib.seams_neighbours_knn(cloud, KNN, CUTOFF, 1, 0, ctypes.byref(handle)), "knn")
            lib.seams_neighbours_destroy(handle)

        def d_q6():
            require(lib, lib.seams_steinhardt(cloud, CUTOFF, 1, 6, q6b, n, b6b, n), "q6")

        def d_chill():
            require(
                lib,
                lib.seams_chill_plus(cloud, CUTOFF, 1, SEAMS_GRAPH_CUTOFF, SEAMS_CHILL_OVITO, lab, cnt),
                "chill+",
            )

        def f_cutoff():
            found = freud.locality.AABBQuery(fbox, pos).query(
                pos, {"r_max": CUTOFF, "exclude_ii": True}
            ).toNeighborList()
            if int(found.neighbor_counts.shape[0]) != n:
                raise SystemExit("freud cutoff")

        def f_knn():
            found = freud.locality.AABBQuery(fbox, pos).query(
                pos,
                {"mode": "nearest", "num_neighbors": KNN, "r_max": CUTOFF + 1.5, "exclude_ii": True},
            ).toNeighborList()
            if int(found.neighbor_counts.shape[0]) != n:
                raise SystemExit("freud knn")

        def f_q6():
            op = freud.order.Steinhardt(l=6, average=True)
            op.compute((fbox, pos), neighbors={"r_max": CUTOFF, "exclude_ii": True})

        def o_cutoff():
            bonds_mod.notify_object_changed()
            bond_pipe.compute()

        def o_knn():
            found, _vecs = NearestNeighborFinder(KNN, data).find_all()
            if found.shape != (n, KNN):
                raise SystemExit("ovito knn")

        def o_q6():
            stein.notify_object_changed()
            stein_pipe.compute()

        def o_chill():
            chill.notify_object_changed()
            chill_pipe.compute()

        def m_cutoff():
            mdapy.NeighborList.build(mdapy_system, cutoff=CUTOFF)

        def m_knn():
            mdapy.NeighborList.build(mdapy_system, num_neighbors=KNN)

        def m_q6():
            mdapy.analysis.SteinhardtParameters((6,), cutoff=CUTOFF, average=False).compute(mdapy_system)
            mdapy.analysis.SteinhardtParameters((6,), cutoff=CUTOFF, average=True).compute(mdapy_system)

        def m_chill():
            mdapy.identify.ChillPlus(cutoff=CUTOFF).compute(mdapy_system)

        jobs = {
            "cutoff": (d_cutoff, {"freud": f_cutoff, "ovito": o_cutoff, "mdapy": m_cutoff}),
            "knn": (d_knn, {"freud": f_knn, "ovito": o_knn, "mdapy": m_knn}),
            "q6": (d_q6, {"freud": f_q6, "ovito": o_q6, "mdapy": m_q6}),
            "chill_plus": (d_chill, {"ovito": o_chill, "mdapy": m_chill}),
        }
        for ours, rivals in jobs.values():
            ours()
            for fn in rivals.values():
                fn()
        for op, (ours, rivals) in jobs.items():
            for rival, fn in rivals.items():
                a, b = [], []
                for _ in range(N):
                    a.append(time_once(ours))
                    b.append(time_once(fn))
                am, ai, amin = stats(a)
                bm, bi, bmin = stats(b)
                win = am < bm and amin < bmin
                print(
                    f"robust op {op} threads {width} molecules {n} "
                    f"dseams_median {am:.4f} dseams_iqr {ai:.4f} dseams_min {amin:.4f} "
                    f"{rival}_median {bm:.4f} {rival}_iqr {bi:.4f} {rival}_min {bmin:.4f} "
                    f"{'WIN' if win else 'LOSE'}",
                    flush=True,
                )
                if not win:
                    losses.append(f"{op} threads {width} molecules {n} vs {rival}")
        lib.seams_system_destroy(cloud)
    print(f"losses {len(losses)}", flush=True)
    for row in losses:
        print(f"loss {row}", flush=True)
    if losses:
        raise SystemExit(1)


if __name__ == "__main__":
    main()
