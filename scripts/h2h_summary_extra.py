#!/usr/bin/env python3
"""Competitors named by the paper, in the summary sitting.

Each timed row asserts the outputs match. An operation the library does
not provide, or a size one sample cannot finish, is recorded as absent.
"""

from __future__ import annotations

import ctypes
import os

import numpy as np

PLUMED_KERNEL = (
    "/cvmfs/software.eessi.io/versions/2026.06/software/linux/x86_64/intel/"
    "sapphirerapids/software/PLUMED/2.10.0-foss-2026.1/lib/libplumedKernel.so"
)
Q6_SCALE = (4.0 * 3.14159265358979323846 / 13.0) ** 0.5


def absent(row_of, competitor, operation, molecules, threads, cores, model, version, why):
    made = row_of(competitor, operation, molecules, threads, cores, [1.0], [1.0], why, version, model)
    made["n"] = 0
    made["dseams_median_ms"] = ""
    made["dseams_min_ms"] = ""
    made["other_median_ms"] = ""
    made["other_min_ms"] = ""
    made["ratio"] = ""
    made["verdict"] = "absent"
    return made


def bind(lib):
    if not hasattr(lib, "seams_system_symmetry") or lib.seams_system_symmetry.argtypes is None:
        lib.seams_system_symmetry.argtypes = [
            ctypes.c_void_p, ctypes.c_double, ctypes.POINTER(ctypes.c_int), ctypes.c_int,
        ]
        lib.seams_system_symmetry.restype = ctypes.c_int
    lib.seams_rings_find.argtypes = [
        ctypes.c_void_p, ctypes.c_double, ctypes.c_int, ctypes.c_int, ctypes.c_int,
        ctypes.POINTER(ctypes.c_void_p),
    ]
    lib.seams_rings_find.restype = ctypes.c_int
    lib.seams_rings_destroy.argtypes = [ctypes.c_void_p]
    lib.seams_rings_histogram.argtypes = [
        ctypes.c_void_p, ctypes.POINTER(ctypes.c_int), ctypes.c_int,
    ]
    lib.seams_rings_histogram.restype = ctypes.c_int


def ring_hist(lib, held, h):
    rings = ctypes.c_void_p()
    h.require(lib, lib.seams_rings_find(held, h.CUTOFF, 1, 0, 6, ctypes.byref(rings)), "rings")
    hist = (ctypes.c_int * 7)()
    h.require(lib, lib.seams_rings_histogram(rings, hist, 7), "histogram")
    lib.seams_rings_destroy(rings)
    return [int(hist[k]) for k in range(7)]


def extra_rows(lib, held, pos, lengths, threads, cores, model, seams, pin, abab, row_of, h):
    n = int(pos.shape[0])
    rows = []
    version = seams

    # OVITO structure modifiers and rings. Adaptive CNA has no diamond
    # type, so it cannot match the ice labels.
    rows.append(absent(
        row_of, "OVITO", "cna", n, threads, cores, model, version,
        "adaptive CNA has no diamond structure and does not match the ice labels",
    ))
    rows.append(absent(
        row_of, "PTM", "ptm", n, threads, cores, model, version,
        "the standalone polyhedral-template-matching library is not linked",
    ))
    rows.append(absent(
        row_of, "ClavisCage", "chill_plus", n, threads, cores, model, version,
        "ClavisCage is not installed; the clavis package here is unrelated",
    ))
    rows.append(absent(
        row_of, "freud", "rings", n, threads, cores, model, version,
        "freud has no primitive-ring finder",
    ))
    rows.append(absent(
        row_of, "SQQ", "cages", n, threads, cores, model, "0.5.8",
        "SQQ counts hydrate cages and does not finish an ice frame of this size",
    ))
    if n > 4096:
        rows.append(absent(
            row_of, "cycless", "rings", n, threads, cores, model, "0.7",
            "cycles_iter takes 12 s on 4096 molecules; this size is not in the sitting",
        ))
        rows.append(absent(
            row_of, "IRA", "shape", n, threads, cores, model, "2.2.0",
            "one IRA match is 0.37 ms per atom; this frame is outside the sitting",
        ))

    from ovito.data import DataCollection
    from ovito.modifiers import (
        CreateBondsModifier,
        FindRingsModifier,
        IdentifyDiamondModifier,
        PolyhedralTemplateMatchingModifier,
    )
    from ovito.pipeline import Pipeline, StaticSource
    from PySide6.QtCore import QThreadPool

    QThreadPool.globalInstance().setMaxThreadCount(threads)
    data = DataCollection()
    data.create_cell(np.diag(lengths), pbc=(True, True, True))
    parts = data.create_particles(count=n)
    parts.create_property("Position", data=np.ascontiguousarray(pos))

    sym = (ctypes.c_int * n)()
    h.require(lib, lib.seams_set_num_threads(threads), "threads")
    h.require(lib, lib.seams_system_symmetry(held, 0.0, sym, n), "symmetry")
    sym_np = np.ctypeslib.as_array(sym).copy()

    def fresh_pipe(*mods):
        pipe = Pipeline(source=StaticSource(data=data))
        for mod in mods:
            pipe.modifiers.append(mod)
        return pipe

    # A serial PTM or ring sample on 884736 molecules does not return
    # inside this sitting. The smaller cells are timed.
    if n > 110592:
        rows.append(absent(
            row_of, "OVITO", "rings", n, threads, cores, model, version,
            "one FindRings or PTM sample at 884736 molecules does not finish in the sitting",
        ))
        rows.append(absent(
            row_of, "OVITO", "ptm", n, threads, cores, model, version,
            "one FindRings or PTM sample at 884736 molecules does not finish in the sitting",
        ))
        bonds = None
        ring_pipe = None
        ptm_pipe = None
        ours_hist = ring_hist(lib, held, h)
        ovito_hex = ours_hist[6]
    else:
        bonds = CreateBondsModifier(cutoff=h.CUTOFF)
        rings_mod = FindRingsModifier(maximum_ring_size=6)
        ring_pipe = fresh_pipe(bonds, rings_mod)
        table = ring_pipe.compute().tables["ring-size-histogram"]
        ovito_hex = 0
        for size, count in zip(table["Ring Size"][...], table["Count"][...]):
            if int(size) == 6:
                ovito_hex = int(count)
        ours_hist = ring_hist(lib, held, h)
        if ours_hist[6] != ovito_hex:
            raise SystemExit(f"OVITO rings {ovito_hex} differ from d-SEAMS {ours_hist[6]} at {n}")
        ptm = PolyhedralTemplateMatchingModifier(rmsd_cutoff=0.1)
        ptm.structures[PolyhedralTemplateMatchingModifier.Type.CUBIC_DIAMOND].enabled = True
        ptm.structures[PolyhedralTemplateMatchingModifier.Type.HEX_DIAMOND].enabled = True
        ptm_pipe = fresh_pipe(ptm)
        ptm_lab = np.asarray(ptm_pipe.compute().particles["Structure Type"])
        # d-SEAMS symmetry 0 is hexagonal. OVITO PTM 7 is hexagonal diamond.
        if not np.all(sym_np == 0) or not np.all(ptm_lab == 7):
            raise SystemExit(f"PTM labels do not match hexagonal ice at {n}")

    dia = IdentifyDiamondModifier()
    dia_pipe = fresh_pipe(dia)
    dia_lab = np.asarray(dia_pipe.compute().particles["Structure Type"])
    # Identify diamond: 4 hexagonal, 5 first shell, 6 second shell.
    if not np.all(np.isin(dia_lab, (4, 5, 6))):
        raise SystemExit(f"identify-diamond is not hexagonal at {n}")

    def dseams_rings():
        pin(threads)
        h.require(lib, lib.seams_set_num_threads(threads), "threads")
        ring_hist(lib, held, h)

    def ovito_rings():
        pin(threads)
        QThreadPool.globalInstance().setMaxThreadCount(threads)
        for mod in ring_pipe.modifiers:
            mod.notify_object_changed()
        ring_pipe.compute()

    def dseams_sym():
        pin(threads)
        h.require(lib, lib.seams_set_num_threads(threads), "threads")
        h.require(lib, lib.seams_system_symmetry(held, 0.0, sym, n), "symmetry")

    def ovito_ptm():
        pin(threads)
        QThreadPool.globalInstance().setMaxThreadCount(threads)
        ptm_pipe.modifiers[0].notify_object_changed()
        ptm_pipe.compute()

    def ovito_dia():
        pin(threads)
        QThreadPool.globalInstance().setMaxThreadCount(threads)
        dia_pipe.modifiers[0].notify_object_changed()
        dia_pipe.compute()

    import ovito
    ov = ".".join(str(v) for v in ovito.version)
    jobs = [("OVITO", "identify_diamond", dseams_sym, ovito_dia, f"hexagonal {n} equal")]
    if n <= 110592:
        jobs = [
            ("OVITO", "rings", dseams_rings, ovito_rings, f"hexagons {ovito_hex} equal"),
            ("OVITO", "ptm", dseams_sym, ovito_ptm, f"hexagonal {n} equal"),
        ] + jobs
    for name, op, ours, other, note in jobs:
        left, right = abab(ours, other)
        rows.append(row_of(
            name, op, n, threads, cores, left, right, note, f"{version}/{ov}", model,
        ))

    # PLUMED Q6, scaled to the local Steinhardt parameter. One serial
    # sample at 884736 molecules does not return inside the sitting.
    if n > 110592:
        rows.append(absent(
            row_of, "PLUMED", "q6", n, threads, cores, model, version,
            "one serial Q6 sample at 884736 molecules does not finish in the sitting",
        ))
    elif os.path.isfile(PLUMED_KERNEL):
        os.environ.setdefault("PLUMED_KERNEL", PLUMED_KERNEL)
        import plumed

        q = (ctypes.c_double * n)()
        bar = (ctypes.c_double * n)()
        h.dseams_q6(lib, held, q, bar)
        local = np.ctypeslib.as_array(q).copy()
        line = f"obs: Q6 SPECIES=1-{n} SWITCH={{CUSTOM FUNC=step(3.5-x) R_0=1 D_MAX=4.2}}"
        calc = plumed.Plumed()
        calc.cmd("setMDEngine", "python")
        calc.cmd("setNatoms", n)
        calc.cmd("setTimestep", 1.0)
        calc.cmd("setKbT", 1.0)
        calc.cmd("setLogFile", "/dev/null")
        calc.cmd("init")
        calc.cmd("readInputLine", line)
        rank = np.zeros(1, dtype=np.int64)
        calc.cmd("getDataRank obs", rank)
        shape = np.zeros(int(rank[0]), dtype=np.int64)
        calc.cmd("getDataShape obs", shape)
        data = np.zeros(tuple(int(v) for v in shape), dtype=np.float64)
        calc.cmd("setMemoryForData obs", data)
        box = np.diag(np.asarray(lengths, dtype=np.float64))
        mass = np.ones(n, dtype=np.float64)
        charge = np.zeros(n, dtype=np.float64)
        forces = np.zeros((n, 3), dtype=np.float64)
        virial = np.zeros((3, 3), dtype=np.float64)
        step = {"n": 0}

        def plumed_q():
            step["n"] += 1
            calc.cmd("setStep", step["n"])
            calc.cmd("setBox", box)
            calc.cmd("setPositions", pos)
            calc.cmd("setMasses", mass)
            calc.cmd("setCharges", charge)
            calc.cmd("setForces", forces)
            calc.cmd("setVirial", virial)
            calc.cmd("calc")

        plumed_q()
        got = np.array(data, copy=True).reshape(-1) * Q6_SCALE
        diff = float(np.max(np.abs(local - got)))
        if diff > 1e-8:
            raise SystemExit(f"PLUMED q6 max_abs {diff:.3e} at {n}")

        def dseams_local():
            pin(threads)
            h.require(lib, lib.seams_set_num_threads(threads), "threads")
            h.dseams_q6(lib, held, q, bar)

        left, right = abab(dseams_local, plumed_q)
        rows.append(row_of(
            "PLUMED", "q6", n, threads, cores, left, right,
            f"q6_max_abs {diff:.3e}", f"{version}/2.10.0", model,
        ))
    else:
        rows.append(absent(
            row_of, "PLUMED", "q6", n, threads, cores, model, version,
            "libplumedKernel.so was not found",
        ))

    # rings-statistics-matter, the Rust shortest-path census.
    import rings_statistics_matter
    from ase import Atoms

    atoms = Atoms(symbols=["O"] * n, positions=pos, cell=list(lengths), pbc=True)

    def rsm_counts():
        _counts, blocks = rings_statistics_matter.ring_statistics(
            atoms, cutoff=h.CUTOFF, maxlvl=8, supercell=False,
            deduplicate=True, backend="rust", num_threads=threads,
        )
        found = [0] * 7
        for block in blocks:
            arr = np.asarray(block)
            if arr.size == 0 or arr.ndim != 2:
                continue
            size = int(arr.shape[1])
            if 0 <= size <= 6:
                found[size] = int(arr.shape[0])
        return found

    rsm = rsm_counts()
    if rsm[6] != ours_hist[6]:
        raise SystemExit(f"rings-statistics-matter {rsm[6]} hexagons vs d-SEAMS {ours_hist[6]}")

    def dseams_r():
        pin(threads)
        h.require(lib, lib.seams_set_num_threads(threads), "threads")
        ring_hist(lib, held, h)

    def rsm_call():
        pin(threads)
        rsm_counts()

    left, right = abab(dseams_r, rsm_call)
    rows.append(row_of(
        "rings-statistics-matter", "rings", n, threads, cores, left, right,
        f"hexagons {rsm[6]} equal", f"{version}/0.3.0", model,
    ))

    if n <= 4096:
        import networkx as nx
        from cycless.cycles import cycles_iter
        import freud

        box = freud.box.Box(*lengths)
        query = freud.locality.AABBQuery(box, pos).query(
            pos, {"r_max": h.CUTOFF, "exclude_ii": True}
        ).toNeighborList()
        graph = nx.Graph()
        graph.add_nodes_from(range(n))
        ii = np.asarray(query[:, 0])
        jj = np.asarray(query[:, 1])
        graph.add_edges_from(
            (int(a), int(b)) for a, b in zip(ii, jj) if int(a) < int(b)
        )
        cyc = sum(1 for _ in cycles_iter(graph, 6))
        if cyc != ours_hist[6]:
            raise SystemExit(f"cycless {cyc} rings vs d-SEAMS {ours_hist[6]}")

        def cyc_call():
            pin(threads)
            sum(1 for _ in cycles_iter(graph, 6))

        left, right = abab(dseams_r, cyc_call)
        rows.append(row_of(
            "cycless", "rings", n, threads, cores, left, right,
            f"hexagons {cyc} equal; cycles_iter has no thread argument",
            f"{version}/0.7", model,
        ))

        import ira_mod

        ira = ira_mod.IRA()
        from h2h_ira import classify_serial, cluster

        ih = cluster(pos, np.asarray(lengths), 0)
        ic_pos, side = __import__("gen_ice_ic", fromlist=["ice_ic"]).ice_ic(4)
        ic = cluster(
            np.ascontiguousarray(ic_pos),
            np.asarray((side, side, side)),
            0,
        )
        labels_ira = classify_serial(ira, ih, ic, pos, np.asarray(lengths))
        # Ideal Ih is hexagonal under symmetry. IRA's nearer template is Ih.
        if int(np.count_nonzero(labels_ira == 0)) != n:
            raise SystemExit(f"IRA labelled {int(np.count_nonzero(labels_ira == 0))} of {n} hexagonal")

        def ira_call():
            pin(threads)
            classify_serial(ira, ih, ic, pos, np.asarray(lengths))

        left, right = abab(dseams_sym, ira_call)
        ver = ira.get_version()
        rows.append(row_of(
            "IRA", "shape", n, threads, cores, left, right,
            f"hexagonal {n} equal", f"{version}/{ver[0]}", model,
        ))
    return rows
