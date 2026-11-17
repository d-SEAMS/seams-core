#!/usr/bin/env python3
"""Equal-hardware head-to-head: d-SEAMS 2.14.0 (seams-core, the in-tree
Python front and the seams CLI) against OVITO 3.17.2, freud 3.5.0,
R.I.N.G.S. 1.3.5 and GenIce2 2.2.13.3 on one pinned core.

    h2h.py prep | chill | q6 | cna | rings | gring | cage | t512 | env

Every stage prints key=value lines; run.sh tees each to raw/<stage>.txt.
Timers: one warm-up then best of REPS; prep work runs outside the timer.
"""
from __future__ import annotations

import sys

sys.dont_write_bytecode = True

import json
import os
import pathlib
import shutil
import statistics
import subprocess
import time
from collections import Counter

import numpy as np

# Every directory is relative to this repository, or set by an environment
# variable. Nothing is pinned to a temporary directory or a home directory.
_ROOT = pathlib.Path(__file__).resolve().parents[1]
_H2H = pathlib.Path(os.environ.get("H2H_ROOT", str(_ROOT / "h2h")))
H = _H2H
RAW = pathlib.Path(os.environ.get("H2H_RAW", str(_H2H / "raw")))
DATA = pathlib.Path(os.environ.get("H2H_DATA", str(_H2H / "data")))
WORK = pathlib.Path(os.environ.get("H2H_WORK", str(_H2H / "work")))
SEAMS = os.environ.get("H2H_SEAMS", str(_ROOT / "build" / "seams"))
RINGS = os.environ.get("H2H_RINGS", str(_H2H / "rings" / "bin" / "rings"))
GI = os.environ.get("H2H_GI", str(_H2H / "env" / "gi"))
GENICE = os.environ.get("H2H_GENICE", str(pathlib.Path(GI) / "bin" / "genice2"))
ANALICE = os.environ.get("H2H_ANALICE", str(pathlib.Path(GI) / "bin" / "analice2"))
TRAJ = pathlib.Path(os.environ.get("H2H_TRAJ", str(_ROOT / "input" / "traj")))
NIU = pathlib.Path(
    os.environ.get(
        "H2H_NIU",
        str(_ROOT / "data" / "public-ice" / "nucleation"
            / "Ice_nucleation_trajectory" / "Nucleation_trajectory-3"),
    )
)
REPS = int(os.environ.get("H2H_REPS", "7"))
CUT = 3.5

# scripts/sota_compare.py lives next to this file, or wherever H2H_SCRIPTS points.
sys.path.insert(0, os.environ.get("H2H_SCRIPTS", str(pathlib.Path(__file__).resolve().parent)))


def emit(**kv):
    print(" ".join(f"{k}={v}" for k, v in kv.items()), flush=True)


# ------------------------------------------------------------------ inputs
def read_lammps_frame(path, frame=1, otype=None):
    """Frame `frame` (1-based; -1 = last) of a LAMMPS dump, wrapped into
    [0, L) with the box origin at zero."""
    lines = pathlib.Path(path).read_text().splitlines()
    starts = [i for i, l in enumerate(lines) if l.startswith("ITEM: TIMESTEP")]
    s = starts[frame - 1] if frame > 0 else starts[frame]
    step = int(lines[s + 1])
    n = int(lines[s + 3])
    lo = np.array([float(lines[s + 5 + k].split()[0]) for k in range(3)])
    hi = np.array([float(lines[s + 5 + k].split()[1]) for k in range(3)])
    cols = lines[s + 8].split()[2:]
    ix = [cols.index(c) for c in ("x", "y", "z")]
    it = cols.index("type")
    arr = [lines[s + 9 + k].split() for k in range(n)]
    typ = np.array([int(r[it]) for r in arr])
    pos = np.array([[float(r[j]) for j in ix] for r in arr])
    if otype is not None:
        pos = pos[typ == otype]
    box = hi - lo
    return (pos - lo) % box, box, dict(step=step, nframes=len(starts), lo=lo)


def read_gro_oxygens(path):
    lines = pathlib.Path(path).read_text().splitlines()
    n = int(lines[1])
    pos = []
    for line in lines[2 : 2 + n]:
        if line[10:15].strip().startswith("O"):
            pos.append((float(line[20:28]), float(line[28:36]), float(line[36:44])))
    cell = [float(x) for x in lines[2 + n].split()]
    assert len(cell) == 3, cell
    box = np.asarray(cell) * 10.0
    return np.asarray(pos) * 10.0 % box, box


def write_lammpstrj(path, pos, box):
    with open(path, "w") as fh:
        fh.write(f"ITEM: TIMESTEP\n0\nITEM: NUMBER OF ATOMS\n{len(pos)}\n")
        fh.write("ITEM: BOX BOUNDS pp pp pp\n")
        for L in box:
            fh.write(f"0.0 {L:.6f}\n")
        fh.write("ITEM: ATOMS id type x y z\n")
        for i, (x, y, z) in enumerate(pos):
            fh.write(f"{i + 1} 1 {x:.6f} {y:.6f} {z:.6f}\n")


def write_xyz(path, pos, box):
    with open(path, "w") as fh:
        fh.write(f"{len(pos)}\nbox {box[0]:.6f} {box[1]:.6f} {box[2]:.6f}\n")
        for x, y, z in pos:
            fh.write(f"O {x:.6f} {y:.6f} {z:.6f}\n")


def write_gro(path, pos, box):
    with open(path, "w") as fh:
        fh.write(f"h2h oxygen-only\n{len(pos)}\n")
        for i, (x, y, z) in enumerate(pos):
            r = (i + 1) % 100000
            fh.write(f"{r:5d}{'SOL':<5s}{'O':>5s}{r:5d}{x / 10:8.3f}{y / 10:8.3f}{z / 10:8.3f}\n")
        fh.write(f"{box[0] / 10:10.5f}{box[1] / 10:10.5f}{box[2] / 10:10.5f}\n")


def load(name):
    """In-memory (pos, box) for a named system, from the files prep wrote."""
    if name == "mw":
        pos, box, _ = read_lammps_frame(TRAJ / "mW_cubic.lammpstrj", 1, 1)
        return pos, box
    if name == "niu":
        pos, box, _ = read_lammps_frame(DATA / "niu_final.lammpstrj", 1)
        return pos, box
    if name == "crn":
        return read_gro_oxygens(DATA / "CRN1.genice.gro")
    if name == "si":
        return read_gro_oxygens(DATA / "sI_2x2x2.genice.gro")
    raise KeyError(name)


def stage_prep():
    DATA.mkdir(parents=True, exist_ok=True)
    WORK.mkdir(parents=True, exist_ok=True)
    if not NIU.is_file():
        fetch = pathlib.Path(__file__).resolve().parent / "fetch_niu_trajectory.py"
        subprocess.run(
            [sys.executable, str(fetch), "--dest", str(NIU.parents[2])],
            check=True,
        )
    env = dict(os.environ, PYTHONPATH=GI)
    for kind, rep, out in (("CS1", (2, 2, 2), "sI_2x2x2"), ("CRN1", (1, 1, 1), "CRN1")):
        cmd = [GENICE, kind, "--rep", *map(str, rep), "--format", "gromacs", "--seed", "1"]
        with open(DATA / f"{out}.genice.gro", "w") as fh:
            subprocess.run(cmd, check=True, stdout=fh, stderr=subprocess.DEVNULL, env=env)
        emit(stage="prep", cmd="'" + " ".join(cmd) + "'", out=f"data/{out}.genice.gro")
    pos, box, meta = read_lammps_frame(NIU, -1)
    write_lammpstrj(DATA / "niu_final.lammpstrj", pos, box)
    emit(stage="prep", system="niu", src=NIU.name, nframes=meta["nframes"],
         final_step=meta["step"], n=len(pos), box=",".join(f"{b:.4f}" for b in box))
    for name in ("mw", "niu", "crn", "si"):
        pos, box = load(name)
        write_xyz(DATA / f"{name}.xyz", pos, box)
        write_gro(DATA / f"{name}.gro", pos, box)
        write_lammpstrj(DATA / f"{name}.lammpstrj", pos, box)
        d = pos[:, None, :] - pos[None, :, :]
        d -= box * np.round(d / box)
        r = np.sqrt((d * d).sum(-1))
        np.fill_diagonal(r, np.inf)
        emit(stage="prep", system=name, n=len(pos), box=",".join(f"{b:.4f}" for b in box),
             edges_3p0=int((r < 3.0).sum() // 2), edges_3p5=int((r < CUT).sum() // 2),
             pairs_3p4_3p6=int(((r > 3.4) & (r < 3.6)).sum() // 2),
             dmin=f"{r.min():.4f}")


# ------------------------------------------------------------------ timing
def bench(label, prep, run, reps=REPS):
    """One warm-up then `reps` timed calls of run(prep()); prep is untimed."""
    times = []
    out = None
    for i in range(reps + 1):
        arg = prep()
        t0 = time.perf_counter()
        out = run(arg)
        dt = time.perf_counter() - t0
        if i:
            times.append(dt * 1e3)
    emit(timer=label, best_ms=f"{min(times):.3f}", median_ms=f"{statistics.median(times):.3f}",
         reps=reps, warmup=1, all_ms=",".join(f"{t:.3f}" for t in times))
    return out


def bench_proc(label, cmd, cwd=None, env=None, reps=REPS, clean=None):
    """Whole-process wall time of `cmd` (fork+exec to exit)."""
    times = []
    r = None
    for i in range(reps + 1):
        if clean:
            clean()
        t0 = time.perf_counter()
        r = subprocess.run(cmd, cwd=cwd, env=env, capture_output=True, text=True)
        dt = time.perf_counter() - t0
        if i:
            times.append(dt * 1e3)
    emit(timer=label, kind="whole-process", best_ms=f"{min(times):.3f}",
         median_ms=f"{statistics.median(times):.3f}", reps=reps, warmup=1,
         exit=r.returncode, cmd="'" + " ".join(map(str, cmd)) + "'",
         all_ms=",".join(f"{t:.3f}" for t in times))
    return r


# ------------------------------------------------------------------ d-SEAMS
def cloud_of(pos, box):
    from sota_compare import make_cloud

    return make_cloud(pos, box)


# ------------------------------------------------------------------ R.I.N.G.S.
RINGS_OPTIONS = """#######################################
       R.I.N.G.S. options file        #
#######################################
 PBC             .true.
 Frac            .false.
 g(r)            .false.
 S(q)            .false.
 S(k)            .false.
 gfft(r)         .false.
 MSD             .false.
 atMSD           .false.
 Bonds           .false.
 Angles          .false.
 Chains          .false.
---- ! chain statistics options ! -----
 Species           0
 AAAA              .false.
 ABAB              .false.
 01                .false.
---------------------------------------
 Rings           .true.
---- ! ring statistics options ! -----
 Species           0
 ABAB              .false.
 Rings0            .false.
 Rings1            .false.
 Rings2            .false.
 Rings3            .false.
 Rings4            .false.
 Prim_Rings        .true.
 Str_Rings         .false.
 BarycRings        .false.
 Prop-1            .false.
 Prop-2            .false.
 Prop-3            .false.
 Prop-4            .false.
 Prop-5            .false.
---------------------------------------
 Vacuum          .false.
#######################################
         ! OUTPUT OPTIONS !           #
#######################################
 Evol            .false.
 Dxout           .false.
! OPTION !       ! DX output options ! 
 RadOut          .false.
 RingsOut        .false.
 DRngOut         .false.
 VoidsOut        .false.
 TetraOut        .false.
 TrajOut         .false.
! OPTION !       ! Output options !   
 Output        h2h.out
#######################################
"""


def rings_workdir(name, pos, box, half_depth=3):
    """R.I.N.G.S. input: one species O, orthorhombic box, O-O cutoff CUT,
    max ring size 2*half_depth, primitive rings only."""
    wd = WORK / f"rings_{name}"
    if wd.exists():
        shutil.rmtree(wd)
    (wd / "data").mkdir(parents=True)
    write_xyz(wd / "data" / "conf.xyz", pos, box)
    (wd / "options").write_text(RINGS_OPTIONS)
    (wd / "input").write_text(
        "#######################################\n"
        "       R.I.N.G.S. input file          #\n"
        "#######################################\n"
        f"h2h-{name}\n{len(pos)}\n1\nO\n1\n0\n"
        f"{box[0]:.6f} {box[1]:.6f} {box[2]:.6f}\n90.0 90.0 90.0\n"
        "1.0\nANI\nconf.xyz\n200\n500\n25.0\n0.125\n90\n20\n"
        f"{half_depth}\n10\n"
        "#######################################\n"
        f"O O {CUT:.3f}\nGrtot {CUT:.3f}\n"
        "#######################################\n"
    )
    return wd


def rings_clean(wd):
    def clean():
        for d in ("rstat", "tmp", "gr", "sq", "cstat", "bonds", "angles"):
            shutil.rmtree(wd / d, ignore_errors=True)
    return clean


def rings_census(wd):
    """Primitive-ring census (absolute counts) from rstat/evol-RINGS-5.dat,
    the 'Number' column after 'Total number of rings'; RINGS-res-5.dat
    only holds counts per node."""
    out = {}
    seen = False
    for line in (wd / "rstat" / "evol-RINGS-5.dat").read_text().splitlines():
        if "Total number of rings" in line:
            seen = True
            continue
        cols = line.split()
        if seen and len(cols) >= 2:
            try:
                out[int(float(cols[0]))] = int(cols[1])
            except ValueError:
                continue
    return out


def stage_rings_mw_warm():
    """R.I.N.G.S. on mW cubic (4096) as the CRN1 row: whole process, one
    warm-up, then one timed run, on the core rings_mw.sh pins. Its 5x5x5
    superlattice makes the O(N^2) distance matrix dominate each run, so
    there is one timed run instead of the best of REPS."""
    pos, box = load("mw")
    wd = rings_workdir("mw", pos, box)
    r = bench_proc("rings-mw-RINGS", [RINGS, "input"], cwd=wd, env=env1(), reps=1,
                   clean=rings_clean(wd))
    (wd / "stdout.txt").write_text(r.stdout)
    emit(row="rings-mw", RINGS_primitive_D6=census_str(rings_census(wd)))


def env1():
    return dict(os.environ, OMP_NUM_THREADS="1", OPENBLAS_NUM_THREADS="1",
                RAYON_NUM_THREADS="1", OVITO_THREAD_COUNT="1")


# ------------------------------------------------------------------ CHILL+
# d-SEAMS AtomStateType name -> common class; OVITO ChillPlusModifier.Type
DS_CLASS = {"cubic": "cubic", "reCubic": "cubic", "hexagonal": "hexagonal",
            "reHex": "hexagonal", "interfacial": "interfacial",
            "clathrate": "hydrate", "interClathrate": "interfacial-hydrate",
            "water": "other", "unknown": "other"}
OV_CLASS = {0: "other", 1: "hexagonal", 2: "cubic", 3: "interfacial", 4: "hydrate",
            5: "interfacial-hydrate"}
OV_NAME = {0: "OTHER", 1: "HEXAGONAL_ICE", 2: "CUBIC_ICE", 3: "INTERFACIAL_ICE",
           4: "HYDRATE", 5: "INTERFACIAL_HYDRATE"}


def ds_chill(cloud, graph="cutoff"):
    from pydseams import yoda

    if graph == "cutoff":
        nl = yoda.neighListO(CUT, cloud, 1)
    else:  # the CLI's waterGraph: mutual 4-NN, candidates within CUT + 1.5
        nl = yoda.kNearestNeighbourList(cloud, 4, CUT + 1.5, 1, True)
    yoda.getCorrelPlus(cloud, nl, False)
    yoda.getIceTypePlusNoPrint(cloud, nl, False)
    return [pt.iceType.name for pt in cloud.pts]


def ovito_data(pos, box):
    from ovito.data import DataCollection

    data = DataCollection()
    data.create_cell(matrix=[[box[0], 0, 0, 0], [0, box[1], 0, 0], [0, 0, box[2], 0]],
                     pbc=(True, True, True))
    p = data.create_particles(count=len(pos))
    p.create_property("Position", data=pos)
    return data


def ov_chill(data):
    from ovito.modifiers import ChillPlusModifier

    data.apply(ChillPlusModifier(cutoff=CUT))
    return np.asarray(data.particles["Structure Type"]).copy()


def hist(labels):
    return ",".join(f"{k}:{v}" for k, v in sorted(Counter(labels).items()))


def stage_chill():
    for name in ("mw", "niu"):
        pos, box = load(name)
        ds = bench(f"chill-{name}-dseams", lambda: cloud_of(pos, box), ds_chill)
        ov = bench(f"chill-{name}-ovito", lambda: ovito_data(pos, box), ov_chill)
        a = np.array([DS_CLASS[x] for x in ds])
        b = np.array([OV_CLASS[int(t)] for t in ov])
        agree = (a == b)
        emit(row=f"chill-{name}", n=len(pos), dseams_native=hist(ds),
             ovito_native=hist(OV_NAME[int(t)] for t in ov))
        emit(row=f"chill-{name}", dseams_mapped=hist(a), ovito_mapped=hist(b),
             same=int(agree.sum()), agree=f"{agree.mean():.6f}",
             confusion=hist(f"{x}>{y}" for x, y in zip(a[~agree], b[~agree])) or "none")
        cli_graph = ds_chill(cloud_of(pos, box), graph="knn")
        c = np.array([DS_CLASS[x] for x in cli_graph])
        emit(row=f"chill-{name}-cligraph", dseams_native=hist(cli_graph),
             agree_vs_ovito=f"{(c == b).mean():.6f}",
             agree_vs_dseams_cutoff=f"{(c == a).mean():.6f}")
        f = TRAJ / "mW_cubic.lammpstrj" if name == "mw" else DATA / "niu_final.lammpstrj"
        r = subprocess.run([SEAMS, "chill-plus", str(f), "--cutoff", "3.5", "--type", "1"],
                           capture_output=True, text=True, env=env1())
        emit(row=f"chill-{name}-cli", cmd=f"'seams chill-plus {f} --cutoff 3.5 --type 1'",
             stdout="'" + r.stdout.strip() + "'")
    # check (a): exampleTraj, oxygen type 2, both graphs
    pos, box, _ = read_lammps_frame(TRAJ / "exampleTraj.lammpstrj", 1, 2)
    for graph in ("cutoff", "knn"):
        lab = ds_chill(cloud_of(pos, box), graph=graph)
        emit(row="check-a", file="exampleTraj.lammpstrj", n=len(pos),
             graph=("neighListO(3.5)" if graph == "cutoff" else "kNN(4,mutual,cand=5.0)"),
             dseams_native=hist(lab))


# ------------------------------------------------------------------ q6
def import_freud():
    import freud

    freud.parallel.set_num_threads(1)
    return freud


def ds_q6(cloud, rows=None):
    from pydseams import yoda

    if rows is None:
        rows = yoda.kNearestNeighbourList(cloud, 4, CUT + 1.5, 1, True)
    r = yoda.steinhardtQl(cloud, rows, 6)
    return np.asarray(r.ql), np.asarray(r.qlBar), rows


def fr_prep(pos, box):
    freud = import_freud()
    return freud.box.Box(Lx=box[0], Ly=box[1], Lz=box[2]), pos - box / 2.0


def fr_q6(arg, nlist=None):
    """q6 and q6bar need two Steinhardt objects: with average=True freud's
    .ql attribute is the averaged value, not the local q6."""
    freud = import_freud()
    fbox, pts = arg
    nb = nlist if nlist is not None else {"num_neighbors": 4}
    op = freud.order.Steinhardt(l=6, average=False)
    op.compute((fbox, pts), neighbors=nb)
    opb = freud.order.Steinhardt(l=6, average=True)
    opb.compute((fbox, pts), neighbors=nb)
    return np.asarray(op.particle_order).copy(), np.asarray(opb.particle_order).copy()


def maxdiff(x, y):
    ok = np.isfinite(x) & np.isfinite(y)
    return float(np.max(np.abs(x[ok] - y[ok]))), int((~ok).sum())


def stage_q6():
    freud = import_freud()
    for name in ("mw", "niu"):
        pos, box = load(name)
        q_ds, qb_ds, rows = bench(f"q6-{name}-dseams", lambda: cloud_of(pos, box), ds_q6)
        q_fr, qb_fr = bench(f"q6-{name}-freud", lambda: fr_prep(pos, box), fr_q6)
        fbox, pts = fr_prep(pos, box)
        op1 = freud.order.Steinhardt(l=6, average=True)
        op1.compute((fbox, pts), neighbors={"num_neighbors": 4})
        d_ql_avgflag = float(np.max(np.abs(np.asarray(op1.ql) - q_fr)))
        d6, nan6 = maxdiff(q_ds, q_fr)
        d6b, nan6b = maxdiff(qb_ds, qb_fr)
        deg = Counter(len(r) - 1 for r in rows)
        # freud's directed 4-NN rows handed to the d-SEAMS kernel
        nl = freud.locality.AABBQuery(fbox, pts).query(
            pts, dict(num_neighbors=4, exclude_ii=True)).toNeighborList()
        frows = [[i + 1] for i in range(len(pos))]
        for qi, pj in zip(nl.query_point_indices, nl.point_indices):
            frows[int(qi)].append(int(pj) + 1)
        same_nbrs = sum(sorted(a[1:]) == sorted(b[1:]) for a, b in zip(rows, frows))
        emit(row=f"q6-{name}", freud_points_dtype=str(
            freud.locality.AABBQuery(fbox, pts).points.dtype))
        q_s, qb_s, _ = ds_q6(cloud_of(pos, box), rows=frows)
        s6, _ = maxdiff(q_s, q_fr)
        s6b, _ = maxdiff(qb_s, qb_fr)
        emit(row=f"q6-{name}", n=len(pos),
             dseams_graph="kNearestNeighbourList(k=4,mutual=True)",
             freud_graph="num_neighbors=4", dseams_degree=",".join(
                 f"{k}:{v}" for k, v in sorted(deg.items())),
             atoms_same_neighbour_set=same_nbrs,
             mean_q6_ds=f"{np.nanmean(q_ds):.6f}", mean_q6_fr=f"{q_fr.mean():.6f}",
             mean_q6bar_ds=f"{np.nanmean(qb_ds):.6f}", mean_q6bar_fr=f"{qb_fr.mean():.6f}")
        emit(row=f"q6-{name}", maxabs_q6=f"{d6:.3e}", maxabs_q6bar=f"{d6b:.3e}",
             nonfinite_q6=nan6, nonfinite_q6bar=nan6b,
             samegraph_maxabs_q6=f"{s6:.3e}", samegraph_maxabs_q6bar=f"{s6b:.3e}",
             freud_ql_attr_under_average_true_vs_q6=f"{d_ql_avgflag:.3e}")


def shared_nlist(freud, fbox, pts, rows):
    """d-SEAMS ID-keyed rows (leading self entry) as a freud NeighborList."""
    qi, pj = [], []
    for r in rows:
        for j in r[1:]:
            qi.append(r[0] - 1)
            pj.append(j - 1)
    qi = np.asarray(qi, dtype=np.uint32)
    pj = np.asarray(pj, dtype=np.uint32)
    vec = fbox.wrap(pts[pj] - pts[qi])
    return freud.locality.NeighborList.from_arrays(len(pts), len(pts), qi, pj, vec)


def stage_q6shared():
    """q6 and q6bar of both codes on one neighbour list: the mutual
    four-nearest list, built once outside both timers."""
    freud = import_freud()
    from pydseams import yoda

    for name in ("mw", "niu"):
        pos, box = load(name)
        cloud = cloud_of(pos, box)
        rows = yoda.kNearestNeighbourList(cloud, 4, CUT + 1.5, 1, True)
        fbox, pts = fr_prep(pos, box)
        nl = shared_nlist(freud, fbox, pts, rows)
        q_ds, qb_ds, _ = bench(f"q6shared-{name}-dseams", lambda: cloud_of(pos, box),
                               lambda c: ds_q6(c, rows=rows))
        q_fr, qb_fr = bench(f"q6shared-{name}-freud", lambda: fr_prep(pos, box),
                            lambda a: fr_q6(a, nlist=nl))
        d6, nan6 = maxdiff(q_ds, q_fr)
        d6b, nan6b = maxdiff(qb_ds, qb_fr)
        deg = Counter(len(r) - 1 for r in rows)
        emit(row=f"q6shared-{name}", n=len(pos),
             graph="kNearestNeighbourList(k=4,mutual=True,cand=5.0), both codes",
             bonds=int(len(nl.point_indices)), degree=",".join(
                 f"{k}:{v}" for k, v in sorted(deg.items())),
             maxabs_q6=f"{d6:.3e}", maxabs_q6bar=f"{d6b:.3e}",
             nonfinite_q6=nan6, nonfinite_q6bar=nan6b,
             mean_q6_ds=f"{np.nanmean(q_ds):.6f}", mean_q6_fr=f"{np.nanmean(q_fr):.6f}")


# ------------------------------------------------------------------ CNA
# OVITO IdentifyDiamondModifier is common-neighbour analysis of the second
# shell; CommonNeighborAnalysisModifier is adaptive CNA of the first shell.
DIAMOND_CORE = {1: "cubic", 4: "hexagonal"}
DIAMOND_SHELL = {1: "cubic", 2: "cubic", 3: "cubic", 4: "hexagonal", 5: "hexagonal",
                 6: "hexagonal"}
DIAMOND_NAME = {0: "OTHER", 1: "CUBIC_DIAMOND", 2: "CUBIC_DIAMOND_FIRST_NEIGHBOR",
                3: "CUBIC_DIAMOND_SECOND_NEIGHBOR", 4: "HEX_DIAMOND",
                5: "HEX_DIAMOND_FIRST_NEIGHBOR", 6: "HEX_DIAMOND_SECOND_NEIGHBOR"}
CNA_NAME = {0: "OTHER", 1: "FCC", 2: "HCP", 3: "BCC", 4: "ICO"}


def ov_diamond(data):
    from ovito.modifiers import IdentifyDiamondModifier

    data.apply(IdentifyDiamondModifier())
    return np.asarray(data.particles["Structure Type"]).copy()


def ov_cna(data):
    from ovito.modifiers import CommonNeighborAnalysisModifier

    data.apply(CommonNeighborAnalysisModifier())
    return np.asarray(data.particles["Structure Type"]).copy()


def ds_seeded(cloud):
    """seams cages --graph seeded: mutual four-nearest seeds, union completion."""
    from pydseams import yoda

    mutual, union = yoda.kNearestNeighbourPair(cloud, 4, CUT + 1.5, 1)
    strict = yoda.neighbourListByIndex(cloud, mutual)
    perm = yoda.neighbourListByIndex(cloud, union)
    six_s = [r for r in yoda.ringNetwork(strict, 6) if len(r) == 6]
    six_u = [r for r in yoda.ringNetwork(perm, 6) if len(r) == 6]
    hc, ddc = yoda.seededCageAffiliation(six_s, strict, six_u, perm, False)
    return ["mixed" if h and d else "hexagonal" if h else "cubic" if d else "water"
            for h, d in zip(hc, ddc)]


def stage_cna():
    import sota_compare as sc

    pos0, box0 = sc.cubic_diamond(4)
    rng = np.random.default_rng(sc.SEED + 300 + 7919 * 5)
    systems = [("mw",) + tuple(load("mw")), ("niu",) + tuple(load("niu")),
               ("t512", sc.jitter(pos0, box0, 0.30, rng), box0)]
    for name, pos, box in systems:
        ds = bench(f"cna-{name}-dseams-seeded", lambda: cloud_of(pos, box), ds_seeded)
        dia = bench(f"cna-{name}-ovito-identify-diamond", lambda: ovito_data(pos, box),
                    ov_diamond)
        cna = bench(f"cna-{name}-ovito-adaptive-cna", lambda: ovito_data(pos, box), ov_cna)
        a = np.array([x if x in ("cubic", "hexagonal") else "other" for x in ds])
        core = np.array([DIAMOND_CORE.get(int(t), "other") for t in dia])
        shell = np.array([DIAMOND_SHELL.get(int(t), "other") for t in dia])
        emit(row=f"cna-{name}", n=len(pos), dseams_native=hist(ds),
             identify_diamond_native=hist(DIAMOND_NAME[int(t)] for t in dia),
             adaptive_cna_native=hist(CNA_NAME[int(t)] for t in cna))
        emit(row=f"cna-{name}", agree_core=f"{(a == core).mean():.6f}",
             agree_core_and_shells=f"{(a == shell).mean():.6f}",
             same_core=int((a == core).sum()), same_core_and_shells=int((a == shell).sum()))


# ------------------------------------------------------------------ rings
def ds_rings(cloud, depth=6):
    from pydseams import yoda

    nl = yoda.neighListO(CUT, cloud, 1)
    idx = yoda.neighbourListByIndex(cloud, nl)
    return Counter(len(r) for r in yoda.ringNetwork(idx, depth))


def census_str(c, sizes=range(3, 7)):
    return ",".join(f"{s}:{int(c.get(s, 0))}" for s in sizes)


def fingerprint_census(stdout):
    first = stdout.splitlines()[0].split()
    k = first.index("rings")
    return {int(t.split(":")[0]): int(t.split(":")[1]) for t in first[k + 1:]
            if ":" in t and t.split(":")[0].isdigit()}


def stage_rings():
    for name in ("mw", "crn"):
        pos, box = load(name)
        c_ds = bench(f"rings-{name}-dseams", lambda: cloud_of(pos, box), ds_rings)
        r = bench_proc(f"rings-{name}-seams-fingerprint",
                       [SEAMS, "fingerprint", str(DATA / f"{name}.lammpstrj"),
                        "--cutoff", "3.5", "--type", "1"], env=env1())
        c_fp = fingerprint_census(r.stdout)
        emit(row=f"rings-{name}", n=len(pos), dseams_D6=census_str(c_ds),
             seams_fingerprint_D7=census_str(c_fp, range(3, 8)),
             fingerprint_stdout_head="'" + r.stdout.splitlines()[0] + "'")
        if name == "crn":
            wd = rings_workdir(name, pos, box)
            r = bench_proc(f"rings-{name}-RINGS", [RINGS, "input"], cwd=wd, env=env1(),
                           clean=rings_clean(wd))
            (wd / "stdout.txt").write_text(r.stdout)
            c_x = rings_census(wd)
            same = all(int(c_ds.get(s, 0)) == int(c_x.get(s, 0)) for s in range(3, 7))
            emit(row=f"rings-{name}", RINGS_primitive_D6=census_str(c_x),
                 identical=("yes" if same else "no"))


# ------------------------------------------------------------------ GenIce
def gi_path():
    if GI not in sys.path:
        sys.path.insert(0, GI)


def gi_graph(pos, box):
    gi_path()
    import networkx as nx
    import pairlist as pl

    rel = pos / box
    g = nx.Graph()
    g.add_nodes_from(range(len(pos)))
    g.add_edges_from(pl.pairs_iter(rel, CUT, np.diag(box), distance=False))
    return g, rel


def gi_rings(arg, maxsize=6):
    gi_path()
    from cycless.cycles import cycles_iter

    pos, box = arg
    g, rel = gi_graph(pos, box)
    return [tuple(int(x) for x in ring) for ring in cycles_iter(g, maxsize, pos=rel)]


def yaplot_census(text):
    c = Counter()
    for line in text.splitlines():
        if line.startswith("p "):
            c[int(line.split()[1])] += 1
    return c


def stage_gring():
    env = dict(env1(), PYTHONPATH=GI)
    for name in ("mw", "crn"):
        pos, box = load(name)
        c_ds = bench(f"gring-{name}-dseams", lambda: cloud_of(pos, box), ds_rings)
        rings = bench(f"gring-{name}-genice-cycless", lambda: (pos, box), gi_rings)
        c_gi = Counter(len(r) for r in rings)
        ds_sets = None
        from pydseams import yoda
        cl = cloud_of(pos, box)
        idx = yoda.neighbourListByIndex(cl, yoda.neighListO(CUT, cl, 1))
        ds_sets = {frozenset(r) for r in yoda.ringNetwork(idx, 6)}
        gi_sets = {frozenset(r) for r in rings}
        emit(row=f"gring-{name}", n=len(pos), dseams_D6=census_str(c_ds),
             genice_cycless_D6=census_str(c_gi),
             same_vertex_sets=("yes" if ds_sets == gi_sets else "no"),
             only_dseams=len(ds_sets - gi_sets), only_genice=len(gi_sets - ds_sets))


def gi_cages(arg):
    """genice2_cage Hook2 with cage[12,14:ring=5,6]: rings up to 6 by
    cycless, polyhedra up to 14 faces, every face a 5- or 6-ring."""
    gi_path()
    from cycless.polyhed import polyhedra_iter

    ringlist = gi_rings(arg, 6)
    cages = []
    for cage in polyhedra_iter(ringlist, 14):
        if len(cage) in (12, 14) and all(len(ringlist[i]) in (5, 6) for i in cage):
            cages.append(list(cage))
    return Counter(cage_label(c, ringlist) for c in cages)


def cage_label(cage, ringlist):
    c = Counter(len(ringlist[i]) for i in cage)
    return "".join(f"{s}^{c[s]}" for s in sorted(c))


def ds_cages(cloud):
    from pydseams import yoda

    nl = yoda.neighListO(CUT, cloud, 1)
    idx = yoda.neighbourListByIndex(cloud, nl)
    rings = yoda.ringNetwork(idx, 6)
    return {sig: len(yoda.findBySignature(rings, idx, sig)) for sig in ("512", "51262")}


def stage_cage():
    pos, box = load("si")
    c_ds = bench("cage-si-dseams", lambda: cloud_of(pos, box), ds_cages)
    c_gi = bench("cage-si-genice-cycless", lambda: (pos, box), gi_cages)
    from pydseams.frame import Frame

    fr = Frame.from_arrays(pos, list(box), cutoff=CUT)
    paper = {sig: len(fr.cages_by_signature(sig)) for sig in ("512", "51262")}
    # analice2 2.2.13.3 crashes in Stage2 (raw/analice2_crash.txt). The
    # plugin's documented route regenerates the sI lattice and orients its
    # hydrogens before Hook2, so it runs once, untimed, for its cage counts;
    # the timed GenIce2 side is gi_cages on the frame both codes read.
    env = dict(env1(), PYTHONPATH=GI)
    r = subprocess.run([GENICE, "CS1", "--rep", "2", "2", "2", "--seed", "1",
                        "-f", "cage[12,14:ring=5,6:geom]"], env=env, capture_output=True,
                       text=True, check=True)
    geom = json.loads(r.stdout)
    c_an = Counter(cage_label(c, geom["rings"]) for c in geom["cages"])
    cli = {}
    for sig in ("512", "51262"):
        rr = subprocess.run([SEAMS, "cages", str(DATA / "si.lammpstrj"), "--cutoff", "3.5",
                             "--type", "1", "--graph", "cutoff", "--signature", sig],
                            env=env1(), capture_output=True, text=True, check=True)
        cli[sig] = rr.stdout.strip()
    emit(row="cage-si", n=len(pos), dseams_512=c_ds["512"], dseams_51262=c_ds["51262"],
         paper_frame_cages_by_signature=f"512:{paper['512']},51262:{paper['51262']}",
         genice_cycless=hist(f"{k}={v}" for k, v in c_gi.items()),
         genice2_cage_cli=hist(f"{k}={v}" for k, v in c_an.items()))
    for sig, out in cli.items():
        emit(row="cage-si-cli", signature=sig, stdout=f"'{out}'")


# ------------------------------------------------------------------ t512
def stage_t512():
    import tempfile

    import sota_compare as sc

    pos0, box = sc.cubic_diamond(4)
    sc.neighbour_check(pos0, box)
    rng = np.random.default_rng(sc.SEED + 300 + 7919 * 5)
    pos = sc.jitter(pos0, box, 0.30, rng)
    WORK.mkdir(parents=True, exist_ok=True)
    scratch = tempfile.mkdtemp(prefix="h2h-chill-", dir=str(WORK)) + os.sep
    for method, run in (
        ("dseams-seeded-adj", lambda _: sc.dseams_topo_seeded(pos, box, ring_adjacent=True)),
        ("ovito-ptm", lambda _: sc.ovito_ptm(pos, box)),
        ("dseams-chill+", lambda _: sc.chill_plus(pos, box, scratch)),
    ):
        labels = bench(f"t512-{method}", lambda: None, run)
        s = sc.score(labels, "cubic")
        emit(row="t512", method=method, n=len(pos), sigma="0.30", seed=5,
             rng_seed=sc.SEED + 300 + 7919 * 5, acc_cubic=f"{s['acc']:.4f}",
             hex=f"{s['hex']:.4f}", crystal=f"{s['crystal']:.4f}", labels=hist(labels))


# ------------------------------------------------------------------ env
def stage_env():
    import importlib.metadata as md
    import platform

    import freud
    import ovito
    from pydseams import yoda  # noqa: F401

    emit(python=platform.python_version(), numpy=np.__version__, numpy_path=np.__file__)
    emit(pydseamslib=md.version("pydseamslib"), ovito=ovito.version_string,
         freud=freud.__version__)
    gi_path()
    for pkg in ("genice2", "genice2-cage", "cycless", "pairlist", "networkx"):
        try:
            emit(package=pkg, version=md.version(pkg))
        except md.PackageNotFoundError:
            emit(package=pkg, version="not-found")


if __name__ == "__main__":
    stage = sys.argv[1]
    globals()["stage_" + stage]()
