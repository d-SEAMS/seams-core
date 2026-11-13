#!/usr/bin/env python3
"""Held-out check of the thermal symmetry gates.

The gates 0.85, 1.20 and 0.15 were chosen on the TIP4P/2005 slab at
252 K. This run builds other systems with the LAMMPS decks already in
the tree and scores the old gates (0.94, 1.15, 0.12) and the new gates
against the same quench-majority and 4 ps-average majority. A reference
atom is crystal when at least two of the new-gate classifier, PTM
diamond and CHILL+ hexagonal-or-cubic agree. Frames stay under
data/holdout/ and are not committed.

Two LAMMPS jobs run at a time. lmp is taken from LMP, then LAMMPS_BIN,
then PATH.

On the runs that produced benchmarks/paper-cells/symmetry-holdout.csv
the new gates raise recall on every held-out slab and keep precision at
or above 0.99. The extra ice labels sit in the first neighbour shell of
majority crystal. Bulk TIP4P/2005 at 240 K with no seed is entirely
other, including after a minimise. The gates were left as they are.
No CHILL+-labelled brine coordinate frame is stored in the tree; the
brine row is the 230 K start built by the growth deck.
"""

from __future__ import annotations

import argparse
import csv
import importlib.util
import os
import random
import shutil
import subprocess
import sys
from concurrent.futures import ThreadPoolExecutor, as_completed
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "scripts"))

import numpy as np
import paper_ptm_templates as templates
from paper_symmetry_truth import average_positions, read_oxygen_frames


def _load(name: str, path: Path):
    spec = importlib.util.spec_from_file_location(name, path)
    module = importlib.util.module_from_spec(spec)
    sys.modules[name] = module
    spec.loader.exec_module(module)
    return module


brine = _load("regenerate", ROOT / "input" / "brine" / "regenerate.py")
tip4p = _load("regenerate_tip4p", ROOT / "input" / "hop" / "regenerate_tip4p.py")

HOLD = ROOT / "data" / "holdout"
OLD = (0.12, 0.94, 1.15)
NEW = (0.15, 0.85, 1.20)


def lmp_bin() -> str:
    for key in ("LMP", "LAMMPS_BIN"):
        env = os.environ.get(key)
        if env:
            return env
    found = shutil.which("lmp")
    if found:
        return found
    raise SystemExit("lmp not found; set LMP")


def lmp_env(threads: int) -> dict[str, str]:
    env = os.environ.copy()
    env["OMP_NUM_THREADS"] = str(threads)
    env.setdefault("UCX_TLS", "tcp,self,sm")
    env.setdefault("OMPI_MCA_btl", "^openib")
    env.setdefault("OMPI_MCA_pml", "ob1")
    return env


def run_lammps(work: Path, deck: str, threads: int = 4) -> None:
    cmd = [lmp_bin(), "-sf", "omp", "-pk", "omp", str(threads), "-nocite",
           "-in", deck, "-log", deck + ".log", "-screen", "none"]
    subprocess.run(cmd, cwd=work, env=lmp_env(threads), check=True)


def write(path: Path, text: str) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(text)


def prepare_tip4p(temp: float) -> None:
    work = HOLD / f"tip4p-{temp:.0f}"
    if (work / "quench.lammpstrj").is_file() and (work / "average.lammpstrj").is_file():
        return
    work.mkdir(parents=True, exist_ok=True)
    waters, box, n_ice, _ice_box, _height = tip4p.build_slab(2522005)
    tip4p.write_data(work / "data.tip4p", waters, box)
    steps = {"melt": 2500, "cool": 1000, "equil": 5000}
    write(work / "in.equil", f"""
variable temp index {temp:.2f}
units real
atom_style full
boundary p p p
atom_modify map array
read_data data.tip4p
bond_style harmonic
bond_coeff 1 450.0 0.9572
angle_style harmonic
angle_coeff 1 55.0 104.52
special_bonds lj/coul 0.0 0.0 0.0
reset_atoms image all
pair_style lj/cut/tip4p/long 1 2 1 1 0.1546 10.0
pair_modify mix arithmetic
kspace_style pppm/tip4p 1.0e-4
pair_coeff 1 1 0.1852 3.1589
pair_coeff 2 2 0.0 1.0
pair_coeff 1 2 0.0 1.0
neighbor 2.0 bin
neigh_modify every 1 delay 0 check yes
timestep 2.0
thermo 1000
thermo_style custom step temp pe press density
group ice molecule 1:{n_ice}
group mobile subtract all ice
fix freeze ice setforce 0.0 0.0 0.0
minimize 1.0e-3 1.0e-4 50 200
fix shake all shake 0.0001 20 0 b 1 a 1
velocity mobile create 300.0 240260 mom yes rot yes dist gaussian
velocity ice set 0.0 0.0 0.0
fix melt mobile nvt temp 300.0 300.0 200.0
run {steps['melt']}
unfix melt
fix cool mobile nvt temp ${{temp}} ${{temp}} 200.0
run {steps['cool']}
unfix cool
unfix freeze
velocity ice create ${{temp}} 17 mom yes rot yes dist gaussian
reset_timestep 0
fix bath all npt temp ${{temp}} ${{temp}} 200.0 z 1.0 1.0 2000.0
run {steps['equil']}
write_restart equil.restart
group oxygen type 1
reset_timestep 0
dump 1 oxygen custom 50 average.lammpstrj id type x y z
dump_modify 1 sort id
run 2000
""")
    run_lammps(work, "in.equil")
    write(work / "in.quench", """
units real
atom_style full
boundary p p p
atom_modify map array
read_restart equil.restart
kspace_style pppm/tip4p 1.0e-4
pair_style lj/cut/tip4p/long 1 2 1 1 0.1546 10.0
pair_modify mix arithmetic
pair_coeff 1 1 0.1852 3.1589
pair_coeff 2 2 0.0 1.0
pair_coeff 1 2 0.0 1.0
bond_style harmonic
bond_coeff 1 450.0 0.9572
angle_style harmonic
angle_coeff 1 55.0 104.52
special_bonds lj/coul 0.0 0.0 0.0
reset_atoms image all
neighbor 2.0 bin
minimize 1.0e-4 1.0e-6 80 300
fix shake all shake 0.0001 20 0 b 1 a 1
minimize 1.0e-4 1.0e-6 80 300
group oxygen type 1
write_dump oxygen custom quench.lammpstrj id type x y z modify sort id
""")
    run_lammps(work, "in.quench")


def prepare_mw() -> None:
    work = HOLD / "mw-274"
    if (work / "quench.lammpstrj").is_file() and (work / "average.lammpstrj").is_file():
        return
    work.mkdir(parents=True, exist_ok=True)
    shutil.copy(ROOT / "input" / "mw215" / "mW.sw", work / "mW.sw")
    shutil.copy(ROOT / "input" / "traj" / "interface_f0.lammpstrj", work / "interface_f0.lammpstrj")
    write(work / "in.equil", """
units real
atom_style atomic
boundary p p p
atom_modify map array
region box block 0 1 0 1 0 1
create_box 1 box
read_dump interface_f0.lammpstrj 0 x y z box yes add yes
mass 1 18.015
pair_style sw
pair_coeff * * mW.sw mW
neighbor 2.0 bin
neigh_modify every 1 delay 0 check yes
timestep 10.0
thermo 500
velocity all create 274.0 274 mom yes rot yes dist gaussian
fix bath all nvt temp 274.0 274.0 100.0
run 2000
write_restart equil.restart
reset_timestep 0
dump 1 all custom 10 average.lammpstrj id type x y z
dump_modify 1 sort id
run 400
""")
    run_lammps(work, "in.equil", threads=4)
    write(work / "in.quench", """
units real
atom_style atomic
boundary p p p
atom_modify map array
read_restart equil.restart
pair_style sw
pair_coeff * * mW.sw mW
neighbor 2.0 bin
minimize 1.0e-4 1.0e-6 200 800
write_dump all custom quench.lammpstrj id type x y z modify sort id
""")
    run_lammps(work, "in.quench", threads=4)


def prepare_brine() -> None:
    work = HOLD / "brine-230"
    if (work / "quench.lammpstrj").is_file() and (work / "average.lammpstrj").is_file():
        return
    work.mkdir(parents=True, exist_ok=True)
    import grow
    waters, na, cl, box, n_ice = grow.build_growth(3, 3, 2, 1.0, 1000 + 230 + 100)
    brine.write_data(work / "data.brine", waters, na, cl, box)
    write(work / "in.run", f"""
units real
atom_style full
boundary p p p
atom_modify map array
read_data data.brine
pair_style lj/cut/tip4p/long 1 2 1 1 0.1546 10.0
pair_modify mix arithmetic
kspace_style pppm/tip4p 1.0e-4
pair_coeff 1 1 0.1852 3.1589
pair_coeff 2 2 0.0 1.0
pair_coeff 3 3 0.351853250 2.21737
pair_coeff 4 4 0.018384799 4.69906
pair_coeff 1 3 0.189624283 2.60838
pair_coeff 1 4 0.014814293 4.23867
pair_coeff 3 4 0.343856119 3.00512
pair_coeff 1 2 0.0 1.0
pair_coeff 2 3 0.0 1.0
pair_coeff 2 4 0.0 1.0
bond_style harmonic
bond_coeff 1 450.0 0.9572
angle_style harmonic
angle_coeff 1 55.0 104.52
special_bonds lj/coul 0.0 0.0 0.0
reset_atoms image all
neighbor 2.0 bin
timestep 2.0
thermo 500
group ice molecule 1:{n_ice}
fix shake all shake 0.0001 20 0 b 1 a 1
velocity all create 230.0 230 mom yes rot yes dist gaussian
fix bath all nvt temp 230.0 230.0 200.0
reset_timestep 0
group oxygen type 1
dump 1 oxygen custom 50 average.lammpstrj id type x y z
dump_modify 1 sort id
run 2000
""")
    run_lammps(work, "in.run", threads=4)
    write(work / "in.quench", """
units real
atom_style full
boundary p p p
atom_modify map array
read_data data.brine
pair_style lj/cut/tip4p/long 1 2 1 1 0.1546 10.0
pair_modify mix arithmetic
kspace_style pppm/tip4p 1.0e-4
pair_coeff 1 1 0.1852 3.1589
pair_coeff 2 2 0.0 1.0
pair_coeff 3 3 0.351853250 2.21737
pair_coeff 4 4 0.018384799 4.69906
pair_coeff 1 3 0.189624283 2.60838
pair_coeff 1 4 0.014814293 4.23867
pair_coeff 3 4 0.343856119 3.00512
pair_coeff 1 2 0.0 1.0
pair_coeff 2 3 0.0 1.0
pair_coeff 2 4 0.0 1.0
bond_style harmonic
bond_coeff 1 450.0 0.9572
angle_style harmonic
angle_coeff 1 55.0 104.52
special_bonds lj/coul 0.0 0.0 0.0
reset_atoms image all
neighbor 2.0 bin
minimize 1.0e-4 1.0e-6 80 300
fix shake all shake 0.0001 20 0 b 1 a 1
minimize 1.0e-4 1.0e-6 40 200
group oxygen type 1
write_dump oxygen custom quench.lammpstrj id type x y z modify sort id
""")
    run_lammps(work, "in.quench", threads=4)


def prepare_liquid() -> None:
    work = HOLD / "tip4p-liquid-240"
    if (work / "equil.lammpstrj").is_file():
        return
    work.mkdir(parents=True, exist_ok=True)
    rng = random.Random(240)
    origins, box = brine.liquid_origins(512, rng)
    waters = [brine.tip4p_at(origin, rng) for origin in origins]
    brine.write_data(work / "data.brine", waters, [], [], box)
    write(work / "in.liquid", """
units real
atom_style full
boundary p p p
atom_modify map array
read_data data.brine
pair_style lj/cut/tip4p/long 1 2 1 1 0.1546 10.0
pair_modify mix arithmetic
kspace_style pppm/tip4p 1.0e-4
pair_coeff 1 1 0.1852 3.1589
pair_coeff 2 2 0.0 1.0
pair_coeff 3 3 0.351853250 2.21737
pair_coeff 4 4 0.018384799 4.69906
pair_coeff 1 3 0.189624283 2.60838
pair_coeff 1 4 0.014814293 4.23867
pair_coeff 3 4 0.343856119 3.00512
pair_coeff 1 2 0.0 1.0
pair_coeff 2 3 0.0 1.0
pair_coeff 2 4 0.0 1.0
bond_style harmonic
bond_coeff 1 450.0 0.9572
angle_style harmonic
angle_coeff 1 55.0 104.52
special_bonds lj/coul 0.0 0.0 0.0
reset_atoms image all
neighbor 2.0 bin
timestep 2.0
thermo 1000
fix shake all shake 0.0001 20 0 b 1 a 1
velocity all create 240.0 240 mom yes rot yes dist gaussian
fix bath all nvt temp 240.0 240.0 200.0
run 5000
group oxygen type 1
write_dump oxygen custom equil.lammpstrj id type x y z modify sort id
minimize 1.0e-4 1.0e-6 40 150
write_dump oxygen custom quench.lammpstrj id type x y z modify sort id
""")
    run_lammps(work, "in.liquid", threads=4)


def open_lib():
    import ctypes
    lib = templates.open_core(str(ROOT / "bbdir-release/src/libyodaLib.so"))
    lib.seams_symmetry_classify_gates.restype = ctypes.c_int
    lib.seams_symmetry_classify_gates.argtypes = [
        ctypes.POINTER(ctypes.c_double), ctypes.c_int, ctypes.POINTER(ctypes.c_double),
        ctypes.c_double, ctypes.c_double, ctypes.c_double, ctypes.POINTER(ctypes.c_int),
    ]
    return lib


def classify_gates(lib, pos, box, gates):
    import ctypes
    rmsd, thermal_q, bond_ratio = gates
    xyz = np.ascontiguousarray(pos, dtype=np.float64)
    lat = np.array(templates_lattice(box), dtype=np.float64).reshape(9)
    n = len(xyz)
    labels = np.empty(n, dtype=np.int32)
    rc = lib.seams_symmetry_classify_gates(
        xyz.ctypes.data_as(ctypes.POINTER(ctypes.c_double)), n,
        lat.ctypes.data_as(ctypes.POINTER(ctypes.c_double)),
        rmsd, thermal_q, bond_ratio,
        labels.ctypes.data_as(ctypes.POINTER(ctypes.c_int)))
    if rc != 0:
        raise RuntimeError("classify failed")
    return labels


def templates_lattice(box):
    return [[box[0], 0.0, 0.0], [0.0, box[1], 0.0], [0.0, 0.0, box[2]]]


def external(pos, box):
    lat = templates_lattice(box)
    ptm = templates.ovito_labels(pos, lat, "ice")
    chill = templates.ovito_labels(pos, lat, "chill")
    return (ptm == 6) | (ptm == 7), (chill == 1) | (chill == 2)


def align(frames):
    ids0 = frames[0][0]
    out = []
    for ids, pos, box in frames:
        if np.array_equal(ids, ids0):
            out.append((ids, pos, box))
            continue
        order = np.argsort(ids)
        ids = ids[order]
        pos = pos[order]
        if not np.array_equal(ids, ids0):
            raise RuntimeError("oxygen ids differ")
        out.append((ids, pos, box))
    return out


def majority(our, ptm_ice, chill_ice):
    return (our.astype(int) + ptm_ice.astype(int) + chill_ice.astype(int)) >= 2


def rates(pred, truth):
    tp = int(np.count_nonzero(pred & truth))
    fp = int(np.count_nonzero(pred & ~truth))
    fn = int(np.count_nonzero(~pred & truth))
    tn = int(np.count_nonzero(~pred & ~truth))
    precision = tp / (tp + fp) if tp + fp else float("nan")
    recall = tp / (tp + fn) if tp + fn else float("nan")
    return tp, fp, fn, tn, precision, recall


def matched(frame, ids0):
    ids, pos, box = frame
    if not np.array_equal(ids, ids0):
        order = np.argsort(ids)
        ids, pos = ids[order], pos[order]
        if not np.array_equal(ids, ids0):
            raise RuntimeError("oxygen ids differ")
    return pos, box


def score_system(name, role, temperature, thermal, quench, average, lib):
    ids0 = thermal[0]
    tpos, tbox = matched(thermal, ids0)
    qpos, qbox = matched(quench, ids0)
    apos, abox = matched(average, ids0)
    rows = []
    for source, pos, box in (("quench", qpos, qbox), ("average", apos, abox)):
        new = classify_gates(lib, pos, box, NEW) >= 0
        ptm_ice, chill_ice = external(pos, box)
        truth = majority(new, ptm_ice, chill_ice)
        old = classify_gates(lib, tpos, tbox, OLD) >= 0
        new_t = classify_gates(lib, tpos, tbox, NEW) >= 0
        ptm_t, chill_t = external(tpos, tbox)
        for method, pred in (
            ("old", old), ("new", new_t), ("PTM", ptm_t), ("CHILL+", chill_t),
        ):
            tp, fp, fn, tn, prec, rec = rates(pred, truth)
            rows.append({
                "system": name,
                "role": role,
                "temperature_K": temperature,
                "truth": source + "-majority",
                "method": method,
                "molecules": len(tpos),
                "truth_ice": int(np.count_nonzero(truth)),
                "predicted_ice": int(np.count_nonzero(pred)),
                "tp": tp, "fp": fp, "fn": fn, "tn": tn,
                "precision": f"{prec:.6f}" if prec == prec else "",
                "recall": f"{rec:.6f}" if rec == rec else "",
            })
    return rows


def score_liquid(lib):
    frames = read_oxygen_frames(HOLD / "tip4p-liquid-240" / "equil.lammpstrj")
    ids, pos, box = frames[0]
    rows = []
    for method, pred in (
        ("old", classify_gates(lib, pos, box, OLD) >= 0),
        ("new", classify_gates(lib, pos, box, NEW) >= 0),
    ):
        ice = int(np.count_nonzero(pred))
        rows.append({
            "system": "tip4p-liquid-240",
            "role": "bulk-liquid",
            "temperature_K": 240,
            "truth": "no-ice-seed",
            "method": method,
            "molecules": len(pos),
            "truth_ice": 0,
            "predicted_ice": ice,
            "tp": 0, "fp": ice, "fn": 0, "tn": len(pos) - ice,
            "precision": "",
            "recall": "",
        })
    ptm_ice, chill_ice = external(pos, box)
    for method, pred in (("PTM", ptm_ice), ("CHILL+", chill_ice)):
        ice = int(np.count_nonzero(pred))
        rows.append({
            "system": "tip4p-liquid-240",
            "role": "bulk-liquid",
            "temperature_K": 240,
            "truth": "no-ice-seed",
            "method": method,
            "molecules": len(pos),
            "truth_ice": 0,
            "predicted_ice": ice,
            "tp": 0, "fp": ice, "fn": 0, "tn": len(pos) - ice,
            "precision": "",
            "recall": "",
        })
    q = read_oxygen_frames(HOLD / "tip4p-liquid-240" / "quench.lammpstrj")[0]
    qnew = classify_gates(lib, q[1], q[2], NEW) >= 0
    rows.append({
        "system": "tip4p-liquid-240",
        "role": "bulk-liquid-quench",
        "temperature_K": 240,
        "truth": "no-ice-seed",
        "method": "new",
        "molecules": len(q[1]),
        "truth_ice": 0,
        "predicted_ice": int(np.count_nonzero(qnew)),
        "tp": 0, "fp": int(np.count_nonzero(qnew)), "fn": 0,
        "tn": len(q[1]) - int(np.count_nonzero(qnew)),
        "precision": "",
        "recall": "",
    })
    return rows


def thermal_frame(work: Path):
    frames = align(read_oxygen_frames(work / "average.lammpstrj"))
    return frames[0]


def load_ref(work: Path, name: str):
    frames = align(read_oxygen_frames(work / name))
    return frames[0]


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--out", type=Path, default=ROOT / "benchmarks/paper-cells")
    parser.add_argument("--skip-lammps", action="store_true")
    args = parser.parse_args()
    os.environ.setdefault("QT_QPA_PLATFORM", "offscreen")
    os.environ.setdefault("OVITO_THREAD_COUNT", "8")
    os.environ.setdefault("OMP_NUM_THREADS", "8")
    if not args.skip_lammps:
        jobs = {
            "mw-274": prepare_mw,
            "tip4p-240": lambda: prepare_tip4p(240),
            "tip4p-260": lambda: prepare_tip4p(260),
            "brine-230": prepare_brine,
            "tip4p-liquid-240": prepare_liquid,
        }
        with ThreadPoolExecutor(max_workers=2) as pool:
            futures = {pool.submit(fn): name for name, fn in jobs.items()}
            for fut in as_completed(futures):
                name = futures[fut]
                fut.result()
                print(f"prepared {name}", flush=True)
    lib = open_lib()
    rows = []
    catalog = [
        ("mw-274", "holdout", 274, HOLD / "mw-274"),
        ("tip4p-240", "holdout", 240, HOLD / "tip4p-240"),
        ("tip4p-260", "holdout", 260, HOLD / "tip4p-260"),
        ("brine-230", "holdout", 230, HOLD / "brine-230"),
        ("tip4p-252", "training", 252, ROOT / "data/hop/truth"),
    ]
    for name, role, temp, work in catalog:
        thermal = thermal_frame(work)
        quench = load_ref(work, "quench.lammpstrj")
        average = average_positions(align(read_oxygen_frames(work / "average.lammpstrj")))
        # average_positions returns ids, mean, box
        average = (average[0], average[1], average[2])
        rows.extend(score_system(name, role, temp, thermal, quench, average, lib))
        print(f"scored {name}", flush=True)
    rows.extend(score_liquid(lib))
    args.out.mkdir(parents=True, exist_ok=True)
    with (args.out / "symmetry-holdout.csv").open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(rows[0]))
        writer.writeheader()
        writer.writerows(rows)
    print(f"wrote {args.out / 'symmetry-holdout.csv'}", flush=True)
    return 0


if __name__ == "__main__":
    sys.exit(main())
