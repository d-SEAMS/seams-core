#!/usr/bin/env python3
# Copyright (c) d-SEAMS developers

"""Build an ice Ih basal slab in contact with liquid TIP4P/Ice.

The oxygen lattice is the GenIce ice Ih cell shipped in
input/traj/ice1h.lammpstrj. The basal normal is the axis whose layer
gaps match the ice Ih bilayer spacing, and that axis is rotated to z.
A liquid block of TIP4P/Ice is placed against the slab. Ice oxygens are
harmonically restrained and the liquid is equilibrated, then production
at SCIENCE_TEMP_K (default 230 K) is stored as an oxygen-only LAMMPS dump.
SCIENCE_PROD_PS sets the production length. Unset, the length is the
wall-clock budget SCIENCE_BUDGET_S, capped at 2 ns. SCIENCE_DUMP_PS sets
the spacing. SCIENCE_THREADS is the OpenMM CPU thread count.

TIP4P/Ice parameters are those of Abascal, Sanz, Fernandez and Vega,
J. Chem. Phys. 122, 234511 (2005). CHILL+ and q6 read the oxygens.

The dump is not committed. A file larger than 1 MB is listed in
paper/tex/figures/data/science/large-frames.csv.
"""

from __future__ import annotations

import math
import os
import random
from pathlib import Path

import numpy as np

# Abascal, Sanz, Fernandez and Vega, J. Chem. Phys. 122 (2005) 234511.
OH_A = 0.9572
HOH_DEG = 104.52
OM_A = 0.1577
Q_H = 0.5897
SIGMA_NM = 0.31668
EPSILON_KJ = 0.8822
CUTOFF_A = 3.5


def repo_root() -> Path:
    return Path(__file__).resolve().parents[1]


def read_ice(path: Path):
    lines = path.read_text().splitlines()
    box = []
    coords = []
    i = 0
    while i < len(lines):
        if lines[i].startswith("ITEM: BOX BOUNDS"):
            for k in range(3):
                lo, hi = lines[i + 1 + k].split()[:2]
                box.append(float(hi) - float(lo))
            i += 4
            continue
        if lines[i].startswith("ITEM: ATOMS"):
            i += 1
            while i < len(lines) and not lines[i].startswith("ITEM:"):
                parts = lines[i].split()
                coords.append((float(parts[2]), float(parts[3]), float(parts[4])))
                i += 1
            break
        i += 1
    if len(coords) == 0 or len(box) != 3:
        raise SystemExit(f"no ice coordinates in {path}")
    return np.asarray(box, dtype=float), np.asarray(coords, dtype=float)


def basal_axis(box: np.ndarray) -> int:
    """Box edge that is an integer multiple of the ice Ih c axis.

    The basal face is perpendicular to c. For the GenIce orthogonal cell
    c is 7.356 angstrom. A puckered bilayer does not put a 3.7 angstrom
    gap between atoms, so the lattice length is the test.
    """
    c_axis = 7.356
    best_axis = 0
    best_err = 1.0e9
    for axis in range(3):
        repeats = max(1, int(round(float(box[axis]) / c_axis)))
        err = abs(float(box[axis]) / repeats - c_axis)
        if err < best_err:
            best_err = err
            best_axis = axis
    if best_err > 0.05:
        raise SystemExit(f"no ice Ih c axis (nearest error {best_err:.3f} A)")
    return best_axis


def rotate_basal_to_z(box: np.ndarray, coords: np.ndarray):
    axis = basal_axis(box)
    if axis == 2:
        return box.copy(), coords.copy(), axis
    order = [i for i in range(3) if i != axis] + [axis]
    return box[order].copy(), coords[:, order].copy(), axis


def tile_axis(box: np.ndarray, coords: np.ndarray, axis: int, repeats: int):
    if repeats <= 1:
        return box, coords
    pieces = []
    for rep in range(repeats):
        shift = coords.copy()
        shift[:, axis] += rep * box[axis]
        pieces.append(shift)
    out = np.vstack(pieces)
    new_box = box.copy()
    new_box[axis] *= repeats
    return new_box, out


def liquid_grid(box_xy, z0: float, z1: float, spacing: float, rng: random.Random):
    count_x = max(1, int(box_xy[0] / spacing))
    count_y = max(1, int(box_xy[1] / spacing))
    count_z = max(1, int((z1 - z0) / spacing))
    dx = box_xy[0] / count_x
    dy = box_xy[1] / count_y
    dz = (z1 - z0) / count_z
    points = []
    for i in range(count_x):
        for j in range(count_y):
            for k in range(count_z):
                points.append(
                    (
                        (i + 0.5) * dx + rng.uniform(-0.15, 0.15),
                        (j + 0.5) * dy + rng.uniform(-0.15, 0.15),
                        z0 + (k + 0.5) * dz + rng.uniform(-0.15, 0.15),
                    )
                )
    return np.asarray(points, dtype=float)


def random_rotation(rng: random.Random) -> np.ndarray:
    # Shoemake's uniform quaternion.
    u1 = rng.random()
    u2 = rng.uniform(0, 2 * math.pi)
    u3 = rng.uniform(0, 2 * math.pi)
    s = math.sqrt(1.0 - u1)
    c = math.sqrt(u1)
    w = c * math.sin(u2)
    x = c * math.cos(u2)
    y = s * math.sin(u3)
    z = s * math.cos(u3)
    return np.array(
        [
            [1 - 2 * (y * y + z * z), 2 * (x * y - z * w), 2 * (x * z + y * w)],
            [2 * (x * y + z * w), 1 - 2 * (x * x + z * z), 2 * (y * z - x * w)],
            [2 * (x * z - y * w), 2 * (y * z + x * w), 1 - 2 * (x * x + y * y)],
        ]
    )


def tip4p_hydrogens(oxygen: np.ndarray, rng: random.Random) -> np.ndarray:
    half = math.radians(HOH_DEG / 2.0)
    template = np.array(
        [
            [math.sin(half) * OH_A, 0.0, math.cos(half) * OH_A],
            [-math.sin(half) * OH_A, 0.0, math.cos(half) * OH_A],
        ]
    )
    rotated = template @ random_rotation(rng).T
    return oxygen + rotated


def virtual_weights() -> tuple[float, float]:
    beta = OM_A / (2.0 * OH_A * math.cos(math.radians(HOH_DEG / 2.0)))
    return 1.0 - 2.0 * beta, beta


def forcefield_xml() -> str:
    weight_o, weight_h = virtual_weights()
    # OpenMM XML lengths are nanometres and energies are kJ/mol.
    return f"""<?xml version="1.0"?>
<ForceField>
 <AtomTypes>
  <Type name="tip4pice-O" class="OW" element="O" mass="15.99943"/>
  <Type name="tip4pice-H" class="HW" element="H" mass="1.007947"/>
  <Type name="tip4pice-M" class="MW" mass="0.0"/>
 </AtomTypes>
 <Residues>
  <Residue name="HOH">
   <Atom name="O" type="tip4pice-O"/>
   <Atom name="H1" type="tip4pice-H"/>
   <Atom name="H2" type="tip4pice-H"/>
   <Atom name="M" type="tip4pice-M"/>
   <VirtualSite type="average3" siteName="M" atomName1="O" atomName2="H1" atomName3="H2" weight1="{weight_o:.8f}" weight2="{weight_h:.8f}" weight3="{weight_h:.8f}"/>
   <Bond atomName1="O" atomName2="H1"/>
   <Bond atomName1="O" atomName2="H2"/>
  </Residue>
 </Residues>
 <HarmonicBondForce>
  <Bond class1="OW" class2="HW" length="0.09572" k="462750.4"/>
 </HarmonicBondForce>
 <HarmonicAngleForce>
  <Angle class1="HW" class2="OW" class3="HW" angle="1.82421813418" k="836.8"/>
 </HarmonicAngleForce>
 <NonbondedForce coulomb14scale="0.833333" lj14scale="0.5">
  <Atom type="tip4pice-O" charge="0.0" sigma="{SIGMA_NM}" epsilon="{EPSILON_KJ}"/>
  <Atom type="tip4pice-H" charge="{Q_H}" sigma="1.0" epsilon="0.0"/>
  <Atom type="tip4pice-M" charge="{-2 * Q_H}" sigma="1.0" epsilon="0.0"/>
 </NonbondedForce>
</ForceField>
"""


def record_large(path: Path, root: Path) -> None:
    csv_path = root / "paper" / "tex" / "figures" / "data" / "science" / "large-frames.csv"
    csv_path.parent.mkdir(parents=True, exist_ok=True)
    rel = path.relative_to(root).as_posix()
    kept = []
    if csv_path.is_file():
        for line in csv_path.read_text().splitlines():
            if not line or line.startswith("path ") or line.startswith(rel + " "):
                continue
            kept.append(line)
    lines = ["path bytes", *kept]
    size = path.stat().st_size
    if size > 1_000_000:
        lines.append(f"{rel} {size}")
    csv_path.write_text("\n".join(lines) + "\n")


def write_frame(handle, step: int, oxygens: np.ndarray, box: np.ndarray) -> None:
    n = oxygens.shape[0]
    handle.write("ITEM: TIMESTEP\n")
    handle.write(f"{step}\n")
    handle.write("ITEM: NUMBER OF ATOMS\n")
    handle.write(f"{n}\n")
    handle.write("ITEM: BOX BOUNDS pp pp pp\n")
    for length in box:
        handle.write(f"0.0 {length:.8f}\n")
    handle.write("ITEM: ATOMS id type x y z\n")
    for i, (x, y, z) in enumerate(oxygens, start=1):
        handle.write(f"{i} 1 {x:.8f} {y:.8f} {z:.8f}\n")


def write_dump(path: Path, frames: list[tuple[int, np.ndarray]], box: np.ndarray) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w") as handle:
        for step, oxygens in frames:
            write_frame(handle, step, oxygens, box)


def wrap_coordinates(positions_nm: np.ndarray, box_nm: np.ndarray) -> np.ndarray:
    out = positions_nm.copy()
    for axis in range(3):
        length = box_nm[axis]
        out[:, axis] = np.mod(out[:, axis], length)
    return out


def run_openmm(ice: np.ndarray, liquid: np.ndarray, box_a: np.ndarray):
    import openmm
    import openmm.app
    import openmm.unit as unit

    rng = random.Random(27011)
    n_ice = ice.shape[0]
    oxygens = np.vstack([ice, liquid])
    n_mol = oxygens.shape[0]
    top = openmm.app.Topology()
    chain = top.addChain()
    positions = np.zeros((n_mol * 4, 3), dtype=float)
    for i, oxygen in enumerate(oxygens):
        residue = top.addResidue("HOH", chain)
        oxygen_atom = top.addAtom("O", openmm.app.element.oxygen, residue)
        hydrogen_a = top.addAtom("H1", openmm.app.element.hydrogen, residue)
        hydrogen_b = top.addAtom("H2", openmm.app.element.hydrogen, residue)
        top.addAtom("M", None, residue)
        top.addBond(oxygen_atom, hydrogen_a)
        top.addBond(oxygen_atom, hydrogen_b)
        hydrogens = tip4p_hydrogens(oxygen, rng)
        base = i * 4
        positions[base] = oxygen
        positions[base + 1] = hydrogens[0]
        positions[base + 2] = hydrogens[1]
        positions[base + 3] = oxygen
    box_nm = box_a / 10.0
    vec_a = openmm.Vec3(float(box_nm[0]), 0.0, 0.0) * unit.nanometer
    vec_b = openmm.Vec3(0.0, float(box_nm[1]), 0.0) * unit.nanometer
    vec_c = openmm.Vec3(0.0, 0.0, float(box_nm[2])) * unit.nanometer
    top.setPeriodicBoxVectors((vec_a, vec_b, vec_c))
    xml = forcefield_xml()
    xml_path = repo_root() / "data" / "science-frames" / "tip4pice.xml"
    xml_path.parent.mkdir(parents=True, exist_ok=True)
    xml_path.write_text(xml)
    ff = openmm.app.ForceField(str(xml_path))
    shortest = float(min(box_nm))
    cutoff_nm = min(0.9, 0.49 * shortest)
    system = ff.createSystem(
        top,
        nonbondedMethod=openmm.app.PME,
        nonbondedCutoff=cutoff_nm * unit.nanometer,
        constraints=openmm.app.HBonds,
        rigidWater=True,
    )
    restraint = openmm.CustomExternalForce("k*((x-x0)^2+(y-y0)^2+(z-z0)^2)")
    restraint.addGlobalParameter("k", 2000.0)
    restraint.addPerParticleParameter("x0")
    restraint.addPerParticleParameter("y0")
    restraint.addPerParticleParameter("z0")
    for i in range(n_ice):
        oxygen = ice[i] / 10.0
        restraint.addParticle(i * 4, oxygen.tolist())
    system.addForce(restraint)
    equil_ps = float(os.environ.get("SCIENCE_EQUIL_PS", "15"))
    budget_s = float(os.environ.get("SCIENCE_BUDGET_S", "900"))
    temperature = float(os.environ.get("SCIENCE_TEMP_K", "230"))
    dt_ps = 0.002
    integrator = openmm.LangevinMiddleIntegrator(
        320 * unit.kelvin, 1.0 / unit.picosecond, dt_ps * unit.picoseconds
    )
    threads = str(int(os.environ.get("SCIENCE_THREADS", "4")))
    platform = openmm.Platform.getPlatformByName("CPU")
    simulation = openmm.app.Simulation(
        top, system, integrator, platform, {"Threads": threads}
    )
    simulation.context.setPositions(positions / 10.0 * unit.nanometer)
    simulation.context.setPeriodicBoxVectors(vec_a, vec_b, vec_c)
    simulation.context.applyConstraints(1e-5)
    simulation.context.computeVirtualSites()
    simulation.minimizeEnergy(maxIterations=400)
    equil_steps = int(round(equil_ps / dt_ps))
    simulation.step(max(equil_steps, 1))
    state = simulation.context.getState(getEnergy=True)
    energy = state.getPotentialEnergy().value_in_unit(unit.kilojoule_per_mole)
    if not math.isfinite(energy):
        raise SystemExit("equilibration energy is not finite")
    integrator.setTemperature(temperature * unit.kelvin)
    cool_steps = int(round(5.0 / dt_ps))
    simulation.step(cool_steps)
    probe = 2000
    import time

    t0 = time.perf_counter()
    simulation.step(probe)
    elapsed = max(time.perf_counter() - t0, 1e-6)
    steps_per_s = probe / elapsed
    prod_ps_env = os.environ.get("SCIENCE_PROD_PS")
    if prod_ps_env:
        prod_steps = int(round(float(prod_ps_env) / dt_ps))
    else:
        prod_steps = int(budget_s * steps_per_s)
        prod_steps = max(prod_steps, int(20.0 / dt_ps))
        prod_steps = min(prod_steps, int(2000.0 / dt_ps))
    dump_ps_env = os.environ.get("SCIENCE_DUMP_PS")
    if dump_ps_env:
        stride = max(1, int(round(float(dump_ps_env) / dt_ps)))
    else:
        stride = max(1, prod_steps // 16)
    n_dump = max(1, prod_steps // stride)
    prod_steps = stride * n_dump
    dump_ps = stride * dt_ps
    print(
        f"production {prod_steps * dt_ps:.1f} ps at {steps_per_s:.0f} steps/s "
        f"on {threads} threads, dump every {dump_ps:.3f} ps ({n_dump} frames)",
        flush=True,
    )
    if os.environ.get("SCIENCE_PROBE_ONLY") == "1":
        meta = {
            "model": "TIP4P/Ice",
            "temperature_K": f"{temperature:.1f}",
            "prod_ps": f"{prod_steps * dt_ps:.3f}",
            "dump_ps": f"{dump_ps:.6f}",
            "openmm_threads": threads,
            "steps_per_s": f"{steps_per_s:.1f}",
            "n_frames": "0",
        }
        return [], meta
    frames = []
    stream = os.environ.get("SCIENCE_STREAM", "").strip()
    stream_path = Path(stream) if stream else None
    if stream_path is not None and not stream_path.is_absolute():
        stream_path = repo_root() / stream_path
    clock = time.perf_counter()

    def grab(step_index: int) -> None:
        state = simulation.context.getState(getPositions=True, enforcePeriodicBox=True)
        pos = state.getPositions(asNumpy=True).value_in_unit(unit.angstrom)
        oxy = np.asarray(pos[0::4], dtype=float)
        frames.append((step_index, wrap_coordinates(oxy, box_a)))
        if stream_path is not None:
            mode = "w" if len(frames) == 1 else "a"
            with stream_path.open(mode) as handle:
                write_frame(handle, step_index, frames[-1][1], box_a)
        done = len(frames) - 1
        elapsed_prod = time.perf_counter() - clock
        rate = (done * stride / elapsed_prod) if done and elapsed_prod > 0 else steps_per_s
        left = ((n_dump - done) * stride / rate) if rate else 0
        print(
            f"dump {done}/{n_dump} step {step_index} waters {oxy.shape[0]} "
            f"eta_s {left:.0f}",
            flush=True,
        )

    step = equil_steps + cool_steps + probe
    grab(step)
    for _ in range(n_dump):
        simulation.step(stride)
        step += stride
        grab(step)
    meta = {
        "model": "TIP4P/Ice",
        "citation": "Abascal 2005 JCP 122 234511",
        "temperature_K": f"{temperature:.1f}",
        "equil_K": "320",
        "equil_ps": f"{equil_ps:.3f}",
        "prod_ps": f"{prod_steps * dt_ps:.3f}",
        "dump_ps": f"{dump_ps:.6f}",
        "restraint_kJ_mol_nm2": "2000",
        "n_ice": str(n_ice),
        "n_liquid": str(liquid.shape[0]),
        "n_frames": str(len(frames)),
        "cutoff_A": f"{CUTOFF_A:.1f}",
        "equil_energy_kJ_mol": f"{energy:.6f}",
        "openmm_threads": threads,
        "steps_per_s": f"{steps_per_s:.1f}",
    }
    return frames, meta


def main() -> int:
    root = repo_root()
    box, coords = read_ice(root / "input" / "traj" / "ice1h.lammpstrj")
    box, coords, source_axis = rotate_basal_to_z(box, coords)
    coords = coords - coords.min(axis=0)
    # PME needs a cutoff below half the shortest lateral edge.
    for axis in (0, 1):
        repeats = 1
        while box[axis] * repeats < 24.0:
            repeats += 1
        box, coords = tile_axis(box, coords, axis, repeats)
    cell_z = float(box[2])
    box, ice = tile_axis(box, coords, 2, 2)
    ice_z0 = 2.5
    ice[:, 2] += ice_z0
    # The crystal occupies the periodic cells, not only the outermost atoms.
    ice_z1 = ice_z0 + float(box[2])
    gap = 2.6
    liquid_length = cell_z
    liquid = liquid_grid(
        box[:2], ice_z1 + gap, ice_z1 + gap + liquid_length, 3.05, random.Random(27011)
    )
    box = box.copy()
    box[2] = ice_z1 + gap + liquid_length + gap
    print(
        f"ice {ice.shape[0]} liquid {liquid.shape[0]} box {box[0]:.2f} {box[1]:.2f} {box[2]:.2f}",
        flush=True,
    )
    tag = os.environ.get("SCIENCE_TAG", "").strip()
    name = f"tip4p-basal-{tag}.lammpstrj" if tag else "tip4p-basal.lammpstrj"
    dest = root / "data" / "science-frames" / name
    dest.parent.mkdir(parents=True, exist_ok=True)
    if os.environ.get("SCIENCE_PROBE_ONLY") != "1":
        os.environ["SCIENCE_STREAM"] = str(dest)
    frames, meta = run_openmm(ice, liquid, box)
    if frames and not os.environ.get("SCIENCE_STREAM"):
        write_dump(dest, frames, box)
    if frames:
        record_large(dest, root)
    meta_name = os.environ.get("SCIENCE_META", "").strip()
    if meta_name:
        meta_path = Path(meta_name)
        if not meta_path.is_absolute():
            meta_path = root / meta_path
    else:
        meta_path = root / "paper" / "tex" / "figures" / "data" / "science" / "interface-meta.csv"
    meta_path.parent.mkdir(parents=True, exist_ok=True)
    rows = [
        "key value",
        "basal_source_axis " + "xyz"[source_axis],
        "basal_axis z",
        f"box_A {box[0]:.6f} {box[1]:.6f} {box[2]:.6f}",
        "ice_oxygens_restrained yes",
        "analyzer oxygen",
    ]
    for key, value in meta.items():
        rows.append(f"{key} {value}")
    meta_path.write_text("\n".join(rows) + "\n")
    print(
        f"ok frames={len(frames)} dump_ps={meta['dump_ps']} "
        f"path={dest.relative_to(root).as_posix()}"
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
