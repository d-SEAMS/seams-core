#!/usr/bin/env python3
# Copyright (c) d-SEAMS developers
"""Run one engine how-to on input/traj/ice1h.lammpstrj. Exit 77 if it is absent."""

import ctypes
import importlib.util
import os
import shutil
import subprocess
import sys
import tempfile
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[2] / "scripts"))
from shared_library import built_library, require_library  # noqa: E402


def _load_e2e():
    """engines/e2e_common.py, without putting engines/ on sys.path.

    engines/openmm, engines/lammps and the other engine folders would then
    import as namespace packages, and a missing OpenMM would pass for an
    installed one.
    """
    import importlib.util

    path = Path(__file__).resolve().parents[2] / "engines" / "e2e_common.py"
    spec = importlib.util.spec_from_file_location("e2e_common", path)
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


e2e = _load_e2e()

HEXAGONAL = 128
SIX_RINGS = 256
Q6 = 0.625855


def built(build, name):
    for path in (Path(build) / name, Path(build) / "engines" / name):
        if path.is_file():
            return path
    return None


def skip(message):
    print(message)
    raise SystemExit(77)


def env_libs(build):
    env = os.environ.copy()
    parts = [str(Path(build) / "src"), str(Path(build) / "bindings")]
    prefix = env.get("CONDA_PREFIX")
    if prefix:
        root = Path(prefix)
        parts.extend(str(root / name) for name in ("lib", "lib.AVX2_256", "lib.AVX_256"))
    current = env.get("LD_LIBRARY_PATH")
    if current:
        parts.append(current)
    env["LD_LIBRARY_PATH"] = os.pathsep.join(parts)
    return env


def read_dump(path):
    lines = Path(path).read_text().splitlines()
    count = None
    bounds = None
    atoms = []
    index = 0
    while index < len(lines):
        line = lines[index]
        if line.startswith("ITEM: NUMBER OF ATOMS"):
            count = int(lines[index + 1])
            index += 2
            continue
        if line.startswith("ITEM: BOX BOUNDS"):
            bounds = []
            for offset in range(1, 4):
                lo, hi = lines[index + offset].split()[:2]
                bounds.append((float(lo), float(hi)))
            index += 4
            continue
        if line.startswith("ITEM: ATOMS"):
            header = line.split()[2:]
            xcol = header.index("x")
            ycol = header.index("y")
            zcol = header.index("z")
            for row in lines[index + 1 : index + 1 + count]:
                cols = row.split()
                atoms.append((float(cols[xcol]), float(cols[ycol]), float(cols[zcol])))
            break
        index += 1
    if count is None or bounds is None or len(atoms) != count:
        raise SystemExit("unreadable dump {}".format(path))
    return atoms, bounds


def require_census(text):
    token = "hexagonal {} six-rings {} q6 ".format(HEXAGONAL, SIX_RINGS)
    if token not in text:
        raise SystemExit(text)
    start = text.index(token) + len(token)
    raw = text[start:].split()[0]
    if abs(float(raw) - Q6) > 5e-7:
        raise SystemExit(text)
    print(text.strip())


def run_gromacs(root, build, work):
    exe = built(build, "gmx-seams")
    if exe is None:
        skip("gromacs skipped: gmx-seams was not built")
    atoms, bounds = read_dump(root / "input" / "traj" / "ice1h.lammpstrj")
    gro = work / "ice1h.gro"
    lines = ["ice Ih", str(len(atoms))]
    for index, (x, y, z) in enumerate(atoms, 1):
        lines.append(
            "{:5d}{:<5s}{:>5s}{:5d}{:8.4f}{:8.4f}{:8.4f}".format(
                1, "ICE", "O", index, x / 10.0, y / 10.0, z / 10.0
            )
        )
    lengths = [hi - lo for lo, hi in bounds]
    lines.append("{:10.5f}{:10.5f}{:10.5f}".format(*(length / 10.0 for length in lengths)))
    gro.write_text("\n".join(lines) + "\n")
    text = subprocess.check_output(
        [str(exe), "-f", str(gro), "-cutoff", "3.5"],
        env=e2e.with_thread_runtime(env_libs(build)),
        text=True,
    )
    require_census(text)


def run_lammps(root, build, work):
    exe = built(build, "run_hook")
    if exe is None:
        skip("lammps skipped: run_hook was not built")
    atoms, bounds = read_dump(root / "input" / "traj" / "ice1h.lammpstrj")
    data = work / "ice1h.data"
    rows = [
        "LAMMPS data file from ice Ih",
        "",
        "{} atoms".format(len(atoms)),
        "1 atom types",
        "",
    ]
    for axis, (lo, hi) in zip(("x", "y", "z"), bounds):
        rows.append("{} {} {}lo {}hi".format(lo, hi, axis, axis))
    rows.extend(["", "Atoms", ""])
    for index, (x, y, z) in enumerate(atoms, 1):
        rows.append("{} 1 {:.10f} {:.10f} {:.10f}".format(index, x, y, z))
    data.write_text("\n".join(rows) + "\n")
    deck = work / "in.ice1h"
    deck.write_text(
        "\n".join(
            [
                "units metal",
                "atom_style atomic",
                "boundary p p p",
                "read_data ice1h.data",
                "mass 1 16.0",
                "pair_style lj/cut 8.0",
                "pair_coeff * * 0.1 2.75",
                "velocity all set 0.0 0.0 0.0",
                "fix freeze all setforce 0.0 0.0 0.0",
                "fix ext all external pf/callback 1 1",
                "run 0",
                "",
            ]
        )
    )
    text = subprocess.check_output(
        [str(exe), str(deck)],
        cwd=work,
        env=e2e.drop_thread_runtime(env_libs(build)),
        text=True,
    )
    require_census(text)


def mean_q6(build, positions, box9):
    os.environ["LD_LIBRARY_PATH"] = env_libs(build)["LD_LIBRARY_PATH"]
    lib = ctypes.CDLL(str(require_library(Path(build) / "bindings", "seams_lang", "SEAMS_LANG_LIB")))
    lib.seams_frame_from_arrays.argtypes = [
        ctypes.POINTER(ctypes.c_double),
        ctypes.c_int,
        ctypes.POINTER(ctypes.c_double),
        ctypes.POINTER(ctypes.c_int),
        ctypes.POINTER(ctypes.c_int),
    ]
    lib.seams_frame_from_arrays.restype = ctypes.c_void_p
    lib.seams_frame_steinhardt_graph.argtypes = [
        ctypes.c_void_p,
        ctypes.c_double,
        ctypes.c_int,
        ctypes.c_int,
        ctypes.c_char_p,
        ctypes.POINTER(ctypes.c_double),
        ctypes.POINTER(ctypes.c_double),
    ]
    lib.seams_frame_steinhardt_graph.restype = ctypes.c_int
    lib.seams_frame_free.argtypes = [ctypes.c_void_p]
    lib.seams_last_error.restype = ctypes.c_char_p
    count = len(positions)
    flat = []
    for point in positions:
        flat.extend(point)
    xyz = (ctypes.c_double * len(flat))(*flat)
    box = (ctypes.c_double * 9)(*box9)
    types = (ctypes.c_int * count)(*([1] * count))
    ids = (ctypes.c_int * count)(*range(1, count + 1))
    frame = lib.seams_frame_from_arrays(xyz, count, box, types, ids)
    if not frame:
        raise SystemExit("openmm frame was rejected")
    ql = (ctypes.c_double * count)()
    code = lib.seams_frame_steinhardt_graph(frame, 3.5, 1, 6, b"cutoff", ql, None)
    lib.seams_frame_free(frame)
    if code != 0:
        raise SystemExit("openmm q6 failed")
    return sum(ql) / count


def run_openmm(root, build, work):
    if importlib.util.find_spec("openmm") is None:
        skip("openmm skipped: the OpenMM package was not found")
    import numpy as np
    import openmm as mm
    import openmm.app as app
    import openmm.unit as unit

    sys.path.insert(0, str(root / "adapters" / "python"))
    from dseams_insitu.engine import InSitu
    from dseams_insitu.openmm import DseamsReporter

    os.environ["SEAMS_LANG_LIB"] = str(require_library(Path(build) / "bindings", "seams_lang", "SEAMS_LANG_LIB"))
    os.environ["LD_LIBRARY_PATH"] = env_libs(build)["LD_LIBRARY_PATH"]
    atoms, bounds = read_dump(root / "input" / "traj" / "ice1h.lammpstrj")
    pos = np.asarray(atoms, dtype=np.float64)
    lengths = [hi - lo for lo, hi in bounds]
    system = mm.System()
    for _ in range(len(pos)):
        system.addParticle(16.0 * unit.dalton)
    system.setDefaultPeriodicBoxVectors(
        (lengths[0], 0, 0) * unit.angstrom,
        (0, lengths[1], 0) * unit.angstrom,
        (0, 0, lengths[2]) * unit.angstrom,
    )
    force = mm.CustomExternalForce("0")
    for index in range(len(pos)):
        force.addParticle(index, [])
    system.addForce(force)
    topology = app.Topology()
    chain = topology.addChain()
    residue = topology.addResidue("ICE", chain)
    element = app.element.oxygen
    for _ in range(len(pos)):
        topology.addAtom("O", element, residue)
    integrator = mm.VerletIntegrator(0.001 * unit.picoseconds)
    platform = mm.Platform.getPlatformByName("Reference")
    simulation = app.Simulation(topology, system, integrator, platform)
    simulation.context.setPositions(pos * unit.angstrom)
    engine = InSitu()
    simulation.reporters.append(DseamsReporter(engine, interval=1))
    simulation.step(1)
    if len(engine.frames) != 1:
        raise SystemExit("openmm frames {}".format(len(engine.frames)))
    frame = engine.frames[0]
    q6 = mean_q6(build, frame["positions"], frame["box"])
    text = "openmm hexagonal {} six-rings {} q6 {:.6f}\n".format(
        int(frame["counts"][1]), int(frame["rings"][6]), q6
    )
    require_census(text)


def run_plumed(root, build, work):
    if shutil.which("plumed") is None:
        skip("plumed skipped: plumed was not found")
    plugin = e2e.external_library(build, "libSeamsPlumed.so")
    if plugin is None:
        skip("plumed skipped: libSeamsPlumed.so was not built")
    atoms, bounds = read_dump(root / "input" / "traj" / "ice1h.lammpstrj")
    lengths = [(hi - lo) / 10.0 for lo, hi in bounds]
    xyz = work / "ice1h.xyz"
    rows = [
        str(len(atoms)),
        "{:.8f} 0 0 0 {:.8f} 0 0 0 {:.8f}".format(*lengths),
    ]
    for x, y, z in atoms:
        rows.append("O {:.8f} {:.8f} {:.8f}".format(x / 10.0, y / 10.0, z / 10.0))
    xyz.write_text("\n".join(rows) + "\n")
    deck = work / "plumed.dat"
    colvar = work / "colvar"
    deck.write_text(
        "LOAD FILE={}\n"
        "seams: SEAMS ATOMS=1-{} CUTOFF=3.5\n"
        "PRINT ARG=seams.q6,seams.hexagonal,seams.six FILE={} STRIDE=1 FMT=%12.6f\n".format(
            plugin, len(atoms), colvar
        )
    )
    env = e2e.drop_thread_runtime(env_libs(build))
    env.pop("LD_PRELOAD", None)
    env["PLUMED_NUM_THREADS"] = "1"
    done = subprocess.run(
        ["plumed", "driver", "--plumed", str(deck), "--ixyz", str(xyz)],
        check=False,
        capture_output=True,
        text=True,
        env=env,
    )
    if done.returncode != 0 or not colvar.is_file():
        sys.stderr.write(done.stdout)
        sys.stderr.write(done.stderr)
        raise SystemExit(done.returncode or 1)
    fields = None
    for line in colvar.read_text().splitlines():
        if not line or line.startswith("#") or line.startswith("#!"):
            continue
        fields = line.split()
    if fields is None or len(fields) < 4:
        raise SystemExit(colvar.read_text())
    text = "plumed hexagonal {:.0f} six-rings {:.0f} q6 {}\n".format(
        float(fields[2]), float(fields[3]), fields[1]
    )
    require_census(text)


def run_insitu(root, build, work):
    libpath = built_library(Path(build) / "bindings", "seams_lang", "SEAMS_LANG_LIB")
    if libpath is None:
        skip("insitu skipped: libseams_lang was not built")
    os.environ["LD_LIBRARY_PATH"] = env_libs(build)["LD_LIBRARY_PATH"]
    lib = ctypes.CDLL(str(libpath))
    lib.seams_frame_read.argtypes = [ctypes.c_char_p, ctypes.c_int, ctypes.c_int]
    lib.seams_frame_read.restype = ctypes.c_void_p
    lib.seams_frame_nop.argtypes = [ctypes.c_void_p]
    lib.seams_frame_nop.restype = ctypes.c_int
    lib.seams_frame_chill_plus.argtypes = [
        ctypes.c_void_p,
        ctypes.c_double,
        ctypes.c_int,
        ctypes.c_char_p,
        ctypes.POINTER(ctypes.c_int),
        ctypes.POINTER(ctypes.c_int),
    ]
    lib.seams_frame_chill_plus.restype = ctypes.c_int
    lib.seams_frame_rings.argtypes = [
        ctypes.c_void_p,
        ctypes.c_double,
        ctypes.c_int,
        ctypes.c_char_p,
        ctypes.c_int,
        ctypes.POINTER(ctypes.c_int),
        ctypes.c_int,
    ]
    lib.seams_frame_rings.restype = ctypes.c_int
    lib.seams_frame_free.argtypes = [ctypes.c_void_p]
    lib.seams_last_error.restype = ctypes.c_char_p
    dump = root / "input" / "traj" / "ice1h.lammpstrj"
    frame = lib.seams_frame_read(str(dump).encode(), 1, 0)
    if not frame:
        raise SystemExit("unreadable frame")
    count = lib.seams_frame_nop(frame)
    labels = (ctypes.c_int * count)()
    counts = (ctypes.c_int * 9)()
    rings = (ctypes.c_int * 7)()
    chill = lib.seams_frame_chill_plus(frame, 3.5, 1, b"cutoff", labels, counts)
    ring = lib.seams_frame_rings(frame, 3.5, 1, b"cutoff", 6, rings, 7)
    lib.seams_frame_free(frame)
    if chill != 0 or ring != 0 or count != 128 or counts[1] != 128 or rings[6] != 256:
        raise SystemExit("insitu n {} counts {} rings {}".format(count, list(counts), list(rings)))
    print("insitu hexagonal {} six-rings {}".format(counts[1], rings[6]))


def main():
    engine, root, build = sys.argv[1], Path(sys.argv[2]), Path(sys.argv[3])
    run = {
        "gromacs": run_gromacs,
        "lammps": run_lammps,
        "openmm": run_openmm,
        "plumed": run_plumed,
        "insitu": run_insitu,
    }[engine]
    # docs_engine_<engine> and the video tape run the same engine at once,
    # and PLUMED keeps a numbered backup of every colvar it finds, 101 at
    # most, so each run writes into a directory of its own.
    base = build / "docs-examples"
    base.mkdir(parents=True, exist_ok=True)
    with tempfile.TemporaryDirectory(prefix=engine + "-", dir=base) as work:
        run(root, build, Path(work))


if __name__ == "__main__":
    main()
