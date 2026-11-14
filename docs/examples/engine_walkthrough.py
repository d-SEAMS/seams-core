#!/usr/bin/env python3
"""One walkthrough for each engine hook, on that engine's own input.

LAMMPS: fix external callback in engines/lammps/in.dseams (stride 2, run 8).
GROMACS: analyzeFrame nanometre-to-angstrom read of engines/gromacs/ice.gro.
OpenMM: DseamsReporter for two steps of the reporter's 8-atom cell.
PLUMED: SEAMS action components on engines/plumed/ice.xyz.

Usage: engine_walkthrough.py ROOT BUILD
"""

from __future__ import annotations

import ctypes
import importlib.util
import os
import shutil
import subprocess
import sys
from pathlib import Path


def load(libdir: Path):
    lib = ctypes.CDLL(str(libdir / "libseams.so"))
    lib.seams_system_create.restype = ctypes.c_void_p
    lib.seams_system_destroy.argtypes = [ctypes.c_void_p]
    lib.seams_system_set_frame.argtypes = [
        ctypes.c_void_p,
        ctypes.POINTER(ctypes.c_double),
        ctypes.c_int,
        ctypes.POINTER(ctypes.c_double),
        ctypes.POINTER(ctypes.c_int),
        ctypes.c_void_p,
    ]
    lib.seams_system_set_frame.restype = ctypes.c_int
    lib.seams_chill_plus.argtypes = [
        ctypes.c_void_p,
        ctypes.c_double,
        ctypes.c_int,
        ctypes.c_int,
        ctypes.c_int,
        ctypes.POINTER(ctypes.c_int),
        ctypes.POINTER(ctypes.c_int),
    ]
    lib.seams_chill_plus.restype = ctypes.c_int
    lib.seams_rings_find.argtypes = [
        ctypes.c_void_p,
        ctypes.c_double,
        ctypes.c_int,
        ctypes.c_int,
        ctypes.c_int,
        ctypes.POINTER(ctypes.c_void_p),
    ]
    lib.seams_rings_find.restype = ctypes.c_int
    lib.seams_rings_histogram.argtypes = [
        ctypes.c_void_p,
        ctypes.POINTER(ctypes.c_int),
        ctypes.c_int,
    ]
    lib.seams_rings_histogram.restype = ctypes.c_int
    lib.seams_rings_fraction.argtypes = [
        ctypes.c_void_p,
        ctypes.c_int,
        ctypes.POINTER(ctypes.c_double),
    ]
    lib.seams_rings_fraction.restype = ctypes.c_int
    lib.seams_rings_destroy.argtypes = [ctypes.c_void_p]
    lib.seams_q4q6.argtypes = [
        ctypes.c_void_p,
        ctypes.c_double,
        ctypes.c_int,
        ctypes.POINTER(ctypes.c_double),
        ctypes.POINTER(ctypes.c_double),
    ]
    lib.seams_q4q6.restype = ctypes.c_int
    lib.seams_set_num_threads.argtypes = [ctypes.c_int]
    lib.seams_set_num_threads.restype = ctypes.c_int
    lib.seams_last_error.restype = ctypes.c_char_p
    if lib.seams_set_num_threads(1) != 0:
        raise SystemExit("threads")
    return lib


def census(lib, xyz, box):
    count = len(xyz) // 3
    points = (ctypes.c_double * len(xyz))(*xyz)
    cell = (ctypes.c_double * 9)(*box)
    types = (ctypes.c_int * count)(*([1] * count))
    system = lib.seams_system_create()
    if lib.seams_system_set_frame(system, points, count, cell, types, None) != 0:
        raise SystemExit(lib.seams_last_error())
    counts = (ctypes.c_int * 9)()
    labels = (ctypes.c_int * count)()
    if lib.seams_chill_plus(system, 3.5, 0, 0, 0, labels, counts) != 0:
        raise SystemExit(lib.seams_last_error())
    rings = ctypes.c_void_p()
    if lib.seams_rings_find(system, 3.5, 0, 0, 6, ctypes.byref(rings)) != 0:
        raise SystemExit(lib.seams_last_error())
    hist = (ctypes.c_int * 13)()
    fraction = ctypes.c_double()
    if lib.seams_rings_histogram(rings, hist, 13) != 0:
        raise SystemExit(lib.seams_last_error())
    if lib.seams_rings_fraction(rings, 6, ctypes.byref(fraction)) != 0:
        raise SystemExit(lib.seams_last_error())
    lib.seams_rings_destroy(rings)
    q4 = (ctypes.c_double * count)()
    q6 = (ctypes.c_double * count)()
    if lib.seams_q4q6(system, 3.5, 0, q4, q6) != 0:
        raise SystemExit(lib.seams_last_error())
    lib.seams_system_destroy(system)
    finite = [q6[i] for i in range(count) if q6[i] == q6[i]]
    return {
        "n": count,
        "cubic": int(counts[0]),
        "hexagonal": int(counts[1]),
        "six": int(hist[6]),
        "ring": float(fraction.value),
        "q6": sum(finite) / len(finite),
    }


def read_lammps(path: Path):
    xyz = []
    for line in path.read_text().splitlines():
        parts = line.split()
        if len(parts) == 5 and parts[0].isdigit() and parts[1].isdigit():
            xyz.extend(float(item) for item in parts[2:5])
    box = [0.0, 0.0, 0.0, 12.70170592, 12.70170592, 12.70170592, 0.0, 0.0, 0.0]
    return xyz, box


def read_gro(path: Path):
    lines = path.read_text().splitlines()
    count = int(lines[1])
    xyz = []
    raw = []
    for line in lines[2 : 2 + count]:
        coords = [float(line[20 + 8 * i : 28 + 8 * i]) for i in range(3)]
        raw.append(coords)
        xyz.extend(value * 10.0 for value in coords)
    lengths = [float(item) * 10.0 for item in lines[2 + count].split()]
    box = [0.0, 0.0, 0.0, lengths[0], lengths[1], lengths[2], 0.0, 0.0, 0.0]
    return xyz, box, raw


def read_xyz_nm(path: Path):
    lines = path.read_text().splitlines()
    count = int(lines[0])
    lattice = [float(item) for item in lines[1].split()]
    xyz = []
    for line in lines[2 : 2 + count]:
        _species, x, y, z = line.split()
        xyz.extend((float(x) * 10.0, float(y) * 10.0, float(z) * 10.0))
    box = [0.0, 0.0, 0.0, lattice[0] * 10.0, lattice[4] * 10.0, lattice[8] * 10.0, 0.0, 0.0, 0.0]
    return xyz, box


def hook_deck(path: Path):
    stride = None
    steps = None
    for line in path.read_text().splitlines():
        stripped = line.strip()
        if stripped.startswith("fix ") and "pf/callback" in stripped:
            parts = stripped.split()
            stride = int(parts[parts.index("pf/callback") + 1])
        if stripped.startswith("run "):
            steps = int(stripped.split()[1])
    if stride is None or steps is None:
        raise SystemExit(f"deck {path} has no fix callback")
    # run_hook sees step 0 and then every stride through the run.
    calls = steps // stride + 1
    return stride, steps, calls


def built(build: Path, name: str):
    for path in (build / name, build / "engines" / name, Path(name)):
        if path.is_file():
            return path
    found = shutil.which(name)
    return Path(found) if found else None


def figure(path: Path, lines):
    path.parent.mkdir(parents=True, exist_ok=True)
    text = "".join(
        f'<text x="16" y="{28 + 22 * i}" font-size="14">{line}</text>' for i, line in enumerate(lines)
    )
    path.write_text(
        f'<svg xmlns="http://www.w3.org/2000/svg" width="520" height="{40 + 22 * len(lines)}">'
        '<rect width="100%" height="100%" fill="white"/>'
        + text
        + "</svg>\n"
    )


def require(text, token):
    if token not in text:
        raise SystemExit(text)


def env_libs(build: Path):
    env = os.environ.copy()
    parts = [str(build / "src"), str(build / "bindings"), str(build / "engines")]
    current = env.get("LD_LIBRARY_PATH")
    if current:
        parts.append(current)
    env["LD_LIBRARY_PATH"] = os.pathsep.join(parts)
    return env


def main():
    root = Path(sys.argv[1])
    build = Path(sys.argv[2])
    libdir = Path(os.environ.get("SEAMS_LIBDIR", str(build / "src")))
    lib = load(libdir)
    figures = root / "docs" / "source" / "walkthrough" / "figures"

    deck = root / "engines" / "lammps" / "in.dseams"
    stride, steps, calls = hook_deck(deck)
    xyz, box = read_lammps(root / "engines" / "lammps" / "ice.data")
    got = census(lib, xyz, box)
    ice = got["cubic"] * calls
    if not (stride == 2 and steps == 8 and calls == 5 and got["cubic"] == 64 and ice == 320):
        raise SystemExit(f"lammps stride {stride} calls {calls} ice {ice} cubic {got['cubic']}")
    line = f"lammps calls {calls} stride {stride} ice {ice}"
    hook = built(build, "run_hook")
    if hook is not None:
        try:
            text = subprocess.check_output([str(hook), str(deck)], cwd=deck.parent, env=env_libs(build), text=True)
        except subprocess.CalledProcessError as exc:
            print(f"lammps hook binary was not runnable ({exc.returncode})")
        else:
            require(text, line)
            print("lammps hook binary")
    print(line)
    figure(figures / "lammps.svg", [line, "fix ext pf/callback 2 2"])

    xyz, box, raw = read_gro(root / "engines" / "gromacs" / "ice.gro")
    scaled = raw[1][0] * 10.0
    got = census(lib, xyz, box)
    if abs(raw[1][0] - 0.1588) > 1e-4 or abs(scaled - 1.588) > 1e-3 or got["n"] != 64 or got["cubic"] != 64:
        raise SystemExit(f"gromacs raw {raw[1]} cubic {got}")
    line = f"gromacs atoms {got['n']} scaled_A {scaled:.4f} ice {got['cubic']}"
    print(line)
    gmx = built(build, "gmx-seams")
    if gmx is not None:
        try:
            text = subprocess.check_output(
                [str(gmx), "-f", str(root / "engines" / "gromacs" / "ice.gro"), "-cutoff", "3.5"],
                cwd=root / "engines" / "gromacs",
                env=env_libs(build),
                text=True,
            )
        except subprocess.CalledProcessError as exc:
            print(f"gromacs binary was not runnable ({exc.returncode})")
        else:
            require(text, "gromacs frames 2 ice 128 atoms 128")
            print("gromacs frames 2 ice 128 atoms 128")
    figure(figures / "gromacs.svg", [line, "analyzeFrame nm to angstrom"])

    if importlib.util.find_spec("openmm") is None:
        print("openmm tutorial skipped: OpenMM was not found", file=sys.stderr)
        if os.environ.get("DOCS_CI") == "1":
            raise SystemExit(1)
    else:
        import numpy as np
        import openmm as mm
        import openmm.app as app
        import openmm.unit as unit

        sys.path.insert(0, str(root / "adapters" / "python"))
        from dseams_insitu.engine import InSitu
        from dseams_insitu.openmm import DseamsReporter

        os.environ["SEAMS_LANG_LIB"] = str(build / "bindings" / "libseams_lang.so")
        os.environ["LD_LIBRARY_PATH"] = env_libs(build)["LD_LIBRARY_PATH"]
        a = 2.75 * 4.0 / np.sqrt(3.0)
        fcc = ((0, 0, 0), (0.5, 0.5, 0), (0.5, 0, 0.5), (0, 0.5, 0.5))
        basis = ((0, 0, 0), (0.25, 0.25, 0.25))
        pos = np.asarray(
            [((shift[i] + site[i]) * a) for shift in fcc for site in basis for i in range(3)],
            dtype=np.float64,
        ).reshape(-1, 3)
        n = len(pos)
        system = mm.System()
        for _ in range(n):
            system.addParticle(16.0 * unit.dalton)
        box_len = a * unit.angstrom
        system.setDefaultPeriodicBoxVectors((box_len, 0, 0), (0, box_len, 0), (0, 0, box_len))
        force = mm.CustomExternalForce("0")
        for index in range(n):
            force.addParticle(index, [])
        system.addForce(force)
        topology = app.Topology()
        residue = topology.addResidue("ICE", topology.addChain())
        for _ in range(n):
            topology.addAtom("O", app.element.oxygen, residue)
        simulation = app.Simulation(
            topology,
            system,
            mm.VerletIntegrator(0.001 * unit.picoseconds),
            mm.Platform.getPlatformByName("Reference"),
        )
        simulation.context.setPositions(pos * unit.angstrom)
        engine = InSitu()
        simulation.reporters.append(DseamsReporter(engine, interval=1))
        simulation.step(2)
        if len(engine.frames) != 2 or n != 8:
            raise SystemExit(f"openmm frames {len(engine.frames)} atoms {n}")
        for frame in engine.frames:
            if not np.all(frame["labels"] == 0):
                raise SystemExit("openmm labels are not all cubic")
        line = f"openmm calls {len(engine.frames)} atoms {n} cubic {n}"
        print(line)
        figure(figures / "openmm.svg", [line, "DseamsReporter interval 1"])

    xyz, box = read_xyz_nm(root / "engines" / "plumed" / "ice.xyz")
    got = census(lib, xyz, box)
    if not (
        got["n"] == 64
        and got["cubic"] == 64
        and got["hexagonal"] == 0
        and got["six"] == 128
        and abs(got["ring"] - 1.0) < 1e-6
    ):
        raise SystemExit(got)
    line = (
        f"plumed ice {got['cubic'] / got['n']:.6f} ring {got['ring']:.6f} "
        f"hexagonal {got['hexagonal']} six {got['six']}"
    )
    print(line)
    plugin = built(build, "libSeamsPlumed.so")
    plumed = shutil.which("plumed")
    if plugin is not None and plumed is not None:
        work = build / "docs-engines" / "plumed"
        work.mkdir(parents=True, exist_ok=True)
        xyz_path = root / "engines" / "plumed" / "ice.xyz"
        deck = work / "plumed.dat"
        colvar = work / "colvar"
        deck.write_text(
            f"LOAD FILE={plugin}\n"
            "seams: SEAMS ATOMS=1-64 CUTOFF=3.5\n"
            "RESTRAINT ARG=seams.q6 AT=0.5 KAPPA=10.0\n"
            f"PRINT ARG=seams.q6,seams.ice,seams.ring,seams.hexagonal,seams.six FILE={colvar} STRIDE=1\n"
        )
        done = subprocess.run(
            [plumed, "driver", "--plumed", str(deck), "--ixyz", str(xyz_path)],
            cwd=work,
            env=env_libs(build),
            capture_output=True,
            text=True,
        )
        if done.returncode != 0:
            print(f"plumed driver was not runnable ({done.returncode})")
        else:
            require(colvar.read_text(), "1.000000 1.000000")
            print("plumed driver ice 1 ring 1")
    figure(figures / "plumed.svg", [line, "SEAMS ice and ring, q6 restrained"])
    print("engine walkthroughs ok")


if __name__ == "__main__":
    main()
