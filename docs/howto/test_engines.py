#!/usr/bin/env python3
"""Run one engine how-to when that engine is built. Exit 77 otherwise."""

import ctypes
import importlib.util
import os
import shutil
import subprocess
import sys
from pathlib import Path


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
    current = env.get("LD_LIBRARY_PATH")
    if current:
        parts.append(current)
    env["LD_LIBRARY_PATH"] = os.pathsep.join(parts)
    return env


def run_gromacs(root, build):
    exe = built(build, "gmx-seams")
    if exe is None:
        skip("gromacs skipped: gmx-seams was not built")
    gro = root / "engines" / "gromacs" / "ice.gro"
    text = subprocess.check_output(
        [str(exe), "-f", str(gro), "-cutoff", "3.5"],
        env=env_libs(build),
        text=True,
    )
    if "gromacs frames" not in text or " ice " not in text:
        raise SystemExit(text)
    print(text.strip())


def run_lammps(root, build):
    exe = built(build, "run_hook")
    if exe is None:
        skip("lammps skipped: run_hook was not built")
    deck = root / "engines" / "lammps" / "in.dseams"
    text = subprocess.check_output(
        [str(exe), str(deck)],
        cwd=root / "engines" / "lammps",
        env=env_libs(build),
        text=True,
    )
    parts = text.split()
    calls = int(parts[parts.index("calls") + 1])
    ice = int(parts[parts.index("ice") + 1])
    if calls < 1 or ice < calls:
        raise SystemExit(text)
    print(text.strip())


def run_openmm(root, build):
    if importlib.util.find_spec("openmm") is None:
        skip("openmm skipped: the OpenMM package was not found")
    script = root / "engines" / "openmm" / "test_reporter.py"
    env = env_libs(build)
    env["PYTHONPATH"] = os.pathsep.join(
        [str(root / "adapters" / "python"), env.get("PYTHONPATH", "")]
    )
    env["SEAMS_LANG_LIB"] = str(Path(build) / "bindings" / "libseams_lang.so")
    text = subprocess.check_output([sys.executable, str(script)], env=env, text=True)
    if "openmm calls" not in text or "ice" not in text:
        raise SystemExit(text)
    print(text.strip())


def run_plumed(root, build):
    if shutil.which("plumed") is None:
        skip("plumed skipped: plumed was not found")
    plugin = built(build, "libSeamsPlumed.so")
    if plugin is None:
        skip("plumed skipped: libSeamsPlumed.so was not built")
    work = Path(build) / "docs-examples"
    work.mkdir(parents=True, exist_ok=True)
    deck = work / "plumed.dat"
    colvar = work / "colvar"
    xyz = root / "engines" / "plumed" / "ice.xyz"
    deck.write_text(
        "LOAD FILE={}\n"
        "seams: SEAMS ATOMS=1-64 CUTOFF=3.5\n"
        "PRINT ARG=seams.q6,seams.ice,seams.ring FILE={} STRIDE=1\n".format(plugin, colvar)
    )
    done = subprocess.run(
        ["plumed", "driver", "--plumed", str(deck), "--ixyz", str(xyz)],
        check=False,
        capture_output=True,
        text=True,
        env=env_libs(build),
    )
    if done.returncode != 0:
        sys.stderr.write(done.stdout)
        sys.stderr.write(done.stderr)
        raise SystemExit(done.returncode)
    ice = None
    for line in colvar.read_text().splitlines():
        if not line or line.startswith("#"):
            continue
        fields = line.split()
        if len(fields) >= 3:
            ice = float(fields[2])
    if ice is None or abs(ice - 1.0) > 1e-6:
        raise SystemExit("seams.ice is {}".format(ice))
    print("plumed seams.ice", ice)


def run_insitu(root, build):
    libpath = Path(build) / "bindings" / "libseams_lang.so"
    if not libpath.is_file():
        skip("insitu skipped: libseams_lang.so was not built")
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
    n = lib.seams_frame_nop(frame)
    labels = (ctypes.c_int * n)()
    counts = (ctypes.c_int * 9)()
    rings = (ctypes.c_int * 7)()
    chill = lib.seams_frame_chill_plus(frame, 3.5, 1, b"cutoff", labels, counts)
    ring = lib.seams_frame_rings(frame, 3.5, 1, b"cutoff", 6, rings, 7)
    lib.seams_frame_free(frame)
    if chill != 0 or ring != 0 or n != 128 or counts[1] != 128 or rings[6] != 256:
        raise SystemExit("insitu n {} counts {} rings {}".format(n, list(counts), list(rings)))
    print("insitu hexagonal {} six-rings {}".format(counts[1], rings[6]))


def main():
    engine, root, build = sys.argv[1], Path(sys.argv[2]), Path(sys.argv[3])
    {
        "gromacs": run_gromacs,
        "lammps": run_lammps,
        "openmm": run_openmm,
        "plumed": run_plumed,
        "insitu": run_insitu,
    }[engine](root, build)


if __name__ == "__main__":
    main()
