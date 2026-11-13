#!/usr/bin/env python3
"""Time seams against mdapy, OVITO, freud and rings-statistics-matter.

One ice Ic frame. The CLI is a fresh process per sample. The library and
the other packages time the call after that frame is loaded. The first
call is a warm-up and is not recorded. Eight samples follow. The report
is the median and the interquartile range. taskset pins the process to
the first N cores and the thread count is that same N.
"""

import math
import os
import statistics
import subprocess
import sys
import tempfile
import time
from pathlib import Path


def ice_ic(ncell, path):
    bond = 2.75
    lattice = 4.0 * bond / math.sqrt(3.0)
    side = ncell * lattice
    n = ncell * ncell * ncell * 8
    fcc = ((0, 0, 0), (0.5, 0.5, 0), (0.5, 0, 0.5), (0, 0.5, 0.5))
    with path.open("w") as out:
        out.write("ITEM: TIMESTEP\n0\n")
        out.write(f"ITEM: NUMBER OF ATOMS\n{n}\n")
        out.write("ITEM: BOX BOUNDS pp pp pp\n")
        out.write(f"0 {side}\n0 {side}\n0 {side}\n")
        out.write("ITEM: ATOMS id type x y z\n")
        aid = 1
        for i in range(ncell):
            for j in range(ncell):
                for k in range(ncell):
                    for b in range(8):
                        fx, fy, fz = fcc[b % 4]
                        shift = 0.0 if b < 4 else 0.25
                        x = ((i + fx + shift) * lattice) % side
                        y = ((j + fy + shift) * lattice) % side
                        z = ((k + fz + shift) * lattice) % side
                        out.write(f"{aid} 1 {x:.8f} {y:.8f} {z:.8f}\n")
                        aid += 1
    return n


def percentile(xs, p):
    k = (len(xs) - 1) * p
    lo = int(k)
    hi = min(lo + 1, len(xs) - 1)
    return xs[lo] * (1 - (k - lo)) + xs[hi] * (k - lo)


def summary(xs):
    xs = sorted(xs)
    return statistics.median(xs), percentile(xs, 0.75) - percentile(xs, 0.25)


def gl_dir():
    proc = subprocess.run(["ldconfig", "-p"], capture_output=True, text=True, check=False)
    for line in proc.stdout.splitlines():
        if "libGL.so.1" not in line or "=>" not in line:
            continue
        found = Path(line.split("=>", 1)[1].strip())
        if found.is_file():
            return str(found.parent)
    return ""


def pin_env(cores, with_gl=False):
    env = os.environ.copy()
    env["OMP_NUM_THREADS"] = str(cores)
    env["OMP_PROC_BIND"] = "close"
    env["OMP_WAIT_POLICY"] = "active"
    env["OMP_DYNAMIC"] = "FALSE"
    env["OPENBLAS_NUM_THREADS"] = str(cores)
    env["MKL_NUM_THREADS"] = str(cores)
    env["NUMEXPR_NUM_THREADS"] = str(cores)
    env["OVITO_THREAD_COUNT"] = str(cores)
    env["QT_QPA_PLATFORM"] = "offscreen"
    if with_gl:
        extra = gl_dir()
        if extra:
            current = env.get("LD_LIBRARY_PATH", "")
            env["LD_LIBRARY_PATH"] = extra if not current else extra + ":" + current
    else:
        env.pop("LD_LIBRARY_PATH", None)
    return env


def pinned(cores, argv, env):
    cpus = ",".join(str(i) for i in range(cores))
    return ["taskset", "-c", cpus, *argv], env


def samples(seams, args, dump, cores, n=8):
    cmd, env = pinned(
        cores, [str(seams), *args, str(dump), "--type", "1"], pin_env(cores, with_gl=False)
    )
    subprocess.run(cmd, env=env, check=True, stdout=subprocess.DEVNULL)
    got = []
    for _ in range(n):
        t0 = time.perf_counter()
        subprocess.run(
            cmd, env=env, check=True, stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL
        )
        got.append((time.perf_counter() - t0) * 1000.0)
    return summary(got)


def load_dump(path):
    import numpy as np

    lines = Path(path).read_text().splitlines()
    n = int(lines[lines.index("ITEM: NUMBER OF ATOMS") + 1])
    box_at = lines.index("ITEM: BOX BOUNDS pp pp pp")
    side = float(lines[box_at + 1].split()[1])
    atom_at = lines.index("ITEM: ATOMS id type x y z")
    pts = np.empty((n, 3), dtype=float)
    for i, line in enumerate(lines[atom_at + 1 : atom_at + 1 + n]):
        parts = line.split()
        pts[i, 0] = float(parts[2])
        pts[i, 1] = float(parts[3])
        pts[i, 2] = float(parts[4])
    return pts, side


def open_lib(lib_path):
    import ctypes

    lib = ctypes.CDLL(str(lib_path))
    lib.seams_frame_read.restype = ctypes.c_void_p
    lib.seams_frame_read.argtypes = [ctypes.c_char_p, ctypes.c_int, ctypes.c_int]
    lib.seams_frame_free.argtypes = [ctypes.c_void_p]
    lib.seams_frame_nop.restype = ctypes.c_int
    lib.seams_frame_nop.argtypes = [ctypes.c_void_p]
    lib.seams_frame_chill_plus.restype = ctypes.c_int
    lib.seams_frame_chill_plus.argtypes = [
        ctypes.c_void_p,
        ctypes.c_double,
        ctypes.c_int,
        ctypes.c_char_p,
        ctypes.POINTER(ctypes.c_int),
        ctypes.POINTER(ctypes.c_int),
    ]
    lib.seams_frame_q4_q6.restype = ctypes.c_int
    lib.seams_frame_q4_q6.argtypes = [
        ctypes.c_void_p,
        ctypes.c_double,
        ctypes.c_int,
        ctypes.POINTER(ctypes.c_double),
        ctypes.POINTER(ctypes.c_double),
        ctypes.c_void_p,
        ctypes.c_void_p,
    ]
    lib.seams_frame_cages.restype = ctypes.c_int
    lib.seams_frame_cages.argtypes = [
        ctypes.c_void_p,
        ctypes.c_double,
        ctypes.c_int,
        ctypes.c_int,
        ctypes.c_char_p,
        ctypes.c_int,
        ctypes.c_void_p,
        ctypes.POINTER(ctypes.c_int),
        ctypes.POINTER(ctypes.c_int),
        ctypes.POINTER(ctypes.c_int),
        ctypes.POINTER(ctypes.c_int),
    ]
    lib.seams_frame_rings.restype = ctypes.c_int
    lib.seams_frame_rings.argtypes = [
        ctypes.c_void_p,
        ctypes.c_double,
        ctypes.c_int,
        ctypes.c_char_p,
        ctypes.c_int,
        ctypes.POINTER(ctypes.c_int),
        ctypes.c_int,
    ]
    return lib


def read_frame(lib, dump):
    frame = lib.seams_frame_read(str(dump).encode(), 1, 1)
    if not frame:
        raise RuntimeError("frame")
    return frame


def time_calls(fn, n=8):
    fn()
    got = []
    sink = None
    for _ in range(n):
        t0 = time.perf_counter()
        sink = fn()
        got.append((time.perf_counter() - t0) * 1000.0)
    return got, sink


def worker(kind, cores, dump, lib_path):
    import numpy as np

    cores = int(cores)
    os.environ["OMP_NUM_THREADS"] = str(cores)
    pts, side = load_dump(dump)
    n = len(pts)

    if kind == "lib-chill":
        import ctypes

        lib = open_lib(lib_path)
        frame = read_frame(lib, dump)
        labels = (ctypes.c_int * n)()
        counts = (ctypes.c_int * 9)()

        def once():
            if lib.seams_frame_chill_plus(frame, 3.5, 1, b"cutoff", labels, counts) != 0:
                raise RuntimeError("chill")
            return int(np.count_nonzero(np.ctypeslib.as_array(labels) == 0))

        got, sink = time_calls(once)
        lib.seams_frame_free(frame)
    elif kind == "lib-q":
        import ctypes

        lib = open_lib(lib_path)
        frame = read_frame(lib, dump)
        q4 = (ctypes.c_double * n)()
        q6 = (ctypes.c_double * n)()

        def once():
            if lib.seams_frame_q4_q6(frame, 3.5, 1, q4, q6, None, None) != 0:
                raise RuntimeError("q")
            return float(q4[0])

        got, sink = time_calls(once)
        lib.seams_frame_free(frame)
    elif kind == "lib-cages":
        import ctypes

        lib = open_lib(lib_path)
        frame = read_frame(lib, dump)
        hexes = ctypes.c_int()
        cubic = ctypes.c_int()
        water = ctypes.c_int()
        both = ctypes.c_int()

        def once():
            if (
                lib.seams_frame_cages(
                    frame, 3.5, 1, 4, b"seeded", 0, None, hexes, cubic, water, both
                )
                != 0
            ):
                raise RuntimeError("cages")
            return int(cubic.value)

        got, sink = time_calls(once)
        lib.seams_frame_free(frame)
    elif kind == "lib-rings":
        import ctypes

        lib = open_lib(lib_path)
        frame = read_frame(lib, dump)
        counts = (ctypes.c_int * 9)()

        def once():
            if lib.seams_frame_rings(frame, 3.5, 1, b"cutoff", 6, counts, 9) != 0:
                raise RuntimeError("rings")
            return int(counts[6])

        got, sink = time_calls(once)
        lib.seams_frame_free(frame)
    elif kind == "mdapy-chill":
        from mdapy import System, set_num_threads
        from mdapy.identify import ChillPlus

        set_num_threads(cores)
        system = System(positions=pts, box=np.diag([side, side, side]), elements="O")
        tool = ChillPlus(3.5)

        def once():
            result = tool.compute(system)
            return int(np.count_nonzero(np.asarray(result.structure) == 2))

        got, sink = time_calls(once)
    elif kind == "ovito-chill":
        from ovito.modifiers import ChillPlusModifier

        def once():
            return ovito_codes(pts, side, ChillPlusModifier(cutoff=3.5), 2)

        got, sink = time_calls(once)
    elif kind == "ovito-ptm":
        from ovito.modifiers import PolyhedralTemplateMatchingModifier

        def once():
            # PTM cubic diamond is type 6 (ptm_template_dcub and alt1).
            # Hexagonal diamond is type 7 (ptm_template_dhex and three
            # alternates). Neither flag is in PTM_CHECK_DEFAULT, and OVITO
            # leaves both off, so a bare modifier labels this ice Other.
            mod = PolyhedralTemplateMatchingModifier(rmsd_cutoff=0.1)
            for structure in mod.structures:
                structure.enabled = structure.id in (6, 7)
            return ovito_codes(pts, side, mod, 6)

        got, sink = time_calls(once)
    elif kind == "ovito-diamond":
        from ovito.modifiers import IdentifyDiamondModifier

        def once():
            return ovito_codes(pts, side, IdentifyDiamondModifier(), 1)

        got, sink = time_calls(once)
    elif kind in ("freud-q4", "freud-q6"):
        import freud

        freud.parallel.set_num_threads(cores)
        box = freud.box.Box.cube(side)
        order = 4 if kind.endswith("4") else 6
        op = freud.order.Steinhardt(l=order, average=False)

        def once():
            op.compute((box, pts), neighbors={"r_max": 3.5})
            return float(op.ql[0])

        got, sink = time_calls(once)
    elif kind == "rings-rsm":
        from ase import Atoms
        import rings_statistics_matter as rsm

        atoms = Atoms(f"O{n}", positions=pts, cell=[side, side, side], pbc=True)

        def once():
            counts, _rings = rsm.ring_statistics(
                atoms, cutoff=3.5, maxlvl=8, num_threads=cores
            )
            return int(counts[5])

        got, sink = time_calls(once)
    elif kind == "agree":
        print_agreement(dump, lib_path, pts, side)
        return
    else:
        raise SystemExit(f"unknown worker {kind}")
    print(
        "samples="
        + ",".join(f"{v:.6f}" for v in got)
        + f" sink={sink}",
        flush=True,
    )


def ovito_codes(pts, side, modifier, code):
    import numpy as np
    from ovito.data import DataCollection, Particles, SimulationCell
    from ovito.pipeline import Pipeline, StaticSource

    data = DataCollection()
    cell = SimulationCell(pbc=(True, True, True))
    matrix = np.zeros((3, 4))
    matrix[0, 0] = side
    matrix[1, 1] = side
    matrix[2, 2] = side
    cell[:, :] = matrix
    parts = Particles()
    parts.create_property("Position", data=pts)
    data.objects.append(cell)
    data.objects.append(parts)
    pipe = Pipeline(source=StaticSource(data=data))
    pipe.modifiers.append(modifier)
    out = pipe.compute()
    labels = np.asarray(out.particles["Structure Type"])
    return int(np.count_nonzero(labels == code))


def print_agreement(dump, lib_path, pts, side):
    import ctypes
    import numpy as np
    import freud
    from ase import Atoms
    from mdapy import System, set_num_threads
    from mdapy.identify import ChillPlus
    from ovito.modifiers import (
        ChillPlusModifier,
        IdentifyDiamondModifier,
        PolyhedralTemplateMatchingModifier,
    )
    import rings_statistics_matter as rsm

    n = len(pts)
    lib = open_lib(lib_path)
    frame = read_frame(lib, dump)
    labels = (ctypes.c_int * n)()
    counts9 = (ctypes.c_int * 9)()
    if lib.seams_frame_chill_plus(frame, 3.5, 1, b"cutoff", labels, counts9) != 0:
        raise RuntimeError("chill")
    q4 = (ctypes.c_double * n)()
    q6 = (ctypes.c_double * n)()
    if lib.seams_frame_q4_q6(frame, 3.5, 1, q4, q6, None, None) != 0:
        raise RuntimeError("q")
    rings = (ctypes.c_int * 9)()
    if lib.seams_frame_rings(frame, 3.5, 1, b"cutoff", 6, rings, 9) != 0:
        raise RuntimeError("rings")
    hexes = ctypes.c_int()
    cubic = ctypes.c_int()
    water = ctypes.c_int()
    both = ctypes.c_int()
    if (
        lib.seams_frame_cages(
            frame, 3.5, 1, 4, b"seeded", 0, None, hexes, cubic, water, both
        )
        != 0
    ):
        raise RuntimeError("cages")
    lib.seams_frame_free(frame)
    chill = np.ctypeslib.as_array(labels).copy()
    set_num_threads(1)
    system = System(positions=pts, box=np.diag([side, side, side]), elements="O")
    mdapy_labels = np.asarray(ChillPlus(3.5).compute(system).structure)
    box = freud.box.Box.cube(side)
    freud.parallel.set_num_threads(1)
    gaps = []
    for order, buf in ((4, q4), (6, q6)):
        op = freud.order.Steinhardt(l=order, average=False)
        op.compute((box, pts), neighbors={"r_max": 3.5})
        ours = np.ctypeslib.as_array(buf)
        gaps.append(float(np.max(np.abs(ours - np.asarray(op.ql, dtype=float)))))
    atoms = Atoms(f"O{n}", positions=pts, cell=[side, side, side], pbc=True)
    rsm_counts, _rings = rsm.ring_statistics(atoms, cutoff=3.5, maxlvl=8, num_threads=1)
    print(f"n={n}")
    print(f"lib-cubic={int(np.count_nonzero(chill == 0))}")
    print(f"lib-other={int(np.count_nonzero(chill != 0))}")
    print(f"mdapy-cubic={int(np.count_nonzero(mdapy_labels == 2))}")
    print(f"mdapy-other={int(np.count_nonzero(mdapy_labels != 2))}")
    print(f"ovito-chill={ovito_codes(pts, side, ChillPlusModifier(cutoff=3.5), 2)}")
    ptm = PolyhedralTemplateMatchingModifier(rmsd_cutoff=0.1)
    for structure in ptm.structures:
        structure.enabled = structure.id in (6, 7)
    print(f"ovito-ptm-cubic={ovito_codes(pts, side, ptm, 6)}")
    print(f"ovito-diamond={ovito_codes(pts, side, IdentifyDiamondModifier(), 1)}")
    print(f"q-gap={max(gaps):.3e}")
    print(f"lib-rings6={int(rings[6])}")
    print(f"rsm-rings6={int(rsm_counts[5])}")
    print(f"lib-cages-cubic={int(cubic.value)}")
    print(f"lib-cages-hex={int(hexes.value)}")
    print(f"lib-cages-water={int(water.value)}")
    print(f"lib-cages-both={int(both.value)}")


def run_worker(kind, cores, dump, lib):
    cmd, env = pinned(
        cores,
        [sys.executable, str(Path(__file__).resolve()), "worker", kind, str(cores), str(dump), str(lib)],
        pin_env(cores, with_gl=True),
    )
    proc = subprocess.run(cmd, env=env, check=True, capture_output=True, text=True)
    line = proc.stdout.strip().splitlines()[-1]
    if not line.startswith("samples="):
        raise RuntimeError(line)
    body, sink = line.split(" sink=")
    got = [float(v) for v in body.split("=", 1)[1].split(",") if v]
    return summary(got), sink


def run_agree(cores, dump, lib):
    cmd, env = pinned(
        cores,
        [sys.executable, str(Path(__file__).resolve()), "worker", "agree", str(cores), str(dump), str(lib)],
        pin_env(cores, with_gl=True),
    )
    proc = subprocess.run(cmd, env=env, check=True, capture_output=True, text=True)
    out = {}
    for line in proc.stdout.splitlines():
        if "=" in line:
            key, value = line.split("=", 1)
            out[key] = value
    return out


def require_agreement(seams, dump, agreed):
    n = agreed["n"]
    env = pin_env(1, with_gl=False)
    chill = subprocess.run(
        ["taskset", "-c", "0", str(seams), "chill-plus", str(dump), "--type", "1"],
        env=env,
        check=True,
        capture_output=True,
        text=True,
    )
    sym = subprocess.run(
        ["taskset", "-c", "0", str(seams), "symmetry", str(dump), "--type", "1"],
        env=env,
        check=True,
        capture_output=True,
        text=True,
    )
    chill_line = chill.stdout.strip().splitlines()[-1]
    if f"cubic {n}" not in chill_line or "hexagonal" in chill_line or "water" in chill_line:
        raise SystemExit(f"CLI CHILL+ is not all cubic: {chill_line}")
    shells = [line for line in sym.stdout.splitlines() if line.startswith("shells ")]
    if not shells or f"ice-Ic {n}" not in shells[-1] or "unmatched 0" not in shells[-1]:
        raise SystemExit(f"CLI symmetry is not all ice Ic: {shells}")
    checks = {
        "lib-cubic": n,
        "lib-other": "0",
        "mdapy-cubic": n,
        "mdapy-other": "0",
        "ovito-chill": n,
        "ovito-ptm-cubic": n,
        "ovito-diamond": n,
        "lib-cages-cubic": n,
        "lib-cages-hex": "0",
        "lib-cages-water": "0",
        "lib-cages-both": "0",
    }
    for key, expect in checks.items():
        if agreed.get(key) != expect:
            raise SystemExit(f"{key}={agreed.get(key)} expected {expect}")
    if agreed["lib-rings6"] != agreed["rsm-rings6"]:
        raise SystemExit(
            f"rings differ {agreed['lib-rings6']} {agreed['rsm-rings6']}"
        )
    if float(agreed["q-gap"]) > 1e-5:
        raise SystemExit(f"q4/q6 gap {agreed['q-gap']}")
    return chill_line, shells[-1]


def main():
    root = Path(__file__).resolve().parents[1]
    seams = root / "build" / "seams"
    lib = root / "build" / "bindings" / "libseams_lang.so"
    if not seams.is_file() or not lib.is_file():
        raise SystemExit("build/seams or the language library is missing")
    import importlib.metadata as meta

    versions = {
        "mdapy": meta.version("mdapy"),
        "ovito": meta.version("ovito"),
        "freud": meta.version("freud-analysis"),
        "rings": meta.version("rings-statistics-matter"),
    }
    cores_available = os.cpu_count() or 1
    with tempfile.TemporaryDirectory(prefix="seams-cli-") as tmp:
        big = Path(tmp) / "iceIc.lammpstrj"
        n_big = ice_ic(24, big)
        agreed = run_agree(1, big, lib)
        chill_line, shell_line = require_agreement(seams, big, agreed)
        print(
            f"agree n={agreed['n']} q-gap={agreed['q-gap']} "
            f"rings6={agreed['lib-rings6']} {chill_line} {shell_line}",
            flush=True,
        )
        rows = []
        plan = [
            ("seams read", ["read"], "cli"),
            ("seams chill-plus", ["chill-plus"], "cli"),
            ("seams symmetry", ["symmetry"], "cli"),
            ("seams steinhardt --l 4,6", ["steinhardt", "--l", "4,6"], "cli"),
            ("seams cages", ["cages"], "cli"),
            ("library chill-plus", None, "lib-chill"),
            ("library q4 q6", None, "lib-q"),
            ("library cages", None, "lib-cages"),
            ("library rings", None, "lib-rings"),
            ("mdapy ChillPlus", None, "mdapy-chill"),
            ("OVITO ChillPlus", None, "ovito-chill"),
            ("OVITO PTM", None, "ovito-ptm"),
            ("OVITO identify diamond", None, "ovito-diamond"),
            ("freud q4", None, "freud-q4"),
            ("freud q6", None, "freud-q6"),
            ("rings-statistics-matter", None, "rings-rsm"),
        ]
        for cores in (1, 8):
            if cores > cores_available:
                continue
            for name, args, kind in plan:
                if kind == "cli":
                    med, iqr = samples(seams, args, big, cores)
                    sink = ""
                else:
                    (med, iqr), sink = run_worker(kind, cores, big, lib)
                rows.append((name, n_big, cores, med, iqr))
                print(
                    f"{name} n={n_big} cores={cores} median_ms={med:.2f} "
                    f"iqr_ms={iqr:.2f} sink={sink}",
                    flush=True,
                )
    out = root / "markdown" / "cli-fast-paths.md"
    out.write_text(render(rows, cores_available, versions, agreed))
    plot(rows, root / "markdown" / "cli-fast-plots")
    print(f"wrote {out.relative_to(root)}")


def render(rows, cores_available, versions, agreed):
    lines = [
        "# CLI times",
        "",
        "Ice Ic, orthorhombic, one oxygen per molecule, one frame of "
        f"{agreed['n']} molecules. The machine reported {cores_available} cores. "
        "Each row is one warm-up and eight samples, median and interquartile range, "
        "with the process pinned to that many cores and the thread count set to the same number. "
        "CLI rows include process start. Library and package rows time the call after the frame is loaded.",
        "",
        f"mdapy {versions['mdapy']}, OVITO {versions['ovito']}, "
        f"freud {versions['freud']}, rings-statistics-matter {versions['rings']}. "
        "Measured in this session on this frame.",
        "",
        "CHILL+ labels every molecule cubic in `seams chill-plus`, "
        "`seams_frame_chill_plus`, mdapy and OVITO. "
        "`seams symmetry` labels every molecule ice Ic, and OVITO's identify-diamond "
        "modifier labels every molecule cubic diamond. "
        "Polyhedral template matching with cubic and hexagonal diamond "
        "(PTM types 6 and 7, the 16-neighbour templates and their alternates) "
        "labels every molecule cubic diamond. "
        f"Local q4 and q6 from `seams_frame_q4_q6` and freud differ by at most {agreed['q-gap']}. "
        f"Primitive six-rings: the library and rings-statistics-matter both report {agreed['lib-rings6']}. "
        "`seams cages` and `seams_frame_cages` label every molecule cubic.",
        "",
        "| command | molecules | cores | median (ms) | IQR (ms) |",
        "| --- | ---: | ---: | ---: | ---: |",
    ]
    for name, n, cores, med, iqr in rows:
        lines.append(f"| `{name}` | {n} | {cores} | {med:.2f} | {iqr:.2f} |")
    lines.append("")
    return "\n".join(lines)


def plot(rows, directory):
    try:
        import matplotlib

        matplotlib.use("Agg")
        import matplotlib.pyplot as plt
    except Exception:
        return
    directory.mkdir(parents=True, exist_ok=True)
    names = []
    for name, _n, cores, _med, _iqr in rows:
        if cores == 1 and name not in names:
            names.append(name)
    for cores, marker in ((1, "o"), (8, "s")):
        xs = []
        ys = []
        for name in names:
            for row_name, _n, row_cores, med, _iqr in rows:
                if row_name == name and row_cores == cores:
                    xs.append(name)
                    ys.append(med)
        if ys:
            plt.plot(xs, ys, marker=marker, label=f"{cores} cores")
    plt.ylabel("median (ms)")
    plt.xticks(rotation=20, ha="right")
    plt.legend()
    plt.tight_layout()
    plt.savefig(directory / "cli-fast.png", dpi=120)
    plt.close()


if __name__ == "__main__":
    if len(sys.argv) > 1 and sys.argv[1] == "worker":
        worker(*sys.argv[2:])
    else:
        main()
