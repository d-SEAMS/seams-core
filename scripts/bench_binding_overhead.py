#!/usr/bin/env python3
# Copyright (c) d-SEAMS developers
"""Per-language overhead of the C API on ice Ih.

One process per binding stays up for the whole run. The serial phase sets
OMP_NUM_THREADS and PLUMED_NUM_THREADS to 1. The threaded phase sets both
to the core count and calls seams_set_num_threads. Eight repeats follow one
warmup. The median of an even sample is the mean of the two central values.
Quartiles are Hyndman-Fan type 7. The log keeps every sample. The CSV keeps
the median and the interquartile range.

GenIce writes an ice Ih cell when that size is not already on disk. Oxygen
sites are kept and converted from nanometres to angstroms, the same reduction
as input/traj/ice1h.lammpstrj. Raw logs and those cells stay out of git.
"""

from __future__ import annotations

import argparse
import csv
import os
import shutil
import subprocess
import sys
import threading
import time
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
SIZES = (4096, 32768, 110592)
OPS = ("read", "knn", "steinhardt", "chill", "rings")
POINTER_BINDINGS = ("c", "cxx", "rust", "fortran", "python", "julia", "r", "lua")
GUEST_OPS = {
    "freud": ("knn", "steinhardt"),
    "ovito": ("read", "knn", "steinhardt", "rings"),
}
CUTOFF = 3.5


def machine():
    cores = os.cpu_count() or 1
    model = "unknown"
    cpuinfo = Path("/proc/cpuinfo")
    if cpuinfo.is_file():
        for line in cpuinfo.read_text(errors="replace").splitlines():
            if line.lower().startswith("model name"):
                model = line.split(":", 1)[1].strip()
                break
    return cores, model


def quantile_linear(ordered, q):
    if len(ordered) == 1:
        return ordered[0]
    pos = q * (len(ordered) - 1)
    lo = int(pos)
    frac = pos - lo
    if lo + 1 >= len(ordered):
        return ordered[-1]
    return ordered[lo] * (1.0 - frac) + ordered[lo + 1] * frac


def summarize(samples):
    ordered = sorted(float(v) for v in samples)
    count = len(ordered)
    mid = count // 2
    if count % 2:
        median = ordered[mid]
    else:
        median = 0.5 * (ordered[mid - 1] + ordered[mid])
    q1 = quantile_linear(ordered, 0.25)
    q3 = quantile_linear(ordered, 0.75)
    return median, q3 - q1, q1, q3


def find_library():
    env = os.environ.get("SEAMS_LIB")
    if env and Path(env).is_file():
        return Path(env)
    for path in (ROOT / "build" / "src" / "libseams.so", ROOT / "bbdir" / "src" / "libseams.so"):
        if path.is_file():
            return path
    raise SystemExit("libseams.so is not built")


def find_lua_module():
    found = []
    for base in (ROOT / "build", ROOT / "bbdir"):
        if base.is_dir():
            found.extend(base.rglob("dseams_core.so"))
    jit = [path for path in found if "luajit" in path.as_posix()]
    if jit:
        exe = shutil.which("luajit")
        if exe:
            return Path(exe), jit[0]
    lua54 = [path for path in found if "lua54" in path.as_posix()]
    exe = shutil.which("lua5.4") or shutil.which("lua")
    if lua54 and exe:
        return Path(exe), lua54[0]
    raise SystemExit("the Lua module dseams_core.so is not built")


def _run(cmd, env=None, cwd=None):
    proc = subprocess.run(cmd, cwd=cwd or ROOT, env=env, text=True, capture_output=True)
    if proc.returncode != 0:
        sys.stderr.write(proc.stdout[-2000:])
        sys.stderr.write(proc.stderr[-4000:])
        raise SystemExit("command failed: " + " ".join(cmd[:6]))


def compile_drivers(lib):
    include = ROOT / "include"
    version = lib.parent
    out = lib.parent.parent
    out.mkdir(parents=True, exist_ok=True)
    rpath = "-Wl,-rpath," + str(lib.parent)
    cc = shutil.which("cc") or "cc"
    cxx = shutil.which("c++") or shutil.which("g++") or "c++"
    _run([cc, "-O2", "-I" + str(include), "-I" + str(version),
          str(ROOT / "benchmarks" / "binding_overhead_c.c"),
          "-L" + str(lib.parent), "-lseams", rpath,
          "-o", str(out / "bench_binding_c")])
    _run([cxx, "-std=c++20", "-O2", "-I" + str(include), "-I" + str(version),
          str(ROOT / "benchmarks" / "binding_overhead_cpp.cpp"),
          "-L" + str(lib.parent), "-lseams", rpath,
          "-o", str(out / "bench_binding_cpp")])
    gfortran = shutil.which("gfortran")
    if not gfortran:
        raise SystemExit("gfortran is not on PATH")
    mod = out / "mod"
    mod.mkdir(exist_ok=True)
    _run([gfortran, "-O2", "-J", str(mod), "-c", str(ROOT / "bindings" / "fortran" / "seams.f90"),
          "-o", str(mod / "seams.o")])
    _run([gfortran, "-O2", "-J", str(mod), str(ROOT / "benchmarks" / "binding_overhead.f90"),
          str(mod / "seams.o"), "-L" + str(lib.parent), "-lseams", rpath,
          "-o", str(out / "bench_binding_fortran")])
    rscript = shutil.which("Rscript")
    if not rscript:
        raise SystemExit("Rscript is not on PATH")
    env = os.environ.copy()
    env["PKG_CPPFLAGS"] = "-I" + str(include) + " -I" + str(version)
    env["PKG_LIBS"] = "-L" + str(lib.parent) + " -lseams " + rpath
    so = out / "seams_call.so"
    _run([str(Path(rscript).with_name("R")), "CMD", "SHLIB", "-o", str(so),
          str(ROOT / "bindings" / "r" / "dseams" / "src" / "seams_call.c")], env=env)
    cargo = shutil.which("cargo")
    if not cargo:
        raise SystemExit("cargo is not on PATH")
    env = os.environ.copy()
    env["SEAMS_LIB_DIR"] = str(lib.parent)
    _run([cargo, "build", "--release", "--example", "binding_overhead",
          "--manifest-path", str(ROOT / "bindings" / "rust" / "seams" / "Cargo.toml")], env=env)
    return {
        "c": out / "bench_binding_c",
        "cxx": out / "bench_binding_cpp",
        "fortran": out / "bench_binding_fortran",
        "rso": so,
        "rust": ROOT / "bindings" / "rust" / "target" / "release" / "examples" / "binding_overhead",
    }


class Worker:
    def __init__(self, argv, env):
        self.proc = subprocess.Popen(
            argv, cwd=ROOT, env=env, stdin=subprocess.PIPE, stdout=subprocess.PIPE,
            stderr=subprocess.PIPE, text=True, bufsize=1,
        )
        self.err = []
        self._thread = threading.Thread(target=self._drain, daemon=True)
        self._thread.start()

    def _drain(self):
        for line in self.proc.stderr:
            self.err.append(line)

    def ask(self, threads, op, path):
        try:
            self.proc.stdin.write("{} {} {}\n".format(threads, op, path))
            self.proc.stdin.flush()
        except BrokenPipeError as exc:
            tail = "".join(self.err)[-2000:]
            raise RuntimeError("worker closed before the reply\n" + tail) from exc
        line = self.proc.stdout.readline()
        if not line:
            tail = "".join(self.err)[-2000:]
            raise RuntimeError("worker closed\n" + tail)
        if line.startswith("error"):
            raise RuntimeError(line.strip())
        return float(line)

    def close(self):
        if self.proc.poll() is None:
            try:
                self.proc.stdin.write("quit\n")
                self.proc.stdin.flush()
            except BrokenPipeError:
                pass
            self.proc.wait(timeout=30)


def child_env(lib, extra=None):
    env = os.environ.copy()
    env["OMP_NUM_THREADS"] = "1"
    env["PLUMED_NUM_THREADS"] = "1"
    env["OMP_DYNAMIC"] = "FALSE"
    env["PYTHONDONTWRITEBYTECODE"] = "1"
    libdir = str(lib.parent)
    current = [part for part in env.get("LD_LIBRARY_PATH", "").split(":") if part and part != libdir]
    env["LD_LIBRARY_PATH"] = ":".join([libdir, *current])
    env["SEAMS_LIB"] = str(lib)
    if extra:
        env.update(extra)
    return env


def _gro_oxygens(text):
    lines = [line for line in text.splitlines() if line.strip()]
    count = int(lines[1])
    body = lines[2:2 + count]
    box = [float(part) for part in lines[2 + count].split()[:3]]
    oxy = []
    for line in body:
        parts = line.split()
        if len(parts) < 5:
            continue
        name = parts[1]
        # Past atom 9999 the fixed gro columns glue the element to the index.
        if name not in ("O", "OW") and not (name.startswith("O") and name[1:].isdigit()):
            continue
        oxy.append(tuple(float(part) * 10.0 for part in parts[-3:]))
    lengths = tuple(edge * 10.0 for edge in box)
    return oxy, lengths


def _write_dump(path, oxy, lengths):
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w") as handle:
        handle.write("ITEM: TIMESTEP\n0\nITEM: NUMBER OF ATOMS\n{}\n".format(len(oxy)))
        handle.write("ITEM: BOX BOUNDS pp pp pp\n")
        for edge in lengths:
            handle.write("0.0 {:.10f}\n".format(edge))
        handle.write("ITEM: ATOMS id type x y z\n")
        for index, (x, y, z) in enumerate(oxy, start=1):
            handle.write("{} 1 {:.8f} {:.8f} {:.8f}\n".format(index, x, y, z))


def _factor(cells):
    best = None
    for nx in range(1, cells + 1):
        if cells % nx:
            continue
        rest = cells // nx
        for ny in range(1, rest + 1):
            if rest % ny:
                continue
            nz = rest // ny
            score = max(nx, ny, nz) / min(nx, ny, nz)
            if best is None or score < best[0]:
                best = (score, nx, ny, nz)
    return best[1:]


def ensure_frames(directory):
    directory.mkdir(parents=True, exist_ok=True)
    missing = [n for n in SIZES if not (directory / "ih-{}.lammpstrj".format(n)).is_file()]
    if not missing:
        return {n: directory / "ih-{}.lammpstrj".format(n) for n in SIZES}
    genice = shutil.which("genice2")
    if not genice:
        raise SystemExit("genice2 is not on PATH")
    unit = subprocess.run(
        [genice, "ice1h", "--rep", "1", "1", "1", "--format", "gromacs"],
        check=True, text=True, capture_output=True,
    )
    oxy, _lengths = _gro_oxygens(unit.stdout)
    per_cell = len(oxy)
    if per_cell <= 0:
        raise SystemExit("GenIce ice Ih cell has no oxygen")
    paths = {}
    for molecules in SIZES:
        path = directory / "ih-{}.lammpstrj".format(molecules)
        paths[molecules] = path
        if path.is_file():
            continue
        if molecules % per_cell:
            raise SystemExit("ice Ih cell of {} does not divide {}".format(per_cell, molecules))
        nx, ny, nz = _factor(molecules // per_cell)
        proc = subprocess.run(
            [genice, "ice1h", "--rep", str(nx), str(ny), str(nz), "--format", "gromacs"],
            check=True, text=True, capture_output=True,
        )
        sites, lengths = _gro_oxygens(proc.stdout)
        if len(sites) != molecules:
            raise SystemExit("GenIce wrote {} oxygens, want {}".format(len(sites), molecules))
        _write_dump(path, sites, lengths)
    return paths


def python_worker(lib):
    import ctypes

    import numpy as np

    library = ctypes.CDLL(str(lib))
    library.seams_set_num_threads.argtypes = [ctypes.c_int]
    library.seams_set_num_threads.restype = ctypes.c_int
    library.seams_system_create.restype = ctypes.c_void_p
    library.seams_system_read.argtypes = [ctypes.c_void_p, ctypes.c_char_p, ctypes.c_int, ctypes.c_int]
    library.seams_system_read.restype = ctypes.c_int
    library.seams_system_count.argtypes = [ctypes.c_void_p]
    library.seams_system_count.restype = ctypes.c_int
    library.seams_handle_error.argtypes = [ctypes.c_void_p]
    library.seams_handle_error.restype = ctypes.c_char_p
    library.seams_neighbours_knn.argtypes = [
        ctypes.c_void_p, ctypes.c_int, ctypes.c_double, ctypes.c_int, ctypes.c_int,
        ctypes.POINTER(ctypes.c_void_p),
    ]
    library.seams_neighbours_knn.restype = ctypes.c_int
    library.seams_neighbours_count.argtypes = [ctypes.c_void_p]
    library.seams_neighbours_count.restype = ctypes.c_int
    library.seams_neighbours_destroy.argtypes = [ctypes.c_void_p]
    library.seams_steinhardt.argtypes = [
        ctypes.c_void_p, ctypes.c_double, ctypes.c_int, ctypes.c_int,
        ctypes.POINTER(ctypes.c_double), ctypes.c_size_t,
        ctypes.POINTER(ctypes.c_double), ctypes.c_size_t,
    ]
    library.seams_steinhardt.restype = ctypes.c_int
    library.seams_chill.argtypes = [
        ctypes.c_void_p, ctypes.c_double, ctypes.c_int, ctypes.c_int, ctypes.c_int,
        ctypes.POINTER(ctypes.c_int), ctypes.c_size_t, ctypes.POINTER(ctypes.c_int),
    ]
    library.seams_chill.restype = ctypes.c_int
    library.seams_rings_find.argtypes = [
        ctypes.c_void_p, ctypes.c_double, ctypes.c_int, ctypes.c_int, ctypes.c_int,
        ctypes.POINTER(ctypes.c_void_p),
    ]
    library.seams_rings_find.restype = ctypes.c_int
    library.seams_rings_count.argtypes = [ctypes.c_void_p]
    library.seams_rings_count.restype = ctypes.c_int
    library.seams_rings_destroy.argtypes = [ctypes.c_void_p]

    import freud

    frames = {}

    def fail(handle):
        text = library.seams_handle_error(handle)
        raise RuntimeError(text.decode() if text else "seams error")

    def slot(path):
        held = frames.get(path)
        if held is None:
            handle = library.seams_system_create()
            held = {"sys": handle, "n": 0, "buf": None, "pts": None, "box": None}
            frames[path] = held
        return held

    def buffers(held):
        if held["buf"] is not None:
            return
        n = held["n"]
        held["buf"] = {
            "ql": np.empty(n, dtype=np.float64),
            "bar": np.empty(n, dtype=np.float64),
            "labels": np.empty(n, dtype=np.int32),
            "counts": np.empty(9, dtype=np.int32),
        }

    def load_points(path):
        ids = []
        coords = []
        box = []
        with open(path) as handle:
            for line in handle:
                if line.startswith("ITEM: BOX"):
                    for _ in range(3):
                        lo, hi = handle.readline().split()[:2]
                        box.append(float(hi) - float(lo))
                elif line.startswith("ITEM: ATOMS"):
                    break
            for line in handle:
                parts = line.split()
                if len(parts) < 5:
                    continue
                coords.append((float(parts[-3]), float(parts[-2]), float(parts[-1])))
        pts = np.asarray(coords, dtype=np.float64)
        return pts, tuple(box)

    def ensure(held, path):
        if held["n"] <= 0:
            if library.seams_system_read(held["sys"], path.encode(), 1, -1) != 0:
                fail(held["sys"])
            held["n"] = library.seams_system_count(held["sys"])
            held["pts"], held["box"] = load_points(path)
            buffers(held)

    sink = 0.0
    for line in sys.stdin:
        line = line.strip()
        if not line or line.startswith("quit"):
            break
        threads_s, op, path = line.split()
        threads = int(threads_s)
        os.environ["OMP_NUM_THREADS"] = str(threads)
        os.environ["PLUMED_NUM_THREADS"] = str(threads)
        if library.seams_set_num_threads(threads) != 0:
            raise RuntimeError("threads")
        freud.parallel.set_num_threads(threads)
        held = slot(path)
        if op.startswith("freud") and held["pts"] is None:
            held["pts"], held["box"] = load_points(path)
        elif op != "read":
            ensure(held, path)
            if op in ("steinhardt", "chill"):
                buffers(held)
        t0 = time.perf_counter()
        if op == "read":
            if library.seams_system_read(held["sys"], path.encode(), 1, -1) != 0:
                fail(held["sys"])
            held["n"] = library.seams_system_count(held["sys"])
            sink += held["n"]
        elif op == "knn":
            nb = ctypes.c_void_p()
            if library.seams_neighbours_knn(held["sys"], 4, CUTOFF, 0, 1, ctypes.byref(nb)) != 0:
                fail(held["sys"])
            sink += library.seams_neighbours_count(nb)
            library.seams_neighbours_destroy(nb)
        elif op == "steinhardt":
            buf = held["buf"]
            n = held["n"]
            if library.seams_steinhardt(
                held["sys"], CUTOFF, 0, 6,
                buf["ql"].ctypes.data_as(ctypes.POINTER(ctypes.c_double)), n,
                buf["bar"].ctypes.data_as(ctypes.POINTER(ctypes.c_double)), n,
            ) != 0:
                fail(held["sys"])
            sink += float(buf["ql"][0] + buf["bar"][0])
        elif op == "chill":
            buf = held["buf"]
            n = held["n"]
            if library.seams_chill(
                held["sys"], CUTOFF, 0, 0, 0,
                buf["labels"].ctypes.data_as(ctypes.POINTER(ctypes.c_int)), n,
                buf["counts"].ctypes.data_as(ctypes.POINTER(ctypes.c_int)),
            ) != 0:
                fail(held["sys"])
            sink += int(buf["counts"][1])
        elif op == "rings":
            rings = ctypes.c_void_p()
            if library.seams_rings_find(held["sys"], CUTOFF, 0, 0, 6, ctypes.byref(rings)) != 0:
                fail(held["sys"])
            sink += library.seams_rings_count(rings)
            library.seams_rings_destroy(rings)
        elif op == "freud-knn":
            lx, ly, lz = held["box"]
            box = freud.box.Box(lx, ly, lz)
            query = freud.locality.AABBQuery(box, held["pts"])
            nlist = query.query(
                held["pts"], {"num_neighbors": 4, "r_max": CUTOFF + 1.5, "exclude_ii": True}
            ).toNeighborList()
            sink += len(nlist)
        elif op == "freud-steinhardt":
            lx, ly, lz = held["box"]
            box = freud.box.Box(lx, ly, lz)
            order = freud.order.Steinhardt(l=6, average=True)
            order.compute((box, held["pts"]), neighbors={"r_max": CUTOFF})
            sink += float(order.ql[0] + order.particle_order[0])
        else:
            raise RuntimeError(op)
        sys.stdout.write("{:.9f}\n".format((time.perf_counter() - t0) * 1e3))
        sys.stdout.flush()
    sys.stderr.write("sink {:.6f}\n".format(sink))


def ovito_worker():
    import numpy as np
    from ovito.data import NearestNeighborFinder
    from ovito.io import import_file
    from ovito.modifiers import CreateBondsModifier, FindRingsModifier, SteinhardtParametersModifier
    from ovito.pipeline import Pipeline, StaticSource

    expected = int(os.environ.get("OVITO_THREAD_COUNT", "1"))
    cache = {}
    sink = 0.0

    def data_of(path):
        held = cache.get(path)
        if held is not None:
            return held
        data = import_file(path).compute()
        cache[path] = data
        return data

    for line in sys.stdin:
        line = line.strip()
        if not line or line.startswith("quit"):
            break
        threads_s, op, path = line.split()
        if int(threads_s) != expected:
            raise RuntimeError("ovito pool is {}".format(expected))
        if op != "read":
            data_of(path)
        t0 = time.perf_counter()
        if op == "read":
            data = import_file(path).compute()
            cache[path] = data
            sink += data.particles.count
        elif op == "knn":
            data = cache[path]
            finder = NearestNeighborFinder(4, data)
            neigh = finder.find_all()
            sink += int(np.size(neigh[0]))
        elif op == "steinhardt":
            data = cache[path]
            pipe = Pipeline(source=StaticSource(data=data))
            mod = SteinhardtParametersModifier(cutoff=CUTOFF, max_l=6, compute_averaged=True)
            mod.neighbor_mode = SteinhardtParametersModifier.NeighborMode.Cutoff
            pipe.modifiers.append(mod)
            result = pipe.compute()
            ql = np.asarray(result.particles["Steinhardt"])
            bar = np.asarray(result.particles["Smooth Steinhardt"])
            sink += float(ql[0, -1] + bar[0, -1])
        elif op == "rings":
            data = cache[path]
            pipe = Pipeline(source=StaticSource(data=data))
            bonds = CreateBondsModifier(cutoff=CUTOFF)
            pipe.modifiers.append(bonds)
            pipe.modifiers.append(FindRingsModifier(maximum_ring_size=6))
            result = pipe.compute()
            sink += result.particles.count
        else:
            raise RuntimeError(op)
        sys.stdout.write("{:.9f}\n".format((time.perf_counter() - t0) * 1e3))
        sys.stdout.flush()
    sys.stderr.write("sink {:.6f}\n".format(sink))


def campaign(reps, out_csv, log_path):
    lib = find_library()
    cores, model = machine()
    frames = ensure_frames(ROOT / "benchmarks" / "logs" / "binding-overhead" / "frames")
    bins = compile_drivers(lib)
    lua_exe, lua_so = find_lua_module()
    julia = shutil.which("julia")
    if not julia:
        raise SystemExit("julia is not on PATH")
    rscript = shutil.which("Rscript")
    base = child_env(lib)
    gl = Path("/usr/lib/x86_64-linux-gnu/libGL.so.1")
    ovito_env = dict(base)
    if gl.is_file():
        gl_dir = lib.parent.parent / "gl-link"
        gl_dir.mkdir(exist_ok=True)
        link = gl_dir / "libGL.so.1"
        if not link.exists():
            link.symlink_to(gl)
        ovito_env["LD_LIBRARY_PATH"] = str(gl_dir) + ":" + base["LD_LIBRARY_PATH"]
    workers = {
        "c": Worker([str(bins["c"])], base),
        "cxx": Worker([str(bins["cxx"])], base),
        "fortran": Worker([str(bins["fortran"])], base),
        "rust": Worker([str(bins["rust"])], base),
        "julia": Worker(
            [julia, "--project=" + str(ROOT / "bindings" / "julia" / "Seams"),
             str(ROOT / "benchmarks" / "binding_overhead.jl")],
            base,
        ),
        "r": Worker([rscript, str(ROOT / "benchmarks" / "binding_overhead.R")],
                    {**base, "SEAMS_R_SO": str(bins["rso"])}),
        "lua": Worker([str(lua_exe), str(ROOT / "benchmarks" / "binding_overhead.lua")],
                      {**base, "LUA_CPATH": str(lua_so.parent / "?.so") + ";;"}),
        "python": Worker([sys.executable, str(Path(__file__).resolve()), "--worker", "python"], base),
        "ovito-1": Worker(
            [sys.executable, str(Path(__file__).resolve()), "--worker", "ovito"],
            {**ovito_env, "OVITO_THREAD_COUNT": "1", "QT_QPA_PLATFORM": "offscreen"},
        ),
        "ovito-n": Worker(
            [sys.executable, str(Path(__file__).resolve()), "--worker", "ovito"],
            {**ovito_env, "OVITO_THREAD_COUNT": str(cores), "QT_QPA_PLATFORM": "offscreen"},
        ),
    }
    log_path.parent.mkdir(parents=True, exist_ok=True)
    rows = []
    try:
        with log_path.open("w") as log:
            log.write("cpu_model\t{}\ncpu_cores\t{}\n".format(model, cores))
            for threads in (1, cores):
                plan = []
                for binding in POINTER_BINDINGS:
                    for molecules in SIZES:
                        for op in OPS:
                            if op == "rings":
                                continue
                            plan.append((binding, binding, op, molecules))
                # Prism classification inside seams_rings_find dominates above
                # 4096. Every binding is timed on the 4096-molecule cell. The
                # larger cells are the C search those wrappers call.
                for binding in POINTER_BINDINGS:
                    for molecules in SIZES:
                        # findPrisms inside seams_rings_find is superlinear.
                        # 4096 molecules is about one second. 32768 was about
                        # a minute in a probe, and 110592 had not returned
                        # after several minutes. Every binding is repeated at
                        # 4096, where the search still dominates the wrapper.
                        if molecules != 4096:
                            continue
                        plan.append((binding, binding, "rings", molecules))
                for molecules in SIZES:
                    for op in GUEST_OPS["freud"]:
                        plan.append(("freud", "python", "freud-" + op, molecules))
                    for op in GUEST_OPS["ovito"]:
                        if op == "rings" and molecules > 4096:
                            continue
                        which = "ovito-1" if threads == 1 else "ovito-n"
                        plan.append(("ovito", which, op, molecules))
                samples = {}
                for binding, worker, op, molecules in plan:
                    public = op.split("-", 1)[-1] if op.startswith("freud-") else op
                    rel = frames[molecules].relative_to(ROOT).as_posix()
                    sys.stderr.write("warmup {} {} {} {}\n".format(worker, op, molecules, threads))
                    sys.stderr.flush()
                    workers[worker].ask(threads, op, rel)
                    values = []
                    for _rep in range(reps):
                        ms = workers[worker].ask(threads, op, rel)
                        values.append(ms)
                        log.write("{}\t{}\t{}\t{}\t{:.9f}\n".format(
                            binding, public, molecules, threads, ms))
                        log.flush()
                    samples[(binding, public, molecules)] = values
                for (binding, op, molecules), values in samples.items():
                    median, iqr, q1, q3 = summarize(values)
                    rows.append({
                        "cpu_model": model,
                        "cpu_cores": cores,
                        "binding": binding,
                        "operation": op,
                        "molecules": molecules,
                        "threads": threads,
                        "n": len(values),
                        "median_ms": "{:.6f}".format(median),
                        "iqr_ms": "{:.6f}".format(iqr),
                        "q1_ms": "{:.6f}".format(q1),
                        "q3_ms": "{:.6f}".format(q3),
                    })
    finally:
        for worker in workers.values():
            worker.close()
    out_csv.parent.mkdir(parents=True, exist_ok=True)
    with out_csv.open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(rows[0]))
        writer.writeheader()
        writer.writerows(rows)
    return rows


def main(argv=None):
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--out", type=Path, default=ROOT / "benchmarks" / "binding-overhead.csv")
    parser.add_argument("--log", type=Path,
                        default=ROOT / "benchmarks" / "logs" / "binding-overhead" / "samples.txt")
    parser.add_argument("--reps", type=int, default=8)
    parser.add_argument("--worker", choices=("python", "ovito"))
    args = parser.parse_args(argv)
    if args.worker == "python":
        python_worker(find_library())
        return 0
    if args.worker == "ovito":
        ovito_worker()
        return 0
    if args.reps < 8:
        raise SystemExit("reps must be at least 8")
    campaign(args.reps, args.out, args.log)
    render = ROOT / "docs" / "render_binding_overhead.py"
    subprocess.run([sys.executable, str(render)], cwd=ROOT, check=True)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
