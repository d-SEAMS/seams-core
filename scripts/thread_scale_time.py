#!/usr/bin/env python3
"""Time CHILL+, q6, rings and neighbours at 1, 2, 4 and 8 threads.

Release -O3 only. One warm-up, six timed calls, median and type-7 IQR.
d-SEAMS takes its width from OMP_NUM_THREADS and gives that width to
Rayon. mdapy uses its own pool. OVITO reads OVITO_THREAD_COUNT before
import. freud uses freud.parallel.set_num_threads.
"""

from __future__ import annotations

import csv
import os
import platform
import statistics
import subprocess
import sys
import tempfile
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
_TMP = Path(tempfile.gettempdir())
H2H_DSEAMS = str(_TMP / "h2h_dseams")
MDAPY_PYTHON = str(_TMP / "mdapy-venv" / "bin" / "python")
CHILL_REL_LIB = str(
    _TMP / "builds" / "chill-rel" / "bindings" / "libseams_lang.so"
)
sys.path.insert(0, str(ROOT / "scripts"))

import release_build

THREADS = (1, 2, 4, 8)
SIZES = (11, 24)  # 10648 and 110592
WARMUP = 1
REPEATS = 6
FLEXI = (
    "/cvmfs/software.eessi.io/versions/2026.06/software/linux/x86_64/intel/"
    "sapphirerapids/software/FlexiBLAS/3.5.0-GCC-15.2.0/lib"
)


def lib_env(lib: str) -> str:
    build = Path(lib).resolve().parents[1]
    parts = [
        str(build / "bindings"),
        str(build / "src"),
        str(build / "subprojects/readcon-core"),
        str(build / "subprojects/readcon-db"),
        str(build / "subprojects/linkcell"),
        str(build / "subprojects/vesin"),
        str(build / "subprojects/minimage"),
        FLEXI,
    ]
    return ":".join(parts)


def summarize(samples):
    ordered = sorted(samples)
    med = statistics.median(ordered)
    q1, _, q3 = statistics.quantiles(ordered, n=4, method="inclusive")
    return med, q1, q3, q3 - q1


def hardware():
    nproc = subprocess.check_output(["nproc"], text=True).strip()
    nproc_all = subprocess.check_output(["nproc", "--all"], text=True).strip()
    cpu = platform.processor() or platform.machine()
    try:
        for line in Path("/proc/cpuinfo").read_text().splitlines():
            if line.startswith("model name"):
                cpu = line.split(":", 1)[1].strip()
                break
    except OSError:
        pass
    allowed = "unknown"
    try:
        for line in Path("/proc/self/status").read_text().splitlines():
            if line.startswith("Cpus_allowed_list:"):
                allowed = line.split(":", 1)[1].strip()
    except OSError:
        pass
    return f"{cpu}; nproc={nproc}; nproc_all={nproc_all}; affinity={allowed}"


def seams_env(threads, lib):
    env = os.environ.copy()
    env.pop("LD_PRELOAD", None)
    env["LD_LIBRARY_PATH"] = lib_env(lib)
    env["OMP_NUM_THREADS"] = str(threads)
    env["OMP_DYNAMIC"] = "FALSE"
    env["OMP_PROC_BIND"] = "close"
    env["RAYON_NUM_THREADS"] = str(threads)
    env["OPENBLAS_NUM_THREADS"] = "1"
    env["MKL_NUM_THREADS"] = "1"
    env["SEAMS_LIB"] = lib
    return env


def run_seams(threads, molecules_m, mode, lib):
    # h2h_dseams warms up once, then prints one line per timed call.
    cmd = ["taskset", "-c", "0-7", H2H_DSEAMS, str(molecules_m), str(REPEATS), mode]
    done = subprocess.run(cmd, env=seams_env(threads, lib), text=True, capture_output=True, check=True)
    return [float(line) for line in done.stdout.split() if line.strip()]


def py_env():
    env = os.environ.copy()
    env.pop("LD_PRELOAD", None)
    env["LD_LIBRARY_PATH"] = "/usr/lib/x86_64-linux-gnu"
    env["QT_QPA_PLATFORM"] = "offscreen"
    env["FASTATOMSTRUCT_NO_TELEMETRY"] = "1"
    env["OPENBLAS_NUM_THREADS"] = "1"
    env["MKL_NUM_THREADS"] = "1"
    env["SEAMS_SCRIPTS"] = str(ROOT / "scripts")
    return env


OTHER = r"""
import json, math, os, sys, time
sys.path.insert(0, os.environ["SEAMS_SCRIPTS"])
import h2h_protocol as proto
import numpy as np
kind, m, threads, repeats = sys.argv[1], int(sys.argv[2]), int(sys.argv[3]), int(sys.argv[4])
pos, side = proto.ice_ic(m)

def samples(prep, call):
    call(prep())
    out = []
    for _ in range(repeats):
        frame = prep()
        t0 = time.perf_counter()
        call(frame)
        out.append((time.perf_counter() - t0) * 1000.0)
    return out

if kind == "mdapy":
    import mdapy
    from mdapy import System
    from mdapy.identify import ChillPlus
    mdapy.set_num_threads(threads)
    if mdapy.get_num_threads() != threads:
        raise SystemExit(f"mdapy threads {mdapy.get_num_threads()} != {threads}")
    def prep():
        return System(pos, box=side)
    def call(system):
        return ChillPlus(cutoff=3.5).compute(system)
    ver = mdapy.__version__
elif kind == "ovito":
    os.environ["OVITO_THREAD_COUNT"] = str(threads)
    import ovito
    from ovito.data import DataCollection
    from ovito.modifiers import ChillPlusModifier
    from ovito.pipeline import Pipeline, StaticSource
    if tuple(ovito.version)[:3] != (3, 17, 2):
        raise SystemExit(f"ovito {ovito.version}")
    def prep():
        data = DataCollection()
        data.create_cell(np.diag([side, side, side]), pbc=(True, True, True))
        parts = data.create_particles(count=len(pos))
        parts.create_property("Position", data=pos)
        return data
    def call(data):
        pipe = Pipeline(source=StaticSource(data=data))
        pipe.modifiers.append(ChillPlusModifier(cutoff=3.5))
        return pipe.compute()
    ver = ".".join(str(v) for v in tuple(ovito.version)[:3])
elif kind == "freud":
    import freud
    freud.parallel.set_num_threads(threads)
    box = freud.box.Box.cube(side)
    pts = pos - side / 2.0
    query = freud.locality.AABBQuery(box, pts)
    def prep():
        return query
    def call(q):
        return q.query(pts, {"r_max": 3.5, "exclude_ii": True}).toNeighborList()
    ver = freud.__version__
else:
    raise SystemExit(kind)
print(ver)
for v in samples(prep, call):
    print(f"{v:.6f}")
"""


def run_other(kind, threads, molecules_m):
    # OVITO reads its thread count at import, so each width is a new process.
    env = py_env()
    if kind == "ovito":
        env["OVITO_THREAD_COUNT"] = str(threads)
    done = subprocess.run(
        ["taskset", "-c", "0-7", MDAPY_PYTHON, "-c", OTHER,
         kind, str(molecules_m), str(threads), str(REPEATS)],
        env=env, text=True, capture_output=True, check=False,
    )
    if done.returncode != 0:
        raise RuntimeError(done.stderr)
    lines = [ln for ln in done.stdout.splitlines() if ln.strip()]
    version = lines[0]
    samples = [float(ln) for ln in lines[1:]]
    return version, samples


def _lang_env(lib: str, threads: int) -> None:
    os.environ.pop("LD_PRELOAD", None)
    os.environ["LD_LIBRARY_PATH"] = lib_env(lib)
    os.environ["OMP_NUM_THREADS"] = str(threads)
    os.environ["OMP_DYNAMIC"] = "FALSE"
    os.environ["OMP_PROC_BIND"] = "close"
    os.environ["RAYON_NUM_THREADS"] = str(threads)
    os.environ["OPENBLAS_NUM_THREADS"] = "1"
    os.environ["MKL_NUM_THREADS"] = "1"
    try:
        os.sched_setaffinity(0, set(range(threads)))
    except (AttributeError, OSError):
        pass


def _ice_frame(lib, cells: int):
    import ctypes
    import numpy as np
    import h2h_protocol as proto

    pos, side = proto.ice_ic(cells)
    n = len(pos)
    xyz = np.ascontiguousarray(pos, dtype=np.float64)
    box = np.array([0, 0, 0, side, side, side, 0, 0, 0], dtype=np.float64)
    types = np.ones(n, dtype=np.int32)
    ids = np.arange(1, n + 1, dtype=np.int32)
    handle = lib.seams_frame_from_arrays(
        xyz.ctypes.data_as(ctypes.POINTER(ctypes.c_double)), n,
        box.ctypes.data_as(ctypes.POINTER(ctypes.c_double)),
        types.ctypes.data_as(ctypes.POINTER(ctypes.c_int)),
        ids.ctypes.data_as(ctypes.POINTER(ctypes.c_int)))
    if not handle:
        raise RuntimeError("could not build the ice frame")
    return handle, n, xyz, box


def _bind_lang(lib):
    import ctypes
    lib.seams_last_error.restype = ctypes.c_char_p
    lib.seams_frame_from_arrays.restype = ctypes.c_void_p
    lib.seams_frame_from_arrays.argtypes = [
        ctypes.POINTER(ctypes.c_double), ctypes.c_int,
        ctypes.POINTER(ctypes.c_double),
        ctypes.POINTER(ctypes.c_int), ctypes.POINTER(ctypes.c_int)]
    lib.seams_frame_free.argtypes = [ctypes.c_void_p]
    lib.seams_frame_steinhardt.restype = ctypes.c_int
    lib.seams_frame_steinhardt.argtypes = [
        ctypes.c_void_p, ctypes.c_double, ctypes.c_int, ctypes.c_int,
        ctypes.POINTER(ctypes.c_double), ctypes.POINTER(ctypes.c_double)]
    lib.seams_frame_q4_q6.restype = ctypes.c_int
    lib.seams_frame_q4_q6.argtypes = [
        ctypes.c_void_p, ctypes.c_double, ctypes.c_int,
        ctypes.POINTER(ctypes.c_double), ctypes.POINTER(ctypes.c_double),
        ctypes.c_void_p, ctypes.c_void_p]
    return lib


def _timed(call, warmup: int, repeats: int) -> list[float]:
    import time
    samples = []
    for index in range(warmup + repeats):
        t0 = time.perf_counter()
        call()
        elapsed = (time.perf_counter() - t0) * 1000.0
        if index >= warmup:
            samples.append(elapsed)
    return samples


def write_paper_q(out_dir: Path, repeats: int, threads: int, lang_lib: str,
                  core_lib: str, protocol: str) -> None:
    """Three q6 protocols. They are not one number.

    tab:ql-eq is seams_frame_steinhardt of order 6 with the averaged
    buffer filled, neighbour build inside the call. That is the driver
    behind the filled 4.005 ms q4 cell. tab:cli is seams_frame_q4_q6
    after the frame is loaded, q4 and q6 together, averaged buffers
    null. tab:lit is seams_ql_cutoff on the core library, the call whose
    q4 median on this size is 25.0 ms.
    """
    import ctypes
    import numpy as np
    sys.path.insert(0, str(ROOT / "scripts"))
    import paper_stats as stats
    import release_build

    if repeats < 8:
        raise SystemExit("paper cells need n >= 8")
    cells = 24
    molecules = 8 * cells ** 3
    out_dir.mkdir(parents=True, exist_ok=True)
    # Record the mask before libgomp initialises. OpenMP's proc bind can
    # narrow the main thread to one place while the team still uses the rest.
    cores, affinity = stats.cores_now()
    _lang_env(lang_lib, threads)
    _label = release_build.assert_release(Path(lang_lib).resolve().parents[1])
    lang = _bind_lang(ctypes.CDLL(lang_lib))
    frame, n, xyz, box = _ice_frame(lang, cells)
    if n != molecules:
        raise SystemExit(f"ice Ic has {n} molecules")
    if cores < threads:
        raise SystemExit(f"affinity has {cores} cores, need {threads}")
    log = [
        f"cores {cores}",
        f"affinity {affinity}",
        f"threads {threads}",
        f"molecules {molecules}",
        "warmup 1",
        f"n {repeats}",
        f"build {_label}",
        "ql-eq seams_frame_steinhardt order 6, ql and qlbar, neighbour inside",
        "cli seams_frame_q4_q6 after the frame is loaded, q4 and q6 together",
        "lit seams_ql_cutoff order 6, neighbour inside the core call",
    ]

    def write_csv(path: Path, rows: list[dict]) -> None:
        with path.open("w", newline="") as handle:
            writer = csv.DictWriter(handle, fieldnames=list(rows[0]))
            writer.writeheader()
            writer.writerows(rows)

    def one(cell: str, code: str, samples: list[float]) -> dict:
        summary = stats.summarize(samples)
        log.append(
            f"{cell} median_ms {summary['median_ms']:.6f} "
            f"iqr_ms {summary['iqr_ms']:.6f} samples {stats.fmt_samples(samples)}"
        )
        print(f"{cell} {summary['median_ms']:.4f} ms iqr {summary['iqr_ms']:.4f}", flush=True)
        return {
            "cell": cell,
            "code": code,
            "molecules": molecules,
            "threads": threads,
            "cores": cores,
            "affinity": affinity,
            "warmup": WARMUP,
            "n": summary["n"],
            "median_ms": f"{summary['median_ms']:.6f}",
            "iqr_ms": f"{summary['iqr_ms']:.6f}",
            "q1_ms": f"{summary['q1_ms']:.6f}",
            "q3_ms": f"{summary['q3_ms']:.6f}",
            "samples_ms": stats.fmt_samples(samples),
        }

    if protocol == "ql-eq":
        ql = (ctypes.c_double * n)()
        bar = (ctypes.c_double * n)()

        def ql_eq():
            if lang.seams_frame_steinhardt(frame, 3.5, 1, 6, ql, bar) != 0:
                raise RuntimeError(lang.seams_last_error().decode())

        eq = one("q6 eight threads", "d-SEAMS", _timed(ql_eq, WARMUP, repeats))
        write_csv(out_dir / "tab-ql-eq.csv", [eq])
    if protocol == "cli":
        q4 = (ctypes.c_double * n)()
        q6 = (ctypes.c_double * n)()

        def both():
            if lang.seams_frame_q4_q6(frame, 3.5, 1, q4, q6, None, None) != 0:
                raise RuntimeError(lang.seams_last_error().decode())

        cli = one("library q4, q6", "d-SEAMS", _timed(both, WARMUP, repeats))
        write_csv(out_dir / "frag-cli-q.csv", [{
            **cli,
            "detail": "library call after the frame is loaded; q4 and q6 together",
        }])
    lang.seams_frame_free(frame)
    if protocol != "lit":
        (out_dir / "raw-q6.txt").write_text(
            ((out_dir / "raw-q6.txt").read_text(encoding="utf-8") if (out_dir / "raw-q6.txt").is_file() else "")
            + "\n".join(log) + "\n",
            encoding="utf-8",
        )
        return

    core = ctypes.CDLL(core_lib)
    core.seams_ql_cutoff.restype = ctypes.c_int
    core.seams_ql_cutoff.argtypes = [
        ctypes.POINTER(ctypes.c_double), ctypes.c_int,
        ctypes.POINTER(ctypes.c_double), ctypes.c_double,
        ctypes.c_int, ctypes.POINTER(ctypes.c_double)]
    values = np.empty(n, dtype=np.float64)
    xyz_ptr = xyz.ctypes.data_as(ctypes.POINTER(ctypes.c_double))
    box_ptr = box.ctypes.data_as(ctypes.POINTER(ctypes.c_double))
    ql_ptr = values.ctypes.data_as(ctypes.POINTER(ctypes.c_double))

    def lit():
        code = core.seams_ql_cutoff(xyz_ptr, n, box_ptr, 3.5, 6, ql_ptr)
        if code != 0:
            raise RuntimeError(f"seams_ql_cutoff returned {code}")

    lit_row = one("q6 eight threads", "d-SEAMS", _timed(lit, WARMUP, repeats))
    lit_row["detail"] = "seams_ql_cutoff; the q4 cell of this table is the 25.0 ms protocol"
    write_csv(out_dir / "frag-lit-q6.csv", [lit_row])
    previous = (out_dir / "raw-q6.txt").read_text(encoding="utf-8") if (out_dir / "raw-q6.txt").is_file() else ""
    (out_dir / "raw-q6.txt").write_text(previous + "\n".join(log) + "\n", encoding="utf-8")


def main():
    global REPEATS
    import argparse
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--repeats", type=int, default=REPEATS)
    parser.add_argument("--paper-cells", type=Path, default=None)
    parser.add_argument("--threads", type=int, default=8)
    parser.add_argument(
        "--lang-lib",
        default=str(ROOT / "bbdir-release/bindings/libseams_lang.so"))
    parser.add_argument(
        "--core-lib",
        default=str(ROOT / "bbdir-release/src/libyodaLib.so"))
    parser.add_argument("--protocol", choices=("ql-eq", "cli", "lit"))
    args, _unknown = parser.parse_known_args()
    REPEATS = args.repeats
    if args.paper_cells is not None and args.protocol is None:
        raw = args.paper_cells / "raw-q6.txt"
        args.paper_cells.mkdir(parents=True, exist_ok=True)
        raw.write_text("", encoding="utf-8")
        cores = f"0-{args.threads - 1}" if args.threads > 1 else "0"
        for protocol in ("ql-eq", "cli", "lit"):
            cmd = [
                "taskset", "-c", cores, sys.executable, str(Path(__file__).resolve()),
                "--paper-cells", str(args.paper_cells),
                "--repeats", str(args.repeats),
                "--threads", str(args.threads),
                "--lang-lib", args.lang_lib,
                "--core-lib", args.core_lib,
                "--protocol", protocol,
            ]
            subprocess.run(cmd, check=True)
        return
    if args.paper_cells is not None:
        write_paper_q(
            args.paper_cells, args.repeats, args.threads,
            args.lang_lib, args.core_lib, args.protocol)
        return
    lib, build = release_build.ensure_release(
        CHILL_REL_LIB)
    hw = hardware()
    out = ROOT / "benchmarks" / "thread-scale.csv"
    fields = [
        "hardware", "build", "code", "version", "threads", "operation",
        "molecules", "warmup", "n", "median_ms", "q1_ms", "q3_ms", "iqr_ms",
        "samples_ms",
    ]
    rows = []
    jobs = []
    for m in SIZES:
        molecules = 8 * m ** 3
        for threads in THREADS:
            for mode, op in (
                ("chill", "published"),
                ("ovito", "ovito-chill+"),
                ("q6", "q6"),
                ("rings", "rings"),
                ("neigh", "neighbours"),
            ):
                jobs.append(("dseams", threads, m, molecules, mode, op))
            jobs.append(("mdapy", threads, m, molecules, "mdapy", "chill+"))
            jobs.append(("ovito", threads, m, molecules, "ovito", "chill+"))
            jobs.append(("freud", threads, m, molecules, "freud", "neighbours"))
    for threads in (1, 8):
        molecules = 8 * 48 ** 3
        for mode, op in (
            ("chill", "published"),
            ("ovito", "ovito-chill+"),
            ("q6", "q6"),
            ("rings", "rings"),
            ("neigh", "neighbours"),
        ):
            jobs.append(("dseams", threads, 48, molecules, mode, op))
        jobs.append(("mdapy", threads, 48, molecules, "mdapy", "chill+"))
        jobs.append(("ovito", threads, 48, molecules, "ovito", "chill+"))
        jobs.append(("freud", threads, 48, molecules, "freud", "neighbours"))
    with out.open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields)
        writer.writeheader()
        for kind, threads, m, molecules, mode, op in jobs:
            print(f"{kind} threads={threads} m={m} {op}", flush=True)
            if kind == "dseams":
                samples = run_seams(threads, m, mode, lib)
                version = "2.14.0"
                code = "d-SEAMS"
            else:
                version, samples = run_other(mode, threads, m)
                code = {"mdapy": "mdapy", "ovito": "OVITO", "freud": "freud"}[kind]
            if len(samples) < REPEATS:
                raise SystemExit(f"short sample {kind} {threads} {m}: {samples}")
            med, q1, q3, iqr = summarize(samples[:REPEATS])
            writer.writerow({
                "hardware": hw,
                "build": build,
                "code": code,
                "version": version,
                "threads": threads,
                "operation": op,
                "molecules": molecules,
                "warmup": WARMUP,
                "n": REPEATS,
                "median_ms": f"{med:.6f}",
                "q1_ms": f"{q1:.6f}",
                "q3_ms": f"{q3:.6f}",
                "iqr_ms": f"{iqr:.6f}",
                "samples_ms": ",".join(f"{v:.6f}" for v in samples[:REPEATS]),
            })
            handle.flush()
            print(f"  {med:.3f} ms  iqr {iqr:.3f}", flush=True)
    print(f"wrote {out}", flush=True)


if __name__ == "__main__":
    main()
