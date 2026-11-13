#!/usr/bin/env python3
"""rings-statistics-matter 0.3.0 against d-SEAMS at 110592 molecules.

tab:rings-214 times the library call in process: two warm-ups, then
eight calls, at one thread and at eight. Each width is its own process
so OpenMP and Rayon see the width at startup. The eight-thread
rings-statistics-matter median is that cell. tab:cli times a fresh
process per sample, so the sample includes process start: one warm-up,
then eight processes. Both must count 221184 six-rings. Cutoff 3.5 A,
maximum ring size 6. The rival uses maxlvl 8 because a six-ring
midpoint sits at distance 3.
"""

from __future__ import annotations

import argparse
import csv
import subprocess
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "scripts"))
import paper_stats as stats

MOLECULES = 110592
SIX = 221184
CELLS = 24


def ice():
    import h2h_protocol as proto
    pos, side = proto.ice_ic(CELLS)
    if len(pos) != MOLECULES:
        raise RuntimeError(f"ice Ic has {len(pos)} molecules")
    return pos, side


def six_rsm(counts) -> int:
    import numpy as np
    values = np.asarray(counts, dtype=np.float64)
    if values.shape[0] < 6:
        raise RuntimeError(f"short ring vector {values.shape}")
    six = float(values[5])
    if abs(six - round(six)) > 1e-6:
        raise RuntimeError(f"non-integer six-ring count {six}")
    return int(round(six))


def atoms_of():
    import numpy as np
    from ase import Atoms
    pos, side = ice()
    return Atoms(
        symbols=["O"] * MOLECULES,
        positions=np.ascontiguousarray(pos),
        cell=[side, side, side],
        pbc=True,
    )


def call_rsm(atoms, threads: int) -> int:
    from rings_statistics_matter import ring_statistics
    counts, _rings = ring_statistics(
        atoms, cutoff=3.5, maxlvl=8, supercell=False, deduplicate=True,
        backend="rust", num_threads=threads)
    return six_rsm(counts)


def samples_rsm(threads: int, warmup: int, repeats: int) -> list[float]:
    import time
    from importlib.metadata import version
    found = version("rings-statistics-matter")
    if found != "0.3.0":
        raise RuntimeError(f"rings-statistics-matter {found}, need 0.3.0")
    atoms = atoms_of()
    samples = []
    for index in range(warmup + repeats):
        t0 = time.perf_counter()
        six = call_rsm(atoms, threads)
        elapsed = (time.perf_counter() - t0) * 1000.0
        if six != SIX:
            raise RuntimeError(f"rsm six-rings {six} != {SIX}")
        if index >= warmup:
            samples.append(elapsed)
    return samples


def samples_dseams(lib_path: str, threads: int, warmup: int, repeats: int) -> list[float]:
    import ctypes
    import time
    import numpy as np

    lib = ctypes.CDLL(lib_path)
    lib.seams_last_error.restype = ctypes.c_char_p
    lib.seams_frame_from_arrays.restype = ctypes.c_void_p
    lib.seams_frame_from_arrays.argtypes = [
        ctypes.POINTER(ctypes.c_double), ctypes.c_int,
        ctypes.POINTER(ctypes.c_double),
        ctypes.POINTER(ctypes.c_int), ctypes.POINTER(ctypes.c_int)]
    lib.seams_frame_free.argtypes = [ctypes.c_void_p]
    lib.seams_frame_rings.restype = ctypes.c_int
    lib.seams_frame_rings.argtypes = [
        ctypes.c_void_p, ctypes.c_double, ctypes.c_int, ctypes.c_char_p,
        ctypes.c_int, ctypes.POINTER(ctypes.c_int), ctypes.c_int]
    pos, side = ice()
    xyz = np.ascontiguousarray(pos, dtype=np.float64)
    n = len(xyz)
    box = np.array([0, 0, 0, side, side, side, 0, 0, 0], dtype=np.float64)
    types = np.ones(n, dtype=np.int32)
    ids = np.arange(1, n + 1, dtype=np.int32)
    handle = lib.seams_frame_from_arrays(
        xyz.ctypes.data_as(ctypes.POINTER(ctypes.c_double)), n,
        box.ctypes.data_as(ctypes.POINTER(ctypes.c_double)),
        types.ctypes.data_as(ctypes.POINTER(ctypes.c_int)),
        ids.ctypes.data_as(ctypes.POINTER(ctypes.c_int)))
    if not handle:
        raise RuntimeError(lib.seams_last_error().decode())
    counts = (ctypes.c_int * 7)()
    samples = []
    try:
        for index in range(warmup + repeats):
            t0 = time.perf_counter()
            status = lib.seams_frame_rings(handle, 3.5, 1, b"cutoff", 6, counts, 7)
            elapsed = (time.perf_counter() - t0) * 1000.0
            if status != 0:
                raise RuntimeError(lib.seams_last_error().decode())
            if int(counts[6]) != SIX:
                raise RuntimeError(f"d-SEAMS six-rings {int(counts[6])} != {SIX}")
            if index >= warmup:
                samples.append(elapsed)
    finally:
        lib.seams_frame_free(handle)
    return samples


def emit_samples(samples: list[float]) -> None:
    print(" ".join(f"{value:.6f}" for value in samples))


def lang_library_path(lib: str) -> str:
    import os
    build = Path(lib).resolve().parents[1]
    parts = [
        build / "bindings",
        build / "src",
        build / "subprojects/readcon-core",
        build / "subprojects/readcon-db",
        build / "subprojects/linkcell",
        build / "subprojects/vesin",
        build / "subprojects/minimage",
    ]
    return os.pathsep.join(str(path) for path in parts)


def run_worker(kind: str, threads: int, lib: str, warmup: int, repeats: int) -> list[float]:
    cmd = [
        "taskset", "-c", stats.taskset_cores(threads),
        sys.executable, str(Path(__file__).resolve()),
        "--worker", kind, "--threads", str(threads),
        "--lib", lib, "--warmup", str(warmup), "--repeats", str(repeats),
    ]
    env = stats.pin_env(threads)
    if kind == "dseams":
        env["LD_LIBRARY_PATH"] = lang_library_path(lib)
    done = subprocess.run(
        cmd, cwd=ROOT, env=env, text=True, capture_output=True)
    if done.returncode != 0:
        tail = (done.stderr or done.stdout or "").strip().splitlines()
        raise RuntimeError(tail[-1] if tail else f"exit {done.returncode}")
    line = [row for row in done.stdout.splitlines() if row.strip()][-1]
    samples = [float(piece) for piece in line.split()]
    if len(samples) != repeats:
        raise RuntimeError(f"{kind} threads {threads} returned {len(samples)} samples")
    return samples


def command_samples(threads: int, warmup: int, repeats: int, lib: str) -> list[float]:
    import time
    samples = []
    cmd = [
        "taskset", "-c", stats.taskset_cores(threads),
        sys.executable, str(Path(__file__).resolve()),
        "--worker", "rsm-once", "--threads", str(threads), "--lib", lib,
    ]
    env = stats.pin_env(threads)
    for index in range(warmup + repeats):
        t0 = time.perf_counter()
        done = subprocess.run(cmd, cwd=ROOT, env=env, text=True, capture_output=True)
        elapsed = (time.perf_counter() - t0) * 1000.0
        if done.returncode != 0:
            tail = (done.stderr or done.stdout or "").strip().splitlines()
            raise RuntimeError(tail[-1] if tail else f"exit {done.returncode}")
        if done.stdout.strip() != str(SIX):
            raise RuntimeError(f"command census {done.stdout.strip()!r}")
        if index >= warmup:
            samples.append(elapsed)
    return samples


def write_table(path: Path, rows: list[dict]) -> None:
    with path.open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(rows[0]))
        writer.writeheader()
        writer.writerows(rows)


def record(rows, log, cell, code, threads, warmup, samples, protocol, cores, affinity):
    summary = stats.summarize(samples)
    rows.append({
        "cell": cell,
        "code": code,
        "protocol": protocol,
        "molecules": MOLECULES,
        "threads": threads,
        "cores": cores,
        "affinity": affinity,
        "warmup": warmup,
        "n": summary["n"],
        "median_ms": f"{summary['median_ms']:.6f}",
        "median_3sig": f"{summary['median_ms']:.3g}",
        "iqr_ms": f"{summary['iqr_ms']:.6f}",
        "q1_ms": f"{summary['q1_ms']:.6f}",
        "q3_ms": f"{summary['q3_ms']:.6f}",
        "six_rings": SIX,
        "samples_ms": stats.fmt_samples(samples),
    })
    log.append(
        f"{protocol} {code} threads {threads} median_ms {summary['median_ms']:.6f} "
        f"iqr_ms {summary['iqr_ms']:.6f} six {SIX} samples {stats.fmt_samples(samples)}"
    )
    print(f"{protocol} {code} threads {threads} {summary['median_ms']:.4f} ms", flush=True)


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--out", type=Path, default=ROOT / "benchmarks/paper-cells")
    parser.add_argument("--lib", default=str(ROOT / "bbdir-release/bindings/libseams_lang.so"))
    parser.add_argument("--worker", choices=("rsm", "dseams", "rsm-once"))
    parser.add_argument("--threads", type=int, default=8)
    parser.add_argument("--warmup", type=int, default=2)
    parser.add_argument("--repeats", type=int, default=8)
    args = parser.parse_args()
    if args.worker == "rsm-once":
        import os
        os.environ["RAYON_NUM_THREADS"] = str(args.threads)
        six = call_rsm(atoms_of(), args.threads)
        if six != SIX:
            raise SystemExit(f"rsm six-rings {six} != {SIX}")
        print(six)
        return 0
    if args.worker == "rsm":
        emit_samples(samples_rsm(args.threads, args.warmup, args.repeats))
        return 0
    if args.worker == "dseams":
        emit_samples(samples_dseams(args.lib, args.threads, args.warmup, args.repeats))
        return 0
    if args.repeats < 8:
        raise SystemExit("n must be at least 8")
    if not Path(args.lib).is_file():
        raise SystemExit(f"missing {args.lib}")
    args.out.mkdir(parents=True, exist_ok=True)
    cores, affinity = stats.cores_now()
    log = [
        *stats.machine_lines(ROOT),
        f"cores {cores}",
        f"affinity {affinity}",
        "cutoff_A 3.5",
        "max_ring 6",
        "rsm 0.3.0 maxlvl 8",
        f"six_rings {SIX}",
        "in-process warmup 2; command warmup 1; n 8",
        "command time includes process start",
    ]
    in_rows = []
    command_rows = []
    for threads in (1, 8):
        samples = run_worker("dseams", threads, args.lib, 2, args.repeats)
        record(in_rows, log, "rings-statistics-matter 0.3.0", "d-SEAMS",
               threads, 2, samples, "in-process", cores, affinity)
        samples = run_worker("rsm", threads, args.lib, 2, args.repeats)
        record(in_rows, log, "rings-statistics-matter 0.3.0",
               "rings-statistics-matter", threads, 2, samples, "in-process",
               cores, affinity)
        samples = command_samples(threads, 1, args.repeats, args.lib)
        record(command_rows, log, "rings-statistics-matter 0.3.0",
               "rings-statistics-matter", threads, 1, samples, "command",
               cores, affinity)
    (args.out / "raw-rings.txt").write_text("\n".join(log) + "\n", encoding="utf-8")
    write_table(args.out / "tab-rings-214.csv", in_rows)
    cli = []
    for row in command_rows:
        cli.append({
            "cell": "rings-statistics-matter 0.3.0" if int(row["threads"]) == 8
            else "rings-statistics-matter 0.3.0, 1 thread",
            "code": row["code"],
            "molecules": row["molecules"],
            "threads": row["threads"],
            "cores": row["cores"],
            "affinity": row["affinity"],
            "warmup": row["warmup"],
            "n": row["n"],
            "median_ms": row["median_ms"],
            "iqr_ms": row["iqr_ms"],
            "samples_ms": row["samples_ms"],
            "detail": "command time includes process start; 221184 six-rings",
        })
    write_table(args.out / "frag-cli-rsm.csv", cli)
    return 0


if __name__ == "__main__":
    sys.exit(main())
