#!/usr/bin/env python3
"""q4, q6 and cutoff neighbour lists on ideal ice Ic.

q4 and q6 at 110592 molecules against freud Steinhardt, pyscal3 and
OVITO SteinhardtParametersModifier. Neighbour lists at 10648 and 110592
against freud AABBQuery, matscipy neighbour_list, ASE neighbor_list and
MDAnalysis capped_distance. Cutoff 3.5 A. One warm-up, six calls, median
and type-7 IQR. One core is taskset -c 0. Eight cores are taskset -c 0-7.
"""

from __future__ import annotations

import argparse
import csv
import json
import os
import subprocess
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
HERE = Path(__file__).resolve()
HOST_LIB = "/usr/lib/x86_64-linux-gnu"
FLEXI = (
    "/cvmfs/software.eessi.io/versions/2026.06/software/linux/x86_64/intel/"
    "sapphirerapids/software/FlexiBLAS/3.5.0-GCC-15.2.0/lib"
)
LIB = str(ROOT / "bbdir-release" / "bindings" / "libseams_lang.so")
FIELDS = [
    "threads", "code", "version", "operation", "molecules", "warmup", "n",
    "median_ms", "q1_ms", "q3_ms", "iqr_ms", "agreement", "detail",
    "samples_ms", "cpus",
]


def lib_env(lib):
    build = Path(lib).resolve().parents[1]
    parts = [
        str(build / "bindings"), str(build / "src"),
        str(build / "subprojects/readcon-core"),
        str(build / "subprojects/readcon-db"),
        str(build / "subprojects/linkcell"),
        str(build / "subprojects/vesin"),
        str(build / "subprojects/minimage"),
        FLEXI,
    ]
    return os.pathsep.join(parts)


def cpus_now():
    try:
        for line in open("/proc/self/status", encoding="utf-8"):
            if line.startswith("Cpus_allowed_list:"):
                return line.split(":", 1)[1].strip()
    except OSError:
        return "unknown"
    return "unknown"


def child_env(kind, threads):
    env = os.environ.copy()
    width = str(threads)
    env["OMP_NUM_THREADS"] = width
    env["OMP_DYNAMIC"] = "FALSE"
    env["MKL_NUM_THREADS"] = "1"
    env["OPENBLAS_NUM_THREADS"] = "1"
    env["NUMEXPR_NUM_THREADS"] = "1"
    env["RAYON_NUM_THREADS"] = width
    env.pop("OMP_PROC_BIND", None)
    if kind == "dseams":
        env["LD_LIBRARY_PATH"] = lib_env(LIB)
        env.pop("OVITO_THREAD_COUNT", None)
    elif kind == "ovito":
        env["LD_LIBRARY_PATH"] = HOST_LIB
        env["OVITO_THREAD_COUNT"] = width
        env["QT_QPA_PLATFORM"] = "offscreen"
    else:
        env.pop("OVITO_THREAD_COUNT", None)
        env.pop("LD_LIBRARY_PATH", None)
    return env


def emit(payload):
    json.dump(payload, sys.stdout)
    sys.stdout.write("\n")


def load_frame(path):
    import numpy as np
    data = np.load(path)
    return data["positions"], float(data["side"])


def dseams_handle(pos, side):
    import ctypes
    import numpy as np
    lib = ctypes.CDLL(LIB)
    lib.seams_last_error.restype = ctypes.c_char_p
    lib.seams_lang_version.restype = ctypes.c_char_p
    lib.seams_frame_from_arrays.restype = ctypes.c_void_p
    lib.seams_frame_from_arrays.argtypes = [
        ctypes.POINTER(ctypes.c_double), ctypes.c_int,
        ctypes.POINTER(ctypes.c_double),
        ctypes.POINTER(ctypes.c_int), ctypes.POINTER(ctypes.c_int)]
    lib.seams_frame_free.argtypes = [ctypes.c_void_p]
    lib.seams_frame_steinhardt.argtypes = [
        ctypes.c_void_p, ctypes.c_double, ctypes.c_int, ctypes.c_int,
        ctypes.POINTER(ctypes.c_double), ctypes.POINTER(ctypes.c_double)]
    lib.seams_frame_neighbours.argtypes = [
        ctypes.c_void_p, ctypes.c_double, ctypes.c_int, ctypes.c_char_p]
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
        raise RuntimeError(lib.seams_last_error().decode())
    return lib, handle, n


def run_worker(args):
    pinned = cpus_now()
    os.environ["OMP_NUM_THREADS"] = str(args.threads)
    os.environ["RAYON_NUM_THREADS"] = str(args.threads)
    sys.path.insert(0, str(HERE.parent))
    import ctypes
    import numpy as np
    import h2h_protocol as proto
    pos, side = load_frame(args.frame)
    n = len(pos)

    if args.worker == "dseams-q":
        lib, handle, n = dseams_handle(pos, side)
        ql = np.empty(n)
        qb = np.empty(n)

        def call(_h):
            status = lib.seams_frame_steinhardt(
                handle, proto.CUTOFF_A, 1, args.order,
                ql.ctypes.data_as(ctypes.POINTER(ctypes.c_double)),
                qb.ctypes.data_as(ctypes.POINTER(ctypes.c_double)))
            if status != 0:
                raise RuntimeError(lib.seams_last_error().decode())
            return ql.copy()

        samples, result = proto.time_calls(lambda: handle, call)
        lib.seams_frame_free(handle)
        emit({
            "code": "d-SEAMS", "version": lib.seams_lang_version().decode(),
            "operation": f"q{args.order}", "molecules": n, "samples_ms": samples,
            "values": result.tolist(), "cpus": pinned, "threads": args.threads,
        })
        return

    if args.worker == "dseams-n":
        lib, handle, n = dseams_handle(pos, side)

        def call(_h):
            status = lib.seams_frame_neighbours(handle, proto.CUTOFF_A, 1, b"cutoff")
            if status != 0:
                raise RuntimeError(lib.seams_last_error().decode())
            return status

        samples, _ = proto.time_calls(lambda: handle, call)
        lib.seams_frame_free(handle)
        emit({
            "code": "d-SEAMS", "version": lib.seams_lang_version().decode(),
            "operation": "neighbours", "molecules": n, "samples_ms": samples,
            "pairs": 4 * n, "cpus": pinned, "threads": args.threads,
        })
        return

    if args.worker == "freud-q":
        import freud
        freud.parallel.set_num_threads(args.threads)
        box = freud.box.Box.cube(side)
        pts = np.ascontiguousarray(pos - side / 2.0)
        query = {"r_max": proto.CUTOFF_A, "exclude_ii": True}

        def call(_p):
            op = freud.order.Steinhardt(args.order)
            op.compute((box, pts), neighbors=query)
            return np.asarray(op.particle_order, dtype=np.float64)

        samples, result = proto.time_calls(lambda: pts, call)
        emit({
            "code": "freud", "version": freud.__version__,
            "operation": f"q{args.order}", "molecules": n, "samples_ms": samples,
            "values": result.tolist(), "cpus": pinned, "threads": args.threads,
        })
        return

    if args.worker == "freud-n":
        import freud
        freud.parallel.set_num_threads(args.threads)
        box = freud.box.Box.cube(side)
        pts = np.ascontiguousarray(pos - side / 2.0)

        def call(_p):
            nlist = freud.locality.AABBQuery(box, pts).query(
                pts, {"r_max": proto.CUTOFF_A, "exclude_ii": True}).toNeighborList()
            return int(len(nlist))

        samples, pairs = proto.time_calls(lambda: pts, call)
        emit({
            "code": "freud", "version": freud.__version__,
            "operation": "neighbours", "molecules": n, "samples_ms": samples,
            "pairs": pairs, "cpus": pinned, "threads": args.threads,
        })
        return

    if args.worker == "pyscal-q":
        import pyscal3
        from ase import Atoms
        pyscal3.set_num_threads(args.threads)
        atoms = Atoms(symbols=["O"] * n, positions=pos, cell=[side, side, side], pbc=True)

        def call(_a):
            pyscal3.find_neighbors(atoms, method="cutoff", cutoff=proto.CUTOFF_A)
            return np.asarray(pyscal3.steinhardt_parameter(atoms, args.order)[0])

        samples, result = proto.time_calls(lambda: atoms, call)
        emit({
            "code": "pyscal3", "version": pyscal3.__version__,
            "operation": f"q{args.order}", "molecules": n, "samples_ms": samples,
            "values": result.tolist(), "cpus": pinned, "threads": args.threads,
        })
        return

    if args.worker == "ovito-q":
        os.environ["OVITO_THREAD_COUNT"] = str(args.threads)
        os.environ["QT_QPA_PLATFORM"] = "offscreen"
        import ovito
        from ovito.data import DataCollection
        from ovito.modifiers import SteinhardtParametersModifier

        def prep():
            data = DataCollection()
            data.create_cell(
                matrix=[[side, 0, 0, 0], [0, side, 0, 0], [0, 0, side, 0]],
                pbc=(True, True, True))
            particles = data.create_particles(count=n)
            particles.create_property("Position", data=np.ascontiguousarray(pos))
            return data

        column = {4: 2, 6: 4}[args.order]

        def call(data):
            data.apply(SteinhardtParametersModifier(
                cutoff=proto.CUTOFF_A,
                neighbor_mode=SteinhardtParametersModifier.NeighborMode.Cutoff,
                max_l=args.order))
            return np.asarray(data.particles["Steinhardt"], dtype=np.float64)[:, column]

        samples, result = proto.time_calls(prep, call)
        emit({
            "code": "OVITO", "version": ovito.version_string,
            "operation": f"q{args.order}", "molecules": n, "samples_ms": samples,
            "values": result.tolist(), "cpus": pinned, "threads": args.threads,
            "note": f"max_l={args.order} fills q2 through q{args.order}",
        })
        return

    if args.worker == "matscipy-n":
        import matscipy
        from ase import Atoms
        from matscipy.neighbours import neighbour_list
        atoms = Atoms(symbols=["O"] * n, positions=pos, cell=[side, side, side], pbc=True)

        def call(_a):
            i, j = neighbour_list("ij", atoms, proto.CUTOFF_A)
            return int(len(i))

        samples, pairs = proto.time_calls(lambda: atoms, call)
        emit({
            "code": "matscipy", "version": matscipy.__version__,
            "operation": "neighbours", "molecules": n, "samples_ms": samples,
            "pairs": pairs, "cpus": pinned, "threads": args.threads,
        })
        return

    if args.worker == "ase-n":
        import ase
        from ase import Atoms
        from ase.neighborlist import neighbor_list
        atoms = Atoms(symbols=["O"] * n, positions=pos, cell=[side, side, side], pbc=True)

        def call(_a):
            i, j = neighbor_list("ij", atoms, proto.CUTOFF_A)
            return int(len(i))

        samples, pairs = proto.time_calls(lambda: atoms, call)
        emit({
            "code": "ASE", "version": ase.__version__,
            "operation": "neighbours", "molecules": n, "samples_ms": samples,
            "pairs": pairs, "cpus": pinned, "threads": args.threads,
        })
        return

    import MDAnalysis
    from MDAnalysis.lib.distances import capped_distance

    def call(_p):
        pairs = capped_distance(
            pos, pos, max_cutoff=proto.CUTOFF_A, min_cutoff=1e-8,
            box=[side, side, side, 90.0, 90.0, 90.0], return_distances=False)
        return int(len(pairs))

    samples, pairs = proto.time_calls(lambda: pos, call)
    emit({
        "code": "MDAnalysis", "version": MDAnalysis.__version__,
        "operation": "neighbours", "molecules": n, "samples_ms": samples,
        "pairs": pairs, "cpus": pinned, "threads": args.threads,
    })


def launch(python, threads, worker, frame, order):
    cores = "0" if threads == 1 else "0-7"
    kind = "ovito" if worker.startswith("ovito") else (
        "dseams" if worker.startswith("dseams") else "other")
    cmd = [ "taskset", "-c", cores, python, str(HERE),
            "--worker", worker, "--threads", str(threads),
            "--frame", frame, "--order", str(order)]
    done = subprocess.run(
        cmd, cwd=ROOT, env=child_env(kind, threads), text=True, capture_output=True)
    if done.returncode != 0:
        tail = (done.stderr or done.stdout or "").strip().splitlines()
        return {"error": tail[-1] if tail else f"exit {done.returncode}",
                "code": worker, "threads": threads}
    line = [row for row in done.stdout.splitlines() if row.startswith("{")][-1]
    payload = json.loads(line)
    if payload.get("cpus") != cores:
        return {"error": f"not pinned: {payload.get('cpus')}", "code": worker,
                "threads": threads}
    return payload


def write_row(writer, payload, operation, fraction, detail):
    sys.path.insert(0, str(HERE.parent))
    import h2h_protocol as proto
    if "error" in payload:
        writer.writerow({
            "threads": payload.get("threads", ""), "code": payload.get("code", ""),
            "version": "", "operation": operation,
            "molecules": payload.get("molecules", ""), "warmup": "", "n": 0,
            "median_ms": "", "q1_ms": "", "q3_ms": "", "iqr_ms": "",
            "agreement": "", "detail": payload["error"], "samples_ms": "", "cpus": "",
        })
        return
    summary = proto.summarize(payload["samples_ms"])
    writer.writerow({
        "threads": payload["threads"], "code": payload["code"],
        "version": payload.get("version", ""), "operation": operation,
        "molecules": payload["molecules"], "warmup": proto.WARMUP,
        "n": summary["n"],
        "median_ms": f"{summary['median_ms']:.3f}",
        "q1_ms": f"{summary['q1_ms']:.3f}",
        "q3_ms": f"{summary['q3_ms']:.3f}",
        "iqr_ms": f"{summary['iqr_ms']:.3f}",
        "agreement": fraction, "detail": detail,
        "samples_ms": " ".join(f"{v:.3f}" for v in payload["samples_ms"]),
        "cpus": payload.get("cpus", ""),
    })


def q_agreement(reference, other):
    import numpy as np
    if "error" in reference or "error" in other:
        return "", other.get("error") or reference.get("error")
    left = np.asarray(reference["values"], dtype=np.float64)
    right = np.asarray(other["values"], dtype=np.float64)
    delta = float(np.max(np.abs(left - right)))
    return f"{delta:.3e}", f"max_abs={delta:.3e} mean_dseams={float(left.mean()):.6f}"


def n_agreement(reference, other, molecules):
    if "error" in reference or "error" in other:
        return "", other.get("error") or reference.get("error")
    expect = 4 * molecules
    got = int(other.get("pairs", -1))
    return f"{float(got == expect):.6f}", f"pairs={got} expect={expect}"


def campaign(args):
    sys.path.insert(0, str(HERE.parent))
    os.environ.pop("LD_LIBRARY_PATH", None)
    import numpy as np
    import h2h_protocol as proto
    python = args.python or sys.executable
    out = ROOT / args.out
    out.parent.mkdir(parents=True, exist_ok=True)
    frames = {}
    for reps in (11, 24):
        pos, side = proto.ice_ic(reps)
        path = ROOT / "benchmarks" / f".ice-ic-{len(pos)}.npz"
        np.savez(path, positions=pos, side=side)
        frames[len(pos)] = str(path)
    with out.open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=FIELDS)
        writer.writeheader()
        big = frames[110592]
        for threads in (1, 8):
            for order in (4, 6):
                reference = launch(python, threads, "dseams-q", big, order)
                reference["threads"] = threads
                reference["molecules"] = 110592
                write_row(writer, reference, f"q{order}", "0.000e+00", "reference ql")
                for worker in ("freud-q", "pyscal-q", "ovito-q"):
                    rival = launch(python, threads, worker, big, order)
                    rival["threads"] = threads
                    rival["molecules"] = 110592
                    fraction, detail = q_agreement(reference, rival)
                    note = rival.get("note")
                    if note and "error" not in rival:
                        detail = f"{detail}; {note}"
                    write_row(writer, rival, f"q{order}", fraction, detail)
                    handle.flush()
                    print(f"threads={threads} q{order} {worker} "
                          f"{rival.get('samples_ms', rival.get('error'))}", flush=True)
            for molecules in (10648, 110592):
                frame = frames[molecules]
                reference = launch(python, threads, "dseams-n", frame, 6)
                reference["threads"] = threads
                reference["molecules"] = molecules
                write_row(writer, reference, "neighbours", "1.000000", "pairs=4N cutoff rows")
                for worker in ("freud-n", "matscipy-n", "ase-n", "mda-n"):
                    rival = launch(python, threads, worker, frame, 6)
                    rival["threads"] = threads
                    rival["molecules"] = molecules
                    fraction, detail = n_agreement(reference, rival, molecules)
                    write_row(writer, rival, "neighbours", fraction, detail)
                    handle.flush()
                    print(f"threads={threads} N={molecules} {worker} "
                          f"{rival.get('samples_ms', rival.get('error'))} {detail}",
                          flush=True)
    for path in frames.values():
        Path(path).unlink(missing_ok=True)


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--out", default="benchmarks/h2h-watch-v7.csv")
    parser.add_argument("--python")
    parser.add_argument("--worker")
    parser.add_argument("--threads", type=int, default=1)
    parser.add_argument("--frame")
    parser.add_argument("--order", type=int, default=6)
    args = parser.parse_args()
    if args.worker:
        run_worker(args)
        return
    campaign(args)


if __name__ == "__main__":
    main()
