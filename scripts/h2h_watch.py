#!/usr/bin/env python3
"""Head-to-head of d-SEAMS against mdapy, fastatomstruct and vesin.

Protocol matches the comparator campaign: taskset on the allowed cores,
one warm-up, six timed calls, median and type-7 interquartile range.
Frames are ideal ice Ic. Scratch frames stay under bbdir and are removed
before exit.

    H2H_PYTHON=/path/to/venv/bin/python scripts/h2h_watch.py \
        --lib bbdir-release/bindings/libseams_lang.so \
        --out benchmarks/h2h-watch.csv

The venv needs mdapy==2.0.0, fastatomstruct and vesin. d-SEAMS is libseams_lang.
"""

from __future__ import annotations

import argparse
import csv
import ctypes
import json
import os
import shutil
import subprocess
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))

import h2h_protocol as proto
import release_build

proto.apply_pin_env()

import numpy as np

ROOT = Path(__file__).resolve().parents[1]
Q6_ATOL = 1e-6


class Seams:
    def __init__(self, path):
        lib = ctypes.CDLL(path)
        lib.seams_lang_version.restype = ctypes.c_char_p
        lib.seams_last_error.restype = ctypes.c_char_p
        lib.seams_frame_from_arrays.restype = ctypes.c_void_p
        lib.seams_frame_from_arrays.argtypes = [
            ctypes.POINTER(ctypes.c_double), ctypes.c_int,
            ctypes.POINTER(ctypes.c_double),
            ctypes.POINTER(ctypes.c_int), ctypes.POINTER(ctypes.c_int)]
        lib.seams_frame_free.argtypes = [ctypes.c_void_p]
        lib.seams_frame_chill_plus.argtypes = [
            ctypes.c_void_p, ctypes.c_double, ctypes.c_int, ctypes.c_char_p,
            ctypes.POINTER(ctypes.c_int), ctypes.POINTER(ctypes.c_int)]
        lib.seams_frame_rings.argtypes = [
            ctypes.c_void_p, ctypes.c_double, ctypes.c_int, ctypes.c_char_p,
            ctypes.c_int, ctypes.POINTER(ctypes.c_int), ctypes.c_int]
        lib.seams_frame_steinhardt.argtypes = [
            ctypes.c_void_p, ctypes.c_double, ctypes.c_int, ctypes.c_int,
            ctypes.POINTER(ctypes.c_double), ctypes.POINTER(ctypes.c_double)]
        lib.seams_frame_neighbours.argtypes = [
            ctypes.c_void_p, ctypes.c_double, ctypes.c_int, ctypes.c_char_p]
        lib.seams_frame_neighbour_pairs.argtypes = [ctypes.c_void_p]
        lib.seams_frame_neighbour_pairs.restype = ctypes.c_int
        self.lib = lib
        self.version = lib.seams_lang_version().decode()

    def frame(self, pos, side):
        n = len(pos)
        xyz = np.ascontiguousarray(pos, dtype=np.float64)
        box = np.array([0, 0, 0, side, side, side, 0, 0, 0], dtype=np.float64)
        types = np.ones(n, dtype=np.int32)
        ids = np.arange(1, n + 1, dtype=np.int32)
        handle = self.lib.seams_frame_from_arrays(
            xyz.ctypes.data_as(ctypes.POINTER(ctypes.c_double)), n,
            box.ctypes.data_as(ctypes.POINTER(ctypes.c_double)),
            types.ctypes.data_as(ctypes.POINTER(ctypes.c_int)),
            ids.ctypes.data_as(ctypes.POINTER(ctypes.c_int)))
        if not handle:
            raise RuntimeError(self.lib.seams_last_error().decode())
        # The library copies the arrays. Keep them alive until the call returns.
        del xyz, box, types, ids
        return handle

    def chill(self, handle, n):
        labels = np.empty(n, dtype=np.int32)
        counts = (ctypes.c_int * 9)()
        status = self.lib.seams_frame_chill_plus(
            handle, proto.CUTOFF_A, 1, b"cutoff",
            labels.ctypes.data_as(ctypes.POINTER(ctypes.c_int)), counts)
        if status != 0:
            raise RuntimeError(self.lib.seams_last_error().decode())
        return labels

    def rings(self, handle):
        counts = (ctypes.c_int * (proto.MAX_RING + 1))()
        status = self.lib.seams_frame_rings(
            handle, proto.CUTOFF_A, 1, b"cutoff", proto.MAX_RING, counts,
            proto.MAX_RING + 1)
        if status != 0:
            raise RuntimeError(self.lib.seams_last_error().decode())
        return [int(counts[i]) for i in range(3, proto.MAX_RING + 1)]

    def neighbours(self, handle):
        status = self.lib.seams_frame_neighbours(
            handle, proto.CUTOFF_A, 1, b"cutoff")
        if status != 0:
            raise RuntimeError(self.lib.seams_last_error().decode())
        return int(self.lib.seams_frame_neighbour_pairs(handle))

    def q6(self, handle, n):
        ql = np.empty(n, dtype=np.float64)
        qlbar = np.empty(n, dtype=np.float64)
        status = self.lib.seams_frame_steinhardt(
            handle, proto.CUTOFF_A, 1, 6,
            ql.ctypes.data_as(ctypes.POINTER(ctypes.c_double)),
            qlbar.ctypes.data_as(ctypes.POINTER(ctypes.c_double)))
        if status != 0:
            raise RuntimeError(self.lib.seams_last_error().decode())
        return ql


def worker(args):
    data = np.load(args.frame)
    pos = data["positions"]
    side = float(data["side"])
    seams = Seams(args.lib)
    handle = seams.frame(pos, side)
    n = len(pos)

    def prep():
        return handle

    if args.operation == "chill+":
        def call(frame):
            labels = seams.chill(frame, n)
            return [proto.DSEAMS_CLASS[int(v)] for v in labels]
    elif args.operation == "rings":
        def call(frame):
            return seams.rings(frame)
    elif args.operation == "q6":
        def call(frame):
            return seams.q6(frame, n)
    elif args.operation == "neighbours":
        def call(frame):
            return seams.neighbours(frame)
    else:
        raise SystemExit(args.operation)
    if args.once:
        call(prep())
        seams.lib.seams_frame_free(handle)
        return
    samples, result = proto.time_calls(prep, call, args.warmup, args.repeats)
    payload = {
        "code": "d-SEAMS",
        "version": seams.version,
        "operation": args.operation,
        "molecules": int(n),
        "samples_ms": samples,
        "cpus": proto.cpus_allowed(),
    }
    if args.operation == "chill+":
        payload["labels"] = result
    elif args.operation == "rings":
        payload["census"] = result
    elif args.operation == "neighbours":
        payload["pairs"] = int(result)
    else:
        payload["q6"] = result.tolist()
    payload["threads"] = int(os.environ.get("OMP_NUM_THREADS", "1"))
    seams.lib.seams_frame_free(handle)
    json.dump(payload, sys.stdout)
    sys.stdout.write("\n")


def core_list(threads):
    if threads <= 1:
        return "0"
    return "0-" + str(threads - 1)


def launch(python, script, extra, frame, env_extra=None, threads=1, timeout=None):
    cores = core_list(threads)
    cmd = ["taskset", "-c", cores, python, str(script), "--frame", str(frame), *extra]
    env = proto.pin_env()
    env["OMP_NUM_THREADS"] = str(threads)
    env["RAYON_NUM_THREADS"] = str(threads)
    if env_extra:
        env.update(env_extra)
    try:
        done = subprocess.run(
            cmd, cwd=ROOT, env=env, text=True, capture_output=True, timeout=timeout)
    except subprocess.TimeoutExpired as exc:
        raise TimeoutError(f"{' '.join(cmd)} exceeded {timeout}s") from exc
    if done.returncode != 0:
        raise RuntimeError(
            f"{' '.join(cmd)}\n{done.stdout}\n{done.stderr}")
    if done.stderr.strip():
        print(done.stderr, file=sys.stderr, end="" if done.stderr.endswith("\n") else "\n")
    line = done.stdout.strip().splitlines()[-1]
    payload = json.loads(line)
    if payload.get("cpus") != cores:
        raise RuntimeError(f"worker not pinned to {cores}: {payload.get('cpus')}")
    return payload


def agreement(reference, other, operation):
    if operation == "chill+":
        left = reference["labels"]
        right = other["labels"]
        if len(left) != len(right):
            raise RuntimeError("label length mismatch")
        matched = sum(a == b for a, b in zip(left, right))
        return matched / len(left), f"matched={matched}/{len(left)}"
    if operation == "rings":
        left = reference["census"]
        right = other["census"]["counts"]
        same = left == right
        return float(same), f"dseams={left};other={right}"
    if operation == "neighbours":
        left = int(reference["pairs"])
        right = int(other["pairs"])
        return float(left == right), f"dseams={left};other={right}"
    left = np.asarray(reference["q6"], dtype=np.float64)
    right = np.asarray(other["q6"], dtype=np.float64)
    delta = np.abs(left - right)
    matched = int(np.count_nonzero(delta <= Q6_ATOL))
    return matched / len(left), f"matched={matched}/{len(left)};max_abs={float(delta.max()):.3e}"


def write_row(handle, reference, other, operation, repetitions, build):
    summary_ref = proto.summarize(reference["samples_ms"])
    summary_other = proto.summarize(other["samples_ms"])
    fraction, detail = agreement(reference, other, operation)
    molecules = molecules_of(reference)
    for summary, payload in ((summary_ref, reference), (summary_other, other)):
        handle.writerow({
            "code": payload["code"],
            "version": payload["version"],
            "build": build,
            "operation": operation,
            "molecules": molecules,
            "repetitions": repetitions,
            "warmup": proto.WARMUP,
            "n": summary["n"],
            "median_ms": f"{summary['median_ms']:.6f}",
            "q1_ms": f"{summary['q1_ms']:.6f}",
            "q3_ms": f"{summary['q3_ms']:.6f}",
            "iqr_ms": f"{summary['iqr_ms']:.6f}",
            "agreement": f"{fraction:.6f}",
            "detail": detail,
            "samples_ms": ",".join(f"{v:.6f}" for v in payload["samples_ms"]),
        })


def molecules_of(payload):
    return int(payload["molecules"])


def bind_library(args):
    path, label = release_build.ensure_release(args.lib or None)
    args.lib = path
    args.build = label
    return path, label


def campaign(args):
    bind_library(args)
    python = args.python or sys.executable
    out = Path(args.out)
    out.parent.mkdir(parents=True, exist_ok=True)
    fields = [
        "code", "version", "build", "operation", "molecules", "repetitions",
        "warmup", "n", "median_ms", "q1_ms", "q3_ms", "iqr_ms",
        "agreement", "detail", "samples_ms",
    ]
    sizes = [int(v) for v in args.sizes.split(",")] if args.sizes else list(proto.REPETITIONS)
    scratch = ROOT / "bbdir" / "h2h-scratch"
    if scratch.exists():
        shutil.rmtree(scratch)
    scratch.mkdir(parents=True)
    try:
        _campaign_body(args, python, out, fields, sizes, scratch)
    finally:
        shutil.rmtree(scratch, ignore_errors=True)


def _campaign_body(args, python, out, fields, sizes, scratch):
    with out.open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields)
        writer.writeheader()
        for repetitions in sizes:
            pos, side = proto.ice_ic(repetitions)
            frame = scratch / f"ic-{len(pos)}.npz"
            np.savez(frame, positions=pos, side=np.float64(side))
            print(f"frame molecules={len(pos)} side_A={side:.6f}", flush=True)
            for operation, script, extra in (
                ("chill+", ROOT / "scripts/h2h_mdapy.py", []),
                ("q6", ROOT / "scripts/h2h_fastatomstruct.py", ["--operation", "q6"]),
                ("rings", ROOT / "scripts/h2h_fastatomstruct.py", ["--operation", "rings"]),
            ):
                reference = launch(python, Path(__file__), [
                    "--worker", "--lib", args.lib, "--operation", operation,
                    "--warmup", str(proto.WARMUP), "--repeats", str(proto.REPEATS),
                ], frame)
                other = launch(python, script, [
                    *extra, "--warmup", str(proto.WARMUP), "--repeats", str(proto.REPEATS),
                ], frame)
                write_row(writer, reference, other, operation, repetitions, args.build)
                if operation == "chill+":
                    ovito = launch(python, ROOT / "scripts/h2h_ovito.py", [
                        "--warmup", str(proto.WARMUP), "--repeats", str(proto.REPEATS),
                    ], frame, {
                        "OVITO_THREAD_COUNT": "1",
                        "QT_QPA_PLATFORM": "offscreen",
                        "LD_LIBRARY_PATH": "/usr/lib/x86_64-linux-gnu",
                    })
                    fraction_o, detail_o = agreement(reference, ovito, "chill+")
                    summary_o = proto.summarize(ovito["samples_ms"])
                    writer.writerow({
                        "code": ovito["code"],
                        "version": ovito["version"],
                        "build": args.build,
                        "operation": "chill+",
                        "molecules": molecules_of(reference),
                        "repetitions": repetitions,
                        "warmup": proto.WARMUP,
                        "n": summary_o["n"],
                        "median_ms": f"{summary_o['median_ms']:.6f}",
                        "q1_ms": f"{summary_o['q1_ms']:.6f}",
                        "q3_ms": f"{summary_o['q3_ms']:.6f}",
                        "iqr_ms": f"{summary_o['iqr_ms']:.6f}",
                        "agreement": f"{fraction_o:.6f}",
                        "detail": detail_o,
                        "samples_ms": ",".join(f"{v:.6f}" for v in ovito["samples_ms"]),
                    })
                    print(
                        f"  chill+ OVITO {summary_o['median_ms']:.3f} ms "
                        f"(iqr {summary_o['iqr_ms']:.3f}) agreement {fraction_o:.4f} {detail_o}",
                        flush=True,
                    )
                handle.flush()
                fraction, detail = agreement(reference, other, operation)
                ref_s = proto.summarize(reference["samples_ms"])
                oth_s = proto.summarize(other["samples_ms"])
                print(
                    f"  {operation} d-SEAMS {ref_s['median_ms']:.3f} ms "
                    f"({reference['code']} iqr {ref_s['iqr_ms']:.3f}) "
                    f"{other['code']} {oth_s['median_ms']:.3f} ms "
                    f"(iqr {oth_s['iqr_ms']:.3f}) agreement {fraction:.4f} {detail}",
                    flush=True,
                )
    print(f"wrote {out}", flush=True)


def once(args):
    """One untimed-setup call for the Airspeed Velocity arm. No CSV."""
    if not getattr(args, "arm", ""):
        args.arm = args.once
    bind_library(args)
    repetitions = {proto.molecules_for(r): r for r in proto.REPETITIONS}[args.molecules]
    pos, side = proto.ice_ic(repetitions)
    scratch = ROOT / "bbdir" / "h2h-once"
    if scratch.exists():
        shutil.rmtree(scratch)
    scratch.mkdir(parents=True)
    try:
        frame = scratch / "frame.npz"
        np.savez(frame, positions=pos, side=np.float64(side))
        if args.arm == "dseams-chill":
            worker(argparse.Namespace(
                frame=frame, lib=args.lib, operation="chill+", once=True,
                warmup=0, repeats=1))
        elif args.arm == "mdapy-chill":
            import importlib.util
            spec = importlib.util.spec_from_file_location(
                "h2h_mdapy", ROOT / "scripts/h2h_mdapy.py")
            module = importlib.util.module_from_spec(spec)
            spec.loader.exec_module(module)
            module.run(frame, 0, 1, True)
        elif args.arm == "mdapy-rings":
            import importlib.util
            spec = importlib.util.spec_from_file_location(
                "h2h_mdapy_rings", ROOT / "scripts/h2h_mdapy_rings.py")
            module = importlib.util.module_from_spec(spec)
            spec.loader.exec_module(module)
            module.run(frame, 0, 1, True)
        elif args.arm == "dseams-neighbours":
            worker(argparse.Namespace(
                frame=frame, lib=args.lib, operation="neighbours", once=True,
                warmup=0, repeats=1))
        elif args.arm == "vesin-neighbours":
            import importlib.util
            spec = importlib.util.spec_from_file_location(
                "h2h_vesin", ROOT / "scripts/h2h_vesin.py")
            module = importlib.util.module_from_spec(spec)
            spec.loader.exec_module(module)
            module.run(frame, 0, 1, True)
        elif args.arm == "dseams-rings":
            worker(argparse.Namespace(
                frame=frame, lib=args.lib, operation="rings", once=True,
                warmup=0, repeats=1))
        elif args.arm == "fastatomstruct-rings":
            import importlib.util
            spec = importlib.util.spec_from_file_location(
                "h2h_fastatomstruct", ROOT / "scripts/h2h_fastatomstruct.py")
            module = importlib.util.module_from_spec(spec)
            spec.loader.exec_module(module)
            module.run(frame, "rings", 0, 1, True)
        elif args.arm == "dseams-q6":
            worker(argparse.Namespace(
                frame=frame, lib=args.lib, operation="q6", once=True,
                warmup=0, repeats=1))
        elif args.arm == "fastatomstruct-q6":
            import importlib.util
            spec = importlib.util.spec_from_file_location(
                "h2h_fastatomstruct", ROOT / "scripts/h2h_fastatomstruct.py")
            module = importlib.util.module_from_spec(spec)
            spec.loader.exec_module(module)
            module.run(frame, "q6", 0, 1, True)
        else:
            raise SystemExit(args.arm)
    finally:
        shutil.rmtree(scratch, ignore_errors=True)


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--lib", default="", help="meson release -O3 libseams_lang.so")
    parser.add_argument("--out", default=str(ROOT / "benchmarks/h2h-watch.csv"))
    parser.add_argument("--python", default="")
    parser.add_argument("--sizes", default="", help="comma-separated cell repetitions")
    parser.add_argument("--worker", action="store_true")
    parser.add_argument("--frame", default="")
    parser.add_argument("--operation", default="chill+")
    parser.add_argument("--warmup", type=int, default=proto.WARMUP)
    parser.add_argument("--repeats", type=int, default=proto.REPEATS)
    parser.add_argument("--once", default="", help="asv arm name")
    parser.add_argument("--molecules", type=int, default=0)
    args = parser.parse_args()
    if args.worker:
        bind_library(args)
        worker(args)
    elif args.once:
        once(args)
    else:
        campaign(args)


if __name__ == "__main__":
    main()
