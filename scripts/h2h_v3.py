#!/usr/bin/env python3
"""Release-build head-to-head.

d-SEAMS is tests/bench_watch from a meson -Dbuildtype=release tree,
which is -O3 and refuses to compile otherwise. One core, one warm-up,
six timed calls. Frames are written under a temporary directory.

CHILL+ is seams_chill_plus (four mutual neighbours). On ideal ice Ic
that neighbour set is the cutoff shell inside 3.5 A. Rings are the
cutoff primitive rings of at most six members.
"""

from __future__ import annotations

import argparse
import csv
import json
import os
import struct
import subprocess
import sys
import tempfile
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))

import h2h_protocol as proto

proto.apply_pin_env()

import numpy as np

ROOT = Path(__file__).resolve().parents[1]
HERE = Path(__file__).resolve()
CLASS = {
    0: "cubic",
    1: "hexagonal",
    2: "other",
    3: "interfacial",
    4: "hydrate",
    5: "interfacial_hydrate",
    6: "other",
    7: "cubic",
    8: "hexagonal",
}


def launch_env():
    env = proto.pin_env()
    env["OVITO_THREAD_COUNT"] = "1"
    lib = "/usr/lib/x86_64-linux-gnu"
    # Only the OVITO worker needs the host libGL. The release binary
    # must not see that directory, or its dynamic linker picks up the
    # wrong libc. Callers that spawn bench_watch drop this path.
    env["LD_LIBRARY_PATH_OVITO"] = lib
    return env


def run(cmd, env):
    return subprocess.run(cmd, cwd=ROOT, env=env, text=True, capture_output=True)


def require_release(bench):
    env = proto.pin_env()
    env.pop("LD_LIBRARY_PATH", None)
    done = run([str(bench), "--build"], env)
    if done.returncode != 0 or done.stdout.strip() != "release":
        raise SystemExit(
            f"{bench} is not a release -O3 build ({done.stdout!r} {done.stderr!r})")


def write_frame(path, pos, side):
    blob = struct.pack("<Id", len(pos), float(side))
    blob += np.ascontiguousarray(pos, dtype="<f8").tobytes()
    Path(path).write_bytes(blob)


def dseams(bench, operation, frame):
    env = proto.pin_env()
    env.pop("LD_LIBRARY_PATH", None)
    op = "chill" if operation == "chill+" else "rings"
    cmd = ["taskset", "-c", proto.CORE, str(bench), op, frame,
           str(proto.WARMUP), str(proto.REPEATS)]
    done = run(cmd, env)
    if done.returncode != 0:
        tail = (done.stderr or done.stdout or "").strip().splitlines()
        return {"error": tail[-1] if tail else f"exit {done.returncode}",
                "code": "d-SEAMS"}
    payload = json.loads(done.stdout.strip().splitlines()[-1])
    payload["code"] = "d-SEAMS"
    payload["version"] = "release-O3"
    if operation == "chill+":
        payload["labels"] = [CLASS[int(v)] for v in payload["labels"]]
    return payload


def other(python, script, extra, npz):
    env = launch_env()
    env["LD_LIBRARY_PATH"] = env.pop("LD_LIBRARY_PATH_OVITO") + os.pathsep + env.get("LD_LIBRARY_PATH", "")
    cmd = ["taskset", "-c", proto.CORE, python, str(script),
           "--frame", npz, "--warmup", str(proto.WARMUP),
           "--repeats", str(proto.REPEATS), *extra]
    done = run(cmd, env)
    if done.returncode != 0:
        tail = (done.stderr or done.stdout or "").strip().splitlines()
        return {"error": tail[-1] if tail else f"exit {done.returncode}"}
    payload = json.loads(done.stdout.strip().splitlines()[-1])
    if payload.get("cpus") not in (None, proto.CORE):
        return {"error": f"not pinned: {payload.get('cpus')}", "code": payload.get("code")}
    return payload


def agreement(reference, other_payload, operation):
    if "error" in reference or "error" in other_payload:
        return "", other_payload.get("error") or reference.get("error")
    if operation == "chill+":
        left, right = reference["labels"], other_payload["labels"]
        matched = sum(a == b for a, b in zip(left, right))
        return f"{matched / len(left):.6f}", f"matched={matched}/{len(left)}"
    same = list(reference["census"]) == list(other_payload["census"])
    return f"{float(same):.6f}", f"dseams={reference['census']};other={other_payload['census']}"


def write_pair(writer, reference, other_payload, operation, repetitions):
    fraction, detail = agreement(reference, other_payload, operation)
    molecules = int(reference.get("molecules") or other_payload.get("molecules") or 0)
    failed = "error" in reference or "error" in other_payload
    for payload in (reference, other_payload):
        if "error" in payload:
            writer.writerow({
                "code": payload.get("code", "unknown"),
                "version": payload.get("version", ""),
                "operation": operation,
                "molecules": molecules,
                "repetitions": repetitions,
                "warmup": proto.WARMUP,
                "n": 0,
                "median_ms": "",
                "q1_ms": "",
                "q3_ms": "",
                "iqr_ms": "",
                "agreement": "",
                "detail": payload["error"],
                "samples_ms": "",
                "build": payload.get("build", ""),
            })
            continue
        summary = proto.summarize(payload["samples_ms"])
        row_detail = detail
        row_fraction = fraction
        if failed:
            row_fraction = ""
            row_detail = "timed; the other code aborted: " + detail
        writer.writerow({
            "code": payload["code"],
            "version": payload.get("version", ""),
            "operation": operation,
            "molecules": payload.get("molecules", molecules),
            "repetitions": repetitions,
            "warmup": proto.WARMUP,
            "n": summary["n"],
            "median_ms": f"{summary['median_ms']:.6f}",
            "q1_ms": f"{summary['q1_ms']:.6f}",
            "q3_ms": f"{summary['q3_ms']:.6f}",
            "iqr_ms": f"{summary['iqr_ms']:.6f}",
            "agreement": row_fraction,
            "detail": row_detail,
            "samples_ms": ",".join(f"{v:.6f}" for v in payload["samples_ms"]),
            "build": payload.get("build", ""),
        })


def campaign(args):
    bench = Path(args.bench)
    require_release(bench)
    python = args.python or sys.executable
    out = Path(args.out)
    log_path = Path(args.log)
    out.parent.mkdir(parents=True, exist_ok=True)
    log_path.parent.mkdir(parents=True, exist_ok=True)
    fields = [
        "code", "version", "operation", "molecules", "repetitions",
        "warmup", "n", "median_ms", "q1_ms", "q3_ms", "iqr_ms",
        "agreement", "detail", "samples_ms", "build",
    ]
    sizes = [int(v) for v in args.sizes.split(",")] if args.sizes else list(proto.REPETITIONS)
    jobs = [
        ("chill+", "OVITO", ROOT / "scripts/h2h_v2.py", ["--worker", "ovito", "--operation", "chill+"]),
        ("chill+", "mdapy", ROOT / "scripts/h2h_mdapy.py", []),
        ("rings", "OVITO", ROOT / "scripts/h2h_v2.py", ["--worker", "ovito", "--operation", "rings"]),
        ("rings", "matscipy", ROOT / "scripts/h2h_v2.py", ["--worker", "matscipy"]),
        ("rings", "cycless", ROOT / "scripts/h2h_v2.py", ["--worker", "cycless"]),
    ]
    with out.open("w", newline="") as handle, log_path.open("w") as log, tempfile.TemporaryDirectory(prefix="h2h-v3-") as tmp:
        writer = csv.DictWriter(handle, fieldnames=fields)
        writer.writeheader()
        def say(text):
            print(text, flush=True)
            log.write(text + "\n")
            log.flush()
        say(f"bench {bench} build release")
        for repetitions in sizes:
            pos, side = proto.ice_ic(repetitions)
            directory = Path(tmp)
            raw = directory / f"ic-{len(pos)}.bin"
            npz = directory / f"ic-{len(pos)}.npz"
            write_frame(raw, pos, side)
            np.savez(npz, positions=pos, side=np.float64(side))
            say(f"frame molecules={len(pos)}")
            reference = {}
            for operation, name, script, extra in jobs:
                if operation not in reference:
                    reference[operation] = dseams(bench, operation, str(raw))
                payload = other(python, script, extra, str(npz))
                payload.setdefault("code", name)
                write_pair(writer, reference[operation], payload, operation, repetitions)
                handle.flush()
                fraction, detail = agreement(reference[operation], payload, operation)
                if "samples_ms" in payload and "samples_ms" in reference[operation]:
                    left = proto.summarize(reference[operation]["samples_ms"])
                    right = proto.summarize(payload["samples_ms"])
                    say(f"  {operation} {name} d-SEAMS {left['median_ms']:.3f} "
                        f"iqr {left['iqr_ms']:.3f} other {right['median_ms']:.3f} "
                        f"iqr {right['iqr_ms']:.3f} agreement {fraction} {detail}")
                else:
                    say(f"  {operation} {name} FAILED {detail}")
    print(f"wrote {out}", flush=True)


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--bench", default=str(ROOT / "bbdir-release/tests/bench_watch"))
    parser.add_argument("--out", default=str(ROOT / "benchmarks/h2h-watch-v3.csv"))
    parser.add_argument("--log", default=str(ROOT / "benchmarks/h2h-watch-v3.log"))
    parser.add_argument("--python", default="")
    parser.add_argument("--sizes", default="")
    args = parser.parse_args()
    campaign(args)


if __name__ == "__main__":
    main()
