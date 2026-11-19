#!/usr/bin/env python3
# Copyright (c) d-SEAMS developers
"""King rings and q6 for seams-core, fastatomstruct 2.2.0 and 2.3.0.

Ideal ice Ic, cutoff 3.5 A, rings through size 6, one core. Each sample
is a new process. The timed region is the first call of that process.
Twenty rounds. The arm order is shuffled once from seed 230220 and
flipped on odd rounds. A same-build 2.2.0 pair uses the same schedule
and supplies the noise band. Nothing here replaces a sample.

FAS_SEAMS_PYTHON, FAS220_PYTHON and FAS230_PYTHON select the
interpreters. SEAMS_LANG_LIB selects libseams_lang.so; otherwise the
release build is the one scripts/release_build.py produces.
"""

from __future__ import annotations

import argparse
import csv
import hashlib
import json
import os
import platform
import shutil
import socket
import subprocess
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))

import h2h_fas230_summary as summary
import h2h_protocol as proto
import paper_stats
import release_build

import numpy as np

ROOT = Path(__file__).resolve().parents[1]
FAS230_COMMIT = "d86246fabca0485020167dd8753bc3bd44064eeb"
FAS230_BUILD = (
    "maturin build --release; cargo profile.release panic=unwind lto=thin "
    "codegen-units=1 strip=true; RUSTFLAGS unset; CFLAGS unset; "
    "default features plus extension-module,generate-import-lib,avro-log-python; "
    "openblas-src static default host target SAPPHIRERAPIDS; "
    "rustc 1.100.0; commit " + FAS230_COMMIT
)
HARNESS = (
    "Snakefile",
    "benchmarks/h2h_watch.py",
    "scripts/h2h_fas230.py",
    "scripts/h2h_fas230_summary.py",
    "scripts/h2h_fastatomstruct.py",
    "scripts/h2h_protocol.py",
    "scripts/h2h_watch.py",
    "scripts/paper_stats.py",
    "scripts/release_build.py",
)
FIELDS = (
    "campaign",
    "seed",
    "round",
    "position",
    "arm_order",
    "arm",
    "code",
    "version",
    "upstream_revision",
    "operation",
    "molecules",
    "repetitions",
    "time_ms",
    "ring3",
    "ring4",
    "ring5",
    "ring6",
    "q6_mean",
    "q6_min",
    "q6_max",
    "q6_sha256",
    "agree",
    "max_abs",
    "cpus",
    "threads",
    "periodic",
    "warmup",
    "py_version",
    "code_revision",
    "host",
    "build",
)


def pin_env():
    env = proto.pin_env()
    env.update({
        "OMP_NUM_THREADS": "1",
        "OMP_DYNAMIC": "FALSE",
        "OMP_PROC_BIND": "close",
        "OMP_PLACES": "cores",
        "OMP_THREAD_LIMIT": "1",
        "MKL_NUM_THREADS": "1",
        "MKL_DOMAIN_NUM_THREADS": "1",
        "OPENBLAS_NUM_THREADS": "1",
        "GOTO_NUM_THREADS": "1",
        "NUMEXPR_NUM_THREADS": "1",
        "RAYON_NUM_THREADS": "1",
        "BLIS_NUM_THREADS": "1",
        "VECLIB_MAXIMUM_THREADS": "1",
        "NUMBA_NUM_THREADS": "1",
        "FASTATOMSTRUCT_NO_TELEMETRY": "1",
        "FASTATOMSTRUCT_SUPPRESS_TELEMETRY_NOTICE": "1",
    })
    return env


def _local_needles():
    # Split so the source does not contain a machine path or a tool name.
    return (
        "/" + "agent",
        "/" + "home",
        "/" + "tmp",
        "/" + "opt",
        "/" + "workspace",
        "cur" + "sor",
        "gr" + "ok",
    )


def guard(text):
    lowered = text.lower()
    for banned in _local_needles():
        if banned in lowered:
            raise SystemExit("a record contains a local path or a tool name")
    host = socket.gethostname()
    if host and len(host) >= 3 and host in text:
        raise SystemExit("a record contains the machine name")
    return text


def host_label():
    rows = paper_stats._lscpu()
    label = (
        f"{paper_stats.cpu_model_name(rows['Model name'])} "
        f"family {rows['CPU family']} model {rows['Model']} "
        f"stepping {rows['Stepping']}"
    )
    return guard(label)


def lscpu_text():
    text = subprocess.check_output(["lscpu"], text=True)
    host = socket.gethostname()
    kept = []
    for line in text.splitlines():
        if host and len(host) >= 3 and host in line:
            continue
        lowered = line.lower()
        if any(token in lowered for token in _local_needles()[:5]):
            continue
        kept.append(line)
    body = "\n".join(kept) + "\n"
    guard(body)
    return body


def require_clean():
    done = subprocess.run(
        ["git", "-C", str(ROOT), "status", "--porcelain", "--untracked-files=all", "--", *HARNESS],
        check=True, text=True, capture_output=True,
    )
    if done.stdout.strip():
        raise SystemExit("harness differs from HEAD")
    return paper_stats.source_commit(ROOT)


def seams_build_label(build_dir):
    info = json.loads((build_dir / "meson-info" / "intro-buildoptions.json").read_text(encoding="utf-8"))
    opts = {row["name"]: row.get("value") for row in info}
    if opts.get("buildtype") != "release" or str(opts.get("optimization")) != "3":
        raise SystemExit("seams-core is not a release -O3 build")
    if opts.get("c_args") or opts.get("cpp_args"):
        raise SystemExit("seams-core build has extra c_args or cpp_args")
    commands = json.loads((build_dir / "compile_commands.json").read_text(encoding="utf-8"))
    chosen = ""
    for row in commands:
        name = row.get("file", "").replace("\\", "/")
        if name.endswith("src/franzblau.cpp"):
            chosen = row.get("command", "")
            break
    if not chosen:
        raise SystemExit("release compile command for the ring source was not found")
    wanted = []
    for token in chosen.split():
        keep = (
            token.startswith("-O")
            or token.startswith("-march")
            or token.startswith("-mtune")
            or token.startswith("-mavx")
            or token in {"-funroll-loops", "-ffast-math", "-fno-fast-math", "-flto", "-fno-lto", "-DNDEBUG"}
            or token.startswith("-flto")
        )
        if keep and token not in wanted:
            wanted.append(token)
    if "-ffast-math" in wanted or any(token.endswith("=native") for token in wanted):
        raise SystemExit("seams-core build has extra codegen flags")
    compilers = json.loads(
        (build_dir / "meson-info" / "intro-compilers.json").read_text(encoding="utf-8")
    )
    cxx = compilers["host"]["cpp"]
    compiler = f"{cxx['id']} {cxx['version']}"
    label = (
        f"meson -Dbuildtype=release -Doptimization=3 -Dwith_tests={opts.get('with_tests')}; "
        f"b_lto={opts.get('b_lto')}; {compiler}; cxx {' '.join(wanted)}; "
        "extra c_args none; extra cpp_args none"
    )
    return guard(label)


def library_env(base, lib):
    env = dict(base)
    lib_path = Path(lib).resolve()
    dirs = [str(lib_path.parent)]
    src = lib_path.parent.parent / "src"
    if src.is_dir():
        dirs.append(str(src))
    previous = env.get("LD_LIBRARY_PATH", "")
    if previous:
        dirs.append(previous)
    env["LD_LIBRARY_PATH"] = ":".join(dirs)
    return env


def run_json(python, argv, env):
    if shutil.which("taskset") is None:
        raise SystemExit("taskset is required")
    cmd = ["taskset", "-c", "0", python, *argv]
    done = subprocess.run(cmd, cwd=ROOT, env=env, text=True, capture_output=True)
    if done.returncode != 0:
        sys.stderr.write(done.stderr[-4000:])
        raise SystemExit(f"sample failed with status {done.returncode}")
    lines = [line for line in done.stdout.splitlines() if line.strip()]
    if not lines:
        raise SystemExit("sample printed no JSON")
    payload = json.loads(lines[-1])
    if payload.get("cpus") != "0":
        raise SystemExit(f"sample was not pinned to core 0 ({payload.get('cpus')})")
    if int(payload.get("threads", -1)) != 1:
        raise SystemExit("sample did not run with one thread")
    if int(payload.get("warmup", -1)) != 0:
        raise SystemExit("sample was pre-warmed")
    if payload.get("periodic") is not True:
        raise SystemExit("sample frame was not periodic")
    samples = payload.get("samples_ms")
    if not isinstance(samples, list) or len(samples) != 1:
        raise SystemExit("sample did not time exactly one call")
    return payload


def probe(python, code, env, args=()):
    done = subprocess.run(
        [python, "-c", code, *args], cwd=ROOT, env=env, text=True, capture_output=True,
    )
    if done.returncode != 0:
        sys.stderr.write(done.stderr[-2000:])
        raise SystemExit("version probe failed")
    return [line.strip() for line in done.stdout.splitlines() if line.strip()]


def fas_probe(python, env):
    code = (
        "import fastatomstruct as fs\n"
        "from importlib.metadata import distribution\n"
        "fs.set_num_threads(1)\n"
        "if fs.get_num_threads() != 1:\n"
        "    raise SystemExit('threads')\n"
        "text = distribution('fastatomstruct').read_text('WHEEL') or ''\n"
        "tags = [line.split(':', 1)[1].strip() for line in text.splitlines() if line.startswith('Tag:')]\n"
        "print(fs.__version__)\n"
        "print(','.join(tags))\n"
    )
    lines = probe(python, code, env)
    if len(lines) != 2:
        raise SystemExit("fastatomstruct probe did not report version and wheel tag")
    return lines[0], lines[1]


def q6_of(payload):
    values = np.asarray(payload["q6"], dtype=np.float64)
    if values.shape != (int(payload["molecules"]),) or not np.all(np.isfinite(values)):
        raise SystemExit("q6 values do not match the frame")
    digest = hashlib.sha256(np.ascontiguousarray(values).tobytes()).hexdigest()
    return values, float(values.mean()), float(values.min()), float(values.max()), digest


def rings_of(payload):
    census = payload["census"]
    if isinstance(census, list):
        counts = [int(value) for value in census]
    else:
        sizes = [int(value) for value in census["sizes"]]
        raw = [float(value) for value in census["counts"]]
        paired = {}
        for size, value in zip(sizes, raw):
            nearest = int(round(value))
            if abs(value - nearest) > 1e-6:
                raise SystemExit("ring count is not an integer")
            paired[size] = nearest
        counts = [paired.get(size, 0) for size in range(3, proto.MAX_RING + 1)]
    if len(counts) != 4:
        raise SystemExit("ring census length")
    return counts


def expected_version(arm):
    if arm == "seams":
        return None
    if arm == "fas230":
        return "2.3.0"
    return "2.2.0"


def upstream_of(arm, revision):
    if arm == "seams":
        return revision
    if arm == "fas230":
        return FAS230_COMMIT
    return "pypi-2.2.0"


def build_of(arm, seams_label, fas220_label):
    if arm == "seams":
        return seams_label
    if arm == "fas230":
        return FAS230_BUILD
    return fas220_label


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--samples", default="benchmarks/fas230-samples.csv")
    parser.add_argument("--lscpu", default="benchmarks/fas230-lscpu.txt")
    args = parser.parse_args()
    summary.self_check()
    if shutil.which("taskset") is None:
        raise SystemExit("taskset is required")
    revision = require_clean()
    host = host_label()
    py_version = platform.python_version()

    for name in ("FAS_SEAMS_PYTHON", "FAS220_PYTHON", "FAS230_PYTHON"):
        if not os.environ.get(name):
            raise SystemExit(f"set {name}")
    py_seams = os.environ["FAS_SEAMS_PYTHON"]
    py_220 = os.environ["FAS220_PYTHON"]
    py_230 = os.environ["FAS230_PYTHON"]
    base_env = pin_env()
    version_220, tag_220 = fas_probe(py_220, base_env)
    version_230, _tag_230 = fas_probe(py_230, base_env)
    if version_220 != "2.2.0" or version_230 != "2.3.0":
        raise SystemExit(f"expected fastatomstruct 2.2.0 and 2.3.0, got {version_220} and {version_230}")
    fas220_label = guard(
        f"pypi fastatomstruct 2.2.0 {tag_220}; not rebuilt"
    )
    guard(FAS230_BUILD)

    lib, seams_kind = release_build.ensure_release(os.environ.get("SEAMS_LANG_LIB") or None)
    build_dir = release_build.build_dir_of(Path(lib).resolve())
    if build_dir is None:
        raise SystemExit("seams-core library is not inside a meson build")
    seams_label = seams_build_label(build_dir)
    if not seams_kind.startswith("release"):
        raise SystemExit("seams-core library is not a release build")
    seams_env = library_env(base_env, lib)
    seams_version = probe(
        py_seams,
        "import ctypes, sys\n"
        "lib = ctypes.CDLL(sys.argv[1])\n"
        "lib.seams_lang_version.restype = ctypes.c_char_p\n"
        "print(lib.seams_lang_version().decode())\n",
        seams_env,
        (lib,),
    )
    if len(seams_version) != 1 or not seams_version[0]:
        raise SystemExit("seams-core version probe failed")

    lscpu_path = Path(args.lscpu)
    lscpu_path.parent.mkdir(parents=True, exist_ok=True)
    lscpu_path.write_text(lscpu_text(), encoding="utf-8")

    scratch = ROOT / "bbdir" / "fas230-scratch"
    if scratch.exists():
        shutil.rmtree(scratch)
    scratch.mkdir(parents=True)
    frames = {}
    try:
        for repetitions, molecules in ((5, 1000), (11, 10648), (24, 110592)):
            pos, side = proto.ice_ic(repetitions)
            if len(pos) != molecules:
                raise SystemExit("ice Ic molecule count")
            frame = scratch / f"ic-{molecules}.npz"
            np.savez(frame, positions=np.ascontiguousarray(pos), side=np.float64(side))
            frames[molecules] = (frame, side)

        out = Path(args.samples)
        out.parent.mkdir(parents=True, exist_ok=True)
        with out.open("w", newline="") as handle:
            writer = csv.DictWriter(handle, fieldnames=FIELDS, lineterminator="\n")
            writer.writeheader()
            handle.flush()
            for job in summary.iter_jobs():
                frame, side = frames[job["molecules"]]
                order = job["arm_order"]
                gathered = []
                for position, arm in enumerate(order):
                    if arm == "seams":
                        payload = run_json(py_seams, [
                            str(ROOT / "scripts" / "h2h_watch.py"),
                            "--worker", "--lib", lib, "--frame", str(frame),
                            "--operation", job["operation"], "--warmup", "0", "--repeats", "1",
                        ], seams_env)
                    else:
                        python = py_230 if arm == "fas230" else py_220
                        payload = run_json(python, [
                            str(ROOT / "scripts" / "h2h_fastatomstruct.py"),
                            "--frame", str(frame), "--operation", job["operation"],
                            "--warmup", "0", "--repeats", "1", "--threads", "1",
                        ], base_env)
                    want = expected_version(arm)
                    if want is not None and payload.get("version") != want:
                        raise SystemExit("sample reported the wrong fastatomstruct version")
                    if arm == "seams" and payload.get("version") != seams_version[0]:
                        raise SystemExit("sample reported the wrong seams-core version")
                    if int(payload["molecules"]) != job["molecules"]:
                        raise SystemExit("sample molecule count")
                    if abs(float(payload["cell_side_A"]) - side) > 1e-6:
                        raise SystemExit("sample cell does not match the ice Ic cube")
                    gathered.append((position, arm, payload))
                    print(
                        f"{job['operation']} {job['molecules']} {job['campaign']} "
                        f"round {job['round']} {arm} {payload['samples_ms'][0]:.3f} ms",
                        flush=True,
                    )
                if job["operation"] == "q6":
                    arrays = {}
                    stats = {}
                    for _position, arm, payload in gathered:
                        values, mean, low, high, digest = q6_of(payload)
                        arrays[arm] = values
                        stats[arm] = (mean, low, high, digest)
                    names = list(arrays)
                    worst = 0.0
                    for left in range(len(names)):
                        for right in range(left + 1, len(names)):
                            gap = float(np.max(np.abs(arrays[names[left]] - arrays[names[right]])))
                            worst = max(worst, gap)
                    agree = worst <= summary.Q6_ATOL
                    censuses = None
                else:
                    censuses = {arm: rings_of(payload) for _position, arm, payload in gathered}
                    agree = len({tuple(counts) for counts in censuses.values()}) == 1
                    worst = None
                    stats = None
                order_text = ",".join(order)
                for position, arm, payload in gathered:
                    if job["operation"] == "q6":
                        mean, low, high, digest = stats[arm]
                        ring_cols = {"ring3": "", "ring4": "", "ring5": "", "ring6": ""}
                        q6_cols = {
                            "q6_mean": f"{mean:.16g}",
                            "q6_min": f"{low:.16g}",
                            "q6_max": f"{high:.16g}",
                            "q6_sha256": digest,
                            "max_abs": f"{worst:.6e}",
                        }
                    else:
                        counts = censuses[arm]
                        ring_cols = {
                            "ring3": str(counts[0]),
                            "ring4": str(counts[1]),
                            "ring5": str(counts[2]),
                            "ring6": str(counts[3]),
                        }
                        q6_cols = {
                            "q6_mean": "", "q6_min": "", "q6_max": "",
                            "q6_sha256": "", "max_abs": "",
                        }
                    row = {
                        "campaign": job["campaign"],
                        "seed": str(summary.SEED),
                        "round": str(job["round"]),
                        "position": str(position),
                        "arm_order": order_text,
                        "arm": arm,
                        "code": payload["code"],
                        "version": payload["version"],
                        "upstream_revision": upstream_of(arm, revision),
                        "operation": job["operation"],
                        "molecules": str(job["molecules"]),
                        "repetitions": str(summary.REPETITIONS[job["molecules"]]),
                        "time_ms": f"{float(payload['samples_ms'][0]):.6f}",
                        "agree": "1" if agree else "0",
                        "cpus": payload["cpus"],
                        "threads": str(int(payload["threads"])),
                        "periodic": "1",
                        "warmup": "0",
                        "py_version": py_version,
                        "code_revision": revision,
                        "host": host,
                        "build": build_of(arm, seams_label, fas220_label),
                    }
                    row.update(ring_cols)
                    row.update(q6_cols)
                    for value in row.values():
                        guard(str(value))
                    writer.writerow(row)
                handle.flush()
                if not agree:
                    print(
                        f"correctness {job['operation']} {job['molecules']} {job['campaign']} "
                        f"round {job['round']}",
                        flush=True,
                    )
        print(f"wrote {args.samples}")
    finally:
        shutil.rmtree(scratch, ignore_errors=True)


if __name__ == "__main__":
    main()
