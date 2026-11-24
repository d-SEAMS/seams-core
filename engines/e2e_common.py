# Copyright (c) d-SEAMS developers
"""Shared pieces of the engine end-to-end checks.

Each engine writes LAMMPS-style dump frames that carry the coordinates it
handed to d-SEAMS in situ and the label it got back. The offline CLI reads
the same frames and writes its own labels with --per-atom. A check passes
only when every label of every frame is the same.

Exit code 77 is a skip: Meson reports it as SKIP.
"""

import json
import os
import shutil
import subprocess
import sys
from pathlib import Path

SKIP = 77


def skip(why):
    print(f"SKIP: {why}")
    sys.exit(SKIP)


def find_build_file(build, name):
    """The built file, not a copy under an install destdir."""
    hits = sorted(
        p for p in Path(build).rglob(name)
        if p.is_file() and not any("destdir" in part for part in p.parts)
    )
    return hits[0] if hits else None


def _thread_runtime():
    runtime = os.environ.get("TSAN_RUNTIME", "")
    return runtime if runtime and Path(runtime).is_file() else ""


def needs_thread_runtime(build):
    """True when this build's libyodaLib references the thread runtime."""
    lib = Path(build) / "src" / "libyodaLib.so"
    if not lib.is_file():
        found = find_build_file(build, "libyodaLib.so")
        if found is None:
            return False
        lib = found
    try:
        return b"__tsan_" in lib.read_bytes()
    except OSError:
        return False


def drop_thread_runtime(env):
    """Environment for a process that must not map libclang_rt.tsan.

    plumed, gmx and lmp abort if that runtime is loaded at startup, and
    loading it later leaves OpenMP workers unregistered. An executable
    already linked with -fsanitize=thread aborts if the shared runtime
    is preloaded as well.
    """
    env = dict(env)
    runtime = _thread_runtime()
    preload = env.get("LD_PRELOAD", "")
    if not runtime or runtime not in preload.split(":"):
        return env
    kept = [part for part in preload.split(":") if part and part != runtime]
    if kept:
        env["LD_PRELOAD"] = ":".join(kept)
    else:
        env.pop("LD_PRELOAD", None)
    return env


def with_thread_runtime(env):
    """Environment for an unsanitized binary that dlopens the instrumented library.

    gmx-seams is that binary. The runtime has to be mapped before the
    process starts; linking it into gmx-seams makes gmx_fio_open jump
    to address 0, while a preload does not.
    """
    env = drop_thread_runtime(env)
    runtime = _thread_runtime()
    if not runtime:
        return env
    env["LD_PRELOAD"] = runtime
    return env


def external_library(build, name):
    """Shared library an external host can dlopen.

    An AddressSanitizer build returns its own plugin. plumed, gmx and lmp
    cannot map libclang_rt.tsan (benchmarks/logs/tsan-engine-hosts.txt),
    so only a ThreadSanitizer build hands them the unsanitized library.
    """
    own = find_build_file(build, name)
    if own is None or not needs_thread_runtime(build):
        return own
    unsan = os.environ.get("SEAMS_UNSANITIZED_PREFIX", "")
    if not unsan:
        return own
    alt = find_build_file(unsan, name)
    return alt if alt is not None else own


def seams_cli(build):
    cli = os.environ.get("SEAMS_CLI") or Path(build) / "seams"
    if not Path(cli).is_file():
        sys.exit(f"the seams CLI was not built under {build}")
    return str(cli)


def launch_env(argv):
    """Process environment for an engine run.

    Hosted runners advertise Azure MANA devices that cannot open a UCX
    queue. Open MPI then aborts inside MPI_Bcast. The TCP BTL does not
    touch that device. A one-rank run leaves the caller's transport alone.
    """
    env = drop_thread_runtime(dict(os.environ, OMP_NUM_THREADS="1"))
    names = {Path(arg).name for arg in argv[:2]}
    if names & {"mpirun", "mpiexec"}:
        env.setdefault("UCX_TLS", "tcp,self,sm")
        env.setdefault("OMPI_MCA_pml", "ob1")
        env.setdefault("OMPI_MCA_btl", "^openib")
        env.setdefault("OMPI_MCA_rmaps_base_oversubscribe", "1")
    return env


def which(name, env_var=None):
    """Engine executable from env_var, else PATH. None when absent."""
    path = os.environ.get(env_var, "") if env_var else ""
    if path:
        return path if Path(path).is_file() or shutil.which(path) else None
    return shutil.which(name)


def plumed_kernel():
    """PLUMED_KERNEL, else the libplumedKernel.so beside plumed on PATH. None when absent."""
    kernel = os.environ.get("PLUMED_KERNEL", "")
    if kernel:
        return kernel if Path(kernel).is_file() else None
    plumed = shutil.which("plumed")
    if plumed is None:
        return None
    lib = Path(plumed).resolve().parent.parent / "lib" / "libplumedKernel.so"
    return str(lib) if lib.is_file() else None


def read_dump(path, column):
    """Frames of a LAMMPS text dump as (timestep, {id: int(column)})."""
    frames = []
    with open(path) as fh:
        while True:
            line = fh.readline()
            if not line:
                break
            if not line.startswith("ITEM: TIMESTEP"):
                continue
            step = int(fh.readline())
            fh.readline()
            n = int(fh.readline())
            fh.readline()
            for _ in range(3):
                fh.readline()
            names = fh.readline().split()[2:]
            ids = names.index("id")
            col = names.index(column)
            labels = {}
            for _ in range(n):
                parts = fh.readline().split()
                labels[int(parts[ids])] = int(float(parts[col]))
            frames.append((step, labels))
    return frames


def offline_labels(cli, dump, nframes, cutoff, workdir, extra=()):
    """seams chill-plus over every frame of dump, one label per atom."""
    out = Path(workdir) / "offline-labels.dump"
    out.unlink(missing_ok=True)
    cmd = [cli, "chill-plus", "--frame", "1", "--last", str(nframes),
           "--cutoff", str(cutoff), "--jobs", "1", "--per-atom", str(out),
           *extra, str(dump)]
    # The CLI is linked with the thread runtime. A preload of the same
    # runtime, inherited from a Python host, aborts it.
    done = subprocess.run(cmd, capture_output=True, text=True, check=False,
                          env=drop_thread_runtime(os.environ))
    if done.returncode != 0:
        sys.stderr.write(done.stdout + done.stderr)
        sys.exit(f"offline CLI failed: {' '.join(cmd)}")
    return read_dump(out, "label")


def compare(insitu, offline):
    """Frames, atoms, and the mismatching (frame, id, in situ, offline)."""
    if len(insitu) != len(offline):
        sys.exit(f"{len(insitu)} in-situ frames but {len(offline)} offline frames")
    atoms = 0
    bad = []
    for k, ((_, mine), (_, theirs)) in enumerate(zip(insitu, offline)):
        if mine.keys() != theirs.keys():
            sys.exit(f"frame {k + 1}: the atom ids differ")
        atoms += len(mine)
        for atom, label in mine.items():
            if theirs[atom] != label:
                bad.append((k + 1, atom, label, theirs[atom]))
    return len(insitu), atoms, bad


def census(frames):
    counts = {}
    for _, labels in frames:
        for label in labels.values():
            counts[label] = counts.get(label, 0) + 1
    return dict(sorted(counts.items()))


def finish(record_path, record, bad):
    """Print and store the record. Any label mismatch fails the check."""
    record["label_agreement"] = f"{record['atom_frames'] - len(bad)}/{record['atom_frames']}"
    print(json.dumps(record, sort_keys=True))
    if record_path:
        Path(record_path).parent.mkdir(parents=True, exist_ok=True)
        Path(record_path).write_text(json.dumps(record, sort_keys=True) + "\n")
    if bad:
        for frame, atom, mine, theirs in bad[:20]:
            print(f"frame {frame} id {atom}: in situ {mine} offline {theirs}", file=sys.stderr)
        sys.exit(f"{len(bad)} labels differ between the engine and the CLI")
    if record.get("native_dump_mismatches"):
        sys.exit(f"{record['native_dump_mismatches']} labels differ on the engine's own dump")
    if record["frames"] < 2 or record["steps"] < 100:
        sys.exit("the run is shorter than 100 steps or has fewer than two frames")
