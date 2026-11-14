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
    hits = sorted(p for p in Path(build).rglob(name) if p.is_file())
    return hits[0] if hits else None


def seams_cli(build):
    cli = os.environ.get("SEAMS_CLI") or Path(build) / "seams"
    if not Path(cli).is_file():
        sys.exit(f"the seams CLI was not built under {build}")
    return str(cli)


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
    done = subprocess.run(cmd, capture_output=True, text=True, check=False)
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
