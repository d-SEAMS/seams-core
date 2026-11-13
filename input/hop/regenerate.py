#!/usr/bin/env python3
# Copyright (c) d-SEAMS developers

"""Regenerate production-stride mW trajectories from the published growth window.

The published dump is the 148-frame nucleation trajectory of d-SEAMS 1
(figshare 10.6084/m9.figshare.11448702.v1, md5
c28fc209a2854d2ba49e8f9b55accbe8). Its frames are 1000 steps apart.
At a 10 fs step that is 10 ps, which is coarser than a production
stride of the same integrator.

This script takes one frame from the growth regime of that file, assigns
velocities at 208 K, equilibrates them under the same Nose-Hoover NPT
barostat, and writes a new dump every 10, 100 and 1000 steps. The dumps
are large and stay outside git. The decks are in.equilibrate and in.stride.
mW.sw is the Stillinger-Weber file in input/mw215.

Frame 104 (timestep 15780000) is the half-rise of the published core
cluster, inside the growth window. The velocity seed is that timestep.
Equilibration is 2000 steps (20 ps) and is not written into the
production dumps. Each production run then saves the equilibrated frame
and one frame per stride for 200 transitions.

Usage:
  python3 input/hop/regenerate.py
  python3 input/hop/regenerate.py --lmp lmp --strides 10,100,1000 --transitions 200
"""

from __future__ import annotations

import argparse
import hashlib
import os
import shutil
import subprocess
import urllib.request
from pathlib import Path

ROOT = Path(__file__).resolve().parents[2]
DECK = Path(__file__).resolve().parent
SW = ROOT / "input" / "mw215" / "mW.sw"
OUT = ROOT / "data" / "hop"

FIGSHARE = "https://ndownloader.figshare.com/files/20418309"
FIGSHARE_MD5 = "c28fc209a2854d2ba49e8f9b55accbe8"
# 1-based frame in the published growth window. The half-rise of the
# published core cluster is frame 104 (docs/nucleation-pathways.md).
FRAME = 104
# Step index of that frame. The file runs 15677000..15824000 by 1000.
STEP = 15780000
HEADER = 9
ATOMS = 4096
FRAME_LINES = HEADER + ATOMS


def md5_file(path: Path) -> str:
    digest = hashlib.md5()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1 << 20), b""):
            digest.update(block)
    return digest.hexdigest()


def fetch_published(dest: Path) -> None:
    if dest.is_file() and md5_file(dest) == FIGSHARE_MD5:
        return
    dest.parent.mkdir(parents=True, exist_ok=True)
    partial = dest.with_suffix(dest.suffix + ".partial")
    urllib.request.urlretrieve(FIGSHARE, partial)
    digest = md5_file(partial)
    if digest != FIGSHARE_MD5:
        partial.unlink(missing_ok=True)
        raise SystemExit(f"published trajectory md5 {digest} != {FIGSHARE_MD5}")
    partial.replace(dest)


def extract_frame(published: Path, dest: Path, frame: int) -> int:
    """Write one 1-based frame. Returns the timestep stored in that frame."""
    start = (frame - 1) * FRAME_LINES
    with published.open() as handle:
        for _ in range(start):
            if handle.readline() == "":
                raise SystemExit(f"frame {frame} is past the end of {published.name}")
        chunk = [handle.readline() for _ in range(FRAME_LINES)]
    if any(line == "" for line in chunk):
        raise SystemExit(f"frame {frame} is short")
    if not chunk[0].startswith("ITEM: TIMESTEP"):
        raise SystemExit("frame does not start at ITEM: TIMESTEP; line count drifted")
    step = int(chunk[1])
    nat = int(chunk[3])
    if nat != ATOMS:
        raise SystemExit(f"frame {frame} has {nat} atoms, expected {ATOMS}")
    dest.write_text("".join(chunk))
    return step


def run_lammps(lmp: str, deck: Path, variables: dict[str, str], log_name: str,
               threads: int) -> None:
    cmd = [lmp, "-in", deck.name, "-log", log_name, "-screen", "none"]
    for key, value in variables.items():
        cmd.extend(["-var", key, value])
    env = os.environ.copy()
    env["OMP_NUM_THREADS"] = str(threads)
    # A single rank. The conda-forge build is MPI, and an unset UCX
    # transport hangs on a host without RDMA. OpenMP is the pair style.
    env.setdefault("UCX_TLS", "tcp,self,sm")
    env.setdefault("OMPI_MCA_btl", "^openib")
    env.setdefault("OMPI_MCA_pml", "ob1")
    omp = [lmp, "-sf", "omp", "-pk", "omp", str(threads)]
    cmd = omp + cmd[1:]
    subprocess.run(cmd, cwd=OUT, env=env, check=True)


def count_frames(path: Path) -> int:
    n = 0
    with path.open() as handle:
        for line in handle:
            if line.startswith("ITEM: TIMESTEP"):
                n += 1
    return n


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--lmp", default=os.environ.get("LMP", "lmp"))
    parser.add_argument("--frame", type=int, default=FRAME)
    parser.add_argument("--strides", default="10,100,1000")
    parser.add_argument("--transitions", type=int, default=200)
    parser.add_argument("--nequil", type=int, default=2000)
    parser.add_argument("--seed", type=int, default=STEP)
    parser.add_argument("--threads", type=int, default=os.cpu_count() or 1)
    parser.add_argument("--published", type=Path, default=None)
    args = parser.parse_args()
    if args.transitions < 1:
        raise SystemExit("need at least one transition")
    strides = [int(piece) for piece in args.strides.split(",") if piece.strip()]
    if not strides or any(stride < 1 for stride in strides):
        raise SystemExit(f"bad strides {args.strides}")
    lmp = args.lmp
    if shutil.which(lmp) is None and not Path(lmp).is_file():
        raise SystemExit(f"lmp not found: {lmp}")

    OUT.mkdir(parents=True, exist_ok=True)
    published = args.published or (OUT / "nucleation.lammpstrj")
    fetch_published(published)
    start = OUT / "start.lammpstrj"
    step = extract_frame(published, start, args.frame)
    if args.frame == FRAME and step != STEP:
        raise SystemExit(f"frame {FRAME} timestep is {step}, expected {STEP}")

    shutil.copy(SW, OUT / "mW.sw")
    shutil.copy(DECK / "in.equilibrate", OUT / "in.equilibrate")
    shutil.copy(DECK / "in.stride", OUT / "in.stride")

    run_lammps(
        lmp,
        OUT / "in.equilibrate",
        {"step": str(step), "nequil": str(args.nequil), "seed": str(args.seed)},
        "log.equil",
        args.threads,
    )
    lammps_version = "LAMMPS"
    equil_log = (OUT / "log.equil").read_text(errors="replace").splitlines()
    if equil_log:
        lammps_version = equil_log[0].strip()
    print(
        f"equilibrated frame {args.frame} step {step} "
        f"seed {args.seed} nequil {args.nequil} threads {args.threads}"
    )
    rows = []
    for stride in strides:
        run_lammps(
            lmp,
            OUT / "in.stride",
            {"stride": str(stride), "ntrans": str(args.transitions)},
            f"log.stride-{stride}",
            args.threads,
        )
        dest = OUT / f"stride-{stride}.lammpstrj"
        produced = OUT / "prod.lammpstrj"
        if not produced.is_file():
            raise SystemExit("LAMMPS did not write prod.lammpstrj")
        produced.replace(dest)
        frames = count_frames(dest)
        size = dest.stat().st_size
        print(
            f"stride {stride} frames {frames} transitions {frames - 1} "
            f"bytes {size} path data/hop/stride-{stride}.lammpstrj"
        )
        if frames < args.transitions + 1:
            raise SystemExit(
                f"stride {stride} wrote {frames} frames, need {args.transitions + 1}"
            )
        rows.append((stride, frames, size))
    command = (
        f"python3 input/hop/regenerate.py --frame {args.frame} "
        f"--strides {args.strides} --transitions {args.transitions} "
        f"--nequil {args.nequil} --seed {args.seed} --threads {args.threads}"
    )
    manifest = [
        f"command {command}",
        f"lammps {lammps_version}",
        "source figshare 10.6084/m9.figshare.11448702.v1",
        f"source_md5 {FIGSHARE_MD5}",
        f"frame {args.frame}",
        f"step {step}",
        f"seed {args.seed}",
        f"nequil {args.nequil}",
        "temperature_K 208",
        "ensemble NPT",
        "pressure_atm 1",
        "thermostat Nose-Hoover",
        "barostat Nose-Hoover iso",
        "tdamp_ps 1",
        "pdamp_ps 5",
        "timestep_fs 10",
        "epsilon_kcal_mol 6.189",
        "sigma_A 2.3925",
        "lambda 23.15",
        f"threads {args.threads}",
    ]
    for stride, frames, size in rows:
        manifest.append(
            f"trajectory stride {stride} frames {frames} "
            f"transitions {frames - 1} bytes {size} "
            f"data/hop/stride-{stride}.lammpstrj"
        )
    (OUT / "manifest.txt").write_text("\n".join(manifest) + "\n")
    print(command)


if __name__ == "__main__":
    main()
