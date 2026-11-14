# Copyright (c) d-SEAMS developers
"""mW ice/liquid in LAMMPS with compute dseams, checked against the CLI.

usage: e2e_mw.py BUILD [--ranks N] [--cells C] [--steps S] [--record FILE]

LMP names the lmp executable (default: lmp on PATH). The run loads the
dseamsplugin.so built under BUILD with `plugin load`.
"""

import argparse
import os
import re
import shutil
import subprocess
import sys
import tempfile
from pathlib import Path

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE.parent))
import e2e_common as e2e  # noqa: E402


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("build")
    ap.add_argument("--ranks", type=int, default=1)
    ap.add_argument("--cells", type=int, default=8)
    ap.add_argument("--steps", type=int, default=200)
    ap.add_argument("--every", type=int, default=10)
    ap.add_argument("--melt", type=int, default=2000)
    ap.add_argument("--record", default="")
    args = ap.parse_args()

    lmp = e2e.which("lmp", "LMP")
    if lmp is None:
        e2e.skip("LAMMPS (lmp or $LMP) was not found")
    plugin = e2e.find_build_file(args.build, "dseamsplugin.so")
    if plugin is None:
        e2e.skip("dseamsplugin.so was not built (LAMMPS class headers missing)")
    cli = e2e.seams_cli(args.build)
    launch = [lmp]
    if args.ranks > 1:
        mpirun = shutil.which("mpirun")
        if mpirun is None:
            e2e.skip("mpirun was not found")
        launch = [mpirun, "-np", str(args.ranks), lmp]

    with tempfile.TemporaryDirectory(prefix="dseams-lammps-") as tmp:
        work = Path(tmp)
        for name in ("in.mw.lmp", "mw_setup.lmp", "mW.sw"):
            shutil.copy(HERE / name, work)
        env = dict(os.environ, OMP_NUM_THREADS="1")
        cmd = launch + [
            "-var", "plugin", str(plugin), "-var", "cells", str(args.cells),
            "-var", "steps", str(args.steps), "-var", "every", str(args.every),
            "-var", "melt", str(args.melt), "-var", "dump", "insitu.lammpstrj",
            "-in", "in.mw.lmp", "-log", "log.lammps", "-screen", "none",
        ]
        done = subprocess.run(cmd, cwd=work, env=env, capture_output=True, text=True)
        log = (work / "log.lammps").read_text() if (work / "log.lammps").is_file() else ""
        if done.returncode != 0:
            sys.stderr.write(done.stdout + done.stderr + log[-4000:])
            sys.exit("lmp failed")
        version = re.search(r"LAMMPS \((.+?)\)", log)
        timing = re.search(r"compute dseams calls (\d+) seconds (\S+) per_call (\S+)", log)
        if timing is None:
            sys.exit("compute dseams did not report its timing")
        insitu = e2e.read_dump(work / "insitu.lammpstrj", "c_ice")
        offline = e2e.offline_labels(cli, work / "insitu.lammpstrj", len(insitu), 3.5, work,
                                     extra=("--type", "1"))
        frames, atom_frames, bad = e2e.compare(insitu, offline)
        if int(timing.group(1)) != frames:
            sys.exit(f"compute dseams ran {timing.group(1)} times for {frames} frames")
        record = {
            "engine": "lammps" if args.ranks == 1 else f"lammps-np{args.ranks}",
            "version": version.group(1) if version else "unknown",
            "system": f"mW ice Ic/liquid, {len(insitu[0][1])} molecules",
            "steps": args.steps,
            "frames": frames,
            "atom_frames": atom_frames,
            "labels": e2e.census(insitu),
            "insitu_s_per_frame": float(timing.group(3)),
        }
        e2e.finish(args.record, record, bad)


if __name__ == "__main__":
    main()
