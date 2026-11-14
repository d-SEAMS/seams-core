# Copyright (c) d-SEAMS developers
"""mW ice/liquid in LAMMPS with fix plumed loading the SEAMS action.

usage: e2e_mw.py BUILD [--cells C] [--steps S] [--record FILE]

Needs lmp (or $LMP) built with the PLUMED package. PLUMED_KERNEL names the
libplumedKernel.so that libSeamsPlumed.so was built against; unset, it is
the one beside plumed on PATH. The
SEAMS labels are checked against the CLI on SEAMS's own dump and on the
LAMMPS dump of the same steps.
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
LAMMPS = HERE.parent / "lammps"
sys.path.insert(0, str(HERE.parent))
import e2e_common as e2e  # noqa: E402


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("build")
    ap.add_argument("--cells", type=int, default=8)
    ap.add_argument("--steps", type=int, default=200)
    ap.add_argument("--every", type=int, default=10)
    ap.add_argument("--melt", type=int, default=2000)
    ap.add_argument("--record", default="")
    args = ap.parse_args()

    kernel = e2e.plumed_kernel()
    if kernel is None:
        e2e.skip("no libplumedKernel.so in PLUMED_KERNEL or beside plumed on PATH")
    lmp = e2e.which("lmp", "LMP")
    if lmp is None:
        e2e.skip("LAMMPS (lmp or $LMP), the PLUMED host, was not found")
    styles = subprocess.run([lmp, "-h"], capture_output=True, text=True).stdout
    if not re.search(r"\bPLUMED\b", styles):
        e2e.skip("this lmp has no PLUMED package")
    plugin = e2e.find_build_file(args.build, "libSeamsPlumed.so")
    if plugin is None:
        e2e.skip("libSeamsPlumed.so was not built")
    cli = e2e.seams_cli(args.build)

    natoms = 8 * args.cells ** 3
    with tempfile.TemporaryDirectory(prefix="dseams-plumed-") as tmp:
        work = Path(tmp)
        shutil.copy(HERE / "in.mw.lmp", work)
        for name in ("mw_setup.lmp", "mW.sw"):
            shutil.copy(LAMMPS / name, work)
        deck = (HERE / "plumed.mw.dat").read_text()
        deck = deck.replace("@PLUGIN@", str(plugin)).replace("@NATOMS@", str(natoms))
        (work / "plumed.dat").write_text(deck.replace("@EVERY@", str(args.every)))
        env = dict(os.environ, OMP_NUM_THREADS="1", PLUMED_NUM_THREADS="1",
                   PLUMED_KERNEL=kernel)
        cmd = [lmp, "-var", "cells", str(args.cells), "-var", "steps", str(args.steps),
               "-var", "every", str(args.every), "-var", "melt", str(args.melt),
               "-var", "dump", "native.lammpstrj", "-in", "in.mw.lmp", "-log", "log.lammps",
               "-screen", "none"]
        done = subprocess.run(cmd, cwd=work, env=env, capture_output=True, text=True)
        plog = (work / "plumed.log").read_text() if (work / "plumed.log").is_file() else ""
        if done.returncode != 0:
            log = (work / "log.lammps").read_text() if (work / "log.lammps").is_file() else ""
            if "Incompatible API version for PLUMED" in log:
                e2e.skip("this lmp's fix plumed does not accept the PLUMED in PLUMED_KERNEL")
            sys.stderr.write(done.stdout + done.stderr + log[-3000:] + plog[-3000:])
            sys.exit("lmp with fix plumed failed")
        version = re.search(r"PLUMED: Version: (\S+)", plog)
        timing = re.search(r"SEAMS calls (\d+) seconds (\S+) per_call (\S+)", plog)
        if timing is None:
            sys.exit("the SEAMS action did not report its timing")
        insitu = e2e.read_dump(work / "insitu.lammpstrj", "label")
        mine = e2e.offline_labels(cli, work / "insitu.lammpstrj", len(insitu), 3.5, work,
                                  extra=("--type", "1"))
        frames, atom_frames, bad = e2e.compare(insitu, mine)
        native = e2e.offline_labels(cli, work / "native.lammpstrj", len(insitu), 3.5, work,
                                    extra=("--type", "1"))
        _, _, bad_native = e2e.compare(insitu, native)
        steps_in = [step for step, _ in insitu]
        steps_native = [step for step, _ in e2e.read_dump(work / "native.lammpstrj", "id")]
        if steps_in != steps_native:
            sys.exit(f"SEAMS ran at steps {steps_in}, LAMMPS dumped {steps_native}")
        record = {
            "engine": "plumed",
            "version": version.group(1) if version else "unknown",
            "system": f"mW ice Ic/liquid, {natoms} molecules, LAMMPS fix plumed",
            "steps": args.steps,
            "frames": frames,
            "atom_frames": atom_frames,
            "labels": e2e.census(insitu),
            "native_dump_mismatches": len(bad_native),
            "insitu_s_per_frame": float(timing.group(3)),
        }
        e2e.finish(args.record, record, bad)


if __name__ == "__main__":
    main()
