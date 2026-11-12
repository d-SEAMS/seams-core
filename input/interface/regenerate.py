# Copyright (c) d-SEAMS developers

"""Regenerate a 1 ns mW trajectory from the shipped interface frame.

The deck is in.growth. LAMMPS writes growth.lammpstrj next to that deck
and this script copies it to data/interface/, which is not committed.
One OpenMP thread keeps the thermostat on a single core.
"""

from __future__ import annotations

import os
import shutil
import subprocess
import tempfile
from pathlib import Path

ROOT = Path(__file__).resolve().parents[2]
DECK = Path(__file__).resolve().parent
OUT = ROOT / "data" / "interface" / "growth.lammpstrj"


def lammps() -> str:
    for name in ("lmp", "lmp_serial", "lammps"):
        found = shutil.which(name)
        if found:
            return found
    raise SystemExit("LAMMPS was not found")


def main() -> None:
    OUT.parent.mkdir(parents=True, exist_ok=True)
    env = os.environ.copy()
    env["OMP_NUM_THREADS"] = "1"
    with tempfile.TemporaryDirectory(prefix="iface-") as work:
        work_path = Path(work)
        shutil.copy(DECK / "in.growth", work_path / "in.growth")
        shutil.copy(ROOT / "input" / "mw215" / "mW.sw", work_path / "mW.sw")
        shutil.copy(ROOT / "input" / "traj" / "interface_f0.lammpstrj",
                    work_path / "interface_f0.lammpstrj")
        subprocess.run(
            [lammps(), "-in", "in.growth", "-log", "log.lammps"],
            cwd=work_path,
            env=env,
            check=True,
        )
        shutil.copy(work_path / "growth.lammpstrj", OUT)
    print(f"wrote {OUT}")


if __name__ == "__main__":
    main()
