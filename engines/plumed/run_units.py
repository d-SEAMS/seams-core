# Copyright (c) d-SEAMS developers
"""SEAMS under UNITS LENGTH=A must match SEAMS in PLUMED's default nm.

Regression: SeamsCV scaled PLUMED positions by 10 whatever the length
unit, so with UNITS LENGTH=A it saw a lattice a hundred times too large
and found no ice. Skip only when plumed itself is absent.
"""

import os
import shutil
import subprocess
import sys
import tempfile
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
from e2e_common import find_build_file  # noqa: E402

if shutil.which("plumed") is None:
    print("plumed is not installed")
    sys.exit(77)

build = Path(sys.argv[1]) if len(sys.argv) > 1 else Path("build")
plugin = find_build_file(build, "libSeamsPlumed.so")
if plugin is None:
    sys.exit("plumed is installed but libSeamsPlumed.so was not built")

here = Path(__file__).resolve().parent
DECK = """{units}LOAD FILE={plugin}
seams: SEAMS ATOMS=1-64 CUTOFF=3.5
PRINT ARG=seams.q6,seams.ice,seams.ring FILE=colvar STRIDE=1 FMT=%.15g
"""


def run(units):
    with tempfile.TemporaryDirectory(prefix="dseams-units-") as tmp:
        work = Path(tmp)
        (work / "plumed.dat").write_text(DECK.format(units=units, plugin=plugin))
        env = dict(os.environ, PLUMED_NUM_THREADS="1")
        done = subprocess.run(
            ["plumed", "driver", "--plumed", "plumed.dat", "--ixyz", str(here / "ice.xyz")],
            cwd=work, env=env, capture_output=True, text=True)
        if done.returncode != 0:
            sys.stderr.write(done.stdout + done.stderr)
            sys.exit(f"plumed driver failed with {units.strip() or 'nm'}")
        rows = [line.split()[1:] for line in (work / "colvar").read_text().splitlines()
                if line and not line.startswith("#")]
        return rows


nm_rows = run("")
a_rows = run("UNITS LENGTH=A\n")
print("nm", nm_rows, "A", a_rows)
if float(nm_rows[0][1]) != 1.0:
    sys.exit(f"seams.ice is {nm_rows[0][1]} in nm, expected 1")
if nm_rows != a_rows:
    sys.exit("seams.q6, seams.ice or seams.ring changed with UNITS LENGTH=A")
print("UNITS LENGTH=A matches nm")
