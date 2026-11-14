# Copyright (c) d-SEAMS developers
"""Run plumed driver on ice.xyz. Skip only when plumed itself is absent."""

import os
import shutil
import subprocess
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
from e2e_common import find_build_file  # noqa: E402

if shutil.which("plumed") is None:
    print("plumed is not installed")
    sys.exit(77)

build = Path(sys.argv[1]) if len(sys.argv) > 1 else Path("build")
plugin = find_build_file(build, "libSeamsPlumed.so")
if plugin is None:
    print("plumed is installed but libSeamsPlumed.so was not built")
    sys.exit(1)

here = Path(__file__).resolve().parent
deck = here / "plumed.dat"
xyz = here / "ice.xyz"
text = deck.read_text().replace("./libSeamsPlumed.so", str(plugin))
work = here / "plumed.driver.dat"
colvar = here / "colvar"
work.write_text(text)
try:
    env = os.environ.copy()
    env["PLUMED_NUM_THREADS"] = env.get("PLUMED_NUM_THREADS", "1")
    done = subprocess.run(
        ["plumed", "driver", "--plumed", str(work), "--ixyz", str(xyz)],
        check=False,
        capture_output=True,
        text=True,
        cwd=here,
        env=env,
    )
    if done.returncode != 0:
        sys.stderr.write(done.stdout)
        sys.stderr.write(done.stderr)
        sys.exit(done.returncode)
    if not colvar.is_file():
        sys.stderr.write(done.stdout)
        sys.stderr.write(done.stderr)
        sys.exit("plumed driver wrote no colvar")
    ice = None
    for line in colvar.read_text().splitlines():
        if not line or line.startswith("#"):
            continue
        fields = line.split()
        if len(fields) < 3:
            continue
        ice = float(fields[2])
    if ice is None or abs(ice - 1.0) > 1e-6:
        sys.stderr.write(colvar.read_text())
        sys.exit(f"seams.ice is {ice}, expected 1")
    print(done.stdout)
    print("seams.ice", ice)
finally:
    work.unlink(missing_ok=True)
    colvar.unlink(missing_ok=True)
