# Copyright (c) d-SEAMS developers
"""Run plumed driver on ice.xyz, or skip when plumed or the plugin is absent."""

import shutil
import subprocess
import sys
from pathlib import Path

if shutil.which("plumed") is None:
    print("plumed is not installed")
    sys.exit(77)

build = Path(sys.argv[1]) if len(sys.argv) > 1 else Path("build")
plugins = list(build.rglob("libSeamsPlumed.so"))
if not plugins:
    print("libSeamsPlumed.so was not built")
    sys.exit(77)

here = Path(__file__).resolve().parent
deck = here / "plumed.dat"
xyz = here / "ice.xyz"
text = deck.read_text().replace("./libSeamsPlumed.so", str(plugins[0]))
work = here / "plumed.driver.dat"
work.write_text(text)
try:
    done = subprocess.run(
        ["plumed", "driver", "--plumed", str(work), "--ixyz", str(xyz), "--length-units", "A"],
        check=False,
        capture_output=True,
        text=True,
    )
finally:
    work.unlink(missing_ok=True)
if done.returncode != 0:
    sys.stderr.write(done.stdout)
    sys.stderr.write(done.stderr)
    sys.exit(done.returncode)
print(done.stdout)
