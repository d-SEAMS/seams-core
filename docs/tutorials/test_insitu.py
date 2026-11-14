#!/usr/bin/env python3
"""Run the in-situ hook inside LAMMPS or OpenMM.

Exit 77 when neither library is available, which meson records as a
skip. A library that is present has to label every atom of the
eight-molecule cubic cell as cubic.
"""

import importlib.util
import os
import shutil
import subprocess
import sys
from pathlib import Path


def have_lammps():
    pkg = shutil.which("pkg-config")
    if pkg and subprocess.call([pkg, "--exists", "lammps"]) == 0:
        return True
    return bool(shutil.which("lmp") or shutil.which("lmp_serial"))


def have_openmm():
    return importlib.util.find_spec("openmm") is not None


def ice_atoms(root):
    for line in (root / "engines" / "lammps" / "ice.data").read_text().splitlines():
        if line.endswith(" atoms"):
            return int(line.split()[0])
    raise SystemExit("ice.data has no atom count")


def run_lammps(root, build):
    hook = Path(build) / "engines" / "run_hook"
    if not hook.is_file():
        hook = Path(build) / "run_hook"
    if not hook.is_file():
        return False
    deck = root / "engines" / "lammps" / "in.dseams"
    env = os.environ.copy()
    lib = str(Path(build) / "src")
    env["LD_LIBRARY_PATH"] = lib + os.pathsep + env.get("LD_LIBRARY_PATH", "")
    try:
        text = subprocess.check_output(
            [str(hook), str(deck)], env=env, text=True, cwd=root / "engines" / "lammps",
            stderr=subprocess.STDOUT,
        )
    except subprocess.CalledProcessError as exc:
        detail = exc.output or ""
        if exc.returncode == 127 or "cannot open shared object file" in detail:
            print("lammps skipped: {}".format(detail.strip().splitlines()[-1] if detail.strip() else "library not found"))
            return False
        raise
    # "lammps calls N stride 2 ice M"
    parts = text.split()
    calls = int(parts[parts.index("calls") + 1])
    ice = int(parts[parts.index("ice") + 1])
    if calls < 1 or ice != calls * ice_atoms(root):
        raise SystemExit(text)
    print(text.strip())
    return True


def run_openmm(root, build):
    script = root / "engines" / "openmm" / "test_reporter.py"
    env = os.environ.copy()
    src = str(Path(build) / "src")
    lang = str(Path(build) / "bindings")
    env["LD_LIBRARY_PATH"] = os.pathsep.join([src, lang, env.get("LD_LIBRARY_PATH", "")])
    env["PYTHONPATH"] = os.pathsep.join(
        [str(root / "adapters" / "python"), env.get("PYTHONPATH", "")]
    )
    env["SEAMS_LANG_LIB"] = str(Path(lang) / "libseams_lang.so")
    text = subprocess.check_output([sys.executable, str(script)], env=env, text=True)
    if "openmm calls" not in text or "ice" not in text:
        raise SystemExit(text)
    print(text.strip())
    return True


def main():
    root = Path(sys.argv[1])
    build = sys.argv[2]
    ran = False
    if have_lammps():
        ran = run_lammps(root, build) or ran
    if have_openmm():
        ran = run_openmm(root, build) or ran
    if not ran:
        print("insitu skipped: LAMMPS and OpenMM are not available")
        raise SystemExit(77)
    print("insitu ok")


if __name__ == "__main__":
    main()
