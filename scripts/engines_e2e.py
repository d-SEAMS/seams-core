# Copyright (c) d-SEAMS developers
"""Run one engine end-to-end check, or tabulate their records.

usage:
  engines_e2e.py run ENGINE BUILD RECORD
  engines_e2e.py table RECORD... --out data/engines-e2e.csv

run writes the record the engine's check prints. A skipped engine (exit
77) gets a record that says so and why, so the table shows it.
"""

import argparse
import csv
import json
import subprocess
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
CHECKS = {
    "lammps": ["engines/lammps/e2e_mw.py"],
    "lammps-np2": ["engines/lammps/e2e_mw.py", "--ranks", "2"],
    "plumed": ["engines/plumed/e2e_mw.py"],
    "gromacs": ["engines/gromacs/e2e_tip4pice.py"],
    "gromacs-plumed": ["engines/gromacs/e2e_tip4pice.py", "--plumed"],
    "openmm": ["engines/openmm/e2e_tip4pice.py"],
}
COLUMNS = ["engine", "version", "steps", "frames", "label_agreement", "insitu_s_per_frame"]


def run(engine, build, record):
    script, *extra = CHECKS[engine]
    Path(record).parent.mkdir(parents=True, exist_ok=True)
    Path(record).unlink(missing_ok=True)
    cmd = [sys.executable, str(ROOT / script), build, *extra, "--record", record]
    done = subprocess.run(cmd, capture_output=True, text=True)
    sys.stdout.write(done.stdout)
    sys.stderr.write(done.stderr)
    if done.returncode == 77:
        why = done.stdout.strip().splitlines()[-1] if done.stdout.strip() else "SKIP"
        Path(record).write_text(json.dumps({"engine": engine, "skipped": why}) + "\n")
        return 0
    return done.returncode


def table(records, out):
    rows = [json.loads(Path(p).read_text()) for p in records]
    Path(out).parent.mkdir(parents=True, exist_ok=True)
    with open(out, "w", newline="") as fh:
        writer = csv.writer(fh, lineterminator="\n")
        writer.writerow(COLUMNS)
        for row in rows:
            if "skipped" in row:
                writer.writerow([row["engine"], row["skipped"], "", "", "", ""])
                continue
            writer.writerow([row["engine"], row["version"], row["steps"], row["frames"],
                             row["label_agreement"], "%.6e" % row["insitu_s_per_frame"]])


def main():
    ap = argparse.ArgumentParser()
    sub = ap.add_subparsers(dest="cmd", required=True)
    one = sub.add_parser("run")
    one.add_argument("engine", choices=sorted(CHECKS))
    one.add_argument("build")
    one.add_argument("record")
    many = sub.add_parser("table")
    many.add_argument("records", nargs="+")
    many.add_argument("--out", required=True)
    args = ap.parse_args()
    if args.cmd == "run":
        sys.exit(run(args.engine, args.build, args.record))
    table(args.records, args.out)


if __name__ == "__main__":
    main()
