#!/usr/bin/env python3
"""After the 245/255 K pair exits, write their CSVs and start the next pass.

Next pass: 235 K, 250 K, and a second 245 K seed. Frames stay in the run directory.
The CSV is the only file this script writes into the source tree.
"""

from __future__ import annotations

import argparse
import os
import subprocess
import sys
import time
from pathlib import Path

ROOT = Path(__file__).resolve().parents[2]
HERE = Path(__file__).resolve().parent
RUNS = Path(os.environ.get("BRINE_DIR", "brine-runs"))
CSV_DIR = ROOT / "paper" / "tex" / "figures" / "data" / "gate"


def lmp_alive() -> bool:
    proc = subprocess.run(["pgrep", "-f", "lmp -sf omp"], capture_output=True, text=True)
    return proc.returncode == 0 and bool(proc.stdout.strip())


def write_csv(tag: str, brine_front: str) -> Path | None:
    dump = RUNS / tag / "dump.lammpstrj"
    if not dump.is_file():
        print(f"no dump for {tag}", flush=True)
        return None
    out = CSV_DIR / f"{tag}.csv"
    subprocess.run(
        [
            sys.executable,
            str(HERE / "front_csv.py"),
            str(dump),
            "--brine-front",
            brine_front,
            "--out",
            str(out),
        ],
        check=True,
    )
    return out


def launch(tag: str, temp: str, seed: str, threads: str, brine_front: str) -> subprocess.Popen:
    log = open(RUNS / f"{tag}.driver.log", "w")
    cmd = [
        sys.executable,
        str(HERE / "coexist.py"),
        "--run",
        "--threads",
        threads,
        "--melt-ps",
        "20",
        "--prod-ps",
        "3000",
        "--dump-ps",
        "20",
        "--molality",
        "0",
        "--tag",
        tag,
        "--temp",
        temp,
        "--seed",
        seed,
    ]
    print("start", " ".join(cmd), flush=True)
    return subprocess.Popen(cmd, stdout=log, stderr=subprocess.STDOUT)


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--brine-front", required=True)
    parser.add_argument("--no-wait", action="store_true")
    args = parser.parse_args()
    if not args.no_wait:
        while lmp_alive():
            time.sleep(30)
    for tag in ("pure-245", "pure-255"):
        write_csv(tag, args.brine_front)
    procs = [
        launch("pure-245-b", "245", "2452", "4", args.brine_front),
        launch("pure-235", "235", "2351", "2", args.brine_front),
        launch("pure-250", "250", "2501", "2", args.brine_front),
    ]
    codes = [p.wait() for p in procs]
    for tag in ("pure-245-b", "pure-235", "pure-250"):
        write_csv(tag, args.brine_front)
    print("codes", codes, flush=True)


if __name__ == "__main__":
    main()
