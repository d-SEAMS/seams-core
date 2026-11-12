# Copyright (c) d-SEAMS developers

"""Regenerate the 215 K mW liquids and compare CHILL+ with OVITO.

The velocity seed is 215 in in.equilibrate. One OpenMP thread keeps
the trajectory fixed. Frames land outside the tree, under MW215_DIR
or `mw215-frames`.
"""

from __future__ import annotations

import os
import shutil
import subprocess
import sys
import tempfile
import time
from collections import Counter
from pathlib import Path

import numpy as np

ROOT = Path(__file__).resolve().parents[2]
DECK = Path(__file__).resolve().parent
OUT = Path(os.environ.get("MW215_DIR", "mw215-frames"))
CELLS = (8, 16, 32)
CUTOFF = 3.5
RUNS = 7
TO_OVITO = {0: 2, 1: 1, 2: 0, 3: 3, 4: 4, 5: 5, 6: 0, 7: 2, 8: 1}
OVITO_NAME = {
    0: "other",
    1: "hexagonal",
    2: "cubic",
    3: "interfacial",
    4: "hydrate",
    5: "interfacial-hydrate",
}


def frame_path(cells: int) -> Path:
    return OUT / f"mw-liquid-{cells ** 3}.lammpstrj"


def generate(cells: int) -> Path:
    OUT.mkdir(parents=True, exist_ok=True)
    dest = frame_path(cells)
    env = os.environ.copy()
    env["OMP_NUM_THREADS"] = "1"
    with tempfile.TemporaryDirectory(prefix="mw215-") as work:
        shutil.copy(DECK / "mW.sw", work)
        shutil.copy(DECK / "in.equilibrate", work)
        subprocess.run(
            ["lmp", "-in", "in.equilibrate", "-var", "n", str(cells), "-log", "log.lammps"],
            cwd=work,
            env=env,
            check=True,
        )
        shutil.copy(Path(work) / "liquid.lammpstrj", dest)
    print(f"wrote {dest.name}")
    return dest


def read_dump(path: Path):
    lines = path.read_text().splitlines()
    n = int(lines[3])
    box = []
    for line in lines[5:8]:
        lo, hi = (float(part) for part in line.split()[:2])
        box.append(hi - lo)
    pos = np.empty((n, 3), dtype=np.float64)
    for i, line in enumerate(lines[9 : 9 + n]):
        parts = line.split()
        pos[i] = (float(parts[2]), float(parts[3]), float(parts[4]))
    return pos, np.asarray(box, dtype=np.float64)


def shell_of(pos, box, index: int) -> int:
    delta = pos[index] - pos
    delta -= box * np.round(delta / box)
    d2 = np.einsum("ij,ij->i", delta, delta)
    return int(np.count_nonzero((d2 > 0.0) & (d2 < CUTOFF * CUTOFF)))


def median(values):
    ordered = sorted(values)
    mid = len(ordered) // 2
    if len(ordered) % 2:
        return ordered[mid]
    return 0.5 * (ordered[mid - 1] + ordered[mid])


def time_rule(pos, box, rule: str):
    sys.path.insert(0, str(ROOT / "bindings" / "python"))
    from chill_plus import chill_plus

    labels = None
    samples = []
    for _ in range(RUNS):
        t0 = time.perf_counter()
        labels, _counts = chill_plus(pos, box, rule=rule, cutoff=CUTOFF)
        samples.append(time.perf_counter() - t0)
    return labels, samples


def time_cli(path: Path, command: str):
    seams = os.environ.get("SEAMS_BIN", "seams")
    env = os.environ.copy()
    env["OMP_NUM_THREADS"] = "1"
    samples = []
    text = ""
    for _ in range(RUNS):
        t0 = time.perf_counter()
        proc = subprocess.run(
            [seams, command, str(path), "--type", "1", "--cutoff", "3.5"],
            check=True,
            capture_output=True,
            text=True,
            env=env,
        )
        samples.append(time.perf_counter() - t0)
        text = proc.stdout.strip()
    return samples, text


def time_ovito(path: Path):
    from ovito.io import import_file
    from ovito.modifiers import ChillPlusModifier

    base = import_file(str(path)).compute()
    mod = ChillPlusModifier(cutoff=CUTOFF)
    samples = []
    last = None
    for _ in range(RUNS):
        data = base.clone()
        t0 = time.perf_counter()
        data.apply(mod)
        samples.append(time.perf_counter() - t0)
        last = data
    return np.asarray(last.particles["Structure Type"]), samples


def report(name: str, seams, ovito, pos, box) -> None:
    mapped = np.asarray([TO_OVITO[int(v)] for v in seams], dtype=np.int32)
    n = len(mapped)
    agree = int(np.count_nonzero(mapped == ovito))
    print(f"rule {name}")
    print(f"overall {agree} {n} {agree / n:.6f}")
    print("label both seams ovito recall")
    for code, label in OVITO_NAME.items():
        got = int(np.count_nonzero(mapped == code))
        ref = int(np.count_nonzero(ovito == code))
        both = int(np.count_nonzero((mapped == code) & (ovito == code)))
        recall = both / ref if ref else 1.0
        print(f"{label} {both} {got} {ref} {recall:.6f}")
    bad = np.flatnonzero(mapped != ovito)
    shells = Counter(shell_of(pos, box, int(i)) for i in bad[:4000])
    print(f"disagree {len(bad)}")
    print("disagree_by_shell " + " ".join(f"{k}:{shells[k]}" for k in sorted(shells)))


def compare(path: Path) -> None:
    pos, box = read_dump(path)
    print(f"frame {path.name} n {len(pos)}")
    ovito, ovito_times = time_ovito(path)
    print("ovito_median_s {:.6f}".format(median(ovito_times)))
    for rule, command in (("chill+", "chill-plus"), ("ovito-chill+", "ovito-chill+")):
        labels, samples = time_rule(pos, box, rule)
        cli, text = time_cli(path, command)
        print(f"seams_median_s {rule} {median(samples):.6f}")
        print(f"cli_median_s {command} {median(cli):.6f}")
        print("cli_last " + " ".join(text.split()))
        report(rule, labels, ovito, pos, box)


def main(argv: list[str]) -> int:
    cells = CELLS
    if len(argv) > 1:
        cells = tuple(int(item) for item in argv[1:])
    for n in cells:
        path = generate(n)
        try:
            compare(path)
        except ImportError as exc:
            print(f"comparison skipped: {exc}", file=sys.stderr)
    return 0


if __name__ == "__main__":
    raise SystemExit(main(sys.argv))
