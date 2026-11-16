#!/usr/bin/env python3
"""Fixture driver for a 64-molecule cell. Not the paper committor.

The paper lane is input/nucleation/in.homogeneous at 4096 molecules.
This driver exists so a short LAMMPS call can be checked. Its table
must not be written over data/committor-table.csv.
"""

from __future__ import annotations

import argparse
import csv
import math
import os
import shutil
import subprocess
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "scripts"))
import committor_basin as basin

BRINE_LMP = (
    "/cvmfs/software.eessi.io/versions/2023.06/software/linux/x86_64/"
    "intel/sapphirerapids/software/LAMMPS/"
    "2Aug2023_update2-foss-2023a-kokkos/bin/lmp"
)
CENSUS_LIB = (
    "/cvmfs/software.eessi.io/versions/2026.06/compat/linux/x86_64/lib64"
)
Z95 = 1.959963984540054
SEEDS = list(range(20261111, 20261131))


def find_lmp(explicit: str) -> str:
    if explicit:
        return explicit
    env = os.environ.get("LAMMPS_BIN", "")
    if env:
        return env
    found = shutil.which("lmp")
    if found:
        return found
    if Path(BRINE_LMP).is_file():
        return BRINE_LMP
    raise SystemExit("lmp not found; the brine-lane binary is absent")


def wilson(successes: int, n: int, z: float = Z95) -> tuple[float, float]:
    if n <= 0:
        return float("nan"), float("nan")
    p = successes / n
    den = 1.0 + z * z / n
    center = (p + z * z / (2.0 * n)) / den
    margin = z * math.sqrt(p * (1.0 - p) / n + z * z / (4.0 * n * n)) / den
    return max(0.0, center - margin), min(1.0, center + margin)


def lmp_env() -> dict[str, str]:
    env = os.environ.copy()
    env.pop("LD_LIBRARY_PATH", None)
    env["OMP_NUM_THREADS"] = "1"
    return env


def census_env(pathway: Path) -> dict[str, str]:
    env = os.environ.copy()
    origin = pathway.resolve().parent
    libs = [
        Path(CENSUS_LIB),
        origin / "src",
        origin / "chemfiles-prefix" / "lib",
        origin / "subprojects" / "readcon-core",
        origin / "subprojects" / "readcon-db",
        Path("/lib/x86_64-linux-gnu"),
    ]
    parts = [str(path) for path in libs if path.is_dir()]
    old = env.get("LD_LIBRARY_PATH", "")
    if old:
        parts.append(old)
    env["LD_LIBRARY_PATH"] = ":".join(parts)
    env["OMP_NUM_THREADS"] = "1"
    return env


def run_lmp(lmp: str, work: Path, args: list[str]) -> None:
    log = work / "lammps.log"
    with log.open("w") as handle:
        proc = subprocess.run(
            [lmp, *args, "-log", "log.lammps", "-screen", "none"],
            cwd=work,
            env=lmp_env(),
            stdout=handle,
            stderr=subprocess.STDOUT,
        )
    if proc.returncode != 0:
        tail = log.read_text(errors="replace")[-1500:]
        raise SystemExit(f"LAMMPS failed in {work}\n{tail}")


def last_frame(path: Path) -> tuple[list[float], list[tuple[int, int, float, float, float]]]:
    box = [0.0, 0.0, 0.0]
    atoms: list[tuple[int, int, float, float, float]] = []
    lines = path.read_text().splitlines()
    i = 0
    while i < len(lines):
        if lines[i].startswith("ITEM: BOX BOUNDS"):
            spans = []
            for k in range(3):
                lo, hi = lines[i + 1 + k].split()[:2]
                spans.append(float(hi) - float(lo))
            box = spans
            i += 4
            continue
        if lines[i].startswith("ITEM: ATOMS"):
            header = lines[i].split()[2:]
            idx = {name: n for n, name in enumerate(header)}
            atoms = []
            i += 1
            while i < len(lines) and not lines[i].startswith("ITEM:"):
                parts = lines[i].split()
                atoms.append(
                    (
                        int(parts[idx["id"]]),
                        int(parts[idx["type"]]),
                        float(parts[idx["x"]]),
                        float(parts[idx["y"]]),
                        float(parts[idx["z"]]),
                    )
                )
                i += 1
            continue
        i += 1
    if not atoms:
        raise SystemExit(f"no atoms in {path}")
    return box, atoms


def diamond8() -> list[tuple[float, float, float]]:
    a = 2.75 * 4.0 / math.sqrt(3.0)
    basis = ((0.0, 0.0, 0.0), (0.25, 0.25, 0.25))
    fcc = ((0.0, 0.0, 0.0), (0.5, 0.5, 0.0), (0.5, 0.0, 0.5), (0.0, 0.5, 0.5))
    return [
        tuple((shift[k] + site[k]) * a for k in range(3))
        for shift in fcc
        for site in basis
    ]


def imprint(box: list[float], atoms: list[tuple[int, int, float, float, float]]) -> list[tuple[int, int, float, float, float]]:
    center = [0.5 * length for length in box]
    ranked = sorted(
        atoms,
        key=lambda atom: sum((atom[2 + k] - center[k]) ** 2 for k in range(3)),
    )
    chosen = ranked[:8]
    crystal = diamond8()
    com_c = [sum(point[k] for point in crystal) / 8.0 for k in range(3)]
    com_a = [sum(atom[2 + k] for atom in chosen) / 8.0 for k in range(3)]
    placed = {
        atom[0]: tuple(crystal[i][k] - com_c[k] + com_a[k] for k in range(3))
        for i, atom in enumerate(chosen)
    }
    out = []
    for atom in atoms:
        if atom[0] in placed:
            x, y, z = placed[atom[0]]
            out.append((atom[0], atom[1], x, y, z))
        else:
            out.append(atom)
    return out


def write_data(path: Path, box: list[float], atoms: list[tuple[int, int, float, float, float]]) -> None:
    lines = [
        "Shared committor configuration",
        "",
        f"{len(atoms)} atoms",
        "1 atom types",
        "",
        f"0.0 {box[0]:.8f} xlo xhi",
        f"0.0 {box[1]:.8f} ylo yhi",
        f"0.0 {box[2]:.8f} zlo zhi",
        "",
        "Masses",
        "",
        "1 18.015",
        "",
        "Atoms",
        "",
    ]
    for atom in sorted(atoms):
        lines.append(f"{atom[0]} {atom[1]} {atom[2]:.8f} {atom[3]:.8f} {atom[4]:.8f}")
    path.write_text("\n".join(lines) + "\n")


def series_of(frames: Path) -> list[float]:
    rows = list(csv.DictReader(frames.open(), delimiter="\t"))
    return [float(row["n_chill_core"]) for row in rows]


def write_tables(out: Path, summary: Path, rows: list[dict], commit: str, command: str) -> None:
    note = f"# commit {commit}\n# command {command}\n"
    fields = ["replica", "molecules", "barrier", "final_nmax", "fate", "cut_at_2n", "commit"]
    out.parent.mkdir(parents=True, exist_ok=True)
    with out.open("w", newline="") as handle:
        handle.write(note)
        writer = csv.DictWriter(handle, fieldnames=fields, lineterminator="\n")
        writer.writeheader()
        writer.writerows(rows)
    used = [row for row in rows if row["fate"] != "truncated"]
    committed = sum(row["fate"] == "committed" for row in used)
    lo, hi = wilson(committed, len(used))
    fraction = committed / len(used) if used else float("nan")
    with summary.open("w", newline="") as handle:
        handle.write(note)
        writer = csv.DictWriter(
            handle,
            fieldnames=[
                "n_replicas",
                "n_used",
                "n_committed",
                "p_commit",
                "p_lo",
                "p_hi",
                "commit",
            ],
            lineterminator="\n",
        )
        writer.writeheader()
        writer.writerow(
            {
                "n_replicas": str(len(rows)),
                "n_used": str(len(used)),
                "n_committed": str(committed),
                "p_commit": f"{fraction:.6f}",
                "p_lo": f"{lo:.6f}",
                "p_hi": f"{hi:.6f}",
                "commit": commit,
            }
        )


def self_test() -> None:
    lo, hi = wilson(0, 20)
    if not (lo == 0.0 and 0.0 < hi < 0.2):
        raise SystemExit(f"wilson zero {lo} {hi}")
    lo, hi = wilson(20, 20)
    if not (0.8 < lo < 1.0 and hi == 1.0):
        raise SystemExit(f"wilson one {lo} {hi}")
    deck = (ROOT / "input" / "nucleation" / "in.committor").read_text()
    body = deck.split("read_data", 1)[-1]
    if "2N" in body or "halt" in body:
        raise SystemExit("the run stops on a cluster cap")
    if "${seed}" not in deck or "dump 1" not in deck:
        raise SystemExit("the deck does not refresh velocities or dump frames")
    print("self-test ok")


def head_commit() -> str:
    done = subprocess.run(
        ["git", "rev-parse", "HEAD"], cwd=ROOT, text=True, capture_output=True, check=True
    )
    return done.stdout.strip()


def main() -> None:
    if "--self-test" in sys.argv:
        self_test()
        return
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--lmp", default="")
    parser.add_argument("--pathway", type=Path, required=True)
    parser.add_argument("--out", type=Path, required=True)
    parser.add_argument("--summary", type=Path, required=True)
    parser.add_argument("--replicas", type=int, default=20)
    parser.add_argument("--nrun", type=int, default=5000)
    parser.add_argument("--ndump", type=int, default=500)
    parser.add_argument("--temperature", type=float, default=220.0)
    parser.add_argument("--work", type=Path, default=ROOT / "data" / "science-frames" / "committor-md")
    parser.add_argument("--commit", default="")
    parser.add_argument("--reuse", action="store_true")
    args = parser.parse_args()
    if args.replicas < 20:
        raise SystemExit("replicas must be at least 20")
    if not args.pathway.is_file():
        raise SystemExit(f"missing {args.pathway}")
    lmp = find_lmp(args.lmp)
    seeds = SEEDS[: args.replicas]
    work = args.work
    work.mkdir(parents=True, exist_ok=True)
    shutil.copy(ROOT / "input" / "mw215" / "mW.sw", work / "mW.sw")
    shutil.copy(ROOT / "input" / "nucleation" / "in.homogeneous", work / "in.homogeneous")
    shutil.copy(ROOT / "input" / "nucleation" / "in.committor", work / "in.committor")
    melt = work / "melt"
    melt.mkdir(exist_ok=True)
    shutil.copy(work / "mW.sw", melt / "mW.sw")
    shutil.copy(work / "in.homogeneous", melt / "in.homogeneous")
    if not (args.reuse and (work / "data.committor").is_file()):
        run_lmp(
        lmp,
        melt,
        [
            "-in",
            "in.homogeneous",
            "-var",
            "n",
            "4",
            "-var",
            "nmelt",
            "2000",
            "-var",
            "nquench",
            "0",
            "-var",
            "ndump",
            "1000",
            "-var",
            "nskip",
            "0",
            "-var",
            "Tq",
            "180",
            "-var",
            "seed",
            "20261111",
        ],
    )
    if args.reuse and (work / "data.committor").is_file():
        n_atoms = 0
        for line in (work / "data.committor").read_text().splitlines():
            if line.endswith(" atoms"):
                n_atoms = int(line.split()[0])
        if n_atoms <= 0:
            raise SystemExit("the shared configuration has no atom count")
    else:
        box, atoms = last_frame(melt / "liquid.lammpstrj")
        atoms = imprint(box, atoms)
        n_atoms = len(atoms)
        write_data(work / "data.committor", box, atoms)
    command = " ".join(
        [
            "python3 scripts/committor_md.py",
            f"--lmp {lmp}",
            f"--pathway {args.pathway}",
            f"--replicas {args.replicas}",
            f"--nrun {args.nrun}",
            f"--ndump {args.ndump}",
            f"--temperature {args.temperature}",
        ]
    )
    commit = args.commit or head_commit()
    rows = []
    for seed in seeds:
        replica = work / f"seed-{seed}"
        replica.mkdir(exist_ok=True)
        shutil.copy(work / "mW.sw", replica / "mW.sw")
        shutil.copy(work / "in.committor", replica / "in.committor")
        shutil.copy(work / "data.committor", replica / "data.committor")
        frames = replica / "frames.tsv"
        if not (args.reuse and frames.is_file()):
            run_lmp(
            lmp,
            replica,
            [
                "-in",
                "in.committor",
                "-var",
                "datafile",
                "data.committor",
                "-var",
                "nrun",
                str(args.nrun),
                "-var",
                "ndump",
                str(args.ndump),
                "-var",
                "Tq",
                str(args.temperature),
                "-var",
                "seed",
                str(seed),
            ],
        )
        if not (args.reuse and frames.is_file()):
            subprocess.run(
            [
                str(args.pathway),
                str(replica / "quench.lammpstrj"),
                "--out",
                str(frames),
                "--summary",
                str(replica / "summary.tsv"),
                "--type",
                "1",
                "--threads",
                "1",
            ],
            cwd=ROOT,
            env=census_env(args.pathway),
            check=True,
        )
        sizes = series_of(frames)
        call = basin.committor(sizes, n_atoms)
        rows.append(
            {
                "replica": str(seed),
                "molecules": str(n_atoms),
                **call,
                "commit": commit,
            }
        )
        print(f"seed {seed} fate {call['fate']} final {call['final_nmax']}", flush=True)
    write_tables(args.out, args.summary, rows, commit, command)
    print(args.summary.read_text())


if __name__ == "__main__":
    main()
