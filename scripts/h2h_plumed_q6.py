# Copyright (c) d-SEAMS developers
"""PLUMED Q6, LOCAL_AVERAGE and LOCAL_Q6 against d-SEAMS q6 and q6bar.

usage: h2h_plumed_q6.py BUILD [--cells 8 24] [--repeats 8] [--threads 1 8]
                        [--out benchmarks/h2h-plumed-q6.csv]

Needs plumed, lmp (or $LMP) and libSeamsPlumed.so under BUILD. LAMMPS
writes mW ice/liquid frames (8 cells is 4096 molecules, 24 is 110592).
plumed driver then reads the same frames for every method.

PLUMED's neighbour weight is the switching function
CUSTOM FUNC=step(c-x) with D_MAX=c and NOSTRETCH, a hard cutoff at the
d-SEAMS cutoff c. PLUMED's Q6 is sqrt(sum_m |q6m|^2); d-SEAMS includes
the Steinhardt factor 4 pi / 13, so PLUMED values are multiplied by
sqrt(4 pi / 13) before the per-atom comparison. LOCAL_AVERAGE is the
Lechner-Dellago average, d-SEAMS q6bar. d-SEAMS has no l=6 bond
correlation, so LOCAL_Q6 is timed and not compared.

Timing is PLUMED's DETAILED_TIMERS for the actions of one method per
driver run, divided by the frames. Each method runs alone; the methods
are interleaved within each repeat.
"""

import argparse
import csv
import math
import os
import re
import shutil
import statistics
import subprocess
import sys
import tempfile
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
SCALE = math.sqrt(4 * math.pi / 13)
SWITCH = "SWITCH={{CUSTOM FUNC=step({c}-x) R_0=1 D_MAX={c} NOSTRETCH}}"

# label: (deck lines, labels whose timers make up the method)
METHODS = {
    "plumed Q6": (["q6: Q6 SPECIES=1-{n} {sw} MEAN"], ["q6"], "q6.mean"),
    "plumed LOCAL_AVERAGE": (["q6: Q6 SPECIES=1-{n} {sw}",
                              "la: LOCAL_AVERAGE SPECIES=q6 {sw} MEAN"], ["q6", "la"], "la.mean"),
    "plumed LOCAL_Q6": (["q6: Q6 SPECIES=1-{n} {sw}",
                         "lq: LOCAL_Q6 SPECIES=q6 {sw} MEAN"], ["q6", "lq"], "lq.mean"),
    "d-SEAMS q6": (["s6: SEAMS ATOMS=1-{n} CUTOFF={c} MODE=q6"], ["s6"], "s6.q6"),
    "d-SEAMS q6bar": (["sb: SEAMS ATOMS=1-{n} CUTOFF={c} MODE=q6bar"], ["sb"], "sb.q6"),
}
PAIRS = {"plumed Q6": "d-SEAMS q6", "plumed LOCAL_AVERAGE": "d-SEAMS q6bar"}
COMPARED = {**PAIRS, **{b: a for a, b in PAIRS.items()}}

FRAMES_DECK = """include mw_setup.lmp
dump d all custom ${every} frames.lammpstrj id type x y z
dump_modify d sort id format float %.17g
run ${steps}
"""


def run(cmd, cwd, env, what):
    done = subprocess.run(cmd, cwd=cwd, env=env, capture_output=True, text=True)
    if done.returncode != 0:
        sys.stderr.write(done.stdout[-4000:] + done.stderr[-4000:])
        sys.exit(f"{what} failed")
    return done.stdout


def make_frames(lmp, cells, frames, every, melt, ranks, work):
    for name in ("mw_setup.lmp", "mW.sw"):
        shutil.copy(ROOT / "engines" / "lammps" / name, work)
    (work / "in.frames").write_text(FRAMES_DECK)
    cmd = [lmp, "-var", "cells", str(cells), "-var", "melt", str(melt), "-var", "every",
           str(every), "-var", "steps", str(every * (frames - 1)), "-in", "in.frames",
           "-log", "log.lammps", "-screen", "none"]
    if ranks > 1:
        cmd = ["mpirun", "-np", str(ranks), *cmd]
    run(cmd, work, dict(os.environ, OMP_NUM_THREADS="1"), "lmp (frames)")
    lines = (work / "frames.lammpstrj").read_text().splitlines()
    k = 0
    with open(work / "frames.xyz", "w") as out:
        while k < len(lines):
            n = int(lines[k + 3])
            lo_hi = [lines[k + 5 + d].split() for d in range(3)]
            out.write(f"{n}\n" + " ".join(repr(float(b) - float(a)) for a, b in lo_hi) + "\n")
            for line in lines[k + 9:k + 9 + n]:
                _, _, x, y, z = line.split()[:5]
                out.write(f"X {x} {y} {z}\n")
            k += 9 + n
    return n


def driver(plugin, deck_lines, work, env, frames_xyz, tag):
    deck = ["UNITS LENGTH=A", f"LOAD FILE={plugin}", *deck_lines, "DEBUG DETAILED_TIMERS"]
    (work / f"{tag}.dat").write_text("\n".join(deck) + "\n")
    log = run(["plumed", "driver", "--plumed", f"{tag}.dat", "--ixyz", str(frames_xyz),
               "--length-units", "A"], work, env, f"plumed driver ({tag})")
    return log


def action_seconds(log, labels):
    total = 0.0
    for label in labels:
        hit = re.search(rf"PLUMED: 4A\s+\d+ {re.escape(label)}\s+(\d+)\s+(\S+)", log)
        if hit is None:
            sys.exit(f"no timer for action {label}")
        total += float(hit.group(2))
    return total


def read_multicolvar(path):
    """First value column of every DUMPMULTICOLVAR row, frame after frame."""
    values = []
    with open(path) as fh:
        while True:
            head = fh.readline()
            if not head:
                break
            n = int(head)
            fh.readline()
            for _ in range(n):
                values.append(float(fh.readline().split()[4]))
    return values


def read_seams_dump(path, column):
    values = []
    with open(path) as fh:
        lines = iter(fh)
        for line in lines:
            if line.startswith("ITEM: NUMBER OF ATOMS"):
                n = int(next(lines))
            elif line.startswith("ITEM: ATOMS"):
                at = line.split()[2:].index(column)
                for _ in range(n):
                    values.append(float(next(lines).split()[at]))
    return values


def compare_values(plugin, n, cutoff, work, env, frames_xyz):
    sw = SWITCH.format(c=cutoff)
    deck = [
        f"q6: Q6 SPECIES=1-{n} {sw}",
        f"la: LOCAL_AVERAGE SPECIES=q6 {sw}",
        f"s6: SEAMS ATOMS=1-{n} CUTOFF={cutoff} MODE=q6 DUMP=s6.dump",
        f"sb: SEAMS ATOMS=1-{n} CUTOFF={cutoff} MODE=q6bar DUMP=sb.dump",
        "DUMPMULTICOLVAR DATA=q6 FILE=q6.xyz PRECISION=13",
        "DUMPMULTICOLVAR DATA=la FILE=la.xyz PRECISION=13",
        "PRINT ARG=s6.q6,sb.q6 FILE=colvar",
    ]
    for name in ("s6.dump", "sb.dump", "q6.xyz", "la.xyz"):
        (work / name).unlink(missing_ok=True)
    log = driver(plugin, deck, work, env, frames_xyz, "values")
    out = {}
    for plumed_file, seams_file, column, method in (("q6.xyz", "s6.dump", "q6", "plumed Q6"),
                                                    ("la.xyz", "sb.dump", "q6bar",
                                                     "plumed LOCAL_AVERAGE")):
        mine = read_multicolvar(work / plumed_file)
        theirs = read_seams_dump(work / seams_file, column)
        if len(mine) != len(theirs):
            sys.exit(f"{method}: {len(mine)} PLUMED values but {len(theirs)} d-SEAMS values")
        out[method] = max(abs(SCALE * a - b) for a, b in zip(mine, theirs))
        out[PAIRS[method]] = out[method]
    version = re.search(r"PLUMED: Version: (\S+)", log)
    seams = re.search(r"d-SEAMS (\S+), cutoff", log)
    return out, version.group(1) if version else "", seams.group(1) if seams else ""


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawTextHelpFormatter)
    ap.add_argument("build")
    ap.add_argument("--cells", type=int, nargs="+", default=[8, 24])
    ap.add_argument("--frames", type=int, default=4)
    ap.add_argument("--every", type=int, default=50)
    ap.add_argument("--melt", type=int, default=2000)
    ap.add_argument("--ranks", type=int, default=4, help="MPI ranks for the LAMMPS frames")
    ap.add_argument("--cutoff", type=float, default=3.5)
    ap.add_argument("--repeats", type=int, default=8)
    ap.add_argument("--threads", type=int, nargs="+", default=[1, 8])
    ap.add_argument("--work", default="", help="keep frames and decks here")
    ap.add_argument("--out", default=str(ROOT / "benchmarks" / "h2h-plumed-q6.csv"))
    args = ap.parse_args()

    if shutil.which("plumed") is None:
        sys.exit("plumed was not found")
    lmp = os.environ.get("LMP") or shutil.which("lmp")
    if not lmp:
        sys.exit("lmp (or $LMP) was not found")
    plugins = sorted(Path(args.build).rglob("libSeamsPlumed.so"))
    if not plugins:
        sys.exit(f"libSeamsPlumed.so was not built under {args.build}")
    plugin = plugins[0]
    cores = len(os.sched_getaffinity(0))

    rows = []
    with tempfile.TemporaryDirectory(prefix="dseams-h2h-q6-") as tmp:
        base = Path(args.work) if args.work else Path(tmp)
        for cells in args.cells:
            work = base / f"cells{cells}"
            work.mkdir(parents=True, exist_ok=True)
            frames_xyz = work / "frames.xyz"
            if frames_xyz.is_file():
                n = int(frames_xyz.open().readline())
            else:
                n = make_frames(lmp, cells, args.frames, args.every, args.melt, args.ranks, work)
            print(f"{n} molecules, {args.frames} frames", flush=True)
            sw = SWITCH.format(c=args.cutoff)
            for threads in args.threads:
                env = dict(os.environ, PLUMED_NUM_THREADS=str(threads),
                           OMP_NUM_THREADS=str(threads))
                diffs, plumed_version, seams_version = compare_values(
                    plugin, n, args.cutoff, work, env, frames_xyz)
                times = {name: [] for name in METHODS}
                for rep in range(args.repeats):
                    for name, (lines, labels, arg) in METHODS.items():
                        deck = [line.format(n=n, sw=sw, c=args.cutoff) for line in lines]
                        deck.append(f"PRINT ARG={arg} FILE=colvar")
                        log = driver(plugin, deck, work, env, frames_xyz, "time")
                        times[name].append(action_seconds(log, labels) / args.frames)
                    print(f"  threads {threads} repeat {rep + 1}", flush=True)
                for name in METHODS:
                    t = times[name]
                    rows.append({
                        "molecules": n,
                        "frames": args.frames,
                        "method": name,
                        "compared_with": COMPARED.get(name, ""),
                        "plumed_scale": "%.15f" % SCALE if name in PAIRS else "",
                        "max_abs_diff": "%.3e" % diffs[name] if name in diffs else "",
                        "threads": threads,
                        "repeats": len(t),
                        "s_per_frame_median": "%.6e" % statistics.median(t),
                        "s_per_frame_min": "%.6e" % min(t),
                        "s_per_frame_max": "%.6e" % max(t),
                        "cores": cores,
                        "plumed": plumed_version,
                        "dseams": seams_version,
                    })
    Path(args.out).parent.mkdir(parents=True, exist_ok=True)
    with open(args.out, "w", newline="") as fh:
        writer = csv.DictWriter(fh, fieldnames=list(rows[0]), lineterminator="\n")
        writer.writeheader()
        writer.writerows(rows)
    print(Path(args.out).read_text())


if __name__ == "__main__":
    main()
