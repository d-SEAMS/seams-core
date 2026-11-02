"""seams prints the same frame range whether a one-worker walk carries rings
and cage affiliation across frames or recomputes them (--rebuild,
SEAMS_CARRY=0, or two frame workers).

usage: cli_carry.py SEAMS TRAJ_DIR
"""

import os
import subprocess
import sys
import tempfile

seams, traj = sys.argv[1:3]
env = {k: v for k, v in os.environ.items() if k != "SEAMS_CARRY"}
failed = []


def run(args, extra=None):
    return subprocess.run([seams, "--format", "json"] + args, capture_output=True,
                          text=True, env=dict(env, **(extra or {})), timeout=600)


def same(args, frames):
    carried = run(args)
    rebuilt = run(["--rebuild"] + args)
    unset = run(args, {"SEAMS_CARRY": "0"})
    parallel = run(["--jobs", "2"] + args)
    ok = (carried.returncode == 0 and
          carried.stdout.count('{"schema"') == frames and
          carried.stdout == rebuilt.stdout == unset.stdout == parallel.stdout)
    if not ok:
        failed.append(" ".join(args))
        print("differs:", " ".join(args), carried.returncode, carried.stderr[-500:])


def perturbed(path, out, frames=6):
    """Frame 1 of the mW crystal, then frames in which three atoms drift
    along x, so each step breaks and makes a few bonds."""
    with open(path, encoding="utf-8") as f:
        lines = f.read().splitlines()
    n = int(lines[3])
    head, atoms = lines[:9], [line.split() for line in lines[9:9 + n]]
    lo, hi = (float(v) for v in lines[5].split())
    with open(out, "w", encoding="utf-8") as f:
        for step in range(frames):
            f.write("\n".join(["ITEM: TIMESTEP", str(step)] + head[2:]) + "\n")
            for a in atoms:
                x = float(a[2])
                if a[0] in ("1", "777", "2049"):
                    x = lo + (x + 0.6 * step - lo) % (hi - lo)
                f.write(f"{a[0]} {a[1]} {x} {a[3]} {a[4]}\n")


mw = os.path.join(traj, "mW_cubic.lammpstrj")
liquid = os.path.join(traj, "dump-l.lammpstrj")
for graph in ("seeded", "cutoff", "knn"):
    same(["--frame", "1", "--last", "11", "cages", "--graph", graph, mw], 11)
same(["--frame", "1", "--last", "11", "cages", "--complete", mw], 11)
same(["--frame", "1", "--last", "11", "cages", "--graph", "cutoff",
      "--signature", "hc", mw], 11)
same(["--frame", "1", "--last", "11", "chill-plus", "--tum-layers", mw], 11)
for graph in ("seeded", "knn"):
    same(["--frame", "1", "--last", "30", "-t", "1", "cages", "--graph", graph,
          liquid], 30)
same(["--frame", "1", "--last", "15", "-t", "1", "cages", "--graph", "knn",
      "--signature", "4:6,6:8", liquid], 15)
same(["--frame", "1", "--last", "8", "-t", "1", "fingerprint", "--graph", "knn",
      "--by-atom", liquid], 8)
with tempfile.TemporaryDirectory() as tmp:
    drift = os.path.join(tmp, "drift.lammpstrj")
    perturbed(mw, drift)
    for graph in ("seeded", "cutoff"):
        same(["--frame", "1", "--last", "6", "cages", "--graph", graph, drift], 6)
    same(["--frame", "1", "--last", "6", "cages", "--graph", "cutoff",
          "--signature", "ddc", drift], 6)
sys.exit(1 if failed else 0)
