#!/usr/bin/env python3
"""Compare primitive ring counts with OVITO FindRings at sizes 4 through 10.

Exit 77 when OVITO cannot be imported. The bond cutoff is the one the
ring benchmark uses.
"""

from __future__ import annotations

import os
import subprocess
import sys

CUTOFF = 3.5
DEPTH = 10
SIZES = range(4, 11)
DUMP = os.path.join("input", "traj", "interface_f0.lammpstrj")


def prepare_ovito() -> bool:
    site = os.environ.get("SEAMS_OVITO_SITE", "/tmp/rings-venv/lib/python3.13/site-packages")
    if os.path.isdir(site) and site not in sys.path:
        sys.path.insert(0, site)
    extra = []
    gl = "/usr/lib/x86_64-linux-gnu"
    qt = os.path.join(site, "PySide6", "Qt", "lib")
    for path in (gl, qt):
        if os.path.isdir(path):
            extra.append(path)
    # The dynamic loader reads LD_LIBRARY_PATH at process start.
    if extra and os.environ.get("SEAMS_OVITO_LIBS") != "1":
        current = os.environ.get("LD_LIBRARY_PATH", "")
        os.environ["LD_LIBRARY_PATH"] = ":".join(extra + ([current] if current else []))
        os.environ["SEAMS_OVITO_LIBS"] = "1"
        os.environ.setdefault("QT_QPA_PLATFORM", "offscreen")
        os.execv(sys.executable, [sys.executable, *sys.argv])
    os.environ.setdefault("QT_QPA_PLATFORM", "offscreen")
    try:
        import ovito  # noqa: F401

        return True
    except ImportError:
        return False


def load_frame(path: str):
    import numpy as np

    lines = open(path).read().splitlines()
    i = 0
    while "NUMBER OF ATOMS" not in lines[i]:
        i += 1
    n = int(lines[i + 1])
    i += 2
    while not lines[i].startswith("ITEM: BOX"):
        i += 1
    i += 1
    lo, hi = [], []
    for _ in range(3):
        a, b = map(float, lines[i].split()[:2])
        lo.append(a)
        hi.append(b)
        i += 1
    header = lines[i].split()[2:]
    i += 1
    cols = {name: k for k, name in enumerate(header)}
    pts = np.empty((n, 3))
    for k in range(n):
        tok = lines[i].split()
        i += 1
        pts[k] = [float(tok[cols[a]]) for a in "xyz"]
    return pts, lo, hi


def dseams_counts(bench: str, dump: str) -> dict[int, int]:
    env = os.environ.copy()
    root = os.path.abspath(os.path.join(os.path.dirname(bench), ".."))
    env["LD_LIBRARY_PATH"] = os.path.join(root, "src")
    env["OMP_NUM_THREADS"] = "1"
    env.pop("SEAMS_RING_PROFILE", None)
    proc = subprocess.run(
        [bench, "--dump", dump, "--frame", "1", "--depth", str(DEPTH), "--threads", "1"],
        env=env,
        check=True,
        capture_output=True,
        text=True,
    )
    counts = {s: 0 for s in SIZES}
    for line in proc.stdout.splitlines():
        if not line.startswith("dump "):
            continue
        parts = line.split()
        for s in SIZES:
            key = f"n{s}"
            if key in parts:
                counts[s] = int(parts[parts.index(key) + 1])
    return counts


def ovito_counts(pts, lo, hi) -> dict[int, int]:
    import numpy as np
    from ovito.data import DataCollection
    from ovito.modifiers import CreateBondsModifier, FindRingsModifier
    from ovito.pipeline import Pipeline, StaticSource

    os.environ["QT_QPA_PLATFORM"] = "offscreen"
    os.environ.setdefault("OVITO_THREAD_COUNT", "1")
    side = [hi[a] - lo[a] for a in range(3)]
    mat = np.zeros((3, 4))
    mat[:, :3] = np.diag(side)
    mat[:, 3] = lo
    data = DataCollection()
    data.create_cell(mat, pbc=(True, True, True))
    parts = data.create_particles(count=len(pts))
    parts.create_property("Position", data=np.asarray(pts))
    pipe = Pipeline(source=StaticSource(data=data))
    pipe.modifiers.append(CreateBondsModifier(cutoff=CUTOFF))
    pipe.modifiers.append(FindRingsModifier(maximum_ring_size=DEPTH))
    out = pipe.compute()
    hist = out.tables["ring-size-histogram"]
    found = {int(s): int(c) for s, c in zip(hist["Ring Size"][...], hist["Count"][...])}
    return {s: found.get(s, 0) for s in SIZES}


def main() -> None:
    if not prepare_ovito():
        print("ovito is not importable", file=sys.stderr)
        raise SystemExit(77)
    if len(sys.argv) != 2:
        print("usage: test_ovito_ring_counts.py BENCH", file=sys.stderr)
        raise SystemExit(2)
    bench = sys.argv[1]
    dump = os.path.abspath(DUMP)
    pts, lo, hi = load_frame(dump)
    ours = dseams_counts(bench, dump)
    theirs = ovito_counts(pts, lo, hi)
    if ours != theirs:
        print(f"dseams {ours}", file=sys.stderr)
        print(f"ovito  {theirs}", file=sys.stderr)
        raise SystemExit(1)
    print("sizes 4-10", " ".join(f"n{s}={ours[s]}" for s in SIZES))


if __name__ == "__main__":
    main()
