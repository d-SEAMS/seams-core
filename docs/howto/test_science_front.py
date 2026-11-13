#!/usr/bin/env python3
"""Run science-front on the small interface frame and the mW census.

Exit 77 when src/science_front.cpp is not in the tree or the
science-front binary was not built. Meson records 77 as a skip.
"""

import os
import subprocess
import sys
from pathlib import Path


def sources_present(root):
    return (root / "src" / "science_front.cpp").is_file() and (
        root / "src" / "include" / "internal" / "science_front.hpp"
    ).is_file()


def find_exe(build):
    candidates = [
        Path(build) / "science-front",
        Path(build) / "src" / "science-front",
    ]
    for path in candidates:
        if path.is_file() and os.access(path, os.X_OK):
            return path
    return None


def run(exe, args, build):
    env = os.environ.copy()
    lib = str(Path(build) / "src")
    env["LD_LIBRARY_PATH"] = lib + os.pathsep + env.get("LD_LIBRARY_PATH", "")
    done = subprocess.run(
        [str(exe), *args],
        cwd=build,
        env=env,
        text=True,
        stdout=subprocess.PIPE,
        stderr=subprocess.PIPE,
    )
    if done.returncode != 0:
        raise SystemExit(done.stderr or done.stdout)
    return done.stderr


def relative_to_build(path, build):
    return os.path.relpath(path, build)


def read_fronts(path):
    rows = []
    with path.open() as handle:
        header = handle.readline().split()
        for line in handle:
            if not line.strip():
                continue
            rows.append(dict(zip(header, line.split())))
    return rows


def read_summary(path):
    rows = {}
    with path.open() as handle:
        header = handle.readline().split()
        for line in handle:
            if not line.strip():
                continue
            fields = dict(zip(header, line.split()))
            rows[fields["definition"]] = fields
    return rows


def check_interface(exe, root, build):
    dump = root / "input" / "traj" / "interface_f0.lammpstrj"
    out = "science-out/interface"
    (Path(build) / out).mkdir(parents=True, exist_ok=True)
    err = run(
        exe,
        ["interface", relative_to_build(dump, build), "--bin", "1", "--out-dir", out],
        build,
    )
    if "interface frame 1" not in err:
        raise SystemExit(err)
    fronts = read_fronts(Path(build) / out / "interface-fronts.csv")
    if len(fronts) != 1 or int(fronts[0]["n"]) != 864:
        raise SystemExit("interface frame: {}".format(fronts))
    print("basal profile on interface_f0: 864 molecules")


def cluster_sizes(err):
    found = {}
    for line in err.splitlines():
        parts = line.split()
        if len(parts) >= 8 and parts[0] == "nucleus" and parts[1] == "frame":
            found[int(parts[2])] = (int(parts[4]), int(parts[6]))
    return found


def check_nucleation(exe, root, build):
    full = root / "data" / "science-frames" / "nucleation.lammpstrj"
    three = root / "input" / "traj" / "nucleation_three.lammpstrj"
    dump = full if full.is_file() else three
    out = "science-out/nucleation"
    (Path(build) / out).mkdir(parents=True, exist_ok=True)
    err = run(
        exe,
        ["nucleation", relative_to_build(dump, build), "--out-dir", out],
        build,
    )
    sizes = cluster_sizes(err)
    # First three frames of the 148-frame mW trajectory.
    expect = {1: (177, 59), 2: (146, 64), 3: (150, 70)}
    for frame, pair in expect.items():
        if frame not in sizes:
            continue
        if sizes[frame] != pair:
            raise SystemExit("frame {} clusters {} != {}".format(frame, sizes[frame], pair))
    if dump == full:
        summary = read_summary(Path(build) / out / "nucleation-summary.csv")
        extended = summary["chill_extended"]
        if extended["upper_bound"] != "1":
            raise SystemExit("critical size was not marked as an upper bound")
        if abs(float(extended["n_star"]) - 153.333333) > 0.02:
            raise SystemExit("n_star {}".format(extended["n_star"]))
        print("census on 148 frames: n* {} upper bound".format(extended["n_star"]))
    else:
        if 1 not in sizes or 2 not in sizes:
            raise SystemExit(err)
        print("census on the first frames of the 148-frame trajectory")


def main():
    root = Path(sys.argv[1])
    build = sys.argv[2]
    if not sources_present(root):
        print("science-front skipped: science sources are absent")
        raise SystemExit(77)
    exe = find_exe(build)
    if exe is None:
        print("science-front skipped: science-front was not built")
        raise SystemExit(77)
    check_interface(exe, root, build)
    check_nucleation(exe, root, build)
    print("science-front ok")


if __name__ == "__main__":
    main()
