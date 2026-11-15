#!/usr/bin/env python3
"""Cage census: d-SEAMS cmdCages against SQQ 0.5.8.

ClavisCage's clavis_engine module is not on PyPI (the package named
clavis is a different client). SQQ 0.5.8 reports 5^12 and 5^12 6^2
cages. sI, sII and sH are replicated to about 120000 molecules.
cmdCages is the seams cages command. SQQ's python cage search does not
finish a 120000-molecule cell in a sample, so its repeated timing stays
on the unit cell, where the counts are compared.
"""

from __future__ import annotations

import os
import subprocess
import sys
import time
from pathlib import Path

import numpy as np

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "scripts"))
from h2h_cells import replicate_sH  # noqa: E402

SEAMS = ROOT / "bbdir-release" / "seams"
N = 20
SIG = "512+51262+51264+51268+sh"


def pin_cores(width: int) -> str:
    return "0" if width == 1 else "0-7"


def read_dump(path: Path):
    lines = path.read_text().splitlines()
    n = int(lines[3])
    bounds = []
    tilts = []
    for row in lines[5:8]:
        parts = row.split()
        bounds.append((float(parts[0]), float(parts[1])))
        tilts.append(float(parts[2]) if len(parts) > 2 else 0.0)
    pos = []
    for line in lines[9 : 9 + n]:
        parts = line.split()
        pos.append((float(parts[-3]), float(parts[-2]), float(parts[-1])))
    return np.asarray(pos, dtype=np.float64), bounds, tilts


def write_ortho(path: Path, pos: np.ndarray, lengths):
    n = len(pos)
    lines = [
        "ITEM: TIMESTEP", "0",
        "ITEM: NUMBER OF ATOMS", str(n),
        "ITEM: BOX BOUNDS pp pp pp",
        f"0 {lengths[0]}", f"0 {lengths[1]}", f"0 {lengths[2]}",
        "ITEM: ATOMS id type x y z",
    ]
    for i, (x, y, z) in enumerate(pos, start=1):
        lines.append(f"{i} 1 {x:.8f} {y:.8f} {z:.8f}")
    path.write_text("\n".join(lines) + "\n")


def replicate_ortho(src: Path, nx, ny, nz, dest: Path):
    pos, bounds, _tilts = read_dump(src)
    lengths = [hi - lo for lo, hi in bounds]
    blocks = []
    for ix in range(nx):
        for iy in range(ny):
            for iz in range(nz):
                shift = np.array([ix * lengths[0], iy * lengths[1], iz * lengths[2]])
                blocks.append(pos + shift)
    big = np.vstack(blocks)
    write_ortho(dest, big, [nx * lengths[0], ny * lengths[1], nz * lengths[2]])
    return len(big)


def parse_cages(text: str) -> dict:
    found = {}
    for line in text.splitlines():
        if "signature" not in line or "cages" not in line:
            continue
        parts = line.split()
        sig = parts[parts.index("signature") + 1]
        n = int(parts[parts.index("cages") + 1])
        found[sig] = n
    return found


def seams_cages(dump: Path, width: int) -> dict:
    env = os.environ.copy()
    env["LD_LIBRARY_PATH"] = str(ROOT / "bbdir-release" / "src")
    env.pop("PLUMED_KERNEL", None)
    cmd = [
        "taskset", "-c", pin_cores(width),
        str(SEAMS), "--jobs", str(width), "--cutoff", "3.5",
        "--graph", "cutoff", "--signature", SIG, "cages", str(dump),
    ]
    done = subprocess.run(cmd, cwd=ROOT, env=env, text=True, capture_output=True)
    if done.returncode != 0:
        raise RuntimeError((done.stderr or done.stdout)[-400:])
    return parse_cages(done.stdout)


def sqq_cages(pos: np.ndarray, lengths, width: int) -> dict:
    from sqq.api import analyze_frame, load_config
    from sqq.models import Atom, Frame

    pos_nm = np.asarray(pos, dtype=np.float64) / 10.0
    box = np.asarray(lengths, dtype=np.float64) / 10.0
    atoms = [
        Atom(index=i, resid=i + 1, resname="SOL", atomname="O", atomid=i + 1,
             xyz=pos_nm[i], molecule_id=i + 1)
        for i in range(len(pos_nm))
    ]
    frame = Frame(name="cell", atoms=atoms, box=box)
    cfg = load_config({
        "engine": "py",
        "graph": {"mode": "oo", "oo_cutoff_nm": 0.35},
        "ring": {"sizes": [4, 5, 6], "definition": "shortest_path", "chordless": True},
        "half_cage": {"enabled": False},
        "quasi_cage": {"enabled": False},
        "cage": {"enabled": True},
        "order_parameter": {"enabled": []},
        "parallel": {"workers": width, "math_threads": 1},
    }, engine="py")
    result = analyze_frame(frame, cfg, engine="py")
    counts = {}
    for cage in result.cages:
        counts[cage.cage_type] = counts.get(cage.cage_type, 0) + 1
    return counts


def stats(samples):
    arr = np.asarray(samples, dtype=np.float64)
    q1, med, q3 = np.percentile(arr, [25, 50, 75])
    return float(med), float(q3 - q1), float(arr.min())


def time_calls(fn, n=N):
    fn()
    out = []
    for _ in range(n):
        t0 = time.perf_counter_ns()
        fn()
        out.append((time.perf_counter_ns() - t0) / 1e6)
    return out


def main():
    width = int(sys.argv[1]) if len(sys.argv) > 1 else 1
    print(f"cores {os.cpu_count()}", flush=True)
    print(f"pinned_threads {width}", flush=True)
    print(f"samples {N}", flush=True)
    print("sqq 0.5.8", flush=True)
    print(
        "note ClavisCage clavis_engine is not installable here; "
        "pip's clavis package is a different client",
        flush=True,
    )
    tmp = Path(os.environ.get("H2H_CAGE_DIR", "bbdir-release/h2h-cages"))
    tmp.mkdir(exist_ok=True)
    units = {
        "sI": ROOT / "input/traj/genice_sI.lammpstrj",
        "sII": ROOT / "input/traj/genice_sII.lammpstrj",
        "sH": ROOT / "input/traj/genice_sH.lammpstrj",
    }
    # Unit-cell agreement, then the large cmdCages timing.
    for name, path in units.items():
        pos, bounds, tilts = read_dump(path)
        lengths = [hi - lo for lo, hi in bounds]
        ours = seams_cages(path, 1)
        # SQQ stores an orthorhombic box. sH is tilted, so the unit
        # comparison for SQQ is sI and sII only.
        if name == "sH":
            print(f"check {name} dseams {ours}", flush=True)
            continue
        theirs = sqq_cages(pos, lengths, 1)
        print(f"check {name} dseams {ours} sqq {theirs}", flush=True)
        expect = {
            "512": ours.get("5:12", 0),
            "51262": ours.get("5:12,6:2", 0),
            "51264": ours.get("5:12,6:4", 0),
        }
        for key, value in expect.items():
            if theirs.get(key, 0) != value:
                raise SystemExit(f"{name} {key} dseams {value} sqq {theirs.get(key, 0)}")

    large = {}
    n = replicate_ortho(units["sI"], 14, 14, 14, tmp / "sI.lammpstrj")
    large["sI"] = (tmp / "sI.lammpstrj", n)
    n = replicate_ortho(units["sII"], 10, 10, 9, tmp / "sII.lammpstrj")
    large["sII"] = (tmp / "sII.lammpstrj", n)
    pos, box9 = replicate_sH(7, 14, 9)
    # Write sH with its tilt so cmdCages sees the triclinic cell.
    sh = tmp / "sH.lammpstrj"
    ox, oy, oz, lx, ly, lz, xy, xz, yz = box9
    lines = [
        "ITEM: TIMESTEP", "0",
        "ITEM: NUMBER OF ATOMS", str(len(pos)),
        "ITEM: BOX BOUNDS xy xz yz pp pp pp",
        f"{ox + min(0, xy, xz, xy + xz)} {ox + lx + max(0, xy, xz, xy + xz)} {xy}",
        f"{oy + min(0, yz)} {oy + ly + max(0, yz)} {xz}",
        f"{oz} {oz + lz} {yz}",
        "ITEM: ATOMS id type x y z",
    ]
    for i, (x, y, z) in enumerate(pos, start=1):
        lines.append(f"{i} 1 {x:.8f} {y:.8f} {z:.8f}")
    sh.write_text("\n".join(lines) + "\n")
    large["sH"] = (sh, len(pos))

    wins = True
    for name, (dump, molecules) in large.items():
        print(f"frame {name} molecules {molecules}", flush=True)

        def run(dump=dump):
            return seams_cages(dump, width)

        samples = time_calls(run)
        med, iqr, amin = stats(samples)
        print(
            f"arm dseams frame {name} threads {width} molecules {molecules} "
            f"median_ms {med:.4f} iqr_ms {iqr:.4f} min_ms {amin:.4f} "
            f"census {run()}",
            flush=True,
        )
    # SQQ repeated timing stays on the 46-molecule sI cell. A cage search
    # there already costs minutes the first time and does not reach 120000.
    pos, bounds, _ = read_dump(units["sI"])
    lengths = [hi - lo for lo, hi in bounds]

    def sqq():
        return sqq_cages(pos, lengths, width)

    # The python cage walk is minutes per 46-molecule cell. Four samples
    # at one thread; the eight-thread repeat is a single sample.
    n_sqq = 4 if width == 1 else 1
    samples = time_calls(sqq, n=n_sqq)
    med, iqr, amin = stats(samples)
    print(
        f"arm sqq frame sI-unit threads {width} molecules {len(pos)} "
        f"median_ms {med:.4f} iqr_ms {iqr:.4f} min_ms {amin:.4f} "
        f"samples {n_sqq} note the python cage walk does not scale to 120000",
        flush=True,
    )
    if not wins:
        raise SystemExit(1)


if __name__ == "__main__":
    main()
