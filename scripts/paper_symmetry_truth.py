#!/usr/bin/env python3
"""Ground truth for thermal misses on the TIP4P/2005 interface.

The instantaneous equil frame defines two miss sets: oxygens the
classifier calls other while CHILL+ calls hexagonal, and oxygens it
calls other while ice-only PTM calls hexagonal diamond. The same atom
ids are scored again on the time average of the following hop and on a
short TIP4P minimisation of the equil restart. A miss is real crystal
when at least two of the three methods call it ice Ih or ice Ic on that
reference frame.
"""

from __future__ import annotations

import os
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "scripts"))
os.environ.setdefault("QT_QPA_PLATFORM", "offscreen")
os.environ.setdefault("OVITO_THREAD_COUNT", "8")
os.environ.setdefault("OMP_NUM_THREADS", "8")

import numpy as np
import paper_ptm_templates as templates


def read_oxygen_frames(path: Path):
    frames = []
    ids = None
    pos = []
    bounds = []
    mode = None
    with path.open() as handle:
        for line in handle:
            text = line.strip()
            if text.startswith("ITEM: TIMESTEP"):
                if pos:
                    frames.append((np.asarray(ids), np.asarray(pos, dtype=float), np.asarray(bounds)))
                ids = []
                pos = []
                bounds = []
                mode = None
                continue
            if text.startswith("ITEM: BOX"):
                mode = "box"
                continue
            if text.startswith("ITEM: ATOMS"):
                mode = "atoms"
                continue
            if text.startswith("ITEM:"):
                mode = None
                continue
            if mode == "box" and len(bounds) < 3:
                lo, hi = text.split()[:2]
                bounds.append(float(hi) - float(lo))
            elif mode == "atoms":
                parts = text.split()
                if len(parts) < 5 or int(parts[1]) != 1:
                    continue
                ids.append(int(parts[0]))
                pos.append((float(parts[2]), float(parts[3]), float(parts[4])))
    if pos:
        frames.append((np.asarray(ids), np.asarray(pos, dtype=float), np.asarray(bounds)))
    return frames


def lattice_of(box):
    return [[box[0], 0.0, 0.0], [0.0, box[1], 0.0], [0.0, 0.0, box[2]]]


def labels_of(lib, pos, box):
    lat = lattice_of(box)
    ours = templates.classify(lib, pos, lat)
    ptm = templates.ovito_labels(pos, lat, "ice")
    chill = templates.ovito_labels(pos, lat, "chill")
    return ours, ptm, chill


def ice_mask(ours, ptm, chill):
    our_ice = ours >= 0
    ptm_ice = (ptm == 6) | (ptm == 7)
    chill_ice = (chill == 1) | (chill == 2)
    votes = our_ice.astype(int) + ptm_ice.astype(int) + chill_ice.astype(int)
    return our_ice, ptm_ice, chill_ice, votes >= 2


def average_positions(frames):
    ids0, pos0, box0 = frames[0]
    acc = pos0.copy()
    cursor = pos0.copy()
    box_acc = box0.astype(float).copy()
    for ids, pos, box in frames[1:]:
        if not np.array_equal(ids, ids0):
            raise RuntimeError("oxygen ids changed between frames")
        delta = pos - cursor
        delta -= box * np.round(delta / box)
        cursor = cursor + delta
        acc += cursor
        box_acc += box
    n = len(frames)
    mean = acc / n
    box = box_acc / n
    mean -= box * np.floor(mean / box)
    return ids0, mean, box


def report(title, mask, our_ice, ptm_ice, chill_ice, majority):
    n = int(np.count_nonzero(mask))
    print(f"{title} n={n}")
    if n == 0:
        return
    print(f"  classifier ice {int(np.count_nonzero(our_ice[mask]))}")
    print(f"  PTM diamond {int(np.count_nonzero(ptm_ice[mask]))}")
    print(f"  CHILL+ hexagonal-or-cubic {int(np.count_nonzero(chill_ice[mask]))}")
    print(f"  majority crystal {int(np.count_nonzero(majority[mask]))}")
    print(f"  all three {int(np.count_nonzero(our_ice[mask] & ptm_ice[mask] & chill_ice[mask]))}")


def main() -> int:
    equil = read_oxygen_frames(Path(sys.argv[1]))[0]
    quench = read_oxygen_frames(Path(sys.argv[2]))[0]
    averaged = average_positions(read_oxygen_frames(Path(sys.argv[3])))
    lib = templates.open_core(str(ROOT / "bbdir-release/src/libyodaLib.so"))
    sets = {}
    for name, (ids, pos, box) in (
        ("equil", equil),
        ("quench", quench),
        ("average", averaged),
    ):
        if not np.array_equal(ids, equil[0]):
            order = np.argsort(ids)
            ids = ids[order]
            pos = pos[order]
            if not np.array_equal(ids, equil[0]):
                raise RuntimeError(f"{name} ids differ")
        ours, ptm, chill = labels_of(lib, pos, box)
        our_ice, ptm_ice, chill_ice, majority = ice_mask(ours, ptm, chill)
        sets[name] = (ours, ptm, chill, our_ice, ptm_ice, chill_ice, majority)
        print(
            f"{name} oxygens {len(pos)} classifier {int(our_ice.sum())} "
            f"PTM {int(ptm_ice.sum())} CHILL+ {int(chill_ice.sum())} "
            f"majority {int(majority.sum())}"
        )
    ours, ptm, chill = sets["equil"][:3]
    chill_miss = (ours == -1) & (chill == 1)
    ptm_miss = (ours == -1) & (ptm == 7)
    for name, mask in (("CHILL+ hexagonal misses", chill_miss), ("PTM hex-diamond misses", ptm_miss)):
        print(f"== {name} defined on the equil frame ==")
        for frame in ("quench", "average"):
            report(frame, mask, *sets[frame][3:])
    return 0


if __name__ == "__main__":
    sys.exit(main())
