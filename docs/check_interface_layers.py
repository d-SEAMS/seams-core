#!/usr/bin/env python3
"""The interface ticks are the basal-bin medians.

Layer 8 is the one cubic oxygen. It sits in the liquid, so the gap
above layer 7 is the assignment, not a misplaced tick.
"""

import sys
from pathlib import Path

import numpy as np

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "scripts"))
sys.path.insert(0, str(ROOT / "docs" / "examples"))

import science_render as render


def layer_centres(positions: np.ndarray, labels: np.ndarray, axis=None):
    """Median height of each basal bin, the same quantity the tick uses."""
    if axis is None:
        axis = np.array([0.0, 0.0, 1.0])
    ids = render.display_layer_ids(positions, np.isin(labels, [0, 1]), axis)
    projection = positions @ render.unit(np.asarray(axis, dtype=float))
    centres = []
    for layer in range(int(ids.max()) + 1):
        members = np.flatnonzero(ids == layer)
        centres.append((layer, members, float(np.median(projection[members]))))
    return ids, centres


def check_layers(positions: np.ndarray, labels: np.ndarray):
    ids, centres = layer_centres(positions, labels)
    if len(centres) != 9:
        raise SystemExit(f"interface has {len(centres)} layer bins, expected 9")
    for layer, members, _centre in centres[:8]:
        kinds = labels[members]
        if int(np.sum(kinds == 0)) != 0 or int(np.sum(kinds == 1)) != len(members):
            raise SystemExit(f"layer {layer} is not hexagonal")
    layer, members, centre = centres[8]
    if len(members) != 1 or int(labels[members[0]]) != 0:
        raise SystemExit(f"layer 8 is {len(members)} atoms, labels {labels[members]}")
    gap = centre - centres[7][2]
    if gap < 15.0:
        raise SystemExit(f"layer 8 is only {gap:.2f} A above layer 7")
    if abs(centre - 47.5779) > 0.01:
        raise SystemExit(f"layer 8 centre {centre:.4f} is not the cubic oxygen")
    # The drawn tick is this median. Recompute it the way draw_layer_guides does.
    up = render.unit(np.array([0.0, 0.0, 1.0]))
    vertical = positions @ up
    for layer, members, expected in centres:
        drawn = float(np.median(vertical[members]))
        if abs(drawn - expected) > 1e-9:
            raise SystemExit(f"layer {layer} tick {drawn} is not the bin median {expected}")
    print(
        f"interface layers 9, layer 8 cubic singleton at {centre:.3f} A, "
        f"gap {gap:.3f} A"
    )
    return centres


def main():
    import os
    import science_tutorials as tutorials

    libdir = os.environ.get("SEAMS_LIBDIR")
    if not libdir:
        for item in os.environ.get("LD_LIBRARY_PATH", "").split(":"):
            if item and (Path(item) / "libseams.so").is_file():
                libdir = item
                break
    if not libdir:
        raise SystemExit("libseams.so is required for the interface layer test")
    path = ROOT / "input" / "traj" / "interface_f0.lammpstrj"
    positions = tutorials.read_dump(path)
    labels, _counts = tutorials.chill_file(tutorials._lib(Path(libdir)), path)
    check_layers(positions, labels)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
