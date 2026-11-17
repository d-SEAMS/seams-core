#!/usr/bin/env python3
# Copyright (c) d-SEAMS developers
"""Interface ticks are crystalline basal bilayers.

The two puckered halves of a bilayer share one index. A bin with fewer
than half the oxygens of the fullest bilayer is not indexed, so the
stray cubic oxygen and the interfacial shell have no layer number.
The numbered centers are about 3.67 angstrom apart.
"""

import sys
from pathlib import Path

import numpy as np

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "scripts"))
sys.path.insert(0, str(ROOT / "docs" / "examples"))
from shared_library import built_library  # noqa: E402

import science_render as render
from render_scenes import basal_bilayers, centres_outside_outline, wrap_into_cell


def check_layers(positions: np.ndarray, labels: np.ndarray):
    basal = np.asarray(render.stacking_axis(positions, labels, np.isin(labels, [1])), dtype=float)
    ids, centers = basal_bilayers(positions, labels, basal)
    if len(centers) < 4:
        raise SystemExit(f"interface has {len(centers)} bilayers")
    gaps = np.diff(centers)
    mean = float(np.mean(gaps))
    if abs(mean - 3.67) > 0.25:
        raise SystemExit(f"bilayer spacing {mean:.3f} A is not about 3.67 A")
    if float(np.max(np.abs(gaps - mean))) > 0.4:
        raise SystemExit(f"bilayer gaps are uneven: {np.round(gaps, 3)}")
    counts = [int(np.sum(ids == layer)) for layer in range(len(centers))]
    fullest = max(counts)
    if min(counts) < 0.5 * fullest:
        raise SystemExit(f"a numbered bilayer is below half of {fullest}: {counts}")
    if np.any(ids[labels != 1] >= 0):
        raise SystemExit("a non-hexagonal oxygen was given a layer index")
    if np.any(ids[labels == 0] >= 0) or np.any(ids[labels == 3] >= 0):
        raise SystemExit("cubic or interfacial oxygen was indexed")
    # The tick is the median of the indexed members.
    projection = positions @ (basal / np.linalg.norm(basal))
    for layer, center in enumerate(centers):
        members = np.flatnonzero(ids == layer)
        drawn = float(np.median(projection[members]))
        if abs(drawn - center) > 1e-9:
            raise SystemExit(f"layer {layer} tick {drawn} is not the bilayer median")
    forward = np.array([-1.0, 0.0, 0.0])
    up = np.array([0.0, 0.0, 1.0])
    _forward, right, _up = render.basis(forward, up)
    horiz = [
        float(np.median(positions[ids == layer] @ right))
        for layer in range(len(centers))
    ]
    if any(horiz[i] >= horiz[i + 1] for i in range(len(horiz) - 1)):
        raise SystemExit(f"bilayer index does not ascend left to right: {np.round(horiz, 2)}")
    print(
        f"interface bilayers {len(centers)} spacing {mean:.3f} A "
        f"counts {counts}"
    )
    return centers


def check_wrapped(raw, origin, lengths):
    before = centres_outside_outline(raw, origin, lengths)
    if before == 0:
        raise SystemExit("the interface dump has no center outside the cell")
    wrapped = wrap_into_cell(raw, origin, lengths)
    after = centres_outside_outline(wrapped, origin, lengths)
    if after:
        raise SystemExit(f"{after} wrapped centers project outside the cell outline")
    print(f"wrapped {before} outside centers into the cell")


def main():
    import os
    import science_tutorials as tutorials

    found = built_library(".", "seams", "SEAMS_LIB") if os.environ.get("SEAMS_LIB") else None
    if found is None:
        for item in os.environ.get("LD_LIBRARY_PATH", "").split(os.pathsep):
            if not item:
                continue
            found = built_library(item, "seams")
            if found is not None:
                break
    if found is None:
        raise SystemExit("libseams is required for the interface layer test")
    libdir = found.parent
    path = ROOT / "input" / "traj" / "interface_f0.lammpstrj"
    positions = tutorials.read_dump(path)
    labels, _counts = tutorials.chill_file(tutorials._lib(Path(libdir)), path)
    check_layers(positions, labels)
    origin = np.zeros(3)
    lengths = []
    with path.open() as handle:
        for line in handle:
            if line.startswith("ITEM: BOX BOUNDS"):
                for _ in range(3):
                    lo, hi = (float(item) for item in handle.readline().split()[:2])
                    lengths.append(hi - lo)
                break
    check_wrapped(positions, origin, np.asarray(lengths))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
