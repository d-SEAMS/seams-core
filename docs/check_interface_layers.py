#!/usr/bin/env python3
"""Interface ticks are crystalline basal bilayers.

The two puckered halves of a bilayer share one index. A bin with fewer
than half the oxygens of the fullest bilayer is not indexed, so the
stray cubic oxygen and the interfacial shell have no layer number.
The numbered centres are about 3.67 angstrom apart.
"""

import sys
from pathlib import Path

import numpy as np

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "scripts"))
sys.path.insert(0, str(ROOT / "docs" / "examples"))

import science_render as render
from render_scenes import basal_bilayers


def check_layers(positions: np.ndarray, labels: np.ndarray):
    basal = np.asarray(render.stacking_axis(positions, labels, np.isin(labels, [1])), dtype=float)
    if float(basal[1]) < 0.0:
        basal = -basal
    ids, centres = basal_bilayers(positions, labels, basal)
    if len(centres) < 4:
        raise SystemExit(f"interface has {len(centres)} bilayers")
    gaps = np.diff(centres)
    mean = float(np.mean(gaps))
    if abs(mean - 3.67) > 0.25:
        raise SystemExit(f"bilayer spacing {mean:.3f} A is not about 3.67 A")
    if float(np.max(np.abs(gaps - mean))) > 0.4:
        raise SystemExit(f"bilayer gaps are uneven: {np.round(gaps, 3)}")
    counts = [int(np.sum(ids == layer)) for layer in range(len(centres))]
    fullest = max(counts)
    if min(counts) < 0.5 * fullest:
        raise SystemExit(f"a numbered bilayer is below half of {fullest}: {counts}")
    if np.any(ids[labels != 1] >= 0):
        raise SystemExit("a non-hexagonal oxygen was given a layer index")
    if np.any(ids[labels == 0] >= 0) or np.any(ids[labels == 3] >= 0):
        raise SystemExit("cubic or interfacial oxygen was indexed")
    # The tick is the median of the indexed members.
    projection = positions @ (basal / np.linalg.norm(basal))
    for layer, centre in enumerate(centres):
        members = np.flatnonzero(ids == layer)
        drawn = float(np.median(projection[members]))
        if abs(drawn - centre) > 1e-9:
            raise SystemExit(f"layer {layer} tick {drawn} is not the bilayer median")
    print(
        f"interface bilayers {len(centres)} spacing {mean:.3f} A "
        f"counts {counts}"
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
