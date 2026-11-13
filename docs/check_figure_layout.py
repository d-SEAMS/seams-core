#!/usr/bin/env python3
"""Annotation boxes do not overlap, and wrapped centres stay inside the cell."""

import sys
import tempfile
from pathlib import Path

import numpy as np
from PIL import Image

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "scripts"))
sys.path.insert(0, str(ROOT / "docs" / "examples"))

import science_render as render
from render_scenes import (
    BRINE_LEGEND,
    _present_interface,
    annotate,
    cell_corners,
    centres_outside_outline,
    display_label,
    marks_overlap,
    padded_frame,
    wrap_into_cell,
)


def check_overlap_predicate():
    separate = [("a", (0, 0, 10, 10)), ("b", (12, 0, 20, 10))]
    if marks_overlap(separate):
        raise SystemExit("separated boxes were reported as overlapping")
    touching = [("a", (0, 0, 10, 10)), ("b", (10, 0, 20, 10))]
    if marks_overlap(touching):
        raise SystemExit("boxes that only share an edge were reported as overlapping")
    stacked = [("a", (0, 0, 10, 10)), ("b", (5, 5, 15, 15))]
    if marks_overlap(stacked) != [("a", "b")]:
        raise SystemExit("an overlapping pair was missed")
    if display_label("Na+") != "Na\u207a" or display_label("Cl-") != "Cl\u207b":
        raise SystemExit("the ion legend does not use superscripts")


def _blank(path, positions, forward, up, width):
    _centre, half_h, half_v, *_rest = render.project(positions, forward, up)
    height = max(int(round(width * half_v / max(half_h, 1e-6))), 32)
    Image.new("RGB", (width, height), (255, 255, 255)).save(path)
    return half_h, half_v


def check_interface_layout():
    origin = np.zeros(3)
    lengths = np.array([12.0, 18.0, 36.0])
    positions = []
    ids = []
    for index, j in enumerate((1.0, 6.0, 11.0, 16.0)):
        for k in (2.0, 18.0, 34.0):
            positions.append((6.0, j, k))
            ids.append(index)
    positions = np.asarray(positions, dtype=float)
    labels = np.ones(len(positions), dtype=int)
    ids = np.asarray(ids, dtype=int)
    forward = np.array([-1.0, 0.0, 0.0])
    up = np.array([0.0, 0.0, 1.0])
    frame = padded_frame(positions, origin, lengths, forward, up, 0.85)
    legend = [item for item in BRINE_LEGEND if item[0] in {"hexagonal", "liquid"}]
    with tempfile.TemporaryDirectory() as tmp:
        path = Path(tmp) / "interface.png"
        _blank(path, frame, forward, up, 900)
        _present_interface(
            render, path, positions, labels, origin, lengths, forward, up, ids, legend,
            frame=frame, pad_radius=0.85,
        )
    print("interface annotation boxes are disjoint")


def check_brine_layout():
    origin = np.array([0.0, 0.0, 1.0])
    lengths = np.array([20.0, 16.0, 40.0])
    positions = []
    for i in (0.4, 19.6):
        for j in (0.4, 15.6):
            for k in (1.4, 40.6):
                positions.append((i, j, k))
    positions = np.asarray(positions, dtype=float)
    forward = np.array([0.0, 1.0, 0.0])
    up = np.array([0.0, 0.0, 1.0])
    radius = 1.16
    frame = padded_frame(positions, origin, lengths, forward, up, radius)
    ice_z = np.array([1.4, 8.0, 16.0, 24.0, 30.0, 51.0])
    legend = [item for item in BRINE_LEGEND if item[0] in {"hexagonal", "liquid", "Na+", "Cl-"}]
    with tempfile.TemporaryDirectory() as tmp:
        path = Path(tmp) / "brine.png"
        _blank(path, frame, forward, up, 900)
        annotate(
            render, path, positions, forward, up, 900, cell_corners(origin, lengths), legend,
            frame=frame, origin=origin, lengths=lengths, ice_z=ice_z,
        )
    print("brine annotation boxes are disjoint")


def check_brine_wrap():
    origin = np.zeros(3)
    lengths = np.array([10.0, 12.0, 20.0])
    forward = np.array([0.0, 1.0, 0.0])
    up = np.array([0.0, 0.0, 1.0])
    raw = np.array([[-0.4, 1.0, 5.0], [10.2, 3.0, 8.0], [1.0, 2.0, 3.0], [4.0, -0.2, 22.0]])
    before = centres_outside_outline(raw, origin, lengths, forward, up)
    if before == 0:
        raise SystemExit("the shifted brine centres were already inside the outline")
    after = centres_outside_outline(wrap_into_cell(raw, origin, lengths), origin, lengths, forward, up)
    if after:
        raise SystemExit(f"{after} wrapped brine centres project outside the cell outline")
    cache = (
        ROOT / "docs" / ".cache" / "brine-front-230K"
        / "4828e33ecf419c70643b2d466bfa4f283b4138dc" / "T230_m20_r0_tend.xyz"
    )
    if cache.is_file():
        import science_tutorials as tutorials

        box, rows = tutorials.read_xyz(cache)
        kept = [row for row in rows if row[0] in ("O", "Na", "Cl")]
        box_origin = np.array(box[:3], dtype=float)
        box_lengths = np.array(box[3:6], dtype=float)
        points = np.array([[row[1], row[2], row[3]] for row in kept], dtype=float)
        outside = centres_outside_outline(
            wrap_into_cell(points, box_origin, box_lengths),
            box_origin,
            box_lengths,
            forward,
            up,
        )
        if outside:
            raise SystemExit(f"{outside} brine centres project outside the cell outline")
        print(f"brine frame centres stay inside the outline ({len(points)} atoms)")
    else:
        print("brine frame cache is absent; the shifted-cell wrap was checked")


def main():
    check_overlap_predicate()
    check_interface_layout()
    check_brine_layout()
    check_brine_wrap()
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
