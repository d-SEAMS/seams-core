#!/usr/bin/env python3
# Copyright (c) d-SEAMS developers
"""Fail when a still leaves the structure small or off centre."""

import sys
from pathlib import Path

import numpy as np
from PIL import Image


def measure(path):
    image = np.asarray(Image.open(path).convert("RGB"))
    mask = np.any(image < 248, axis=2)
    height, width = mask.shape
    ys, xs = np.where(mask)
    if len(xs) == 0:
        return 0.0, 1.0, 1.0, 0
    box_w = int(xs.max() - xs.min() + 1)
    box_h = int(ys.max() - ys.min() + 1)
    fill = (box_w * box_h) / float(width * height)
    # The caption sits in the bottom band. Centroid is the drawing, so a
    # scale bar on the left does not drag the structure off centre.
    band = max(1, int(height * 0.78))
    drawing = mask.copy()
    drawing[band:, :] = False
    dys, dxs = np.where(drawing)
    if len(dxs) == 0:
        dxs, dys = xs, ys
        band = height
    # Box centre, not the ink mean. A liquid on one side must not count as off-centre.
    dx = float((dxs.min() + dxs.max()) / 2.0 / width - 0.5)
    dy = float((dys.min() + dys.max()) / 2.0 / band - 0.5)
    return fill, dx, dy, _clipped_highlights(image)


def _clipped_highlights(image):
    """Near-white pixels sitting on a coloured atom, not the page or a ring hole."""
    white = np.all(image >= 250, axis=2)
    chroma = image.max(axis=2) - image.min(axis=2)
    colored = chroma > 16
    near = np.zeros(colored.shape, dtype=np.int16)
    height, width = colored.shape
    for dy in range(-2, 3):
        y0, y1 = max(0, dy), min(height, height + dy)
        sy0, sy1 = max(0, -dy), min(height, height - dy)
        for dx in range(-2, 3):
            x0, x1 = max(0, dx), min(width, width + dx)
            sx0, sx1 = max(0, -dx), min(width, width - dx)
            near[y0:y1, x0:x1] += colored[sy0:sy1, sx0:sx1]
    # A clipped specular sits inside the sphere. A ring hole is mostly white.
    return int(np.count_nonzero(white & (near >= 16)))


def edge_ink(image):
    """Dark pixels on the outer rim. A clipped label or atom lands here."""
    dark = np.any(image < 80, axis=2)
    rim = 6
    return int(np.count_nonzero(dark[:rim, :]) + np.count_nonzero(dark[-rim:, :])
               + np.count_nonzero(dark[:, :rim]) + np.count_nonzero(dark[:, -rim:]))


def column_point_size(path):
    """Legend line height in points when the PNG is printed at its dpi."""
    image = Image.open(path)
    dpi = 300.0
    info = image.info.get("dpi")
    if info:
        dpi = float(info[0] or 300.0)
    array = np.asarray(image.convert("RGB"))
    height = array.shape[0]
    band = array[int(height * 0.72) :]
    ink = np.any(band < 80, axis=2)
    rows = np.where(ink.any(axis=1))[0]
    if len(rows) == 0:
        return 0.0
    breaks = np.where(np.diff(rows) > 4)[0]
    starts = np.concatenate(([0], breaks + 1))
    ends = np.concatenate((breaks, [len(rows) - 1]))
    line = max(int(rows[end] - rows[start] + 1) for start, end in zip(starts, ends))
    return line / dpi * 72.0


def main(argv):
    failed = False
    for name in argv:
        fill, dx, dy, clipped = measure(name)
        points = column_point_size(name)
        rim = edge_ink(np.asarray(Image.open(name).convert("RGB")))
        print(
            "{} fill={:.3f} dx={:+.3f} dy={:+.3f} clipped={} label={:.1f}pt rim={}".format(
                Path(name).name, fill, dx, dy, clipped, points, rim
            ),
            flush=True,
        )
        minimum = 0.72 if Path(name).parent.name == "nucleation" else 0.60
        if fill < minimum or abs(dx) > 0.05 or abs(dy) > 0.05 or clipped > 12 or points < 5.0 or rim > 0:
            print(
                "frame failed {}: fill {:.3f} (need {:.2f}), centroid ({:+.3f}, {:+.3f}), clipped {}, label {:.1f} pt, rim {}".format(
                    name, fill, minimum, dx, dy, clipped, points, rim
                ),
                file=sys.stderr,
            )
            failed = True
    return 1 if failed else 0


if __name__ == "__main__":
    raise SystemExit(main(sys.argv[1:]))
