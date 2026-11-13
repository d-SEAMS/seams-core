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
    dx = float(xs.mean() / width - 0.5)
    dy = float(ys.mean() / height - 0.5)
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


def main(argv):
    failed = False
    for name in argv:
        fill, dx, dy, clipped = measure(name)
        print(
            "{} fill={:.3f} dx={:+.3f} dy={:+.3f} clipped={}".format(
                Path(name).name, fill, dx, dy, clipped
            ),
            flush=True,
        )
        minimum = 0.8 if Path(name).parent.name == "nucleation" else 0.7
        if fill < minimum or abs(dx) > 0.05 or abs(dy) > 0.05 or clipped > 12:
            print(
                "frame failed {}: fill {:.3f} (need {:.2f}), centroid ({:+.3f}, {:+.3f}), clipped {}".format(
                    name, fill, minimum, dx, dy, clipped
                ),
                file=sys.stderr,
            )
            failed = True
    return 1 if failed else 0


if __name__ == "__main__":
    raise SystemExit(main(sys.argv[1:]))
