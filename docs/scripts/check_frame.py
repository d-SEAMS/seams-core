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
        return 0.0, 1.0, 1.0
    box_w = int(xs.max() - xs.min() + 1)
    box_h = int(ys.max() - ys.min() + 1)
    fill = (box_w * box_h) / float(width * height)
    dx = float(xs.mean() / width - 0.5)
    dy = float(ys.mean() / height - 0.5)
    return fill, dx, dy


def main(argv):
    failed = False
    for name in argv:
        fill, dx, dy = measure(name)
        print(
            "{} fill={:.3f} dx={:+.3f} dy={:+.3f}".format(Path(name).name, fill, dx, dy),
            flush=True,
        )
        if fill < 0.7 or abs(dx) > 0.05 or abs(dy) > 0.05:
            print(
                "frame failed {}: fill {:.3f}, centroid offset ({:+.3f}, {:+.3f})".format(
                    name, fill, dx, dy
                ),
                file=sys.stderr,
            )
            failed = True
    return 1 if failed else 0


if __name__ == "__main__":
    raise SystemExit(main(sys.argv[1:]))
