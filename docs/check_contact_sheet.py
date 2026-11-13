#!/usr/bin/env python3
"""Fail if a docs raster is under 1200 px on its long side."""

import sys
from pathlib import Path

from PIL import Image

NAMES = ("interface.png", "brine.png", "sii.png", "sh.png")


def check(directory: Path):
    for name in NAMES:
        path = directory / name
        if not path.is_file():
            raise SystemExit(f"missing render {path}")
        image = Image.open(path)
        long_side = max(image.size)
        if long_side < 1200:
            raise SystemExit(f"{name} is {image.size[0]}x{image.size[1]}, long side under 1200")
        print(f"contact {name} {image.size[0]}x{image.size[1]}")


def main():
    directory = Path(sys.argv[1]) if len(sys.argv) > 1 else Path("docs/source/science/figures")
    check(directory)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
