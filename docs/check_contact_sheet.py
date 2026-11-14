#!/usr/bin/env python3
"""Fail if a docs raster is under 1200 px on its short side."""

import os
import sys
from pathlib import Path

from PIL import Image

NAMES = ("interface.png", "brine.png", "sii.png", "sh.png")


def check(directory: Path):
    fallbacks = []
    for name in NAMES:
        path = directory / name
        if not path.is_file():
            svg = directory / name.replace(".png", ".svg")
            if svg.is_file() and os.environ.get("DOCS_CI") != "1":
                fallbacks.append(name)
                continue
            raise SystemExit(f"missing render {path}")
        image = Image.open(path)
        short_side = min(image.size)
        if short_side < 1200:
            raise SystemExit(f"{name} is {image.size[0]}x{image.size[1]}, short side under 1200")
        print(f"contact {name} {image.size[0]}x{image.size[1]}")
    if fallbacks:
        print(
            "contact sheet skipped: the ovito vulkan renderer is not available ({})".format(
                ", ".join(fallbacks)
            )
        )
        raise SystemExit(77)


def main():
    directory = Path(sys.argv[1]) if len(sys.argv) > 1 else Path("docs/source/science/figures")
    check(directory)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
