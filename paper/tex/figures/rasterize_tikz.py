"""Rasterize the six TikZ figures to PNG.

Run from the repository root, with tectonic and pdftoppm on PATH:

    python3 paper/tex/figures/rasterize_tikz.py

The script writes a temporary wrapper next to the manuscript and deletes it.
"""

import os
import shutil
import subprocess
import sys
from pathlib import Path

from PIL import Image

NAMES = (
    "fig-cages",
    "fig-update",
    "fig-versions",
    "fig-codes",
    "fig-churn",
    "fig-nucleation",
)
TEX = Path("paper/tex")
WRAPPER = r"""\documentclass{article}
\usepackage{amsmath,amssymb}
\usepackage{xcolor}
\usepackage{siunitx}
\sisetup{exponent-product=\times, group-minimum-digits=4}
\usepackage{tikz}
\usetikzlibrary{arrows.meta,calc}
\usepackage{pgfplots}
\usepgfplotslibrary{groupplots}
\pgfplotsset{compat=1.17}
\pagestyle{empty}
\begin{document}
\input{figures/style}
\input{figures/data/macros}
\input{figures/%s}
\end{document}
"""


def trim(path, dest):
    image = Image.open(path).convert("RGB")
    background = (255, 255, 255)
    pixels = image.load()
    width, height = image.size

    def row_blank(y):
        return all(pixels[x, y] == background for x in range(width))

    def col_blank(x):
        return all(pixels[x, y] == background for y in range(height))

    top = 0
    while top < height and row_blank(top):
        top += 1
    bottom = height - 1
    while bottom > top and row_blank(bottom):
        bottom -= 1
    left = 0
    while left < width and col_blank(left):
        left += 1
    right = width - 1
    while right > left and col_blank(right):
        right -= 1
    pad = 12
    box = (
        max(0, left - pad),
        max(0, top - pad),
        min(width, right + 1 + pad),
        min(height, bottom + 1 + pad),
    )
    image.crop(box).save(dest)


def main():
    tectonic = shutil.which("tectonic")
    pdftoppm = shutil.which("pdftoppm")
    if tectonic is None or pdftoppm is None:
        print("tectonic and pdftoppm must be on PATH", file=sys.stderr)
        return 1
    os.chdir(TEX)
    wrap = Path("_wrap.tex")
    code = 0
    try:
        for name in NAMES:
            wrap.write_text(WRAPPER % name, encoding="utf-8")
            proc = subprocess.run([tectonic, "-X", "compile", str(wrap)], check=False)
            if proc.returncode != 0:
                code = proc.returncode
                break
            subprocess.run(
                [pdftoppm, "-png", "-r", "160", "_wrap.pdf", name],
                check=True,
            )
            trim(Path(f"{name}-1.png"), Path("figures") / f"{name}.png")
            Path(f"{name}-1.png").unlink(missing_ok=True)
    finally:
        for extra in (
            "_wrap.tex",
            "_wrap.pdf",
            "_wrap.aux",
            "_wrap.log",
            "_wrap.out",
            "_wrap.fls",
        ):
            Path(extra).unlink(missing_ok=True)
    return code


if __name__ == "__main__":
    sys.exit(main())
