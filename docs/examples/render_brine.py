#!/usr/bin/env python3
"""Render the brine front for paper Fig 18 and SI Figs 25-26.

The picture is the science-tutorial brine figure: ice seated on the
lower surface, depth in nanometres, and the shared CHILL+ palette.
``--width-mm`` is the atom-panel width at 300 ppi. The file written to
``--out`` is wider, because the depth axis and the legend sit outside
the panel. The legend is raised until it is at least 7 pt when that
whole png is printed at ``--width-mm``.

    python docs/examples/render_brine.py \\
        --xyz docs/.cache/brine-front-230K/<rev>/T230_m20_r0_tend.xyz \\
        --out brine.png \\
        --width-mm 85
"""

from __future__ import annotations

import argparse
import math
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT / "docs" / "examples"))

import science_tutorials as tutorials  # noqa: E402
from render_scenes import render_brine  # noqa: E402

PPI = 300.0
MIN_POINT = 7.0


def panel_width_px(width_mm, ppi=PPI):
    if width_mm <= 0 or ppi <= 0:
        raise SystemExit("--width-mm and --ppi must be positive")
    return int(round(float(width_mm) / 25.4 * float(ppi)))


def printed_points(font_px, image_width_px, width_mm):
    """Point size of ``font_px`` when the whole png is printed at ``width_mm``."""
    return float(font_px) * (float(width_mm) / float(image_width_px)) * (72.0 / 25.4)


def write_figure(render, rows, box, out, width_mm, ppi=PPI):
    """Write ``out`` and return the legend point size at ``width_mm`` print width."""
    import render_scenes
    from PIL import Image

    panel = panel_width_px(width_mm, ppi)
    saved = render_scenes.LEGEND_FONT_PX
    font = saved
    out.parent.mkdir(parents=True, exist_ok=True)
    try:
        points = 0.0
        for _ in range(3):
            render_scenes.LEGEND_FONT_PX = font
            render_brine(render, rows, out, box, width=panel)
            image = Image.open(out)
            points = printed_points(font, image.width, width_mm)
            if points >= MIN_POINT:
                return points
            font = int(math.ceil(font * MIN_POINT / points))
        raise SystemExit(f"legend is {points:.2f} pt at {width_mm:g} mm, under {MIN_POINT:g} pt")
    finally:
        render_scenes.LEGEND_FONT_PX = saved


def main(argv=None):
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("--xyz", type=Path, required=True, help="labelled brine xyz frame")
    parser.add_argument("--out", type=Path, required=True, help="png path to write")
    parser.add_argument("--width-mm", type=float, required=True, help="atom panel width, and the print width used for the 7 pt check")
    parser.add_argument("--ppi", type=float, default=PPI, help="pixels per inch used to size the atom panel")
    args = parser.parse_args(argv)
    render = tutorials._load_science_render(ROOT)
    box, rows = tutorials.read_xyz(args.xyz)
    points = write_figure(render, rows, box, args.out, args.width_mm, args.ppi)
    panel = panel_width_px(args.width_mm, args.ppi)
    print(f"wrote {args.out} panel {panel} px legend {points:.2f} pt")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
