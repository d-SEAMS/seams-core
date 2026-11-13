#!/usr/bin/env python3
"""Render the brine front for paper Fig 18 and SI Figs 25-26.

The picture is the science-tutorial brine figure: ice seated on the
lower surface, depth in nanometres, and the shared CHILL+ palette.
``--width-mm`` is the atom-panel width at 300 ppi. The file written to
``--out`` is wider, because the depth axis and the legend sit outside
the panel.

    python docs/examples/render_brine.py \\
        --xyz docs/.cache/brine-front-230K/<rev>/T230_m20_r0_tend.xyz \\
        --out brine.png \\
        --width-mm 85
"""

from __future__ import annotations

import argparse
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT / "docs" / "examples"))

import science_tutorials as tutorials  # noqa: E402
from render_scenes import render_brine  # noqa: E402

PPI = 300.0


def panel_width_px(width_mm, ppi=PPI):
    if width_mm <= 0 or ppi <= 0:
        raise SystemExit("--width-mm and --ppi must be positive")
    return int(round(float(width_mm) / 25.4 * float(ppi)))


def main(argv=None):
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("--xyz", type=Path, required=True, help="labelled brine xyz frame")
    parser.add_argument("--out", type=Path, required=True, help="png path to write")
    parser.add_argument("--width-mm", type=float, required=True, help="atom panel width in millimetres")
    parser.add_argument("--ppi", type=float, default=PPI, help="pixels per inch for --width-mm")
    args = parser.parse_args(argv)
    width = panel_width_px(args.width_mm, args.ppi)
    render = tutorials._load_science_render(ROOT)
    _box, rows = tutorials.read_xyz(args.xyz)
    args.out.parent.mkdir(parents=True, exist_ok=True)
    render_brine(render, rows, args.out, _box, width=width)
    print(f"wrote {args.out} panel {width} px")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
