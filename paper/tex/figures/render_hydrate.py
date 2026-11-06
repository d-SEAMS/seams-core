#!/usr/bin/env python3
"""Draw the two methane-hydrate cells from their xyz records."""

import math
import sys
from pathlib import Path

from PIL import Image, ImageDraw

sys.path.insert(0, str(Path(__file__).resolve().parent))
from render_finish import (  # noqa: E402
    BOND_A,
    RUHI,
    draw_structure,
    font,
    neighbour_graph,
)

ROOT = Path(__file__).resolve().parents[3]
XYZ = ROOT / "paper" / "tex" / "figures" / "xyzrender_xyz"
FIG = ROOT / "paper" / "tex" / "figures"
CUTOFF = 3.5


def read_xyz(path):
    lines = Path(path).read_text(encoding="utf-8").splitlines()
    comment = lines[1]
    raw = comment.split('Lattice="', 1)[1].split('"', 1)[0].split()
    values = [float(item) for item in raw[:9]]
    lengths = (values[0], values[4], values[8])
    box = [(0.0, lengths[0]), (0.0, lengths[1]), (0.0, lengths[2])]
    points = []
    species = []
    for row in lines[2:]:
        if not row.strip():
            continue
        parts = row.split()
        species.append(parts[0])
        points.append((float(parts[-3]), float(parts[-2]), float(parts[-1])))
    return points, species, box


def main():
    image = Image.new("RGB", (1280, 640), RUHI["white"])
    draw = ImageDraw.Draw(image)
    fonts = (font(22),)
    panels = ("si_hydrate.xyz", "sii_hydrate.xyz")
    for index, name in enumerate(panels):
        points, species, box = read_xyz(XYZ / name)
        _graph, edges = neighbour_graph(points, box, CUTOFF)
        label = "(a) structure I" if index == 0 else "(b) structure II"

        def edge_colour(_a, _b, colour=RUHI["bond"]):
            return colour

        def atom_colour(atom, kinds=species):
            return RUHI["teal"] if kinds[atom] == "C" else RUHI["oxygen"]

        draw_structure(
            draw,
            (index * 640, 0),
            (640, 640),
            points,
            box,
            edges,
            edge_colour,
            atom_colour,
            label,
            fonts,
        )
        out = FIG / name.replace(".xyz", "_xyzrender.png")
        # Each cell is also written on its own, which is what the paper includes.
        panel = Image.new("RGB", (640, 640), RUHI["white"])
        panel_draw = ImageDraw.Draw(panel)
        draw_structure(
            panel_draw,
            (0, 0),
            (640, 640),
            points,
            box,
            edges,
            edge_colour,
            atom_colour,
            label,
            fonts,
        )
        panel.save(out, "PNG")
        print(f"wrote {out.relative_to(ROOT)}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
