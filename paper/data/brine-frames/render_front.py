#!/usr/bin/env python3
"""Two-panel OVITO figure of a 230 K ice/brine front.

Panel (a) is the frame at the start of production, after the liquid half
has been melted and quenched. Panel (b) is the frame at 5 ns. Water
oxygens are coloured by the d-SEAMS CHILL+ label. Na and Cl are drawn as
larger spheres of their own colours. Hydrogens are omitted so the ions
remain visible in the liquid.

The PDF is 300 ppi. OVITO renders the panels; Pillow packs them.

    python3 paper/data/brine-frames/render_front.py \
        --xyz paper/data/brine-frames/out/T230_m20_r0_t0.xyz \
              paper/data/brine-frames/out/T230_m20_r0_tend.xyz \
        --pdf paper/tex/figures/brine_front.pdf
"""

from __future__ import annotations

import argparse
import os
import sys
from pathlib import Path

import numpy as np

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))
import regenerate as deck  # noqa: E402

# Published CHILL+ classes, in the d-SEAMS label order.
COLOUR = {
    "cubic": (0.90, 0.49, 0.13),
    "hexagonal": (0.00, 0.35, 0.64),
    "water": (0.62, 0.70, 0.76),
    "interfacial": (0.00, 0.55, 0.45),
    "clathrate": (0.55, 0.34, 0.64),
    "interClathrate": (0.35, 0.62, 0.80),
    "unclassified": (0.45, 0.45, 0.45),
    "reCubic": (0.90, 0.49, 0.13),
    "reHex": (0.00, 0.35, 0.64),
}
NA = (0.96, 0.75, 0.10)
CL = (0.70, 0.13, 0.18)
RADIUS = {"O": 1.05, "Na": 1.55, "Cl": 1.75}
# Two panels, each matching the cell aspect, at 300 ppi.
PPI = 300
PANEL_H_IN = 4.6


def select(atoms, labels):
    pos, col, rad = [], [], []
    for i, (aid, typ, x, y, z) in enumerate(atoms):
        if typ == 2:
            continue
        if typ == 1:
            name = deck.CHILL[labels[i]]
            colour = COLOUR.get(name, COLOUR["unclassified"])
            radius = RADIUS["O"]
        elif typ == 3:
            colour, radius = NA, RADIUS["Na"]
        elif typ == 4:
            colour, radius = CL, RADIUS["Cl"]
        else:
            continue
        pos.append((x, y, z))
        col.append(colour)
        rad.append(radius)
    return (
        np.asarray(pos, dtype=np.float64),
        np.asarray(col, dtype=np.float32),
        np.asarray(rad, dtype=np.float32),
    )


def render_panel(path, origin, length, pos, col, rad, size):
    # Imported here so a missing OVITO names itself on the command line.
    from ovito.data import DataCollection, Particles, SimulationCell
    from ovito.pipeline import Pipeline, StaticSource
    from ovito.vis import TachyonRenderer, Viewport

    data = DataCollection()
    cell = SimulationCell(pbc=(True, True, True))
    cell[:, 0] = (length[0], 0.0, 0.0)
    cell[:, 1] = (0.0, length[1], 0.0)
    cell[:, 2] = (0.0, 0.0, length[2])
    cell[:, 3] = tuple(origin)
    parts = Particles()
    parts.create_property("Position", data=pos)
    parts.create_property("Color", data=col)
    parts.create_property("Radius", data=rad)
    data.objects.append(cell)
    data.objects.append(parts)
    if getattr(data.cell, "vis", None) is not None:
        data.cell.vis.enabled = False
    # A fresh scene, so a second panel does not keep the first pipeline.
    import ovito

    for existing in list(ovito.scene.pipelines):
        existing.remove_from_scene()
    Pipeline(source=StaticSource(data=data)).add_to_scene()
    view = Viewport(type=Viewport.Type.Ortho)
    view.camera_pos = (
        origin[0] + 0.5 * length[0],
        origin[1] + length[1] + 40.0,
        origin[2] + 0.5 * length[2],
    )
    view.camera_dir = (0.0, -1.0, 0.0)
    view.camera_up = (0.0, 0.0, 1.0)
    # Orthographic fov is the half-height of the view, in angstrom.
    # Half the cell height, plus a small margin, puts the slab on the edges.
    view.fov = length[2] * 0.52
    view.render_image(filename=str(path), size=size, renderer=TachyonRenderer(), background=(1, 1, 1))


def compose(panels, pdf, size):
    from PIL import Image, ImageDraw, ImageFont

    width, height = size
    gap = 24
    canvas = Image.new("RGB", (width * 2 + gap, height), (255, 255, 255))
    tags = ("(a)", "(b)")
    try:
        font = ImageFont.truetype("DejaVuSans-Bold.ttf", 42)
    except OSError:
        font = ImageFont.load_default()
    for i, panel in enumerate(panels):
        image = Image.open(panel).convert("RGB")
        if image.size != size:
            image = image.resize(size, Image.Resampling.LANCZOS)
        canvas.paste(image, (i * (width + gap), 0))
        draw = ImageDraw.Draw(canvas)
        draw.text((i * (width + gap) + 18, 12), tags[i], fill=(20, 20, 20), font=font)
    # 300 ppi is stored in the PDF page size: pixels / 300 = inches.
    canvas.save(pdf, "PDF", resolution=PPI, quality=90)
    return pdf.stat().st_size


def parse_xyz(path):
    """Read one labelled frame written by regenerate.write_xyz."""
    lines = path.read_text().splitlines()
    if len(lines) < 3:
        raise SystemExit(f"empty frame {path}")
    header = lines[1]
    step = int(header.split("step=")[1].split()[0])
    origin = [float(value) for value in header.split('origin="')[1].split('"')[0].split()]
    lattice = [float(value) for value in header.split('Lattice="')[1].split('"')[0].split()]
    length = [lattice[0], lattice[4], lattice[8]]
    species = {"O": 1, "H": 2, "Na": 3, "Cl": 4}
    atoms = []
    labels = []
    for index, line in enumerate(lines[2:], start=1):
        name, x, y, z, chill = line.split()
        atoms.append((index, species[name], float(x), float(y), float(z)))
        labels.append(deck.CHILL.index(chill) if chill in deck.CHILL else 6)
    return step, origin, length, atoms, labels


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--dump", type=Path, nargs="+")
    parser.add_argument("--xyz", type=Path, nargs="+")
    parser.add_argument("--lib", default="")
    parser.add_argument("--pdf", type=Path, required=True)
    parser.add_argument("--end-step", type=int, default=2_500_000)
    args = parser.parse_args()
    if args.xyz:
        labelled = [parse_xyz(path) for path in args.xyz]
        labelled.sort(key=lambda frame: frame[0])
        start = labelled[0]
        end = min(labelled, key=lambda frame: abs(frame[0] - args.end_step))
        if end[0] < args.end_step - 20_000:
            raise SystemExit(f"5 ns frame is missing; last step is {labelled[-1][0]}")
        panels_in = [(start[0], start[1], start[2], start[3], start[4]), (end[0], end[1], end[2], end[3], end[4])]
    elif args.dump:
        if not args.lib:
            raise SystemExit("pass --lib when rendering a dump")
        merged = {}
        for path in args.dump:
            for frame in deck.parse_dump(path):
                merged[frame[0]] = frame
        series = [merged[step] for step in sorted(merged)]
        if not series:
            raise SystemExit("no frames")
        start = series[0]
        end = min(series, key=lambda frame: abs(frame[0] - args.end_step))
        if end[0] < args.end_step - 20_000:
            raise SystemExit(f"5 ns frame is missing; last step is {series[-1][0]}")
        lib = deck.load_chill(args.lib)
        panels_in = []
        for frame in (start, end):
            step, origin, box, atoms = frame
            labels = deck.chill_labels(lib, atoms, origin, box)
            panels_in.append((step, origin, box, atoms, labels))
    else:
        raise SystemExit("pass --xyz or --dump")
    height = int(round(PANEL_H_IN * PPI))
    # Width follows the cell so the slab fills the panel.
    length = panels_in[0][2]
    width = int(round(height * length[0] / length[2]))
    tmp = args.pdf.parent
    tmp.mkdir(parents=True, exist_ok=True)
    panels = []
    for name, (step, origin, box, atoms, labels) in zip(("a", "b"), panels_in):
        pos, col, rad = select(atoms, labels)
        panel = tmp / f".brine-panel-{name}.png"
        render_panel(panel, origin, box, pos, col, rad, (width, height))
        panels.append(panel)
        print(f"panel {name} step {step}", flush=True)
    nbytes = compose(panels, args.pdf, (width, height))
    for panel in panels:
        panel.unlink(missing_ok=True)
    print(f"wrote {args.pdf} bytes {nbytes}")
    if nbytes > 1_000_000:
        raise SystemExit(f"PDF is {nbytes} bytes, over 1 MB")


if __name__ == "__main__":
    main()
