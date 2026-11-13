#!/usr/bin/env python3
"""Figure of the 230 K ice/brine front, at CPC single-column width.

Panel (a) is two OVITO renders of the 20-pair cell: the frame at the
start of production, after the liquid half has been melted and quenched,
and the frame at 5 ns. Water oxygens are coloured by the d-SEAMS CHILL+
label. Na and Cl are drawn as larger spheres of their own colours.
Hydrogens are omitted so the ions remain visible in the liquid. The
legend sits under the renders.

Panel (b) is read from replicas.csv: at 5 ns, the hexagonal and cubic
counts and the ice thickness against the NaCl loading. Each replica is
one marker. A bar spans the replicas when at least two have the frame.
A loading with no frames in the dataset is marked, not plotted.

The CHILL+ colours are the palette of scripts/science_render.py. The
figure is 85 mm wide. OVITO renders the panels at 300 ppi and matplotlib
lays out the PDF.

    python3 paper/data/brine-frames/render_front.py \
        --xyz d-SEAMS/brine-front-230K/data/frames/T230_m20_r0_t0.xyz \
              d-SEAMS/brine-front-230K/data/frames/T230_m20_r0_tend.xyz \
        --replicas paper/data/brine-frames/replicas.csv \
        --pdf paper/tex/figures/brine_front.pdf
"""

from __future__ import annotations

import argparse
import csv
import os
import sys
from pathlib import Path

import numpy as np

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))
import regenerate as deck  # noqa: E402

# scripts/science_render.py, published CHILL+ classes.
CUBIC = (0.000, 0.447, 0.698)
HEXAGONAL = (0.902, 0.624, 0.000)
INTERFACIAL = (0.000, 0.620, 0.451)
LIQUID = (0.780, 0.780, 0.780)
CLATHRATE = (0.337, 0.706, 0.914)
INTER_CLATHRATE = (0.800, 0.475, 0.655)
COLOUR = {
    "cubic": CUBIC,
    "hexagonal": HEXAGONAL,
    "water": LIQUID,
    "interfacial": INTERFACIAL,
    "clathrate": CLATHRATE,
    "interClathrate": INTER_CLATHRATE,
    "unclassified": (0.45, 0.45, 0.45),
    "reCubic": CUBIC,
    "reHex": HEXAGONAL,
}
NA = (0.37, 0.24, 0.60)
CL = (0.12, 0.12, 0.12)
LEGEND = (
    ("hexagonal", HEXAGONAL),
    ("cubic", CUBIC),
    ("interfacial", INTERFACIAL),
    ("liquid", LIQUID),
    ("clathrate", CLATHRATE),
    ("interf. clathrate", INTER_CLATHRATE),
    ("Na$^+$", NA),
    ("Cl$^-$", CL),
)
RADIUS = {"O": 1.05, "Na": 1.55, "Cl": 1.75}
PPI = 300
WIDTH_IN = 85.0 / 25.4
FONT_PT = 8.0
MARKERS = {"r0": "o", "r1": "s"}
PAIRS = (0, 10, 20, 40)


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
    os.environ.setdefault("QT_QPA_PLATFORM", "ovitoheadless")
    # OVITO needs a Vulkan driver even for Tachyon; Mesa lavapipe is the software one.
    lvp = "/usr/share/vulkan/icd.d/lvp_icd.json"
    if Path(lvp).is_file():
        os.environ.setdefault("VK_DRIVER_FILES", lvp)
        os.environ.setdefault("VK_ICD_FILENAMES", lvp)
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


def read_replicas(path):
    with open(path) as handle:
        rows = list(csv.DictReader(line for line in handle if not line.startswith("#")))
    return [row for row in rows if row["frame"] == "tend"]


def panel_b(fig, rect_left, rect_right, legend_y, rows):
    import matplotlib.pyplot as plt

    left = fig.add_axes(rect_left)
    right = fig.add_axes(rect_right)
    have = {int(r["pairs"]) for r in rows if r["status"] == "ok"}
    replicas = sorted({r["replica"] for r in rows if r["status"] == "ok"})
    absent = sorted({r["replica"] for r in rows} - set(replicas))
    offsets = {name: (i - (len(replicas) - 1) / 2.0) * 1.4 for i, name in enumerate(replicas)}

    def series(ax, column, colour, label):
        for pairs in PAIRS:
            values = [
                (r["replica"], float(r[column]))
                for r in rows
                if int(r["pairs"]) == pairs and r["status"] == "ok" and r[column] != ""
            ]
            if len(values) >= 2:
                lo = min(v for _n, v in values)
                hi = max(v for _n, v in values)
                ax.plot([pairs, pairs], [lo, hi], color=colour, lw=2.4, alpha=0.45, solid_capstyle="butt")
            for name, value in values:
                ax.plot(
                    pairs + offsets[name], value, MARKERS.get(name, "D"), ms=4.5,
                    mfc=colour, mec="black", mew=0.5, ls="none",
                )
        if label:
            handles.append(plt.Line2D([], [], marker="s", ms=5, mfc=colour, mec="black", mew=0.5, ls="none",
                                      label=label))

    handles = []
    series(left, "hexagonal", HEXAGONAL, "hexagonal")
    series(left, "cubic", CUBIC, "cubic")
    series(right, "ice_thickness_A", (0.35, 0.35, 0.35), None)
    left.set_ylabel("molecules")
    right.set_ylabel("ice thickness (Å)")
    for ax in (left, right):
        ax.set_xticks(PAIRS)
        ax.set_xlim(-5, 45)
        ax.set_xlabel("NaCl pairs")
        ax.tick_params(length=2.5, pad=1.5)
        for pairs in PAIRS:
            if pairs not in have:
                ax.axvspan(pairs - 3, pairs + 3, color=(0.92, 0.92, 0.92), lw=0)
                ax.text(pairs, 0.5, "not in dataset", transform=ax.get_xaxis_transform(),
                        rotation=90, ha="center", va="center", fontsize=FONT_PT - 1.5, color="0.35")
    _lo, hi = left.get_ylim()
    left.set_ylim(-0.06 * hi, hi * 1.08)
    lo, hi = right.get_ylim()
    right.set_ylim(lo - 0.08 * (hi - lo), hi + 0.08 * (hi - lo))
    for name in replicas:
        handles.append(plt.Line2D([], [], marker=MARKERS.get(name, "D"), ms=4.5, mfc="white", mec="black",
                                  mew=0.6, ls="none", label=name))
    handles.append(plt.Line2D([], [], color="0.4", lw=2.4, alpha=0.45, label="range"))
    fig.legend(handles=handles, loc="lower center", ncol=len(handles), frameon=False, handletextpad=0.25,
               columnspacing=0.9, bbox_to_anchor=(0.5, legend_y))
    for name in absent:
        print(f"panel b: replica {name} has no 5 ns frame in the dataset", flush=True)
    for pairs in PAIRS:
        if pairs not in have:
            print(f"panel b: {pairs} NaCl pairs has no 5 ns frame in the dataset", flush=True)


def figure(pdf, panels, rows):
    import matplotlib

    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    from matplotlib.patches import Patch

    plt.rcParams.update({
        "font.size": FONT_PT,
        "axes.labelsize": FONT_PT,
        "xtick.labelsize": FONT_PT,
        "ytick.labelsize": FONT_PT,
        "legend.fontsize": FONT_PT,
        "axes.linewidth": 0.6,
        "pdf.fonttype": 42,
    })
    (img_a, step_a), (img_b, step_b) = panels
    aspect = img_a.shape[0] / img_a.shape[1]
    gap_in = 0.06
    render_w = (WIDTH_IN - gap_in) / 2.0
    render_h = render_w * aspect
    title_h = 0.20
    legend_h = 0.46
    b_legend_h = 0.30
    b_h = 1.45
    b_pad = 0.36
    height = title_h + render_h + legend_h + b_legend_h + b_h + b_pad
    fig = plt.figure(figsize=(WIDTH_IN, height))

    def frac(x, y, w, h):
        return (x / WIDTH_IN, y / height, w / WIDTH_IN, h / height)

    y_render = height - title_h - render_h
    for i, (image, step) in enumerate(((img_a, step_a), (img_b, step_b))):
        ax = fig.add_axes(frac(i * (render_w + gap_in), y_render, render_w, render_h))
        ax.imshow(image, interpolation="lanczos")
        ax.set_axis_off()
        ax.set_title(f"t = {step * 2e-6:.0f} ns" if step else "t = 0", pad=2)
    fig.text(0.0, 1.0, "(a)", ha="left", va="top", fontweight="bold", fontsize=FONT_PT + 1)
    handles = [Patch(facecolor=c, edgecolor="black", linewidth=0.4, label=n) for n, c in LEGEND]
    fig.legend(
        handles=handles, loc="upper center", ncol=4, frameon=False, handlelength=1.1, handleheight=1.1,
        columnspacing=0.8, handletextpad=0.35, bbox_to_anchor=(0.5, (y_render - 0.02) / height),
    )
    y_b = b_pad
    left_margin, mid_gap, right_margin = 0.46, 0.50, 0.04
    panel_w = (WIDTH_IN - left_margin - mid_gap - right_margin) / 2.0
    panel_b(
        fig,
        frac(left_margin, y_b, panel_w, b_h),
        frac(left_margin + panel_w + mid_gap, y_b, panel_w, b_h),
        (y_b + b_h + 0.02) / height,
        rows,
    )
    fig.text(0.0, (y_b + b_h + b_legend_h) / height, "(b)", ha="left", va="top", fontweight="bold",
             fontsize=FONT_PT + 1)
    fig.savefig(pdf, dpi=PPI)
    plt.close(fig)
    return pdf.stat().st_size


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("--dump", type=Path, nargs="+")
    parser.add_argument("--xyz", type=Path, nargs="+")
    parser.add_argument("--lib", default="")
    parser.add_argument("--replicas", type=Path, required=True)
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
        panels_in = [start, end]
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
    rows = read_replicas(args.replicas)
    if not any(row["status"] == "ok" for row in rows):
        raise SystemExit(f"{args.replicas} has no 5 ns frame")
    length = panels_in[0][2]
    width = int(round((WIDTH_IN / 2.0) * PPI))
    height = int(round(width * length[2] / length[0]))
    out = args.pdf.parent
    out.mkdir(parents=True, exist_ok=True)
    from matplotlib.image import imread

    panels = []
    for name, (step, origin, box, atoms, labels) in zip(("a", "b"), panels_in):
        pos, col, rad = select(atoms, labels)
        path = out / f".brine-panel-{name}.png"
        render_panel(path, origin, box, pos, col, rad, (width, height))
        image = imread(path)[..., :3]
        panels.append(((image * 255).round().astype(np.uint8), step))
        path.unlink(missing_ok=True)
        print(f"panel {name} step {step}", flush=True)
    nbytes = figure(args.pdf, panels, rows)
    print(f"wrote {args.pdf} bytes {nbytes}")
    if nbytes > 1_000_000:
        raise SystemExit(f"PDF is {nbytes} bytes, over 1 MB")


if __name__ == "__main__":
    main()
