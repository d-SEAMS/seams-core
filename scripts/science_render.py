#!/usr/bin/env python3
# Copyright (c) d-SEAMS developers

"""Offscreen OVITO renders of the basal slab and the nucleus.

Tachyon ray-traces the oxygens. CHILL+ colours follow the published
Nguyen-Molinero classes written by science-front (the chill column).
Liquid molecules on the nucleus views are faint, and oxygen-oxygen
contacts inside the largest ice cluster are drawn as bonds.

Writes PNG at the requested pixel size and a PDF that embeds that image
at 300 dpi. Neither file belongs in the repository.
"""

from __future__ import annotations

import argparse
import math
import os
import subprocess
import sys
from pathlib import Path

import numpy as np

# Published CHILL+ classes, molSys::atom_state_type.
CHILL = {
    0: ("cubic", (0.000, 0.447, 0.698)),
    1: ("hexagonal", (0.902, 0.624, 0.000)),
    3: ("interfacial", (0.000, 0.620, 0.451)),
    2: ("liquid", (0.780, 0.780, 0.780)),
    4: ("clathrate", (0.337, 0.706, 0.914)),
    5: ("interfacial clathrate", (0.800, 0.475, 0.655)),
}
ICE = {0, 1, 3, 4, 5}
ROOT = Path(__file__).resolve().parents[1]


def ensure_runtime() -> None:
    os.environ["QT_QPA_PLATFORM"] = "ovitoheadless"
    os.environ.setdefault("LIBGL_ALWAYS_SOFTWARE", "1")
    os.environ.setdefault("OVITO_THREAD_COUNT", "2")
    os.environ.pop("VK_LOADER_DEBUG", None)
    os.environ.pop("QT_DEBUG_PLUGINS", None)
    lvp = "/usr/share/vulkan/icd.d/lvp_icd.json"
    if Path(lvp).is_file():
        os.environ.setdefault("VK_DRIVER_FILES", lvp)
        os.environ.setdefault("VK_ICD_FILENAMES", lvp)
    # science-front needs the toolchain libm (GLIBC_2.43). Qt needs libGL
    # from the system directory. The toolchain directory comes first so the
    # newer libm wins, and the system directory still supplies libGL.
    wanted = [
        "/cvmfs/software.eessi.io/versions/2026.06/compat/linux/x86_64/lib64",
        "/lib/x86_64-linux-gnu",
    ]
    wanted = [path for path in wanted if Path(path).is_dir()]
    current = [path for path in os.environ.get("LD_LIBRARY_PATH", "").split(":") if path]
    merged = wanted + [path for path in current if path not in wanted]
    if merged != current:
        os.environ["LD_LIBRARY_PATH"] = ":".join(merged)
        os.execv(sys.executable, [sys.executable, *sys.argv])


def read_table(path: Path):
    header = None
    box = None
    rows = []
    for line in path.read_text().splitlines():
        if not line:
            continue
        if line.startswith("# box"):
            box = np.array([float(v) for v in line.split()[2:5]])
            continue
        if line.startswith("#"):
            continue
        parts = line.split()
        if header is None:
            header = parts
            continue
        rows.append(parts)
    return header, rows, box


def column(header, rows, name: str) -> np.ndarray:
    index = header.index(name)
    return np.array([float(row[index]) for row in rows])


def read_dump_frame(path: Path, frame: int):
    """Return box and oxygen coordinates for a 1-based frame index."""
    if frame < 1:
        raise SystemExit(f"frame {frame} is not 1-based")
    seen = 0
    box = None
    coords = []
    with path.open() as handle:
        while True:
            line = handle.readline()
            if not line:
                break
            if not line.startswith("ITEM: TIMESTEP"):
                continue
            handle.readline()
            seen += 1
            if handle.readline().strip() != "ITEM: NUMBER OF ATOMS":
                raise SystemExit("dump is missing the atom count")
            n = int(handle.readline())
            bounds = handle.readline()
            if not bounds.startswith("ITEM: BOX BOUNDS"):
                raise SystemExit("dump is missing box bounds")
            lengths = []
            for _ in range(3):
                lo, hi = handle.readline().split()[:2]
                lengths.append(float(hi) - float(lo))
            atoms = handle.readline()
            if not atoms.startswith("ITEM: ATOMS"):
                raise SystemExit("dump is missing the atom table")
            names = atoms.split()[2:]
            xi, yi, zi = names.index("x"), names.index("y"), names.index("z")
            block = []
            for _ in range(n):
                parts = handle.readline().split()
                block.append((float(parts[xi]), float(parts[yi]), float(parts[zi])))
            if seen == frame:
                return np.array(lengths), np.array(block)
            coords.append(block)
            box = lengths
    raise SystemExit(f"{path} has {seen} frames, not {frame}")


def write_one_frame(path: Path, box: np.ndarray, pos: np.ndarray) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w") as handle:
        handle.write("ITEM: TIMESTEP\n0\nITEM: NUMBER OF ATOMS\n")
        handle.write(f"{pos.shape[0]}\nITEM: BOX BOUNDS pp pp pp\n")
        for length in box:
            handle.write(f"0.0 {length:.8f}\n")
        handle.write("ITEM: ATOMS id type x y z\n")
        for i, (x, y, z) in enumerate(pos, start=1):
            handle.write(f"{i} 1 {x:.6f} {y:.6f} {z:.6f}\n")


def run_front(dump: Path, out: Path, command: str) -> None:
    out.mkdir(parents=True, exist_ok=True)
    binary = os.environ.get("SCIENCE_FRONT", str(ROOT / "bbdir" / "science-front"))
    threads = os.environ.get("SCIENCE_THREADS", "1")
    env = os.environ.copy()
    env["OMP_NUM_THREADS"] = threads
    cmd = [binary, command, str(dump), "--out-dir", str(out), "--threads", threads]
    if command == "interface":
        cmd[3:3] = ["--ps", "1", "--bin", "1"]
    subprocess.run(cmd, check=True, env=env)


def label_positions(box: np.ndarray, pos: np.ndarray, work: Path, command: str):
    work.mkdir(parents=True, exist_ok=True)
    dump = work / "one.lammpstrj"
    write_one_frame(dump, box, pos)
    out = work / "labels"
    run_front(dump, out, command)
    name = "nucleus-atoms.csv" if command == "nucleation" else "interface-atoms.csv"
    header, rows, _ = read_table(out / name)
    chill = column(header, rows, "chill").astype(int)
    extended = None
    if "in_extended" in header:
        extended = column(header, rows, "in_extended") > 0.5
    return chill, extended


def load_labeled(path: Path):
    header, rows, box = read_table(path)
    pos = np.column_stack([column(header, rows, name) for name in ("x", "y", "z")])
    chill = column(header, rows, "chill").astype(int)
    extended = None
    if "in_extended" in header:
        extended = column(header, rows, "in_extended") > 0.5
    if box is None:
        raise SystemExit(f"{path} has no box comment")
    return box, pos, chill, extended


def minimum_image_pairs(pos: np.ndarray, box: np.ndarray, mask: np.ndarray, cutoff: float):
    chosen = np.flatnonzero(mask)
    if chosen.size < 2:
        return np.zeros((0, 2), dtype=np.int32)
    pts = pos[chosen]
    pairs = []
    cut2 = cutoff * cutoff
    # A few thousand ice molecules: a cell list keeps this off the quadratic path.
    cell = cutoff
    shape = np.maximum(1, np.floor(box / cell).astype(int))
    bins = {}
    keys = np.floor(pts / box * shape).astype(int) % shape
    for index, key in enumerate(keys):
        bins.setdefault(tuple(key), []).append(index)
    for index, key in enumerate(keys):
        for dx in (-1, 0, 1):
            for dy in (-1, 0, 1):
                for dz in (-1, 0, 1):
                    neighbour = tuple((key + (dx, dy, dz)) % shape)
                    for other in bins.get(neighbour, []):
                        if other <= index:
                            continue
                        delta = pts[other] - pts[index]
                        delta -= box * np.round(delta / box)
                        if float(np.dot(delta, delta)) <= cut2:
                            pairs.append((int(chosen[index]), int(chosen[other])))
    if not pairs:
        return np.zeros((0, 2), dtype=np.int32)
    return np.array(pairs, dtype=np.int32)


def unwrap_cluster(pos: np.ndarray, box: np.ndarray, mask: np.ndarray, cutoff: float):
    """Place the ice cluster in one image and park the liquid around it.

    OVITO draws a bond as a straight segment. A periodic neighbour pair
    otherwise becomes a line across the cell.
    """
    placed = pos.copy()
    index = np.flatnonzero(mask)
    if index.size == 0:
        return placed, box.copy()
    pairs = minimum_image_pairs(pos, box, mask, cutoff)
    neighbours: dict[int, list[int]] = {}
    for a, b in pairs:
        neighbours.setdefault(int(a), []).append(int(b))
        neighbours.setdefault(int(b), []).append(int(a))
    seed = int(index[0])
    placed[seed] = pos[seed]
    seen = {seed}
    stack = [seed]
    while stack:
        current = stack.pop()
        for other in neighbours.get(current, []):
            if other in seen:
                continue
            delta = pos[other] - pos[current]
            delta -= box * np.round(delta / box)
            placed[other] = placed[current] + delta
            seen.add(other)
            stack.append(other)
    for atom in index:
        if int(atom) not in seen:
            placed[int(atom)] = pos[int(atom)]
            seen.add(int(atom))
    centre = placed[index].mean(axis=0)
    placed[index] -= centre
    origin = pos[seed]
    anchor = placed[seed]
    for atom in range(pos.shape[0]):
        if atom in seen:
            continue
        delta = pos[atom] - origin
        delta -= box * np.round(delta / box)
        placed[atom] = anchor + delta
    placed -= placed.min(axis=0) - 4.0
    cell = placed.max(axis=0) + 4.0
    return placed, cell


def scene(box, pos, chill, bond_mask, liquid_alpha: float, ice_radius: float, liquid_radius: float, cutoff: float):
    from ovito.data import DataCollection, Particles, SimulationCell
    from ovito.pipeline import Pipeline, StaticSource
    from ovito.vis import BondsVis

    data = DataCollection()
    periodic = bond_mask is None
    cell = SimulationCell(pbc=(periodic, periodic, periodic))
    cell[:, :3] = np.diag(box.astype(float))
    data.objects.append(cell)
    particles = Particles()
    particles.create_property("Position", data=np.ascontiguousarray(pos, dtype=np.float64))
    colours = np.array([CHILL.get(int(code), ("other", (0.5, 0.5, 0.5)))[1] for code in chill])
    particles.create_property("Color", data=colours)
    liquid = np.array([int(code) == 2 for code in chill])
    transparency = np.where(liquid, 1.0 - liquid_alpha, 0.0)
    particles.create_property("Transparency", data=transparency.astype(np.float64))
    radius = np.where(liquid, liquid_radius, ice_radius).astype(np.float64)
    particles.create_property("Radius", data=radius)
    if bond_mask is not None and np.any(bond_mask):
        topology = minimum_image_pairs(pos, box, bond_mask, cutoff)
        if topology.shape[0]:
            bonds = particles.create_bonds()
            bonds.create_property("Topology", data=topology)
            bonds.vis.width = 0.22
            bonds.vis.color = (0.15, 0.15, 0.15)
            bonds.vis.coloring_mode = BondsVis.ColoringMode.Uniform
    data.objects.append(particles)
    particles.vis.radius = ice_radius
    return Pipeline(source=StaticSource(data=data))


def frame_view(pipeline, box, pos, focus_mask, png: Path, width: int, height: int, perspective: bool):
    from ovito.vis import TachyonRenderer, Viewport

    if focus_mask is not None and np.any(focus_mask):
        centre = pos[focus_mask].mean(axis=0)
        extent = pos[focus_mask].max(axis=0) - pos[focus_mask].min(axis=0)
        span = float(np.max(extent))
    else:
        centre = box * 0.5
        span = float(np.max(box))
    distance = max(span * 1.85, 18.0)
    view = Viewport(type=Viewport.Type.Perspective if perspective else Viewport.Type.Ortho)
    pipeline.add_to_scene()
    try:
        if perspective:
            view.camera_pos = centre + np.array([distance * 0.55, -distance, distance * 0.32])
            view.camera_dir = centre - np.array(view.camera_pos)
            view.camera_up = (0.0, 0.0, 1.0)
            view.fov = math.radians(30.0)
        else:
            view.camera_dir = (0.0, -1.0, 0.0)
            view.camera_up = (0.0, 0.0, 1.0)
            view.zoom_all(size=(width, height))
        renderer = TachyonRenderer(shadows=True, antialiasing_samples=8, direct_light_intensity=0.95)
        png.parent.mkdir(parents=True, exist_ok=True)
        view.render_image(
            filename=str(png),
            size=(width, height),
            renderer=renderer,
            background=(1.0, 1.0, 1.0),
            frame=0,
        )
    finally:
        pipeline.remove_from_scene()
    return view


def write_pdf(png: Path, pdf: Path, caption: str) -> None:
    import matplotlib

    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    from PIL import Image

    image = np.asarray(Image.open(png))
    height, width = image.shape[:2]
    dpi = 300
    legend_in = 0.42
    fig = plt.figure(figsize=(width / dpi, height / dpi + legend_in), dpi=dpi)
    ax = fig.add_axes((0.0, legend_in / (height / dpi + legend_in), 1.0, 1.0 - legend_in / (height / dpi + legend_in)))
    ax.imshow(image, interpolation="nearest")
    ax.set_axis_off()
    handles = [
        plt.Line2D([0], [0], marker="o", color="none", markerfacecolor=colour, markersize=5, label=name)
        for _, (name, colour) in CHILL.items()
    ]
    fig.legend(handles=handles, loc="lower center", ncol=3, frameon=False, fontsize=7)
    fig.text(0.01, 0.012, caption, fontsize=7, color="#222222", ha="left", va="bottom")
    pdf.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(pdf, dpi=dpi)
    plt.close(fig)


def render_labeled(box, pos, chill, extended, png: Path, pdf: Path, caption: str, kind: str, pixels: int) -> None:
    if kind == "nucleus":
        bond_mask = extended if extended is not None else np.isin(chill, list(ICE))
        pos, box = unwrap_cluster(pos, box, bond_mask, 3.5)
        focus = bond_mask
        pipeline = scene(box, pos, chill, bond_mask, 0.18, 1.05, 0.48, 3.5)
        perspective = True
        height = int(round(pixels * 0.8))
        width = pixels
    else:
        bond_mask = None
        focus = None
        pipeline = scene(box, pos, chill, bond_mask, 0.9, 1.05, 0.72, 3.5)
        perspective = False
        # Looking along y: the picture is x against the basal normal z.
        aspect = float(box[2] / max(box[0], 1e-6))
        width = pixels
        height = max(int(round(pixels * aspect)), pixels)
    frame_view(pipeline, box, pos, focus, png, width, height, perspective)
    write_pdf(png, pdf, caption)


def main() -> int:
    ensure_runtime()
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--out", type=Path, required=True)
    parser.add_argument("--pixels", type=int, default=2400)
    parser.add_argument("--tag", default="", help="filename stem, for example 230K")
    parser.add_argument("--interface-dump", type=Path)
    parser.add_argument("--interface-frames", default="1,-1", help="1-based indices; -1 is the last frame")
    parser.add_argument("--nucleus-dump", type=Path)
    parser.add_argument("--nucleus-min-frame", type=int, default=2)
    parser.add_argument("--nucleus-late-frame", type=int, default=-1)
    parser.add_argument("--work", type=Path, default=ROOT / "data" / "science-frames" / "render-work")
    args = parser.parse_args()
    if args.pixels < 900:
        raise SystemExit("publication renders need at least 900 px on the long edge (300 dpi at 3 in)")
    args.out.mkdir(parents=True, exist_ok=True)
    stem = f"{args.tag}-" if args.tag else ""

    if args.interface_dump:
        nframes = count_frames(args.interface_dump)
        tokens = [token for token in args.interface_frames.split(",") if token]
        if len(tokens) == 1 and tokens[0].startswith("-"):
            names = ("end",)
        elif len(tokens) == 1:
            names = ("start",)
        else:
            names = ("start", "end")
        for token, tag in zip(tokens, names):
            index = int(token)
            if index < 0:
                index = nframes + 1 + index
            box, pos = read_dump_frame(args.interface_dump, index)
            chill, _ = label_positions(box, pos, args.work / f"interface-{tag}", "interface")
            render_labeled(
                box,
                pos,
                chill,
                None,
                args.out / f"interface-slab-{stem}{tag}.png",
                args.out / f"interface-slab-{stem}{tag}.pdf",
                f"TIP4P/Ice basal slab, {tag} frame, CHILL+",
                "slab",
                args.pixels,
            )
            print(f"rendered interface {tag} frame {index}", flush=True)

    if args.nucleus_dump:
        nframes = count_frames(args.nucleus_dump)
        late = args.nucleus_late_frame
        if late < 0:
            late = nframes + 1 + late
        for index, tag in ((args.nucleus_min_frame, "minimum"), (late, "late")):
            box, pos = read_dump_frame(args.nucleus_dump, index)
            chill, extended = label_positions(box, pos, args.work / f"nucleus-{tag}", "nucleation")
            if extended is None:
                extended = np.isin(chill, list(ICE))
            render_labeled(
                box,
                pos,
                chill,
                extended,
                args.out / f"nucleus-{stem}{tag}.png",
                args.out / f"nucleus-{stem}{tag}.pdf",
                f"mW nucleus, {tag} frame, CHILL+ with liquid faint",
                "nucleus",
                args.pixels,
            )
            print(f"rendered nucleus {tag} frame {index}", flush=True)
    return 0


def count_frames(path: Path) -> int:
    count = 0
    with path.open() as handle:
        for line in handle:
            if line.startswith("ITEM: TIMESTEP"):
                count += 1
    return count


if __name__ == "__main__":
    raise SystemExit(main())
