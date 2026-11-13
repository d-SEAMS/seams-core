#!/usr/bin/env python3
# Copyright (c) d-SEAMS developers

"""Offscreen OVITO renders of the basal slab and the nucleus.

Tachyon ray-traces the oxygens. CHILL+ colours follow the published
Nguyen-Molinero classes written by science-front (the chill column).
Liquid molecules on the nucleus views are faint, and oxygen-oxygen
contacts inside the largest ice cluster are drawn as bonds.

The picture is one Computer Physics Communications column (85 mm) at
600 dpi. Every atom uses the same radius and every bond the same width.
The frame includes each sphere plus a margin, so no atom is clipped and
no bond runs off the crop. Neither file belongs in the repository.
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
# Okabe-Ito. Liquid is grey so the five categorical colours stay on the ice classes.
# In the 230 K slab the molecules that look pink and sky blue inside the liquid
# are interfacial clathrate (33) and clathrate (10). CHILL+ assigns those
# patterns; the colours are not a drawing error and the cutoffs are not refit.
CHILL = {
    0: ("cubic Ic", (0.000, 0.447, 0.698)),
    1: ("hexagonal Ih", (0.902, 0.624, 0.000)),
    3: ("interfacial", (0.000, 0.620, 0.451)),
    2: ("liquid", (0.780, 0.780, 0.780)),
    4: ("clathrate", (0.337, 0.706, 0.914)),
    5: ("interfacial clathrate", (0.800, 0.475, 0.655)),
}
ICE = {0, 1, 3, 4, 5}
# One radius and one bond width for every panel. 85 mm is the CPC column.
ATOM_RADIUS = 0.90
BOND_WIDTH = 0.16
HBOND_A = 3.2
SLICE_A = 10.0
COLUMN_MM = 85.0
PRINT_DPI = 600
ROOT = Path(__file__).resolve().parents[1]


def column_width_px() -> int:
    return int(round(COLUMN_MM / 25.4 * PRINT_DPI))


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



def unit(vector: np.ndarray) -> np.ndarray:
    vector = np.asarray(vector, dtype=float)
    return vector / np.linalg.norm(vector)


def basis(forward: np.ndarray, up: np.ndarray):
    forward = unit(forward)
    right = np.cross(up, forward)
    right = unit(right)
    up = unit(np.cross(forward, right))
    return forward, right, up


def basal_forward(pos: np.ndarray, chill: np.ndarray, box: np.ndarray) -> np.ndarray:
    """In-plane nearest-neighbour row of hexagonal ice.

    That row is the basal [11-20] direction. Looking exactly along it
    stacks each bilayer into columns instead of a doubled pair.
    """
    chosen = np.flatnonzero(chill == 1)
    pts = pos[chosen]
    cell = 3.3
    keys = np.floor(pts / cell).astype(int)
    bins: dict[tuple[int, int, int], list[int]] = {}
    for index, key in enumerate(map(tuple, keys)):
        bins.setdefault(key, []).append(index)
    total = np.zeros(2)
    count = 0
    for index, key in enumerate(map(tuple, keys)):
        for dx in (-1, 0, 1):
            for dy in (-1, 0, 1):
                for dz in (-1, 0, 1):
                    for other in bins.get((key[0] + dx, key[1] + dy, key[2] + dz), []):
                        if other <= index:
                            continue
                        delta = pts[other] - pts[index]
                        delta -= box * np.round(delta / box)
                        length = float(np.linalg.norm(delta))
                        if not (2.4 < length < 3.15 and abs(delta[2]) < 1.0):
                            continue
                        xy = delta[:2]
                        if xy[0] < 0.0:
                            xy = -xy
                        angle = math.degrees(math.atan2(float(xy[1]), float(xy[0]))) % 180.0
                        if angle < 12.0 or angle > 168.0:
                            total += xy
                            count += 1
    if count < 8:
        return np.array([1.0, 0.0, 0.0])
    direction = np.array([total[0], total[1], 0.0])
    return unit(direction)


def stacking_axis(pos: np.ndarray, mask: np.ndarray) -> np.ndarray:
    """Smallest eigenvector of the ice-ice bond tensor.

    A hexagonal or cubic cluster has three basal bond directions and one
    stacking direction, so the stacking axis is the least populated.
    """
    chosen = np.flatnonzero(mask)
    pts = pos[chosen]
    if pts.shape[0] < 8:
        return np.array([0.0, 0.0, 1.0])
    tensor = np.zeros((3, 3))
    cell = HBOND_A
    keys = np.floor(pts / cell).astype(int)
    bins: dict[tuple[int, int, int], list[int]] = {}
    for index, key in enumerate(map(tuple, keys)):
        bins.setdefault(key, []).append(index)
    bonds = 0
    for index, key in enumerate(map(tuple, keys)):
        for dx in (-1, 0, 1):
            for dy in (-1, 0, 1):
                for dz in (-1, 0, 1):
                    for other in bins.get((key[0] + dx, key[1] + dy, key[2] + dz), []):
                        if other <= index:
                            continue
                        delta = pts[other] - pts[index]
                        length = float(np.linalg.norm(delta))
                        if length <= 0.0 or length > HBOND_A:
                            continue
                        hat = delta / length
                        tensor += np.outer(hat, hat)
                        bonds += 1
    if bonds < 4:
        return np.array([0.0, 0.0, 1.0])
    _values, vectors = np.linalg.eigh(tensor)
    axis = vectors[:, 0]
    if axis[2] < 0.0:
        axis = -axis
    return unit(axis)


def layer_kinds(pos: np.ndarray, chill: np.ndarray, mask: np.ndarray, axis: np.ndarray) -> list[int]:
    """One CHILL+ class per 3.7 A bin along the stacking axis."""
    pts = pos[mask]
    kinds = chill[mask]
    if pts.shape[0] == 0:
        return []
    projection = pts @ axis
    bins = np.floor((projection - projection.min()) / 3.7).astype(int)
    letters = []
    for index in range(int(bins.max()) + 1):
        selected = kinds[bins == index]
        if selected.size == 0:
            continue
        counts = {code: int(np.sum(selected == code)) for code in (1, 0, 3)}
        winner = max(counts, key=counts.get)
        if counts[winner] == 0:
            continue
        letters.append(winner)
    return letters


def slice_mask(pos: np.ndarray, chill: np.ndarray, forward: np.ndarray, crystal: np.ndarray) -> np.ndarray:
    """Crystal in full, liquid-like classes only inside a 1 nm depth slice."""
    depth = pos @ unit(forward)
    centre = float(np.median(depth[crystal])) if np.any(crystal) else float(np.median(depth))
    in_slice = np.abs(depth - centre) <= (0.5 * SLICE_A)
    liquid_like = np.isin(chill, [2, 4, 5])
    return crystal | (liquid_like & in_slice)


def contained_bonds(pos: np.ndarray, box: np.ndarray | None, mask: np.ndarray) -> np.ndarray:
    """H-bonds whose straight segment stays inside the drawn coordinates."""
    if box is None:
        pairs = []
        chosen = np.flatnonzero(mask)
        pts = pos[chosen]
        cut2 = HBOND_A * HBOND_A
        cell = HBOND_A
        keys = np.floor(pts / cell).astype(int)
        bins: dict[tuple[int, int, int], list[int]] = {}
        for index, key in enumerate(map(tuple, keys)):
            bins.setdefault(key, []).append(index)
        for index, key in enumerate(map(tuple, keys)):
            for dx in (-1, 0, 1):
                for dy in (-1, 0, 1):
                    for dz in (-1, 0, 1):
                        for other in bins.get((key[0] + dx, key[1] + dy, key[2] + dz), []):
                            if other <= index:
                                continue
                            delta = pts[other] - pts[index]
                            if float(np.dot(delta, delta)) <= cut2:
                                pairs.append((int(chosen[index]), int(chosen[other])))
        if not pairs:
            return np.zeros((0, 2), dtype=np.int32)
        return np.array(pairs, dtype=np.int32)
    raw = minimum_image_pairs(pos, box, mask, HBOND_A)
    kept = []
    for first, second in raw:
        delta = pos[second] - pos[first]
        if float(np.dot(delta, delta)) <= (HBOND_A * 1.02) ** 2:
            kept.append((int(first), int(second)))
    if not kept:
        return np.zeros((0, 2), dtype=np.int32)
    return np.array(kept, dtype=np.int32)


def depth_cue(colours: np.ndarray, depth: np.ndarray) -> np.ndarray:
    span = float(depth.max() - depth.min())
    near = 1.0 if span < 1e-6 else 1.0 - (depth - depth.min()) / span
    factor = 0.78 + 0.22 * near
    return colours * factor[:, None]


def scene(pos, chill, bonds, liquid_alpha: float, forward: np.ndarray):
    from ovito.data import DataCollection, Particles, SimulationCell
    from ovito.pipeline import Pipeline, StaticSource
    from ovito.vis import BondsVis

    data = DataCollection()
    lo = pos.min(axis=0) - 4.0
    hi = pos.max(axis=0) + 4.0
    cell = SimulationCell(pbc=(False, False, False))
    cell[:, :3] = np.diag((hi - lo).astype(float))
    cell[:, 3] = lo
    cell.vis.enabled = False
    data.objects.append(cell)
    particles = Particles()
    particles.create_property("Position", data=np.ascontiguousarray(pos, dtype=np.float64))
    colours = np.array([CHILL.get(int(code), ("other", (0.5, 0.5, 0.5)))[1] for code in chill])
    colours = depth_cue(colours, pos @ unit(forward))
    particles.create_property("Color", data=np.clip(colours, 0.0, 1.0))
    liquid = np.isin(chill, [2, 4, 5])
    transparency = np.where(liquid, 1.0 - liquid_alpha, 0.0)
    particles.create_property("Transparency", data=transparency.astype(np.float64))
    particles.create_property("Radius", data=np.full(pos.shape[0], ATOM_RADIUS))
    if bonds is not None and len(bonds):
        topology = particles.create_bonds()
        topology.create_property("Topology", data=np.ascontiguousarray(bonds, dtype=np.int32))
        topology.vis.width = BOND_WIDTH
        topology.vis.color = (0.12, 0.12, 0.12)
        topology.vis.coloring_mode = BondsVis.ColoringMode.Uniform
    data.objects.append(particles)
    particles.vis.radius = ATOM_RADIUS
    return Pipeline(source=StaticSource(data=data))


def project(pos: np.ndarray, forward: np.ndarray, up: np.ndarray):
    forward, right, up = basis(forward, up)
    origin = pos.mean(axis=0)
    relative = pos - origin
    horizontal = relative @ right
    vertical = relative @ up
    margin = ATOM_RADIUS + 0.5 * BOND_WIDTH + 2.2
    h_lo, h_hi = float(horizontal.min()) - margin, float(horizontal.max()) + margin
    v_lo, v_hi = float(vertical.min()) - margin, float(vertical.max()) + margin
    centre = origin + right * (0.5 * (h_lo + h_hi)) + up * (0.5 * (v_lo + v_hi))
    return centre, 0.5 * (h_hi - h_lo), 0.5 * (v_hi - v_lo), forward, right, up


def render_view(pos, chill, bonds, forward, up, png: Path, half_h: float, half_v: float, width_px: int):
    from ovito.vis import TachyonRenderer, Viewport

    forward, _right, up = basis(forward, up)
    centre, _own_h, _own_v, forward, _right, up = project(pos, forward, up)
    # The caller passes the shared frame so panels share a scale. Re-centre
    # on this cloud so a smaller cluster is not pushed against the edge.
    height_px = max(int(round(width_px * half_v / half_h)), 32)
    pipeline = scene(pos, chill, bonds, 0.0, forward)
    # liquid alpha is already in the subset; scene uses liquid_alpha 0 here
    # and the caller encoded transparency by passing prefiltered particles.
    view = Viewport(type=Viewport.Type.Ortho)
    pipeline.add_to_scene()
    try:
        view.camera_dir = tuple(float(v) for v in forward)
        view.camera_up = tuple(float(v) for v in up)
        view.camera_pos = centre - forward * (6.0 * max(half_h, half_v))
        view.fov = float(half_v)
        renderer = TachyonRenderer(
            shadows=True,
            ambient_occlusion=True,
            ambient_occlusion_samples=8,
            antialiasing_samples=6,
            direct_light_intensity=0.82,
        )
        png.parent.mkdir(parents=True, exist_ok=True)
        view.render_image(
            filename=str(png),
            size=(width_px, height_px),
            renderer=renderer,
            background=(1.0, 1.0, 1.0),
            frame=0,
        )
    finally:
        pipeline.remove_from_scene()
    return height_px


def prepare_slab(pos, chill, box):
    forward = basal_forward(pos, chill, box)
    up = np.array([0.0, 0.0, 1.0])
    crystal = np.isin(chill, [0, 1, 3])
    keep = slice_mask(pos, chill, forward, crystal)
    # Clathrate classes stay opaque; only true liquid is faint.
    faint = keep & (chill == 2)
    solid = keep & ~faint
    order = np.concatenate([np.flatnonzero(faint), np.flatnonzero(solid)])
    # Draw solid after liquid by listing liquid first (painter). OVITO does
    # not use that order; transparency handles the liquid.
    drawn = pos[keep]
    drawn_chill = chill[keep]
    bond_source = np.isin(chill, [0, 1, 3]) & keep
    bonds = contained_bonds(pos, box, bond_source)
    remap = -np.ones(pos.shape[0], dtype=int)
    remap[np.flatnonzero(keep)] = np.arange(int(keep.sum()))
    if len(bonds):
        bonds = np.stack([remap[bonds[:, 0]], remap[bonds[:, 1]]], axis=1)
    return drawn, drawn_chill, bonds, forward, up


def prepare_nucleus(pos, chill, extended, forward, up):
    crystal = extended if extended is not None else np.isin(chill, [0, 1, 3])
    depth = pos @ unit(forward)
    centre = float(np.median(depth[crystal])) if np.any(crystal) else 0.0
    liquid = (chill == 2) & (np.abs(depth - centre) <= 0.5 * SLICE_A)
    keep = crystal | liquid
    drawn = pos[keep]
    drawn_chill = chill[keep]
    local = np.zeros(pos.shape[0], dtype=bool)
    local[np.flatnonzero(keep)] = crystal[keep]
    bonds = contained_bonds(pos, None, local)
    remap = -np.ones(pos.shape[0], dtype=int)
    remap[np.flatnonzero(keep)] = np.arange(int(keep.sum()))
    if len(bonds):
        bonds = np.stack([remap[bonds[:, 0]], remap[bonds[:, 1]]], axis=1)
    # Faint liquid: rewrite the chill of liquid to stay class 2 and let scene
    # transparency do the work. Mark them by returning a transparency vector.
    return drawn, drawn_chill, bonds


def render_molecule(pos, chill, bonds, forward, up, png, half_h, half_v, width_px, liquid_alpha):
    from ovito.data import DataCollection, Particles, SimulationCell
    from ovito.pipeline import Pipeline, StaticSource
    from ovito.vis import BondsVis, TachyonRenderer, Viewport

    forward, _right, up = basis(forward, up)
    centre, _h, _v, forward, _right, up = project(pos, forward, up)
    height_px = max(int(round(width_px * half_v / max(half_h, 1e-6))), 32)
    data = DataCollection()
    lo = pos.min(axis=0) - 6.0
    hi = pos.max(axis=0) + 6.0
    cell = SimulationCell(pbc=(False, False, False))
    cell[:, :3] = np.diag((hi - lo).astype(float))
    cell[:, 3] = lo
    cell.vis.enabled = False
    data.objects.append(cell)
    particles = Particles()
    particles.create_property("Position", data=np.ascontiguousarray(pos, dtype=np.float64))
    colours = np.array([CHILL.get(int(code), ("other", (0.5, 0.5, 0.5)))[1] for code in chill])
    colours = depth_cue(colours, pos @ forward)
    particles.create_property("Color", data=np.clip(colours, 0.0, 1.0))
    liquid = np.isin(chill, [2])
    guest = np.isin(chill, [4, 5])
    transparency = np.zeros(pos.shape[0])
    transparency[liquid] = 1.0 - liquid_alpha
    transparency[guest] = 0.15
    particles.create_property("Transparency", data=transparency)
    particles.create_property("Radius", data=np.full(pos.shape[0], ATOM_RADIUS))
    if bonds is not None and len(bonds):
        bond_data = particles.create_bonds()
        bond_data.create_property("Topology", data=np.ascontiguousarray(bonds, dtype=np.int32))
        bond_data.vis.width = BOND_WIDTH
        bond_data.vis.color = (0.10, 0.10, 0.10)
        bond_data.vis.coloring_mode = BondsVis.ColoringMode.Uniform
    data.objects.append(particles)
    particles.vis.radius = ATOM_RADIUS
    pipeline = Pipeline(source=StaticSource(data=data))
    view = Viewport(type=Viewport.Type.Ortho)
    pipeline.add_to_scene()
    try:
        view.camera_dir = tuple(float(v) for v in forward)
        view.camera_up = tuple(float(v) for v in up)
        view.camera_pos = centre - forward * (8.0 * max(half_h, half_v))
        view.fov = float(half_v)
        renderer = TachyonRenderer(
            shadows=True,
            ambient_occlusion=True,
            ambient_occlusion_samples=8,
            antialiasing_samples=6,
            direct_light_intensity=0.8,
        )
        png.parent.mkdir(parents=True, exist_ok=True)
        view.render_image(
            filename=str(png),
            size=(width_px, height_px),
            renderer=renderer,
            background=(1.0, 1.0, 1.0),
            frame=0,
        )
    finally:
        pipeline.remove_from_scene()
    return height_px


def compose(slab, nuclei, stacks, width_mm: float, png: Path, pdf: Path) -> None:
    import matplotlib

    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    from PIL import Image

    dpi = PRINT_DPI
    # Right margin holds the z scale bar. Bottom margin holds the legend.
    margin_l, margin_r, margin_t, margin_b = 10.0, 4.0, 5.0, 18.0
    gap = 3.0
    if width_mm >= 150.0:
        slab_w = 70.0
        content_w = width_mm - margin_l - margin_r
        nucleus_w = (content_w - slab_w - gap) / 2.0
        slab_h = slab_w * slab["half_v"] / slab["half_h"]
        nucleus_h = nucleus_w * nuclei["half_v"] / nuclei["half_h"]
        bar_h = 12.0
        right_h = nucleus_h + 6.0 + bar_h
        height_mm = margin_t + margin_b + max(slab_h, right_h)
        slab_x, slab_y = margin_l, margin_b + max(0.0, right_h - slab_h)
        nuc_x = margin_l + slab_w + gap
        nuc_y = margin_b + bar_h + 6.0 + max(0.0, slab_h - right_h)
    else:
        content_w = width_mm - margin_l - margin_r
        slab_w = content_w
        slab_h = slab_w * slab["half_v"] / slab["half_h"]
        nucleus_w = (content_w - gap) / 2.0
        nucleus_h = nucleus_w * nuclei["half_v"] / nuclei["half_h"]
        bar_h = 16.0
        height_mm = margin_t + margin_b + slab_h + gap + nucleus_h + gap + bar_h
        slab_x, slab_y = margin_l, margin_b + bar_h + gap + nucleus_h + gap
        nuc_x, nuc_y = margin_l, margin_b + bar_h + gap
    if height_mm > 230.0:
        raise SystemExit(f"figure is {height_mm:.0f} mm tall, past a CPC page")

    fig = plt.figure(figsize=(width_mm / 25.4, height_mm / 25.4), dpi=dpi)
    fig.patch.set_facecolor("white")

    def place(image_path, x_mm, y_mm, w_mm, h_mm):
        axis = fig.add_axes((x_mm / width_mm, y_mm / height_mm, w_mm / width_mm, h_mm / height_mm))
        axis.imshow(np.asarray(Image.open(image_path)), interpolation="nearest")
        axis.set_axis_off()
        return axis

    def letter(char, x_mm, y_mm):
        fig.text(
            x_mm / width_mm,
            y_mm / height_mm,
            char,
            fontsize=9,
            fontweight="bold",
            ha="left",
            va="top",
            color="black",
        )

    place(slab["png"], slab_x, slab_y, slab_w, slab_h)
    letter("(a)", slab_x, slab_y + slab_h + 0.3)
    angstrom_per_mm = (2.0 * slab["half_v"]) / slab_h
    bar_mm = 10.0 / angstrom_per_mm
    bar_x = 2.2 / width_mm
    bar_y = (slab_y + 0.5 * (slab_h - bar_mm)) / height_mm
    fig.add_artist(
        plt.Line2D(
            [bar_x, bar_x],
            [bar_y, bar_y + bar_mm / height_mm],
            transform=fig.transFigure,
            color="black",
            lw=1.1,
            solid_capstyle="butt",
        )
    )
    fig.text(
        bar_x + 0.9 / width_mm,
        bar_y + 0.5 * bar_mm / height_mm,
        "z\n10 Å",
        fontsize=7,
        ha="left",
        va="center",
        linespacing=1.15,
    )

    for offset, key, tag in ((0.0, "minimum", "minimum"), (nucleus_w + gap, "late", "late")):
        place(nuclei[key], nuc_x + offset, nuc_y, nucleus_w, nucleus_h)
    letter("(b)", nuc_x, nuc_y + nucleus_h + 3.2)
    fig.text(
        (nuc_x + 0.5 * nucleus_w) / width_mm,
        (nuc_y + nucleus_h + 0.15) / height_mm,
        "minimum",
        fontsize=7,
        ha="center",
        va="bottom",
    )
    fig.text(
        (nuc_x + nucleus_w + gap + 0.5 * nucleus_w) / width_mm,
        (nuc_y + nucleus_h + 0.15) / height_mm,
        "late",
        fontsize=7,
        ha="center",
        va="bottom",
    )

    if width_mm >= 150.0:
        bar_left = nuc_x + 16.0
        bar_width = 2 * nucleus_w + gap - 16.0
        bar_bottom = nuc_y - 6.0 - bar_h
    else:
        bar_left = margin_l + 16.0
        bar_width = content_w - 16.0
        bar_bottom = margin_b
    letter("(c)", bar_left - 16.0, bar_bottom + bar_h)
    axis = fig.add_axes(
        (
            bar_left / width_mm,
            (bar_bottom + 1.0) / height_mm,
            bar_width / width_mm,
            (bar_h - 4.0) / height_mm,
        )
    )
    for row, (name, kinds) in enumerate((("minimum", stacks["minimum"]), ("late", stacks["late"]))):
        for column, kind in enumerate(kinds):
            colour = CHILL.get(kind, ("", (0.5, 0.5, 0.5)))[1]
            axis.add_patch(
                plt.Rectangle((column, row), 1.0, 0.72, facecolor=colour, edgecolor="white", lw=0.3)
            )
        fig.text(
            (bar_left - 1.2) / width_mm,
            (bar_bottom + 1.0 + (bar_h - 4.0) * (row + 0.36) / 1.9) / height_mm,
            name,
            fontsize=7,
            ha="right",
            va="center",
        )
    axis.set_xlim(0, max(len(stacks["minimum"]), len(stacks["late"]), 1))
    axis.set_ylim(-0.15, 1.9)
    axis.set_axis_off()
    axis.set_aspect("auto")
    fig.text(
        (bar_left + bar_width) / width_mm,
        (bar_bottom + 0.2) / height_mm,
        "stacking",
        fontsize=7,
        ha="right",
        va="bottom",
    )

    handles = [
        plt.Line2D([0], [0], marker="o", color="none", markerfacecolor=colour, markeredgecolor="none", markersize=5, label=name)
        for _code, (name, colour) in CHILL.items()
    ]
    fig.legend(
        handles=handles,
        loc="lower center",
        ncol=3 if width_mm < 150.0 else 6,
        frameon=False,
        fontsize=7,
        bbox_to_anchor=(0.5, 0.028),
    )
    fig.text(
        0.5,
        0.008,
        "Sky blue: clathrate. Purple: interfacial clathrate.",
        fontsize=6.5,
        ha="center",
        va="bottom",
        color="#333333",
    )
    png.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(png, dpi=dpi)
    fig.savefig(pdf)
    plt.close(fig)


def build(interface_dump: Path, interface_frame: int, nucleus_dump: Path, min_frame: int, late_frame: int, out: Path, work: Path) -> None:
    work.mkdir(parents=True, exist_ok=True)
    box, pos = read_dump_frame(interface_dump, interface_frame)
    chill, _extended = label_positions(box, pos, work / "slab", "interface")
    slab_pos, slab_chill, slab_bonds, forward, up = prepare_slab(pos, chill, box)
    _centre, half_h, half_v, forward, _right, up = project(slab_pos, forward, up)
    slab_px = int(round(72.0 / 25.4 * PRINT_DPI))
    slab_png = work / "slab.png"
    render_molecule(slab_pos, slab_chill, slab_bonds, forward, up, slab_png, half_h, half_v, slab_px, 0.42)
    slab = {"png": slab_png, "half_h": half_h, "half_v": half_v}

    n_frames = count_frames(nucleus_dump)
    if late_frame < 0:
        late_frame = n_frames + 1 + late_frame
    clouds = {}
    axes = {}
    for label, frame in (("minimum", min_frame), ("late", late_frame)):
        nbox, npos = read_dump_frame(nucleus_dump, frame)
        nchill, extended = label_positions(nbox, npos, work / label, "nucleation")
        if extended is None:
            extended = np.isin(nchill, [0, 1, 3])
        npos, _nbox = unwrap_cluster(npos, nbox, extended, HBOND_A)
        axis = stacking_axis(npos, extended)
        hint = np.array([1.0, 0.0, 0.0])
        if abs(float(np.dot(hint, axis))) > 0.85:
            hint = np.array([0.0, 1.0, 0.0])
        view = unit(np.cross(axis, hint))
        clouds[label] = (npos, nchill, extended, view, axis)
        axes[label] = axis
    # Shared orthographic scale from the late cluster, which is the larger one.
    halves = []
    prepared = {}
    for label, (npos, nchill, extended, view, axis) in clouds.items():
        drawn, drawn_chill, bonds = prepare_nucleus(npos, nchill, extended, view, axis)
        _c, half_h, half_v, view, _r, axis = project(drawn, view, axis)
        prepared[label] = (drawn, drawn_chill, bonds, view, axis, half_h, half_v)
        halves.append((half_h, half_v))
    share_h = max(item[0] for item in halves)
    share_v = max(item[1] for item in halves)
    nucleus_px = int(round(48.0 / 25.4 * PRINT_DPI))
    images = {}
    stacks = {}
    for label, (drawn, drawn_chill, bonds, view, axis, _h, _v) in prepared.items():
        images[label] = work / f"nucleus-{label}.png"
        render_molecule(drawn, drawn_chill, bonds, view, axis, images[label], share_h, share_v, nucleus_px, 0.12)
        npos, nchill, extended, _view, axis = clouds[label]
        stacks[label] = layer_kinds(npos, nchill, extended, axis)
    nuclei = {"minimum": images["minimum"], "late": images["late"], "half_h": share_h, "half_v": share_v}
    for width, stem in ((170.0, "figure-170mm"), (85.0, "figure-85mm")):
        compose(slab, nuclei, stacks, width, out / f"{stem}.png", out / f"{stem}.pdf")
        print(f"wrote {stem}", flush=True)


def main() -> int:
    ensure_runtime()
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--out", type=Path, required=True)
    parser.add_argument("--interface-dump", type=Path, required=True)
    parser.add_argument("--interface-frame", type=int, default=1)
    parser.add_argument("--nucleus-dump", type=Path, required=True)
    parser.add_argument("--nucleus-min-frame", type=int, default=2)
    parser.add_argument("--nucleus-late-frame", type=int, default=-1)
    parser.add_argument("--work", type=Path, default=ROOT / "data" / "science-frames" / "render-work" / "journal")
    args = parser.parse_args()
    args.out.mkdir(parents=True, exist_ok=True)
    frame = args.interface_frame
    if frame < 0:
        frame = count_frames(args.interface_dump) + 1 + frame
    build(
        args.interface_dump,
        frame,
        args.nucleus_dump,
        args.nucleus_min_frame,
        args.nucleus_late_frame,
        args.out,
        args.work,
    )
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
