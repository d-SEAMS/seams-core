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
# Okabe-Ito for the ice classes. Clathrate is Paul Tol green and interfacial
# clathrate is wine, so neither sits on the cubic blue after a deuteranopia
# simulation. In the 230 K slab the green and wine molecules in the liquid are
# clathrate and interfacial clathrate. The cutoffs are not refit.
CHILL = {
    0: ("cubic Ic", (0.000, 0.447, 0.698)),
    1: ("hexagonal Ih", (0.902, 0.624, 0.000)),
    3: ("interfacial", (0.941, 0.894, 0.259)),
    2: ("liquid", (0.940, 0.940, 0.940)),
    4: ("clathrate", (0.133, 0.533, 0.200)),
    5: ("interfacial clathrate", (0.533, 0.133, 0.333)),
}
ICE = {0, 1, 3, 4, 5}
# One radius and one bond width for every panel. 85 mm is the CPC column.
ATOM_RADIUS = 0.90
BOND_WIDTH = 0.25
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
    name = "nucleus-atoms.csv" if command == "nucleation" else "interface-atoms.csv"
    out = work / "labels"
    cached = out / name
    if os.environ.get("SCIENCE_REUSE_LABELS") == "1" and cached.is_file():
        header, rows, _ = read_table(cached)
        chill = column(header, rows, "chill").astype(int)
        if chill.shape[0] != pos.shape[0]:
            raise SystemExit(f"{cached} has {chill.shape[0]} rows for {pos.shape[0]} atoms")
        extended = None
        if "in_extended" in header:
            extended = column(header, rows, "in_extended") > 0.5
        print(f"reuse labels {cached}", flush=True)
        return chill, extended
    dump = work / "one.lammpstrj"
    write_one_frame(dump, box, pos)
    run_front(dump, out, command)
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


def _ice_points(pos: np.ndarray, chill: np.ndarray) -> np.ndarray:
    return pos[np.isin(chill, [0, 1])]


def basal_forward(pos: np.ndarray, chill: np.ndarray, box: np.ndarray) -> np.ndarray:
    """[11-20]: the in-plane nearest-neighbour direction of hexagonal ice.

    The three bond azimuths are 60 degrees apart. Their mean, not the box
    axis, is the view. Front and back oxygens of a column then coincide.
    """
    chosen = np.flatnonzero(chill == 1)
    pts = pos[chosen]
    cell = 3.3
    keys = np.floor(pts / cell).astype(int)
    bins: dict[tuple[int, int, int], list[int]] = {}
    for index, key in enumerate(map(tuple, keys)):
        bins.setdefault(key, []).append(index)
    angles = []
    vectors = []
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
                        if not (2.4 < length < 3.15 and abs(float(delta[2])) < 1.0):
                            continue
                        xy = delta[:2]
                        if xy[0] < 0.0:
                            xy = -xy
                        angles.append(math.degrees(math.atan2(float(xy[1]), float(xy[0]))) % 180.0)
                        vectors.append(xy)
    if len(angles) < 8:
        return np.array([1.0, 0.0, 0.0])
    folded = np.array(angles) % 60.0
    hist, edges = np.histogram(folded, bins=12, range=(0.0, 60.0))
    centre = float(edges[int(np.argmax(hist))] + 2.5)
    # The three [11-20] axes are this peak and the same peak plus 60 and 120.
    members = []
    for angle, xy in zip(angles, vectors):
        delta = (angle - centre + 30.0) % 60.0 - 30.0
        if abs(delta) <= 8.0:
            members.append(xy)
    # Average only one of the three [11-20] axes. Folding all three
    # together points the camera between lattice rows.
    best_count = -1
    best_mean = np.array([1.0, 0.0])
    for shift in (0.0, 60.0, 120.0):
        target = (centre + shift) % 180.0
        picked = []
        for angle, xy in zip(angles, vectors):
            delta = (angle - target + 90.0) % 180.0 - 90.0
            if abs(delta) <= 8.0:
                picked.append(xy if float(xy[0]) >= 0.0 else -xy)
        if len(picked) > best_count:
            best_count = len(picked)
            best_mean = np.mean(np.asarray(picked), axis=0)
    return unit(np.array([best_mean[0], best_mean[1], 0.0]))


def _column_offsets(pts: np.ndarray, direction: np.ndarray) -> np.ndarray:
    """In-plane distance from each core oxygen to its partner one lattice step along the view.

    Ice Ih repeats every a ≈ 4.5 Å along [11-20]. That partner is the other
    oxygen of the column.
    """
    direction = unit(direction)
    _forward, right, up = basis(direction, np.array([0.0, 0.0, 1.0]))
    horizontal = pts @ right
    vertical = pts @ up
    depth = pts @ direction
    offsets = []
    for index in range(len(pts)):
        delta_depth = np.abs(depth - depth[index])
        plane = np.hypot(horizontal - horizontal[index], vertical - vertical[index])
        chosen = (delta_depth > 3.8) & (delta_depth < 5.1) & (plane < 1.2)
        if np.any(chosen):
            offsets.append(float(np.min(plane[chosen])))
    if not offsets:
        return np.array([9.0])
    return np.asarray(offsets, dtype=float)


def column_axis(pos: np.ndarray, chill: np.ndarray, box: np.ndarray) -> tuple[np.ndarray, float]:
    """[11-20], the basal row, from the oxygen columns.

    The camera looks along the lattice repeat, which lies in the basal
    plane. A nearest-neighbour bond is the wrong axis: it is tilted by
    the tetrahedral angle, and the two oxygens of a column then sit about
    half a diameter apart. The returned spread is the median in-plane
    offset of the column partner.
    """
    core = np.flatnonzero(np.isin(chill, [0, 1]) & (pos[:, 2] > 8.0) & (pos[:, 2] < 26.0))
    if core.size < 30:
        core = np.flatnonzero(np.isin(chill, [0, 1]))
    pts = pos[core]
    if pts.shape[0] < 10:
        return basal_forward(pos, chill, box), 9.0

    def candidate(azimuth: float, tilt: float) -> tuple[np.ndarray, np.ndarray]:
        theta = math.radians(azimuth)
        phi = math.radians(tilt)
        in_plane = np.array([math.cos(theta), math.sin(theta), 0.0])
        direction = math.cos(phi) * in_plane + math.sin(phi) * np.array([0.0, 0.0, 1.0])
        return unit(direction), _column_offsets(pts, direction)

    best: tuple | None = None
    for azimuth in np.linspace(0.0, 180.0, 180, endpoint=False):
        direction, offsets = candidate(float(azimuth), 0.0)
        median = float(np.median(offsets))
        # More matched columns, then a tighter median.
        rank = (-int(np.sum(offsets < 1.0)), median)
        if best is None or rank < best[0]:
            best = (rank, float(azimuth), direction, offsets)
    assert best is not None
    azimuth = best[1]
    refined: tuple | None = None
    for delta in np.linspace(-2.0, 2.0, 17):
        for tilt in np.linspace(-2.0, 2.0, 9):
            direction, offsets = candidate(azimuth + float(delta), float(tilt))
            median = float(np.median(offsets))
            rank = (float(np.percentile(offsets, 90)), median, -int(np.sum(offsets < 1.0)))
            if refined is None or rank < refined[0]:
                refined = (rank, direction, offsets)
    assert refined is not None
    direction, offsets = refined[1], refined[2]
    median = float(np.median(offsets))
    print(
        "column axis "
        f"{np.round(direction, 3)} median {median:.3f} A "
        f"p90 {float(np.percentile(offsets, 90)):.3f} A "
        f"max {float(np.max(offsets)):.3f} A",
        flush=True,
    )
    return direction, median


def layer_period(projection: np.ndarray) -> float:
    """Basal repeat along a stacking axis, from the strongest density wave."""
    best_score = -1.0
    best = 3.6
    centred = projection - float(np.mean(projection))
    for period in np.linspace(3.30, 4.05, 16):
        wave = np.exp(1j * (2.0 * np.pi / float(period)) * centred)
        score = float(np.abs(np.mean(wave)))
        if score > best_score:
            best_score = score
            best = float(period)
    return best


def layer_kinds(pos: np.ndarray, chill: np.ndarray, mask: np.ndarray, axis: np.ndarray) -> list[int]:
    """Ih (1) or Ic (0) for each basal layer. Interfacial molecules do not colour a layer."""
    chosen = np.flatnonzero(mask & np.isin(chill, [0, 1]))
    if chosen.size == 0:
        return []
    projection = pos[chosen] @ unit(axis)
    kinds = chill[chosen]
    bins = np.floor((projection - projection.min()) / layer_period(projection)).astype(int)
    letters = []
    for index in range(int(bins.max()) + 1):
        selected = kinds[bins == index]
        hexagonal = int(np.sum(selected == 1))
        cubic = int(np.sum(selected == 0))
        if hexagonal + cubic == 0:
            continue
        letters.append(1 if hexagonal >= cubic else 0)
    return letters


def hexagonality(chill: np.ndarray, mask: np.ndarray) -> float:
    hexagonal = int(np.sum(mask & (chill == 1)))
    cubic = int(np.sum(mask & (chill == 0)))
    if hexagonal + cubic == 0:
        return float("nan")
    return hexagonal / float(hexagonal + cubic)


def cubicity(kinds: list[int]) -> float:
    """Cubic share of the basal layers, one majority class per layer."""
    if not kinds:
        return float("nan")
    return sum(1 for kind in kinds if kind == 0) / float(len(kinds))


def stacking_axis(pos: np.ndarray, chill: np.ndarray, mask: np.ndarray) -> np.ndarray:
    """Normal of the basal planes, from the density wave of the Ih/Ic lattice.

    A purity score on fixed bins can sit a few degrees off the planes. The
    layers then smear into each other and the edge-on view is not a stripe.
    """
    lattice = mask & np.isin(chill, [0, 1])
    pts = pos[lattice]
    if pts.shape[0] < 8:
        return np.array([0.0, 0.0, 1.0])

    def amplitude(axis: np.ndarray, period: float) -> float:
        projection = pts @ axis
        wave = np.exp(1j * (2.0 * np.pi / period) * (projection - projection.mean()))
        return float(np.abs(np.mean(wave)))

    best_score = -1.0
    best_axis = np.array([0.0, 0.0, 1.0])
    best_period = 3.6
    for theta in np.linspace(0.0, np.pi, 24, endpoint=False):
        for phi in np.linspace(0.0, np.pi, 24, endpoint=False):
            axis = np.array(
                [np.sin(theta) * np.cos(phi), np.sin(theta) * np.sin(phi), np.cos(theta)]
            )
            for period in (3.40, 3.60, 3.80):
                score = amplitude(axis, period)
                if score > best_score:
                    best_score = score
                    best_axis = axis
                    best_period = period
    for scale in (0.12, 0.04, 0.015):
        tangent = np.array([1.0, 0.0, 0.0]) if abs(float(best_axis[2])) > 0.9 else np.array([0.0, 0.0, 1.0])
        east = unit(np.cross(tangent, best_axis))
        north = np.cross(best_axis, east)
        for step_east in np.linspace(-scale, scale, 7):
            for step_north in np.linspace(-scale, scale, 7):
                axis = unit(best_axis + step_east * east + step_north * north)
                for period in np.linspace(best_period - 0.12, best_period + 0.12, 5):
                    score = amplitude(axis, float(period))
                    if score > best_score:
                        best_score = score
                        best_axis = axis
                        best_period = float(period)
    if best_axis[2] < 0.0:
        best_axis = -best_axis
    return best_axis


def bilayer_sheet(pos: np.ndarray, crystal: np.ndarray, forward: np.ndarray, box: np.ndarray) -> np.ndarray:
    """Both oxygens of each puckered basal pair, about 2.8 A along [11-20].

    Hydrogen-bond partners just outside that window are kept as well, so a
    bond that crosses the slice boundary is not left dangling.
    """
    depth = pos @ unit(forward)
    ice_depth = depth[crystal]
    window = 2.8
    best_lo = float(ice_depth.min()) if ice_depth.size else 0.0
    best_n = -1
    if ice_depth.size:
        for lo in np.linspace(float(ice_depth.min()), float(ice_depth.max()) - window, 60):
            count = int(np.sum((ice_depth >= lo) & (ice_depth < lo + window)))
            if count > best_n:
                best_n = count
                best_lo = float(lo)
    inside = crystal & (depth >= best_lo) & (depth < best_lo + window)
    partners = minimum_image_pairs(pos, box, crystal, HBOND_A)
    extra = np.zeros(pos.shape[0], dtype=bool)
    for first, second in partners:
        if inside[first]:
            extra[second] = True
        if inside[second]:
            extra[first] = True
    extra &= ~inside
    # The oxygen one lattice step along [11-20] projects onto a column
    # that is already drawn. Keeping it paints a second sphere on the
    # first. A partner that lands on a new vertex stays, so a bond that
    # crosses the slice still has both ends.
    if np.any(inside) and np.any(extra):
        _forward, right, up = basis(unit(forward), np.array([0.0, 0.0, 1.0]))
        horizontal = pos @ right
        vertical = pos @ up
        centre = float(np.median(depth[inside]))
        order = np.flatnonzero(extra)[np.argsort(np.abs(depth[extra] - centre))]
        kept = list(np.flatnonzero(inside))
        for atom in order:
            plane = np.hypot(horizontal[kept] - horizontal[atom], vertical[kept] - vertical[atom])
            if float(np.min(plane)) < 1.05:
                extra[atom] = False
            else:
                kept.append(int(atom))
    return inside | extra


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


def depth_cue(
    colours: np.ndarray,
    depth: np.ndarray,
    fade: float = 0.22,
    tint: np.ndarray | None = None,
) -> tuple[np.ndarray, np.ndarray]:
    """Report nearness. Far atoms darken, or blend toward a light tint.

    The tint keeps the hue and lifts the far side toward the background.
    It does not push the colour through black.
    """
    span = float(depth.max() - depth.min())
    near = np.ones(len(depth)) if span < 1e-6 else 1.0 - (depth - depth.min()) / span
    if tint is None:
        factor = (1.0 - fade) + fade * near
        return colours * factor[:, None], near
    amount = 0.55 + 0.45 * near
    lifted = amount[:, None] * colours + (1.0 - amount[:, None]) * np.asarray(tint, dtype=float)
    return lifted, near


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
    colours, _near = depth_cue(colours, pos @ unit(forward))
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
    forward, spread = column_axis(pos, chill, box)
    if spread > 0.22:
        raise SystemExit(f"column spread {spread:.3f} A is above 0.22 A")
    up = np.array([0.0, 0.0, 1.0])
    crystal = np.isin(chill, [0, 1, 3])
    sheet = bilayer_sheet(pos, crystal, forward, box)
    ice_z = pos[np.isin(chill, [0, 1]), 2]
    z_mid = float(np.median(ice_z)) if ice_z.size else 0.0
    liquid = np.isin(chill, [2, 4, 5])
    # The drawn ice is one [11-20] slice. The liquid has to be that same
    # slice, minimum-imaged, or the upper contact layer — which sits in
    # this slice and not at the middle of the block — is left out and the
    # picture shows an empty band.
    anchor = pos[sheet].mean(axis=0) if np.any(sheet) else pos[crystal].mean(axis=0)
    delta = pos - anchor
    delta -= box * np.round(delta / box)
    along = delta @ unit(forward)
    in_slice = np.abs(along) <= 0.5 * SLICE_A
    keep = sheet | (liquid & (pos[:, 2] > z_mid) & in_slice) | (liquid & (pos[:, 2] < z_mid) & in_slice)
    # Clathrate classes stay opaque; only true liquid is faint.
    faint = keep & (chill == 2)
    solid = keep & ~faint
    order = np.concatenate([np.flatnonzero(faint), np.flatnonzero(solid)])
    # Draw solid after liquid by listing liquid first (painter). OVITO does
    # not use that order; transparency handles the liquid.
    drawn = pos[keep]
    drawn_chill = chill[keep]
    bonds = contained_bonds(pos, box, sheet)
    remap = -np.ones(pos.shape[0], dtype=int)
    remap[np.flatnonzero(keep)] = np.arange(int(keep.sum()))
    if len(bonds):
        bonds = np.stack([remap[bonds[:, 0]], remap[bonds[:, 1]]], axis=1)
    return drawn, drawn_chill, bonds, forward, up


def ice_bonds(pos: np.ndarray, box: np.ndarray, mask: np.ndarray, cutoff: float = 3.5) -> np.ndarray:
    """Minimum-image oxygen pairs inside the cutoff, at most four per water.

    The four shortest contacts are kept. A fifth neighbour inside the
    cutoff is dropped, and the result is checked.
    """
    chosen = np.flatnonzero(mask)
    if chosen.size < 2:
        return np.zeros((0, 2), dtype=np.int32)
    pts = pos[chosen]
    pairs: list[tuple[int, int, float]] = []
    cut2 = cutoff * cutoff
    cell = cutoff
    shape = np.maximum(1, np.floor(box / cell).astype(int))
    bins: dict[tuple[int, int, int], list[int]] = {}
    keys = np.floor(pts / box * shape).astype(int) % shape
    for index, key in enumerate(keys):
        bins.setdefault(tuple(int(v) for v in key), []).append(index)
    for index, key in enumerate(keys):
        for dx in (-1, 0, 1):
            for dy in (-1, 0, 1):
                for dz in (-1, 0, 1):
                    neighbour = (
                        int((key[0] + dx) % shape[0]),
                        int((key[1] + dy) % shape[1]),
                        int((key[2] + dz) % shape[2]),
                    )
                    for other in bins.get(neighbour, []):
                        if other <= index:
                            continue
                        delta = pts[other] - pts[index]
                        delta -= box * np.round(delta / box)
                        length2 = float(np.dot(delta, delta))
                        if length2 <= cut2:
                            pairs.append((index, other, length2))
    pairs.sort(key=lambda item: item[2])
    degree = np.zeros(chosen.size, dtype=int)
    kept: list[tuple[int, int]] = []
    for first, second, _length2 in pairs:
        if degree[first] >= 4 or degree[second] >= 4:
            continue
        degree[first] += 1
        degree[second] += 1
        kept.append((int(chosen[first]), int(chosen[second])))
    if degree.size and int(degree.max()) > 4:
        raise SystemExit("a water molecule was given more than 4 bonds")
    if not kept:
        return np.zeros((0, 2), dtype=np.int32)
    return np.array(kept, dtype=np.int32)


def column_view(pos: np.ndarray, mask: np.ndarray, stacking: np.ndarray) -> np.ndarray:
    """[11-20] in the basal plane: the direction the columns hide behind."""
    stacking = unit(stacking)
    hint = np.array([1.0, 0.0, 0.0]) if abs(float(stacking[0])) < 0.9 else np.array([0.0, 1.0, 0.0])
    east = unit(np.cross(stacking, hint))
    north = unit(np.cross(stacking, east))
    pts = pos[mask]
    centre = pts.mean(axis=0)
    relative = pts - centre
    best: tuple | None = None
    for degree in np.linspace(0.0, 180.0, 180, endpoint=False):
        theta = math.radians(float(degree))
        view = math.cos(theta) * east + math.sin(theta) * north
        depth = relative @ view
        right = unit(np.cross(stacking, view))
        horizontal = relative @ right
        vertical = relative @ stacking
        offsets = []
        separations = []
        for index in range(0, len(pts), 2):
            delta_depth = np.abs(depth - depth[index])
            plane = np.hypot(horizontal - horizontal[index], vertical - vertical[index])
            near = (delta_depth > 2.2) & (delta_depth < 6.0)
            if not np.any(near):
                continue
            partner = int(np.argmin(np.where(near, plane, 1e9)))
            if plane[partner] < 1.3:
                offsets.append(float(plane[partner]))
                separations.append(float(delta_depth[partner]))
        if len(offsets) < 12:
            continue
        median = float(np.median(offsets))
        rank = (median, -len(offsets))
        if best is None or rank < best[0]:
            best = (rank, view, float(np.median(separations)))
    if best is None:
        return east
    return unit(best[1])


def layer_rows(pos: np.ndarray, chill: np.ndarray, mask: np.ndarray, axis: np.ndarray) -> list[tuple[int, float, int]]:
    """Majority class, purity and count of each basal layer along axis.

    Bins follow the same repeat as the stacking bar. A layer shorter than
    four molecules is an edge fragment and is left out of the purity check.
    """
    axis = unit(axis)
    chosen = np.flatnonzero(mask & np.isin(chill, [0, 1]))
    if chosen.size == 0:
        return []
    projection = pos[chosen] @ axis
    kinds = chill[chosen]
    centred = projection - float(projection.mean())
    best: tuple | None = None
    for period in np.linspace(3.30, 4.05, 16):
        for shift in np.linspace(0.0, float(period), 8, endpoint=False):
            bins = np.floor((centred - shift) / float(period)).astype(int)
            rows = []
            for index in range(int(bins.min()), int(bins.max()) + 1):
                selected = kinds[bins == index]
                if selected.size < 4:
                    continue
                cubic = int(np.sum(selected == 0))
                hexagonal = int(np.sum(selected == 1))
                purity = max(cubic, hexagonal) / float(selected.size)
                rows.append((1 if hexagonal >= cubic else 0, purity, int(selected.size)))
            if len(rows) < 2:
                continue
            score = min(row[1] for row in rows)
            if best is None or score > best[0]:
                best = (score, rows)
    return [] if best is None else best[1]


def outlined_layers(pos: np.ndarray, mask: np.ndarray, axis: np.ndarray, shown: int) -> np.ndarray:
    """Atoms in the first `shown` basal layers, the run outlined on panel (c)."""
    chosen = np.flatnonzero(mask)
    if shown <= 0 or chosen.size == 0:
        return np.ones(pos.shape[0], dtype=bool)
    projection = pos[chosen] @ unit(axis)
    period = layer_period(projection)
    bins = np.floor((projection - projection.min()) / period).astype(int)
    order = [index for index in range(int(bins.min()), int(bins.max()) + 1) if np.any(bins == index)]
    keep = set(order[:shown])
    labels = -np.ones(pos.shape[0], dtype=int)
    labels[chosen] = bins
    return np.isin(labels, list(keep))


def edge_on_view(pos: np.ndarray, mask: np.ndarray, axis: np.ndarray) -> np.ndarray:
    """In-plane direction that lays the c-axis bands across the picture.

    The camera looks along the short basal direction, so each stacking
    layer is a horizontal band and the c-axis of panel (c) is vertical.
    """
    axis = unit(axis)
    hint = np.array([1.0, 0.0, 0.0]) if abs(float(axis[0])) < 0.9 else np.array([0.0, 1.0, 0.0])
    east = unit(np.cross(axis, hint))
    north = unit(np.cross(axis, east))
    pts = pos[mask]
    relative = pts - pts.mean(axis=0)
    best_width = -1.0
    best = east
    for degree in np.linspace(0.0, 180.0, 180, endpoint=False):
        theta = math.radians(float(degree))
        view = math.cos(theta) * east + math.sin(theta) * north
        right = unit(np.cross(axis, view))
        width = float(np.ptp(relative @ right))
        if width > best_width:
            best_width = width
            best = view
    return unit(best)


def orient_nucleus(pos: np.ndarray, chill: np.ndarray, mask: np.ndarray):
    """Local c-axis vertical, view along [11-20], two layers through the centre.

    The axis is this nucleus's own, chosen so the sheet's horizontal layers
    are as pure as the CHILL+ labels allow.
    """
    lattice = mask & np.isin(chill, [0, 1])
    centre = pos[lattice].mean(axis=0)

    def pack(axis: np.ndarray, view: np.ndarray):
        axis = unit(axis)
        view = view - axis * float(view @ axis)
        view = unit(view)
        depth = (pos - centre) @ view
        sheet = lattice & (np.abs(depth) <= 4.6)
        rows = layer_rows(pos, chill, sheet, axis)
        if len(rows) < 2:
            return None
        return (min(row[1] for row in rows), rows, axis, view, sheet)

    seeds = [stacking_axis(pos, chill, mask)]
    best = None
    for theta in np.linspace(0.2, np.pi - 0.2, 14):
        for phi in np.linspace(0.0, np.pi, 14, endpoint=False):
            seeds.append(
                np.array([np.sin(theta) * np.cos(phi), np.sin(theta) * np.sin(phi), np.cos(theta)])
            )
    for seed in seeds:
        axis = unit(seed)
        hint = np.array([1.0, 0.0, 0.0]) if abs(float(axis[0])) < 0.9 else np.array([0.0, 1.0, 0.0])
        east = unit(np.cross(axis, hint))
        north = unit(np.cross(axis, east))
        for degree in np.linspace(0.0, 180.0, 24, endpoint=False):
            theta = math.radians(float(degree))
            view = math.cos(theta) * east + math.sin(theta) * north
            packed = pack(axis, view)
            if packed is not None and (best is None or packed[0] > best[0]):
                best = packed
    if best is None:
        axis = stacking_axis(pos, chill, mask)
        view = column_view(pos, lattice, axis)
        packed = pack(axis, view)
        best = packed if packed is not None else (0.0, [], unit(axis), unit(view), lattice)
    _score, rows, axis, view, sheet = best
    # Stay perpendicular to the c-axis. Prefer the in-plane direction that
    # keeps the sheet's layer order aligned with panel (c).
    tightened = column_view(pos, sheet if int(sheet.sum()) > 12 else lattice, axis)
    repacked = pack(axis, tightened)
    if repacked is not None and repacked[0] + 0.02 >= best[0]:
        _score, rows, axis, view, sheet = repacked
    else:
        # The purity-maximising view is already perpendicular to the c-axis,
        # so the layers of (c) are horizontal bands.
        view = view - axis * float(view @ axis)
        view = unit(view)
    if axis[2] < 0.0:
        axis = -axis
    print("nucleus c-axis", np.round(axis, 3), flush=True)
    for kind, purity, count in rows:
        print(f"  layer class {kind} purity {purity:.3f} n {count}", flush=True)
    weak = [row for row in rows if row[2] >= 4 and row[1] < 0.80]
    if weak:
        print(
            "  CHILL+ mixes both classes inside "
            f"{len(weak)} basal layer(s); purity stays the label, not a recolouring",
            flush=True,
        )
    return sheet, view, axis, rows


def nucleus_sheet(pos: np.ndarray, chill: np.ndarray, mask: np.ndarray, stacking: np.ndarray):
    """Two [11-20] layers through the centre of mass. Liquid and interfacial water are left out."""
    lattice = mask & np.isin(chill, [0, 1])
    if int(lattice.sum()) < 8:
        lattice = np.isin(chill, [0, 1])
    view = column_view(pos, lattice, stacking)
    centre = pos[lattice].mean(axis=0)
    depth = (pos - centre) @ view
    # The column repeat is the partner separation, about one lattice step.
    relative = pos[lattice] - centre
    stacking_u = unit(stacking)
    right = unit(np.cross(stacking_u, view))
    horizontal = relative @ right
    vertical = relative @ stacking_u
    along = relative @ view
    separations = []
    step = max(1, len(relative) // 80)
    for index in range(0, len(relative), step):
        delta_depth = np.abs(along - along[index])
        plane = np.hypot(horizontal - horizontal[index], vertical - vertical[index])
        near = (delta_depth > 2.2) & (delta_depth < 6.0) & (plane < 1.3)
        if np.any(near):
            separations.append(float(np.min(delta_depth[near])))
    period = float(np.median(separations)) if separations else 4.5
    period = float(np.clip(period, 3.6, 5.2))
    layer = np.round(depth / period)
    present, counts = np.unique(layer[lattice], return_counts=True)
    order = np.argsort(np.abs(present) * 1000 - counts)
    keep_layers = set(int(item) for item in present[order[:2]])
    sheet = lattice & np.isin(layer, list(keep_layers))
    print(
        f"nucleus sheet layers {sorted(keep_layers)} period {period:.2f} A "
        f"n {int(sheet.sum())} of {int(lattice.sum())}",
        flush=True,
    )
    return sheet, view


def display_layer_ids(pos: np.ndarray, mask: np.ndarray, axis: np.ndarray) -> np.ndarray:
    """Layer index of each atom, numbered as on the stacking bar."""
    ids = -np.ones(pos.shape[0], dtype=int)
    chosen = np.flatnonzero(mask)
    if chosen.size == 0:
        return ids
    projection = pos[chosen] @ unit(axis)
    period = layer_period(projection)
    raw = np.floor((projection - projection.min()) / period).astype(int)
    order = [index for index in range(int(raw.min()), int(raw.max()) + 1) if np.any(raw == index)]
    remap = {bin_id: number for number, bin_id in enumerate(order)}
    ids[chosen] = [remap[int(bin_id)] for bin_id in raw]
    return ids


def draw_layer_guides(png: Path, pos: np.ndarray, layer_ids: np.ndarray, forward, up, half_v: float) -> None:
    """Faint horizontal lines at each drawn layer, numbered as in panel (c)."""
    from PIL import Image, ImageDraw, ImageFont

    image = Image.open(png).convert("RGBA")
    forward, _right, up = basis(forward, up)
    centre, _h, _v, forward, _right, up = project(pos, forward, up)
    height = image.height
    draw = ImageDraw.Draw(image)
    # The nucleus image is already at print scale. 7 pt at 600 dpi is 58 px.
    font_px = max(28, int(round(7.0 * PRINT_DPI / 72.0)))
    try:
        font = ImageFont.truetype("/usr/share/fonts/truetype/dejavu/DejaVuSans.ttf", font_px)
    except OSError:
        font = ImageFont.load_default()
        font_px = 14
    present = sorted({int(item) for item in layer_ids if item >= 0})
    vertical = (pos - centre) @ up
    span = 2.0 * half_v
    for layer in present:
        members = np.flatnonzero(layer_ids == layer)
        if members.size == 0:
            continue
        height_world = float(np.median(vertical[members]))
        y = (half_v - height_world) / span * height
        if not 0 <= y < height:
            continue
        draw.line([(font_px + 8, y), (image.width - 8, y)], fill=(168, 168, 168, 255), width=2)
        draw.text((4, max(0, y - 0.55 * font_px)), str(layer), font=font, fill=(70, 70, 70, 255))
    image.convert("RGB").save(png)


def _paint(pos, chill, bonds, forward, up, png, half_h, half_v, width_px, liquid_alpha, nucleus: bool, frame_pos=None):
    from ovito.data import DataCollection, Particles, SimulationCell
    from ovito.pipeline import Pipeline, StaticSource
    from ovito.vis import BondsVis, TachyonRenderer, Viewport

    forward, _right, up = basis(forward, up)
    focus = pos if frame_pos is None else frame_pos
    if frame_pos is None and liquid_alpha < 0.3:
        solid = ~np.isin(chill, [2])
        if np.any(solid):
            focus = pos[solid]
    centre, _h, _v, forward, _right, up = project(focus, forward, up)
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
    if nucleus:
        colours, near = depth_cue(colours, pos @ forward, tint=np.array([0.94, 0.94, 0.95]))
    else:
        colours, near = depth_cue(colours, pos @ forward, 0.22)
    particles.create_property("Color", data=np.clip(colours, 0.0, 1.0))
    liquid = np.isin(chill, [2])
    guest = np.isin(chill, [4, 5])
    transparency = np.zeros(pos.shape[0])
    # Tachyon multiplies overlapping transparent spheres, so a crowd of
    # liquid becomes opaque. The caller composites one liquid layer instead.
    transparency[liquid] = 0.0
    transparency[guest] = 0.0
    radius = np.full(pos.shape[0], ATOM_RADIUS)
    if nucleus:
        radius *= 0.58 + 0.42 * near
    particles.create_property("Transparency", data=transparency)
    particles.create_property("Radius", data=radius)
    if bonds is not None and len(bonds):
        bond_data = particles.create_bonds()
        bond_data.create_property("Topology", data=np.ascontiguousarray(bonds, dtype=np.int32))
        bond_data.vis.width = 0.10 if nucleus else BOND_WIDTH
        bond_data.vis.color = (0.62, 0.62, 0.62) if nucleus else (0.10, 0.10, 0.10)
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


def render_molecule(pos, chill, bonds, forward, up, png, half_h, half_v, width_px, liquid_alpha, nucleus: bool = False):
    """Ray-trace the molecules.

    Liquid is painted on its own and blended at ``liquid_alpha``. Overlapping
    transparent spheres in one Tachyon pass multiply and turn the liquid solid.
    """
    from PIL import Image

    liquid = np.isin(chill, [2])
    if (not nucleus) and 0.0 <= liquid_alpha < 0.999 and np.any(liquid) and np.any(~liquid):
        solid = ~liquid
        ice_png = png.with_name(png.stem + "-ice.png")
        liq_png = png.with_name(png.stem + "-liq.png")
        solid_bonds = bonds
        if bonds is not None and len(bonds):
            keep = solid[bonds[:, 0]] & solid[bonds[:, 1]]
            solid_bonds = bonds[keep]
            remap = -np.ones(pos.shape[0], dtype=int)
            remap[np.flatnonzero(solid)] = np.arange(int(solid.sum()))
            if len(solid_bonds):
                solid_bonds = np.stack([remap[solid_bonds[:, 0]], remap[solid_bonds[:, 1]]], axis=1)
        _paint(pos[solid], chill[solid], solid_bonds, forward, up, ice_png, half_h, half_v, width_px, 1.0, False, frame_pos=pos)
        _paint(pos[liquid], chill[liquid], None, forward, up, liq_png, half_h, half_v, width_px, 1.0, False, frame_pos=pos)
        ice = np.asarray(Image.open(ice_png).convert("RGB"), dtype=np.float32)
        liq = np.asarray(Image.open(liq_png).convert("RGB"), dtype=np.float32)
        veil = np.any(liq < 250, axis=2)
        blended = ice.copy()
        blended[veil] = liquid_alpha * liq[veil] + (1.0 - liquid_alpha) * ice[veil]
        Image.fromarray(np.clip(blended, 0, 255).astype(np.uint8)).save(png)
        covered = blended[veil]
        print(
            f"liquid alpha {liquid_alpha:.2f} pixels {int(veil.sum())} mean {covered.mean():.1f}",
            flush=True,
        )
        return ice.shape[0]
    return _paint(pos, chill, bonds, forward, up, png, half_h, half_v, width_px, liquid_alpha, nucleus)


def compose(slab, nuclei, stacks, width_mm: float, png: Path, pdf: Path) -> None:
    import matplotlib

    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    from PIL import Image

    dpi = PRINT_DPI
    margin_l, margin_r, margin_t = 8.0, 3.5, 3.2
    legend_h, caption_h = 7.5, 0.0
    margin_b = legend_h + 1.2
    gap = 2.2
    narrow = width_mm < 150.0
    content_w = width_mm - margin_l - margin_r
    slab_w = content_w if narrow else 70.0
    right_w = content_w if narrow else content_w - slab_w - gap
    slab_h = slab_w * slab["half_v"] / slab["half_h"]
    angstrom_per_mm = (2.0 * slab["half_h"]) / slab_w

    def panel_mm(entry: dict) -> tuple[float, float]:
        image = Image.open(entry["png"])
        # The render used the slab scale. Keep that scale on the wide figure.
        # The narrow figure fits the column, so both nuclei still share one scale.
        width = entry["w_mm"]
        height = entry["h_mm"]
        if narrow or width > right_w:
            scale = right_w / width
            width *= scale
            height *= scale
        return width, height

    late_w, late_h = panel_mm(nuclei["late"])
    min_w, min_h = panel_mm(nuclei["minimum"])
    # On the wide figure the two sheets share the slab scale, so the 10 A
    # bars have one length. The narrow figure scales both by the same factor.
    if not narrow:
        scale = 1.0
    else:
        scale = late_w / nuclei["late"]["w_mm"]
    bar_mm = 10.0 / angstrom_per_mm * (scale if narrow else 1.0)
    title_h = 4.6
    pair_gap = 7.2
    bar_row = 7.2
    stack_h = 11.5
    right_h = title_h + late_h + pair_gap + min_h + bar_row + stack_h
    if narrow:
        height_mm = margin_t + slab_h + gap + right_h + margin_b
        slab_x, slab_y = margin_l, margin_b + right_h + gap
        right_x, right_top = margin_l, margin_b + right_h
    else:
        height_mm = margin_t + max(slab_h, right_h) + margin_b
        slab_x = margin_l
        slab_y = margin_b + max(0.0, right_h - slab_h)
        right_x = margin_l + slab_w + gap
        right_top = margin_b + right_h
    if (not narrow) and height_mm > 230.0:
        raise SystemExit(f"figure is {height_mm:.0f} mm tall, past a CPC page")
    print(f"figure {width_mm:.0f} x {height_mm:.1f} mm", flush=True)

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

    def scale_bar(x_mm, y_mm, vertical: bool, label: str):
        if vertical:
            xs = [x_mm / width_mm, x_mm / width_mm]
            ys = [y_mm / height_mm, (y_mm + bar_mm) / height_mm]
            text_x = (x_mm + 1.0) / width_mm
            text_y = (y_mm + 0.5 * bar_mm) / height_mm
            ha, va = "left", "center"
        else:
            xs = [x_mm / width_mm, (x_mm + bar_mm) / width_mm]
            ys = [y_mm / height_mm, y_mm / height_mm]
            text_x = (x_mm + bar_mm + 0.8) / width_mm
            text_y = y_mm / height_mm
            ha, va = "left", "center"
        fig.add_artist(
            plt.Line2D(xs, ys, transform=fig.transFigure, color="black", lw=1.1, solid_capstyle="butt")
        )
        fig.text(text_x, text_y, label, fontsize=7, ha=ha, va=va, linespacing=1.15)

    place(slab["png"], slab_x, slab_y, slab_w, slab_h)
    letter("(a)", slab_x, slab_y + slab_h + 0.2)
    scale_bar(1.1, slab_y + 0.5 * slab_h - 0.5 * bar_mm, True, "z\n10 Å")

    # (b) stacked above (c), both in the right column. One 10 A bar for (b).
    late_y = right_top - title_h - late_h
    min_y = late_y - pair_gap - min_h
    place(nuclei["late"]["png"], right_x, late_y, late_w, late_h)
    place(nuclei["minimum"]["png"], right_x, min_y, min_w, min_h)
    letter("(b)", right_x, right_top + 0.2)
    fig.text(
        (right_x + 0.5 * late_w) / width_mm,
        (late_y + late_h + 0.15) / height_mm,
        "late",
        fontsize=7,
        ha="center",
        va="bottom",
    )
    fig.text(
        (right_x + 0.5 * min_w) / width_mm,
        (min_y + min_h + 0.15) / height_mm,
        "minimum",
        fontsize=7,
        ha="center",
        va="bottom",
    )
    scale_bar(right_x, late_y - 2.8, False, "10 Å")
    scale_bar(right_x, min_y - 2.8, False, "10 Å")

    bar_h = stack_h
    bar_bottom = min_y - bar_row - bar_h
    letter("(c)", right_x, bar_bottom + bar_h + 0.4)
    label_w = 16.0
    axis = fig.add_axes(
        (
            (right_x + label_w) / width_mm,
            (bar_bottom + 1.2) / height_mm,
            (right_w - label_w) / width_mm,
            (bar_h - 3.2) / height_mm,
        )
    )
    for row, name in enumerate(("minimum", "late")):
        kinds = stacks[name]
        shown = int(stacks["shown"][name])
        for column, kind in enumerate(kinds):
            colour = CHILL.get(kind, ("", (0.5, 0.5, 0.5)))[1]
            axis.add_patch(plt.Rectangle((column, row), 1.0, 0.72, facecolor=colour, edgecolor="white", lw=0.3))
        if shown > 0:
            axis.add_patch(
                plt.Rectangle((0.0, row), shown, 0.72, fill=False, edgecolor="black", lw=0.9, zorder=3)
            )
        fraction = stacks["cubicity"][name]
        fig.text(
            (right_x + label_w - 0.6) / width_mm,
            (bar_bottom + 1.2 + (bar_h - 3.2) * (row + 0.36) / 1.9) / height_mm,
            f"{name}  {fraction:.2f}",
            fontsize=7,
            ha="right",
            va="center",
        )
    axis.set_xlim(0, max(len(stacks["minimum"]), len(stacks["late"]), 1))
    axis.set_ylim(-0.15, 1.9)
    axis.set_axis_off()
    axis.set_aspect("auto")
    fig.text(
        (right_x + right_w) / width_mm,
        (bar_bottom - 0.3) / height_mm,
        "layer index along [0001]",
        fontsize=7,
        ha="right",
        va="top",
    )

    handles = [
        plt.Line2D([0], [0], marker="o", color="none", markerfacecolor=colour, markeredgecolor="none", markersize=5, label=name)
        for _code, (name, colour) in CHILL.items()
    ]
    fig.legend(
        handles=handles,
        loc="lower center",
        ncol=3 if narrow else 6,
        frameon=False,
        fontsize=7,
        bbox_to_anchor=(0.5, (caption_h + 0.4) / height_mm),
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
    render_molecule(slab_pos, slab_chill, slab_bonds, forward, up, slab_png, half_h, half_v, slab_px, 0.40)
    slab = {"png": slab_png, "half_h": half_h, "half_v": half_v}

    n_frames = count_frames(nucleus_dump)
    if late_frame < 0:
        late_frame = n_frames + 1 + late_frame
    clouds = {}
    loaded = {}
    for label, frame in (("late", late_frame), ("minimum", min_frame)):
        nbox, npos = read_dump_frame(nucleus_dump, frame)
        nchill, extended = label_positions(nbox, npos, work / label, "nucleation")
        if extended is None:
            extended = np.isin(nchill, [0, 1, 3])
        raw = npos.copy()
        npos, _nbox = unwrap_cluster(npos, nbox, extended, HBOND_A)
        loaded[label] = (npos, nchill, extended, raw, nbox)
    # Each nucleus keeps its own c-axis. That axis is vertical and the
    # camera looks along the local [11-20], through two layers at the centre.
    prepared = {}
    orient_rows = {}
    for label, (npos, nchill, extended, raw, raw_box) in loaded.items():
        sheet, view, axis, orient_rows[label] = orient_nucleus(npos, nchill, extended)
        # Bonds only inside the outlined run of panel (c). Atoms with no
        # bond in that run are strays and are not drawn.
        kinds_now = layer_kinds(npos, nchill, extended, axis)
        sheet_kinds = [kind for kind, _purity, count in orient_rows[label] if count >= 4]
        shown_now = 0
        for full_kind, sheet_kind in zip(kinds_now, sheet_kinds):
            if full_kind != sheet_kind:
                break
            shown_now += 1
        outlined = outlined_layers(npos, extended & np.isin(nchill, [0, 1]), axis, shown_now)
        sheet = sheet & outlined
        bonds = ice_bonds(raw, raw_box, sheet, 3.5)
        if len(bonds):
            straight = np.linalg.norm(npos[bonds[:, 1]] - npos[bonds[:, 0]], axis=1)
            bonds = bonds[straight < 3.5]
        degree = np.bincount(bonds.ravel(), minlength=npos.shape[0]) if len(bonds) else np.zeros(npos.shape[0], dtype=int)
        if degree.size and int(degree.max()) > 4:
            raise SystemExit("a water molecule was given more than 4 bonds")
        linked = degree > 0
        sheet = sheet & linked
        if len(bonds):
            keep_bond = linked[bonds[:, 0]] & linked[bonds[:, 1]]
            bonds = bonds[keep_bond]
        drawn = npos[sheet]
        drawn_chill = nchill[sheet]
        # Number the guides with the stacking-bar index, not a fresh count
        # of whichever layers survived the sheet cut.
        layer_ids = display_layer_ids(npos, extended & np.isin(nchill, [0, 1]), axis)[sheet]
        remap = -np.ones(npos.shape[0], dtype=int)
        remap[np.flatnonzero(sheet)] = np.arange(int(sheet.sum()))
        if len(bonds):
            bonds = np.stack([remap[bonds[:, 0]], remap[bonds[:, 1]]], axis=1)
        print(f"  {label} drawn {int(sheet.sum())} bonds {len(bonds)}", flush=True)
        prepared[label] = (drawn, drawn_chill, bonds, view, axis, layer_ids)
        clouds[label] = (npos, nchill, extended, view, axis)
    def raw_half(label: str) -> tuple[float, float]:
        drawn, drawn_chill, _bonds, view_i, axis_i, _ids = prepared[label]
        cluster = drawn[np.isin(drawn_chill, [0, 1])]
        _forward_i, right_i, up_i = basis(view_i, axis_i)
        centred = cluster - cluster.mean(axis=0)
        half_h = float(np.max(np.abs(centred @ right_i))) + ATOM_RADIUS
        half_v = float(np.max(np.abs(centred @ up_i))) + ATOM_RADIUS
        return half_h, half_v

    # Both sheets use the slab's angstrom scale, so one 10 A bar serves (b).
    angstrom_per_mm = (2.0 * slab["half_h"]) / 70.0
    angstrom_per_px = angstrom_per_mm / (PRINT_DPI / 25.4)
    images = {}
    stacks = {"shown": {}, "cubicity": {}}
    for label, (drawn, drawn_chill, bonds, view, axis, layer_ids) in prepared.items():
        images[label] = work / f"nucleus-{label}.png"
        half_h, half_v = raw_half(label)
        width_px = max(32, int(round(2.0 * half_h / angstrom_per_px)))
        render_molecule(
            drawn, drawn_chill, bonds, view, axis, images[label], half_h, half_v, width_px, 0.18, nucleus=True
        )
        draw_layer_guides(images[label], drawn, layer_ids, view, axis, half_v)
        print(f"  {label} guides {sorted({int(i) for i in layer_ids if i >= 0})}", flush=True)
        npos, nchill, extended, _view, axis = clouds[label]
        kinds = layer_kinds(npos, nchill, extended, axis)
        stacks[label] = kinds
        stacks["cubicity"][label] = cubicity(kinds)
        sheet_kinds = [kind for kind, _purity, count in orient_rows[label] if count >= 4]
        shown = 0
        for full_kind, sheet_kind in zip(kinds, sheet_kinds):
            if full_kind != sheet_kind:
                break
            shown += 1
        stacks["shown"][label] = shown
    shared = stacking_axis(loaded["late"][0], loaded["late"][1], loaded["late"][2])
    stacks["cubicity_shared"] = {
        label: cubicity(layer_kinds(npos, nchill, extended, shared))
        for label, (npos, nchill, extended, _view, _axis) in clouds.items()
    }
    px_per_mm = PRINT_DPI / 25.4
    nuclei = {}
    for label, path in images.items():
        from PIL import Image

        image = Image.open(path)
        nuclei[label] = {"png": path, "w_mm": image.width / px_per_mm, "h_mm": image.height / px_per_mm}
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
