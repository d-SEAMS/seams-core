#!/usr/bin/env python3
# Copyright (c) d-SEAMS developers
"""Draw molecule stills.

chemparseplot is the entry point for the oxygen-oxygen bond list.
It does not path-trace, it does not colour atoms by a classifier,
and it does not draw translucent faces. The picture is Blender
Cycles when ``bpy`` imports: glossy spheres, cylinder bonds, glass
cage faces, a three-point studio and a linear Standard view.
Tachyon is the fallback, with antialiasing of 64 samples.

    python docs/scripts/render_molecules.py --out <dir> [--slug <slug>]

Each still is rendered above 2048 px and stored at 2048 by 1536,
plus a web PNG and a dark pair.
The files are build products and are not stored in git.
"""

from __future__ import annotations

import argparse
import os
import sys
import warnings
from pathlib import Path

ROOT = Path(__file__).resolve().parents[2]


def _loader_path():
    """Keep the system OpenGL libraries visible to OVITO's Qt build."""
    folders = []
    for candidate in ("/usr/lib/x86_64-linux-gnu", "/usr/lib64", "/usr/lib"):
        if os.path.isdir(candidate):
            folders.append(candidate)
    current = os.environ.get("LD_LIBRARY_PATH", "")
    if current:
        folders.append(current)
    if folders:
        os.environ["LD_LIBRARY_PATH"] = os.pathsep.join(folders)


_loader_path()
WIDTH = 2048
HEIGHT = 1536
RENDER_SCALE = max(1, int(os.environ.get("MOLECULES_SCALE", "2")))
WEB = (960, 720)

# Ruhi palette. Every slug reads colours from this dict and from nowhere else.
PALETTE = {
    "hexagonal": "#1E88E5",
    "cubic": "#FF655D",
    "interfacial": "#F1DB4B",
    "clathrate": "#D81B60",
    "interClathrate": "#004D40",
    "water": "#B0BEC5",
    "unclassified": "#78909C",
    "reCubic": "#FF655D",
    "reHex": "#1E88E5",
    "sodium": "#AA00FF",
    "chloride": "#00C853",
    "methane": "#F9A825",
    "THF": "#EF6C00",
    "hydrogen": "#90A4AE",
    # ColorBrewer Dark2 purple. It stays clear of the interfacial yellow.
    "mixed": "#7570B3",
}
CHILL_NAME = {
    0: "cubic",
    1: "hexagonal",
    2: "water",
    3: "interfacial",
    4: "clathrate",
    5: "interClathrate",
    6: "unclassified",
    7: "reCubic",
    8: "reHex",
}


def palette_rgb(token):
    text = PALETTE[token].lstrip("#")
    return tuple(int(text[i:i + 2], 16) / 255.0 for i in (0, 2, 4))


def srgb_to_linear(colour):
    """sRGB display colour to the linear values Blender materials expect."""

    def channel(value):
        value = float(value)
        if value <= 0.04045:
            return value / 12.92
        return ((value + 0.055) / 1.055) ** 2.4

    return tuple(channel(value) for value in colour)


def _lab(colour):
    """sRGB 0-1 tuple to CIE Lab (D65)."""
    linear = srgb_to_linear(colour)
    matrix = (
        (0.4124564, 0.3575761, 0.1804375),
        (0.2126729, 0.7151522, 0.0721750),
        (0.0193339, 0.1191920, 0.9503041),
    )
    xyz = [sum(row[i] * linear[i] for i in range(3)) for row in matrix]
    white = (0.95047, 1.0, 1.08883)

    def f_lab(value):
        if value > 0.008856:
            return value ** (1.0 / 3.0)
        return 7.787 * value + 16.0 / 116.0

    f = [f_lab(xyz[i] / white[i]) for i in range(3)]
    return (116.0 * f[1] - 16.0, 500.0 * (f[0] - f[1]), 200.0 * (f[1] - f[2]))


def delta_e(left, right):
    """CIE76 distance between two sRGB colours in 0-1."""
    lab_left = _lab(left)
    lab_right = _lab(right)
    return sum((lab_left[i] - lab_right[i]) ** 2 for i in range(3)) ** 0.5


def atom_style(token):
    """Colour, radius and transparency. Oxygen is a 0.6 angstrom sphere."""
    colour = palette_rgb(token)
    if token == "water":
        return colour, 0.28, 0.78
    if token == "hydrogen":
        return colour, 0.22, 0.0
    if token == "sodium":
        return colour, 1.05, 0.0
    if token == "chloride":
        return colour, 1.30, 0.0
    if token in ("methane", "THF"):
        return colour, 0.56, 0.0
    return colour, 0.60, 0.0

# Video slug -> stills, in order. The first file is overlay.png.
VIDEO_SLUGS = {
    "c": ("ice-ih", "ice-ic"),
    "cpp": ("ice-ih", "ice-ic"),
    "python": ("ice-ih", "ice-ic"),
    "fortran": ("ice-ih", "ice-ic"),
    "rust": ("ice-ih", "ice-ic"),
    "julia": ("ice-ih", "ice-ic"),
    "r": ("ice-ih", "ice-ic"),
    "cli": ("ice-ih", "ice-ic"),
    "lua": ("ice-ih", "ice-ic"),
    "migrate": ("ice-ih",),
    "lammps": ("ice-ih", "ice-ic"),
    "gromacs": ("ice-ih",),
    "plumed": ("ice-ih",),
    "insitu": ("ice-ih", "ice-ic"),
    "interface": ("interface",),
    "rings-cages": ("rings", "prisms", "cage-si", "cage-sii", "cage-sh", "cage-guest"),
    "brine": ("brine",),
    "nucleation": ("nucleation",),
}


def fail(message):
    print(message, file=sys.stderr)
    raise SystemExit(1)


def stack_ok():
    try:
        import ase  # noqa: F401
        import chemparseplot  # noqa: F401
        import ovito  # noqa: F401

        warnings.filterwarnings("ignore", message=".*OVITO.*PyPI")
        return True
    except Exception:
        return False


def relaunch(argv):
    """Run this file under another interpreter. argv is the script's own arguments."""
    import subprocess

    if os.environ.get("MOLECULES_PY") == "1":
        fail("install chemparseplot, ovito and ase")
    for folder in os.environ.get("PATH", "").split(os.pathsep):
        program = os.path.join(folder, "python3")
        if not os.path.isfile(program):
            continue
        try:
            if os.path.samefile(program, sys.executable):
                continue
        except OSError:
            continue
        probe = subprocess.run(
            [program, "-c", "import chemparseplot, ovito, ase"],
            check=False,
            capture_output=True,
        )
        if probe.returncode != 0:
            continue
        env = os.environ.copy()
        env["MOLECULES_PY"] = "1"
        done = subprocess.run([program, str(Path(__file__).resolve()), *argv], env=env)
        return done.returncode
    fail("install chemparseplot, ovito and ase")


def load_renderer():
    """Return chemparseplot's structure renderer.

    chemparseplot 1.9 evaluates a forward reference while the class body
    is still running, so the module is loaded with postponed annotations
    when the plain import fails.
    """
    try:
        from chemparseplot.plot.neb import _render_atoms

        return _render_atoms
    except Exception:
        import chemparseplot

        path = Path(chemparseplot.__file__).resolve().parent / "plot" / "neb.py"
        source = "from __future__ import annotations\n" + path.read_text()
        module = type(sys)("chemparseplot.plot.neb")
        module.__file__ = str(path)
        module.__package__ = "chemparseplot.plot"
        sys.modules["chemparseplot.plot.neb"] = module
        exec(compile(source, str(path), "exec"), module.__dict__)
        return module._render_atoms


def find_lib(explicit):
    if explicit:
        return Path(explicit)
    env = os.environ.get("SEAMS_LIBDIR")
    if env:
        return Path(env)
    for relative in ("bbdir/src/libseams.so", "build/src/libseams.so"):
        candidate = ROOT / relative
        if candidate.is_file():
            return candidate.parent
    fail("set SEAMS_LIBDIR to the directory that holds libseams.so")


class Frame:
    def __init__(self, xyz, types, box9, symbols=None):
        import numpy as np

        self.xyz = np.asarray(xyz, dtype=float).reshape(-1, 3)
        self.types = np.asarray(types, dtype=int)
        self.box9 = list(box9)
        self.symbols = symbols

    def subset(self, index):
        import numpy as np

        index = np.asarray(index, dtype=int)
        symbols = None if self.symbols is None else [self.symbols[int(i)] for i in index]
        return Frame(self.xyz[index], self.types[index], self.box9, symbols)


def read_dump(path, frame_index=0):
    lines = Path(path).read_text().splitlines()
    frames = []
    i = 0
    while i < len(lines):
        if not lines[i].startswith("ITEM: NUMBER OF ATOMS"):
            i += 1
            continue
        count = int(lines[i + 1])
        i += 2
        if not lines[i].startswith("ITEM: BOX BOUNDS"):
            fail("box bounds missing in {}".format(path))
        tilted = "xy" in lines[i]
        bounds = []
        for _ in range(3):
            i += 1
            bounds.append([float(item) for item in lines[i].split()])
        i += 1
        if tilted:
            xy, xz, yz = bounds[0][2], bounds[1][2], bounds[2][2]
            xlo, xhi = bounds[0][0], bounds[0][1]
            ylo, yhi = bounds[1][0], bounds[1][1]
            zlo, zhi = bounds[2][0], bounds[2][1]
        else:
            xy = xz = yz = 0.0
            (xlo, xhi), (ylo, yhi), (zlo, zhi) = bounds
        box9 = [xlo, ylo, zlo, xhi - xlo, yhi - ylo, zhi - zlo, xy, xz, yz]
        if not lines[i].startswith("ITEM: ATOMS"):
            fail("atom header missing in {}".format(path))
        columns = lines[i].split()[2:]
        rows = []
        for _ in range(count):
            i += 1
            rows.append(dict(zip(columns, lines[i].split())))
        i += 1
        xyz = [[float(row["x"]), float(row["y"]), float(row["z"])] for row in rows]
        types = [int(row.get("type", "1")) for row in rows]
        frames.append(Frame(xyz, types, box9))
        if len(frames) > frame_index:
            break
    if frame_index >= len(frames):
        fail("{} has no frame {}".format(path, frame_index))
    return frames[frame_index]


def load_api(libdir):
    from ctypes import CDLL, POINTER, c_double, c_int, c_size_t, c_void_p

    previous = os.environ.get("LD_LIBRARY_PATH", "")
    os.environ["LD_LIBRARY_PATH"] = os.pathsep.join([str(libdir), previous])
    library = CDLL(str(Path(libdir) / "libseams.so"))
    library.seams_system_create.restype = c_void_p
    library.seams_system_destroy.argtypes = [c_void_p]
    library.seams_system_set_frame.argtypes = [
        c_void_p, POINTER(c_double), c_int, POINTER(c_double), POINTER(c_int), POINTER(c_int),
    ]
    library.seams_system_set_frame.restype = c_int
    library.seams_chill_plus.argtypes = [
        c_void_p, c_double, c_int, c_int, c_int, POINTER(c_int), POINTER(c_int),
    ]
    library.seams_chill_plus.restype = c_int
    library.seams_rings_find.argtypes = [c_void_p, c_double, c_int, c_int, c_int, POINTER(c_void_p)]
    library.seams_rings_find.restype = c_int
    library.seams_rings_count.argtypes = [c_void_p]
    library.seams_rings_count.restype = c_int
    library.seams_rings_vertices.argtypes = [
        c_void_p, c_int, POINTER(c_int), c_size_t, POINTER(c_size_t),
    ]
    library.seams_rings_vertices.restype = c_int
    library.seams_rings_destroy.argtypes = [c_void_p]
    return library


def _as_atoms(frame, symbols):
    from ase import Atoms

    return Atoms(symbols=symbols, positions=frame.xyz.tolist())


def chill_labels(library, frame):
    from ctypes import POINTER, c_double, c_int

    count = len(frame.xyz)
    xyz = (c_double * (count * 3))(*frame.xyz.reshape(-1).tolist())
    box = (c_double * 9)(*frame.box9)
    types = (c_int * count)(*frame.types.tolist())
    system = library.seams_system_create()
    if library.seams_system_set_frame(system, xyz, count, box, types, None) != 0:
        library.seams_system_destroy(system)
        fail("seams_system_set_frame failed")
    labels = (c_int * count)()
    counts = (c_int * 9)()
    if library.seams_chill_plus(system, 3.5, 0, 0, 0, labels, counts) != 0:
        library.seams_system_destroy(system)
        fail("seams_chill_plus failed")
    library.seams_system_destroy(system)
    return [int(labels[i]) for i in range(count)]


def ring_members(library, frame, wanted):
    from ctypes import POINTER, byref, c_double, c_int, c_size_t, c_void_p

    count = len(frame.xyz)
    xyz = (c_double * (count * 3))(*frame.xyz.reshape(-1).tolist())
    box = (c_double * 9)(*frame.box9)
    types = (c_int * count)(*frame.types.tolist())
    system = library.seams_system_create()
    if library.seams_system_set_frame(system, xyz, count, box, types, None) != 0:
        library.seams_system_destroy(system)
        fail("seams_system_set_frame failed")
    rings = c_void_p()
    if library.seams_rings_find(system, 3.5, 0, 0, 6, byref(rings)) != 0:
        library.seams_system_destroy(system)
        fail("seams_rings_find failed")
    chosen = set()
    total = library.seams_rings_count(rings)
    buffer = (c_int * 16)()
    written = c_size_t()
    for index in range(total):
        if library.seams_rings_vertices(rings, index, buffer, 16, written) != 0:
            continue
        size = int(written.value)
        if size not in wanted:
            continue
        for slot in range(size):
            chosen.add(int(buffer[slot]))
    library.seams_rings_destroy(rings)
    library.seams_system_destroy(system)
    return chosen


def central_slice(frame, labels, limit):
    import numpy as np

    if len(frame.xyz) <= limit:
        return frame, labels
    ice = [i for i, label in enumerate(labels) if label in (0, 1, 7, 8)]
    anchor = ice if ice else list(range(len(labels)))
    centre = frame.xyz[anchor].mean(axis=0)
    distance = np.linalg.norm(frame.xyz - centre, axis=1)
    keep = [int(i) for i in np.argsort(distance)[:limit]]
    keep.sort()
    return frame.subset(keep), [labels[i] for i in keep]


def legend_rows(names, colors):
    rows = []
    seen = set()
    for name, rgb in zip(names, colors):
        if not name or str(name).startswith("_") or name in seen:
            continue
        seen.add(name)
        rows.append((name, tuple(rgb)))
    return rows


def oxygen_bonds(positions, heavy):
    """Oxygen-oxygen contacts between 2.30 and 3.30 angstroms.

    chemparseplot's pair list is used when that module is present.
    """
    import numpy as np

    symbols = ["O" if flag else "H" for flag in heavy]
    edges = None
    try:
        from chemparseplot.plot.rings import bonds_within_cutoff

        edges = bonds_within_cutoff(symbols, positions, 3.30)
    except Exception:
        edges = None
    if edges is None:
        xyz = np.asarray(positions, dtype=float)
        index = [i for i, flag in enumerate(heavy) if flag]
        edges = []
        for cursor, left in enumerate(index):
            for right in index[cursor + 1 :]:
                gap = float(np.linalg.norm(xyz[left] - xyz[right]))
                if 2.30 <= gap <= 3.30:
                    edges.append((left, right))
        return edges
    kept = []
    xyz = np.asarray(positions, dtype=float)
    for left, right in edges:
        gap = float(np.linalg.norm(xyz[left] - xyz[right]))
        if gap >= 2.30:
            kept.append((left, right))
    return kept


def place_guests(positions, box9, limit=24):
    """Cage centres: grid points far from every oxygen, then separated."""
    import numpy as np

    lengths = np.asarray(box9[3:6], dtype=float)
    origin = np.asarray(box9[:3], dtype=float)
    local = np.asarray(positions, dtype=float) - origin
    step = 0.7
    axes = [np.arange(step / 2.0, length, step) for length in lengths]
    grid = np.stack(np.meshgrid(*axes, indexing="ij"), axis=-1).reshape(-1, 3)
    nearest = np.full(len(grid), 1e9)
    block = 4000
    for start in range(0, len(grid), block):
        piece = grid[start : start + block]
        delta = piece[:, None, :] - local[None, :, :]
        delta -= lengths * np.round(delta / lengths)
        nearest[start : start + block] = np.linalg.norm(delta, axis=2).min(axis=1)
    picked = []
    for index in np.argsort(-nearest):
        if nearest[index] < 3.05:
            break
        point = grid[index]
        if picked:
            delta = point - np.asarray(picked)
            delta -= lengths * np.round(delta / lengths)
            if np.linalg.norm(delta, axis=1).min() < 3.4:
                continue
        local_point = point
        if np.any(local_point < 1.6) or np.any(lengths - local_point < 1.6):
            continue
        picked.append(local_point)
        if len(picked) >= limit:
            break
    if not picked:
        return np.zeros((0, 3))
    return np.asarray(picked) + origin


def _spin(pos, axis, degrees):
    import numpy as np

    angle = np.deg2rad(degrees)
    cosine, sine = np.cos(angle), np.sin(angle)
    out = np.array(pos, copy=True)
    if axis == "x":
        y_coord, z_coord = out[:, 1].copy(), out[:, 2].copy()
        out[:, 1] = cosine * y_coord - sine * z_coord
        out[:, 2] = sine * y_coord + cosine * z_coord
    else:
        x_coord, z_coord = out[:, 0].copy(), out[:, 2].copy()
        out[:, 0] = cosine * x_coord + sine * z_coord
        out[:, 2] = -sine * x_coord + cosine * z_coord
    return out


def _basal_flat(points, center):
    """Basal coordinates before the small tilt. The c axis lies on z."""
    import numpy as np

    pos = np.asarray(points, dtype=float) - np.asarray(center, dtype=float)
    return np.column_stack((pos[:, 0], pos[:, 2], pos[:, 1]))


def _seat(canvas, light):
    """Shift the ink so its centroid sits on the middle of the frame."""
    import numpy as np
    from PIL import Image

    image = canvas.convert("RGB")
    arr = np.asarray(image)
    if light:
        mask = np.any(arr < 248, axis=2)
    else:
        mask = np.any(arr > 40, axis=2)
    ys, xs = np.where(mask)
    if len(xs) == 0:
        return image
    height, width = mask.shape
    dx = float(xs.mean() / width - 0.5)
    dy = float(ys.mean() / height - 0.5)
    if abs(dx) <= 0.02 and abs(dy) <= 0.02:
        return image
    pad = 2
    left = max(0, int(xs.min()) - pad)
    right = min(width, int(xs.max()) + 1 + pad)
    top = max(0, int(ys.min()) - pad)
    bottom = min(height, int(ys.max()) + 1 + pad)
    cx = float(xs.mean() - left)
    cy = float(ys.mean() - top)
    crop = image.crop((left, top, right, bottom))
    room_x = width / 2.0 - 8.0
    room_y = height / 2.0 - 8.0
    scale = min(
        room_x / max(cx, crop.size[0] - cx),
        room_y / max(cy, crop.size[1] - cy),
        1.0,
    )
    if scale < 0.999:
        sized = (
            max(1, int(round(crop.size[0] * scale))),
            max(1, int(round(crop.size[1] * scale))),
        )
        crop = crop.resize(sized, Image.Resampling.LANCZOS)
        cx *= scale
        cy *= scale
    ox = int(round(width / 2.0 - cx))
    oy = int(round(height / 2.0 - cy))
    ox = max(0, min(width - crop.size[0], ox))
    oy = max(0, min(height - crop.size[1], oy))
    out = Image.new("RGB", (width, height), image.getpixel((1, 1)))
    out.paste(crop, (ox, oy))
    return out


def _principal_basis(positions, mask):
    """Centre and axes. Columns are major, middle, then the short axis."""
    import numpy as np

    chosen = np.asarray(positions, dtype=float)[mask]
    centre = chosen.mean(axis=0)
    delta = chosen - centre
    _values, vectors = np.linalg.eigh(np.cov(delta.T))
    order = np.argsort(_values)[::-1]
    basis = vectors[:, order].copy()
    if np.linalg.det(basis) < 0:
        basis[:, 2] *= -1
    return centre, basis


def _apply_basis(points, centre, basis):
    import numpy as np

    return (np.asarray(points, dtype=float) - centre) @ basis


def _basal_fog(positions, colors, transparency, radii):
    """Fade stacked ice layers so the back bilayers recede.

    The front bilayer stays solid, which keeps the chair hexagon.
    Each bilayer behind it shrinks and fades toward the background.
    """
    import numpy as np

    z = np.asarray(positions, dtype=float)[:, 2]
    colors = np.asarray(colors, dtype=float)
    transparency = np.asarray(transparency, dtype=float)
    radii = np.asarray(radii, dtype=float)
    layers = []
    for idx in np.argsort(z):
        if not layers or float(z[idx] - z[layers[-1][0]]) > 0.55:
            layers.append([int(idx)])
        else:
            layers[-1].append(int(idx))
    # The camera looks down minus z, so the last cluster is nearest.
    for rank, members in enumerate(reversed(layers)):
        bilayer = rank // 2
        rear = rank % 2 == 1
        fog = min(0.42, 0.16 * bilayer + (0.04 if rear else 0.0))
        if fog <= 0.0:
            continue
        idx = np.asarray(members, dtype=int)
        # Stay blue. Mixing toward white is what blew the back spheres out.
        colors[idx] = colors[idx] * (1.0 - 0.18 * fog)
        transparency[idx] = np.clip(transparency[idx] + 0.22 * fog, 0.0, 0.28)
        radii[idx] *= 1.0 - 0.06 * fog
    return colors, transparency, radii


def orient(positions, view, center):
    """Place the structure for the camera, which looks down minus z.

    In the ice Ih dump the c axis is y: the layers alternate a short
    gap and a long one. Basal view looks down that axis, four degrees
    off, so the chair rings stay readable.
    """
    import numpy as np

    pos = np.asarray(positions, dtype=float) - np.asarray(center, dtype=float)
    if view == "front":
        pos = np.column_stack((pos[:, 2], pos[:, 1], -pos[:, 0]))
        return _spin(_spin(pos, "y", 18.0), "x", 24.0)
    if view == "cage":
        return _spin(_spin(pos, "y", 34.0), "x", 26.0)
    if view == "basal":
        pos = np.column_stack((pos[:, 0], pos[:, 2], pos[:, 1]))
        return _spin(pos, "x", 4.0)
    if view == "channel":
        # The ice c axis is y. The long cell axis z runs across the page.
        return np.column_stack((pos[:, 2], pos[:, 0], pos[:, 1]))
    return pos


CYCLES_SAMPLES = 256
RENDER_ENGINE = "cycles"
PANEL_NAME = "panel"


def cycles_python():
    """A Python that can import bpy. None when Cycles is not installed."""
    import subprocess

    programs = []
    chosen = os.environ.get("MOLECULES_BPY")
    if chosen:
        programs.append(chosen)
    programs.append(sys.executable)
    for folder in os.environ.get("PATH", "").split(os.pathsep):
        for name in ("python3.13", "python3"):
            programs.append(os.path.join(folder, name))
    seen = set()
    for program in programs:
        if not program or program in seen or not os.path.isfile(program):
            continue
        seen.add(program)
        probe = subprocess.run(
            [program, "-c", "import bpy"],
            check=False,
            capture_output=True,
        )
        if probe.returncode == 0:
            return program
    return None


def frame_bounds(positions, radii, transparency):
    """World box of every drawn atom, plus one radius so spheres stay inside."""
    import numpy as np

    points = np.asarray(positions, dtype=float)
    radius = np.asarray(radii, dtype=float)
    pad = float(radius.max()) + 0.35 if len(radius) else 0.6
    low = points.min(axis=0) - pad
    high = points.max(axis=0) + pad
    return low, high


def cycles_job():
    """Headless Cycles. The parent passes one numpy archive on argv."""
    import numpy as np
    import bpy

    job = np.load(sys.argv[1], allow_pickle=False)
    positions = np.asarray(job["positions"], dtype=float)
    colors = np.asarray(job["colors"], dtype=float)
    radii = np.asarray(job["radii"], dtype=float)
    transparency = np.asarray(job["transparency"], dtype=float)
    bonds = np.asarray(job["bonds"], dtype=int)
    background = tuple(float(channel) for channel in job["background"])
    mesh_count = int(job["mesh_count"])
    light = background[0] > 0.5
    framed = [positions]
    for slot in range(mesh_count):
        key = "mesh_vertices_{}".format(slot)
        if key in job.files:
            framed.append(np.asarray(job[key], dtype=float))
    framed = np.vstack(framed)
    com = framed.mean(axis=0)
    pad = float(radii.max()) + 1.1 if len(radii) else 1.1
    low = framed.min(axis=0) - pad
    high = framed.max(axis=0) + pad
    half = np.maximum(com - low, high - com)
    # Look at the atom centroid and fit the farther side, so the mass sits
    # in the middle and nothing is clipped.
    centre = np.array([float(com[0]), float(com[1]), float((low[2] + high[2]) / 2.0)])
    width = float(2.0 * half[0])
    height = float(2.0 * half[1])
    look = centre.copy()
    vertical = max(height, width * HEIGHT / float(WIDTH))

    for obj in list(bpy.data.objects):
        bpy.data.objects.remove(obj, do_unlink=True)

    def socket(node, *names):
        for name in names:
            if name in node.inputs:
                return node.inputs[name]
        return None

    def glossy(name, rgb, alpha, transmission):
        """Palette colour is sRGB, stored as a linear Principled base.

        Opaque spheres get ambient occlusion at 2 angstroms, a soft coat
        and a little emission so the centre pixel stays near the hex.
        """
        material = bpy.data.materials.new(name)
        material.use_nodes = True
        tree = material.node_tree
        tree.nodes.clear()
        output = tree.nodes.new("ShaderNodeOutputMaterial")
        linear = srgb_to_linear(rgb)
        principled = tree.nodes.new("ShaderNodeBsdfPrincipled")
        colour_in = socket(principled, "Base Color")
        if colour_in is not None:
            colour_in.default_value = (*linear, 1.0)
        rough = socket(principled, "Roughness")
        if rough is not None:
            rough.default_value = 0.05 if transmission > 0.5 else 0.32
        spec = socket(principled, "Specular IOR Level", "Specular")
        if spec is not None:
            spec.default_value = 0.5
        coat = socket(principled, "Coat Weight")
        if coat is not None and transmission < 0.5:
            coat.default_value = 0.25
        trans = socket(principled, "Transmission Weight", "Transmission")
        if trans is not None:
            trans.default_value = float(transmission)
        ior = socket(principled, "IOR")
        if ior is not None and transmission > 0.5:
            ior.default_value = 1.45
        fade = socket(principled, "Alpha")
        if fade is not None:
            fade.default_value = float(alpha)
        if transmission > 0.5 or alpha < 0.85:
            tree.links.new(principled.outputs["BSDF"], output.inputs["Surface"])
            return material
        ao = tree.nodes.new("ShaderNodeAmbientOcclusion")
        ao.inputs["Distance"].default_value = 2.0
        ao.inputs["Color"].default_value = (1.0, 1.0, 1.0, 1.0)
        tint = tree.nodes.new("ShaderNodeMixRGB")
        tint.blend_type = "MULTIPLY"
        tint.inputs["Fac"].default_value = 0.75
        tint.inputs["Color1"].default_value = (*linear, 1.0)
        tree.links.new(ao.outputs["Color"], tint.inputs["Color2"])
        tree.links.new(tint.outputs["Color"], colour_in)
        if rough is not None:
            rough.default_value = 0.46
        if spec is not None:
            spec.default_value = 0.16
        if coat is not None:
            coat.default_value = 0.05
        emission = tree.nodes.new("ShaderNodeEmission")
        emission.inputs["Color"].default_value = (*linear, 1.0)
        emission.inputs["Strength"].default_value = 0.18
        added = tree.nodes.new("ShaderNodeAddShader")
        tree.links.new(principled.outputs["BSDF"], added.inputs[0])
        tree.links.new(emission.outputs["Emission"], added.inputs[1])
        tree.links.new(added.outputs["Shader"], output.inputs["Surface"])
        return material

    def add_mesh(name, mesh, material):
        obj = bpy.data.objects.new(name, mesh)
        bpy.context.scene.collection.objects.link(obj)
        obj.data.materials.append(material)
        return obj

    def sphere_proto(name, radius):
        import bmesh

        mesh = bpy.data.meshes.new(name)
        pen = bmesh.new()
        bmesh.ops.create_uvsphere(
            pen, u_segments=24, v_segments=16, radius=float(radius),
        )
        pen.to_mesh(mesh)
        pen.free()
        for poly in mesh.polygons:
            poly.use_smooth = True
        return mesh

    def bond_proto(name):
        import bmesh

        mesh = bpy.data.meshes.new(name)
        pen = bmesh.new()
        bmesh.ops.create_cone(
            pen, cap_ends=True, segments=10,
            radius1=1.0, radius2=1.0, depth=1.0,
        )
        pen.to_mesh(mesh)
        pen.free()
        return mesh

    def link_instances(name, mesh, material, matrices):
        mesh.materials.append(material)
        for matrix in matrices:
            obj = bpy.data.objects.new(name, mesh)
            obj.matrix_world = matrix
            bpy.context.scene.collection.objects.link(obj)

    groups = {}
    for index, colour in enumerate(colors):
        key = tuple(round(float(channel), 4) for channel in colour)
        groups.setdefault(key, {"index": [], "alpha": [], "radius": []})
        groups[key]["index"].append(index)
        groups[key]["alpha"].append(float(1.0 - transparency[index]))
        groups[key]["radius"].append(float(radii[index]))
    for key, group in groups.items():
        alpha = min(group["alpha"])
        # Faint context (the liquid shell) stays a flat ghost. Mid alpha is glass.
        transmission = 0.55 if 0.22 < alpha < 0.75 else 0.0
        material = glossy("atom", key, alpha, transmission)
        radius = float(np.median(group["radius"]))
        from mathutils import Matrix

        centres = positions[group["index"]]
        matrices = [
            Matrix.Translation(tuple(float(v) for v in centre))
            for centre in centres
        ]
        link_instances("spheres", sphere_proto("spheres", radius), material, matrices)
    def add_bond_set(pairs, radius, name):
        if len(pairs) == 0:
            return
        from mathutils import Matrix, Vector

        ink = tuple(float(v) for v in job["bond_rgb"]) if "bond_rgb" in job.files else ((0.42, 0.44, 0.47) if light else (0.78, 0.80, 0.82))
        material = glossy(name, ink, 1.0, 0.0)
        rough = material.node_tree.nodes.get("Principled BSDF")
        if rough is not None:
            bond_rough = socket(rough, "Roughness")
            if bond_rough is not None:
                bond_rough.default_value = 0.45
        up = Vector((0.0, 0.0, 1.0))
        matrices = []
        for left, right in pairs:
            start = positions[int(left)]
            stop = positions[int(right)]
            span = stop - start
            length = float(np.linalg.norm(span))
            if length < 1e-6:
                continue
            direction = span / length
            gap = float(radii[int(left)] + radii[int(right)])
            visible = length - gap * 0.92
            if visible < 0.08:
                continue
            mid = (start + stop) / 2.0
            quat = up.rotation_difference(Vector(tuple(float(v) for v in direction)))
            scale = Matrix.Diagonal((float(radius), float(radius), float(visible), 1.0))
            matrices.append(Matrix.Translation(tuple(float(v) for v in mid)) @ quat.to_matrix().to_4x4() @ scale)
        if matrices:
            link_instances(name, bond_proto(name), material, matrices)

    bond_radius = float(job["bond_radius"]) if "bond_radius" in job.files else 0.045
    bold = np.asarray(job["bold_bonds"], dtype=int) if "bold_bonds" in job.files else np.zeros((0, 2), dtype=int)
    if len(bonds) and len(bold):
        bold_key = {tuple(sorted((int(pair[0]), int(pair[1])))) for pair in bold}
        thin_pairs = [
            pair for pair in bonds
            if tuple(sorted((int(pair[0]), int(pair[1])))) not in bold_key
        ]
        add_bond_set(thin_pairs, bond_radius, "bonds")
        add_bond_set(bold, bond_radius * 1.2, "cage-bonds")
    elif len(bonds):
        add_bond_set(bonds, bond_radius, "bonds")
    for slot in range(mesh_count):
        vertices = np.asarray(job["mesh_vertices_{}".format(slot)], dtype=float)
        faces = [tuple(int(v) for v in face) for face in job["mesh_faces_{}".format(slot)]]
        colour = tuple(float(v) for v in job["mesh_color_{}".format(slot)])
        material = glossy("cage", colour, 0.22, 0.12)
        mesh = bpy.data.meshes.new("cage")
        mesh.from_pydata(vertices.tolist(), [], faces)
        mesh.update()
        for poly in mesh.polygons:
            poly.use_smooth = False
        add_mesh("cage", mesh, material)

    from mathutils import Vector

    ortho = bool(int(job["ortho"])) if "ortho" in job.files else False
    camera_data = bpy.data.cameras.new("camera")
    camera_data.sensor_fit = "VERTICAL"
    camera = bpy.data.objects.new("camera", camera_data)
    bpy.context.scene.collection.objects.link(camera)
    target = Vector((float(centre[0]), float(centre[1]), float(centre[2])))
    if ortho:
        camera_data.type = "ORTHO"
        # The atom box fills about 92 percent of the limiting frame edge.
        vertical = vertical / 0.95
        camera_data.ortho_scale = vertical
        camera.location = target + Vector((0.0, 0.0, 40.0))
        camera_data.clip_start = 0.1
        camera_data.clip_end = 400.0
    else:
        camera_data.type = "PERSP"
        camera_data.lens = 65.0
        camera_data.sensor_width = 36.0
        half = 0.5 * vertical
        fov = 2.0 * float(np.arctan(18.0 / 65.0))
        distance = half / np.tan(fov / 2.0) / 0.90
        camera.location = target + Vector((0.0, 0.0, distance))
        vertical = 2.0 * distance * float(np.tan(fov / 2.0))
        camera_data.clip_start = 0.1
        camera_data.clip_end = max(800.0, distance * 4.0)
    camera.rotation_euler = (target - camera.location).to_track_quat("-Z", "Y").to_euler()
    bpy.context.scene.camera = camera

    def area(name, location, energy, size):
        from mathutils import Vector

        data = bpy.data.lights.new(name, "AREA")
        data.energy = energy
        data.size = size
        data.shape = "DISK"
        obj = bpy.data.objects.new(name, data)
        bpy.context.scene.collection.objects.link(obj)
        obj.location = location
        direction = Vector(look) - Vector(location)
        obj.rotation_euler = direction.to_track_quat("-Z", "Y").to_euler()
        return obj

    area("key", (float(centre[0]) - 9.0, float(centre[1]) + 11.0, float(high[2]) + 14.0), 1800.0, 18.0)
    area("fill", (float(centre[0]) + 14.0, float(centre[1]) + 4.0, float(high[2]) + 18.0), 480.0, 36.0)
    area("rim", (float(centre[0]) + 1.0, float(centre[1]) - 16.0, float(low[2]) + 5.0), 520.0, 10.0)

    world = bpy.data.worlds.new("studio")
    bpy.context.scene.world = world
    world.use_nodes = True
    backdrop = world.node_tree.nodes.get("Background")
    if backdrop is not None:
        backdrop.inputs["Color"].default_value = (*background, 1.0)
        backdrop.inputs["Strength"].default_value = 0.34 if light else 0.08

    panel = str(job["panel"]) if "panel" in job.files else "witness"
    scene = bpy.context.scene
    scene.render.engine = "CYCLES"
    scene.cycles.device = "CPU"
    scene.cycles.samples = int(job["samples"])
    scene.cycles.use_denoising = True
    try:
        scene.cycles.denoiser = "OPENIMAGEDENOISE"
    except TypeError:
        pass
    scene.cycles.use_adaptive_sampling = False
    try:
        scene.cycles.use_fast_gi = True
    except AttributeError:
        pass
    scene.cycles.adaptive_min_samples = 96
    try:
        scene.cycles.adaptive_threshold = 0.01
    except (AttributeError, TypeError):
        pass
    scene.cycles.max_bounces = 6
    scene.cycles.transmission_bounces = 6
    scene.cycles.glossy_bounces = 2
    for attr, value in (("sample_clamp_direct", 8.0), ("sample_clamp_indirect", 4.0)):
        if hasattr(scene.cycles, attr):
            setattr(scene.cycles, attr, value)
    if str(panel).startswith("witness"):
        scene.render.resolution_x = 640
        scene.render.resolution_y = 480
    else:
        scene.render.resolution_x = WIDTH * RENDER_SCALE
        scene.render.resolution_y = HEIGHT * RENDER_SCALE
    scene.render.resolution_percentage = 100
    scene.render.film_transparent = True
    scene.render.image_settings.file_format = "PNG"
    scene.render.image_settings.color_mode = "RGBA"
    scene.view_settings.view_transform = "Standard"
    scene.view_settings.look = "None"
    scene.view_settings.exposure = 0.0
    scene.view_settings.gamma = 1.0
    scene.cycles.use_adaptive_sampling = False
    scene.render.filepath = sys.argv[2]
    print("engine CYCLES samples {} panel {}".format(int(job["samples"]), panel), flush=True)
    bpy.ops.render.render(write_still=True)
    print("WORLD {:.8f}".format(vertical))
    return 0


def cycles_launch(positions, colors, radii, transparency, bonds, meshes, background, panel="panel", samples=None, ortho=False, bold=None, bond_radius=0.045, bond_rgb=None):
    """Render with Cycles in a fresh interpreter. None when bpy is absent."""
    import subprocess
    import tempfile

    import numpy as np

    program = cycles_python()
    if program is None:
        return None
    handle = tempfile.NamedTemporaryFile(suffix=".npz", delete=False)
    picture = tempfile.NamedTemporaryFile(suffix=".png", delete=False)
    handle.close()
    picture.close()
    payload = {
        "positions": np.asarray(positions, dtype=np.float64),
        "colors": np.asarray(colors, dtype=np.float64),
        "radii": np.asarray(radii, dtype=np.float64),
        "transparency": np.asarray(transparency, dtype=np.float64),
        "bonds": np.asarray(bonds if len(bonds) else np.zeros((0, 2)), dtype=np.int32),
        "background": np.asarray(background, dtype=np.float64),
        "mesh_count": np.int32(len(meshes)),
        "samples": np.int32(int(os.environ.get("MOLECULES_SAMPLES", str(CYCLES_SAMPLES if samples is None else samples)))),
        "panel": np.array(str(panel)),
        "ortho": np.int32(1 if ortho else 0),
        "bond_radius": np.float64(bond_radius),
        "bond_rgb": np.asarray(bond_rgb if bond_rgb is not None else (0.42, 0.44, 0.47), dtype=np.float64),
    }
    if bold is not None and len(bold):
        payload["bold_bonds"] = np.asarray(bold, dtype=np.int32)
    for slot, mesh in enumerate(meshes):
        payload["mesh_vertices_{}".format(slot)] = np.asarray(mesh["vertices"], dtype=np.float64)
        payload["mesh_faces_{}".format(slot)] = np.asarray(mesh["faces"], dtype=np.int32)
        payload["mesh_color_{}".format(slot)] = np.asarray(mesh["color"], dtype=np.float64)
    np.savez(handle.name, **payload)
    env = os.environ.copy()
    env["MOLECULES_CYCLES"] = "1"
    env.pop("OMP_NUM_THREADS", None)
    try:
        done = subprocess.run(
            [program, str(Path(__file__).resolve()), handle.name, picture.name],
            check=False,
            capture_output=True,
            text=True,
            env=env,
            timeout=1800,
        )
        world = 0.0
        for line in (done.stdout or "").splitlines():
            if line.startswith("WORLD "):
                world = float(line.split()[1])
            elif line.startswith("engine "):
                print(line, flush=True)
        if done.returncode != 0 or not os.path.isfile(picture.name) or os.path.getsize(picture.name) < 1000:
            detail = (done.stderr or done.stdout or "cycles failed")[-500:]
            print(detail, file=sys.stderr)
            return None
        from PIL import Image

        image = Image.open(picture.name).convert("RGBA")
        plate = tuple(int(round(channel * 255)) for channel in background)
        canvas = Image.new("RGBA", image.size, plate + (255,))
        image = Image.alpha_composite(canvas, image).convert("RGB")
        image = image.resize((WIDTH, HEIGHT), Image.Resampling.LANCZOS)
        image.load()
        return image, world
    except (OSError, subprocess.TimeoutExpired, ValueError) as exc:
        print(str(exc), file=sys.stderr)
        return None
    finally:
        for path in (handle.name, picture.name):
            try:
                os.unlink(path)
            except OSError:
                pass


def raytrace(positions, colors, radii, transparency, bonds, meshes, background, perspective, ortho=False, bold=None, bond_radius=0.045, bond_rgb=None):
    """Cycles when requested. Tachyon is only used when the engine is tachyon."""
    if RENDER_ENGINE == "cycles":
        drawn = cycles_launch(
            positions, colors, radii, transparency, bonds, meshes, background,
            panel=PANEL_NAME, ortho=ortho, bold=bold,
            bond_radius=bond_radius, bond_rgb=bond_rgb,
        )
        if drawn is None:
            fail("Cycles was requested (--engine cycles) and did not produce a picture")
        return drawn
    print("engine TACHYON samples 64 panel {}".format(PANEL_NAME), flush=True)
    return tachyon_raytrace(
        positions, colors, radii, transparency, bonds, meshes, background, perspective
    )


def tachyon_raytrace(positions, colors, radii, transparency, bonds, meshes, background, perspective):
    """Tachyon fallback. Antialiasing stays at 64 samples. Depth of field stays off."""
    import tempfile

    import numpy as np
    import ovito
    from ovito.data import DataCollection, Particles, TriangleMesh
    from ovito.pipeline import Pipeline, StaticSource
    from ovito.vis import TachyonRenderer, Viewport
    from PIL import Image

    warnings.filterwarnings("ignore", message=".*OVITO.*PyPI")
    while ovito.scene.pipelines:
        ovito.scene.pipelines[0].remove_from_scene()
    data = DataCollection()
    particles = Particles()
    particles.create_property("Position", data=np.asarray(positions, dtype=float))
    particles.create_property("Color", data=np.asarray(colors, dtype=float))
    particles.create_property("Radius", data=np.asarray(radii, dtype=float))
    particles.create_property("Transparency", data=np.asarray(transparency, dtype=float))
    data.objects.append(particles)
    light = background[0] > 0.5
    if len(bonds):
        particles.create_bonds().create_property("Topology", data=np.asarray(bonds, dtype=int))
        particles.bonds.vis.width = 0.040
        particles.bonds.vis.color = (0.28, 0.28, 0.30) if light else (0.78, 0.79, 0.82)
    particles.vis.radius = 0.38
    for mesh in meshes:
        tri = TriangleMesh()
        tri.set_vertices(np.asarray(mesh["vertices"], dtype=float))
        tri.set_faces(np.asarray(mesh["faces"], dtype=int))
        tri.vis.color = tuple(float(channel) for channel in mesh["color"])
        tri.vis.transparency = float(mesh["transparency"])
        tri.vis.backface_culling = False
        tri.vis.highlight_edges = True
        tri.vis.wireframe_color = (0.20, 0.20, 0.22) if light else (0.86, 0.87, 0.90)
        data.objects.append(tri)
    pipeline = Pipeline(source=StaticSource(data=data))
    pipeline.add_to_scene()
    handle = tempfile.NamedTemporaryFile(suffix=".png", delete=False)
    handle.close()
    try:
        view = Viewport(type=Viewport.Type.Perspective if perspective else Viewport.Type.Ortho)
        # zoom_all keeps the side the camera already occupies, so start on +z.
        view.camera_pos = (0.0, 0.0, 40.0)
        view.camera_dir = (0.0, 0.0, -1.0)
        view.camera_up = (0.0, 1.0, 0.0)
        if perspective:
            view.fov = 28.0
        low, high = frame_bounds(positions, radii, transparency)
        centre = (low + high) / 2.0
        span_x = float(high[0] - low[0])
        span_y = float(high[1] - low[1])
        aspect = WIDTH / float(HEIGHT)
        look_x = centre[0] + 0.04 * span_x
        look_y = centre[1] - 0.02 * span_y
        view.camera_dir = (0.0, 0.0, -1.0)
        view.camera_up = (0.0, 1.0, 0.0)
        view.camera_pos = (look_x, look_y, float(high[2]) + 40.0)
        if perspective:
            view.fov = 28.0
            distance = max(span_y, span_x / aspect) / (2.0 * float(np.tan(np.deg2rad(14.0))))
            view.camera_pos = (look_x, look_y, float(high[2]) + distance)
            world = 2.0 * distance * float(np.tan(np.deg2rad(float(view.fov) / 2.0)))
        else:
            world = max(span_y, span_x / aspect) * 1.12
            # Orthographic fov is the half-height of the view.
            view.fov = world / 2.0
        renderer = TachyonRenderer(
            ambient_occlusion=True,
            ambient_occlusion_brightness=0.55,
            ambient_occlusion_samples=16,
            antialiasing=True,
            antialiasing_samples=64,
            shadows=True,
            depth_of_field=False,
            aperture=0.0,
            direct_light=True,
            direct_light_intensity=0.95,
            max_ray_recursion=2,
        )
        view.render_image(
            filename=handle.name,
            size=(2400, 1800),
            background=tuple(background),
            renderer=renderer,
            frame=0,
        )
        picture = Image.open(handle.name).convert("RGB")
        return picture, world
    finally:
        pipeline.remove_from_scene()
        try:
            os.unlink(handle.name)
        except OSError:
            pass


def crop_margin(picture, light):
    import numpy as np
    from PIL import Image

    array = np.asarray(picture)
    if light:
        mask = np.any(array < 246, axis=2)
    else:
        mask = np.any(array > 40, axis=2)
    rows = np.where(mask.any(axis=1))[0]
    cols = np.where(mask.any(axis=0))[0]
    if len(rows) == 0:
        return picture, picture.size[1], picture.size[1]
    pad = 28
    top = max(0, int(rows[0]) - pad)
    bottom = min(array.shape[0], int(rows[-1]) + pad)
    left = max(0, int(cols[0]) - pad)
    right = min(array.shape[1], int(cols[-1]) + pad)
    cropped = Image.fromarray(array[top:bottom, left:right])
    return cropped, bottom - top, array.shape[0]


def _font(size):
    from PIL import ImageFont

    for path in (
        "/usr/share/fonts/truetype/dejavu/DejaVuSans.ttf",
        "/usr/share/fonts/truetype/liberation/LiberationSans-Regular.ttf",
    ):
        if os.path.isfile(path):
            return ImageFont.truetype(path, size)
    return ImageFont.load_default()


def _contact_shadow(picture, light):
    """A soft silhouette under the atoms. The page behind it stays white."""
    import numpy as np
    from PIL import Image, ImageFilter

    array = np.asarray(picture.convert("RGB"))
    if light:
        mask = np.any(array < 246, axis=2)
    else:
        mask = np.any(array > 40, axis=2)
    plate = Image.fromarray((mask.astype(np.uint8) * 255))
    plate = plate.filter(ImageFilter.GaussianBlur(radius=7))
    shade = (186, 188, 190) if light else (12, 12, 14)
    shadow = Image.new("RGB", picture.size, shade)
    return shadow, plate


def compose(structure, legend, world_height, light):
    """Centre the structure. The legend is two rows across the panel."""
    from PIL import Image, ImageDraw

    bg = (255, 255, 255) if light else (22, 24, 28)
    ink = (25, 25, 28) if light else (236, 237, 239)
    picture, crop_rows, render_rows = crop_margin(structure, light)
    font = _font(40)
    row_h = 52
    caption = 16 + row_h * 2
    avail_w = int(round(WIDTH * 0.94))
    avail_h = int(round(HEIGHT * 0.96)) - caption
    scale = min(avail_w / float(picture.size[0]), avail_h / float(picture.size[1]))
    sized = (
        max(1, int(round(picture.size[0] * scale))),
        max(1, int(round(picture.size[1] * scale))),
    )
    picture = picture.resize(sized, Image.Resampling.LANCZOS)
    canvas = Image.new("RGB", (WIDTH, HEIGHT), bg)
    left = (WIDTH - sized[0]) // 2
    top = max(12, (avail_h - sized[1]) // 2)
    shadow, plate = _contact_shadow(picture, light)
    shadow_layer = Image.new("RGBA", (WIDTH, HEIGHT), (0, 0, 0, 0))
    tint = Image.new("RGBA", picture.size, (186, 188, 190, 70) if light else (8, 8, 10, 90))
    tint.putalpha(plate.point(lambda value: int(value * 0.45)))
    shadow_layer.paste(tint, (left, top + 10), tint)
    canvas = Image.alpha_composite(canvas.convert("RGBA"), shadow_layer).convert("RGB")
    import numpy as np

    opaque = np.any(np.asarray(picture.convert("RGB")) < 246, axis=2)
    sprite = picture.convert("RGBA")
    sprite.putalpha(Image.fromarray((opaque.astype(np.uint8) * 255)))
    canvas.paste(sprite, (left, top), sprite)
    draw = ImageDraw.Draw(canvas)
    bar = 0.0
    if world_height > 0 and crop_rows > 0 and render_rows > 0:
        world_crop = world_height * (crop_rows / float(render_rows))
        bar = 10.0 * sized[1] / world_crop
    items = []
    for name, colour_rgb in legend:
        colour = tuple(int(max(0, min(255, round(channel * 255)))) for channel in colour_rgb)
        items.append((name, colour, draw.textlength(name, font=font)))
    # Two rows, split as evenly as the text allows, each entry on one line.
    if len(items) <= 1:
        cuts = [items]
    else:
        best = 1
        best_gap = 1.0e9
        widths = [44 + text_w for _name, _colour, text_w in items]
        for cut in range(1, len(items)):
            gap = abs(sum(widths[:cut]) - sum(widths[cut:]))
            if gap < best_gap:
                best, best_gap = cut, gap
        cuts = [items[:best], items[best:]]
    y0 = top + sized[1] + 8
    margin = 28
    if bar > 8:
        bar_y = y0 + 18
        draw.line((margin, bar_y, margin + bar, bar_y), fill=ink, width=8)
        draw.line((margin, bar_y - 10, margin, bar_y + 10), fill=ink, width=4)
        draw.line((margin + bar, bar_y - 10, margin + bar, bar_y + 10), fill=ink, width=4)
        draw.text((margin + bar + 10, bar_y - 22), "10 \u00c5", fill=ink, font=font)
        legend_left = margin + bar + 10 + draw.textlength("10 \u00c5", font=font) + 24
    else:
        legend_left = margin
    span = WIDTH - margin - legend_left
    for row_index, row in enumerate(cuts):
        if not row:
            continue
        y = y0 + row_index * row_h
        need = [44 + text_w for _name, _colour, text_w in row]
        spare = max(0.0, span - sum(need))
        gap = spare / (len(row) - 1) if len(row) > 1 else 0.0
        x = float(legend_left)
        if len(row) == 1:
            x = legend_left + spare / 2.0
        for (name, colour, _text_w), width in zip(row, need):
            draw.rounded_rectangle((x, y + 4, x + 28, y + 32), radius=4, fill=colour, outline=ink, width=2)
            draw.text((x + 36, y), name, fill=ink, font=font)
            x += width + gap
    return _seat(canvas, light)


def save_pair(image, directory, stem):
    from PIL import Image

    directory.mkdir(parents=True, exist_ok=True)
    light, dark = image if isinstance(image, tuple) else (image, None)
    light.save(directory / "{}.png".format(stem))
    light.resize(WEB, Image.Resampling.LANCZOS).save(directory / "{}-web.png".format(stem))
    if dark is not None:
        dark.save(directory / "{}-dark.png".format(stem))
        dark.resize(WEB, Image.Resampling.LANCZOS).save(directory / "{}-web-dark.png".format(stem))
    return directory / "{}.png".format(stem)


def largest_cluster(frame, labels, shell=3.6):
    """The biggest crystal component, plus waters that touch it."""
    import numpy as np
    from scipy.spatial import cKDTree

    crystal = np.array([i for i, label in enumerate(labels) if label in (0, 1, 3, 7, 8)], dtype=int)
    if len(crystal) == 0:
        return central_slice(frame, labels, 400)
    tree = cKDTree(frame.xyz[crystal])
    parent = np.arange(len(crystal))

    def find(node):
        while parent[node] != node:
            parent[node] = parent[parent[node]]
            node = parent[node]
        return node

    for left, right in tree.query_pairs(3.45):
        root_left, root_right = find(left), find(right)
        if root_left != root_right:
            parent[root_right] = root_left
    roots = [find(i) for i in range(len(crystal))]
    best = max(set(roots), key=roots.count)
    members = crystal[[i for i, root in enumerate(roots) if root == best]]
    near, _ = cKDTree(frame.xyz[members]).query(frame.xyz, distance_upper_bound=shell)
    shell_ids = np.where(np.isfinite(near) & (near <= shell))[0]
    keep = sorted(set(int(i) for i in members) | set(int(i) for i in shell_ids))
    return frame.subset(keep), [labels[i] for i in keep], frame.xyz[members].mean(axis=0)


def genice_atoms(lattice, rep, guests):
    """Real water and guest molecules from GenIce. Virtual sites are dropped."""
    import subprocess

    command = [
        "genice2", lattice, "-r", str(rep[0]), str(rep[1]), str(rep[2]),
        "-w", "tip4p", "-f", "gromacs", "--quiet",
    ]
    for guest in guests:
        command.extend(["-g", guest])
    done = subprocess.run(command, check=False, capture_output=True, text=True)
    if done.returncode != 0:
        fail(done.stderr[-600:])
    lines = [line for line in done.stdout.splitlines() if line.strip()]
    count = int(lines[1])
    atoms = []
    for line in lines[2 : 2 + count]:
        residue = line[5:10].strip()
        atom = line[10:15].strip()
        if atom == "MW":
            continue
        xyz = [float(value) * 10.0 for value in line[20:].split()[:3]]
        atoms.append((residue, atom, xyz, int(line[0:5])))
    box = [float(value) * 10.0 for value in lines[2 + count].split()[:3]]
    return atoms, box


def cage_meshes(oxygens, guests, box):
    """One small cage and one large cage, each wrapped around its guest.

    Waters across the periodic boundary are unfolded next to the guest
    so the polyhedron and the molecule stay in one piece.
    """
    import numpy as np
    from scipy.spatial import ConvexHull

    if len(guests) == 0 or len(oxygens) < 8:
        return [], np.zeros((0, 3)), []
    lengths = np.asarray(box[:3], dtype=float)
    origin = np.asarray(oxygens, dtype=float).mean(axis=0)
    found = []
    for kind, heavy in guests:
        centre = np.asarray(heavy, dtype=float).mean(axis=0)
        delta = np.asarray(oxygens, dtype=float) - centre
        delta -= lengths * np.round(delta / lengths)
        distance = np.linalg.norm(delta, axis=1)
        limit = 4.45 if kind == "CH4" else 5.45
        chosen = (distance > 2.4) & (distance < limit)
        if int(np.count_nonzero(chosen)) < (16 if kind == "CH4" else 24):
            continue
        shell = centre + delta[chosen]
        try:
            hull = ConvexHull(shell, qhull_options="QJ")
        except Exception:
            continue
        token = "hexagonal" if kind == "CH4" else "clathrate"
        found.append((float(np.linalg.norm(centre - origin)), kind, shell, np.where(chosen)[0], {
            "vertices": shell,
            "faces": hull.simplices.astype(int),
            "color": palette_rgb(token),
            "transparency": 0.82,
            "legend": "small cage" if kind == "CH4" else "large cage",
        }))
    meshes = []
    extra = []
    highlighted = []
    oxygens = np.asarray(oxygens, dtype=float)
    for kind in ("CH4", "THF"):
        ranked = sorted(item for item in found if item[1] == kind)
        if not ranked:
            continue
        _distance, _kind, shell, indices, mesh = ranked[0]
        meshes.append(mesh)
        highlighted.append(np.asarray(indices, dtype=int))
        if len(oxygens):
            nearest = np.linalg.norm(shell[:, None, :] - oxygens[None, :, :], axis=2).min(axis=1)
            extra.append(shell[nearest > 1.2])
    if extra:
        extra = np.vstack(extra)
    else:
        extra = np.zeros((0, 3))
    return meshes, extra, highlighted


def ring_meshes(positions, rings):
    """Filled polygons. Six-rings are orange and every other ring is blue."""
    import numpy as np

    buckets = {
        "hexagonal ring": {"token": "hexagonal", "vertices": [], "faces": []},
        "other ring": {"token": "cubic", "vertices": [], "faces": []},
    }
    for size, indices in rings:
        if len(indices) < 3:
            continue
        key = "hexagonal ring" if size == 6 else "other ring"
        bucket = buckets[key]
        points = np.asarray(positions)[list(indices)]
        centre = points.mean(axis=0)
        normal = np.cross(points[1] - points[0], points[min(2, len(points) - 1)] - points[0])
        length = float(np.linalg.norm(normal)) or 1.0
        centre = centre + normal / length * 0.12
        points = points + normal / length * 0.12
        base = len(bucket["vertices"])
        bucket["vertices"].append(centre)
        bucket["vertices"].extend(points.tolist())
        count = len(points)
        for step in range(count):
            bucket["faces"].append((base, base + 1 + step, base + 1 + ((step + 1) % count)))
    meshes = []
    for name, bucket in buckets.items():
        if not bucket["faces"]:
            continue
        meshes.append({
            "vertices": np.asarray(bucket["vertices"], dtype=float),
            "faces": np.asarray(bucket["faces"], dtype=int),
            "color": palette_rgb(bucket["token"]),
            "transparency": 0.55,
            "legend": name,
        })
    return meshes


def finish(positions, colors, radii, transparency, names, bonds, meshes, view, center, perspective, bold=None, align=False, bond_radius=0.045, bond_rgb=None):
    import numpy as np

    positions = np.asarray(positions, dtype=float)
    colors = np.asarray(colors, dtype=float)
    transparency = np.asarray(transparency, dtype=float)
    radii = np.asarray(radii, dtype=float)
    if align and len(positions) > 8:
        solid = transparency < 0.5
        if int(np.count_nonzero(solid)) >= 8:
            origin, basis = _principal_basis(positions, solid)
            positions = _spin(_apply_basis(positions, origin, basis), "x", 12.0)
            posed = []
            for mesh in meshes:
                posed.append({
                    "vertices": _spin(_apply_basis(mesh["vertices"], origin, basis), "x", 12.0),
                    "faces": mesh["faces"],
                    "color": mesh["color"],
                    "transparency": mesh["transparency"],
                    "legend": mesh.get("legend"),
                })
            meshes = posed
        else:
            positions = orient(positions, view, center)
            meshes = [
                {
                    "vertices": orient(mesh["vertices"], view, center),
                    "faces": mesh["faces"],
                    "color": mesh["color"],
                    "transparency": mesh["transparency"],
                    "legend": mesh.get("legend"),
                }
                for mesh in meshes
            ]
    else:
        if view == "basal":
            positions = _basal_flat(positions, center)
            meshes = [
                {
                    "vertices": _basal_flat(mesh["vertices"], center),
                    "faces": mesh["faces"],
                    "color": mesh["color"],
                    "transparency": mesh["transparency"],
                    "legend": mesh.get("legend"),
                }
                for mesh in meshes
            ]
            if len(positions) > 1:
                colors, transparency, radii = _basal_fog(positions, colors, transparency, radii)
            positions = _spin(positions, "x", 4.0)
            for mesh in meshes:
                mesh["vertices"] = _spin(mesh["vertices"], "x", 4.0)
        else:
            positions = orient(positions, view, center)
            meshes = [
                {
                    "vertices": orient(mesh["vertices"], view, center),
                    "faces": mesh["faces"],
                    "color": mesh["color"],
                    "transparency": mesh["transparency"],
                    "legend": mesh.get("legend"),
                }
                for mesh in meshes
            ]
    legend = legend_rows(names, colors)
    seen = {name for name, _colour in legend}
    for mesh in meshes:
        label = mesh.get("legend")
        if label and label not in seen:
            legend.append((label, mesh["color"]))
            seen.add(label)
    positions = np.asarray(positions, dtype=float)
    radii = np.asarray(radii, dtype=float)
    for index, name in enumerate(names):
        if name in ("methane", "THF"):
            positions[index, 2] += 4.6 if name == "methane" else 2.2
            radii[index] = max(float(radii[index]), 0.62 if name == "methane" else 0.50)
        elif name == "_rim":
            positions[index, 2] -= 0.25
    ortho = view in ("basal", "channel") or bool(align)
    light_picture, world = raytrace(
        positions, colors, radii, transparency, bonds, meshes, (1, 1, 1), perspective,
        ortho=ortho, bold=bold, bond_radius=bond_radius, bond_rgb=bond_rgb,
    )
    if os.environ.get("MOLECULES_THEMES") == "light":
        dark_picture = light_picture
    else:
        dark_picture, _world = raytrace(
            positions, colors, radii, transparency, bonds, meshes, (0.09, 0.10, 0.12), perspective,
            ortho=ortho, bold=bold, bond_radius=bond_radius, bond_rgb=bond_rgb,
        )
    return (
        compose(light_picture, legend, world, True),
        compose(dark_picture, legend, world, False),
    )


def labelled_picture(frame, labels, view, liquid_clear, ions=False, perspective=False, meshes=None, center=None, cluster=False):
    import numpy as np

    colors, names, radii, transparency, heavy = [], [], [], [], []
    kept_xyz = []
    kept_labels = []
    for index in range(len(frame.xyz)):
        kind = int(frame.types[index])
        if ions and kind in (3, 4):
            token = "sodium" if kind == 3 else "chloride"
            colour, radius, alpha = atom_style(token)
            rim = (0.12, 0.12, 0.14)
            colors.append(rim)
            names.append("_rim")
            radii.append(radius * 1.18)
            transparency.append(0.0)
            heavy.append(False)
            kept_xyz.append(frame.xyz[index])
            kept_labels.append(int(labels[index]))
            colors.append(colour)
            names.append(token)
            radii.append(radius)
            transparency.append(alpha)
            heavy.append(False)
            kept_xyz.append(frame.xyz[index])
            kept_labels.append(int(labels[index]))
            continue
        label = int(labels[index])
        if cluster and label == 2:
            continue
        if cluster and label not in (0, 1, 3, 7, 8):
            continue
        token = "mixed" if cluster and label == 3 else CHILL_NAME.get(label, "unclassified")
        colour, radius, alpha = atom_style(token)
        colors.append(colour)
        names.append(token)
        radii.append(radius)
        transparency.append(alpha)
        heavy.append(True)
        kept_xyz.append(frame.xyz[index])
        kept_labels.append(label)
    xyz = np.asarray(kept_xyz, dtype=float)
    bonds = oxygen_bonds(xyz, heavy)
    if cluster:
        crystal = {0, 1, 3, 7, 8}
        bonds = [
            pair for pair in bonds
            if int(kept_labels[pair[0]]) in crystal and int(kept_labels[pair[1]]) in crystal
        ]
    if cluster and bonds:
        keep = set(range(len(xyz)))
        for _pass in range(4):
            degree = {index: 0 for index in keep}
            for left, right in bonds:
                if left in keep and right in keep:
                    degree[left] += 1
                    degree[right] += 1
            drop = {index for index, count in degree.items() if count < 2}
            if not drop or len(keep - drop) < 12:
                break
            keep -= drop
        order = sorted(keep)
        remap = {old: new for new, old in enumerate(order)}
        xyz = xyz[order]
        colors = [colors[index] for index in order]
        names = [names[index] for index in order]
        radii = [radii[index] for index in order]
        transparency = [transparency[index] for index in order]
        kept_labels = [kept_labels[index] for index in order]
        bonds = [
            (remap[left], remap[right])
            for left, right in bonds
            if left in remap and right in remap
        ]
    elif liquid_clear:
        bonds = [
            pair for pair in bonds
            if not (int(kept_labels[pair[0]]) == 2 and int(kept_labels[pair[1]]) == 2)
        ]
    if view == "channel":
        radii = [value * (0.72 if name not in ("sodium", "chloride", "_rim") else 1.0) for value, name in zip(radii, names)]
    bond_radius = 0.038 if view == "channel" else 0.050
    return finish(
        xyz, colors, radii, transparency, names, bonds, meshes or [],
        view, frame.xyz.mean(axis=0) if center is None else center, perspective,
        align=cluster, bond_radius=bond_radius,
    )


def ring_lists(library, frame):
    from ctypes import POINTER, byref, c_double, c_int, c_size_t, c_void_p

    count = len(frame.xyz)
    xyz = (c_double * (count * 3))(*frame.xyz.reshape(-1).tolist())
    box = (c_double * 9)(*frame.box9)
    types = (c_int * count)(*frame.types.tolist())
    system = library.seams_system_create()
    if library.seams_system_set_frame(system, xyz, count, box, types, None) != 0:
        library.seams_system_destroy(system)
        fail("seams_system_set_frame failed")
    rings = c_void_p()
    if library.seams_rings_find(system, 3.5, 0, 0, 6, byref(rings)) != 0:
        library.seams_system_destroy(system)
        fail("seams_rings_find failed")
    found = []
    buffer = (c_int * 16)()
    written = c_size_t()
    for index in range(library.seams_rings_count(rings)):
        if library.seams_rings_vertices(rings, index, buffer, 16, written) != 0:
            continue
        size = int(written.value)
        found.append((size, [int(buffer[slot]) for slot in range(size)]))
    library.seams_rings_destroy(rings)
    library.seams_system_destroy(system)
    return found


def reimage_atoms(atoms, box):
    """Keep every molecule whole and put each guest inside one cage.

    Each guest moves by one lattice vector, chosen so the periodic
    image has the most cage waters around it and does not land on
    another guest. Every atom of the residue shares that vector.
    """
    import numpy as np
    from collections import defaultdict

    lengths = np.asarray(box[:3], dtype=float)
    grouped = defaultdict(list)
    for residue, atom, xyz, resid in atoms:
        grouped[resid].append((residue, atom, np.asarray(xyz, dtype=float), resid))
    molecules = []
    kinds = []
    for items in grouped.values():
        heavy_at = [xyz for _residue, atom, xyz, _resid in items if not str(atom).startswith("H")]
        base = heavy_at[0] if heavy_at else items[0][2]
        local = []
        for residue, atom, xyz, resid in items:
            delta = xyz - base
            delta -= lengths * np.round(delta / lengths)
            local.append((residue, atom, base + delta, resid))
        molecules.append(local)
        kinds.append(local[0][0])
    anchors = []
    for mol in molecules:
        heavy = [xyz for _residue, atom, xyz, _resid in mol if not str(atom).startswith("H")]
        anchors.append(np.mean(heavy if heavy else [mol[0][2]], axis=0))
    anchors = np.asarray(anchors, dtype=float)
    oxygen = np.asarray([
        xyz for mol in molecules for _residue, atom, xyz, _resid in mol if atom == "OW"
    ], dtype=float)
    guests = [index for index, kind in enumerate(kinds) if kind in ("CH4", "THF")]

    def images_of(index):
        limit = 4.55 if kinds[index] == "CH4" else 5.45
        found = []
        for image in np.ndindex(3, 3, 3):
            shift = (np.asarray(image, dtype=float) - 1.0) * lengths
            trial = anchors[index] + shift
            count = int(np.sum(np.linalg.norm(oxygen - trial, axis=1) < limit))
            found.append((count, shift))
        found.sort(key=lambda item: item[0], reverse=True)
        return found

    occupied = []
    for index in sorted(guests, key=lambda item: images_of(item)[0][0], reverse=True):
        chosen = np.zeros(3)
        for count, shift in images_of(index):
            trial = anchors[index] + shift
            if occupied and float(np.linalg.norm(trial - np.asarray(occupied), axis=1).min()) < 4.5:
                continue
            chosen = shift
            break
        occupied.append(anchors[index] + chosen)
        molecules[index] = [
            (residue, atom, xyz + chosen, resid)
            for residue, atom, xyz, resid in molecules[index]
        ]
    seated = []
    for mol in molecules:
        seated.extend(mol)
    return seated


def host_token(lattice):
    if lattice in ("1h", "ice1h"):
        return "hexagonal"
    if lattice in ("1c", "ice1c"):
        return "cubic"
    return "clathrate"


def genice_picture(lattice, rep, guests, cage_faces, hydrogens=True):
    import numpy as np

    atoms, box = genice_atoms(lattice, rep, guests)
    atoms = reimage_atoms(atoms, box)
    host = host_token(lattice)
    oxygen_now = np.asarray(
        [xyz for residue, atom, xyz, _resid in atoms if atom == "OW"], dtype=float,
    )
    weak_guests = set()
    guest_heavy = {}
    for residue, atom, xyz, resid in atoms:
        if residue in ("CH4", "THF") and not str(atom).startswith("H"):
            guest_heavy.setdefault(resid, {"kind": residue, "xyz": []})["xyz"].append(xyz)
    for resid, group in guest_heavy.items():
        centre = np.mean(group["xyz"], axis=0)
        limit = 4.55 if group["kind"] == "CH4" else 5.45
        if len(oxygen_now) and int(np.sum(np.linalg.norm(oxygen_now - centre, axis=1) < limit)) < 8:
            weak_guests.add(resid)
    positions, colors, names, radii, transparency = [], [], [], [], []
    residue_index = []
    oxygens = []
    guest_groups = {}
    for residue, atom, xyz, resid in atoms:
        if resid in weak_guests:
            continue
        if str(atom).startswith("H") and not hydrogens:
            continue
        positions.append(xyz)
        residue_index.append(resid)
        if residue in ("CH4", "THF"):
            kind = "CH4" if residue == "CH4" else "THF"
            group = guest_groups.setdefault(resid, {"kind": kind, "heavy": []})
            token = "methane" if kind == "CH4" else "THF"
            if str(atom).startswith("H"):
                colour, radius, alpha = atom_style("hydrogen")
                names.append("hydrogen")
            else:
                group["heavy"].append(xyz)
                colour, radius, alpha = atom_style(token)
                names.append(token)
            colors.append(colour)
            radii.append(radius)
            transparency.append(alpha)
        elif str(atom).startswith("H"):
            colour, radius, alpha = atom_style("hydrogen")
            colors.append(colour)
            names.append("hydrogen")
            radii.append(radius)
            transparency.append(alpha)
        else:
            oxygens.append(xyz)
            colour, radius, alpha = atom_style(host)
            colors.append(colour)
            names.append(host)
            radii.append(radius)
            transparency.append(alpha)
    xyz = np.asarray(positions, dtype=float)
    bonds = []
    for left in range(len(xyz)):
        for right in range(left + 1, len(xyz)):
            if residue_index[left] != residue_index[right]:
                continue
            gap = float(np.linalg.norm(xyz[left] - xyz[right]))
            if gap < 1.70:
                bonds.append((left, right))
    oxygen_ids = [i for i, name in enumerate(names) if name == host]
    for cursor, left in enumerate(oxygen_ids):
        for right in oxygen_ids[cursor + 1 :]:
            gap = float(np.linalg.norm(xyz[left] - xyz[right]))
            if 2.40 <= gap <= 3.20:
                bonds.append((left, right))
    meshes = []
    bold = []
    if cage_faces:
        meshes, extra, highlighted = cage_meshes(
            np.asarray(oxygens, dtype=float),
            [
                (group["kind"], np.asarray(group["heavy"], dtype=float))
                for group in guest_groups.values()
                if group["heavy"]
            ],
            box,
        )
        cage_atoms = set()
        for indices in highlighted:
            for slot in indices:
                cage_atoms.add(oxygen_ids[int(slot)])
        # Draw the methane and the THF that sit in the two glass cages.
        # The other guests sit on the cell boundary and read as clutter.
        featured = set()
        for mesh in meshes:
            centre = np.asarray(mesh["vertices"], dtype=float).mean(axis=0)
            best_resid, best_gap = None, 1.0e9
            for resid, group in guest_groups.items():
                gap = float(np.linalg.norm(np.mean(group["heavy"], axis=0) - centre))
                if gap < best_gap:
                    best_resid, best_gap = resid, gap
            if best_resid is not None and best_gap < 3.0:
                featured.add(best_resid)
        guest_resids = set(guest_groups)
        if not featured:
            keep = list(range(len(residue_index)))
        else:
            keep = [
                index for index, resid in enumerate(residue_index)
                if resid not in guest_resids or resid in featured
            ]
        remap = {old: new for new, old in enumerate(keep)}
        xyz = xyz[keep]
        colors = [colors[index] for index in keep]
        names = [names[index] for index in keep]
        radii = [radii[index] for index in keep]
        transparency = [transparency[index] for index in keep]
        bonds = [
            (remap[pair[0]], remap[pair[1]])
            for pair in bonds
            if pair[0] in remap and pair[1] in remap
        ]
        cage_atoms = {remap[index] for index in cage_atoms if index in remap}
        oxygen_ids = [remap[index] for index in oxygen_ids if index in remap]
        water_colour = palette_rgb("water")
        centres = [np.asarray(mesh["vertices"], dtype=float).mean(axis=0) for mesh in meshes]
        near_cage = set()
        if centres and oxygen_ids:
            stacked = np.vstack(centres)
            distance = np.linalg.norm(xyz[oxygen_ids, None, :] - stacked[None, :, :], axis=2).min(axis=1)
            for atom_index, gap in zip(oxygen_ids, distance):
                if float(gap) <= 6.4:
                    near_cage.add(int(atom_index))
        cage_atoms = {index for index in cage_atoms if index in near_cage}
        guest_names = {"methane", "THF", "hydrogen"}
        keep = [
            index for index, name in enumerate(names)
            if name in guest_names or index in near_cage
        ]
        remap = {old: new for new, old in enumerate(keep)}
        xyz = xyz[keep]
        colors = [colors[index] for index in keep]
        names = [names[index] for index in keep]
        radii = [radii[index] for index in keep]
        transparency = [transparency[index] for index in keep]
        cage_atoms = {remap[index] for index in cage_atoms if index in remap}
        bonds = [
            (remap[pair[0]], remap[pair[1]])
            for pair in bonds
            if pair[0] in remap and pair[1] in remap and remap[pair[0]] in cage_atoms and remap[pair[1]] in cage_atoms
        ]
        for index, name in enumerate(names):
            if name in guest_names:
                continue
            if index in cage_atoms:
                radii[index] = 0.32
                transparency[index] = 0.0
            else:
                radii[index] = 0.24
                transparency[index] = 0.35
            colors[index] = water_colour
            names[index] = "water"
        if len(extra):
            xyz = np.vstack([xyz, np.asarray(extra, dtype=float)])
            for _point in extra:
                colors.append(water_colour)
                names.append("water")
                radii.append(0.32)
                transparency.append(0.0)
        bold = list(bonds)
    return finish(
        xyz, colors, radii, transparency, names, bonds, meshes,
        "cage", xyz.mean(axis=0), False, bold=bold,
        bond_radius=0.036, bond_rgb=(0.33, 0.35, 0.38),
    )


def whole_hexagon_index(points):
    """Atoms that sit on a closed six-ring. Edge fragments are left out."""
    import numpy as np
    from scipy.spatial import cKDTree

    pts = np.asarray(points, dtype=float)
    if len(pts) < 6:
        return np.arange(len(pts))
    tree = cKDTree(pts)
    adj = [[] for _ in range(len(pts))]
    for left, right in tree.query_pairs(3.25):
        gap = float(np.linalg.norm(pts[left] - pts[right]))
        if gap < 2.30:
            continue
        adj[left].append(int(right))
        adj[right].append(int(left))
    keep = np.zeros(len(pts), dtype=bool)
    for start, firsts in enumerate(adj):
        for second in firsts:
            if second <= start:
                continue
            for third in adj[second]:
                if third == start:
                    continue
                for fourth in adj[third]:
                    if fourth in (start, second):
                        continue
                    for fifth in adj[fourth]:
                        if fifth in (start, second, third):
                            continue
                        for sixth in adj[fifth]:
                            if sixth in (second, third, fourth, fifth):
                                continue
                            if start in adj[sixth] and len({start, second, third, fourth, fifth, sixth}) == 6:
                                keep[[start, second, third, fourth, fifth, sixth]] = True
    if int(keep.sum()) < 12:
        return np.arange(len(pts))
    return np.where(keep)[0]


def tile_frame(frame, labels, nx, ny, nz):
    """Repeat a crystal by its box so the basal view holds more channels."""
    import numpy as np

    lengths = np.asarray(frame.box9[3:6], dtype=float)
    pieces, kinds, repeated = [], [], []
    for ix in range(nx):
        for iy in range(ny):
            for iz in range(nz):
                shift = lengths * np.array([ix, iy, iz], dtype=float)
                pieces.append(np.asarray(frame.xyz, dtype=float) + shift)
                kinds.append(np.asarray(frame.types))
                repeated.extend(labels)
    box = list(frame.box9)
    box[3] *= nx
    box[4] *= ny
    box[5] *= nz
    return Frame(np.vstack(pieces), np.concatenate(kinds), box), repeated


def clip_axis(frame, labels, axis, low, high):
    """Keep the slice of a long cell that still shows the lattice."""
    import numpy as np

    mask = (frame.xyz[:, axis] >= low) & (frame.xyz[:, axis] <= high)
    index = np.where(mask)[0]
    return frame.subset(index), [labels[int(i)] for i in index]


def draw(library, scene):
    global PANEL_NAME
    PANEL_NAME = scene
    if scene == "ice-ih":
        frame = read_dump(ROOT / "input/traj/ice1h.lammpstrj")
        labels = chill_labels(library, frame)
        frame, labels = tile_frame(frame, labels, 3, 1, 2)
        kept = whole_hexagon_index(frame.xyz)
        frame = frame.subset(kept)
        labels = [labels[int(index)] for index in kept]
        return labelled_picture(frame, labels, "basal", False)
    if scene == "ice-ih-h":
        return genice_picture("1h", (2, 2, 1), [], False)
    if scene == "ice-ic":
        frame = read_dump(ROOT / "input/traj/mW_cubic.lammpstrj")
        labels = chill_labels(library, frame)
        frame, labels = central_slice(frame, labels, 700)
        return labelled_picture(frame, labels, "basal", False)
    if scene == "interface":
        frame = read_dump(ROOT / "input/traj/interface_f0.lammpstrj")
        labels = chill_labels(library, frame)
        frame, labels = clip_axis(frame, labels, 2, 0.0, 34.0)
        return labelled_picture(frame, labels, "channel", True)
    if scene == "nucleation":
        frame = read_dump(ROOT / "input/traj/nucleation_three.lammpstrj", frame_index=1)
        labels = chill_labels(library, frame)
        frame, labels, centre = largest_cluster(frame, labels, shell=0.0)
        return labelled_picture(
            frame, labels, "cage", True, perspective=False, center=centre, cluster=True,
        )
    if scene == "brine":
        frame = read_dump(ROOT / "input/traj/nacl_front.lammpstrj")
        oxygen = [i for i, kind in enumerate(frame.types) if int(kind) == 1]
        labelled = chill_labels(library, frame.subset(oxygen))
        labels = [2] * len(frame.xyz)
        for index, label in zip(oxygen, labelled):
            labels[index] = label
        frame, labels = clip_axis(frame, labels, 2, 8.0, 50.0)
        return labelled_picture(frame, labels, "channel", True, ions=True)
    if scene == "rings":
        frame = read_dump(ROOT / "input/traj/genice_sI.lammpstrj")
        meshes = ring_meshes(frame.xyz, ring_lists(library, frame))
        return labelled_picture(
            frame, chill_labels(library, frame), "cage", False, perspective=False, meshes=meshes,
        )
    if scene == "prisms":
        frame = read_dump(ROOT / "input/traj/ice1h.lammpstrj")
        rings = [(size, ids) for size, ids in ring_lists(library, frame) if size == 6]
        upright = []
        for size, ids in rings:
            span = float(frame.xyz[ids, 2].max() - frame.xyz[ids, 2].min())
            if span > 1.4:
                upright.append((size, ids))
        meshes = ring_meshes(frame.xyz, upright)
        for mesh in meshes:
            mesh["legend"] = "prism face"
            mesh["color"] = palette_rgb("cubic")
        return labelled_picture(
            frame, chill_labels(library, frame), "basal", False, meshes=meshes,
        )
    if scene == "cage-si":
        return genice_picture("sI", (1, 1, 1), ["12=ch4", "14=ch4"], True, hydrogens=False)
    if scene == "cage-sii":
        return genice_picture("CS2", (1, 1, 1), ["12=ch4", "16=thf"], True, hydrogens=False)
    if scene == "cage-sh":
        return genice_picture("sH", (1, 1, 1), ["12=ch4", "12_1=ch4", "20=thf"], True, hydrogens=False)
    if scene == "cage-guest":
        return genice_picture("CS2", (1, 1, 1), ["12=ch4", "16=thf"], True, hydrogens=False)
    if scene in ("genice", "solvis"):
        return genice_picture("sI", (1, 1, 1), ["12=ch4", "14=ch4"], True, hydrogens=False)
    fail("unknown scene {}".format(scene))


def guest_cage():
    """One THF guest and the water within 6 angstroms of it."""
    import numpy as np

    frame = read_dump(ROOT / "input/traj/clathrate-thf.lammpstrj")
    guest = np.where(np.isin(frame.types, (3, 4, 5)))[0]
    if len(guest) == 0:
        fail("THF dump has no guest types")
    centre = frame.xyz[guest[:8]].mean(axis=0)
    distance = np.linalg.norm(frame.xyz - centre, axis=1)
    keep = np.where(distance <= 6.0)[0]
    element = {1: "O", 2: "H", 3: "C", 4: "C", 5: "O"}
    symbols = [element.get(int(frame.types[i]), "C") for i in keep]
    return frame.subset(keep).with_symbols(symbols)


def with_symbols(self, symbols):
    self.symbols = list(symbols)
    return self


Frame.with_symbols = with_symbols


def check_colour(token="hexagonal"):
    """Render one sphere and require the centre pixel to match the palette."""
    colour = palette_rgb(token)
    drawn = cycles_launch(
        [[0.0, 0.0, 0.0]],
        [colour],
        [0.60],
        [0.0],
        [],
        [],
        (1.0, 1.0, 1.0),
        panel="witness-{}".format(token),
        samples=32,
    )
    if drawn is None:
        fail("Cycles was requested (--engine cycles) and did not produce a picture")
    image, _world = drawn
    pixel = image.getpixel((image.size[0] // 2, image.size[1] // 2))
    got = tuple(channel / 255.0 for channel in pixel[:3])
    distance = delta_e(got, colour)
    print(
        "colour {} pixel {} deltaE {:.2f}".format(token, tuple(pixel[:3]), distance),
        flush=True,
    )
    if distance >= 10.0:
        fail("deltaE {:.2f} for {} is not under 10".format(distance, token))
    return distance


def fresh(directory, stem):
    path = directory / "{}.png".format(stem)
    if not path.is_file():
        return True
    from PIL import Image

    width, height = Image.open(path).size
    return width < WIDTH or height < HEIGHT


def write_pages(out):
    """Sphinx figure list. Paths are relative to docs/source/tutorials."""
    lines = [
        "Molecule stills",
        "---------------",
        "",
        "Drawn by ``python docs/scripts/render_molecules.py`` with chemparseplot.",
        "The repository stores the script. The PNG files are produced at build time.",
        "",
    ]
    order = (
        "ice-ih", "ice-ih-h", "ice-ic", "rings", "prisms", "cage-si", "cage-sii", "cage-sh",
        "cage-guest", "brine", "interface", "nucleation", "genice", "solvis",
    )
    captions = {
        "ice-ih": "Ice Ih coloured by published CHILL+.",
        "ice-ih-h": "Ice Ih with hydrogen atoms and the hydrogen-bond network.",
        "ice-ic": "Ice Ic coloured by published CHILL+.",
        "rings": "Five-rings and six-rings drawn as translucent polygons.",
        "prisms": "Upright six-ring faces of ice Ih.",
        "cage-si": "GenIce sI with a methane molecule in each cage.",
        "cage-sii": "GenIce sII with methane in the small cages and THF in the large cages.",
        "cage-sh": "GenIce sH with methane and THF guests.",
        "cage-guest": "THF guest and the surrounding water.",
        "brine": "NaCl front. Water, sodium and chloride.",
        "interface": "Basal interface coloured by published CHILL+.",
        "nucleation": "Nucleation frame coloured by published CHILL+.",
        "genice": "GenIce FAU host.",
        "solvis": "sI oxygen cell drawn with the solvis backend.",
    }
    for name in order:
        web = out / name / "overlay-web.png"
        if not web.is_file():
            continue
        lines.append(".. figure:: molecules/{}/overlay-web.png".format(name))
        lines.append("   :alt: {}".format(captions[name]))
        lines.append("   :width: 480px")
        lines.append("")
        lines.append("   {}".format(captions[name]))
        lines.append("")
    dest = out.parent / "molecules.rst"
    dest.write_text("\n".join(lines) + "\n")


def main(argv):
    if not stack_ok():
        return relaunch(argv)
    parser = argparse.ArgumentParser(description="Render molecule stills with chemparseplot")
    parser.add_argument("--out", help="Directory that receives <slug>/ stills")
    parser.add_argument("--slug", help="One video slug or one scene name")
    parser.add_argument("--libdir", help="Directory that holds libseams.so")
    parser.add_argument("--pages", action="store_true", help="Write the Sphinx figure list")
    parser.add_argument("--force", action="store_true", help="Draw again even if a PNG is already large enough")
    parser.add_argument("--engine", choices=("cycles", "tachyon"), default="cycles")
    parser.add_argument("--check-colour", action="store_true", help="Render one sphere and check deltaE")
    args = parser.parse_args(argv)
    global RENDER_ENGINE
    RENDER_ENGINE = args.engine
    if args.check_colour:
        if args.engine != "cycles":
            fail("the colour check requires --engine cycles")
        for token in ("hexagonal", "interfacial", "cubic"):
            check_colour(token)
        return 0
    if not args.out:
        fail("arguments are required: --out")
    out = Path(args.out)
    library = load_api(find_lib(args.libdir))
    if args.slug and args.slug not in VIDEO_SLUGS:
        scenes = (args.slug,)
        plan = {args.slug: (args.slug,)}
    elif args.slug:
        plan = {args.slug: VIDEO_SLUGS[args.slug]}
        scenes = VIDEO_SLUGS[args.slug]
    else:
        plan = dict(VIDEO_SLUGS)
        for name in (
            "ice-ih", "ice-ih-h", "ice-ic", "rings", "prisms", "cage-si", "cage-sii", "cage-sh",
            "cage-guest", "brine", "interface", "nucleation", "genice", "solvis",
        ):
            plan.setdefault(name, (name,))
        scenes = []
        for names in plan.values():
            for name in names:
                if name not in scenes:
                    scenes.append(name)
    cache = {}
    for scene in scenes:
        targets = [slug for slug, names in plan.items() if scene in names]
        pending = []
        for slug in targets:
            names = plan[slug]
            stem = "overlay" if names[0] == scene else scene
            if args.force or fresh(out / slug, stem):
                pending.append((slug, stem))
        if not pending:
            continue
        print("render", scene, flush=True)
        cache[scene] = draw(library, scene)
        for slug, stem in pending:
            save_pair(cache[scene], out / slug, stem)
            print("wrote", out / slug / "{}.png".format(stem), flush=True)
    if args.pages:
        write_pages(out)
    return 0


if __name__ == "__main__":
    if os.environ.get("MOLECULES_CYCLES") == "1":
        raise SystemExit(cycles_job())
    raise SystemExit(main(sys.argv[1:]))
