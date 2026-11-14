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
sys.path.insert(0, str(ROOT / "scripts"))
from lammps_box import cell_from_bounds


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
COLUMN_MM = 85.0
# 8 pt at 85 mm on the 2048-wide still. The scale bar and the legend share it.
CAPTION_PX = int(round(8.0 / 72.0 * (WIDTH / (COLUMN_MM / 25.4))))
CAPTION_MARGIN = 128
RENDER_SCALE = max(1, int(os.environ.get("MOLECULES_SCALE", "2")))
WEB = (960, 720)

# Colours come from scripts/chill_palette.py. Guest colours stay here and
# must not repeat a palette triple.
from chill_palette import COLOUR as _SHARED


def _hex(colour):
    return "#{:02X}{:02X}{:02X}".format(*(int(round(channel * 255)) for channel in colour))


PALETTE = {name: _hex(colour) for name, colour in _SHARED.items()}
PALETTE.update({
    "methane": "#F9A825",
    "THF": "#EF6C00",
    "hydrogen": "#90A4AE",
    # ColorBrewer Dark2 purple. It stays clear of the interfacial yellow.
    "mixed": "#7570B3",
    "ring O": _hex(_SHARED["red"]),
})
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


# One framework radius and one bond width for every still.
ATOM_RADIUS = 0.50
BOND_RADIUS = 0.06


def atom_style(token):
    """Colour, radius and transparency. Framework oxygen is one radius."""
    colour = palette_rgb(token)
    if token == "water":
        return colour, 0.36, 0.72
    if token == "hydrogen":
        return colour, 0.22, 0.0
    if token == "sodium":
        return colour, 0.90, 0.0
    if token == "chloride":
        return colour, 1.10, 0.0
    if token == "ring O":
        return colour, ATOM_RADIUS, 0.0
    return colour, ATOM_RADIUS, 0.0

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
        print("chemparseplot, ovito or ase is not installed", file=sys.stderr)
        raise SystemExit(77)
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
    print("chemparseplot, ovito or ase is not installed", file=sys.stderr)
    raise SystemExit(77)


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
        i += 1
        rows = []
        for _ in range(3):
            rows.append([float(item) for item in lines[i].split()])
            i += 1
        origin, edges, tilts = cell_from_bounds(rows)
        box9 = [*origin, *edges, *tilts]
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


def _lattice(box9):
    """Columns are the LAMMPS cell vectors, including tilt."""
    import numpy as np

    length_x, length_y, length_z = box9[3:6]
    tilt_xy, tilt_xz, tilt_yz = box9[6:9]
    return np.array(
        [
            [length_x, tilt_xy, tilt_xz],
            [0.0, length_y, tilt_yz],
            [0.0, 0.0, length_z],
        ],
        dtype=float,
    )


def unwrap_around(points, centre, box9):
    """Bring every atom to the image nearest the centre of mass.

    A wrapped cluster otherwise projects a neighbour onto the far side of
    the cell, and the bond is drawn across the panel.
    """
    import numpy as np

    lattice = _lattice(box9)
    delta = np.asarray(points, dtype=float) - np.asarray(centre, dtype=float)
    fractional = delta @ np.linalg.inv(lattice).T
    fractional -= np.round(fractional)
    return np.asarray(centre, dtype=float) + fractional @ lattice.T


def crystal_bonds(points, labels):
    """O–O contacts under 3.5 Å. Each water keeps its four shortest.

    The points are already unwrapped, so the cartesian length is the
    minimum-image length. A longer segment is not drawn.
    """
    import numpy as np

    crystal = {0, 1, 7, 8}
    points = np.asarray(points, dtype=float)
    ids = [index for index, label in enumerate(labels) if int(label) in crystal]
    found = []
    for cursor, left in enumerate(ids):
        for right in ids[cursor + 1 :]:
            span = float(np.linalg.norm(points[left] - points[right]))
            if 2.30 <= span < 3.5:
                found.append((span, left, right))
    degree = {index: [] for index in ids}
    for span, left, right in found:
        degree[left].append((span, right))
        degree[right].append((span, left))
    dropped = set()
    for index, neighbours in degree.items():
        if len(neighbours) <= 4:
            continue
        for _span, other in sorted(neighbours)[4:]:
            dropped.add(tuple(sorted((index, other))))
    bonds = [
        (left, right)
        for span, left, right in found
        if tuple(sorted((left, right))) not in dropped
    ]
    count = {index: 0 for index in ids}
    for left, right in bonds:
        span = float(np.linalg.norm(points[left] - points[right]))
        if span >= 3.5:
            fail("bond is {:.2f} Å, longer than 3.5".format(span))
        count[left] += 1
        count[right] += 1
    over = [index for index, number in count.items() if number > 4]
    if over:
        fail("a water has {} bonds; the limit is 4".format(max(count.values())))
    longest = max(
        (float(np.linalg.norm(points[left] - points[right])) for left, right in bonds),
        default=0.0,
    )
    degree_max = max(count.values(), default=0)
    print(
        "nucleus bonds {} longest {:.2f} A degree {}".format(len(bonds), longest, degree_max),
        flush=True,
    )
    return bonds


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


def _rotation_taking(source, target):
    """Matrix that sends the column vector source onto target."""
    import numpy as np

    source = np.asarray(source, dtype=float)
    target = np.asarray(target, dtype=float)
    source = source / (np.linalg.norm(source) or 1.0)
    target = target / (np.linalg.norm(target) or 1.0)
    cross = np.cross(source, target)
    sine = float(np.linalg.norm(cross))
    cosine = float(np.dot(source, target))
    if sine < 1e-8:
        if cosine > 0.0:
            return np.eye(3)
        axis = np.array([1.0, 0.0, 0.0])
        if abs(float(source[0])) > 0.9:
            axis = np.array([0.0, 1.0, 0.0])
        axis = np.cross(source, axis)
        axis = axis / np.linalg.norm(axis)
        skew = np.array([
            [0.0, -axis[2], axis[1]],
            [axis[2], 0.0, -axis[0]],
            [-axis[1], axis[0], 0.0],
        ])
        return np.eye(3) + 2.0 * (skew @ skew)
    axis = cross / sine
    skew = np.array([
        [0.0, -axis[2], axis[1]],
        [axis[2], 0.0, -axis[0]],
        [-axis[1], axis[0], 0.0],
    ])
    return np.eye(3) + sine * skew + (1.0 - cosine) * (skew @ skew)


def _cluster_mapper(positions, names):
    """Put the stacking axis up the page, with a turn so the bands have depth."""
    import numpy as np

    pts = np.asarray(positions, dtype=float)
    crystal = np.array([
        name in ("cubic", "hexagonal", "reCubic", "reHex") for name in names
    ])
    if int(np.count_nonzero(crystal)) < 8:
        crystal = np.ones(len(pts), dtype=bool)
    chosen = pts[crystal]
    centre = chosen.mean(axis=0)
    delta = chosen - centre
    _values, vectors = np.linalg.eigh(np.cov(delta.T))
    order = np.argsort(_values)
    hex_mark = np.array([
        name in ("hexagonal", "reHex")
        for name, keep in zip(names, crystal) if keep
    ], dtype=float)
    best_score = -1.0
    stacking = vectors[:, int(order[0])]
    for col in range(3):
        axis = vectors[:, col]
        seq = hex_mark[np.argsort(delta @ axis)]
        if len(seq) < 6:
            score = 0.0
        else:
            smooth = np.convolve(seq, np.ones(5) / 5.0, mode="valid")
            score = float(np.sum(np.abs(np.diff(smooth)) > 0.2))
        if col == int(order[0]):
            score += 1.0
        if score > best_score:
            best_score = score
            stacking = axis
    stacking = stacking / (np.linalg.norm(stacking) or 1.0)
    perp = delta - np.outer(delta @ stacking, stacking)
    spread, plane = np.linalg.eigh(np.cov(perp.T))
    across = plane[:, int(np.argmax(spread))]
    across = across - stacking * float(np.dot(across, stacking))
    across = across / (np.linalg.norm(across) or 1.0)
    depth = np.cross(across, stacking)
    depth = depth / (np.linalg.norm(depth) or 1.0)
    basis = np.column_stack((across, stacking, depth))
    if np.linalg.det(basis) < 0.0:
        basis[:, 2] *= -1.0
    # row' = row @ spin matches _spin. Yaw, then a small tilt.
    yaw, tilt = np.deg2rad(28.0), np.deg2rad(16.0)
    cy, sy = np.cos(yaw), np.sin(yaw)
    ct, st = np.cos(tilt), np.sin(tilt)
    spin_y = np.array([[cy, 0.0, -sy], [0.0, 1.0, 0.0], [sy, 0.0, cy]])
    spin_x = np.array([[1.0, 0.0, 0.0], [0.0, ct, st], [0.0, -st, ct]])
    matrix = basis @ spin_y @ spin_x

    def apply(points):
        return (np.asarray(points, dtype=float) - centre) @ matrix

    return apply


def _depth_cue(positions, radii, transparency, names):
    """Shrink and fade atoms toward the back of the cluster.

    Cubic and hexagonal stay a little larger than the Ice Ih oxygen.
    Mixed atoms stay smaller. The renderer treats mid alpha as glass,
    so the fade steps from solid to a flat ghost.
    """
    import numpy as np

    z = np.asarray(positions, dtype=float)[:, 2]
    span = float(z.max() - z.min()) or 1.0
    near = (z - float(z.min())) / span
    radii = np.asarray(radii, dtype=float).copy()
    transparency = np.asarray(transparency, dtype=float).copy()
    for index, name in enumerate(names):
        if name in ("cubic", "hexagonal", "reCubic", "reHex"):
            base = ATOM_RADIUS * 1.12
        elif name in ("mixed", "interfacial"):
            base = ATOM_RADIUS * 0.64
        elif name == "water":
            base = ATOM_RADIUS * 0.40
        else:
            continue
        radii[index] = base * (0.58 + 0.42 * float(near[index]))
        transparency[index] = 0.0 if float(near[index]) >= 0.50 else 0.78
    return radii, transparency, near


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
    room_x = width / 2.0 - 48.0
    room_y = height / 2.0 - 48.0
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
    if view == "slab":
        # One depth (x). The interface normal z runs across the page and the
        # c axis y runs up it, so the ice reads as columns.
        return np.column_stack((pos[:, 2], pos[:, 1], -pos[:, 0]))
    if view == "pair":
        # Cages already sit side by side. A small tilt shows the guest inside.
        return _spin(pos, "x", 20.0)
    if view == "row":
        # sH cages are posed per cage, including a pentagon facing the camera.
        return pos
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


def camera_fit(positions, radii, meshes, panel, names=None):
    """Ortho camera for the panels whose framing is checked.

    None leaves the historical fit alone, so Ice Ih, sII, the interface
    and the brine slab are not moved.
    """
    import numpy as np

    panel = str(panel)
    if not (panel.startswith("cage-sh") or panel.startswith("nucleation")):
        return None
    points = np.asarray(positions, dtype=float)
    radius = np.asarray(radii, dtype=float)
    framed = [points]
    for mesh in meshes or []:
        framed.append(np.asarray(mesh["vertices"], dtype=float))
    framed = np.vstack(framed)
    if panel.startswith("cage-sh"):
        # The sphere, and the face that sits just outside it, stay in the box.
        pad = float(radius.max()) + 0.55 if len(radius) else 1.0
        low = framed.min(axis=0) - pad
        high = framed.max(axis=0) + pad
        centre = np.array([
            float(low[0] + high[0]) / 2.0,
            float(low[1] + high[1]) / 2.0,
            float(framed[:, 2].mean()),
        ])
        width = float(high[0] - low[0])
        height = float(high[1] - low[1])
        vertical = max(height, width * HEIGHT / float(WIDTH))
        # The padded box occupies 90% of the view: 5% of the axes on each side.
        return centre, vertical / 0.90
    names = list(names or [])
    crystal = [
        index for index, name in enumerate(names)
        if name in ("cubic", "hexagonal", "reCubic", "reHex")
    ]
    if len(crystal) < 8:
        crystal = list(range(len(points)))
    chosen = points[crystal]
    span_x = float(chosen[:, 0].max() - chosen[:, 0].min())
    span_y = float(chosen[:, 1].max() - chosen[:, 1].min())
    centre = np.array([
        float(chosen[:, 0].mean()),
        float(chosen[:, 1].mean()),
        float(chosen[:, 2].mean()),
    ])
    # The crystal itself fills about 70% of the panel. The caption and the
    # crop take some of the page, so the camera is a little tighter than that.
    ortho = max(span_y / 0.78, (span_x / 0.78) * HEIGHT / float(WIDTH))
    horiz = ortho * WIDTH / float(HEIGHT)
    for _step in range(16):
        left = 0.5 + (points[:, 0] - radius - centre[0]) / horiz
        right = 0.5 + (points[:, 0] + radius - centre[0]) / horiz
        lower = 0.5 + (points[:, 1] - radius - centre[1]) / ortho
        upper = 0.5 + (points[:, 1] + radius - centre[1]) / ortho
        if left.min() >= 0.0 and right.max() <= 1.0 and lower.min() >= 0.0 and upper.max() <= 1.0:
            break
        ortho *= 1.05
        horiz = ortho * WIDTH / float(HEIGHT)
    return centre, ortho


def assert_atoms_inside(positions, radii, centre, ortho_scale, margin=0.05):
    """Fail when a sphere falls outside the camera axes.

    ``margin`` is the fraction of the axes that must stay empty. The sH
    cages use 5%. A nucleus only has to stay inside the axes.
    """
    import numpy as np

    points = np.asarray(positions, dtype=float)
    radius = np.asarray(radii, dtype=float)
    horiz = float(ortho_scale) * WIDTH / float(HEIGHT)
    left = 0.5 + (points[:, 0] - radius - centre[0]) / horiz
    right = 0.5 + (points[:, 0] + radius - centre[0]) / horiz
    lower = 0.5 + (points[:, 1] - radius - centre[1]) / float(ortho_scale)
    upper = 0.5 + (points[:, 1] + radius - centre[1]) / float(ortho_scale)
    if np.any(left < 0.0) or np.any(right > 1.0) or np.any(lower < 0.0) or np.any(upper > 1.0):
        fail("projected atom falls outside the axes")
    edge = min(float(left.min()), float(1.0 - right.max()), float(lower.min()), float(1.0 - upper.max()))
    if edge < float(margin) - 1.0e-6:
        fail("projected atom is {:.1%} from the axes; need {:.0%} margin".format(edge, margin))
    print("axes margin {:.1%}".format(edge), flush=True)


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
    preset = "view_ortho" in job.files
    if preset:
        # The parent already fitted this panel and checked the axes.
        centre = np.asarray(job["view_centre"], dtype=float)
        vertical = float(job["view_ortho"])
        look = centre.copy()
        high = framed.max(axis=0)
    else:
        com = framed.mean(axis=0)
        pad = float(radii.max()) + 1.8 if len(radii) else 1.8
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
            # Matte. A specular peak was blowing the bottom of the spheres out to white.
            rough.default_value = 1.0
        spec = socket(principled, "Specular IOR Level", "Specular")
        if spec is not None:
            spec.default_value = 0.0
        coat = socket(principled, "Coat Weight")
        if coat is not None:
            coat.default_value = 0.0
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
            rough.default_value = 1.0
        if spec is not None:
            spec.default_value = 0.0
        if coat is not None:
            coat.default_value = 0.0
        # Flat palette colour. Lighting was painting a hotspot on a few oxygens.
        emission = tree.nodes.new("ShaderNodeEmission")
        emission.inputs["Color"].default_value = (*linear, 1.0)
        emission.inputs["Strength"].default_value = 1.0
        tree.links.new(emission.outputs["Emission"], output.inputs["Surface"])
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
        # Radius and alpha stay per atom so a depth cue can shrink the back.
        key = (
            tuple(round(float(channel), 4) for channel in colour),
            round(float(radii[index]), 3),
            round(float(1.0 - transparency[index]), 2),
        )
        groups.setdefault(key, []).append(index)
    for key, members in groups.items():
        colour, radius, alpha = key
        # Faint context (the liquid shell) stays a flat ghost. Mid alpha is glass.
        transmission = 0.55 if 0.22 < alpha < 0.75 else 0.0
        material = glossy("atom", colour, alpha, transmission)
        from mathutils import Matrix

        centres = positions[members]
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
                bond_rough.default_value = 1.0
            bond_spec = socket(rough, "Specular IOR Level", "Specular")
            if bond_spec is not None:
                bond_spec.default_value = 0.0
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
            mid = np.array(mid, dtype=float, copy=True)
            # Sit the bond just behind the atoms so a flat sheet does not z-fight.
            mid[2] -= 0.08
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
        add_bond_set(bold, bond_radius, "cage-bonds")
    elif len(bonds):
        add_bond_set(bonds, bond_radius, "bonds")
    for slot in range(mesh_count):
        vertices = np.asarray(job["mesh_vertices_{}".format(slot)], dtype=float)
        faces = [tuple(int(v) for v in face) for face in job["mesh_faces_{}".format(slot)]]
        colour = tuple(float(v) for v in job["mesh_color_{}".format(slot)])
        face = 0.42
        key = "mesh_alpha_{}".format(slot)
        if key in job.files:
            face = float(job[key])
        material = glossy("cage", colour, face, 0.35)
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
        if not preset:
            vertical = vertical / 0.86
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

    # One large lamp on the camera axis. Every sphere then has the same shade.
    area("key", (float(centre[0]), float(centre[1]), float(high[2]) + 28.0), 700.0, 48.0)

    world = bpy.data.worlds.new("studio")
    bpy.context.scene.world = world
    world.use_nodes = True
    backdrop = world.node_tree.nodes.get("Background")
    if backdrop is not None:
        backdrop.inputs["Color"].default_value = (*background, 1.0)
        backdrop.inputs["Strength"].default_value = 0.85 if light else 0.2

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


def cycles_launch(positions, colors, radii, transparency, bonds, meshes, background, panel="panel", samples=None, ortho=False, bold=None, bond_radius=0.045, bond_rgb=None, view_centre=None, view_ortho=None):
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
    if view_centre is not None and view_ortho is not None:
        payload["view_centre"] = np.asarray(view_centre, dtype=np.float64)
        payload["view_ortho"] = np.float64(view_ortho)
    if bold is not None and len(bold):
        payload["bold_bonds"] = np.asarray(bold, dtype=np.int32)
    for slot, mesh in enumerate(meshes):
        payload["mesh_vertices_{}".format(slot)] = np.asarray(mesh["vertices"], dtype=np.float64)
        payload["mesh_faces_{}".format(slot)] = np.asarray(mesh["faces"], dtype=np.int32)
        payload["mesh_color_{}".format(slot)] = np.asarray(mesh["color"], dtype=np.float64)
        payload["mesh_alpha_{}".format(slot)] = np.float64(mesh.get("alpha", 0.42))
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


def raytrace(positions, colors, radii, transparency, bonds, meshes, background, perspective, ortho=False, bold=None, bond_radius=0.045, bond_rgb=None, view_centre=None, view_ortho=None):
    """Cycles when requested. Tachyon is only used when the engine is tachyon."""
    if RENDER_ENGINE == "cycles":
        drawn = cycles_launch(
            positions, colors, radii, transparency, bonds, meshes, background,
            panel=PANEL_NAME, ortho=ortho, bold=bold,
            bond_radius=bond_radius, bond_rgb=bond_rgb,
            view_centre=view_centre, view_ortho=view_ortho,
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


def _ink_mask(picture, light):
    """Atoms, including a highlight. Alpha is the mask when the render has one."""
    import numpy as np

    if picture.mode == "RGBA":
        alpha = np.asarray(picture.split()[-1])
        return alpha > 16
    array = np.asarray(picture.convert("RGB"))
    if light:
        return np.any(array < 246, axis=2)
    return np.any(array > 40, axis=2)


def crop_margin(picture, light):
    import numpy as np
    from PIL import Image

    array = np.asarray(picture)
    mask = _ink_mask(picture, light)
    rows = np.where(mask.any(axis=1))[0]
    cols = np.where(mask.any(axis=0))[0]
    if len(rows) == 0:
        return picture, picture.size[1], picture.size[1]
    # Keep 5% of the render, so a sphere that is inside the axes is not
    # cropped back to the ink and sliced on the page.
    pad_x = max(28, int(round(array.shape[1] * 0.05)))
    pad_y = max(28, int(round(array.shape[0] * 0.05)))
    top = max(0, int(rows[0]) - pad_y)
    bottom = min(array.shape[0], int(rows[-1]) + pad_y)
    left = max(0, int(cols[0]) - pad_x)
    right = min(array.shape[1], int(cols[-1]) + pad_x)
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

    mask = _ink_mask(picture, light)
    plate = Image.fromarray((mask.astype(np.uint8) * 255))
    plate = plate.filter(ImageFilter.GaussianBlur(radius=7))
    shade = (186, 188, 190) if light else (12, 12, 14)
    shadow = Image.new("RGB", picture.size, shade)
    return shadow, plate


def _pack_labels(draw, items, font, span, swatch, gap):
    """Rows that stay inside ``span``. None when one label is wider than the span."""
    rows = []
    row = []
    used = 0.0
    for name, colour in items:
        text_w = float(draw.textlength(name, font=font))
        width = float(swatch + 16 + text_w)
        if width > span:
            return None
        if row and used + gap + width > span:
            rows.append(row)
            row = []
            used = 0.0
        if row:
            used += gap
        row.append((name, colour, text_w, width))
        used += width
    if row:
        rows.append(row)
    return _balance_rows(rows, span, gap)


def _row_width(row, gap):
    if not row:
        return 0.0
    return sum(item[3] for item in row) + gap * (len(row) - 1)


def _balance_rows(rows, span, gap):
    """Keep row lengths within one entry when the shorter row still fits."""
    guard = 0
    while len(rows) >= 2 and guard < 8:
        guard += 1
        longest = max(range(len(rows)), key=lambda index: len(rows[index]))
        shortest = min(range(len(rows)), key=lambda index: len(rows[index]))
        if len(rows[longest]) - len(rows[shortest]) < 2:
            break
        if shortest != longest + 1 and shortest != longest - 1:
            break
        if shortest == longest + 1:
            moved = rows[longest][-1]
            donor = rows[longest][:-1]
            taken = [moved, *rows[shortest]]
        else:
            moved = rows[longest][0]
            donor = rows[longest][1:]
            taken = [*rows[shortest], moved]
        if not donor or _row_width(taken, gap) > span:
            break
        rows[longest] = donor
        rows[shortest] = taken
    return rows


def _legend_layout(items):
    """Caption at 8 pt when the still is printed 85 mm wide.

    The scale bar label uses this same font. Rows wrap so a label is not clipped.
    """
    from PIL import Image, ImageDraw

    probe = Image.new("RGB", (WIDTH, 400), (255, 255, 255))
    draw = ImageDraw.Draw(probe)
    margin = CAPTION_MARGIN
    full = WIDTH - 2 * margin
    size = CAPTION_PX
    font = _font(size)
    swatch = int(round(size * 0.62))
    rows = _pack_labels(draw, items, font, full, swatch, 36.0)
    if rows is None or len(rows) > 3:
        rows = _pack_labels(draw, items, font, full, swatch, 18.0) or [[]]
    return font, size, swatch, rows


def compose(structure, legend, world_height, light):
    """Structure in the drawing band. Scale bar and legend stay in the caption."""
    from PIL import Image, ImageDraw
    import numpy as np

    bg = (255, 255, 255) if light else (22, 24, 28)
    ink = (25, 25, 28) if light else (236, 237, 239)
    picture, crop_rows, render_rows = crop_margin(structure, light)
    items = []
    for name, colour_rgb in legend:
        colour = tuple(int(max(0, min(255, round(channel * 255)))) for channel in colour_rgb)
        items.append((name, colour))
    font, font_size, swatch, rows = _legend_layout(items)
    margin = CAPTION_MARGIN
    bar_band = font_size + 28
    row_h = font_size + 26
    caption = 20 + bar_band + 8 + row_h * max(1, len(rows)) + 16
    top = 28
    gap = 72
    avail_w = WIDTH - 2 * margin
    avail_h = HEIGHT - caption - top - gap
    scale = min(avail_w / float(picture.size[0]), avail_h / float(picture.size[1]))
    sized = (
        max(1, int(round(picture.size[0] * scale))),
        max(1, int(round(picture.size[1] * scale))),
    )
    picture = picture.resize(sized, Image.Resampling.LANCZOS)
    left = (WIDTH - sized[0]) // 2
    region_top = top
    region_bottom = HEIGHT - caption - gap
    struct_y = region_top + max(0, (region_bottom - region_top - sized[1]) // 2)
    # Three cages in a row are short. Sit them higher so the row is not
    # stranded under a band of empty page.
    if str(PANEL_NAME).startswith("cage-sh"):
        struct_y = max(region_top + 12, struct_y - 48)

    def structure_layer(origin_y):
        layer = Image.new("RGBA", (WIDTH, HEIGHT), (0, 0, 0, 0))
        shadow, plate = _contact_shadow(picture, light)
        tint = Image.new("RGBA", picture.size, (186, 188, 190, 70) if light else (8, 8, 10, 90))
        tint.putalpha(plate.point(lambda value: int(value * 0.45)))
        layer.paste(tint, (left, origin_y + 10), tint)
        sprite = picture.convert("RGBA")
        sprite.putalpha(Image.fromarray((_ink_mask(picture, light).astype(np.uint8) * 255)))
        layer.paste(sprite, (left, origin_y), sprite)
        return layer

    canvas = Image.new("RGB", (WIDTH, HEIGHT), bg)
    draw = ImageDraw.Draw(canvas)
    cap_top = HEIGHT - caption
    bar = 0.0
    if world_height > 0 and crop_rows > 0 and render_rows > 0:
        world_crop = world_height * (crop_rows / float(render_rows))
        bar = 10.0 * sized[1] / world_crop
    bar = min(bar, WIDTH * 0.42)
    bar_y = cap_top + 8 + font_size // 2
    if bar > 8:
        draw.line((margin, bar_y, margin + bar, bar_y), fill=ink, width=8)
        draw.line((margin, bar_y - 12, margin, bar_y + 12), fill=ink, width=4)
        draw.line((margin + bar, bar_y - 12, margin + bar, bar_y + 12), fill=ink, width=4)
        draw.text((margin + bar + 14, bar_y - font_size * 0.55), "10 \u00c5", fill=ink, font=font)
    label_y0 = cap_top + bar_band + 4
    span = WIDTH - 2 * margin
    for row_index, row in enumerate(rows):
        if not row:
            continue
        y = label_y0 + row_index * row_h
        used = sum(item[3] for item in row)
        spare = max(0.0, span - used)
        gap_px = spare / (len(row) - 1) if len(row) > 1 else 0.0
        # A short row stays together. Stretching three words across the page
        # leaves the swatches looking lost.
        if gap_px > 72:
            gap_px = 36
        x = float(margin)
        if len(row) == 1:
            x = float(WIDTH - margin - row[0][3])
        for name, colour, _text_w, width in row:
            draw.rounded_rectangle(
                (x, y + 6, x + swatch, y + 6 + swatch),
                radius=6, fill=colour, outline=ink, width=2,
            )
            draw.text((x + swatch + 14, y + 2), name, fill=ink, font=font)
            x += width + gap_px
    # Nudge only the drawing so the ink centroid stays near the middle,
    # and stop before the atoms meet the caption.
    placed = struct_y
    best = structure_layer(placed)
    base = canvas.convert("RGBA")
    merged = Image.alpha_composite(base, best).convert("RGB")
    arr = np.asarray(merged)
    mask = np.any(arr < 248, axis=2)
    ys, xs = np.where(mask)
    if len(ys):
        dy = float(ys.mean() / HEIGHT - 0.5)
        shift = int(round(-dy * HEIGHT))
        floor = region_top
        ceiling = region_bottom - sized[1]
        shifted = min(max(placed + shift, floor), max(floor, ceiling))
        if shifted != placed:
            placed = shifted
            best = structure_layer(placed)
            merged = Image.alpha_composite(base, best).convert("RGB")
    # The row is short. A small lift keeps it from sitting on the caption
    # while the horizontal margin stays put.
    if str(PANEL_NAME).startswith("cage-sh"):
        lifted = max(region_top, placed - 28)
        if lifted != placed:
            best = structure_layer(lifted)
            merged = Image.alpha_composite(base, best).convert("RGB")
    return merged


COLUMN_DPI = 300


def write_vector(directory, stem, geom):
    """Orthographic circles and lines, plus a scale bar and the legend.

    The PDF is one column wide (85 mm) so the caption is the print size.
    """
    import numpy as np

    pos = np.asarray(geom["positions"], dtype=float)
    colours = np.asarray(geom["colors"], dtype=float)
    radii = np.asarray(geom["radii"], dtype=float)
    bonds = geom["bonds"]
    legend = list(geom.get("legend") or [])
    if len(pos) == 0:
        return
    order = np.argsort(pos[:, 2])
    pad = float(radii.max()) + 0.4
    xmin, xmax = float(pos[:, 0].min() - pad), float(pos[:, 0].max() + pad)
    ymin, ymax = float(pos[:, 1].min() - pad), float(pos[:, 1].max() + pad)
    width = max(xmax - xmin, 1.0)
    height = max(ymax - ymin, 1.0)
    # 8 pt at 85 mm, expressed in the angstrom viewBox.
    font_world = 8.0 / 72.0 * 25.4 * width / COLUMN_MM
    caption = font_world * (2.4 + 1.5 * max(1, (len(legend) + 2) // 3))
    svg = [
        '<?xml version="1.0" encoding="UTF-8"?>',
        '<svg xmlns="http://www.w3.org/2000/svg" viewBox="0 0 {:.3f} {:.3f}">'.format(width, height + caption),
        '<rect width="100%" height="100%" fill="white"/>',
    ]
    for left, right in bonds:
        a, b = pos[int(left)], pos[int(right)]
        svg.append(
            '<line x1="{:.3f}" y1="{:.3f}" x2="{:.3f}" y2="{:.3f}" stroke="#555860" stroke-width="{:.3f}" stroke-linecap="round"/>'.format(
                a[0] - xmin, ymax - a[1], b[0] - xmin, ymax - b[1], BOND_RADIUS * 2.0,
            )
        )
    for index in order:
        red, green, blue = [int(round(max(0.0, min(1.0, float(channel))) * 255)) for channel in colours[index]]
        svg.append(
            '<circle cx="{:.3f}" cy="{:.3f}" r="{:.3f}" fill="#{:02x}{:02x}{:02x}"/>'.format(
                pos[index, 0] - xmin, ymax - pos[index, 1], float(radii[index]), red, green, blue,
            )
        )
    bar_y = height + font_world * 1.1
    svg.append(
        '<line x1="{:.3f}" y1="{:.3f}" x2="{:.3f}" y2="{:.3f}" stroke="#19191c" stroke-width="{:.3f}"/>'.format(
            0.15, bar_y, 10.15, bar_y, font_world * 0.08,
        )
    )
    svg.append(
        '<text x="{:.3f}" y="{:.3f}" font-size="{:.3f}" font-family="DejaVu Sans, sans-serif" fill="#19191c">10 &#197;</text>'.format(
            10.4, bar_y + font_world * 0.35, font_world,
        )
    )
    x_cursor = 0.15
    y_text = bar_y + font_world * 1.8
    for name, colour in legend:
        red, green, blue = [int(round(max(0.0, min(1.0, float(channel))) * 255)) for channel in colour]
        svg.append(
            '<rect x="{:.3f}" y="{:.3f}" width="{:.3f}" height="{:.3f}" rx="0.08" fill="#{:02x}{:02x}{:02x}" stroke="#19191c" stroke-width="0.03"/>'.format(
                x_cursor, y_text - font_world * 0.8, font_world * 0.7, font_world * 0.7, red, green, blue,
            )
        )
        svg.append(
            '<text x="{:.3f}" y="{:.3f}" font-size="{:.3f}" font-family="DejaVu Sans, sans-serif" fill="#19191c">{}</text>'.format(
                x_cursor + font_world * 0.9, y_text, font_world, name,
            )
        )
        x_cursor += font_world * (1.6 + 0.55 * len(name))
        if x_cursor > width - font_world * 4:
            x_cursor = 0.15
            y_text += font_world * 1.5
    svg.append("</svg>")
    (directory / "{}.svg".format(stem)).write_text("\n".join(svg) + "\n")
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    from matplotlib.patches import Circle, FancyBboxPatch
    fig_w = COLUMN_MM / 25.4
    fig_h = fig_w * (height + caption) / width
    fig, ax = plt.subplots(figsize=(fig_w, fig_h))
    ax.set_aspect("equal")
    ax.set_xlim(xmin, xmax)
    ax.set_ylim(ymin - caption, ymax)
    ax.axis("off")
    for left, right in bonds:
        a, b = pos[int(left)], pos[int(right)]
        ax.plot([a[0], b[0]], [a[1], b[1]], color="#555860", lw=0.7, solid_capstyle="round", zorder=1)
    for index in order:
        ax.add_patch(Circle(
            (pos[index, 0], pos[index, 1]), float(radii[index]),
            facecolor=tuple(float(c) for c in colours[index]), edgecolor="none", zorder=2,
        ))
    bar_y = ymin - font_world * 1.2
    ax.plot([xmin + 0.15, xmin + 10.15], [bar_y, bar_y], color="#19191c", lw=1.1, solid_capstyle="butt")
    ax.plot([xmin + 0.15, xmin + 0.15], [bar_y - 0.25, bar_y + 0.25], color="#19191c", lw=0.8)
    ax.plot([xmin + 10.15, xmin + 10.15], [bar_y - 0.25, bar_y + 0.25], color="#19191c", lw=0.8)
    ax.text(xmin + 10.4, bar_y, "10 \u00c5", va="center", ha="left", fontsize=8, color="#19191c")
    x_cursor = xmin + 0.15
    y_text = bar_y - font_world * 1.6
    for name, colour in legend:
        ax.add_patch(FancyBboxPatch(
            (x_cursor, y_text - font_world * 0.15), font_world * 0.55, font_world * 0.55,
            boxstyle="round,pad=0.02,rounding_size=0.08",
            facecolor=tuple(float(c) for c in colour), edgecolor="#19191c", linewidth=0.4,
        ))
        ax.text(x_cursor + font_world * 0.75, y_text + font_world * 0.12, name, va="center", ha="left", fontsize=8, color="#19191c")
        x_cursor += font_world * (1.5 + 0.52 * len(name))
        if x_cursor > xmax - font_world * 5:
            x_cursor = xmin + 0.15
            y_text -= font_world * 1.4
    fig.savefig(directory / "{}.pdf".format(stem), bbox_inches="tight")
    plt.close(fig)


def save_pair(image, directory, stem):
    from PIL import Image

    directory.mkdir(parents=True, exist_ok=True)
    geom = None
    if isinstance(image, tuple) and len(image) == 3:
        light, dark, geom = image
    else:
        light, dark = image if isinstance(image, tuple) else (image, None)
    # 2048 px across 85 mm is 612 dpi. The column file is exactly 300 dpi.
    light.save(directory / "{}.png".format(stem), dpi=(612, 612))
    column_w = int(round(COLUMN_MM / 25.4 * COLUMN_DPI))
    column_h = int(round(column_w * light.size[1] / light.size[0]))
    light.resize((column_w, column_h), Image.Resampling.LANCZOS).save(
        directory / "{}-300.png".format(stem), dpi=(COLUMN_DPI, COLUMN_DPI),
    )
    light.resize(WEB, Image.Resampling.LANCZOS).save(directory / "{}-web.png".format(stem))
    if dark is not None:
        dark.save(directory / "{}-dark.png".format(stem), dpi=(612, 612))
        dark.resize(WEB, Image.Resampling.LANCZOS).save(directory / "{}-web-dark.png".format(stem))
    if geom is not None:
        write_vector(directory, stem, geom)
    return directory / "{}.png".format(stem)


def largest_cluster(frame, labels, shell=3.6, core=(0, 1, 3, 7, 8)):
    """The biggest crystal component, plus waters that touch it."""
    import numpy as np
    from scipy.spatial import cKDTree

    crystal = np.array([i for i, label in enumerate(labels) if label in core], dtype=int)
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
        minimum = 16 if kind == "CH4" else 24
        chosen = (distance > 2.4) & (distance < limit)
        # The sH large cage sits further from its THF than an sII cage does.
        if kind != "CH4" and int(np.count_nonzero(chosen)) < minimum:
            limit = 8.8
            minimum = 18
            chosen = (distance > 2.4) & (distance < limit)
        if int(np.count_nonzero(chosen)) < minimum:
            continue
        full = centre + delta[chosen]
        try:
            hull = ConvexHull(full, qhull_options="QJ")
        except Exception:
            continue
        # Interior waters are the next shell. They read as loose balls inside the cage.
        surface = np.asarray(hull.vertices, dtype=int)
        remap = {int(old): new for new, old in enumerate(surface)}
        faces = [[remap[int(vertex)] for vertex in tri] for tri in hull.simplices]
        shell = full[surface]
        token = "hexagonal" if kind == "CH4" else "clathrate"
        found.append((float(np.linalg.norm(centre - origin)), kind, shell, np.where(chosen)[0][surface], {
            "vertices": shell,
            "faces": np.asarray(faces, dtype=int),
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
            mapper = _cluster_mapper(positions, list(names))
            positions = mapper(positions)
            if str(PANEL_NAME).startswith("nucleation"):
                for left, right in bonds:
                    span = float(np.linalg.norm(positions[int(left)] - positions[int(right)]))
                    if span >= 3.5:
                        fail("projected bond is {:.2f} Å, longer than 3.5".format(span))
            posed = []
            for mesh in meshes:
                posed.append({
                    "vertices": mapper(mesh["vertices"]),
                    "faces": mesh["faces"],
                    "color": mesh["color"],
                    "transparency": mesh["transparency"],
                    "legend": mesh.get("legend"),
                })
            meshes = posed
            radii, transparency, _near = _depth_cue(positions, radii, transparency, names)
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
            # True orthographic view down c. One bilayer, flattened so every
            # oxygen takes the same light and a ring stays a whole hexagon.
            positions[:, 2] = 0.0
            for mesh in meshes:
                mesh["vertices"][:, 2] = 0.0
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
        if view == "cage" and name in ("methane", "THF", "ring O"):
            positions[index, 2] += 0.9
        elif name == "_rim":
            positions[index, 2] -= 0.25
    ortho = view in ("basal", "channel", "slab", "pair", "row") or bool(align)
    fit = camera_fit(positions, radii, meshes, PANEL_NAME, names)
    view_centre = view_ortho = None
    if fit is not None:
        view_centre, view_ortho = fit
        # Cages need 5% of the axes clear. The nucleus only has to stay inside.
        margin = 0.05 if str(PANEL_NAME).startswith("cage") else 0.0
        assert_atoms_inside(positions, radii, view_centre, view_ortho, margin)
    light_picture, world = raytrace(
        positions, colors, radii, transparency, bonds, meshes, (1, 1, 1), perspective,
        ortho=ortho, bold=bold, bond_radius=bond_radius, bond_rgb=bond_rgb,
        view_centre=view_centre, view_ortho=view_ortho,
    )
    if os.environ.get("MOLECULES_THEMES") == "light":
        dark_picture = light_picture
    else:
        dark_picture, _world = raytrace(
            positions, colors, radii, transparency, bonds, meshes, (0.09, 0.10, 0.12), perspective,
            ortho=ortho, bold=bold, bond_radius=bond_radius, bond_rgb=bond_rgb,
            view_centre=view_centre, view_ortho=view_ortho,
        )
    return (
        compose(light_picture, legend, world, True),
        compose(dark_picture, legend, world, False),
        {
            "positions": np.asarray(positions, dtype=float),
            "colors": np.asarray(colors, dtype=float),
            "radii": np.asarray(radii, dtype=float),
            "bonds": list(bonds),
            "legend": [(name, tuple(colour)) for name, colour in legend],
        },
    )


def labelled_picture(frame, labels, view, liquid_clear, ions=False, perspective=False, meshes=None, center=None, cluster=False, slab=False):
    import numpy as np

    colors, names, radii, transparency, heavy = [], [], [], [], []
    kept_xyz = []
    kept_labels = []
    for index in range(len(frame.xyz)):
        kind = int(frame.types[index])
        if ions and kind in (3, 4):
            # brine.rst: type 3 is chloride and type 4 is sodium.
            token = "chloride" if kind == 3 else "sodium"
            colour = palette_rgb(token)
            # Shannon radii relative to the oxygen used for the ice: Cl 1.81, Na 1.02.
            radius = ATOM_RADIUS * (1.8 if token == "chloride" else 1.0)
            colors.append(colour)
            names.append(token)
            radii.append(radius)
            transparency.append(0.0)
            heavy.append(False)
            kept_xyz.append(frame.xyz[index])
            kept_labels.append(int(labels[index]))
            continue
        label = int(labels[index])
        if cluster and label not in (0, 1, 2, 3, 7, 8):
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
    bonds = [] if cluster else oxygen_bonds(xyz, heavy)
    if cluster:
        # Cubic and hexagonal only. Mixed stays unbound, and liquid has no bonds.
        # Unwrap before the bonds are chosen, so a periodic image is not drawn
        # as a line across the cluster.
        crystal = {0, 1, 7, 8}
        from scipy.spatial import cKDTree

        core = [index for index, label in enumerate(kept_labels) if int(label) in crystal]
        if core:
            xyz = unwrap_around(xyz, xyz[core].mean(axis=0), frame.box9)
        keep = set(core)
        if core:
            tree = cKDTree(xyz[core])
            for index, label in enumerate(kept_labels):
                # Liquid is left out. A mixed oxygen is kept only when it
                # touches the cubic or hexagonal lattice.
                if int(label) != 3:
                    continue
                dist, _near = tree.query(xyz[index], k=1)
                if float(dist) <= 3.0:
                    keep.add(index)
        if keep and len(keep) < len(xyz):
            order = sorted(keep)
            xyz = xyz[order]
            colors = [colors[index] for index in order]
            names = [names[index] for index in order]
            radii = [radii[index] for index in order]
            transparency = [transparency[index] for index in order]
            kept_labels = [kept_labels[index] for index in order]
        bonds = crystal_bonds(xyz, kept_labels)
    elif slab:
        # Ice and the interfacial layer are bonded. Liquid stays as atoms
        # so a crystal bond does not run off into the solvent.
        bonds = [
            pair for pair in bonds
            if int(kept_labels[pair[0]]) != 2 and int(kept_labels[pair[1]]) != 2
        ]
    elif liquid_clear:
        # An interfacial oxygen must not grow a bond into a liquid that is not drawn.
        bonds = [
            pair for pair in bonds
            if int(kept_labels[pair[0]]) != 2 and int(kept_labels[pair[1]]) != 2
        ]
        keep = set(range(len(xyz)))
        for _pass in range(6):
            degree = {}
            for left, right in bonds:
                if left in keep and right in keep:
                    degree[left] = degree.get(left, 0) + 1
                    degree[right] = degree.get(right, 0) + 1
            drop = {index for index, count in degree.items() if count < 2}
            if not drop:
                break
            keep -= drop
            bonds = [pair for pair in bonds if pair[0] in keep and pair[1] in keep]
        linked = set()
        for left, right in bonds:
            linked.add(left)
            linked.add(right)
        # Liquid stays. A crystal oxygen with no bond left is a dangling end.
        keep = {
            index for index in keep
            if int(kept_labels[index]) == 2 or index in linked
        }
        if len(keep) < len(xyz):
            order = sorted(keep)
            remap = {old: new for new, old in enumerate(order)}
            xyz = xyz[order]
            colors = [colors[index] for index in order]
            names = [names[index] for index in order]
            radii = [radii[index] for index in order]
            transparency = [transparency[index] for index in order]
            kept_labels = [kept_labels[index] for index in order]
            bonds = [(remap[left], remap[right]) for left, right in bonds]
    if slab:
        # One oxygen radius, the same one Ice Ih uses. Phase is colour alone.
        oxygen = {
            "hexagonal", "cubic", "water", "interfacial", "clathrate",
            "interClathrate", "unclassified", "reCubic", "reHex", "mixed",
        }
        for index, name in enumerate(names):
            if name in oxygen:
                radii[index] = ATOM_RADIUS
                transparency[index] = 0.0
    if view == "channel":
        radii = [value * (0.72 if name not in ("sodium", "chloride", "_rim") else 1.0) for value, name in zip(radii, names)]
    bond_radius = BOND_RADIUS * (1.45 if cluster else 1.0)
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


def sh_ring_mesh(oxygens, guest_centre, box, kind="20"):
    """One sH cage, with each face taken from its own ring.

    12 is the 5^12 cage, 12_1 is 4^3 5^6 6^3, and 20 is 5^12 6^8.
    Waters are the ideal cage snapped onto the relaxed oxygens. A convex
    hull of a wide shell pulls in neighbouring cages.
    """
    import numpy as np

    if kind == "12":
        legend, colour = "small cage", palette_rgb("hexagonal")
        homes_frac = ((0.5, 0.5, 0.5),)
        offsets = np.array([
            3.1059, -1.7932, 1.4072, 3.1059, -1.7932, -1.4072, -3.1059, 1.7932, -1.4072,
            -3.1059, 1.7932, 1.4072, -0.6950, -1.2037, -3.6518, 0.6950, 1.2037, -3.6518,
            -0.6950, -1.2037, 3.6538, 0.6950, 1.2037, 3.6538, -0.7819, 3.1335, -2.2607,
            -0.7819, 3.1335, 2.2607, 0.7819, -3.1335, 2.2607, 0.7819, -3.1335, -2.2607,
            -3.1059, -0.8896, -2.2607, 3.1059, 0.8896, -2.2607, 3.1059, 0.8896, 2.2607,
            -3.1059, -0.8896, 2.2607, -3.1047, -2.5720, 0.0000, 3.1047, 2.5720, 0.0000,
            0.6745, 3.9758, 0.0000, -0.6745, -3.9758, 0.0000,
        ], dtype=float).reshape(-1, 3)
        rings = (
            (0, 1, 11, 19, 10), (0, 1, 13, 17, 14), (0, 10, 6, 7, 14), (1, 11, 4, 5, 13),
            (2, 3, 9, 18, 8), (2, 3, 15, 16, 12), (2, 8, 5, 4, 12), (3, 9, 7, 6, 15),
            (4, 11, 19, 16, 12), (5, 8, 18, 17, 13), (6, 10, 19, 16, 15), (7, 9, 18, 17, 14),
        )
    elif kind == "12_1":
        legend, colour = "medium cage", (0.0, 0.537, 0.482)
        homes_frac = ((0.75, 0.3333, 0.75),)
        offsets = np.array([
            2.7650, 0.5895, -1.3790, 2.7650, 0.5895, 1.3810, 1.2875, -1.3414, 2.7721,
            1.2875, -1.3414, -2.7721, 1.3750, 2.9969, -1.3790, 1.3750, 2.9969, 1.3810,
            -1.0347, 0.0000, -3.6256, -1.0347, 0.0000, 3.6256, 0.3546, -3.5853, -1.3790,
            0.3546, -3.5853, 1.3810, -1.0347, 2.6828, 2.7721, -1.0347, 2.6828, -2.7721,
            -2.4252, -3.5853, -1.3790, -2.4252, -3.5853, 1.3810, -3.3587, -1.3403, 2.7721,
            -3.3587, -1.3403, -2.7721, -3.4456, 2.9969, -1.3790, -3.4456, 2.9969, 1.3810,
            -4.8356, 0.5895, -1.3790, -4.8356, 0.5895, 1.3810,
        ], dtype=float).reshape(-1, 3)
        rings = (
            (0, 1, 5, 4), (8, 9, 13, 12), (16, 17, 19, 18),
            (0, 3, 6, 11, 4), (1, 2, 7, 10, 5), (2, 7, 14, 13, 9),
            (3, 6, 15, 12, 8), (6, 11, 16, 18, 15), (7, 10, 17, 19, 14),
            (0, 1, 2, 9, 8, 3), (4, 5, 10, 17, 16, 11), (12, 13, 14, 19, 18, 15),
        )
    else:
        legend, colour = "large cage", palette_rgb("clathrate")
        homes_frac = ((1.0, 0.0, 0.2), (0.0, 0.0, 0.75))
        offsets = np.array([
        3.8878, -2.2439, -1.7655, -3.8878, 2.2439, -1.7655, 3.8872, 2.2450, -1.7655,
        -0.0006, 4.4889, -1.7655, 0.0006, -4.4889, -1.7655, -3.8872, -2.2450, -1.7655,
        -2.4103, 4.1748, -0.3724, 2.4103, 4.1748, -0.3724, -2.4103, -4.1748, -0.3724,
        2.4103, -4.1748, -0.3724, 4.8207, 0.0000, -0.3724, -4.8207, 0.0000, -0.3724,
        -2.4302, 1.4038, -4.0263, 2.4302, -1.4038, -4.0263, -0.0006, -2.8065, -4.0263,
        -2.4308, -1.4027, -4.0263, 2.4308, 1.4027, -4.0263, 0.0006, 2.8065, -4.0263,
        2.4103, 4.1748, 2.3876, -2.4103, 4.1748, 2.3876, -2.4103, -4.1748, 2.3876,
        2.4103, -4.1748, 2.3876, 4.8207, 0.0000, 2.3876, -4.8207, 0.0000, 2.3876,
        3.8872, 2.2450, 3.7786, -3.8878, 2.2439, 3.7786, 3.8878, -2.2439, 3.7786,
        -0.0006, 4.4889, 3.7786, 0.0006, -4.4889, 3.7786, -3.8872, -2.2450, 3.7786,
        -2.4302, 1.4038, 6.0394, 0.0006, 2.8065, 6.0394, -2.4308, -1.4027, 6.0394,
        -0.0006, -2.8065, 6.0394, 2.4302, -1.4038, 6.0394,         2.4308, 1.4027, 6.0394,
        ], dtype=float).reshape(-1, 3)
        rings = (
            (0, 9, 4, 14, 13), (0, 10, 2, 16, 13), (1, 6, 3, 17, 12), (1, 11, 5, 15, 12),
            (2, 7, 3, 17, 16), (4, 8, 5, 15, 14), (18, 24, 35, 31, 27), (19, 25, 30, 31, 27),
            (20, 28, 33, 32, 29), (21, 26, 34, 33, 28), (22, 24, 35, 34, 26), (23, 25, 30, 32, 29),
            (0, 9, 21, 26, 22, 10), (1, 6, 19, 25, 23, 11), (2, 7, 18, 24, 22, 10), (3, 6, 19, 27, 18, 7),
            (4, 8, 20, 28, 21, 9), (5, 8, 20, 29, 23, 11), (12, 15, 14, 13, 16, 17), (30, 31, 35, 34, 33, 32),
        )
    length_a = 12.4212113964
    length_c = 10.0656455142 * 2.0
    angle = np.deg2rad(120.0)

    def cartesian(frac):
        frac = np.atleast_2d(np.asarray(frac, dtype=float))
        return np.column_stack((
            length_a * frac[:, 0] + length_a * np.cos(angle) * frac[:, 1],
            length_a * np.sin(angle) * frac[:, 1],
            length_c * frac[:, 2],
        ))

    lattice = cartesian([[1.0, 0.0, 0.0], [0.0, 1.0, 0.0], [0.0, 0.0, 1.0]])
    homes = cartesian(homes_frac)
    guest = np.asarray(guest_centre, dtype=float)
    images = []
    for ix in (-1, 0, 1):
        for iy in (-1, 0, 1):
            for iz in (-1, 0, 1):
                images.append(homes + ix * lattice[0] + iy * lattice[1] + iz * lattice[2])
    images = np.vstack(images)
    home = images[int(np.linalg.norm(images - guest, axis=1).argmin())]
    ideal = home + offsets
    oxy = np.asarray(oxygens, dtype=float)
    copies = []
    for ix in (-1, 0, 1):
        for iy in (-1, 0, 1):
            for iz in (-1, 0, 1):
                copies.append(oxy + ix * lattice[0] + iy * lattice[1] + iz * lattice[2])
    copies = np.vstack(copies)
    shell = []
    for point in ideal:
        shell.append(copies[int(np.linalg.norm(copies - point, axis=1).argmin())])
    shell = np.asarray(shell, dtype=float)
    vertices = [row for row in shell]
    faces = []
    for ring in rings:
        points = shell[list(ring)]
        centre = points.mean(axis=0)
        normal = np.zeros(3)
        for index, left in enumerate(points):
            right = points[(index + 1) % len(points)]
            normal[0] += (left[1] - right[1]) * (left[2] + right[2])
            normal[1] += (left[2] - right[2]) * (left[0] + right[0])
            normal[2] += (left[0] - right[0]) * (left[1] + right[1])
        length = float(np.linalg.norm(normal)) or 1.0
        normal = normal / length
        if float(np.dot(normal, home - centre)) < 0.0:
            normal = -normal
        centre = centre + normal * 0.12
        base = len(vertices)
        vertices.append(centre)
        for index in range(len(ring)):
            faces.append((base, ring[index], ring[(index + 1) % len(ring)]))
    return {
        "vertices": np.asarray(vertices, dtype=float),
        "shell": shell,
        "faces": np.asarray(faces, dtype=int),
        "color": colour,
        "transparency": 0.55,
        "legend": legend,
        "alpha": 0.42,
    }


def host_token(lattice):
    if lattice in ("1h", "ice1h"):
        return "hexagonal"
    if lattice in ("1c", "ice1c"):
        return "cubic"
    return "clathrate"


def _pentagon_matrix(shell, faces):
    """Rotate a 5^12 cage so one pentagon faces the camera, vertex up."""
    import numpy as np

    shell = np.asarray(shell, dtype=float)
    centre = shell.mean(axis=0)
    groups = {}
    for tri in np.asarray(faces, dtype=int):
        groups.setdefault(int(tri[0]), []).append((int(tri[1]), int(tri[2])))
    best_normal = None
    best_cycle = None
    best_score = -2.0
    for _hub, edges in groups.items():
        nxt = {left: right for left, right in edges}
        if len(nxt) < 5:
            continue
        start = next(iter(nxt))
        cycle = [start]
        seen = {start}
        while nxt.get(cycle[-1]) not in (None, start) and len(cycle) < 8:
            nxt_i = nxt[cycle[-1]]
            if nxt_i in seen:
                break
            cycle.append(nxt_i)
            seen.add(nxt_i)
        if nxt.get(cycle[-1]) != start or len(cycle) != 5:
            continue
        pts = shell[cycle]
        normal = np.zeros(3)
        for index, left in enumerate(pts):
            right = pts[(index + 1) % len(pts)]
            normal[0] += (left[1] - right[1]) * (left[2] + right[2])
            normal[1] += (left[2] - right[2]) * (left[0] + right[0])
            normal[2] += (left[0] - right[0]) * (left[1] + right[1])
        length = float(np.linalg.norm(normal)) or 1.0
        normal = normal / length
        if float(np.dot(normal, pts.mean(axis=0) - centre)) < 0.0:
            normal = -normal
        score = float(normal[2])
        if score > best_score:
            best_score = score
            best_normal = normal
            best_cycle = cycle
    if best_normal is None:
        return np.eye(3)
    aiming = _rotation_taking(best_normal, np.array([0.0, 0.0, 1.0]))
    turned = (shell[best_cycle] - centre) @ aiming.T
    top = turned[int(np.argmax(turned[:, 1]))]
    yaw = float(np.arctan2(top[0], top[1]))
    cosine, sine = np.cos(yaw), np.sin(yaw)
    upright = np.array([
        [cosine, -sine, 0.0],
        [sine, cosine, 0.0],
        [0.0, 0.0, 1.0],
    ])
    return upright @ aiming


def _pose_cage(legend, waters, vertices, guests, faces):
    """Small cage shows a pentagon. The others take a small tilt."""
    import numpy as np

    origin = np.asarray(waters, dtype=float).mean(axis=0)
    if legend == "small cage":
        matrix = _pentagon_matrix(waters, faces)
    else:
        angle = np.deg2rad(18.0)
        cosine, sine = np.cos(angle), np.sin(angle)
        matrix = np.array([
            [1.0, 0.0, 0.0],
            [0.0, cosine, -sine],
            [0.0, sine, cosine],
        ])

    def turn(points):
        return (np.asarray(points, dtype=float) - origin) @ matrix.T + origin

    return turn(waters), turn(vertices), [turn(point) for point in guests]


def side_by_side_cages(meshes, xyz, names, colors, radii, transparency, residue_index):
    """One water radius, two cages, a gap between them.

    Each cage is its own shell. Shared crystal waters are copied so the
    cages do not occupy the same space.
    """
    import numpy as np

    guest_names = {"methane", "THF", "ring O", "hydrogen"}
    water_colour = palette_rgb("water")
    xyz = np.asarray(xyz, dtype=float)
    centres = [
        np.asarray(mesh["shell"] if "shell" in mesh else mesh["vertices"], dtype=float).mean(axis=0)
        for mesh in meshes
    ]
    guests = {}
    for index, name in enumerate(names):
        if name not in guest_names:
            continue
        guests.setdefault(residue_index[index], []).append(index)
    order = sorted(range(len(meshes)), key=lambda slot: 0 if meshes[slot].get("legend") == "small cage" else 1)
    cursor = None
    # Three sH cages sit closer so the large one can match the sII large cage.
    gap = 1.4 if any(mesh.get("legend") == "medium cage" for mesh in meshes) else 4.0
    out_pos, out_col, out_name, out_rad, out_alpha = [], [], [], [], []
    bonds = []
    shifted = []
    for slot in order:
        mesh = meshes[slot]
        waters = np.asarray(mesh["shell"] if "shell" in mesh else mesh["vertices"], dtype=float).copy()
        mine = []
        best_dist = 1.0e9
        for indices in guests.values():
            centre = xyz[indices].mean(axis=0)
            dist = float(np.linalg.norm(centre - centres[slot]))
            guest_kind = {names[index] for index in indices}
            limit = 5.5 if guest_kind & {"THF", "ring O"} else 3.2
            if dist < best_dist and dist <= limit:
                best_dist = dist
                mine = list(indices)
        guest_points = [np.array(xyz[index], dtype=float, copy=True) for index in mine]
        if "shell" in mesh:
            waters, posed_vertices, guest_points = _pose_cage(
                mesh.get("legend"), waters, mesh["vertices"], guest_points, mesh["faces"],
            )
            mesh = dict(mesh)
            mesh["vertices"] = posed_vertices
            mesh["shell"] = np.asarray(waters, dtype=float)
        stacked = waters
        if guest_points:
            stacked = np.vstack([waters, np.asarray(guest_points, dtype=float)])
        xmin = float(stacked[:, 0].min())
        xmax = float(stacked[:, 0].max())
        shift = 0.0 if cursor is None else cursor + gap - xmin
        cursor = xmax + shift
        # Share one vertical centre so a small cage is not pushed off the top.
        anchor = waters.mean(axis=0)
        waters[:, 0] += shift
        waters[:, 1] -= anchor[1]
        waters[:, 2] -= anchor[2]
        moved = {
            "vertices": np.asarray(mesh["vertices"], dtype=float).copy(),
            "faces": mesh["faces"],
            "color": mesh["color"],
            "transparency": mesh.get("transparency", 0.55),
            "legend": mesh.get("legend"),
            "alpha": 0.42,
        }
        if "shell" in mesh:
            moved["shell"] = np.asarray(mesh["shell"], dtype=float).copy()
            moved["shell"][:, 0] += shift
            moved["shell"][:, 1] -= anchor[1]
            moved["shell"][:, 2] -= anchor[2]
        moved["vertices"][:, 0] += shift
        moved["vertices"][:, 1] -= anchor[1]
        moved["vertices"][:, 2] -= anchor[2]
        shifted.append(moved)
        start = len(out_pos)
        for point in waters:
            out_pos.append(point)
            out_col.append(water_colour)
            out_name.append("water")
            out_rad.append(ATOM_RADIUS)
            out_alpha.append(0.0)
        for left in range(len(waters)):
            for right in range(left + 1, len(waters)):
                span = float(np.linalg.norm(waters[left] - waters[right]))
                limit = 3.30 if "shell" in mesh else 3.20
                if 2.40 <= span <= limit:
                    bonds.append((start + left, start + right))
        guest_at = start + len(waters)
        for offset, index in enumerate(mine):
            point = np.array(guest_points[offset], dtype=float, copy=True)
            point[0] += shift
            point[1] -= anchor[1]
            point[2] -= anchor[2]
            out_pos.append(point)
            out_col.append(colors[index])
            out_name.append(names[index])
            out_rad.append(radii[index] if names[index] == "hydrogen" else ATOM_RADIUS)
            out_alpha.append(transparency[index])
            for earlier in range(offset):
                other = mine[earlier]
                if residue_index[other] != residue_index[index]:
                    continue
                if float(np.linalg.norm(out_pos[guest_at + offset] - out_pos[guest_at + earlier])) < 1.70:
                    bonds.append((guest_at + earlier, guest_at + offset))
    return (
        np.asarray(out_pos, dtype=float),
        out_col,
        out_rad,
        out_alpha,
        out_name,
        bonds,
        shifted,
    )


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
            elif kind == "THF" and str(atom).startswith("O"):
                group["heavy"].append(xyz)
                colour, radius, alpha = atom_style("ring O")
                names.append("ring O")
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
        meshes, _extra, _highlighted = cage_meshes(
            np.asarray(oxygens, dtype=float),
            [
                (group["kind"], np.asarray(group["heavy"], dtype=float))
                for group in guest_groups.values()
                if group["heavy"]
            ],
            box,
        )
        if lattice == "sH":
            built = []
            for kind in ("12", "12_1", "20"):
                best_mesh = None
                best_dist = 1.0e9
                want = "THF" if kind == "20" else "CH4"
                for group in guest_groups.values():
                    if group["kind"] != want or not group["heavy"]:
                        continue
                    centre = np.mean(group["heavy"], axis=0)
                    trial = sh_ring_mesh(np.asarray(oxygens, dtype=float), centre, box, kind=kind)
                    dist = float(np.linalg.norm(centre - trial["shell"].mean(axis=0)))
                    if dist < best_dist:
                        best_dist = dist
                        best_mesh = trial
                if best_mesh is not None and best_dist < 3.0:
                    built.append(best_mesh)
            if built:
                meshes = built
        if meshes:
            xyz, colors, radii, transparency, names, bonds, meshes = side_by_side_cages(
                meshes, xyz, names, colors, radii, transparency, residue_index,
            )
            # sH is posed per cage, so the row view must not tilt the pentagon.
            view = "row" if lattice == "sH" else "pair"
            return finish(
                xyz, colors, radii, transparency, names, bonds, meshes,
                view, xyz.mean(axis=0), False, bold=bonds,
                bond_radius=BOND_RADIUS, bond_rgb=(0.33, 0.35, 0.38),
            )
    return finish(
        xyz, colors, radii, transparency, names, bonds, meshes,
        "cage", xyz.mean(axis=0), False,
        bond_radius=BOND_RADIUS, bond_rgb=(0.33, 0.35, 0.38),
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


def front_sheet(frame, labels, axis=1, keep_liquid=False):
    """One ice bilayer. Stacked layers are what drew each bond two or three times."""
    import numpy as np

    y = np.asarray(frame.xyz[:, int(axis)], dtype=float)
    crystal = [i for i, label in enumerate(labels) if int(label) != 2]
    if len(crystal) < 8:
        return frame, labels
    centres = []
    for value in np.sort(y[crystal]):
        if not centres or float(value) - centres[-1] > 0.55:
            centres.append(float(value))
    best = None
    best_n = -1
    for left, right in zip(centres, centres[1:]):
        if right - left > 1.6:
            continue
        members = [
            i for i in crystal
            if abs(float(y[i]) - left) <= 0.45 or abs(float(y[i]) - right) <= 0.45
        ]
        if len(members) > best_n:
            best = members
            best_n = len(members)
    if not best:
        return frame, labels
    chosen = set(best)
    if keep_liquid:
        chosen.update(i for i, label in enumerate(labels) if int(label) == 2)
    index = sorted(chosen)
    return frame.subset(index), [labels[i] for i in index]


def closed_crystal_rings(frame, labels):
    """Keep crystal oxygens that sit on a closed six-ring. Liquid stays."""
    crystal = [i for i, label in enumerate(labels) if int(label) != 2]
    if len(crystal) < 6:
        return frame, labels
    local = whole_hexagon_index(frame.xyz[crystal])
    chosen = {crystal[int(i)] for i in local}
    chosen.update(i for i, label in enumerate(labels) if int(label) == 2)
    index = sorted(chosen)
    return frame.subset(index), [labels[i] for i in index]


def quiet_stray_clathrate(frame, labels):
    """CHILL+ can mark a liquid oxygen interClathrate when every neighbour is water.

    Those two atoms in the interface cell have no clathrate neighbour, so the
    still draws them as liquid instead of as dark green specks.
    """
    import numpy as np
    from scipy.spatial import cKDTree

    labels = list(labels)
    xyz = np.asarray(frame.xyz, dtype=float)
    marked = [i for i, label in enumerate(labels) if int(label) in (4, 5)]
    if not marked:
        return labels
    tree = cKDTree(xyz)
    for index in marked:
        _dist, neigh = tree.query(xyz[index], k=8)
        friends = 0
        for other in np.atleast_1d(neigh):
            other = int(other)
            if other == index:
                continue
            if float(np.linalg.norm(xyz[other] - xyz[index])) > 3.5:
                continue
            if int(labels[other]) in (4, 5):
                friends += 1
        if friends < 2:
            labels[index] = 2
    return labels


def depth_slab(frame, labels, x_lo=6.0, x_hi=12.5, z_lo=4.5, z_hi=26.0, liquid_reach=6.0):
    """One depth through the interface: ice, then the interfacial layer, then liquid.

    x is the depth. The window is a few angstroms, so a liquid ion does not
    land on the ice by projection. Liquid is the shell that touches the ice,
    not the far solvent, which is what opened a gap in the picture.
    """
    import numpy as np
    from scipy.spatial import cKDTree

    labels = list(labels)
    for index, label in enumerate(labels):
        if int(label) in (4, 5):
            labels[index] = 2
    xyz = np.asarray(frame.xyz, dtype=float)
    depth = (xyz[:, 0] >= x_lo) & (xyz[:, 0] < x_hi)
    hexagonal = depth & np.array([int(label) == 1 for label in labels])
    if not np.any(hexagonal):
        index = np.where(depth)[0]
        return frame.subset(index), [labels[int(i)] for i in index]
    front = float(np.percentile(xyz[hexagonal, 2], 75))
    ice_max = float(xyz[hexagonal, 2].max())
    # Ice through the front, plus the interfacial waters just beyond it.
    # The periodic interface on the low-z side is left out.
    chosen = set(
        int(i) for i in np.where(
            depth & (xyz[:, 2] >= z_lo) & (xyz[:, 2] <= ice_max + 0.6)
            & np.array([int(label) != 2 for label in labels])
        )[0]
    )
    chosen.update(
        int(i) for i in np.where(
            depth & np.array([int(label) == 3 for label in labels])
            & (xyz[:, 2] >= front) & (xyz[:, 2] <= ice_max + 6.0)
        )[0]
    )
    liquid = [
        i for i, label in enumerate(labels)
        if depth[i] and int(label) == 2
        and ice_max - 1.5 <= float(xyz[i, 2]) <= ice_max + 14.0
    ]
    chosen.update(liquid)
    # Ions that sit on a drawn water. Type 3 in the ice is a chloride.
    # A liquid sodium (type 4) is added below when the slice would otherwise miss it.
    if chosen and getattr(frame, "types", None) is not None:
        tree = cKDTree(xyz[sorted(chosen)])
        for index, kind in enumerate(frame.types):
            if int(kind) not in (3, 4) or not depth[index]:
                continue
            gap, _near = tree.query(xyz[index], k=1)
            if float(gap) <= 4.0 and z_lo <= float(xyz[index, 2]) <= ice_max + 16.0:
                chosen.add(int(index))
        liquid_sodium = [
            index for index in chosen
            if int(frame.types[index]) == 4 and float(xyz[index, 2]) > ice_max
        ]
        if not liquid_sodium:
            outside = [
                index for index, kind in enumerate(frame.types)
                if int(kind) == 4 and depth[index] and float(xyz[index, 2]) > ice_max
            ]
            if outside:
                pick = min(outside, key=lambda index: float(xyz[index, 2]))
                chosen.add(int(pick))
                for other, label in enumerate(labels):
                    if not depth[other] or int(label) != 2:
                        continue
                    if float(np.linalg.norm(xyz[other] - xyz[pick])) <= 4.0:
                        chosen.add(int(other))
        # The liquid sodium must not sit on a chloride. Keep the sodium
        # whose nearest chloride is furthest away, and drop any chloride
        # still within 4 angstroms of it.
        chlorides = [i for i in chosen if int(frame.types[i]) == 3]
        liquid_sodium = [
            i for i in chosen
            if int(frame.types[i]) == 4 and float(xyz[i, 2]) > ice_max
        ]
        if liquid_sodium and chlorides:
            def clearance(index):
                return min(float(np.linalg.norm(xyz[index] - xyz[other])) for other in chlorides)

            separated = [index for index in liquid_sodium if clearance(index) >= 4.0]
            pick = separated[0] if separated else max(liquid_sodium, key=clearance)
            for other in liquid_sodium:
                if other != pick:
                    chosen.discard(other)
            if clearance(pick) < 4.0:
                for other in chlorides:
                    if float(np.linalg.norm(xyz[pick] - xyz[other])) < 4.0:
                        chosen.discard(other)
    index = sorted(chosen)
    return frame.subset(index), [labels[i] for i in index]


def best_nucleus(library):
    """The largest cubic/hexagonal cluster in the trajectory.

    Frame 1 is the historical still. A later frame replaces it when that
    frame holds a larger crystal.
    """
    path = ROOT / "input/traj/nucleation_three.lammpstrj"
    best = None
    for index in range(3):
        frame = read_dump(path, index)
        labels = chill_labels(library, frame)
        chosen, chosen_labels, centre = largest_cluster(
            frame, labels, shell=3.0, core=(0, 1, 7, 8),
        )
        score = sum(1 for label in chosen_labels if int(label) in (0, 1, 7, 8))
        if best is None or score > best[0] or (score == best[0] and index > best[1]):
            best = (score, index, chosen, chosen_labels, centre)
    print("nucleation frame {} crystal {}".format(best[1], best[0]), flush=True)
    return best[2], best[3], best[4]


def draw(library, scene):
    global PANEL_NAME
    PANEL_NAME = scene
    if scene == "ice-ih":
        frame = read_dump(ROOT / "input/traj/ice1h.lammpstrj")
        labels = chill_labels(library, frame)
        frame, labels = tile_frame(frame, labels, 3, 1, 2)
        frame, labels = front_sheet(frame, labels, axis=1, keep_liquid=False)
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
        frame, labels = depth_slab(frame, labels)
        return labelled_picture(frame, labels, "slab", False, slab=True)
    if scene == "nucleation":
        frame, labels, centre = best_nucleus(library)
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
        # Type 3 at about (11, 10, 12) is the chloride in the crystalline
        # slab. CHILL+ runs on the water alone, so the ion is not ice.
        frame, labels = depth_slab(frame, labels)
        return labelled_picture(frame, labels, "slab", False, ions=True, slab=True)
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
        "brine": "NaCl front. The chloride in the ice slice is not an ice label. CHILL+ runs on the water alone.",
        "interface": "Basal interface, one ice bilayer. Crystal bonds stop at the ice. An interClathrate oxygen with only liquid neighbours is drawn as water.",
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
