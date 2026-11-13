#!/usr/bin/env python3
# Copyright (c) d-SEAMS developers
"""Draw molecule stills with chemparseplot.

Run from a checkout, after libseams is built:

    python docs/scripts/render_molecules.py --out <dir> [--slug <slug>]

Each still is written to <dir>/<slug>/ as a 2400 by 1800 PNG and a
smaller web PNG. The pictures are build products and are not stored
in git.
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
WIDTH = 2400
HEIGHT = 1800
WEB = (960, 720)

# Okabe-Ito. One radius for every water oxygen, so a label is a colour and not a size.
LABEL_COLOR = {
    0: (0.00, 0.45, 0.70),
    1: (0.90, 0.62, 0.00),
    2: (0.35, 0.70, 0.90),
    3: (0.95, 0.90, 0.25),
    4: (0.80, 0.47, 0.65),
    5: (0.00, 0.62, 0.45),
    6: (0.55, 0.55, 0.55),
    7: (0.00, 0.45, 0.70),
    8: (0.90, 0.62, 0.00),
}
ION_COLOR = {"sodium": (0.55, 0.05, 0.55), "chloride": (0.10, 0.55, 0.20)}
GUEST_COLOR = (0.98, 0.78, 0.15)
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
        if name in seen:
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


def orient(positions, view, center):
    import numpy as np

    pos = np.asarray(positions, dtype=float) - np.asarray(center, dtype=float)
    if view == "front":
        pos = np.column_stack((pos[:, 2], pos[:, 1], -pos[:, 0]))
    tilt = np.deg2rad(18.0 if view != "cage" else 32.0)
    cosine, sine = np.cos(tilt), np.sin(tilt)
    y_coord, z_coord = pos[:, 1].copy(), pos[:, 2].copy()
    pos[:, 1] = cosine * y_coord - sine * z_coord
    pos[:, 2] = sine * y_coord + cosine * z_coord
    if view == "cage":
        turn = np.deg2rad(24.0)
        cosine, sine = np.cos(turn), np.sin(turn)
        x_coord, z_coord = pos[:, 0].copy(), pos[:, 2].copy()
        pos[:, 0] = cosine * x_coord + sine * z_coord
        pos[:, 2] = -sine * x_coord + cosine * z_coord
    return pos


def tachyon(positions, colors, radii, bonds):
    """Ray-trace with ambient occlusion. Bonds are the oxygen contacts."""
    import tempfile

    import numpy as np
    import ovito
    from ovito.data import DataCollection, Particles
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
    data.objects.append(particles)
    if len(bonds):
        topology = particles.create_bonds().create_property("Topology", data=np.asarray(bonds, dtype=int))
        del topology
        particles.bonds.vis.width = 0.12
        particles.bonds.vis.color = (0.22, 0.22, 0.22)
    particles.vis.radius = 0.45
    pipeline = Pipeline(source=StaticSource(data=data))
    pipeline.add_to_scene()
    handle = tempfile.NamedTemporaryFile(suffix=".png", delete=False)
    handle.close()
    try:
        view = Viewport(type=Viewport.Type.Ortho)
        view.zoom_all()
        renderer = TachyonRenderer(
            ambient_occlusion=True,
            ambient_occlusion_brightness=0.82,
            ambient_occlusion_samples=12,
            antialiasing=True,
            antialiasing_samples=8,
            shadows=True,
        )
        view.render_image(
            filename=handle.name,
            size=(1880, 1620),
            background=(1, 1, 1),
            renderer=renderer,
            frame=0,
        )
        picture = Image.open(handle.name).convert("RGB")
        world = 2.0 * float(view.fov) if view.fov else 0.0
        return picture, world
    finally:
        pipeline.remove_from_scene()
        try:
            os.unlink(handle.name)
        except OSError:
            pass


def crop_margin(picture):
    import numpy as np
    from PIL import Image

    array = np.asarray(picture)
    mask = np.any(array < 248, axis=2)
    rows = np.where(mask.any(axis=1))[0]
    cols = np.where(mask.any(axis=0))[0]
    if len(rows) == 0:
        return picture
    pad = 18
    top = max(0, int(rows[0]) - pad)
    bottom = min(array.shape[0], int(rows[-1]) + pad)
    left = max(0, int(cols[0]) - pad)
    right = min(array.shape[1], int(cols[-1]) + pad)
    return Image.fromarray(array[top:bottom, left:right])


def compose(structure, legend, world_height):
    from PIL import Image, ImageDraw, ImageFont

    canvas = Image.new("RGB", (WIDTH, HEIGHT), (255, 255, 255))
    picture = crop_margin(structure)
    picture.thumbnail((1760, 1580), Image.Resampling.LANCZOS)
    left, top = 48, 56
    canvas.paste(picture, (left, top))
    draw = ImageDraw.Draw(canvas)
    font = ImageFont.load_default()
    for path in (
        "/usr/share/fonts/truetype/dejavu/DejaVuSans.ttf",
        "/usr/share/fonts/truetype/liberation/LiberationSans-Regular.ttf",
    ):
        if os.path.isfile(path):
            font = ImageFont.truetype(path, 32)
            break
    x = left + picture.size[0] + 28
    y = top
    for name, rgb in legend:
        colour = tuple(int(max(0, min(255, channel * 255))) for channel in rgb)
        draw.rectangle((x, y, x + 40, y + 40), fill=colour, outline=(20, 20, 20), width=2)
        draw.text((x + 54, y + 4), name, fill=(20, 20, 20), font=font)
        y += 58
    if world_height > 0:
        bar = 10.0 * (picture.size[1] / world_height)
        bar = max(80.0, min(bar, picture.size[0] * 0.45))
        bar_y = top + picture.size[1] + 28
        draw.line((left, bar_y, left + bar, bar_y), fill=(20, 20, 20), width=6)
        draw.line((left, bar_y - 12, left, bar_y + 12), fill=(20, 20, 20), width=4)
        draw.line((left + bar, bar_y - 12, left + bar, bar_y + 12), fill=(20, 20, 20), width=4)
        draw.text((left, bar_y + 16), "10 A", fill=(20, 20, 20), font=font)
    return canvas


def save_pair(image, directory, stem):
    directory.mkdir(parents=True, exist_ok=True)
    dest = directory / "{}.png".format(stem)
    web = directory / "{}-web.png".format(stem)
    image.save(dest)
    from PIL import Image

    small = image.resize(WEB, Image.Resampling.LANCZOS)
    small.save(web)
    return dest


def cluster_shell(frame, labels, shell=4.5):
    import numpy as np

    ice = [i for i, label in enumerate(labels) if label in (0, 1, 7, 8)]
    if not ice:
        return central_slice(frame, labels, 700)
    centre = frame.xyz[ice]
    keep = []
    for index, point in enumerate(frame.xyz):
        gap = np.linalg.norm(centre - point, axis=1).min()
        if gap <= shell:
            keep.append(index)
    return frame.subset(keep), [labels[i] for i in keep]


def paint(frame, labels, view, guests=False, ions=False, highlight=None):
    import numpy as np

    count = len(frame.xyz)
    colors = []
    names = []
    radii = []
    heavy = []
    for index in range(count):
        kind = int(frame.types[index])
        if ions and kind == 3:
            colors.append(ION_COLOR["sodium"])
            names.append("sodium")
            radii.append(0.78)
            heavy.append(False)
        elif ions and kind == 4:
            colors.append(ION_COLOR["chloride"])
            names.append("chloride")
            radii.append(0.95)
            heavy.append(False)
        elif labels[index] < 0:
            colors.append(GUEST_COLOR)
            names.append("guest")
            radii.append(0.85)
            heavy.append(False)
        else:
            label = int(labels[index])
            if highlight is not None and index in highlight:
                colors.append((0.85, 0.20, 0.15))
                names.append("highlighted")
            else:
                colors.append(LABEL_COLOR.get(label, (0.5, 0.5, 0.5)))
                names.append(CHILL_NAME.get(label, "other"))
            radii.append(0.45)
            heavy.append(True)
    center = frame.xyz.mean(axis=0)
    positions = orient(frame.xyz, view, center)
    if guests:
        extra = place_guests(frame.xyz, frame.box9)
        if len(extra):
            positions = np.vstack((positions, orient(extra, view, center)))
            colors.extend([GUEST_COLOR] * len(extra))
            names.extend(["guest"] * len(extra))
            radii.extend([0.85] * len(extra))
            heavy.extend([False] * len(extra))
    host = np.asarray(positions[:count], dtype=float)
    bonds = oxygen_bonds(host, heavy[:count]) if any(heavy[:count]) else []
    picture, world = tachyon(positions, colors, radii, bonds)
    return compose(picture, legend_rows(names, colors), world)


def draw(library, scene):
    if scene == "ice-ih":
        frame = read_dump(ROOT / "input/traj/ice1h.lammpstrj")
        labels = chill_labels(library, frame)
        return paint(frame, labels, "basal")
    if scene == "ice-ic":
        frame = read_dump(ROOT / "input/traj/mW_cubic.lammpstrj")
        labels = chill_labels(library, frame)
        frame, labels = central_slice(frame, labels, 700)
        return paint(frame, labels, "basal")
    if scene == "interface":
        frame = read_dump(ROOT / "input/traj/interface_f0.lammpstrj")
        labels = chill_labels(library, frame)
        return paint(frame, labels, "front")
    if scene == "nucleation":
        frame = read_dump(ROOT / "input/traj/nucleation_three.lammpstrj", frame_index=1)
        labels = chill_labels(library, frame)
        frame, labels = cluster_shell(frame, labels)
        return paint(frame, labels, "cage")
    if scene == "brine":
        frame = read_dump(ROOT / "input/traj/nacl_front.lammpstrj")
        oxygen = [i for i, kind in enumerate(frame.types) if int(kind) == 1]
        sub = frame.subset(oxygen)
        labelled = chill_labels(library, sub)
        labels = [-1] * len(frame.xyz)
        for index, label in zip(oxygen, labelled):
            labels[index] = label
        return paint(frame, labels, "front", ions=True)
    if scene == "rings":
        frame = read_dump(ROOT / "input/traj/genice_sI.lammpstrj")
        labels = [2] * len(frame.xyz)
        return paint(frame, labels, "cage", highlight=ring_members(library, frame, {5, 6}))
    if scene == "prisms":
        frame = read_dump(ROOT / "input/traj/ice1h.lammpstrj")
        labels = chill_labels(library, frame)
        return paint(frame, labels, "basal", highlight=ring_members(library, frame, {6}))
    if scene in ("cage-si", "cage-sii", "cage-sh", "genice"):
        path = {
            "cage-si": "genice_sI.lammpstrj",
            "cage-sii": "genice_sII.lammpstrj",
            "cage-sh": "genice_sH.lammpstrj",
            "genice": "genice_fau.lammpstrj",
        }[scene]
        frame = read_dump(ROOT / "input/traj" / path)
        labels = [2] * len(frame.xyz)
        return paint(frame, labels, "cage", guests=(scene == "cage-sii"))
    if scene == "cage-guest":
        frame = guest_cage()
        labels = []
        for symbol in frame.symbols:
            if symbol == "H":
                labels.append(6)
            elif symbol == "O":
                labels.append(2)
            else:
                labels.append(-1)
        return paint(frame, labels, "cage")
    if scene == "solvis":
        frame = read_dump(ROOT / "input/traj/genice_sI.lammpstrj")
        return paint(frame, [2] * len(frame.xyz), "cage")
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
        "ice-ih", "ice-ic", "rings", "prisms", "cage-si", "cage-sii", "cage-sh",
        "cage-guest", "brine", "interface", "nucleation", "genice", "solvis",
    )
    captions = {
        "ice-ih": "Ice Ih coloured by published CHILL+.",
        "ice-ic": "Ice Ic coloured by published CHILL+.",
        "rings": "Five- and six-rings on the sI host.",
        "prisms": "Six-rings on ice Ih.",
        "cage-si": "GenIce sI host.",
        "cage-sii": "GenIce sII host with a guest in each cage centre.",
        "cage-sh": "GenIce sH host.",
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
    parser.add_argument("--out", required=True, help="Directory that receives <slug>/ stills")
    parser.add_argument("--slug", help="One video slug or one scene name")
    parser.add_argument("--libdir", help="Directory that holds libseams.so")
    parser.add_argument("--pages", action="store_true", help="Write the Sphinx figure list")
    parser.add_argument("--force", action="store_true", help="Draw again even if a PNG is already large enough")
    args = parser.parse_args(argv)
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
            "ice-ih", "ice-ic", "rings", "prisms", "cage-si", "cage-sii", "cage-sh",
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
    raise SystemExit(main(sys.argv[1:]))
