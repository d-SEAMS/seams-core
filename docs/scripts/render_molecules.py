#!/usr/bin/env python3
# Copyright (c) d-SEAMS developers
"""Draw molecule stills.

chemparseplot is the entry point for the oxygen-oxygen bond list.
It does not ray-trace, it does not colour atoms by a classifier,
and it does not draw translucent faces. The pixels are OVITO's
Tachyon ray tracer: soft shadows, ambient occlusion, depth of field
off, thin sticks, and triangle meshes for cage faces and rings.
OSPRay is available in the same module, but on this build its Vulkan
device does not initialise and the CPU path does not finish a
translucent cage at print resolution.

    python docs/scripts/render_molecules.py --out <dir> [--slug <slug>]

Each still is a 2400 by 1800 PNG, a web PNG, and a dark pair.
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
WIDTH = 2400
HEIGHT = 1800
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
    "sodium": "#6A1B9A",
    "chloride": "#2E7D32",
    "methane": "#F9A825",
    "THF": "#EF6C00",
    "hydrogen": "#90A4AE",
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


def atom_style(token):
    """Colour, radius and transparency. Liquid water is the small grey."""
    colour = palette_rgb(token)
    if token == "water":
        return colour, 0.20, 0.72
    if token == "hydrogen":
        return colour, 0.14, 0.0
    if token == "sodium":
        return colour, 0.62, 0.0
    if token == "chloride":
        return colour, 0.78, 0.0
    if token in ("methane", "THF"):
        return colour, 0.42, 0.0
    return colour, 0.36, 0.0

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


def raytrace(positions, colors, radii, transparency, bonds, meshes, background, perspective):
    """Tachyon still. Depth of field stays off. The render is twice the panel."""
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
        view.zoom_all()
        camera = np.array(view.camera_pos, dtype=float)
        forward = np.array(view.camera_dir, dtype=float)
        if float(np.dot(forward, -camera)) < 0.0:
            camera = -camera
            view.camera_pos = (float(camera[0]), float(camera[1]), float(camera[2]))
            view.camera_dir = (0.0, 0.0, -1.0)
        if perspective:
            view.camera_pos = (float(camera[0]) * 1.12, float(camera[1]) * 1.12, float(camera[2]) * 1.12)
            distance = abs(float(view.camera_pos[2]))
            world = 2.0 * distance * float(np.tan(np.deg2rad(float(view.fov) / 2.0)))
        else:
            view.fov = float(view.fov) * 1.12
            # Orthographic fov is the half-height of the view.
            world = 2.0 * float(view.fov)
        renderer = TachyonRenderer(
            ambient_occlusion=True,
            ambient_occlusion_brightness=0.55,
            ambient_occlusion_samples=4,
            antialiasing=True,
            antialiasing_samples=4,
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


def compose(structure, legend, world_height, light):
    """One family of frames: molecule, legend, and a 10 angstrom bar."""
    from PIL import Image, ImageDraw

    bg = (255, 255, 255) if light else (22, 24, 28)
    ink = (25, 25, 28) if light else (236, 237, 239)
    canvas = Image.new("RGB", (WIDTH, HEIGHT), bg)
    picture, crop_rows, render_rows = crop_margin(structure, light)
    legend_w = 360
    strip = 86
    pad = 18
    area_w = WIDTH - legend_w - 2 * pad
    area_h = HEIGHT - strip - 2 * pad
    scale = min(area_w / picture.size[0], area_h / picture.size[1])
    sized = (
        max(1, int(round(picture.size[0] * scale))),
        max(1, int(round(picture.size[1] * scale))),
    )
    picture = picture.resize(sized, Image.Resampling.LANCZOS)
    left = pad + (area_w - sized[0]) // 2
    top = pad + (area_h - sized[1]) // 2
    canvas.paste(picture, (left, top))
    draw = ImageDraw.Draw(canvas)
    font = _font(28)
    small = _font(26)
    x = WIDTH - legend_w + 8
    y = pad
    for name, colour_rgb in legend:
        colour = tuple(int(max(0, min(255, round(channel * 255)))) for channel in colour_rgb)
        draw.rounded_rectangle((x, y, x + 32, y + 32), radius=4, fill=colour, outline=ink, width=2)
        draw.text((x + 44, y + 1), name, fill=ink, font=font)
        y += 48
    if world_height > 0 and crop_rows > 0 and render_rows > 0:
        world_crop = world_height * (crop_rows / float(render_rows))
        px_per = sized[1] / world_crop
        bar = 10.0 * px_per
        bar_x = left
        bar_y = HEIGHT - 46
        draw.line((bar_x, bar_y, bar_x + bar, bar_y), fill=ink, width=5)
        draw.line((bar_x, bar_y - 9, bar_x, bar_y + 9), fill=ink, width=3)
        draw.line((bar_x + bar, bar_y - 9, bar_x + bar, bar_y + 9), fill=ink, width=3)
        draw.text((bar_x, bar_y + 12), "10 \u00c5", fill=ink, font=small)
    return canvas


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

    crystal = np.array([i for i, label in enumerate(labels) if label in (0, 1, 7, 8)], dtype=int)
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


def cage_meshes(oxygens, guests):
    """Convex hull of each guest whose water shell sits inside the cell."""
    import numpy as np
    from scipy.spatial import ConvexHull

    if len(guests) == 0 or len(oxygens) < 8:
        return []
    origin = np.asarray(oxygens, dtype=float).mean(axis=0)
    found = []
    for kind, heavy in guests:
        centre = np.asarray(heavy, dtype=float).mean(axis=0)
        distance = np.linalg.norm(oxygens - centre, axis=1)
        limit = 4.40 if kind == "CH4" else 5.40
        chosen = (distance > 2.4) & (distance < limit)
        if int(np.count_nonzero(chosen)) < 12:
            continue
        shell = oxygens[chosen]
        try:
            hull = ConvexHull(shell, qhull_options="QJ")
        except Exception:
            continue
        token = "hexagonal" if kind == "CH4" else "clathrate"
        found.append((float(np.linalg.norm(centre - origin)), kind, {
            "vertices": shell,
            "faces": hull.simplices.astype(int),
            "color": palette_rgb(token),
            "transparency": 0.86,
            "legend": "small cage" if kind == "CH4" else "large cage",
        }))
    meshes = []
    for kind in ("CH4", "THF"):
        ranked = sorted(item for item in found if item[1] == kind)
        meshes.extend(item[2] for item in ranked[:1])
    return meshes


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


def finish(positions, colors, radii, transparency, names, bonds, meshes, view, center, perspective):
    import numpy as np

    positions = orient(positions, view, center)
    depth = positions[:, 2]
    span = float(depth.max() - depth.min()) or 1.0
    radii = np.asarray(radii, dtype=float) * (0.90 + 0.12 * (depth - float(depth.min())) / span)
    turned = []
    for mesh in meshes:
        turned.append({
            "vertices": orient(mesh["vertices"], view, center),
            "faces": mesh["faces"],
            "color": mesh["color"],
            "transparency": mesh["transparency"],
        })
    legend = legend_rows(names, colors)
    seen = {name for name, _colour in legend}
    for mesh in meshes:
        label = mesh.get("legend")
        if label and label not in seen:
            legend.append((label, mesh["color"]))
            seen.add(label)
    light_picture, world = raytrace(
        positions, colors, radii, transparency, bonds, turned, (1, 1, 1), perspective
    )
    dark_picture, _world = raytrace(
        positions, colors, radii, transparency, bonds, turned, (0.09, 0.10, 0.12), perspective
    )
    return (
        compose(light_picture, legend, world, True),
        compose(dark_picture, legend, world, False),
    )


def labelled_picture(frame, labels, view, liquid_clear, ions=False, perspective=False, meshes=None, center=None):
    import numpy as np

    colors, names, radii, transparency, heavy = [], [], [], [], []
    for index in range(len(frame.xyz)):
        kind = int(frame.types[index])
        if ions and kind == 3:
            colour, radius, alpha = atom_style("sodium")
            colors.append(colour)
            names.append("sodium")
            radii.append(radius)
            transparency.append(alpha)
            heavy.append(False)
            continue
        if ions and kind == 4:
            colour, radius, alpha = atom_style("chloride")
            colors.append(colour)
            names.append("chloride")
            radii.append(radius)
            transparency.append(alpha)
            heavy.append(False)
            continue
        label = int(labels[index])
        token = CHILL_NAME.get(label, "unclassified")
        colour, radius, alpha = atom_style(token)
        colors.append(colour)
        names.append(token)
        radii.append(radius)
        transparency.append(alpha)
        heavy.append(True)
    xyz = np.asarray(frame.xyz, dtype=float)
    bonds = oxygen_bonds(xyz, heavy)
    if liquid_clear:
        bonds = [
            pair for pair in bonds
            if not (int(labels[pair[0]]) == 2 and int(labels[pair[1]]) == 2)
        ]
    return finish(
        xyz, colors, radii, transparency, names, bonds, meshes or [],
        view, frame.xyz.mean(axis=0) if center is None else center, perspective,
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
    """Move each molecule into the image nearest the cell centre."""
    import numpy as np
    from collections import defaultdict

    lengths = np.asarray(box[:3], dtype=float)
    origin = lengths / 2.0
    grouped = defaultdict(list)
    for residue, atom, xyz, resid in atoms:
        grouped[resid].append((residue, atom, np.asarray(xyz, dtype=float), resid))
    placed = []
    for items in grouped.values():
        base = items[0][2]
        local = []
        for residue, atom, xyz, resid in items:
            delta = xyz - base
            delta -= lengths * np.round(delta / lengths)
            local.append((residue, atom, base + delta, resid))
        heavy = [xyz for _residue, atom, xyz, _resid in local if not str(atom).startswith("H")]
        anchor = np.mean(heavy if heavy else [local[0][2]], axis=0)
        image = np.round((anchor - origin) / lengths)
        for residue, atom, xyz, resid in local:
            placed.append((residue, atom, xyz - image * lengths, resid))
    return placed


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
    positions, colors, names, radii, transparency = [], [], [], [], []
    residue_index = []
    oxygens = []
    guest_groups = {}
    for residue, atom, xyz, resid in atoms:
        if str(atom).startswith("H") and residue not in ("CH4", "THF") and not hydrogens:
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
    if cage_faces:
        meshes = cage_meshes(
            np.asarray(oxygens, dtype=float),
            [
                (group["kind"], np.asarray(group["heavy"], dtype=float))
                for group in guest_groups.values()
                if group["heavy"]
            ],
        )
    return finish(
        xyz, colors, radii, transparency, names, bonds, meshes,
        "cage", xyz.mean(axis=0), False,
    )


def draw(library, scene):
    if scene == "ice-ih":
        frame = read_dump(ROOT / "input/traj/ice1h.lammpstrj")
        return labelled_picture(frame, chill_labels(library, frame), "basal", False)
    if scene == "ice-ih-h":
        return genice_picture("1h", (2, 2, 1), [], False)
    if scene == "ice-ic":
        frame = read_dump(ROOT / "input/traj/mW_cubic.lammpstrj")
        labels = chill_labels(library, frame)
        frame, labels = central_slice(frame, labels, 700)
        return labelled_picture(frame, labels, "basal", False)
    if scene == "interface":
        frame = read_dump(ROOT / "input/traj/interface_f0.lammpstrj")
        return labelled_picture(frame, chill_labels(library, frame), "front", True)
    if scene == "nucleation":
        frame = read_dump(ROOT / "input/traj/nucleation_three.lammpstrj", frame_index=1)
        labels = chill_labels(library, frame)
        frame, labels, centre = largest_cluster(frame, labels)
        picture = labelled_picture(frame, labels, "cage", True, perspective=False, center=centre)
        return picture
    if scene == "brine":
        frame = read_dump(ROOT / "input/traj/nacl_front.lammpstrj")
        oxygen = [i for i, kind in enumerate(frame.types) if int(kind) == 1]
        labelled = chill_labels(library, frame.subset(oxygen))
        labels = [2] * len(frame.xyz)
        for index, label in zip(oxygen, labelled):
            labels[index] = label
        return labelled_picture(frame, labels, "front", True, ions=True)
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
    raise SystemExit(main(sys.argv[1:]))
