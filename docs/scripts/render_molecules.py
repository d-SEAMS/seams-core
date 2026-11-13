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

# Distinct elements so chemparseplot's CPK palette separates CHILL+ classes.
CHILL_ELEMENT = {
    0: "Ni",
    1: "Fe",
    2: "O",
    3: "S",
    4: "N",
    5: "Cl",
    6: "C",
    7: "P",
    8: "F",
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


def symbols_for(labels, highlight=None):
    if highlight is None:
        return [CHILL_ELEMENT.get(int(label), "C") for label in labels]
    return ["Fe" if i in highlight else "O" for i in range(len(labels))]


def legend_for(symbols, names):
    from ase.data.colors import jmol_colors
    from ase.data import atomic_numbers

    rows = []
    seen = set()
    for symbol, name in zip(symbols, names):
        if name in seen:
            continue
        seen.add(name)
        rgb = tuple(float(channel) for channel in jmol_colors[atomic_numbers[symbol]])
        rows.append((name, rgb))
    return rows


def render_atoms(atoms, renderer):
    _render_atoms = load_renderer()
    image = _render_atoms(
        atoms,
        renderer,
        0.3,
        "15x,25y,5z",
        canvas_size=1800,
        perspective_tilt=8.0,
    )
    return image


def to_rgb(image):
    import numpy as np
    from PIL import Image

    array = np.asarray(image)
    if array.dtype != np.uint8:
        array = np.clip(array, 0.0, 1.0)
        array = (array * 255.0).round().astype(np.uint8)
    if array.shape[-1] == 4:
        colour = array[..., :3].astype(np.float32)
        alpha = array[..., 3:4].astype(np.float32)
        if alpha.max() > 1.0:
            alpha = alpha / 255.0
        array = (colour * alpha + 255.0 * (1.0 - alpha)).round().astype(np.uint8)
    return Image.fromarray(array[..., :3], mode="RGB")


def compose(structure, legend):
    from PIL import Image, ImageDraw, ImageFont

    canvas = Image.new("RGB", (WIDTH, HEIGHT), (255, 255, 255))
    picture = to_rgb(structure)
    picture.thumbnail((1800, 1700), Image.Resampling.LANCZOS)
    canvas.paste(picture, (40, 40))
    draw = ImageDraw.Draw(canvas)
    font = ImageFont.load_default()
    for path in (
        "/usr/share/fonts/truetype/dejavu/DejaVuSans.ttf",
        "/usr/share/fonts/truetype/liberation/LiberationSans-Regular.ttf",
    ):
        if os.path.isfile(path):
            font = ImageFont.truetype(path, 28)
            break
    x = 40 + picture.size[0] + 36
    y = 48
    for name, rgb in legend:
        colour = tuple(int(channel * 255) for channel in rgb)
        draw.rectangle((x, y, x + 36, y + 36), fill=colour, outline=(20, 20, 20))
        draw.text((x + 48, y + 4), name, fill=(20, 20, 20), font=font)
        y += 52
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


def draw(library, scene):
    renderer = "ovito"
    if scene == "ice-ih":
        frame = read_dump(ROOT / "input/traj/ice1h.lammpstrj")
        labels = chill_labels(library, frame)
        symbols = symbols_for(labels)
        names = [CHILL_NAME[int(label)] for label in labels]
    elif scene == "ice-ic":
        frame = read_dump(ROOT / "input/traj/mW_cubic.lammpstrj")
        labels = chill_labels(library, frame)
        frame, labels = central_slice(frame, labels, 900)
        symbols = symbols_for(labels)
        names = [CHILL_NAME[int(label)] for label in labels]
    elif scene == "interface":
        frame = read_dump(ROOT / "input/traj/interface_f0.lammpstrj")
        labels = chill_labels(library, frame)
        symbols = symbols_for(labels)
        names = [CHILL_NAME[int(label)] for label in labels]
    elif scene == "nucleation":
        frame = read_dump(ROOT / "input/traj/nucleation_three.lammpstrj", frame_index=1)
        labels = chill_labels(library, frame)
        frame, labels = central_slice(frame, labels, 900)
        symbols = symbols_for(labels)
        names = [CHILL_NAME[int(label)] for label in labels]
    elif scene == "brine":
        frame = read_dump(ROOT / "input/traj/nacl_front.lammpstrj")
        mapping = {1: "O", 3: "Na", 4: "Cl"}
        naming = {1: "water", 3: "sodium", 4: "chloride"}
        symbols = [mapping.get(int(kind), "C") for kind in frame.types]
        names = [naming.get(int(kind), "other") for kind in frame.types]
    elif scene == "rings":
        frame = read_dump(ROOT / "input/traj/genice_sI.lammpstrj")
        members = ring_members(library, frame, {5, 6})
        symbols = ["Fe" if i in members else "O" for i in range(len(frame.xyz))]
        names = ["ring" if i in members else "other" for i in range(len(frame.xyz))]
    elif scene == "prisms":
        frame = read_dump(ROOT / "input/traj/ice1h.lammpstrj")
        members = ring_members(library, frame, {6})
        symbols = ["Fe" if i in members else "O" for i in range(len(frame.xyz))]
        names = ["six-ring" if i in members else "other" for i in range(len(frame.xyz))]
    elif scene in ("cage-si", "cage-sii", "cage-sh", "genice"):
        path = {
            "cage-si": "genice_sI.lammpstrj",
            "cage-sii": "genice_sII.lammpstrj",
            "cage-sh": "genice_sH.lammpstrj",
            "genice": "genice_fau.lammpstrj",
        }[scene]
        frame = read_dump(ROOT / "input/traj" / path)
        symbols = ["O"] * len(frame.xyz)
        names = ["host oxygen"] * len(frame.xyz)
    elif scene == "cage-guest":
        frame = guest_cage()
        symbols = frame.symbols
        names = ["guest" if symbol != "O" and symbol != "H" else ("hydrogen" if symbol == "H" else "water") for symbol in symbols]
    elif scene == "solvis":
        frame = read_dump(ROOT / "input/traj/genice_sI.lammpstrj")
        symbols = ["O"] * len(frame.xyz)
        names = ["oxygen"] * len(frame.xyz)
        renderer = "solvis"
    else:
        fail("unknown scene {}".format(scene))
    atoms = _as_atoms(frame, symbols)
    try:
        image = render_atoms(atoms, renderer)
    except Exception:
        if renderer != "ovito":
            image = render_atoms(atoms, "ovito")
        else:
            raise
    return compose(image, legend_for(symbols, names))


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
        "cage-sii": "GenIce sII host.",
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
