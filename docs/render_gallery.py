#!/usr/bin/env python3
# Copyright (c) d-SEAMS developers
"""Draw the science gallery from the slow paper tutorials.

Images are written when the docs are built. OVITO colours the atoms
when it imports. Matplotlib's 3D axes are the fallback. Nothing here
is meant to be committed.
"""

import importlib.util
import os
import subprocess
import sys
from ctypes import CDLL, POINTER, byref, c_char_p, c_double, c_int, c_void_p
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
IMAGE_DIR = ROOT / "docs" / "source" / "tutorials" / "gallery"
CAPTIONS = ROOT / "docs" / "source" / "tutorials" / "gallery-captions.rst"
RUNS = ROOT / "docs" / "build" / "gallery-runs"

LABEL_COLOR = {
    0: (0.12, 0.47, 0.71),
    1: (0.90, 0.45, 0.13),
    2: (0.65, 0.82, 0.89),
    3: (0.95, 0.78, 0.20),
    4: (0.55, 0.34, 0.70),
    5: (0.17, 0.62, 0.35),
    6: (0.60, 0.60, 0.60),
    7: (0.20, 0.30, 0.70),
    8: (0.70, 0.25, 0.15),
}
ION_COLOR = {2: (0.85, 0.15, 0.15), 3: (0.15, 0.65, 0.25)}


def have_stack(program):
    probe = subprocess.run(
        [program, "-c", "import openmm, numpy"],
        check=False,
        capture_output=True,
    )
    return probe.returncode == 0


def find_python():
    if have_stack(sys.executable):
        return None
    seen = []
    for folder in os.environ.get("PATH", "").split(os.pathsep):
        program = os.path.join(folder, "python3")
        if program in seen or not os.path.isfile(program):
            continue
        seen.append(program)
        if have_stack(program):
            return program
    return None


def paper_module():
    path = ROOT / "docs" / "tutorials" / "paper_science.py"
    spec = importlib.util.spec_from_file_location("paper_science", path)
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


def locate_tools(argv):
    args = [item for item in argv if not item.startswith("--")]
    if len(args) >= 3:
        return args[0], args[1], args[2]
    seams = os.environ.get("SEAMS", "")
    libdir = os.environ.get("SEAMS_LIBDIR", "")
    langdir = os.environ.get("SEAMS_LANGDIR", "")
    if seams and libdir and langdir:
        return seams, libdir, langdir
    for folder in (ROOT / "bbdir", ROOT / "build"):
        binary = folder / "seams"
        if binary.is_file():
            return str(binary), str(folder / "src"), str(folder / "bindings")
    raise SystemExit("gallery needs SEAMS, SEAMS_LIBDIR and SEAMS_LANGDIR")


def load_lang(libdir, langdir):
    os.environ["LD_LIBRARY_PATH"] = os.pathsep.join(
        [str(libdir), str(langdir), os.environ.get("LD_LIBRARY_PATH", "")]
    )
    lang = CDLL(str(Path(langdir) / "libseams_lang.so"))
    api = CDLL(str(Path(libdir) / "libseams.so"))
    lang.seams_frame_read.restype = c_void_p
    lang.seams_frame_read.argtypes = [c_char_p, c_int, c_int]
    lang.seams_frame_free.argtypes = [c_void_p]
    lang.seams_frame_nop.restype = c_int
    lang.seams_frame_nop.argtypes = [c_void_p]
    lang.seams_frame_positions.restype = POINTER(c_double)
    lang.seams_frame_positions.argtypes = [c_void_p]
    lang.seams_frame_types.restype = POINTER(c_int)
    lang.seams_frame_types.argtypes = [c_void_p]
    lang.seams_frame_chill_plus.restype = c_int
    lang.seams_frame_chill_plus.argtypes = [
        c_void_p, c_double, c_int, c_char_p, POINTER(c_int), POINTER(c_int),
    ]
    lang.seams_frame_steinhardt_graph.restype = c_int
    lang.seams_frame_steinhardt_graph.argtypes = [
        c_void_p, c_double, c_int, c_int, c_char_p, POINTER(c_double), POINTER(c_double),
    ]
    api.seams_system_create.restype = c_void_p
    api.seams_system_destroy.argtypes = [c_void_p]
    api.seams_system_set_frame.restype = c_int
    api.seams_system_set_frame.argtypes = [
        c_void_p, POINTER(c_double), c_int, POINTER(c_double), POINTER(c_int), POINTER(c_int),
    ]
    lang.seams_frame_box.argtypes = [c_void_p, POINTER(c_double)]
    api.seams_neighbours_build.restype = c_int
    api.seams_neighbours_build.argtypes = [c_void_p, c_double, c_int, c_int, POINTER(c_void_p)]
    api.seams_neighbours_destroy.argtypes = [c_void_p]
    api.seams_neighbours_csr.restype = c_int
    api.seams_neighbours_csr.argtypes = [
        c_void_p, POINTER(POINTER(c_int)), POINTER(c_int), POINTER(POINTER(c_int)), POINTER(c_int),
    ]
    return lang, api


def cloud(lang, path, frame, type_filter):
    import numpy as np

    handle = lang.seams_frame_read(str(path).encode(), frame, type_filter)
    if not handle:
        raise SystemExit("could not read {} frame {}".format(path.name, frame))
    count = lang.seams_frame_nop(handle)
    positions = np.ctypeslib.as_array(lang.seams_frame_positions(handle), shape=(count, 3)).copy()
    types = np.ctypeslib.as_array(lang.seams_frame_types(handle), shape=(count,)).copy()
    return handle, positions, types


def chill_labels(lang, handle, count):
    import numpy as np

    labels = (c_int * count)()
    counts = (c_int * 9)()
    if lang.seams_frame_chill_plus(handle, 3.5, 1, b"cutoff", labels, counts) != 0:
        raise SystemExit("CHILL+ failed")
    return np.ctypeslib.as_array(labels, shape=(count,)).copy()


def q6bar_values(lang, handle, count):
    import numpy as np

    local = (c_double * count)()
    averaged = (c_double * count)()
    if lang.seams_frame_steinhardt_graph(handle, 3.5, 1, 6, b"cutoff", local, averaged) != 0:
        raise SystemExit("q6bar failed")
    return np.ctypeslib.as_array(averaged, shape=(count,)).copy()


def undirected_bonds(lang, api, handle):
    count = lang.seams_frame_nop(handle)
    positions = lang.seams_frame_positions(handle)
    types = lang.seams_frame_types(handle)
    box = (c_double * 9)()
    lang.seams_frame_box(handle, box)
    system = api.seams_system_create()
    if api.seams_system_set_frame(system, positions, count, box, types, None) != 0:
        raise SystemExit("neighbour frame failed")
    built = c_void_p()
    if api.seams_neighbours_build(system, 3.5, 0, 0, byref(built)) != 0:
        raise SystemExit("neighbour build failed")
    offsets = POINTER(c_int)()
    columns = POINTER(c_int)()
    n_off = c_int()
    nnz = c_int()
    if api.seams_neighbours_csr(built, byref(offsets), byref(n_off), byref(columns), byref(nnz)) != 0:
        raise SystemExit("neighbour csr failed")
    pairs = []
    width = n_off.value
    for atom in range(max(0, width - 1)):
        for slot in range(offsets[atom], offsets[atom + 1]):
            other = columns[slot]
            if atom < other:
                pairs.append((atom, other))
    api.seams_neighbours_destroy(built)
    api.seams_system_destroy(system)
    return pairs


def color_from_labels(labels):
    import numpy as np

    return np.array([LABEL_COLOR.get(int(label), (0.5, 0.5, 0.5)) for label in labels], dtype=float)


def color_q6(values):
    import numpy as np

    span = np.clip((values - 0.2) / 0.45, 0.0, 1.0)
    colors = np.zeros((len(values), 3))
    colors[:, 0] = span
    colors[:, 1] = 0.25
    colors[:, 2] = 1.0 - span
    return colors


def render_scene(path, positions, colors, radii=None, bonds=None, bond_positions=None):
    paper_module().prepare_vulkan()
    try:
        import warnings

        warnings.filterwarnings("ignore", message=".*OVITO.*PyPI")
        import ovito
        from ovito.data import DataCollection, Particles
        from ovito.pipeline import Pipeline, StaticSource
        from ovito.vis import Viewport

        while ovito.scene.pipelines:
            ovito.scene.pipelines[0].remove_from_scene()
    except ImportError:
        return render_matplotlib(path, positions, colors, radii, bonds, bond_positions)
    try:
        data = DataCollection()
        particles = Particles()
        particles.create_property("Position", data=positions)
        particles.create_property("Color", data=colors)
        if radii is not None:
            particles.create_property("Radius", data=radii)
        if bonds is not None and len(bonds):
            particles.create_bonds().create_property("Topology", data=bonds)
        data.objects.append(particles)
        particles.vis.radius = 0.55
        pipeline = Pipeline(source=StaticSource(data=data))
        pipeline.add_to_scene()
        try:
            view = Viewport(type=Viewport.Type.Ortho)
            view.zoom_all()
            view.render_image(filename=str(path), size=(480, 360), background=(1, 1, 1))
        finally:
            pipeline.remove_from_scene()
    except Exception:
        return render_matplotlib(path, positions, colors, radii, bonds, bond_positions)
    if not path.is_file():
        return render_matplotlib(path, positions, colors, radii, bonds, bond_positions)
    return "ovito"


def render_matplotlib(path, positions, colors, radii, bonds, bond_positions):
    import matplotlib

    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    figure = plt.figure(figsize=(6.4, 4.8))
    axis = figure.add_subplot(111, projection="3d")
    sizes = 8.0 if radii is None else radii * 10.0
    axis.scatter(positions[:, 0], positions[:, 1], positions[:, 2], c=colors, s=sizes, depthshade=False)
    if bonds and bond_positions is not None:
        for left, right in bonds[:: max(1, len(bonds) // 4000)]:
            segment = bond_positions[[left, right]]
            axis.plot(segment[:, 0], segment[:, 1], segment[:, 2], color=(0.85, 0.65, 0.10), linewidth=0.4)
    axis.set_axis_off()
    figure.tight_layout()
    figure.savefig(path, dpi=100)
    plt.close(figure)
    return "matplotlib"


def save(path, positions, colors, radii=None, bonds=None, bond_positions=None):
    import numpy as np

    kind = render_scene(
        path,
        np.asarray(positions, dtype=float),
        np.asarray(colors, dtype=float),
        None if radii is None else np.asarray(radii, dtype=float),
        None if not bonds else np.asarray(bonds, dtype=int),
        None if bond_positions is None else np.asarray(bond_positions, dtype=float),
    )
    if not path.is_file() or path.stat().st_size == 0:
        raise SystemExit("missing figure {}".format(path.name))
    return kind


def brine_colors(lang, dump, frame):
    import numpy as np

    water, water_pos, _types = cloud(lang, dump, frame, 1)
    labels = chill_labels(lang, water, len(water_pos))
    lang.seams_frame_free(water)
    everyone, positions, types = cloud(lang, dump, frame, -1)
    lang.seams_frame_free(everyone)
    colors = np.zeros((len(positions), 3))
    radii = np.full(len(positions), 0.45)
    cursor = 0
    for index, kind in enumerate(types):
        if int(kind) in ION_COLOR:
            colors[index] = ION_COLOR[int(kind)]
            radii[index] = 1.15
        else:
            colors[index] = LABEL_COLOR.get(int(labels[cursor]), (0.7, 0.7, 0.7))
            cursor += 1
    return positions, colors, radii


def ring_scene(lang, api, dump):
    handle, positions, _types = cloud(lang, dump, 1, 1)
    try:
        bonds = undirected_bonds(lang, api, handle)
    finally:
        lang.seams_frame_free(handle)
    import numpy as np

    colors = np.tile((0.75, 0.84, 0.92), (len(positions), 1))
    return positions, colors, bonds


def caption_block(title, image_names, prose, report_text):
    lines = [title, "^" * len(title), ""]
    for name in image_names:
        lines.extend(
            [
                ".. figure:: gallery/{}".format(name),
                "   :alt: {}".format(title),
                "",
            ]
        )
    lines.append(prose)
    lines.append("")
    lines.append("::")
    lines.append("")
    for row in report_text.splitlines():
        if row.strip():
            lines.append("   {}".format(row.strip()))
    lines.append("")
    return "\n".join(lines)


def build(seams, libdir, langdir):
    import numpy as np

    science = paper_module()
    lang, api = load_lang(libdir, langdir)
    IMAGE_DIR.mkdir(parents=True, exist_ok=True)
    RUNS.mkdir(parents=True, exist_ok=True)
    reports = {}
    for topic, function in (
        ("interface", science.interface_full),
        ("brine", science.brine_full),
        ("rings", science.rings_full),
        ("nucleation", science.nucleation_full),
        ("cages", science.cages_full),
    ):
        out = RUNS / topic
        if out.exists():
            import shutil

            shutil.rmtree(out)
        out.mkdir(parents=True, exist_ok=True)
        function(seams, libdir, out)
        reports[topic] = (out / "report.txt").read_text().strip()

    interface = RUNS / "interface" / "interface-full.lammpstrj"
    handle, positions, _types = cloud(lang, interface, 1, 1)
    labels = chill_labels(lang, handle, len(positions))
    lang.seams_frame_free(handle)
    save(IMAGE_DIR / "gallery-interface.png", positions, color_from_labels(labels))

    brine = RUNS / "brine" / "brine-full.lammpstrj"
    for frame, name in ((1, "gallery-brine-start.png"), (3, "gallery-brine-end.png")):
        positions, colors, radii = brine_colors(lang, brine, frame)
        save(IMAGE_DIR / name, positions, colors, radii)

    rings = RUNS / "rings" / "rings-full.lammpstrj"
    positions, colors, bonds = ring_scene(lang, api, rings)
    save(IMAGE_DIR / "gallery-rings.png", positions, colors, bonds=bonds, bond_positions=positions)

    nucleus = RUNS / "nucleation" / "nucleation-full.lammpstrj"
    handle, positions, _types = cloud(lang, nucleus, 2, 1)
    values = q6bar_values(lang, handle, len(positions))
    lang.seams_frame_free(handle)
    save(IMAGE_DIR / "gallery-nucleation.png", positions, color_q6(values))

    for name, dump in (
        ("gallery-cages-sI.png", RUNS / "cages" / "sI-full.lammpstrj"),
        ("gallery-cages-sII.png", RUNS / "cages" / "sII-full.lammpstrj"),
    ):
        handle, positions, _types = cloud(lang, dump, 1, 1)
        try:
            bonds = undirected_bonds(lang, api, handle)
        finally:
            lang.seams_frame_free(handle)
        tint = (0.55, 0.34, 0.70) if "sI" in name else (0.15, 0.45, 0.70)
        colors = np.tile(tint, (len(positions), 1))
        save(IMAGE_DIR / name, positions, colors, bonds=bonds, bond_positions=positions)

    def piece(line, key):
        parts = line.split()
        return parts[parts.index(key) + 1]

    interface_lines = reports["interface"].splitlines()
    published = interface_lines[1]
    ovito = interface_lines[2]
    brine_lines = reports["brine"].splitlines()
    ring_line = reports["rings"].splitlines()[-1]
    q_lines = reports["nucleation"].splitlines()
    cage_lines = reports["cages"].splitlines()
    front_start = piece(brine_lines[0], "front")
    front_end = piece(brine_lines[2], "front")
    blocks = [
        ".. Generated by docs/render_gallery.py from the slow tutorial reports.",
        "",
        caption_block(
            "Interface slab",
            ["gallery-interface.png"],
            " ".join(
                [
                    "The picture is the 2048-water slab, coloured by published CHILL+.",
                    published.replace("interface full published ", "Published CHILL+ reports ") + ".",
                    ovito.replace("interface full ovito ", "The OVITO rule reports ") + ".",
                    "Hexagonal ice is the warm colour and liquid water is pale. The interfacial shell is where the two rules part.",
                ]
            ),
            reports["interface"],
        ),
        caption_block(
            "Brine front at 240 K",
            ["gallery-brine-start.png", "gallery-brine-end.png"],
            " ".join(
                [
                    "The first picture is the start of the 240 K run and the second picture is the end.",
                    "Sodium is red and chloride is green.",
                    brine_lines[0].replace("brine full frame 0 ", "At the start, ") + ".",
                    brine_lines[2].replace("brine full frame 2 ", "At the end, ") + ".",
                    "The front moves from {} Å to {} Å as hexagonal water fills the next lattice planes.".format(
                        front_start, front_end
                    ),
                ]
            ),
            reports["brine"],
        ),
        caption_block(
            "Ice Ih six-rings",
            ["gallery-rings.png"],
            " ".join(
                [
                    "The cell is 4096 molecules of ice Ih.",
                    "The drawn bonds are the cutoff neighbours, and in this cell they close six-rings.",
                    ring_line.replace("rings full histogram ", "The histogram is ") + ".",
                    "Bulk ice Ih has two six-rings per molecule and no other size through 8.",
                ]
            ),
            reports["rings"],
        ),
        caption_block(
            "Nucleation nucleus",
            ["gallery-nucleation.png"],
            " ".join(
                [
                    "Each molecule is coloured by its own q6bar on the frame that restores the basal nucleus.",
                    "Blue is liquid-like and red is ice-like.",
                    q_lines[0].replace("nucleation full frame 0 ", "Frame 0 has ") + ".",
                    q_lines[1].replace("nucleation full frame 1 ", "Frame 1 has ") + ".",
                    q_lines[2].replace("nucleation full frame 2 ", "Frame 2 has ") + ".",
                    "The ice threshold is 0.36, and the last frame is past it.",
                ]
            ),
            reports["nucleation"],
        ),
        caption_block(
            "sI and sII cages",
            ["gallery-cages-sI.png", "gallery-cages-sII.png"],
            " ".join(
                [
                    "The first cluster is sI and the second cluster is sII.",
                    "Bonds draw the water rings that close the cages.",
                    "sI holds {} of the small cage and {} of the large cage.".format(
                        cage_lines[0].split()[3], cage_lines[0].split()[4]
                    ),
                    "sII holds {} of the small cage and {} of the large cage.".format(
                        cage_lines[1].split()[3], cage_lines[1].split()[4]
                    ),
                ]
            ),
            reports["cages"],
        ),
    ]
    CAPTIONS.write_text("\n".join(blocks))
    return reports


def check(reports):
    text = CAPTIONS.read_text()
    errors = []
    for topic, report in reports.items():
        for line in report.splitlines():
            if line.strip() and line.strip() not in text:
                errors.append("caption missing {} line {}".format(topic, line.strip()))
    for name in (
        "gallery-interface.png",
        "gallery-brine-start.png",
        "gallery-brine-end.png",
        "gallery-rings.png",
        "gallery-nucleation.png",
        "gallery-cages-sI.png",
        "gallery-cages-sII.png",
    ):
        image = IMAGE_DIR / name
        if not image.is_file() or image.stat().st_size == 0:
            errors.append("missing {}".format(name))
        elif name not in text:
            errors.append("caption does not show {}".format(name))
    return errors


def main(argv):
    other = find_python()
    if other:
        done = subprocess.run([other, str(Path(__file__).resolve()), *argv], check=False)
        return done.returncode
    seams, libdir, langdir = locate_tools(argv)
    reports = build(seams, libdir, langdir)
    if "--check" in argv:
        errors = check(reports)
        if errors:
            print("\n".join(errors), file=sys.stderr)
            return 1
        print("gallery {} images".format(len(list(IMAGE_DIR.glob('*.png')))))
    return 0


if __name__ == "__main__":
    raise SystemExit(main(sys.argv[1:]))
