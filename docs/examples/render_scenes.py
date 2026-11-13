"""Tachyon scenes for the brine front and the hydrate cages.

Waters use the CHILL+ colours from scripts/science_render.py. Sodium
and chloride are extra spheres. Cage faces are polyhedra in those same
colours: 5^12 is hexagonal orange, 5^12 6^4 is cubic blue, the sH
medium cage is interfacial yellow, and 5^12 6^8 is clathrate green.
"""

from __future__ import annotations

import ctypes
import sys
from collections import defaultdict
from pathlib import Path

import numpy as np

_SCRIPTS = Path(__file__).resolve().parents[2] / "scripts"
if str(_SCRIPTS) not in sys.path:
    sys.path.insert(0, str(_SCRIPTS))
from chill_palette import CHLORIDE, COLOUR, LEGEND, SODIUM
from lammps_box import cell_from_bounds

# Face-size census of each cage, and the CHILL+ colour index it borrows.
CAGE_KINDS = (
    ("512", {5: 12}, 1),
    ("51264", {5: 12, 6: 4}, 0),
    ("sh", {4: 3, 5: 6, 6: 3}, 3),
    ("51268", {5: 12, 6: 8}, 4),
)
NA = np.array(SODIUM)
CL = np.array(CHLORIDE)
BRINE_COLOUR = {name: np.array(colour) for name, colour in COLOUR.items()}
BRINE_LEGEND = tuple((label, np.array(colour)) for label, colour in LEGEND)
# Drawn oxygen is 0.90 for a 1.4 Å oxygen. Sodium and chloride use that
# same scale at 1.0 Å and 1.8 Å. Liquid stays the small veil.
_OXYGEN_A = 1.4
BRINE_RADIUS = {
    "O": 0.90,
    "water": 0.50,
    "Na": 0.90 * (1.0 / _OXYGEN_A),
    "Cl": 0.90 * (1.8 / _OXYGEN_A),
}
# One size for the interface legend and the brine legend.
LEGEND_FONT_PX = 44
MAX_FACE_EDGE = 3.5
MAX_PLANARITY = 0.35


def _rings(lib, path: str):
    lib.seams_rings_find.argtypes = [
        ctypes.c_void_p, ctypes.c_double, ctypes.c_int, ctypes.c_int, ctypes.c_int,
        ctypes.POINTER(ctypes.c_void_p),
    ]
    lib.seams_rings_find.restype = ctypes.c_int
    lib.seams_rings_count.argtypes = [ctypes.c_void_p]
    lib.seams_rings_count.restype = ctypes.c_int
    lib.seams_rings_vertices.argtypes = [
        ctypes.c_void_p, ctypes.c_int, ctypes.POINTER(ctypes.c_int), ctypes.c_size_t,
        ctypes.POINTER(ctypes.c_size_t),
    ]
    lib.seams_rings_vertices.restype = ctypes.c_int
    lib.seams_rings_destroy.argtypes = [ctypes.c_void_p]
    system = lib.seams_system_create()
    if lib.seams_system_read(system, path.encode(), 1, -1) != 0:
        raise RuntimeError(lib.seams_last_error())
    rings = ctypes.c_void_p()
    if lib.seams_rings_find(system, 3.5, 0, 0, 6, ctypes.byref(rings)) != 0:
        raise RuntimeError(lib.seams_last_error())
    count = lib.seams_rings_count(rings)
    cycles = []
    for index in range(count):
        buf = (ctypes.c_int * 16)()
        n_out = ctypes.c_size_t()
        if lib.seams_rings_vertices(rings, index, buf, 16, ctypes.byref(n_out)) != 0:
            raise RuntimeError("ring vertices")
        cycles.append(tuple(int(buf[k]) for k in range(n_out.value)))
    lib.seams_rings_destroy(rings)
    lib.seams_system_destroy(system)
    return cycles


def find_cages(cycles, need, limit):
    """Closed trivalent polyhedra whose face sizes match ``need``."""
    edge_of = {}
    face_edges = []
    edge_faces = defaultdict(list)
    for face, cycle in enumerate(cycles):
        edges = []
        for index in range(len(cycle)):
            edge = tuple(sorted((cycle[index], cycle[(index + 1) % len(cycle)])))
            edges.append(edge_of.setdefault(edge, len(edge_of)))
        face_edges.append(edges)
        for edge in edges:
            edge_faces[edge].append(face)
    n_edges = len(edge_of)
    ecount = [0] * n_edges
    vcount = [0] * (max(max(cycle) for cycle in cycles) + 1)
    used = [False] * len(cycles)
    sizes = [len(cycle) for cycle in cycles]
    target = sum(need.values())
    have = {size: 0 for size in need}
    chosen = []
    found = []
    seen = set()

    def legal(face):
        size = sizes[face]
        if used[face] or need.get(size, 0) <= have.get(size, 0):
            return False
        if any(ecount[edge] >= 2 for edge in face_edges[face]):
            return False
        return all(vcount[atom] < 3 for atom in cycles[face])

    def add(face):
        used[face] = True
        chosen.append(face)
        have[sizes[face]] += 1
        for edge in face_edges[face]:
            ecount[edge] += 1
        for atom in cycles[face]:
            vcount[atom] += 1

    def remove(face):
        used[face] = False
        chosen.pop()
        have[sizes[face]] -= 1
        for edge in face_edges[face]:
            ecount[edge] -= 1
        for atom in cycles[face]:
            vcount[atom] -= 1

    def search():
        if len(found) >= limit:
            return
        if len(chosen) == target:
            closed = all(ecount[edge] == 2 for face in chosen for edge in face_edges[face])
            trivalent = all(vcount[atom] == 3 for face in chosen for atom in cycles[face])
            if closed and trivalent:
                key = frozenset(chosen)
                if key not in seen:
                    seen.add(key)
                    found.append(list(chosen))
            return
        boundary = []
        marked = set()
        for face in chosen:
            for edge in face_edges[face]:
                if ecount[edge] == 1 and edge not in marked:
                    marked.add(edge)
                    boundary.append([other for other in edge_faces[edge] if legal(other)])
        if not boundary or any(len(options) == 0 for options in boundary):
            return
        forced = []
        branch = None
        for options in boundary:
            if len(options) == 1:
                if options[0] not in forced:
                    forced.append(options[0])
            elif branch is None or len(options) < len(branch):
                branch = options
        if forced:
            face = forced[0]
            if legal(face):
                add(face)
                search()
                remove(face)
            return
        for face in branch:
            add(face)
            search()
            remove(face)
            if len(found) >= limit:
                return

    for seed, size in enumerate(sizes):
        if size not in need or len(found) >= limit:
            continue
        add(seed)
        search()
        remove(seed)
    return found


def _box_vectors(path):
    with open(path) as handle:
        for line in handle:
            if not line.startswith("ITEM: BOX"):
                continue
            rows = [handle.readline().split() for _ in range(3)]
            _origin, (lx, ly, lz), (xy, xz, yz) = cell_from_bounds(rows)
            return [
                np.array([lx, 0.0, 0.0]),
                np.array([xy, ly, 0.0]),
                np.array([xz, yz, lz]),
            ]
    raise RuntimeError(f"no box in {path}")


def _shifts(vectors):
    return [
        i * vectors[0] + j * vectors[1] + k * vectors[2]
        for i in (-1, 0, 1)
        for j in (-1, 0, 1)
        for k in (-1, 0, 1)
    ]


def _mic(point, target, shifts):
    return point + min(shifts, key=lambda shift: np.linalg.norm(point + shift - target))


def cage_centre(atom_ids, positions, vectors):
    """Mean of each atom's image nearest a seed atom."""
    shifts = _shifts(vectors)
    seed = positions[atom_ids[0]]
    placed = [_mic(positions[atom], seed, shifts) for atom in atom_ids]
    return np.mean(placed, axis=0)


def face_metrics(coords):
    """Longest edge and the largest distance from the best-fit plane."""
    rolled = np.roll(coords, -1, axis=0)
    longest = float(np.linalg.norm(rolled - coords, axis=1).max())
    centred = coords - coords.mean(axis=0)
    _u, _s, vh = np.linalg.svd(centred, full_matrices=False)
    normal = vh[-1]
    residual = float(np.abs(centred @ normal).max())
    return longest, residual


def drawn_ring(cycle, positions, centre, vectors):
    """Place one ring on the images nearest the cage centre.

    A face whose longest edge reaches 3.5 angstrom, or whose planarity
    residual exceeds 0.35 angstrom, is not drawn. A non-planar ring is
    split into triangles from its first vertex, and only the triangles
    that pass the same tests are kept.
    """
    shifts = _shifts(vectors)
    coords = np.vstack([_mic(positions[atom], centre, shifts) for atom in cycle])
    longest, residual = face_metrics(coords)
    if longest < MAX_FACE_EDGE and residual <= MAX_PLANARITY:
        return [coords], [(longest, residual)]
    kept = []
    stats = []
    for index in range(1, len(coords) - 1):
        triangle = coords[[0, index, index + 1]]
        edge, _flat = face_metrics(triangle)
        if edge < MAX_FACE_EDGE:
            kept.append(triangle)
            stats.append((edge, 0.0))
    return kept, stats


def cage_meshes(cycles, groups, positions, vectors):
    """One vertex array and triangle list per colour index.

    Each ring is unwrapped onto the images nearest the cage centre.
    The returned stats are ``(longest edge, planarity residual)`` for
    every polygon that is actually drawn.
    """
    by_colour = defaultdict(list)
    reference = None
    shifts = _shifts(vectors)
    stats = []
    for faces, colour in groups:
        atoms = sorted({atom for face in faces for atom in cycles[face]})
        centre = cage_centre(atoms, positions, vectors)
        polygons = []
        for face in faces:
            pieces, piece_stats = drawn_ring(cycles[face], positions, centre, vectors)
            for coords in pieces:
                normal = np.cross(coords[1] - coords[0], coords[2] - coords[0])
                if float(np.dot(normal, coords.mean(axis=0) - centre)) < 0.0:
                    coords = coords[::-1]
                polygons.append(coords)
            stats.extend(piece_stats)
        if not polygons:
            continue
        cage_centre_now = np.mean([poly.mean(axis=0) for poly in polygons], axis=0)
        if reference is None:
            reference = cage_centre_now
            delta = np.zeros(3)
        else:
            delta = min(shifts, key=lambda shift: np.linalg.norm(cage_centre_now + shift - reference))
        by_colour[colour].append([poly + delta for poly in polygons])
    meshes = []
    for colour, cages in by_colour.items():
        vertices = []
        triangles = []
        for cage in cages:
            for polygon in cage:
                base = len(vertices)
                vertices.extend(polygon)
                for index in range(1, len(polygon) - 1):
                    triangles.append([base, base + index, base + index + 1])
        meshes.append((np.asarray(vertices, dtype=float), triangles, colour))
    return meshes, stats


def assert_drawn_faces(stats):
    """Every drawn face is a short, nearly planar polygon."""
    if not stats:
        raise RuntimeError("no cage faces were drawn")
    for longest, residual in stats:
        if longest >= MAX_FACE_EDGE or residual > MAX_PLANARITY:
            raise RuntimeError(f"drawn face edge {longest:.3f} A residual {residual:.3f} A")


def _camera(render, positions, forward, up):
    centre, half_h, half_v, forward, right, up = render.project(positions, forward, up)
    return centre, half_h, half_v, forward, right, up


def tachyon_scene(render, path, positions, colours, radii, meshes, forward, up, width, frame=None):
    """One Tachyon frame. Settings match science_render._paint."""
    from ovito.data import DataCollection, Particles, SimulationCell, SurfaceMesh
    from ovito.pipeline import Pipeline, StaticSource
    from ovito.vis import TachyonRenderer, Viewport

    centre, half_h, half_v, forward, _right, up = _camera(render, frame if frame is not None else positions, forward, up)
    height = max(int(round(width * half_v / max(half_h, 1e-6))), 32)
    data = DataCollection()
    low = positions.min(axis=0) - 6.0
    high = positions.max(axis=0) + 6.0
    cell = SimulationCell(pbc=(False, False, False))
    cell[:, :3] = np.diag((high - low).astype(float))
    cell[:, 3] = low
    cell.vis.enabled = False
    data.objects.append(cell)
    tinted, near = render.depth_cue(
        colours, positions @ forward, fade=0.38, tint=np.array((0.94, 0.94, 0.96))
    )
    radii = np.asarray(radii, dtype=np.float64) * (0.84 + 0.16 * near)
    particles = Particles()
    particles.create_property("Position", data=np.ascontiguousarray(positions, dtype=np.float64))
    particles.create_property("Color", data=np.clip(tinted, 0.0, 1.0))
    particles.create_property("Radius", data=radii)
    particles.create_property("Transparency", data=np.zeros(len(positions)))
    data.objects.append(particles)
    particles.vis.radius = float(np.median(radii))
    for vertices, faces, colour in meshes:
        mesh = SurfaceMesh()
        mesh.create_vertices(np.ascontiguousarray(vertices, dtype=np.float64))
        mesh.create_faces(faces)
        red, green, blue = render.CHILL[colour][1]
        mesh.vis.surface_color = (float(red), float(green), float(blue))
        mesh.vis.surface_transparency = 0.08
        mesh.vis.highlight_edges = False
        mesh.vis.cap_color = mesh.vis.surface_color
        data.objects.append(mesh)
    pipeline = Pipeline(source=StaticSource(data=data))
    view = Viewport(type=Viewport.Type.Ortho)
    pipeline.add_to_scene()
    try:
        view.camera_dir = tuple(float(value) for value in forward)
        view.camera_up = tuple(float(value) for value in up)
        view.camera_pos = centre - forward * (8.0 * max(half_h, half_v))
        view.fov = float(half_v)
        renderer = TachyonRenderer(
            shadows=True,
            ambient_occlusion=True,
            ambient_occlusion_samples=2,
            antialiasing_samples=3,
            direct_light_intensity=0.8,
        )
        path.parent.mkdir(parents=True, exist_ok=True)
        view.render_image(
            filename=str(path),
            size=(width, height),
            renderer=renderer,
            background=(1.0, 1.0, 1.0),
            frame=0,
        )
    finally:
        pipeline.remove_from_scene()


def render_brine(render, rows, path, box, width=1400):
    """Oxygens by the paper CHILL+ palette. Na is violet and Cl is black.

    Hydrogens are omitted. Clathrate and interfacial clathrate are the
    light-blue and pink labels from render_front.py, and both appear in
    the legend.
    """
    kept = [row for row in rows if row[0] in ("O", "Na", "Cl")]
    origin = np.array(box[:3], dtype=float)
    lengths = np.array(box[3:6], dtype=float)
    positions = wrap_into_cell(
        np.array([[row[1], row[2], row[3]] for row in kept], dtype=float),
        origin,
        lengths,
    )
    colours = []
    radii = []
    present = set()
    for species, _x, _y, _z, label in kept:
        if species == "Na":
            colours.append(NA)
            radii.append(BRINE_RADIUS["Na"])
            present.add("Na+")
        elif species == "Cl":
            colours.append(CL)
            radii.append(BRINE_RADIUS["Cl"])
            present.add("Cl-")
        else:
            colours.append(BRINE_COLOUR.get(label, BRINE_COLOUR["unclassified"]))
            radii.append(BRINE_RADIUS["O"] if label != "water" else BRINE_RADIUS["water"])
            present.add("liquid" if label == "water" else ("interf. clathrate" if label == "interClathrate" else label))
    colours = np.vstack(colours)
    radii = np.asarray(radii)
    # Liquid is a light veil so the ions stay visible. The interface slab keeps 0.40.
    liquid = np.array([row[0] == "O" and row[4] == "water" for row in kept])
    solid = ~liquid
    ice_png = path.with_name(path.stem + "-ice.png")
    liq_png = path.with_name(path.stem + "-liq.png")
    forward = np.array([0.0, 1.0, 0.0])
    up = np.array([0.0, 0.0, 1.0])
    pad_radius = float(np.max(radii))
    frame = padded_frame(positions, origin, lengths, forward, up, pad_radius)
    tachyon_scene(render, ice_png, positions[solid], colours[solid], radii[solid], [], forward, up, width, frame=frame)
    tachyon_scene(render, liq_png, positions[liquid], colours[liquid], radii[liquid], [], forward, up, width, frame=frame)
    from PIL import Image

    ice = np.asarray(Image.open(ice_png).convert("RGB"), dtype=np.float32)
    liquid_image = np.asarray(Image.open(liq_png).convert("RGB"), dtype=np.float32)
    veil = np.any(liquid_image < 250, axis=2)
    blended = ice.copy()
    blended[veil] = 0.18 * liquid_image[veil] + 0.82 * ice[veil]
    Image.fromarray(np.clip(blended, 0, 255).astype(np.uint8)).save(path)
    corners = cell_corners(origin, lengths)
    legend = [item for item in BRINE_LEGEND if item[0] in present]
    ice_z = np.array(
        [positions[index, 2] for index, row in enumerate(kept) if row[0] == "O" and row[4] == "hexagonal"]
    )
    annotate(
        render, path, positions, forward, up, width, corners, legend,
        frame=frame, origin=origin, lengths=lengths, ice_z=ice_z,
    )


def render_host(render, lib, dump, path, kinds, width=1280):
    cycles = _rings(lib, str(dump))
    positions = _positions(dump)
    groups = []
    expected = []
    for _name, need, colour, limit in kinds:
        facesets = find_cages(cycles, need, limit)
        expected.append(len(facesets))
        for faces in facesets:
            groups.append((faces, colour))
    vectors = _box_vectors(dump)
    meshes, stats = cage_meshes(cycles, groups, positions, vectors)
    assert_drawn_faces(stats)
    corners = np.vstack([vertices for vertices, _faces, _colour in meshes])
    colours = np.tile(np.array((0.96, 0.96, 0.96)), (len(corners), 1))
    radii = np.full(len(corners), 0.22)
    forward = np.array([0.55, 0.35, 0.76])
    up = np.array([0.0, 0.0, 1.0])
    tachyon_scene(render, path, corners, colours, radii, meshes, forward, up, width)
    return expected


def _legend_font(size):
    from PIL import ImageFont

    try:
        return ImageFont.truetype("/usr/share/fonts/truetype/dejavu/DejaVuSans.ttf", size)
    except OSError:
        return ImageFont.load_default()


def display_label(label):
    """Legend text. Sodium and chloride use real superscripts."""
    if label == "Na+":
        return "Na\u207a"
    if label == "Cl-":
        return "Cl\u207b"
    return label


def text_bbox(font, xy, text):
    left, top, right, bottom = font.getbbox(text)
    x, y = xy
    return (x + left, y + top, x + right, y + bottom)


def marks_overlap(marks):
    """Pairs of named boxes whose interiors meet."""
    hits = []
    for index, (name, box) in enumerate(marks):
        for other, obox in marks[index + 1 :]:
            if box[0] < obox[2] and obox[0] < box[2] and box[1] < obox[3] and obox[1] < box[3]:
                hits.append((name, other))
    return hits


def require_disjoint(marks):
    hits = marks_overlap(marks)
    if hits:
        detail = ", ".join(f"{left} overlaps {right}" for left, right in hits)
        raise RuntimeError(detail)


def annotate(render, path, positions, forward, up, width, corners, legend=None, layer_title=None, frame=None, origin=None, lengths=None, ice_z=None):
    """Box outline, a 1 nm bar, and the colour legend under the atoms.

    The bar is below the legend. When ``ice_z`` is given, a bracket marks
    the ice and an arrow marks the growth front.
    """
    from PIL import Image, ImageDraw

    centre, half_h, half_v, forward, right, up = _camera(render, frame if frame is not None else positions, forward, up)
    image = Image.open(path).convert("RGB")
    font = _legend_font(LEGEND_FONT_PX)
    swatch = 32
    row_h = 52
    legend_rows = 1
    labels = [(display_label(label), colour) for label, colour in (legend or ())]
    if labels:
        x = 16
        rows = 1
        for label, _colour in labels:
            width_px = swatch + 16 + float(font.getlength(label))
            if x + width_px > image.width - 8:
                rows += 1
                x = 16
            x += width_px + 18
        legend_rows = rows
    margin = 24 + row_h * legend_rows + LEGEND_FONT_PX + 48
    left = 170 if ice_z is not None else (110 if layer_title else 0)
    right_margin = 500 if ice_z is not None else 0
    top = 8 if ice_z is not None else 0
    canvas = Image.new("RGB", (image.width + left + right_margin, image.height + top + margin), (255, 255, 255))
    canvas.paste(image, (left, top))
    draw = ImageDraw.Draw(canvas)
    marks = []

    def add_text(name, xy, text):
        marks.append((name, text_bbox(font, xy, text)))
        draw.text(xy, text, font=font, fill=(20, 20, 20))

    if layer_title:
        text = Image.new("RGB", (int(font.getlength(layer_title)) + 8, 52), (255, 255, 255))
        ImageDraw.Draw(text).text((0, 4), layer_title, font=font, fill=(20, 20, 20))
        turned = text.rotate(90, expand=True)
        canvas.paste(turned, (6, max(8, (image.height - turned.height) // 2)))

    def project(point):
        offset = np.asarray(point, dtype=float) - centre
        x = left + (float(np.dot(offset, right)) / (2.0 * half_h) + 0.5) * image.width
        y = top + (0.5 - float(np.dot(offset, up)) / (2.0 * half_v)) * image.height
        return (x, y)

    projected = [project(corner) for corner in corners]
    for first in range(8):
        for second in range(first + 1, 8):
            if bin(first ^ second).count("1") != 1:
                continue
            draw.line([projected[first], projected[second]], fill=(70, 70, 70), width=2)
    if ice_z is not None:
        z_cut = connected_slab_top(ice_z)
        z_lo = float(origin[2])
        z_hi = z_lo + float(lengths[2])

        def y_of(z_value):
            point = np.asarray(origin, dtype=float) + np.array((0.5 * lengths[0], 0.5 * lengths[1], z_value - z_lo))
            return project(point)[1]

        x_left = min(point[0] for point in projected)
        x_right = max(point[0] for point in projected)
        brine = "brine"
        brine_w = float(font.getlength(brine))
        add_text("brine", (max(8, x_left - brine_w - 16), y_of(0.5 * (z_cut + z_hi)) - LEGEND_FONT_PX / 2.0), brine)
        bracket_x = x_right + 28
        y_front = y_of(z_cut)
        y_ice = y_of(z_lo)
        draw.line([(bracket_x, y_front), (bracket_x, y_ice)], fill=(20, 20, 20), width=3)
        add_text("ice Ih", (bracket_x + 14, 0.5 * (y_front + y_ice) - LEGEND_FONT_PX / 2.0), "ice Ih")
        draw.polygon([(x_right + 4, y_front), (x_right + 22, y_front - 7), (x_right + 22, y_front + 7)], fill=(20, 20, 20))
        draw.line([(bracket_x - 8, y_front), (bracket_x + 8, y_front)], fill=(20, 20, 20), width=3)
        add_text("growth front", (bracket_x + 14, y_front - LEGEND_FONT_PX - 4), "growth front")
    legend_bottom = image.height + top
    if labels:
        x, y = left + 16, image.height + top + 16
        for label, colour in labels:
            text_w = float(font.getlength(label))
            if x + swatch + 16 + text_w > canvas.width - 8:
                x = left + 16
                y += row_h
            rgb = tuple(int(round(channel * 255)) for channel in colour)
            draw.rectangle([x, y, x + swatch, y + swatch], fill=rgb, outline=(30, 30, 30))
            marks.append((f"swatch {label}", (x, y, x + swatch, y + swatch)))
            add_text(f"legend {label}", (x + swatch + 10, y), label)
            legend_bottom = max(legend_bottom, y + swatch, marks[-1][1][3])
            x += swatch + 16 + text_w + 18
    bar = 10.0 * image.width / (2.0 * half_h)
    bar_y = legend_bottom + 28
    bar_x = left + 16
    draw.line([(bar_x, bar_y), (bar_x + bar, bar_y)], fill=(20, 20, 20), width=6)
    marks.append(("1 nm bar", (bar_x, bar_y - 4, bar_x + bar, bar_y + 4)))
    add_text("1 nm", (bar_x + bar + 14, bar_y - LEGEND_FONT_PX / 2.0), "1 nm")
    require_disjoint(marks)
    canvas.save(path)


def basal_bilayers(positions, labels, axis, gap=0.8):
    """Index crystalline basal bilayers along ``axis``.

    A gap larger than ``gap`` angstrom starts a new sheet. Two sheets
    whose centres are closer than 2.2 angstrom are the puckered halves
    of one bilayer and share an index. A bilayer is numbered only when
    it holds at least half the oxygens of the fullest one. Cubic and
    interfacial oxygens stay at -1.
    """
    ids = np.full(len(positions), -1, dtype=int)
    axis = np.asarray(axis, dtype=float)
    axis = axis / np.linalg.norm(axis)
    chosen = np.flatnonzero(np.asarray(labels) == 1)
    if chosen.size == 0:
        return ids, []
    projection = positions @ axis
    order = np.argsort(projection[chosen])
    ordered = chosen[order]
    values = projection[ordered]
    breaks = np.flatnonzero(np.diff(values) > gap) + 1
    sheets = [sheet for sheet in np.split(ordered, breaks) if len(sheet)]
    merged = []
    for sheet in sheets:
        centre = float(np.median(projection[sheet]))
        if merged and centre - merged[-1][0] < 2.2:
            members = np.concatenate([merged[-1][1], sheet])
            merged[-1] = (float(np.median(projection[members])), members)
        else:
            merged.append((centre, sheet))
    fullest = max(len(members) for _centre, members in merged)
    kept = [item for item in merged if len(item[1]) >= 0.5 * fullest]
    centres = []
    for index, (_centre, members) in enumerate(kept):
        ids[members] = index
        centres.append(float(np.median(projection[members])))
    return ids, centres


def wrap_into_cell(positions, origin, lengths):
    """Fold coordinates into the periodic cell."""
    origin = np.asarray(origin, dtype=float)
    lengths = np.asarray(lengths, dtype=float)
    return origin + np.mod(np.asarray(positions, dtype=float) - origin, lengths)


def cell_corners(origin, lengths):
    origin = np.asarray(origin, dtype=float)
    lengths = np.asarray(lengths, dtype=float)
    return [
        origin + np.array((i * lengths[0], j * lengths[1], k * lengths[2]))
        for i in (0, 1)
        for j in (0, 1)
        for k in (0, 1)
    ]


def padded_frame(positions, origin, lengths, forward, up, radius):
    """Centres plus the cell corners pushed out by one sphere radius."""
    import science_render as science

    _forward, right, up_axis = science.basis(np.asarray(forward, dtype=float), np.asarray(up, dtype=float))
    extra = []
    for corner in cell_corners(origin, lengths):
        extra.append(corner + radius * right)
        extra.append(corner - radius * right)
        extra.append(corner + radius * up_axis)
        extra.append(corner - radius * up_axis)
    return np.vstack([np.asarray(positions, dtype=float), np.asarray(extra, dtype=float)])


def connected_slab_top(z_values, gap=8.0):
    """Top of the ice slab. A periodic image split off by ``gap`` is left out."""
    z = np.sort(np.asarray(z_values, dtype=float))
    if z.size == 0:
        raise ValueError("the ice slab is empty")
    breaks = np.flatnonzero(np.diff(z) > gap)
    end = int(breaks[0]) + 1 if breaks.size else int(z.size)
    return float(z[end - 1])


def centres_outside_outline(positions, origin, lengths, forward=None, up=None):
    """How many atom centres project outside the orthographic cell rectangle."""
    import science_render as science

    if forward is None:
        forward = np.array([-1.0, 0.0, 0.0])
    if up is None:
        up = np.array([0.0, 0.0, 1.0])
    _forward, right, up = science.basis(np.asarray(forward, dtype=float), np.asarray(up, dtype=float))
    origin = np.asarray(origin, dtype=float)
    lengths = np.asarray(lengths, dtype=float)
    face = cell_corners(origin, lengths)
    horiz = np.array([float(np.dot(point - origin, right)) for point in face])
    vert = np.array([float(np.dot(point - origin, up)) for point in face])
    point_h = (np.asarray(positions, dtype=float) - origin) @ right
    point_v = (np.asarray(positions, dtype=float) - origin) @ up
    tol = 1e-4
    outside = (
        (point_h < horiz.min() - tol)
        | (point_h > horiz.max() + tol)
        | (point_v < vert.min() - tol)
        | (point_v > vert.max() + tol)
    )
    return int(np.count_nonzero(outside))


def _present_interface(render, path, positions, labels, origin, lengths, forward, up, ids, legend, frame=None, pad_radius=0.85):
    """Rectangle cell, bottom bilayer ticks, and labels outside the box.

    The camera looks down x with up along z, so growth is vertical, as in
    the brine figure, and an orthorhombic cell is a rectangle. Tick marks
    sit on the bottom axis only.
    """
    from PIL import Image, ImageDraw

    centre, half_h, half_v, forward, right, up = _camera(
        render, frame if frame is not None else positions, forward, up
    )
    raw = Image.open(path).convert("RGB")
    font = _legend_font(LEGEND_FONT_PX)
    swatch = 32
    row_h = 52
    left, right_margin, top = 150, 520, 8
    legend_rows = 1
    if legend:
        x = 16
        rows = 1
        for label, _colour in legend:
            width_px = swatch + 16 + float(font.getlength(label))
            if x + width_px > raw.width - 8:
                rows += 1
                x = 16
            x += width_px + 18
        legend_rows = rows
    def raw_pix(point):
        offset = np.asarray(point, dtype=float) - centre
        x = (float(np.dot(offset, right)) / (2.0 * half_h) + 0.5) * raw.width
        y = (0.5 - float(np.dot(offset, up)) / (2.0 * half_v)) * raw.height
        return x, y

    raw_face = [
        np.asarray(origin, dtype=float) + np.array((0.0, j * lengths[1], k * lengths[2]))
        for j, k in ((0, 0), (1, 0), (1, 1), (0, 1))
    ]
    raw_corners = [raw_pix(point) for point in raw_face]
    pad = int(np.ceil(float(pad_radius) * raw.width / (2.0 * half_h))) + 2
    crop_x0 = int(min(point[0] for point in raw_corners)) - pad
    crop_y0 = int(min(point[1] for point in raw_corners)) - pad
    crop_x1 = int(max(point[0] for point in raw_corners)) + pad + 1
    crop_y1 = int(max(point[1] for point in raw_corners)) + pad + 1
    if crop_x0 < 0 or crop_y0 < 0 or crop_x1 > raw.width or crop_y1 > raw.height:
        raise RuntimeError("the render does not include one sphere radius outside the cell")
    image = raw.crop((crop_x0, crop_y0, crop_x1, crop_y1))
    bottom = 340 + row_h * legend_rows
    canvas = Image.new("RGB", (image.width + left + right_margin, image.height + top + bottom), (255, 255, 255))
    canvas.paste(image, (left, top))
    draw = ImageDraw.Draw(canvas)

    def pix(point):
        x, y = raw_pix(point)
        return x - crop_x0 + left, y - crop_y0 + top

    face = [
        np.asarray(origin, dtype=float) + np.array((0.0, j * lengths[1], k * lengths[2]))
        for j, k in ((0, 0), (1, 0), (1, 1), (0, 1))
    ]
    corners = [pix(point) for point in face]
    for first, second in ((0, 1), (1, 2), (2, 3), (3, 0)):
        draw.line([corners[first], corners[second]], fill=(40, 40, 40), width=3)
    y_bottom = max(point[1] for point in corners)
    x_left = min(point[0] for point in corners)
    x_right = max(point[0] for point in corners)
    y_top = min(point[1] for point in corners)

    marks = []

    def add_text(name, xy, text):
        marks.append((name, text_bbox(font, xy, text)))
        draw.text(xy, text, font=font, fill=(20, 20, 20))

    for layer in sorted({int(item) for item in ids if item >= 0}):
        members = np.flatnonzero(ids == layer)
        horiz = float(np.median((positions[members] - centre) @ right))
        x = left + (horiz / (2.0 * half_h) + 0.5) * raw.width - crop_x0
        draw.line([(x, y_bottom), (x, y_bottom + 16)], fill=(20, 20, 20), width=3)
        text = str(layer)
        text_w = float(font.getlength(text))
        add_text(f"layer {text}", (x - text_w / 2.0, y_bottom + 20), text)
    title = "layer index"
    title_w = float(font.getlength(title))
    title_y = y_bottom + 28 + LEGEND_FONT_PX
    add_text("layer index", ((x_left + x_right - title_w) / 2.0, title_y), title)

    z_ice = positions[np.asarray(labels) == 1, 2]
    z_cut = float(np.percentile(z_ice, 98)) if z_ice.size else float(origin[2] + 0.5 * lengths[2])
    z_lo = float(origin[2])
    z_hi = z_lo + float(lengths[2])

    def y_of(z_value):
        point = np.asarray(origin, dtype=float) + np.array((0.0, 0.5 * lengths[1], z_value - z_lo))
        return pix(point)[1]

    for text, z_mid in (("liquid", 0.5 * (z_cut + z_hi)),):
        text_w = float(font.getlength(text))
        add_text(text, (max(8, x_left - text_w - 16), y_of(z_mid) - LEGEND_FONT_PX / 2.0), text)

    # One bracket over the crystal. The arrows name the two ends of that slab.
    bracket_x = x_right + 28
    y_front = y_of(z_cut)
    y_periodic = y_of(z_lo)
    draw.line([(bracket_x, y_front), (bracket_x, y_periodic)], fill=(20, 20, 20), width=3)
    slab = "ice Ih slab"
    add_text(slab, (bracket_x + 14, 0.5 * (y_front + y_periodic) - LEGEND_FONT_PX / 2.0), slab)
    for y_mark, text in ((y_front, "growth front"), (y_periodic, "periodic boundary")):
        draw.polygon([(x_right + 4, y_mark), (x_right + 22, y_mark - 7), (x_right + 22, y_mark + 7)], fill=(20, 20, 20))
        draw.line([(bracket_x - 8, y_mark), (bracket_x + 8, y_mark)], fill=(20, 20, 20), width=3)
        add_text(text, (bracket_x + 14, y_mark - LEGEND_FONT_PX - 4), text)

    legend_bottom = title_y + LEGEND_FONT_PX
    if legend:
        x, y = left + 16, title_y + LEGEND_FONT_PX + 16
        for label, colour in legend:
            shown = display_label(label)
            text_w = float(font.getlength(shown))
            if x + swatch + 16 + text_w > canvas.width - 8:
                x = left + 16
                y += row_h
            rgb = tuple(int(round(channel * 255)) for channel in colour)
            draw.rectangle([x, y, x + swatch, y + swatch], fill=rgb, outline=(30, 30, 30))
            marks.append((f"swatch {shown}", (x, y, x + swatch, y + swatch)))
            add_text(f"legend {shown}", (x + swatch + 10, y), shown)
            legend_bottom = max(legend_bottom, y + swatch, marks[-1][1][3])
            x += swatch + 16 + text_w + 18
    bar = 10.0 * raw.width / (2.0 * half_h)
    bar_y = legend_bottom + 28
    bar_x = x_left
    draw.line([(bar_x, bar_y), (bar_x + bar, bar_y)], fill=(20, 20, 20), width=6)
    marks.append(("1 nm bar", (bar_x, bar_y - 4, bar_x + bar, bar_y + 4)))
    add_text("1 nm", (bar_x + bar + 14, bar_y - LEGEND_FONT_PX / 2.0), "1 nm")
    require_disjoint(marks)
    canvas.save(path)


def render_interface(render, positions, labels, origin, lengths, path, width=None):
    """Growth vertical, orthographic view straight down x.

    Liquid uses the brine veil: radius 0.50 and 18 percent opacity.
    Bilayer ticks are short marks on the bottom axis.
    """
    positions = wrap_into_cell(positions, origin, lengths)
    shown = np.array(labels, copy=True)
    ice_z = positions[shown == 1, 2]
    front = float(np.percentile(ice_z, 98)) if ice_z.size else float(origin[2])
    # CHILL+ is left as the engine computed it. Only the picture fades
    # classes more than 1.5 nm from the growth front into liquid grey.
    far = np.abs(positions[:, 2] - front) > 15.0
    shown[far & (shown != 1)] = 2
    colours = np.array([render.CHILL.get(int(code), ("", (0.5, 0.5, 0.5)))[1] for code in shown])
    radii = np.where(shown == 2, BRINE_RADIUS["water"], 0.85)
    liquid = shown == 2
    basal = np.asarray(render.stacking_axis(positions, labels, np.isin(labels, [1])), dtype=float)
    # Look along -x so bilayer 0 is on the left. up stays +z, so the
    # camera (right, up, forward) stays right-handed.
    forward = np.array([-1.0, 0.0, 0.0])
    up = np.array([0.0, 0.0, 1.0])
    pad_radius = 0.85
    frame = padded_frame(positions, origin, lengths, forward, up, pad_radius)
    _centre, half_h, half_v, forward, _right, up = _camera(render, frame, forward, up)
    if width is None:
        short = 1700
        width = short if half_h <= half_v else int(np.ceil(short * half_h / max(half_v, 1e-6)))
    ice_png = path.with_name(path.stem + "-ice.png")
    liq_png = path.with_name(path.stem + "-liq.png")
    tachyon_scene(render, ice_png, positions[~liquid], colours[~liquid], radii[~liquid], [], forward, up, width, frame=frame)
    tachyon_scene(render, liq_png, positions[liquid], colours[liquid], radii[liquid], [], forward, up, width, frame=frame)
    from PIL import Image

    ice = np.asarray(Image.open(ice_png).convert("RGB"), dtype=np.float32)
    liquid_image = np.asarray(Image.open(liq_png).convert("RGB"), dtype=np.float32)
    veil = np.any(liquid_image < 250, axis=2)
    blended = ice.copy()
    blended[veil] = 0.18 * liquid_image[veil] + 0.82 * ice[veil]
    Image.fromarray(np.clip(blended, 0, 255).astype(np.uint8)).save(path)
    ids, _centres = basal_bilayers(positions, labels, basal)
    names = {0: "cubic", 1: "hexagonal", 2: "liquid", 3: "interfacial", 4: "clathrate", 5: "interf. clathrate"}
    present = {names[int(code)] for code in shown if int(code) in names}
    legend = [item for item in BRINE_LEGEND if item[0] in present]
    _present_interface(
        render, path, positions, labels, origin, lengths, forward, up, ids, legend,
        frame=frame, pad_radius=pad_radius,
    )


def _positions(path):
    rows = []
    reading = False
    with open(path) as handle:
        for line in handle:
            if line.startswith("ITEM: ATOMS"):
                reading = True
            elif reading:
                parts = line.split()
                if len(parts) < 5:
                    break
                rows.append([float(parts[2]), float(parts[3]), float(parts[4])])
    return np.asarray(rows, dtype=float)
