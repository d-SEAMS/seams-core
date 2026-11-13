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
from chill_palette import CHLORIDE, CLATHRATE, COLOUR, CUBIC, HEXAGONAL, INTERFACIAL, LEGEND, SODIUM
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
    origin, vectors = _box_cell(path)
    del origin
    return vectors


def _box_cell(path):
    with open(path) as handle:
        for line in handle:
            if not line.startswith("ITEM: BOX"):
                continue
            rows = [handle.readline().split() for _ in range(3)]
            origin, (lx, ly, lz), (xy, xz, yz) = cell_from_bounds(rows)
            return np.array(origin, dtype=float), [
                np.array([lx, 0.0, 0.0]),
                np.array([xy, ly, 0.0]),
                np.array([xz, yz, lz]),
            ]
    raise RuntimeError(f"no box in {path}")


# Cage colour index to the shared palette. The symbols match the paper.
CAGE_LEGEND = {
    0: (r"$5^{12}6^{4}$", np.array(CUBIC)),
    1: (r"$5^{12}$", np.array(HEXAGONAL)),
    3: (r"$4^{3}5^{6}6^{3}$", np.array(INTERFACIAL)),
    4: (r"$5^{12}6^{8}$", np.array(CLATHRATE)),
}


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
    raw = np.array([[row[1], row[2], row[3]] for row in kept], dtype=float)
    hex_z = np.array([raw[i, 2] for i, row in enumerate(kept) if row[0] == "O" and row[4] == "hexagonal"])
    shell_z = np.array([raw[i, 2] for i, row in enumerate(kept) if row[0] == "O" and row[4] == "interfacial"])
    origin = origin.copy()
    origin[2] = lower_ice_origin(hex_z, shell_z, float(lengths[2]))
    positions = wrap_into_cell(raw, origin, lengths)
    ice_z = np.array(
        [positions[i, 2] for i, row in enumerate(kept) if row[0] == "O" and row[4] == "hexagonal"]
    )
    front = connected_slab_top(ice_z)
    colours = []
    radii = []
    present = set()
    drawn_liquid = []
    for index, (species, _x, _y, _z, label) in enumerate(kept):
        if species == "Na":
            colours.append(NA)
            radii.append(BRINE_RADIUS["Na"])
            present.add("Na+")
            drawn_liquid.append(False)
        elif species == "Cl":
            colours.append(CL)
            radii.append(BRINE_RADIUS["Cl"])
            present.add("Cl-")
            drawn_liquid.append(False)
        else:
            # CHILL+ is unchanged. Only the picture fades classes more than
            # 1.5 nm from the growth front into liquid grey. Ions stay.
            shown = label
            if shown != "hexagonal" and abs(float(positions[index, 2]) - front) > 15.0:
                shown = "water"
            colours.append(BRINE_COLOUR.get(shown, BRINE_COLOUR["unclassified"]))
            radii.append(BRINE_RADIUS["O"] if shown != "water" else BRINE_RADIUS["water"])
            present.add("liquid" if shown == "water" else ("interf. clathrate" if shown == "interClathrate" else shown))
            drawn_liquid.append(shown == "water")
    colours = np.vstack(colours)
    radii = np.asarray(radii)
    liquid = np.asarray(drawn_liquid, dtype=bool)
    solid = ~liquid
    ice_png = path.with_name(path.stem + "-ice.png")
    liq_png = path.with_name(path.stem + "-liq.png")
    forward = np.array([0.0, 1.0, 0.0])
    up = np.array([0.0, 0.0, 1.0])
    pad_radius = float(np.max(radii))
    # The viewport is wider than the crop so a sphere centred on the cell
    # edge, including the top edge, is drawn whole.
    frame = padded_frame(positions, origin, lengths, forward, up, pad_radius + 1.0)
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
    annotate(
        render, path, positions, forward, up, width, corners, legend,
        frame=frame, origin=origin, lengths=lengths, ice_z=ice_z, pad_radius=pad_radius,
    )


def _cell_corners(origin, vectors):
    return [
        np.asarray(origin, dtype=float) + i * vectors[0] + j * vectors[1] + k * vectors[2]
        for i in (0, 1)
        for j in (0, 1)
        for k in (0, 1)
    ]


def present_cage(render, path, origin, vectors, legend, title, forward, up, frame, cage_points, pad_radius=2.0):
    """Legend and scale. The crop keeps every cage vertex, not only the cell."""
    from PIL import Image, ImageDraw

    corners = _cell_corners(origin, vectors)
    centre, half_h, half_v, forward, right, up = _camera(render, frame, forward, up)
    raw = Image.open(path).convert("RGB")
    font = _legend_font(LEGEND_FONT_PX)
    swatch = 32
    row_h = 52

    def raw_pix(point):
        offset = np.asarray(point, dtype=float) - centre
        x = (float(np.dot(offset, right)) / (2.0 * half_h) + 0.5) * raw.width
        y = (0.5 - float(np.dot(offset, up)) / (2.0 * half_v)) * raw.height
        return x, y

    pad = int(np.ceil(float(pad_radius) * raw.width / (2.0 * half_h))) + 4
    covered = list(corners) + [np.asarray(point, dtype=float) for point in cage_points]
    image, (crop_x0, crop_y0) = clip_to_padded_cell(raw, raw_pix, covered, pad)
    left, right_margin, top = 24, 24, 8
    margin = 24 + row_h * 2 + LEGEND_FONT_PX + 48
    canvas = Image.new("RGB", (image.width + left + right_margin, image.height + top + margin), (255, 255, 255))
    canvas.paste(image, (left, top))
    draw = ImageDraw.Draw(canvas)
    marks = []

    def add_text(name, xy, text):
        marks.append((name, text_bbox(font, xy, text)))
        draw.text(xy, text, font=font, fill=(20, 20, 20))

    def project(point):
        x, y = raw_pix(point)
        return x - crop_x0 + left, y - crop_y0 + top

    projected = [project(corner) for corner in corners]
    for first in range(8):
        for second in range(first + 1, 8):
            if bin(first ^ second).count("1") != 1:
                continue
            draw.line([projected[first], projected[second]], fill=(70, 70, 70), width=2)
    x_left = min(point[0] for point in projected)
    x_right = max(point[0] for point in projected)
    y_top = min(point[1] for point in projected)
    y_bottom = max(point[1] for point in projected)
    legend_bottom = image.height + top
    y = image.height + top + 8
    title_w = float(font.getlength(title))
    add_text(title, (left + 16, y), title)
    marks[-1] = (title, marks[-1][1])
    y += row_h
    x = left + 16
    for label, colour in legend:
        glyph = mathtext_image(label, LEGEND_FONT_PX)
        if x + swatch + 16 + glyph.width > canvas.width - 8:
            x = left + 16
            y += row_h
        rgb = tuple(int(round(channel * 255)) for channel in colour)
        draw.rectangle([x, y, x + swatch, y + swatch], fill=rgb, outline=(30, 30, 30))
        marks.append((f"swatch {label}", (x, y, x + swatch, y + swatch)))
        canvas.paste(glyph, (int(x + swatch + 10), int(y)), glyph)
        marks.append((f"legend {label}", (x + swatch + 10, y, x + swatch + 10 + glyph.width, y + glyph.height)))
        legend_bottom = max(legend_bottom, y + max(swatch, glyph.height))
        x += swatch + 16 + glyph.width + 18
    bar = 10.0 * raw.width / (2.0 * half_h)
    ten = project(centre + 10.0 * right)
    origin_px = project(centre)
    assert_nanometre_bar(bar, abs(ten[0] - origin_px[0]), 10.0)
    bar_y = legend_bottom + 28
    bar_x = left + 16
    draw.line([(bar_x, bar_y), (bar_x + bar, bar_y)], fill=(20, 20, 20), width=6)
    marks.append(("1 nm bar", (bar_x, bar_y - 4, bar_x + bar, bar_y + 4)))
    add_text("1 nm", (bar_x + bar + 14, bar_y - LEGEND_FONT_PX / 2.0), "1 nm")
    require_disjoint(marks)
    canvas.save(path)


def render_host(render, lib, dump, path, kinds, width=1700):
    cycles = _rings(lib, str(dump))
    positions = _positions(dump)
    groups = []
    expected = []
    legend = []
    seen = set()
    for _name, need, colour, limit in kinds:
        facesets = find_cages(cycles, need, limit)
        expected.append(len(facesets))
        for faces in facesets:
            groups.append((faces, colour))
        if colour not in seen and colour in CAGE_LEGEND:
            seen.add(colour)
            legend.append(CAGE_LEGEND[colour])
    origin, vectors = _box_cell(dump)
    meshes, stats = cage_meshes(cycles, groups, positions, vectors)
    assert_drawn_faces(stats)
    corners = np.vstack([vertices for vertices, _faces, _colour in meshes])
    colours = np.tile(np.array((0.96, 0.96, 0.96)), (len(corners), 1))
    radii = np.full(len(corners), 0.22)
    forward = np.array([0.55, 0.35, 0.76])
    up = np.array([0.0, 0.0, 1.0])
    cell = _cell_corners(origin, vectors)
    import science_render as science

    _forward, right, up_axis = science.basis(forward, up)
    extra = []
    for corner in cell:
        extra.append(corner + 3.0 * right)
        extra.append(corner - 3.0 * right)
        extra.append(corner + 3.0 * up_axis)
        extra.append(corner - 3.0 * up_axis)
    frame = np.vstack([corners, np.asarray(extra, dtype=float)])
    tachyon_scene(render, path, corners, colours, radii, meshes, forward, up, width, frame=frame)
    title = {"sii": "sII", "sh": "sH"}.get(path.stem, path.stem)
    present_cage(render, path, origin, vectors, legend, title, forward, up, frame, corners, pad_radius=2.0)
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


def legend_groups(labels, font, swatch, max_width):
    """One row when it fits, otherwise two rows of similar width.

    A greedy wrap leaves the last swatch alone. The split keeps both
    rows non-empty and as close in width as the order allows.
    """
    if not labels:
        return []
    widths = [swatch + 16.0 + float(font.getlength(label)) + 18.0 for label, _colour in labels]
    if sum(widths) <= max_width:
        return [list(labels)]
    best = None
    acc = 0.0
    total = sum(widths)
    for index in range(1, len(labels)):
        acc += widths[index - 1]
        other = total - acc
        if acc <= max_width and other <= max_width:
            score = abs(acc - other)
            if best is None or score < best[0]:
                best = (score, index)
    if best is None:
        groups = [[]]
        used = 0.0
        for item, width in zip(labels, widths):
            if groups[-1] and used + width > max_width:
                groups.append([])
                used = 0.0
            groups[-1].append(item)
            used += width
        return groups
    cut = best[1]
    return [list(labels[:cut]), list(labels[cut:])]


def assert_nanometre_bar(bar_px, span_px, span_angstroms, bar_angstroms=10.0):
    """The drawn 1 nm bar must match the on-screen length of the cell edge."""
    if span_px <= 1.0 or span_angstroms <= 0.0:
        raise RuntimeError("the scale bar has no cell edge to check")
    expected = bar_angstroms * float(span_px) / float(span_angstroms)
    if abs(float(bar_px) - expected) > 1.5:
        raise RuntimeError(
            f"1 nm bar is {bar_px:.2f} px but {span_angstroms:.3f} A of cell is {span_px:.2f} px"
        )


def mathtext_image(expr, height_px):
    """Rasterise one mathtext label at the legend height."""
    from io import BytesIO

    from matplotlib import mathtext
    from PIL import Image

    buf = BytesIO()
    mathtext.math_to_image(expr, buf, dpi=200, format="png")
    buf.seek(0)
    image = Image.open(buf).convert("RGBA")
    scale = float(height_px) / float(image.height)
    width = max(1, int(round(image.width * scale)))
    return image.resize((width, int(height_px)), Image.Resampling.LANCZOS)


def stroke_bracket(draw, x, y0, y1, tick=16):
    """A vertical span with a tick at each end."""
    top, bottom = (y0, y1) if y0 <= y1 else (y1, y0)
    if bottom - top < 4:
        return
    draw.line([(x, top), (x, bottom)], fill=(20, 20, 20), width=3)
    draw.line([(x - tick / 2.0, top), (x + tick / 2.0, top)], fill=(20, 20, 20), width=3)
    draw.line([(x - tick / 2.0, bottom), (x + tick / 2.0, bottom)], fill=(20, 20, 20), width=3)


def dashed_hline(draw, x0, x1, y, dash=18, gap=12, width=4):
    """A horizontal dashed line long enough to read at figure scale."""
    x = float(min(x0, x1))
    end = float(max(x0, x1))
    while x < end:
        draw.line([(x, y), (min(x + dash, end), y)], fill=(20, 20, 20), width=width)
        x += dash + gap


def clip_to_padded_cell(raw, project, corners, pad_px):
    """White out every pixel outside the cell expanded by ``pad_px`` on each side."""
    pts = [project(corner) for corner in corners]
    x0 = int(np.floor(min(point[0] for point in pts) - pad_px))
    y0 = int(np.floor(min(point[1] for point in pts) - pad_px))
    x1 = int(np.ceil(max(point[0] for point in pts) + pad_px))
    y1 = int(np.ceil(max(point[1] for point in pts) + pad_px))
    if x0 < 0 or y0 < 0 or x1 > raw.width or y1 > raw.height:
        raise RuntimeError("the render does not include one sphere radius outside the cell")
    arr = np.array(raw)
    white = np.array([255, 255, 255], dtype=arr.dtype)
    if y0 > 0:
        arr[:y0] = white
    if y1 < arr.shape[0]:
        arr[y1:] = white
    if x0 > 0:
        arr[:, :x0] = white
    if x1 < arr.shape[1]:
        arr[:, x1:] = white
    from PIL import Image

    return Image.fromarray(arr).crop((x0, y0, x1, y1)), (x0, y0)


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


def annotate(render, path, positions, forward, up, width, corners, legend=None, layer_title=None, frame=None, origin=None, lengths=None, ice_z=None, pad_radius=None):
    """Box outline, a 1 nm bar, and the colour legend under the atoms.

    The bar is below the legend. When ``ice_z`` is given, a bracket marks
    the ice, a depth axis runs up the cell, and an arrow marks the growth
    front clear of the boundary spheres.
    """
    from PIL import Image, ImageDraw

    centre, half_h, half_v, forward, right, up = _camera(render, frame if frame is not None else positions, forward, up)
    raw = Image.open(path).convert("RGB")
    font = _legend_font(LEGEND_FONT_PX)
    swatch = 32
    row_h = 52

    def raw_pix(point):
        offset = np.asarray(point, dtype=float) - centre
        x = (float(np.dot(offset, right)) / (2.0 * half_h) + 0.5) * raw.width
        y = (0.5 - float(np.dot(offset, up)) / (2.0 * half_v)) * raw.height
        return x, y

    crop_x0 = crop_y0 = 0
    sphere_pad = 0
    image = raw
    if pad_radius is not None and origin is not None and lengths is not None:
        sphere_pad = int(np.ceil(float(pad_radius) * raw.width / (2.0 * half_h))) + 2
        image, (crop_x0, crop_y0) = clip_to_padded_cell(
            raw, raw_pix, cell_corners(origin, lengths), sphere_pad
        )
    labels = [(display_label(label), colour) for label, colour in (legend or ())]
    left = 220 if ice_z is not None else (110 if layer_title else 0)
    right_margin = 720 if ice_z is not None else 0
    # The legend starts at left+16 and must end before the canvas edge.
    groups = legend_groups(labels, font, swatch, image.width + right_margin - 40)
    legend_rows = max(1, len(groups))
    margin = 24 + row_h * legend_rows + LEGEND_FONT_PX + 48
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
        x, y = raw_pix(point)
        return x - crop_x0 + left, y - crop_y0 + top

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
        # Depth, in nanometres, from the lower ice surface. The interface
        # figure numbers its basal layers; this slab is viewed along that
        # stack, so the visible axis is the growth direction.
        step = 10.0
        depth = 0
        z_tick = z_lo
        while z_tick <= z_hi + 1e-6:
            y_tick = y_of(z_tick)
            draw.line([(x_left - 14, y_tick), (x_left, y_tick)], fill=(20, 20, 20), width=3)
            label = str(depth)
            label_w = float(font.getlength(label))
            add_text(f"depth {label}", (x_left - 22 - label_w, y_tick - LEGEND_FONT_PX / 2.0), label)
            z_tick += step
            depth += 1
        depth_label = "depth from lower ice surface (nm)"
        plate = Image.new("RGB", (int(font.getlength(depth_label)) + 8, LEGEND_FONT_PX + 8), (255, 255, 255))
        ImageDraw.Draw(plate).text((0, 0), depth_label, font=font, fill=(20, 20, 20))
        turned = plate.rotate(90, expand=True)
        depth_x = 8
        depth_y = int(0.5 * (y_of(z_lo) + y_of(z_hi)) - turned.height / 2.0)
        depth_y = max(top, depth_y)
        canvas.paste(turned, (depth_x, depth_y))
        marks.append(("depth", (depth_x, depth_y, depth_x + turned.width, depth_y + turned.height)))
        bracket_x = x_right + max(sphere_pad, 28) + 18
        y_front = y_of(z_cut)
        y_ice = y_of(z_lo)
        y_top = y_of(z_hi)
        gap = 28
        # Screen y grows downward: the brine span is the upper one.
        stroke_bracket(draw, bracket_x, y_top, y_front - gap)
        stroke_bracket(draw, bracket_x, y_front + gap, y_ice)
        add_text("brine", (bracket_x + 14, 0.5 * (y_top + y_front) - LEGEND_FONT_PX / 2.0), "brine")
        add_text("ice Ih", (bracket_x + 14, 0.5 * (y_front + y_ice) - LEGEND_FONT_PX / 2.0), "ice Ih")
        dashed_hline(draw, x_left, x_right, y_front)
        add_text("growth front", (bracket_x + 200, y_front - LEGEND_FONT_PX / 2.0), "growth front")
    legend_bottom = image.height + top
    y = image.height + top + 16
    for group in groups:
        x = left + 16
        for label, colour in group:
            text_w = float(font.getlength(label))
            rgb = tuple(int(round(channel * 255)) for channel in colour)
            draw.rectangle([x, y, x + swatch, y + swatch], fill=rgb, outline=(30, 30, 30))
            marks.append((f"swatch {label}", (x, y, x + swatch, y + swatch)))
            add_text(f"legend {label}", (x + swatch + 10, y), label)
            legend_bottom = max(legend_bottom, y + swatch, marks[-1][1][3])
            x += swatch + 16 + text_w + 18
        y += row_h
    bar = 10.0 * raw.width / (2.0 * half_h)
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


def lower_ice_origin(z_hex, z_shell, length_z, gap=8.0):
    """Cell origin on the lowest hexagonal oxygen.

    The periodic image of the brine can sit just under that oxygen.
    Seating the origin on the ice, not on the brine, wraps that row to
    the top of the cell. ``z_shell`` is accepted so callers can pass the
    interfacial coordinates; the ice surface is hexagonal.
    """
    del z_shell
    z_hex = np.asarray(z_hex, dtype=float)
    order = np.sort(z_hex)
    if order.size == 0:
        raise ValueError("the ice slab is empty")
    diffs = np.diff(order)
    if diffs.size == 0 or float(diffs.max()) <= gap:
        return float(order[0])
    cut = float(order[int(np.argmax(diffs))])
    top = z_hex[z_hex > cut]
    return float(np.min(top) - length_z)


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
    """Rectangle cell and the same depth axis as the brine figure.

    The camera looks down x with up along z, so growth is vertical.
    ``ids`` is retained for the bilayer test and is not drawn.
    """
    del ids
    from PIL import Image, ImageDraw

    centre, half_h, half_v, forward, right, up = _camera(
        render, frame if frame is not None else positions, forward, up
    )
    raw = Image.open(path).convert("RGB")
    font = _legend_font(LEGEND_FONT_PX)
    swatch = 32
    row_h = 52
    left, right_margin, top = 220, 720, 8
    shown_legend = [(display_label(label), colour) for label, colour in (legend or ())]
    groups = [shown_legend] if shown_legend else []
    legend_rows = 1
    def raw_pix(point):
        offset = np.asarray(point, dtype=float) - centre
        x = (float(np.dot(offset, right)) / (2.0 * half_h) + 0.5) * raw.width
        y = (0.5 - float(np.dot(offset, up)) / (2.0 * half_v)) * raw.height
        return x, y

    raw_face = [
        np.asarray(origin, dtype=float) + np.array((0.0, j * lengths[1], k * lengths[2]))
        for j, k in ((0, 0), (1, 0), (1, 1), (0, 1))
    ]
    pad = int(np.ceil(float(pad_radius) * raw.width / (2.0 * half_h))) + 2
    image, (crop_x0, crop_y0) = clip_to_padded_cell(raw, raw_pix, raw_face, pad)
    groups = [shown_legend] if shown_legend else []
    legend_rows = 1
    legend_width = 16.0
    for label, _colour in shown_legend:
        legend_width += swatch + 16.0 + float(font.getlength(label)) + 18.0
    right_margin = max(right_margin, int(legend_width - image.width) + 48)
    bottom = 120 + row_h * legend_rows
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

    z_ice = positions[np.asarray(labels) == 1, 2]
    z_cut = float(np.percentile(z_ice, 98)) if z_ice.size else float(origin[2] + 0.5 * lengths[2])
    z_lo = float(origin[2])
    z_hi = z_lo + float(lengths[2])

    def y_of(z_value):
        point = np.asarray(origin, dtype=float) + np.array((0.0, 0.5 * lengths[1], z_value - z_lo))
        return pix(point)[1]

    step = 10.0
    depth = 0
    z_tick = z_lo
    while z_tick <= z_hi + 1e-6:
        y_tick = y_of(z_tick)
        draw.line([(x_left - 14, y_tick), (x_left, y_tick)], fill=(20, 20, 20), width=3)
        numeral = str(depth)
        numeral_w = float(font.getlength(numeral))
        add_text(f"depth {numeral}", (x_left - 22 - numeral_w, y_tick - LEGEND_FONT_PX / 2.0), numeral)
        z_tick += step
        depth += 1
    depth_label = "depth from lower ice surface (nm)"
    plate = Image.new("RGB", (int(font.getlength(depth_label)) + 8, LEGEND_FONT_PX + 8), (255, 255, 255))
    ImageDraw.Draw(plate).text((0, 0), depth_label, font=font, fill=(20, 20, 20))
    turned = plate.rotate(90, expand=True)
    depth_x = 8
    depth_y = max(top, int(0.5 * (y_of(z_lo) + y_of(z_hi)) - turned.height / 2.0))
    canvas.paste(turned, (depth_x, depth_y))
    marks.append(("depth", (depth_x, depth_y, depth_x + turned.width, depth_y + turned.height)))

    bracket_x = x_right + pad + 18
    y_front = y_of(z_cut)
    y_ice = y_of(z_lo)
    y_top_cell = y_of(z_hi)
    gap = 28
    stroke_bracket(draw, bracket_x, y_top_cell, y_front - gap)
    stroke_bracket(draw, bracket_x, y_front + gap, y_ice)
    add_text("liquid", (bracket_x + 14, 0.5 * (y_top_cell + y_front) - LEGEND_FONT_PX / 2.0), "liquid")
    add_text("ice Ih", (bracket_x + 14, 0.5 * (y_front + y_ice) - LEGEND_FONT_PX / 2.0), "ice Ih")
    dashed_hline(draw, x_left, x_right, y_front)
    add_text("growth front", (bracket_x + 200, y_front - LEGEND_FONT_PX / 2.0), "growth front")

    legend_bottom = image.height + top
    y = image.height + top + 16
    for group in groups:
        x = left + 16
        for label, colour in group:
            text_w = float(font.getlength(label))
            rgb = tuple(int(round(channel * 255)) for channel in colour)
            draw.rectangle([x, y, x + swatch, y + swatch], fill=rgb, outline=(30, 30, 30))
            marks.append((f"swatch {label}", (x, y, x + swatch, y + swatch)))
            add_text(f"legend {label}", (x + swatch + 10, y), label)
            legend_bottom = max(legend_bottom, y + swatch, marks[-1][1][3])
            x += swatch + 16 + text_w + 18
        y += row_h
    bar = 10.0 * raw.width / (2.0 * half_h)
    # The horizontal edge is the cell's y length when the camera looks down x.
    assert_nanometre_bar(bar, x_right - x_left, float(lengths[1]))
    assert_nanometre_bar(bar, abs(y_bottom - y_top), float(lengths[2]))
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
    The depth axis matches the brine figure.
    """
    origin = np.array(origin, dtype=float).copy()
    lengths = np.array(lengths, dtype=float)
    hex_z = positions[np.asarray(labels) == 1, 2]
    origin[2] = lower_ice_origin(hex_z, [], float(lengths[2]))
    positions = wrap_into_cell(positions, origin, lengths)
    if np.any(positions[:, 2] < float(np.min(positions[np.asarray(labels) == 1, 2])) - 1e-6):
        raise RuntimeError("an interface centre sits below the lower ice surface")
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
    frame = padded_frame(positions, origin, lengths, forward, up, pad_radius + 1.0)
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
