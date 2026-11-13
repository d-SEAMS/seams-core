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
# Slightly under the paper radii so the ice and liquid stay readable.
# Liquid is smaller than the ice oxygen. Ions are a step larger than the ice.
BRINE_RADIUS = {"O": 0.90, "water": 0.50, "Na": 1.32, "Cl": 1.45}
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
    positions = np.array([[row[1], row[2], row[3]] for row in kept], dtype=float)
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
    tachyon_scene(render, ice_png, positions[solid], colours[solid], radii[solid], [], forward, up, width, frame=positions)
    tachyon_scene(render, liq_png, positions[liquid], colours[liquid], radii[liquid], [], forward, up, width, frame=positions)
    from PIL import Image

    ice = np.asarray(Image.open(ice_png).convert("RGB"), dtype=np.float32)
    liquid_image = np.asarray(Image.open(liq_png).convert("RGB"), dtype=np.float32)
    veil = np.any(liquid_image < 250, axis=2)
    blended = ice.copy()
    blended[veil] = 0.18 * liquid_image[veil] + 0.82 * ice[veil]
    Image.fromarray(np.clip(blended, 0, 255).astype(np.uint8)).save(path)
    origin = np.array(box[:3], dtype=float)
    lengths = np.array(box[3:6], dtype=float)
    corners = [
        origin + np.array((i * lengths[0], j * lengths[1], k * lengths[2]))
        for i in (0, 1)
        for j in (0, 1)
        for k in (0, 1)
    ]
    legend = [item for item in BRINE_LEGEND if item[0] in present]
    annotate(render, path, positions, forward, up, width, corners, legend)


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


def annotate(render, path, positions, forward, up, width, corners, legend=None, layer_title=None):
    """Box outline, a 1 nm bar, and the colour legend under the atoms."""
    from PIL import Image, ImageDraw

    centre, half_h, half_v, forward, right, up = _camera(render, positions, forward, up)
    image = Image.open(path).convert("RGB")
    font = _legend_font(LEGEND_FONT_PX)
    swatch = 32
    row_h = 52
    legend_rows = 1
    if legend:
        x = 16
        rows = 1
        for label, _colour in legend:
            width_px = swatch + 16 + float(font.getlength(label))
            if x + width_px > image.width - 8:
                rows += 1
                x = 16
            x += width_px + 18
        legend_rows = rows
    margin = 24 + row_h * legend_rows + 64
    left = 110 if layer_title else 0
    canvas = Image.new("RGB", (image.width + left, image.height + margin), (255, 255, 255))
    canvas.paste(image, (left, 0))
    draw = ImageDraw.Draw(canvas)
    if layer_title:
        text = Image.new("RGB", (int(font.getlength(layer_title)) + 8, 52), (255, 255, 255))
        ImageDraw.Draw(text).text((0, 4), layer_title, font=font, fill=(20, 20, 20))
        turned = text.rotate(90, expand=True)
        canvas.paste(turned, (6, max(8, (image.height - turned.height) // 2)))

    def project(point):
        offset = np.asarray(point, dtype=float) - centre
        x = left + (float(np.dot(offset, right)) / (2.0 * half_h) + 0.5) * image.width
        y = (0.5 - float(np.dot(offset, up)) / (2.0 * half_v)) * image.height
        return (x, y)

    for first in range(8):
        for second in range(first + 1, 8):
            if bin(first ^ second).count("1") != 1:
                continue
            draw.line([project(corners[first]), project(corners[second])], fill=(70, 70, 70), width=2)
    if legend:
        x, y = left + 16, image.height + 16
        for label, colour in legend:
            text_w = float(font.getlength(label))
            if x + swatch + 16 + text_w > canvas.width - 8:
                x = left + 16
                y += row_h
            rgb = tuple(int(round(channel * 255)) for channel in colour)
            draw.rectangle([x, y, x + swatch, y + swatch], fill=rgb, outline=(30, 30, 30))
            draw.text((x + swatch + 10, y), label, font=font, fill=(20, 20, 20))
            x += swatch + 16 + text_w + 18
        bar_y = y + row_h + 8
    else:
        bar_y = image.height + 36
    bar = 10.0 * image.width / (2.0 * half_h)
    draw.line([(left + 16, bar_y), (left + 16 + bar, bar_y)], fill=(20, 20, 20), width=6)
    draw.text((left + 28 + bar, bar_y - 22), "1 nm", font=font, fill=(20, 20, 20))
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


def _present_interface(render, path, positions, labels, origin, lengths, forward, up, ids, legend):
    """Rectangle cell, bottom bilayer ticks, and labels outside the box.

    The camera looks down x with up along z, so growth is vertical, as in
    the brine figure, and an orthorhombic cell is a rectangle. Tick marks
    sit on the bottom axis only.
    """
    from PIL import Image, ImageDraw

    centre, half_h, half_v, forward, right, up = _camera(render, positions, forward, up)
    image = Image.open(path).convert("RGB")
    font = _legend_font(LEGEND_FONT_PX)
    swatch = 32
    row_h = 52
    left, right_margin, top = 170, 210, 36
    legend_rows = 1
    if legend:
        x = 16
        rows = 1
        for label, _colour in legend:
            width_px = swatch + 16 + float(font.getlength(label))
            if x + width_px > image.width - 8:
                rows += 1
                x = 16
            x += width_px + 18
        legend_rows = rows
    bottom = 36 + LEGEND_FONT_PX + 8 + LEGEND_FONT_PX + 16 + row_h * legend_rows + 110
    canvas = Image.new("RGB", (image.width + left + right_margin, image.height + top + bottom), (255, 255, 255))
    canvas.paste(image, (left, top))
    draw = ImageDraw.Draw(canvas)

    def pix(point):
        offset = np.asarray(point, dtype=float) - centre
        x = left + (float(np.dot(offset, right)) / (2.0 * half_h) + 0.5) * image.width
        y = top + (0.5 - float(np.dot(offset, up)) / (2.0 * half_v)) * image.height
        return x, y

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

    for layer in sorted({int(item) for item in ids if item >= 0}):
        members = np.flatnonzero(ids == layer)
        horiz = float(np.median((positions[members] - centre) @ right))
        x = left + (horiz / (2.0 * half_h) + 0.5) * image.width
        draw.line([(x, y_bottom), (x, y_bottom + 16)], fill=(20, 20, 20), width=3)
        text = str(layer)
        text_w = float(font.getlength(text))
        draw.text((x - text_w / 2.0, y_bottom + 20), text, font=font, fill=(20, 20, 20))
    title = "layer index"
    title_w = float(font.getlength(title))
    draw.text(((x_left + x_right - title_w) / 2.0, y_bottom + 20 + LEGEND_FONT_PX + 4), title, font=font, fill=(20, 20, 20))

    z_ice = positions[np.asarray(labels) == 1, 2]
    z_cut = float(np.percentile(z_ice, 98)) if z_ice.size else float(origin[2] + 0.5 * lengths[2])
    z_lo = float(origin[2])
    z_hi = z_lo + float(lengths[2])

    def y_of(z_value):
        point = np.asarray(origin, dtype=float) + np.array((0.0, 0.5 * lengths[1], z_value - z_lo))
        return pix(point)[1]

    for text, z_mid in (("ice Ih", 0.5 * (z_lo + z_cut)), ("liquid", 0.5 * (z_cut + z_hi))):
        text_w = float(font.getlength(text))
        draw.text((max(8, x_left - text_w - 16), y_of(z_mid) - LEGEND_FONT_PX / 2.0), text, font=font, fill=(20, 20, 20))

    # Both interfaces, outside the cell: the growth front and its periodic image.
    bracket_x = x_right + 28
    y_front = y_of(z_cut)
    y_periodic = y_of(z_lo)
    draw.line([(bracket_x, y_front), (bracket_x, y_periodic)], fill=(20, 20, 20), width=3)
    for y_mark, text in ((y_front, "interface"), (y_periodic, "periodic")):
        draw.polygon([(x_right + 4, y_mark), (x_right + 22, y_mark - 7), (x_right + 22, y_mark + 7)], fill=(20, 20, 20))
        draw.line([(bracket_x - 8, y_mark), (bracket_x + 8, y_mark)], fill=(20, 20, 20), width=3)
        draw.text((bracket_x + 14, y_mark - LEGEND_FONT_PX / 2.0), text, font=font, fill=(20, 20, 20))

    if legend:
        x, y = left + 16, image.height + top + 20 + 2 * LEGEND_FONT_PX + 24
        for label, colour in legend:
            text_w = float(font.getlength(label))
            if x + swatch + 16 + text_w > canvas.width - 8:
                x = left + 16
                y += row_h
            rgb = tuple(int(round(channel * 255)) for channel in colour)
            draw.rectangle([x, y, x + swatch, y + swatch], fill=rgb, outline=(30, 30, 30))
            draw.text((x + swatch + 10, y), label, font=font, fill=(20, 20, 20))
            x += swatch + 16 + text_w + 18
        bar_y = y + row_h + 16
    else:
        bar_y = image.height + top + 36
    bar = 10.0 * image.width / (2.0 * half_h)
    draw.line([(left + 16, bar_y), (left + 16 + bar, bar_y)], fill=(20, 20, 20), width=6)
    draw.text((left + 16, bar_y + 10), "1 nm", font=font, fill=(20, 20, 20))
    canvas.save(path)


def render_interface(render, positions, labels, origin, lengths, path, width=None):
    """Growth vertical, orthographic view straight down x.

    Liquid uses the brine veil: radius 0.50 and 18 percent opacity.
    Bilayer ticks are short marks on the bottom axis.
    """
    colours = np.array([render.CHILL.get(int(code), ("", (0.5, 0.5, 0.5)))[1] for code in labels])
    radii = np.where(np.asarray(labels) == 2, BRINE_RADIUS["water"], 0.85)
    liquid = np.asarray(labels) == 2
    basal = np.asarray(render.stacking_axis(positions, labels, np.isin(labels, [1])), dtype=float)
    forward = np.array([1.0, 0.0, 0.0])
    up = np.array([0.0, 0.0, 1.0])
    _centre, half_h, half_v, forward, _right, up = _camera(render, positions, forward, up)
    if width is None:
        short = 1700
        width = short if half_h <= half_v else int(np.ceil(short * half_h / max(half_v, 1e-6)))
    ice_png = path.with_name(path.stem + "-ice.png")
    liq_png = path.with_name(path.stem + "-liq.png")
    tachyon_scene(render, ice_png, positions[~liquid], colours[~liquid], radii[~liquid], [], forward, up, width, frame=positions)
    tachyon_scene(render, liq_png, positions[liquid], colours[liquid], radii[liquid], [], forward, up, width, frame=positions)
    from PIL import Image

    ice = np.asarray(Image.open(ice_png).convert("RGB"), dtype=np.float32)
    liquid_image = np.asarray(Image.open(liq_png).convert("RGB"), dtype=np.float32)
    veil = np.any(liquid_image < 250, axis=2)
    blended = ice.copy()
    blended[veil] = 0.18 * liquid_image[veil] + 0.82 * ice[veil]
    Image.fromarray(np.clip(blended, 0, 255).astype(np.uint8)).save(path)
    ids, _centres = basal_bilayers(positions, labels, basal)
    names = {0: "cubic", 1: "hexagonal", 2: "liquid", 3: "interfacial", 4: "clathrate", 5: "interf. clathrate"}
    present = {names[int(code)] for code in labels if int(code) in names}
    legend = [item for item in BRINE_LEGEND if item[0] in present]
    _present_interface(render, path, positions, labels, origin, lengths, forward, up, ids, legend)


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
