"""Tachyon scenes for the brine front and the hydrate cages.

Waters use the CHILL+ colours from scripts/science_render.py. Sodium
and chloride are extra spheres. Cage faces are polyhedra in those same
colours: 5^12 is hexagonal orange, 5^12 6^4 is cubic blue, the sH
medium cage is interfacial yellow, and 5^12 6^8 is clathrate green.
"""

from __future__ import annotations

import ctypes
from collections import defaultdict

import numpy as np

# Face-size census of each cage, and the CHILL+ colour index it borrows.
CAGE_KINDS = (
    ("512", {5: 12}, 1),
    ("51264", {5: 12, 6: 4}, 0),
    ("sh", {4: 3, 5: 6, 6: 3}, 3),
    ("51268", {5: 12, 6: 8}, 4),
)
NA = np.array((0.55, 0.17, 0.89))
CL = np.array((0.15, 0.78, 0.82))


def _rings(lib, path: str):
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
            rows = [[float(item) for item in handle.readline().split()] for _ in range(3)]
            xy = rows[0][2] if len(rows[0]) > 2 else 0.0
            xz = rows[1][2] if len(rows[1]) > 2 else 0.0
            yz = rows[2][2] if len(rows[2]) > 2 else 0.0
            return [
                np.array([rows[0][1] - rows[0][0], 0.0, 0.0]),
                np.array([xy, rows[1][1] - rows[1][0], 0.0]),
                np.array([xz, yz, rows[2][1] - rows[2][0]]),
            ]
    raise RuntimeError(f"no box in {path}")


def _unwrap_cage(faces, positions, vectors):
    """Bring every atom of one cage next to its neighbours."""
    neighbours = defaultdict(set)
    for face in faces:
        for index, atom in enumerate(face):
            other = face[(index + 1) % len(face)]
            neighbours[atom].add(other)
            neighbours[other].add(atom)
    shifts = [
        i * vectors[0] + j * vectors[1] + k * vectors[2]
        for i in (-1, 0, 1)
        for j in (-1, 0, 1)
        for k in (-1, 0, 1)
    ]
    start = faces[0][0]
    placed = {start: positions[start].copy()}
    queue = [start]
    while queue:
        atom = queue.pop()
        for other in neighbours[atom]:
            if other in placed:
                continue
            point = positions[other]
            placed[other] = point + min(shifts, key=lambda shift: np.linalg.norm(point + shift - placed[atom]))
            queue.append(other)
    return placed


def cage_meshes(cycles, groups, positions, vectors):
    """One vertex array and face list per colour index."""
    by_colour = defaultdict(list)
    reference = None
    shifts = [
        i * vectors[0] + j * vectors[1] + k * vectors[2]
        for i in (-1, 0, 1)
        for j in (-1, 0, 1)
        for k in (-1, 0, 1)
    ]
    for faces, colour in groups:
        placed = _unwrap_cage([cycles[face] for face in faces], positions, vectors)
        centre = np.mean(list(placed.values()), axis=0)
        if reference is None:
            reference = centre
        else:
            delta = min(shifts, key=lambda shift: np.linalg.norm(centre + shift - reference))
            placed = {atom: point + delta for atom, point in placed.items()}
            centre = centre + delta
        span = np.ptp(np.vstack(list(placed.values())), axis=0)
        if float(span.max()) > 14.0:
            continue
        oriented = []
        for face in faces:
            cycle = cycles[face]
            coords = np.vstack([placed[atom] for atom in cycle])
            normal = np.cross(coords[1] - coords[0], coords[2] - coords[0])
            if float(np.dot(normal, coords.mean(axis=0) - centre)) < 0.0:
                cycle = tuple(reversed(cycle))
            oriented.append(tuple(placed[atom] for atom in cycle))
        by_colour[colour].append(oriented)
    meshes = []
    for colour, cages in by_colour.items():
        vertices = []
        faces = []
        for cage in cages:
            for polygon in cage:
                base = len(vertices)
                vertices.extend(polygon)
                for index in range(1, len(polygon) - 1):
                    faces.append([base, base + index, base + index + 1])
        meshes.append((np.asarray(vertices, dtype=float), faces, colour))
    return meshes


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
    tinted, _near = render.depth_cue(colours, positions @ forward, 0.22)
    particles = Particles()
    particles.create_property("Position", data=np.ascontiguousarray(positions, dtype=np.float64))
    particles.create_property("Color", data=np.clip(tinted, 0.0, 1.0))
    particles.create_property("Radius", data=np.asarray(radii, dtype=np.float64))
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


def render_brine(render, rows, path, width=560):
    """Oxygens by CHILL+, sodium and chloride as larger spheres."""
    names = {name: index for index, name in enumerate((
        "cubic", "hexagonal", "water", "interfacial", "clathrate",
        "interClathrate", "unclassified", "reCubic", "reHex",
    ))}
    kept = [row for row in rows if row[0] in ("O", "Na", "Cl")]
    positions = np.array([[row[1], row[2], row[3]] for row in kept], dtype=float)
    colours = []
    radii = []
    for species, _x, _y, _z, label in kept:
        if species == "Na":
            colours.append(NA)
            radii.append(1.25)
        elif species == "Cl":
            colours.append(CL)
            radii.append(1.45)
        else:
            colours.append(np.array(render.CHILL.get(names[label], ("", (0.5, 0.5, 0.5)))[1]))
            radii.append(0.72 if label != "water" else 0.62)
    colours = np.vstack(colours)
    radii = np.asarray(radii)
    # Liquid is painted faintly, the same 0.40 blend as the interface.
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
    blended[veil] = 0.40 * liquid_image[veil] + 0.60 * ice[veil]
    Image.fromarray(np.clip(blended, 0, 255).astype(np.uint8)).save(path)


def render_host(render, lib, dump, path, kinds, width=560):
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
    meshes = cage_meshes(cycles, groups, positions, vectors)
    corners = np.vstack([vertices for vertices, _faces, _colour in meshes])
    colours = np.tile(np.array((0.96, 0.96, 0.96)), (len(corners), 1))
    radii = np.full(len(corners), 0.22)
    forward = np.array([0.55, 0.35, 0.76])
    up = np.array([0.0, 0.0, 1.0])
    tachyon_scene(render, path, corners, colours, radii, meshes, forward, up, width)
    return expected


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
