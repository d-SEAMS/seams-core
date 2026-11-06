"""Draw the structure, classification, ring and interface figures.

Run from the repository root:

    python3 paper/tex/figures/render_finish.py

The ideal crystals and tubes use the geometry of worked_example.py.
The other frames are the committed trajectories named in the log.
The ring bars are the type 7 median and interquartile range of the
one-thread samples in figures/data/ring-head2head.json.
Both arms are one thread, and the file records ovito_threads as 1.
"""

import json
import math
import os
import subprocess
import sys
from collections import deque
from pathlib import Path

from PIL import Image, ImageDraw, ImageFont

ROOT = Path(".").resolve()
FIG = Path("paper/tex/figures")
DATA = FIG / "data"
LOG_PATH = FIG / "logs" / "render-finish.log"
CUTOFF_A = 3.5
BOND_A = 2.75
PRISM_HEIGHT_A = 2.845
TUBE_LAYERS = 12
TUBE_BOX_A = 25.0

RUHI = {
    "coral": (255, 101, 93),
    "sunshine": (241, 219, 75),
    "teal": (0, 77, 64),
    "sky": (30, 136, 229),
    "magenta": (216, 27, 96),
    "black": (0, 0, 0),
    "oxygen": (255, 13, 13),
    "white": (255, 255, 255),
    "bond": (40, 40, 40),
}
RING_COLOUR = {4: "sky", 5: "sunshine", 6: "teal"}


def run_text(args):
    return subprocess.check_output(args, text=True).strip()


def cpu_model():
    with open("/proc/cpuinfo", encoding="utf-8") as handle:
        for line in handle:
            if line.startswith("model name"):
                return " ".join(line.split(":", 1)[1].replace("(R)", "").split())
    return ""


def write_log(lines):
    LOG_PATH.parent.mkdir(parents=True, exist_ok=True)
    LOG_PATH.write_text("\n".join(lines) + "\n", encoding="utf-8")


def header(command):
    status = subprocess.check_output(["git", "status", "--short"], text=True).rstrip()
    lines = [
        run_text(["date", "-u"]),
        "uname_m " + run_text(["uname", "-m"]),
        "nproc " + run_text(["nproc"]),
        "cpu_model " + cpu_model(),
        "clone https://github.com/d-SEAMS/seams-core",
        "head " + run_text(["git", "rev-parse", "HEAD"]),
        "status_begin",
        status,
        "status_end",
        "pwd .",
        "TMPDIR .",
        "command " + command,
    ]
    return lines


def font(size):
    return ImageFont.truetype("/usr/share/fonts/truetype/dejavu/DejaVuSans.ttf", size)


def cubic_diamond(repetitions):
    lattice = 4.0 * BOND_A / math.sqrt(3.0)
    fcc = ((0.0, 0.0, 0.0), (0.5, 0.5, 0.0), (0.5, 0.0, 0.5), (0.0, 0.5, 0.5))
    basis = fcc + tuple((x + 0.25, y + 0.25, z + 0.25) for x, y, z in fcc)
    side = repetitions * lattice
    points = [
        (
            ((i + x) * lattice) % side,
            ((j + y) * lattice) % side,
            ((k + z) * lattice) % side,
        )
        for i in range(repetitions)
        for j in range(repetitions)
        for k in range(repetitions)
        for x, y, z in basis
    ]
    return points, [(0.0, side), (0.0, side), (0.0, side)]


def lonsdaleite(repetitions):
    lattice = BOND_A / ((3.0 / 8.0) * math.sqrt(8.0 / 3.0))
    height = lattice * math.sqrt(8.0 / 3.0)
    hexagonal = (
        (1.0 / 3.0, 2.0 / 3.0, 0.0),
        (2.0 / 3.0, 1.0 / 3.0, 0.5),
        (1.0 / 3.0, 2.0 / 3.0, 3.0 / 8.0),
        (2.0 / 3.0, 1.0 / 3.0, 7.0 / 8.0),
    )
    x_length = lattice
    y_length = lattice * math.sqrt(3.0)
    basis = []
    for hx, hy, hz in hexagonal:
        x = (hx - 0.5 * hy) * lattice
        y = hy * math.sqrt(3.0) * lattice / 2.0
        for shift_x, shift_y in ((0.0, 0.0), (x_length / 2.0, y_length / 2.0)):
            basis.append(
                (
                    (x + shift_x) % x_length,
                    (y + shift_y) % y_length,
                    hz * height,
                )
            )
    box = (repetitions * x_length, repetitions * y_length, repetitions * height)
    points = [
        (
            (i * x_length + x) % box[0],
            (j * y_length + y) % box[1],
            (k * height + z) % box[2],
        )
        for i in range(repetitions)
        for j in range(repetitions)
        for k in range(repetitions)
        for x, y, z in basis
    ]
    return points, [(0.0, box[0]), (0.0, box[1]), (0.0, box[2])]


def nanotube(sides):
    radius = BOND_A / (2.0 * math.sin(math.pi / sides))
    centre = TUBE_BOX_A / 2.0
    points = []
    for layer in range(TUBE_LAYERS):
        for vertex in range(sides):
            angle = 2.0 * math.pi * vertex / sides
            points.append(
                (
                    centre + radius * math.cos(angle),
                    centre + radius * math.sin(angle),
                    (layer + 0.5) * PRISM_HEIGHT_A,
                )
            )
    length = TUBE_LAYERS * PRISM_HEIGHT_A
    return points, [(0.0, TUBE_BOX_A), (0.0, TUBE_BOX_A), (0.0, length)]


def read_lammps(path, atom_type=None):
    lines = Path(path).read_text(encoding="utf-8", errors="replace").splitlines()
    index = 0
    while index < len(lines) and "NUMBER OF ATOMS" not in lines[index]:
        index += 1
    count = int(lines[index + 1])
    index += 2
    while not lines[index].startswith("ITEM: BOX"):
        index += 1
    index += 1
    bounds = []
    while not lines[index].startswith("ITEM:"):
        parts = lines[index].split()
        bounds.append((float(parts[0]), float(parts[1])))
        index += 1
    columns = lines[index].split()[2:]
    index += 1
    points = []
    for row in lines[index : index + count]:
        rec = dict(zip(columns, row.split()))
        if atom_type is not None and int(float(rec["type"])) != atom_type:
            continue
        points.append((float(rec["x"]), float(rec["y"]), float(rec["z"])))
    return points, bounds


def read_xyz(path):
    lines = Path(path).read_text(encoding="utf-8").splitlines()
    count = int(lines[0].split()[0])
    comment = lines[1]
    points = []
    for row in lines[2 : 2 + count]:
        parts = row.split()
        points.append((float(parts[-3]), float(parts[-2]), float(parts[-1])))
    box = None
    if "Lattice=" in comment:
        raw = comment.split("Lattice=", 1)[1].strip().strip('"')
        values = []
        for token in raw.replace('"', " ").split():
            try:
                values.append(float(token))
            except ValueError:
                break
            if len(values) == 9:
                break
        lengths = (values[0], values[4], values[8])
        box = [(0.0, lengths[0]), (0.0, lengths[1]), (0.0, lengths[2])]
    if box is None:
        span = []
        for axis in range(3):
            coords = [point[axis] for point in points]
            span.append((min(coords) - 5.0, max(coords) + 5.0))
        box = span
    return points, box


def mic(delta, length):
    return delta - round(delta / length) * length


def pair_vector(a, b, box):
    return tuple(mic(a[axis] - b[axis], box[axis][1] - box[axis][0]) for axis in range(3))


def neighbour_graph(points, box, cutoff):
    count = len(points)
    cell = cutoff
    origins = [box[axis][0] for axis in range(3)]
    lengths = [box[axis][1] - box[axis][0] for axis in range(3)]
    dims = [max(1, int(lengths[axis] / cell)) for axis in range(3)]
    buckets = {}
    keys = []
    for index, point in enumerate(points):
        key = tuple(
            min(dims[axis] - 1, max(0, int((point[axis] - origins[axis]) / cell)))
            for axis in range(3)
        )
        keys.append(key)
        buckets.setdefault(key, []).append(index)
    graph = [[] for _ in range(count)]
    edges = []
    cutoff2 = cutoff * cutoff
    for index, point in enumerate(points):
        key = keys[index]
        seen_cells = set()
        for dx in (-1, 0, 1):
            for dy in (-1, 0, 1):
                for dz in (-1, 0, 1):
                    other = (
                        (key[0] + dx) % dims[0],
                        (key[1] + dy) % dims[1],
                        (key[2] + dz) % dims[2],
                    )
                    if other in seen_cells:
                        continue
                    seen_cells.add(other)
                    for neighbour in buckets.get(other, ()):
                        if neighbour <= index:
                            continue
                        vector = pair_vector(point, points[neighbour], box)
                        dist2 = vector[0] ** 2 + vector[1] ** 2 + vector[2] ** 2
                        if dist2 <= cutoff2 and dist2 > 1e-8:
                            graph[index].append(neighbour)
                            graph[neighbour].append(index)
                            edges.append((index, neighbour))
    return graph, edges


def primitive_rings(graph, max_size):
    found = set()
    for start, _ in enumerate(graph):
        def walk(path):
            head = path[-1]
            if len(path) > max_size:
                return
            for nxt in graph[head]:
                if nxt == start and len(path) >= 3:
                    if path[1] < path[-1]:
                        found.add(tuple(path))
                    continue
                if nxt < start or nxt in path:
                    continue
                walk(path + [nxt])

        walk([start])
    rings = []
    for cycle in found:
        if is_primitive(cycle, graph):
            rings.append(cycle)
    return rings


def distances(graph, source, limit):
    dist = {source: 0}
    queue = deque([source])
    while queue:
        node = queue.popleft()
        if dist[node] == limit:
            continue
        for nxt in graph[node]:
            if nxt not in dist:
                dist[nxt] = dist[node] + 1
                queue.append(nxt)
    return dist


def is_primitive(cycle, graph):
    size = len(cycle)
    for left, start in enumerate(cycle):
        dist = distances(graph, start, size // 2)
        for step in range(1, size // 2 + 1):
            other = cycle[(left + step) % size]
            if dist.get(other, size) < step:
                return False
    return True


def bond_ring_size(edges, rings):
    smallest = {}
    for cycle in rings:
        size = len(cycle)
        members = set(cycle)
        for a, b in edges:
            if a in members and b in members:
                key = (a, b) if a < b else (b, a)
                # An edge of the ring, not merely two vertices that are bonded
                # off the ring. Check adjacency on the cycle.
                ia = cycle.index(a)
                ib = cycle.index(b)
                gap = abs(ia - ib)
                if min(gap, size - gap) != 1:
                    continue
                if key not in smallest or size < smallest[key]:
                    smallest[key] = size
    return smallest


def view_centre(points, box):
    centre = []
    for axis in range(3):
        lo, hi = box[axis]
        coords = [point[axis] for point in points]
        span = max(coords) - min(coords)
        length = hi - lo
        if span < 0.8 * length:
            centre.append(sum(coords) / len(coords))
        else:
            centre.append(0.5 * (lo + hi))
    return centre


def project_point(point, centre, box):
    azimuth = math.radians(32.0)
    elevation = math.radians(18.0)
    # The bond endpoint is the short image, so a second wrap would throw it
    # across the cell. Coordinates stay where the caller placed them.
    del box
    vector = tuple(point[axis] - centre[axis] for axis in range(3))
    x, y, z = vector
    xz = x * math.cos(azimuth) + z * math.sin(azimuth)
    zz = -x * math.sin(azimuth) + z * math.cos(azimuth)
    yy = y * math.cos(elevation) - zz * math.sin(elevation)
    depth = y * math.sin(elevation) + zz * math.cos(elevation)
    return xz, yy, depth


def draw_structure(draw, origin, size, points, box, edges, edge_colour, atom_colour, label, fonts):
    left, top = origin
    width, height = size
    centre = view_centre(points, box)
    projected = [project_point(point, centre, box) for point in points]
    xs = [item[0] for item in projected]
    ys = [item[1] for item in projected]
    span = max(max(xs) - min(xs), max(ys) - min(ys), 1.0)
    scale = 0.78 * min(width, height - 28) / span
    mid_x = 0.5 * (min(xs) + max(xs))
    mid_y = 0.5 * (min(ys) + max(ys))

    def pixel(item):
        return (
            left + width / 2.0 + (item[0] - mid_x) * scale,
            top + (height - 28) / 2.0 - (item[1] - mid_y) * scale,
        )

    order = sorted(range(len(points)), key=lambda index: projected[index][2])
    pixels = [pixel(item) for item in projected]
    radius = max(2.0, min(7.0, 0.22 * BOND_A * scale))
    for a, b in edges:
        colour = edge_colour(a, b)
        delta = pair_vector(points[b], points[a], box)
        end = tuple(points[a][axis] + delta[axis] for axis in range(3))
        end_px = pixel(project_point(end, centre, box))
        draw.line([pixels[a], end_px], fill=colour, width=2 if radius > 3.5 else 1)
    for index in order:
        x, y = pixels[index]
        fill = atom_colour(index)
        draw.ellipse(
            [x - radius, y - radius, x + radius, y + radius],
            fill=fill,
            outline=RUHI["black"],
        )
    draw.text((left + 8, top + height - 24), label, fill=RUHI["black"], font=fonts[0])


def panel_figure(panels, columns, path, panel_size=(460, 400)):
    rows = math.ceil(len(panels) / columns)
    width = columns * panel_size[0]
    height = rows * panel_size[1]
    image = Image.new("RGB", (width, height), RUHI["white"])
    draw = ImageDraw.Draw(image)
    fonts = (font(22),)
    for index, panel in enumerate(panels):
        column = index % columns
        row = index // columns
        origin = (column * panel_size[0], row * panel_size[1])
        draw_structure(draw, origin, panel_size, *panel, fonts)
    image.save(path, "PNG")
    return len(panels)


def slab(points, box, axis, thickness):
    coords = [point[axis] for point in points]
    mid = 0.5 * (min(coords) + max(coords))
    kept = [point for point in points if abs(point[axis] - mid) <= 0.5 * thickness]
    return kept, box


def constant_colour(name):
    colour = RUHI[name]

    def paint(*_args):
        return colour

    return paint


def systems(log):
    frames = []
    built = [
        ("(a) ice Ih", lonsdaleite(2)),
        ("(b) ice Ic", cubic_diamond(2)),
        ("(c) square tube", nanotube(4)),
        ("(d) pentagonal tube", nanotube(5)),
        ("(e) hexagonal tube", nanotube(6)),
    ]
    for label, (points, box) in built:
        frames.append((label, points, box, "built"))
    frames.append(("(f) structure I",) + read_lammps("input/traj/genice_sI.lammpstrj") + ("input/traj/genice_sI.lammpstrj",))
    frames.append(("(g) structure II",) + read_lammps("input/traj/genice_sII.lammpstrj") + ("input/traj/genice_sII.lammpstrj",))
    tip_points, tip_box = read_lammps("input/traj/exampleTraj.lammpstrj", atom_type=2)
    frames.append(
        (
            "(h) TIP4P window",
            *slab(tip_points, tip_box, 2, 36.0),
            "input/traj/exampleTraj.lammpstrj type 2, 36 A along z",
        )
    )
    nuc_points, nuc_box = read_lammps("input/traj/nucleation_f1.lammpstrj")
    frames.append(
        (
            "(i) nucleation slab",
            *slab(nuc_points, nuc_box, 2, 12.0),
            "input/traj/nucleation_f1.lammpstrj, 12 A slab",
        )
    )
    panels = []
    for label, points, box, source in frames:
        graph, edges = neighbour_graph(points, box, CUTOFF_A)
        degree = (2.0 * len(edges) / len(points)) if points else 0.0
        log.append(
            f"system {label} source {source} atoms {len(points)} bonds {len(edges)} "
            f"mean_degree {degree:.3f}"
        )
        panels.append(
            (
                points,
                box,
                edges,
                constant_colour("bond"),
                constant_colour("oxygen"),
                label,
            )
        )
    panel_figure(panels, 3, FIG / "systems.png", (520, 430))
    log.append("wrote paper/tex/figures/systems.png")


def classify(log):
    specs = [
        ("(a) hexagonal cage",) + read_xyz("templates/hc.xyz") + ("templates/hc.xyz",),
        ("(b) double-diamond cage",) + read_xyz("templates/ddc.xyz") + ("templates/ddc.xyz",),
        ("(c) hexagonal tube",) + nanotube(6) + ("built hexagonal tube",),
        ("(d) structure I",) + read_lammps("input/traj/genice_sI.lammpstrj") + ("input/traj/genice_sI.lammpstrj",),
    ]
    panels = []
    for label, points, box, source in specs:
        graph, edges = neighbour_graph(points, box, CUTOFF_A)
        rings = primitive_rings(graph, 6)
        sizes = bond_ring_size(edges, rings)
        counts = {}
        for size in sizes.values():
            counts[size] = counts.get(size, 0) + 1
        log.append(
            f"classify {label} source {source} atoms {len(points)} "
            f"rings {len(rings)} bond_sizes {counts}"
        )

        def edge_colour(a, b, sizes=sizes):
            key = (a, b) if a < b else (b, a)
            return RUHI[RING_COLOUR.get(sizes.get(key, 0), "black")]

        panels.append(
            (points, box, edges, edge_colour, constant_colour("oxygen"), label)
        )
    panel_figure(panels, 2, FIG / "classify.png", (640, 520))
    log.append("wrote paper/tex/figures/classify.png")


def percentile_type7(samples, percent):
    ordered = sorted(samples)
    count = len(ordered)
    if count == 0:
        raise ValueError("empty sample")
    if count == 1:
        return ordered[0]
    position = (count - 1) * percent / 100.0
    low = math.floor(position)
    high = math.ceil(position)
    weight = position - low
    return ordered[low] * (1.0 - weight) + ordered[high] * weight


def ring_record(log):
    record = json.loads((DATA / "ring-head2head.json").read_text(encoding="utf-8"))
    policy = record["thread_policy"]
    ovito_threads = int(policy["ovito_threads"])
    dseams_threads = [int(value) for value in policy["dseams_threads_reported"]]
    if ovito_threads != 1 or any(value != 1 for value in dseams_threads):
        raise SystemExit(
            f"ring record is not one thread: ovito {ovito_threads} dseams {dseams_threads}"
        )
    log.append(
        "ring_record machine uname_m {uname_m} nproc {nproc} cpu_model {cpu_model}".format(
            **record["machine"]
        )
    )
    log.append(
        "ring_record repeats {repeats} warmup_runs {warmup_runs}".format(**record["protocol"])
    )
    log.append(
        "ring_record ovito_threads "
        + str(ovito_threads)
        + " dseams_threads "
        + " ".join(str(value) for value in dseams_threads)
    )
    arms = ("dseams_omp1", "ovito")
    stats = {}
    plot_rows = []
    for name, frame in record["datasets"].items():
        stats[name] = {}
        for arm in arms:
            samples = frame["ring_stage_ms"][arm]["samples_ms"]
            if len(samples) < 6:
                raise SystemExit(f"{name} {arm} has fewer than 6 samples")
            median = percentile_type7(samples, 50)
            other = quantile_type7(samples, 0.5)
            if abs(median - other) > 1e-9:
                raise SystemExit(f"type 7 mismatch {name} {arm} {median} {other}")
            q1 = percentile_type7(samples, 25)
            q3 = percentile_type7(samples, 75)
            stats[name][arm] = (median, q1, q3, q3 - q1)
            log.append(
                f"ring {name} atoms {frame['n_atoms']} arm {arm} "
                f"n {len(samples)} median_ms {median:.6f} q1_ms {q1:.6f} "
                f"q3_ms {q3:.6f} iqr_ms {q3 - q1:.6f} "
                f"printed_median {median:.3f} samples {' '.join(f'{v:.6f}' for v in samples)}"
            )
        dseams_median = stats[name]["dseams_omp1"][0]
        ovito_median = stats[name]["ovito"][0]
        relation = "dseams_slower" if dseams_median > ovito_median else "dseams_faster"
        log.append(
            f"ring {name} atoms {frame['n_atoms']} relation {relation} "
            f"dseams_median_ms {dseams_median:.6f} ovito_median_ms {ovito_median:.6f}"
        )
        plot_rows.append(
            (int(frame["n_atoms"]), dseams_median, ovito_median)
        )
    log.append("type7_check two formulas agree pass")
    macros = ["% Generated by paper/tex/figures/render_finish.py.", "% snapper:off"]
    for name, arms_stats in stats.items():
        token = name.replace("-", "")
        for arm, (median, q1, q3, iqr) in arms_stats.items():
            key = f"rp@{token}@{arm}"
            macros.append(f"\\expandafter\\def\\csname dsv@{key}@median\\endcsname{{{median:.6f}}}")
            macros.append(f"\\expandafter\\def\\csname dsv@{key}@print\\endcsname{{{median:.3f}}}")
            macros.append(f"\\expandafter\\def\\csname dsv@{key}@q1\\endcsname{{{q1:.6f}}}")
            macros.append(f"\\expandafter\\def\\csname dsv@{key}@q3\\endcsname{{{q3:.6f}}}")
            macros.append(f"\\expandafter\\def\\csname dsv@{key}@iqr\\endcsname{{{iqr:.6f}}}")
    macros.append("% snapper:on")
    (DATA / "ringperf_macros.tex").write_text("\n".join(macros) + "\n", encoding="utf-8")
    plot_path = DATA / "plot-serial-rings.dat"
    plot_lines = ["sites dseams ovito"]
    for atoms, dseams_median, ovito_median in sorted(plot_rows):
        plot_lines.append(f"{atoms} {dseams_median:.6f} {ovito_median:.6f}")
    plot_path.write_text("\n".join(plot_lines) + "\n", encoding="utf-8")
    draw_ring_plot(record, stats)
    log.append("wrote paper/tex/figures/ringperf.png")
    log.append("wrote paper/tex/figures/data/ringperf_macros.tex")
    log.append("wrote paper/tex/figures/data/plot-serial-rings.dat")
    return stats


def quantile_type7(values, probability):
    ordered = sorted(float(value) for value in values)
    count = len(ordered)
    position = (count - 1) * probability + 1
    index = int(math.floor(position))
    weight = position - index
    left = ordered[index - 1]
    if index >= count:
        return left
    return left + weight * (ordered[index] - left)


def draw_ring_plot(record, stats):
    image = Image.new("RGB", (1280, 860), RUHI["white"])
    draw = ImageDraw.Draw(image)
    title = font(26)
    label = font(20)
    small = font(16)
    names = ["orth-4k", "tri-136", "orth-1e5", "tri-1e5"]
    arms = [
        ("dseams_omp1", "d-SEAMS, 1 thread", "black"),
        ("ovito", "OVITO, 1 thread", "coral"),
    ]
    left, right, top, bottom = 110, 980, 90, 760
    decades = (0, 5)  # 10^0 to 10^5 ms
    def y_of(value):
        logv = math.log10(max(value, 1.0))
        frac = (logv - decades[0]) / (decades[1] - decades[0])
        return bottom - frac * (bottom - top)

    for decade in range(decades[0], decades[1] + 1):
        y = y_of(10 ** decade)
        draw.line([(left, y), (right, y)], fill=(220, 220, 220), width=1)
        draw.text((20, y - 10), f"{10 ** decade:g}", fill=RUHI["black"], font=small)
    draw.line([(left, top), (left, bottom), (right, bottom)], fill=RUHI["black"], width=2)
    slot = (right - left) / len(names)
    bar_w = 28
    for index, name in enumerate(names):
        centre = left + (index + 0.5) * slot
        atoms = record["datasets"][name]["n_atoms"]
        for arm_index, (arm, _caption, colour) in enumerate(arms):
            median, q1, q3, _iqr = stats[name][arm]
            x = centre + (arm_index - 0.5) * (bar_w + 8)
            y_med = y_of(median)
            draw.rectangle([x - bar_w / 2, y_med, x + bar_w / 2, bottom], fill=RUHI[colour])
            draw.line([(x, y_of(q1)), (x, y_of(q3))], fill=RUHI["black"], width=3)
            draw.line([(x - 8, y_of(q1)), (x + 8, y_of(q1))], fill=RUHI["black"], width=2)
            draw.line([(x - 8, y_of(q3)), (x + 8, y_of(q3))], fill=RUHI["black"], width=2)
        draw.text((centre - 70, bottom + 16), f"{name}", fill=RUHI["black"], font=label)
        draw.text((centre - 70, bottom + 42), f"{atoms} atoms", fill=RUHI["black"], font=small)
    draw.text((left, 28), "Ring stage, median and type 7 interquartile range", fill=RUHI["black"], font=title)
    draw.text((24, 420), "ms", fill=RUHI["black"], font=label)
    legend_x = 1020
    for index, (_arm, caption, colour) in enumerate(arms):
        y = 120 + index * 48
        draw.rectangle([legend_x, y, legend_x + 28, y + 28], fill=RUHI[colour])
        draw.text((legend_x + 40, y + 2), caption, fill=RUHI["black"], font=small)
    image.save(FIG / "ringperf.png", "PNG")


def interface(log):
    image = Image.new("RGB", (1480, 1060), RUHI["white"])
    draw = ImageDraw.Draw(image)
    title = font(28)
    body = font(22)
    small = font(20)

    def font_size_of(face, text):
        bbox = draw.textbbox((0, 0), text, font=face)
        return bbox[2] - bbox[0], bbox[3] - bbox[1]

    def box(xy, text, fill, outline):
        draw.rounded_rectangle(xy, radius=12, fill=fill, outline=outline, width=3)
        x0, y0, x1, y1 = xy
        lines = text.split("\n")
        heights = [font_size_of(body, line)[1] for line in lines]
        block = sum(heights) + 6 * (len(lines) - 1)
        y = y0 + (y1 - y0 - block) / 2
        for line, th in zip(lines, heights):
            tw, _ = font_size_of(body, line)
            draw.text((x0 + (x1 - x0 - tw) / 2, y), line, fill=RUHI["black"], font=body)
            y += th + 6

    def arrow(start, end):
        draw.line([start, end], fill=RUHI["sky"], width=4)
        dx, dy = end[0] - start[0], end[1] - start[1]
        length = math.hypot(dx, dy) or 1.0
        ux, uy = dx / length, dy / length
        px, py = -uy, ux
        tip = end
        left = (tip[0] - 16 * ux + 8 * px, tip[1] - 16 * uy + 8 * py)
        right = (tip[0] - 16 * ux - 8 * px, tip[1] - 16 * uy - 8 * py)
        draw.polygon([tip, left, right], fill=RUHI["sky"])

    draw.text((48, 28), "C interface of d-SEAMS", fill=RUHI["black"], font=title)
    draw.text(
        (860, 36),
        "A LAMMPS compute calls this interface during a run.",
        fill=RUHI["black"],
        font=small,
    )
    box((500, 90, 980, 180), "LAMMPS dump, XYZ, CON", (232, 244, 252), RUHI["sky"])
    box((400, 250, 1080, 360), "seams, Python, R, Julia, Lua", (255, 232, 230), RUHI["coral"])
    box((500, 430, 980, 530), "C interface", (255, 255, 255), RUHI["teal"])
    analyses = [(48, "rings"), (340, "cages"), (632, "CHILL+"), (924, "prisms"), (1216, "clathrates")]
    for x, name in analyses:
        box((x, 620, x + 220, 720), name, (232, 244, 252), RUHI["sky"])
    box((360, 820, 1120, 920), "labels of each molecule and counts of each frame", (255, 255, 255), RUHI["magenta"])
    arrow((740, 180), (740, 250))
    arrow((740, 360), (740, 430))
    draw.line([(158, 575), (1326, 575)], fill=RUHI["sky"], width=4)
    arrow((740, 530), (740, 575))
    for x, _name in analyses:
        arrow((x + 110, 575), (x + 110, 620))
    arrow((740, 720), (740, 820))
    image.save(FIG / "interface.png", "PNG")
    log.append("wrote paper/tex/figures/interface.png")


def main():
    os.chdir(ROOT)
    os.environ["TMPDIR"] = "."
    command = "python3 paper/tex/figures/render_finish.py"
    lines = header(command)
    code = 0
    try:
        if not (DATA / "ring-head2head.json").is_file():
            raise SystemExit("missing paper/tex/figures/data/ring-head2head.json")
        systems(lines)
        classify(lines)
        ring_record(lines)
        interface(lines)
    except SystemExit as exc:
        code = int(exc.code) if isinstance(exc.code, int) else 1
        lines.append(f"error {exc}")
    except Exception as exc:  # noqa: BLE001
        code = 1
        lines.append(f"error {type(exc).__name__}: {exc}")
    lines.append(f"exit_code {code}")
    write_log(lines)
    return code


if __name__ == "__main__":
    sys.exit(main())
