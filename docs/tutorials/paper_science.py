#!/usr/bin/env python3
"""Build the five paper tutorials from GenIce, packmol, or a short OpenMM run.

Usage: paper_science.py TOPIC SEAMS LIBDIR OUTDIR
TOPIC is interface, cages, rings, brine, or nucleation.
The same names with a -full suffix are the longer runs.
Exit 77 when the engine that topic needs is absent.
"""

import math
import os
import random
import shutil
import subprocess
import sys
from pathlib import Path


def skip(message):
    print(message)
    raise SystemExit(77)


def need(program):
    sibling = Path(sys.executable).resolve().parent / program
    if sibling.is_file():
        return str(sibling)
    found = shutil.which(program)
    if found:
        return found
    skip("{} is absent".format(program))


def seams_env(libdir):
    env = os.environ.copy()
    env["NO_COLOR"] = "1"
    env["OMP_NUM_THREADS"] = "2"
    parts = [str(libdir)]
    if env.get("LD_LIBRARY_PATH"):
        parts.append(env["LD_LIBRARY_PATH"])
    env["LD_LIBRARY_PATH"] = os.pathsep.join(parts)
    return env


def prepare_vulkan():
    if os.environ.get("VK_ICD_FILENAMES"):
        return
    folders = [Path("/usr/share/vulkan/icd.d")]
    prefix = os.environ.get("CONDA_PREFIX")
    if prefix:
        folders.insert(0, Path(prefix) / "share" / "vulkan" / "icd.d")
    for folder in folders:
        if not folder.is_dir():
            continue
        matches = sorted(folder.glob("lvp_icd*.json"))
        if matches:
            os.environ["VK_ICD_FILENAMES"] = str(matches[0])
            os.environ["VK_DRIVER_FILES"] = str(matches[0])
            return


def render(dump, png):
    prepare_vulkan()
    try:
        import warnings

        warnings.filterwarnings("ignore", message=".*OVITO.*PyPI")
        from ovito.io import import_file
        from ovito.vis import Viewport
    except ImportError:
        skip("ovito is absent")
    pipeline = import_file(str(dump))
    pipeline.add_to_scene()
    view = Viewport(type=Viewport.Type.Ortho)
    view.zoom_all()
    view.render_image(filename=str(png), size=(480, 360), background=(1, 1, 1))
    if not png.is_file() or png.stat().st_size > 300 * 1024:
        raise SystemExit("figure {} is missing or larger than 300 KB".format(png.name))


def genice_oxygens(lattice, rep):
    genice = need("genice2")
    gro = subprocess.check_output(
        [genice, lattice, "-r", *(str(v) for v in rep), "--format", "gromacs", "--water", "tip4p", "--quiet"],
        text=True,
        stderr=subprocess.DEVNULL,
    )
    xyz = subprocess.check_output(
        [genice, lattice, "-r", *(str(v) for v in rep), "--format", "xyz", "--water", "tip4p", "--quiet"],
        text=True,
        stderr=subprocess.DEVNULL,
    )
    grows = gro.splitlines()
    box_nm = [float(v) for v in grows[int(grows[1]) + 2].split()]
    box = [value * 10.0 for value in box_nm]
    rows = xyz.splitlines()
    index = 0
    while rows[index].startswith("#"):
        index += 1
    count = int(rows[index])
    index += 2
    oxygens = []
    for line in rows[index : index + count]:
        parts = line.split()
        if parts and parts[0] == "OW":
            oxygens.append([float(value) for value in parts[1:4]])
    return oxygens, box


def write_dump(path, frames):
    """frames is a list of (oxygens, box, extras). extras are (type, x, y, z)."""
    chunks = []
    for step, (oxygens, box, extras) in enumerate(frames):
        if len(box) == 3:
            oxygens = [
                (x % box[0], y % box[1], z % box[2]) for x, y, z in oxygens
            ]
        points = [(1, *xyz) for xyz in oxygens] + list(extras)
        if len(box) == 3:
            bounds = [
                "ITEM: BOX BOUNDS pp pp pp",
                "0.0 {:.8f}".format(box[0]),
                "0.0 {:.8f}".format(box[1]),
                "0.0 {:.8f}".format(box[2]),
            ]
        else:
            lx, ly, lz, xy, xz, yz = box
            xlo_b = min(0.0, xy, xz, xy + xz)
            xhi_b = lx + max(0.0, xy, xz, xy + xz)
            ylo_b = min(0.0, yz)
            yhi_b = ly + max(0.0, yz)
            bounds = [
                "ITEM: BOX BOUNDS xy xz yz pp pp pp",
                "{:.8f} {:.8f} {:.8f}".format(xlo_b, xhi_b, xy),
                "{:.8f} {:.8f} {:.8f}".format(ylo_b, yhi_b, xz),
                "0.0 {:.8f} {:.8f}".format(lz, yz),
            ]
        chunks.append("ITEM: TIMESTEP")
        chunks.append(str(step))
        chunks.append("ITEM: NUMBER OF ATOMS")
        chunks.append(str(len(points)))
        chunks.extend(bounds)
        chunks.append("ITEM: ATOMS id type x y z")
        for index, (kind, x, y, z) in enumerate(points, 1):
            chunks.append("{} {} {:.6f} {:.6f} {:.6f}".format(index, kind, x, y, z))
    path.write_text("\n".join(chunks) + "\n")


def box_of(lattice_box):
    if len(lattice_box) <= 3:
        return lattice_box[:3]
    xx, yy, zz, _gxy, _gxz, gyx, gyz, gzx, _gzy = lattice_box
    return [xx, yy, zz, gyx, gzx, gyz]


def run_seams(seams, libdir, args):
    done = subprocess.run(
        [seams, *args],
        check=False,
        capture_output=True,
        text=True,
        env=seams_env(libdir),
    )
    if done.returncode != 0:
        sys.stderr.write(done.stdout)
        sys.stderr.write(done.stderr)
        raise SystemExit(done.returncode)
    return done.stdout


def counts(text):
    found = {}
    for line in text.splitlines():
        parts = line.split()
        if not parts or parts[0] != "nop":
            continue
        for index in range(2, len(parts) - 1, 2):
            found[parts[index]] = int(parts[index + 1])
    return found


def openmm_relax(oxygens, box, steps, restrained, stiffness=1000.0):
    try:
        import openmm as mm
        import openmm.unit as unit
    except ImportError:
        skip("openmm is absent")
    system = mm.System()
    for _ in oxygens:
        system.addParticle(18.015 * unit.dalton)
    if len(box) == 3:
        lx, ly, lz = box
        xy = xz = yz = 0.0
    else:
        lx, ly, lz, xy, xz, yz = box
    system.setDefaultPeriodicBoxVectors(
        (lx, 0, 0) * unit.angstrom,
        (xy, ly, 0) * unit.angstrom,
        (xz, yz, lz) * unit.angstrom,
    )
    force = mm.NonbondedForce()
    force.setNonbondedMethod(mm.NonbondedForce.CutoffPeriodic)
    force.setCutoffDistance(0.35 * unit.nanometer)
    for _ in oxygens:
        force.addParticle(0.0, 0.3165 * unit.nanometer, 0.650 * unit.kilojoule_per_mole)
    system.addForce(force)
    if restrained:
        hold = mm.CustomExternalForce("k*((x-x0)^2+(y-y0)^2+(z-z0)^2)")
        hold.addGlobalParameter("k", float(stiffness))
        hold.addPerParticleParameter("x0")
        hold.addPerParticleParameter("y0")
        hold.addPerParticleParameter("z0")
        for index in restrained:
            x, y, z = oxygens[index]
            hold.addParticle(index, [x * 0.1, y * 0.1, z * 0.1])
        system.addForce(hold)
    integrator = mm.VerletIntegrator(0.002 * unit.picoseconds)
    platform = mm.Platform.getPlatformByName("Reference")
    context = mm.Context(system, integrator, platform)
    context.setPositions([(x * 0.1, y * 0.1, z * 0.1) for x, y, z in oxygens] * unit.nanometer)
    mm.LocalEnergyMinimizer.minimize(context, 10.0, steps)
    state = context.getState(getPositions=True)
    raw = state.getPositions(asNumpy=True).value_in_unit(unit.angstrom)
    return [tuple(float(v) for v in point) for point in raw]


def interface(seams, libdir, out):
    crystal, box = genice_oxygens("1h", (2, 2, 1))
    length = box_of(box)
    liquid = []
    rng = random.Random(7)
    for x, y, z in crystal:
        liquid.append(
            (
                x + rng.uniform(-0.8, 0.8),
                y + rng.uniform(-0.8, 0.8),
                z + length[2] + rng.uniform(-0.4, 0.4),
            )
        )
    stacked = list(crystal) + liquid
    tall = [length[0], length[1], length[2] * 2.0]
    relaxed = openmm_relax(stacked, tall, 8, list(range(len(crystal))), stiffness=1.0e5)
    dump = out / "interface.lammpstrj"
    write_dump(dump, [(relaxed, tall, [])])
    published = counts(run_seams(seams, libdir, ["chill-plus", str(dump), "--type", "1", "--cutoff", "3.5", "--rule", "published"]))
    ovito = counts(run_seams(seams, libdir, ["ovito-chill+", str(dump), "--type", "1", "--cutoff", "3.5"]))
    print(
        "interface published hexagonal {} interfacial {} water {}".format(
            published.get("hexagonal", 0),
            published.get("interfacial", 0),
            published.get("water", 0),
        )
    )
    print(
        "interface ovito hexagonal {} interfacial {} water {}".format(
            ovito.get("hexagonal", 0),
            ovito.get("interfacial", 0),
            ovito.get("water", 0),
        )
    )
    render(dump, out / "interface.png")


def cages(seams, libdir, out):
    jobs = (
        ("CS1", (1, 1, 1), "512+51262", None, "sI"),
        ("CS2", (1, 1, 1), "512+51264", None, "sII"),
        ("HS3", (1, 1, 1), "512+sh+51268", "none", "sH"),
    )
    for lattice, rep, signature, edge, name in jobs:
        oxygens, gro_box = genice_oxygens(lattice, rep)
        cell = box_of(gro_box) if len(gro_box) > 3 else gro_box
        dump = out / "{}.lammpstrj".format(name)
        write_dump(dump, [(oxygens, cell, [])])
        args = ["cages", str(dump), "--type", "1", "--cutoff", "3.5", "--graph", "cutoff", "--signature", signature]
        if edge:
            args.extend(["--cage-edge", edge])
        text = run_seams(seams, libdir, args)
        numbers = []
        for line in text.splitlines():
            if " cages " in line:
                numbers.append(line.split(" cages ")[1].split()[0])
        print("cages {} {}".format(name, " ".join(numbers)))
    render(out / "sI.lammpstrj", out / "cages.png")


def neighbour_graph(oxygens, box, cutoff=3.5):
    lx, ly, lz = box[:3]
    graph = [[] for _ in oxygens]
    for i, (x, y, z) in enumerate(oxygens):
        for j in range(i + 1, len(oxygens)):
            dx = x - oxygens[j][0]
            dy = y - oxygens[j][1]
            dz = z - oxygens[j][2]
            dx -= lx * round(dx / lx)
            dy -= ly * round(dy / ly)
            dz -= lz * round(dz / lz)
            if dx * dx + dy * dy + dz * dz <= cutoff * cutoff:
                graph[i].append(j)
                graph[j].append(i)
    return graph


def simple_cycles(graph, limit=6):
    cycles = set()

    def walk(start, current, path, seen):
        if len(path) > limit:
            return
        for nxt in graph[current]:
            if nxt == start and len(path) >= 3:
                ring = tuple(sorted(path))
                cycles.add(ring)
            elif nxt > start and nxt not in seen:
                seen.add(nxt)
                path.append(nxt)
                walk(start, nxt, path, seen)
                path.pop()
                seen.remove(nxt)

    for start in range(len(graph)):
        walk(start, start, [start], {start})
    return cycles


def graph_distance(graph, source, limit):
    dist = {source: 0}
    layer = [source]
    while layer:
        nxt = []
        for node in layer:
            if dist[node] >= limit:
                continue
            for other in graph[node]:
                if other not in dist:
                    dist[other] = dist[node] + 1
                    nxt.append(other)
        layer = nxt
    return dist


def shortest_path_rings(graph, cycles):
    kept = []
    for ring in cycles:
        ok = True
        members = list(ring)
        index = {atom: i for i, atom in enumerate(members)}
        # The stored ring is a sorted tuple, so rebuild one cyclic order by a walk.
        order = [members[0]]
        seen = {members[0]}
        while len(order) < len(members):
            found = None
            for other in graph[order[-1]]:
                if other in index and other not in seen:
                    found = other
                    break
            if found is None:
                ok = False
                break
            order.append(found)
            seen.add(found)
        if not ok:
            continue
        place = {atom: i for i, atom in enumerate(order)}
        size = len(order)
        for a_i, atom in enumerate(order):
            dist = graph_distance(graph, atom, size // 2)
            for b_i in range(a_i + 1, size):
                other = order[b_i]
                along = min(b_i - a_i, size - (b_i - a_i))
                if dist.get(other, size) < along:
                    ok = False
                    break
            if not ok:
                break
        if ok:
            kept.append(ring)
    return kept


def rings(seams, libdir, out):
    oxygens, gro_box = genice_oxygens("1h", (1, 1, 1))
    cell = gro_box[:3]
    dump = out / "rings.lammpstrj"
    write_dump(dump, [(oxygens, cell, [])])
    text = run_seams(seams, libdir, ["rings", str(dump), "--type", "1", "--cutoff", "3.5", "--max", "6"])
    primitive = {}
    for line in text.splitlines():
        parts = line.split()
        if len(parts) == 2 and parts[0].isdigit():
            primitive[int(parts[0])] = int(parts[1])
    graph = neighbour_graph(oxygens, cell)
    cycles = simple_cycles(graph, 6)
    short = shortest_path_rings(graph, cycles)
    def hist(items):
        found = {}
        for ring in items:
            found[len(ring)] = found.get(len(ring), 0) + 1
        return found
    all_hist = hist(cycles)
    short_hist = hist(short)
    print(
        "rings primitive 4 {} 6 {}".format(primitive.get(4, 0), primitive.get(6, 0))
    )
    print(
        "rings shortest 4 {} 6 {}".format(short_hist.get(4, 0), short_hist.get(6, 0))
    )
    print(
        "rings searched 4 {} 6 {}".format(all_hist.get(4, 0), all_hist.get(6, 0))
    )
    render(dump, out / "rings.png")


def brine(seams, libdir, out):
    packmol = need("packmol")
    crystal, box = genice_oxygens("1h", (2, 2, 1))
    length = box_of(box)
    ice_xyz = out / "ice.xyz"
    brine_xyz = out / "brine.xyz"
    ice_xyz.write_text(
        "{}\nice\n".format(len(crystal))
        + "\n".join("O {:.5f} {:.5f} {:.5f}".format(*point) for point in crystal)
        + "\n"
    )
    # A handful of waters and ions in the region above the slab. Packmol
    # places them. The ice block stays fixed, which is why the front has
    # a crystal to grow from.
    waters = ["O {:.5f} {:.5f} {:.5f}".format(1.5 + i * 1.1, 1.5, length[2] + 2.0) for i in range(6)]
    ions = [
        "Na {:.5f} {:.5f} {:.5f}".format(2.0, 3.0, length[2] + 4.0),
        "Cl {:.5f} {:.5f} {:.5f}".format(4.0, 3.0, length[2] + 4.5),
    ]
    brine_xyz.write_text("{}\nbrine\n".format(len(waters) + len(ions)) + "\n".join(waters + ions) + "\n")
    deck = out / "brine.inp"
    z0 = length[2] + 1.0
    z1 = length[2] * 2.0 - 1.0
    deck.write_text(
        "\n".join(
            [
                "tolerance 2.0",
                "filetype xyz",
                "output {}".format(out / "packed.xyz"),
                "structure {}".format(ice_xyz),
                "  number 1",
                "  fixed 0. 0. 0. 0. 0. 0.",
                "end structure",
                "structure {}".format(brine_xyz),
                "  number 1",
                "  inside box 0. 0. {:.3f} {:.3f} {:.3f} {:.3f}".format(z0, length[0], length[1], z1),
                "end structure",
                "",
            ]
        )
    )
    done = subprocess.run(
        [packmol],
        stdin=deck.open(),
        text=True,
        capture_output=True,
        check=False,
    )
    if done.returncode != 0 or not (out / "packed.xyz").is_file():
        sys.stderr.write(done.stdout[-500:])
        sys.stderr.write(done.stderr[-500:])
        raise SystemExit(done.returncode or 1)
    packed = (out / "packed.xyz").read_text().splitlines()
    n = int(packed[0])
    oxygens = []
    extras = []
    for line in packed[2 : 2 + n]:
        kind, x, y, z = line.split()[:4]
        point = (float(x), float(y), float(z))
        if kind == "O":
            oxygens.append(point)
        else:
            extras.append((2 if kind == "Na" else 3, *point))
    tall = [length[0], length[1], length[2] * 2.0]
    frames = []
    current = list(oxygens)
    ice_index = list(range(len(crystal)))
    for stage in range(3):
        # Move the nearest brine oxygen onto the next lattice plane and
        # restrain it. That is the growth step. The relax lets the ions
        # and the remaining liquid answer it.
        if stage:
            targets = [(x, y, z + 3.6 * stage) for x, y, z in crystal]
            for target in targets:
                best = None
                best_d = 2.4 * 2.4
                for index in range(len(crystal), len(current)):
                    dx = current[index][0] - target[0]
                    dy = current[index][1] - target[1]
                    dz = current[index][2] - target[2]
                    dist = dx * dx + dy * dy + dz * dz
                    if dist < best_d:
                        best_d = dist
                        best = index
                if best is not None:
                    current[best] = target
                    ice_index.append(best)
        ice_index = sorted(set(ice_index))
        current = list(openmm_relax(current, tall, 8, ice_index, stiffness=1.0e5))
        frames.append((list(current), tall, extras))
    dump = out / "brine.lammpstrj"
    write_dump(dump, frames)
    fronts = []
    for frame in range(3):
        text = run_seams(
            seams,
            libdir,
            [
                "chill-plus",
                str(dump),
                "--type",
                "1",
                "--cutoff",
                "3.5",
                "--frame",
                str(frame + 1),
                "--last",
                str(frame + 1),
                "--layers",
                "--axis",
                "2",
            ],
        )
        stack = ""
        for line in text.splitlines():
            if " stack " in line:
                stack = line.split(" stack ", 1)[1].split()[0]
        # H and C are crystal slabs. The front is the first slab that is not.
        front = len(stack) * 3.7
        for index, mark in enumerate(stack):
            if mark not in ("H", "C"):
                front = index * 3.7
                break
        fronts.append("{:.2f}".format(front))
        found = counts(text)
        print(
            "brine frame {} hexagonal {} water {} front {}".format(
                frame, found.get("hexagonal", 0), found.get("water", 0), fronts[-1]
            )
        )
    render(dump, out / "brine.png")


def nucleation(seams, libdir, out):
    crystal, box = genice_oxygens("1h", (2, 2, 1))
    length = box_of(box)
    rng = random.Random(3)
    liquid = [
        (
            min(max(x + rng.uniform(-1.4, 1.4), 0.2), length[0] - 0.2),
            min(max(y + rng.uniform(-1.4, 1.4), 0.2), length[1] - 0.2),
            min(max(z + rng.uniform(-1.4, 1.4), 0.2), length[2] - 0.2),
        )
        for x, y, z in crystal
    ]
    # The nucleus is the crystal copy of the first quarter of the cell.
    # Each saved frame returns more of those molecules to the lattice,
    # then OpenMM relaxes the rest. That is the short nucleation run.
    nucleus = [index for index, (_x, _y, z) in enumerate(crystal) if z < length[2] * 0.5]
    frames = []
    current = list(liquid)
    for stage, count in enumerate((0, len(nucleus) // 2, len(nucleus))):
        for index in nucleus[:count]:
            current[index] = tuple(crystal[index])
        current = list(openmm_relax(current, length, 2, nucleus[:count], stiffness=1.0e5))
        frames.append((current, length, []))
    dump = out / "nucleation.lammpstrj"
    write_dump(dump, frames)
    for frame in range(3):
        text = run_seams(
            seams,
            libdir,
            [
                "steinhardt",
                str(dump),
                "--type",
                "1",
                "--cutoff",
                "3.5",
                "--l",
                "6",
                "--frame",
                str(frame + 1),
                "--last",
                str(frame + 1),
            ],
        )
        q6bar = "nan"
        for line in text.splitlines():
            if " q6bar " in line:
                q6bar = line.split(" q6bar ", 1)[1].split()[0]
        print("nucleation frame {} q6bar {}".format(frame, q6bar))
    render(dump, out / "nucleation.png")


def publish(out, lines, dump, png):
    render(dump, png)
    text = "\n".join(lines) + "\n"
    (out / "report.txt").write_text(text)
    sys.stdout.write(text)


def openmm_modules():
    os.environ["OPENMM_CPU_THREADS"] = "1"
    try:
        import openmm as mm
        import openmm.unit as unit
    except ImportError:
        skip("openmm is absent")
    return mm, unit


def reference_platform(mm):
    return mm.Platform.getPlatformByName("Reference")


def stack_front(text):
    stack = ""
    phi = ""
    for line in text.splitlines():
        if "phi-c" in line:
            phi = line.split("phi-c", 1)[1].split()[0]
        if " stack " in line:
            stack = line.split(" stack ", 1)[1].split()[0]
    front = len(stack) * 3.7
    for index, mark in enumerate(stack):
        if mark not in ("H", "C"):
            front = index * 3.7
            break
    return "{:.2f}".format(front), stack, phi


def label_line(prefix, found, names):
    parts = [prefix]
    for name in names:
        parts.append("{} {}".format(name, found.get(name, 0)))
    return " ".join(parts)


def interface_full(seams, libdir, out):
    mm, unit = openmm_modules()
    crystal, gro_box = genice_oxygens("1h", (4, 4, 4))
    length = box_of(gro_box)
    rng = random.Random(11)
    liquid = [
        (
            rng.uniform(0.5, length[0] - 0.5),
            rng.uniform(0.5, length[1] - 0.5),
            length[2] + rng.uniform(1.0, length[2] - 1.0),
        )
        for _ in crystal
    ]
    pts = list(crystal) + liquid
    tall = [length[0], length[1], length[2] * 2.0]
    system = mm.System()
    for _ in pts:
        system.addParticle(18.015 * unit.dalton)
    system.setDefaultPeriodicBoxVectors(
        (tall[0], 0, 0) * unit.angstrom,
        (0, tall[1], 0) * unit.angstrom,
        (0, 0, tall[2]) * unit.angstrom,
    )
    force = mm.NonbondedForce()
    force.setNonbondedMethod(mm.NonbondedForce.CutoffPeriodic)
    force.setCutoffDistance(0.35 * unit.nanometer)
    for _ in pts:
        force.addParticle(0.0, 0.3165 * unit.nanometer, 2.0 * unit.kilojoule_per_mole)
    system.addForce(force)
    hold = mm.CustomExternalForce("k*((x-x0)^2+(y-y0)^2+(z-z0)^2)")
    hold.addGlobalParameter("k", 5.0e4)
    hold.addPerParticleParameter("x0")
    hold.addPerParticleParameter("y0")
    hold.addPerParticleParameter("z0")
    for index, (x, y, z) in enumerate(crystal):
        hold.addParticle(index, [x * 0.1, y * 0.1, z * 0.1])
    system.addForce(hold)
    integrator = mm.LangevinMiddleIntegrator(
        240 * unit.kelvin, 1 / unit.picosecond, 0.002 * unit.picoseconds
    )
    integrator.setRandomNumberSeed(240)
    context = mm.Context(system, integrator, reference_platform(mm))
    context.setPositions([(x * 0.1, y * 0.1, z * 0.1) for x, y, z in pts] * unit.nanometer)
    mm.LocalEnergyMinimizer.minimize(context, 10.0, 40)
    integrator.step(400)
    raw = context.getState(getPositions=True).getPositions(asNumpy=True).value_in_unit(unit.angstrom)
    relaxed = [tuple(float(v) for v in point) for point in raw]
    dump = out / "interface-full.lammpstrj"
    write_dump(dump, [(relaxed, tall, [])])
    published_text = run_seams(
        seams,
        libdir,
        ["chill-plus", str(dump), "--type", "1", "--cutoff", "3.5", "--rule", "published", "--layers", "--axis", "2"],
    )
    ovito_text = run_seams(
        seams, libdir, ["ovito-chill+", str(dump), "--type", "1", "--cutoff", "3.5"]
    )
    published = counts(published_text)
    ovito = counts(ovito_text)
    _front, stack, phi = stack_front(published_text)
    names = ("hexagonal", "interfacial", "water", "clathrate", "interClathrate")
    lines = [
        "interface full waters {}".format(len(relaxed)),
        label_line("interface full published", published, names),
        label_line("interface full ovito", ovito, names),
        "interface full phi-c {}".format(phi),
        "interface full stack {}".format(stack),
    ]
    publish(out, lines, dump, out / "interface-full.png")


def cages_full(seams, libdir, out):
    jobs = (
        ("CS1", (2, 2, 2), "512+51262", None, "sI"),
        ("CS2", (2, 2, 1), "512+51264", None, "sII"),
        ("HS3", (2, 2, 1), "512+sh+51268", "none", "sH"),
    )
    lines = []
    for lattice, rep, signature, edge, name in jobs:
        oxygens, gro_box = genice_oxygens(lattice, rep)
        cell = box_of(gro_box) if len(gro_box) > 3 else gro_box
        dump = out / "{}-full.lammpstrj".format(name)
        write_dump(dump, [(oxygens, cell, [])])
        args = ["cages", str(dump), "--type", "1", "--cutoff", "3.5", "--graph", "cutoff", "--signature", signature]
        if edge:
            args.extend(["--cage-edge", edge])
        text = run_seams(seams, libdir, args)
        numbers = []
        for line in text.splitlines():
            if " cages " in line:
                numbers.append(line.split(" cages ")[1].split()[0])
        lines.append("cages full {} {}".format(name, " ".join(numbers)))
    publish(out, lines, out / "sI-full.lammpstrj", out / "cages-full.png")


def rings_full(seams, libdir, out):
    oxygens, gro_box = genice_oxygens("1h", (8, 8, 4))
    cell = gro_box[:3]
    dump = out / "rings-full.lammpstrj"
    write_dump(dump, [(oxygens, cell, [])])
    text = run_seams(seams, libdir, ["rings", str(dump), "--type", "1", "--cutoff", "3.5", "--max", "8"])
    nop = "0"
    total = "0"
    histogram = []
    for line in text.splitlines():
        parts = line.split()
        if parts and parts[0] == "nop" and "rings" in parts:
            nop = parts[1]
            total = parts[3]
        elif len(parts) == 2 and parts[0].isdigit():
            histogram.append("{} {}".format(parts[0], parts[1]))
    lines = [
        "rings full nop {} rings {}".format(nop, total),
        "rings full histogram {}".format(" ".join(histogram)),
    ]
    publish(out, lines, dump, out / "rings-full.png")


def brine_full(seams, libdir, out):
    mm, unit = openmm_modules()
    ice, gro_box = genice_oxygens("1h", (2, 2, 2))
    box = box_of(gro_box)
    rng = random.Random(5)
    layer = [point for point in ice if point[2] < box[2] * 0.5]

    def jitter(zshift):
        shifted = []
        for x, y, z in layer:
            shifted.append(
                (
                    x + rng.uniform(-0.6, 0.6),
                    y + rng.uniform(-0.6, 0.6),
                    z + zshift + rng.uniform(-0.4, 0.4),
                )
            )
        return shifted

    liquid = jitter(box[2]) + jitter(1.55 * box[2])
    waters = [tuple(point) for point in ice] + liquid
    ions = [
        (2, 2.2, 2.2, box[2] * 2.05),
        (3, 8.0, 3.0, box[2] * 2.12),
        (2, 4.5, 9.0, box[2] * 2.18),
        (3, 11.0, 6.5, box[2] * 2.22),
    ]
    everyone = waters + [(x, y, z) for _kind, x, y, z in ions]
    tall = [box[0], box[1], box[2] * 2.4]
    system = mm.System()
    for _ in everyone:
        system.addParticle(22.0 * unit.dalton)
    system.setDefaultPeriodicBoxVectors(
        (tall[0], 0, 0) * unit.angstrom,
        (0, tall[1], 0) * unit.angstrom,
        (0, 0, tall[2]) * unit.angstrom,
    )
    nb = mm.NonbondedForce()
    nb.setNonbondedMethod(mm.NonbondedForce.CutoffPeriodic)
    nb.setCutoffDistance(0.4 * unit.nanometer)
    nb.setUseSwitchingFunction(True)
    nb.setSwitchingDistance(0.32 * unit.nanometer)
    for _ in waters:
        nb.addParticle(0.0, 0.31 * unit.nanometer, 1.2 * unit.kilojoule_per_mole)
    for kind, _x, _y, _z in ions:
        charge = 0.3 if kind == 2 else -0.3
        nb.addParticle(charge * unit.elementary_charge, 0.40 * unit.nanometer, 0.15 * unit.kilojoule_per_mole)
    system.addForce(nb)

    def spring(name, pairs):
        hold = mm.CustomExternalForce("{}*((x-x0)^2+(y-y0)^2+(z-z0)^2)".format(name))
        hold.addGlobalParameter(name, 0.0)
        hold.addPerParticleParameter("x0")
        hold.addPerParticleParameter("y0")
        hold.addPerParticleParameter("z0")
        for index, (x, y, z) in pairs:
            hold.addParticle(index, [x * 0.1, y * 0.1, z * 0.1])
        return hold

    system.addForce(spring("kice", list(enumerate(ice))))
    system.addForce(
        spring("k1", [(len(ice) + index, (x, y, z + box[2])) for index, (x, y, z) in enumerate(layer)])
    )
    system.addForce(
        spring(
            "k2",
            [(len(ice) + len(layer) + index, (x, y, z + 1.55 * box[2])) for index, (x, y, z) in enumerate(layer)],
        )
    )
    integrator = mm.LangevinMiddleIntegrator(
        240 * unit.kelvin, 5 / unit.picosecond, 0.001 * unit.picoseconds
    )
    integrator.setRandomNumberSeed(240)
    context = mm.Context(system, integrator, reference_platform(mm))
    context.setPositions([(x * 0.1, y * 0.1, z * 0.1) for x, y, z in everyone] * unit.nanometer)
    context.setParameter("kice", 5.0e4)
    context.setParameter("k1", 0.0)
    context.setParameter("k2", 0.0)

    def settle(name=None, value=0.0):
        if name:
            context.setParameter(name, value)
        mm.LocalEnergyMinimizer.minimize(context, 5.0, 30)
        integrator.step(800)

    saved = []

    def grab():
        state = context.getState(getPositions=True).getPositions(asNumpy=True).value_in_unit(unit.angstrom)
        oxygens = [tuple(float(v) for v in point) for point in state[: len(waters)]]
        extras = []
        for offset, (kind, _x, _y, _z) in enumerate(ions):
            x, y, z = (float(v) for v in state[len(waters) + offset])
            extras.append((kind, x, y, z))
        saved.append((oxygens, tall, extras))

    settle()
    grab()
    settle("k1", 4.0e4)
    grab()
    settle("k2", 4.0e4)
    grab()
    dump = out / "brine-full.lammpstrj"
    write_dump(dump, saved)
    lines = []
    for frame in range(3):
        text = run_seams(
            seams,
            libdir,
            [
                "chill-plus",
                str(dump),
                "--type",
                "1",
                "--cutoff",
                "3.5",
                "--frame",
                str(frame + 1),
                "--last",
                str(frame + 1),
                "--layers",
                "--axis",
                "2",
            ],
        )
        front, stack, _phi = stack_front(text)
        found = counts(text)
        lines.append(
            "brine full frame {} hexagonal {} water {} front {} stack {}".format(
                frame, found.get("hexagonal", 0), found.get("water", 0), front, stack
            )
        )
    figure = out / "brine-full-last.lammpstrj"
    write_dump(figure, [saved[-1]])
    publish(out, lines, figure, out / "brine-full.png")


def nucleation_full(seams, libdir, out):
    crystal, gro_box = genice_oxygens("1h", (4, 4, 4))
    length = box_of(gro_box)
    rng = random.Random(19)
    liquid = []
    for x, y, z in crystal:
        liquid.append(
            (
                min(max(x + rng.uniform(-1.6, 1.6), 0.3), length[0] - 0.3),
                min(max(y + rng.uniform(-1.6, 1.6), 0.3), length[1] - 0.3),
                min(max(z + rng.uniform(-1.6, 1.6), 0.3), length[2] - 0.3),
            )
        )
    basal = [index for index, (_x, _y, z) in enumerate(crystal) if z < length[2] * 0.5]
    current = list(liquid)
    frames = []
    for stage, count in enumerate((0, len(basal) // 2, len(crystal))):
        restored = list(range(len(crystal))) if count == len(crystal) else basal[:count]
        for index in restored:
            current[index] = tuple(crystal[index])
        current = list(openmm_relax(current, length, 2 if stage == 0 else 8, restored, stiffness=1.0e5))
        frames.append((list(current), length, []))
    dump = out / "nucleation-full.lammpstrj"
    write_dump(dump, frames)
    lines = []
    for frame in range(3):
        text = run_seams(
            seams,
            libdir,
            [
                "steinhardt",
                str(dump),
                "--type",
                "1",
                "--cutoff",
                "3.5",
                "--l",
                "6",
                "--frame",
                str(frame + 1),
                "--last",
                str(frame + 1),
            ],
        )
        q6bar = "nan"
        for line in text.splitlines():
            if " q6bar " in line:
                q6bar = line.split(" q6bar ", 1)[1].split()[0]
        lines.append("nucleation full frame {} q6bar {}".format(frame, q6bar))
    figure = out / "nucleation-full-last.lammpstrj"
    write_dump(figure, [frames[-1]])
    publish(out, lines, figure, out / "nucleation-full.png")


def main():
    topic, seams, libdir, out = sys.argv[1], sys.argv[2], sys.argv[3], Path(sys.argv[4])
    out.mkdir(parents=True, exist_ok=True)
    topics = {
        "interface": interface,
        "cages": cages,
        "rings": rings,
        "brine": brine,
        "nucleation": nucleation,
        "interface-full": interface_full,
        "cages-full": cages_full,
        "rings-full": rings_full,
        "brine-full": brine_full,
        "nucleation-full": nucleation_full,
    }
    if topic not in topics:
        raise SystemExit("unknown topic {}".format(topic))
    topics[topic](seams, libdir, out)


if __name__ == "__main__":
    main()
