#!/usr/bin/env python3
"""Science tutorials: interface, brine front, and hydrate cages.

Python and the C program in science_tutorials.c both call libseams.
Layer period, layer kind, and cubicity follow scripts/science_render.py
at 702d85b2454e459dd14d9ec5fa3bd6b331e73862.

Usage: science_tutorials.py ROOT [BUILD]
"""

from __future__ import annotations

import ctypes
import hashlib
import math
import os
import shutil
import subprocess
import sys
import urllib.error
import urllib.request
from pathlib import Path

import numpy as np

NAMES = (
    "cubic",
    "hexagonal",
    "water",
    "interfacial",
    "clathrate",
    "interClathrate",
    "unclassified",
    "reCubic",
    "reHex",
)
BRINE_REV = "4828e33ecf419c70643b2d466bfa4f283b4138dc"
BRINE_URL = (
    "https://huggingface.co/datasets/d-SEAMS/brine-front-230K/resolve/" + BRINE_REV
)
BRINE_FILES = {
    "data/frames/T230_m20_r0_t0.xyz":
        "dc55832fcae664aa5e528283b3be88762efdf1dabf4db1520bb1c806401b6128",
    "data/frames/T230_m20_r0_tend.xyz":
        "dcbc05f81f321b45195d259cd4c22e69340f9c1ae57d93c47a74c995c58dd05d",
}
SKIP = (
    "brine tutorial skipped: cannot fetch d-SEAMS/brine-front-230K "
    "at 4828e33ecf419c70643b2d466bfa4f283b4138dc"
)


def _lib(libdir: Path):
    lib = ctypes.CDLL(str(libdir / "libseams.so"))
    lib.seams_system_create.restype = ctypes.c_void_p
    lib.seams_system_destroy.argtypes = [ctypes.c_void_p]
    lib.seams_system_read.argtypes = [
        ctypes.c_void_p, ctypes.c_char_p, ctypes.c_int, ctypes.c_int
    ]
    lib.seams_system_read.restype = ctypes.c_int
    lib.seams_system_count.argtypes = [ctypes.c_void_p]
    lib.seams_system_count.restype = ctypes.c_int
    lib.seams_system_set_frame.argtypes = [
        ctypes.c_void_p,
        ctypes.POINTER(ctypes.c_double),
        ctypes.c_int,
        ctypes.POINTER(ctypes.c_double),
        ctypes.POINTER(ctypes.c_int),
        ctypes.c_void_p,
    ]
    lib.seams_system_set_frame.restype = ctypes.c_int
    lib.seams_chill_plus.argtypes = [
        ctypes.c_void_p,
        ctypes.c_double,
        ctypes.c_int,
        ctypes.c_int,
        ctypes.c_int,
        ctypes.POINTER(ctypes.c_int),
        ctypes.POINTER(ctypes.c_int),
    ]
    lib.seams_chill_plus.restype = ctypes.c_int
    lib.seams_cages_signature.argtypes = [
        ctypes.c_void_p,
        ctypes.c_double,
        ctypes.c_int,
        ctypes.c_int,
        ctypes.c_char_p,
        ctypes.POINTER(ctypes.c_size_t),
        ctypes.POINTER(ctypes.c_size_t),
    ]
    lib.seams_cages_signature.restype = ctypes.c_int
    lib.seams_brine_analyse.argtypes = [
        ctypes.c_void_p,
        ctypes.POINTER(ctypes.c_int),
        ctypes.c_size_t,
        ctypes.c_int,
        ctypes.c_int,
        ctypes.c_int,
        ctypes.c_int,
        ctypes.POINTER(ctypes.c_void_p),
    ]
    lib.seams_brine_analyse.restype = ctypes.c_int
    lib.seams_brine_destroy.argtypes = [ctypes.c_void_p]
    lib.seams_brine_fronts.argtypes = [
        ctypes.c_void_p,
        ctypes.POINTER(ctypes.c_double),
        ctypes.POINTER(ctypes.c_double),
        ctypes.POINTER(ctypes.c_double),
        ctypes.POINTER(ctypes.c_int),
    ]
    lib.seams_brine_fronts.restype = ctypes.c_int
    lib.seams_brine_rejection.argtypes = [
        ctypes.c_void_p,
        ctypes.POINTER(ctypes.c_double),
        ctypes.POINTER(ctypes.c_double),
        ctypes.POINTER(ctypes.c_double),
        ctypes.POINTER(ctypes.c_int),
    ]
    lib.seams_brine_rejection.restype = ctypes.c_int
    lib.seams_last_error.restype = ctypes.c_char_p
    lib.seams_set_num_threads.argtypes = [ctypes.c_int]
    lib.seams_set_num_threads.restype = ctypes.c_int
    if lib.seams_set_num_threads(1) != 0:
        raise SystemExit("threads")
    return lib


def read_dump(path: Path):
    pos = []
    reading = False
    with path.open() as handle:
        for line in handle:
            if line.startswith("ITEM: ATOMS"):
                reading = True
            elif reading:
                parts = line.split()
                if len(parts) < 5:
                    break
                pos.append([float(parts[2]), float(parts[3]), float(parts[4])])
    return np.asarray(pos, dtype=float)


def chill_file(lib, path: Path):
    system = lib.seams_system_create()
    if lib.seams_system_read(system, str(path).encode(), 1, -1) != 0:
        raise SystemExit(lib.seams_last_error())
    count = lib.seams_system_count(system)
    labels = (ctypes.c_int * count)()
    counts = (ctypes.c_int * 9)()
    if lib.seams_chill_plus(system, 3.5, 0, 0, 0, labels, counts) != 0:
        raise SystemExit(lib.seams_last_error())
    lib.seams_system_destroy(system)
    return np.asarray(list(labels), dtype=int), [int(counts[i]) for i in range(9)]


def layer_period(projection: np.ndarray) -> float:
    """Basal repeat, the same scan as science_render.layer_period."""
    best_score = -1.0
    best = 3.6
    centred = projection - float(np.mean(projection))
    for period in np.linspace(3.30, 4.05, 16):
        wave = np.exp(1j * (2.0 * math.pi / float(period)) * centred)
        score = float(np.abs(np.mean(wave)))
        if score > best_score:
            best_score = score
            best = float(period)
    return best


def layer_score(projection: np.ndarray, period: float) -> float:
    centred = projection - float(np.mean(projection))
    wave = np.exp(1j * (2.0 * math.pi / period) * centred)
    return float(np.abs(np.mean(wave)))


def layer_kinds(pos, chill, axis) -> list[int]:
    """Ih (1) or Ic (0) for each basal layer. Same rule as science_render."""
    chosen = np.flatnonzero(np.isin(chill, [0, 1]))
    axis = np.asarray(axis, dtype=float)
    axis = axis / np.linalg.norm(axis)
    projection = pos[chosen] @ axis
    kinds = chill[chosen]
    bins = np.floor((projection - projection.min()) / layer_period(projection)).astype(int)
    letters = []
    for index in range(int(bins.max()) + 1):
        selected = kinds[bins == index]
        hexagonal = int(np.sum(selected == 1))
        cubic = int(np.sum(selected == 0))
        if hexagonal + cubic == 0:
            continue
        letters.append(1 if hexagonal >= cubic else 0)
    return letters


def cubicity(kinds: list[int]) -> float:
    return sum(1 for kind in kinds if kind == 0) / float(len(kinds))


def strongest_axis(pos, chill):
    chosen = np.flatnonzero(np.isin(chill, [0, 1]))
    best = None
    for name, axis in (
        ("x", np.array([1.0, 0.0, 0.0])),
        ("y", np.array([0.0, 1.0, 0.0])),
        ("z", np.array([0.0, 0.0, 1.0])),
    ):
        projection = pos[chosen] @ axis
        period = layer_period(projection)
        score = layer_score(projection, period)
        if best is None or score > best[0]:
            best = (score, name, axis)
    return best[1], best[2]


def signature(lib, path: Path, spec: str) -> int:
    system = lib.seams_system_create()
    if lib.seams_system_read(system, str(path).encode(), 1, -1) != 0:
        raise SystemExit(lib.seams_last_error())
    cages = ctypes.c_size_t()
    atoms = ctypes.c_size_t()
    if lib.seams_cages_signature(
        system, 3.5, 0, 0, spec.encode(), ctypes.byref(cages), ctypes.byref(atoms)
    ) != 0:
        raise SystemExit(lib.seams_last_error())
    lib.seams_system_destroy(system)
    return int(cages.value)


def read_xyz(path: Path):
    lines = path.read_text().splitlines()
    count = int(lines[0])
    header = lines[1]
    lattice = header.split('Lattice="')[1].split('"')[0].split()
    origin = header.split('origin="')[1].split('"')[0].split()
    box = (
        float(origin[0]),
        float(origin[1]),
        float(origin[2]),
        float(lattice[0]),
        float(lattice[4]),
        float(lattice[8]),
    )
    rows = []
    for line in lines[2 : 2 + count]:
        species, x, y, z, chill = line.split()
        rows.append((species, float(x), float(y), float(z), chill))
    return box, rows


def brine_frame(lib, path: Path):
    box, rows = read_xyz(path)
    type_of = {"O": 1, "H": 2, "Na": 3, "Cl": 4}
    xyz = []
    types = []
    for species, x, y, z, _chill in rows:
        xyz.extend((x, y, z))
        types.append(type_of[species])
    count = len(rows)
    points = (ctypes.c_double * len(xyz))(*xyz)
    box9 = (ctypes.c_double * 9)(box[0], box[1], box[2], box[3], box[4], box[5], 0, 0, 0)
    kinds = (ctypes.c_int * count)(*types)
    system = lib.seams_system_create()
    if lib.seams_system_set_frame(system, points, count, box9, kinds, None) != 0:
        raise SystemExit(lib.seams_last_error())
    labels = (ctypes.c_int * count)()
    counts = (ctypes.c_int * 9)()
    if lib.seams_chill_plus(system, 3.5, 1, 0, 0, labels, counts) != 0:
        raise SystemExit(lib.seams_last_error())
    oxygen = [i for i, row in enumerate(rows) if row[0] == "O"]
    agree = sum(labels[i] == NAMES.index(rows[i][4]) for i in oxygen)
    histogram = {name: 0 for name in NAMES}
    for index in oxygen:
        histogram[NAMES[labels[index]]] += 1
    ice = (ctypes.c_int * count)(
        *[1 if rows[i][0] == "O" and labels[i] in (0, 1) else 0 for i in range(count)]
    )
    front = ctypes.c_void_p()
    # water type 1, sodium type 3, 20 bins, axis z.
    if lib.seams_brine_analyse(system, ice, count, 1, 3, 20, 2, ctypes.byref(front)) != 0:
        raise SystemExit(lib.seams_last_error())
    enter = ctypes.c_double()
    leave = ctypes.c_double()
    thick = ctypes.c_double()
    found = ctypes.c_int()
    if lib.seams_brine_fronts(
        front, ctypes.byref(enter), ctypes.byref(leave), ctypes.byref(thick), ctypes.byref(found)
    ) != 0:
        raise SystemExit("fronts")
    c_ice = ctypes.c_double()
    c_liquid = ctypes.c_double()
    coeff = ctypes.c_double()
    defined = ctypes.c_int()
    if lib.seams_brine_rejection(
        front, ctypes.byref(c_ice), ctypes.byref(c_liquid), ctypes.byref(coeff), ctypes.byref(defined)
    ) != 0:
        raise SystemExit("rejection")
    lib.seams_brine_destroy(front)
    lib.seams_system_destroy(system)
    return {
        "agree": agree,
        "oxygen": len(oxygen),
        "histogram": histogram,
        "found": int(found.value),
        "thickness": float(thick.value),
        "c_ice": float(c_ice.value),
        "k": float(coeff.value),
        "defined": int(defined.value),
    }


def fetch_brine(cache: Path):
    cache.mkdir(parents=True, exist_ok=True)
    paths = []
    for name, digest in BRINE_FILES.items():
        dest = cache / Path(name).name
        if dest.is_file() and hashlib.sha256(dest.read_bytes()).hexdigest() == digest:
            paths.append(dest)
            continue
        try:
            urllib.request.urlretrieve(BRINE_URL + "/" + name, dest)
        except (urllib.error.URLError, TimeoutError, OSError) as exc:
            print(f"{SKIP}: {exc}")
            return None
        got = hashlib.sha256(dest.read_bytes()).hexdigest()
        if got != digest:
            dest.unlink(missing_ok=True)
            print(f"{SKIP}: checksum mismatch for {name}")
            return None
        paths.append(dest)
    return paths


def _bars(title, pairs, svg: Path):
    width = 640
    height = 360
    left = 48
    top = 36
    plot_h = 260
    peak = max(value for _name, value in pairs) or 1
    gap = (width - left - 24) / len(pairs)
    rects = []
    labels = []
    for index, (name, value) in enumerate(pairs):
        bar_h = plot_h * value / peak
        x = left + index * gap + 8
        y = top + plot_h - bar_h
        rects.append(
            f'<rect x="{x:.1f}" y="{y:.1f}" width="{gap - 16:.1f}" height="{bar_h:.1f}" fill="#1f4e79"/>'
        )
        labels.append(
            f'<text x="{x + (gap - 16) / 2:.1f}" y="{top + plot_h + 18}" text-anchor="middle" font-size="12">{name}</text>'
        )
        labels.append(
            f'<text x="{x + (gap - 16) / 2:.1f}" y="{y - 6:.1f}" text-anchor="middle" font-size="12">{value}</text>'
        )
    svg.parent.mkdir(parents=True, exist_ok=True)
    svg.write_text(
        f'<svg xmlns="http://www.w3.org/2000/svg" width="{width}" height="{height}">'
        f'<rect width="100%" height="100%" fill="white"/>'
        f'<text x="24" y="24" font-size="16">{title}</text>'
        + "".join(rects)
        + "".join(labels)
        + "</svg>\n"
    )


def write_figure(stem: str, title: str, pairs, figures: Path):
    _bars(title, pairs, figures / f"{stem}.svg")


def render_interface_v12(pos, chill, box, figures: Path):
    """Paint the slab with the v12 rules from science_render.py.

    Liquid is blended at 0.40. Basal layers are numbered by
    display_layer_ids. The picture has no baked caption. Tachyon is
    used when Vulkan starts; otherwise the same colours and guides are
    written as SVG.
    """
    import science_render as render

    forward, _spread = render.column_axis(pos, chill, box)
    up = np.array([0.0, 0.0, 1.0])
    centre, half_h, half_v, forward, _right, up = render.project(pos, forward, up)
    png = figures / "interface-tachyon.png"
    painted = False
    try:
        render.render_molecule(pos, chill, None, forward, up, png, half_h, half_v, 480, 0.40, False)
        ids = render.display_layer_ids(pos, np.isin(chill, [0, 1]), up)
        render.draw_layer_guides(png, pos, ids, forward, up, half_v)
        painted = png.is_file() and png.stat().st_size > 0
    except Exception as exc:
        print(f"v12 Tachyon render unavailable ({exc}); painting the same blend and guides")
        painted = False
    if painted:
        return
    width, height = 640, 480
    ids = render.display_layer_ids(pos, np.isin(chill, [0, 1]), up)
    vertical = (pos - centre) @ up
    span = 2.0 * half_v

    def xy_of(point):
        offset = point - centre
        x = (float(offset @ _right) / (2.0 * half_h) + 0.5) * (width - 80) + 40
        y = (0.5 - float(offset @ up) / span) * (height - 40) + 20
        return x, y

    shapes = []
    order = np.argsort(np.where(chill == 2, 0, 1))
    for index in order:
        x, y = xy_of(pos[index])
        red, green, blue = render.CHILL.get(int(chill[index]), ("", (0.5, 0.5, 0.5)))[1]
        opacity = 0.40 if int(chill[index]) == 2 else 1.0
        shapes.append(
            f'<circle cx="{x:.1f}" cy="{y:.1f}" r="4.2" fill="rgb({int(red*255)},{int(green*255)},{int(blue*255)})" fill-opacity="{opacity:.2f}"/>'
        )
    guides = []
    for layer in sorted({int(item) for item in ids if item >= 0}):
        members = np.flatnonzero(ids == layer)
        y = (0.5 - float(np.median(vertical[members])) / span) * (height - 40) + 20
        guides.append(f'<line x1="36" y1="{y:.1f}" x2="{width-8}" y2="{y:.1f}" stroke="#a8a8a8" stroke-width="1"/>')
        guides.append(f'<text x="8" y="{y+4:.1f}" font-size="12" fill="#464646">{layer}</text>')
    (figures / "interface.svg").write_text(
        f'<svg xmlns="http://www.w3.org/2000/svg" width="{width}" height="{height}">'
        '<rect width="100%" height="100%" fill="white"/>'
        + "".join(shapes)
        + "".join(guides)
        + "</svg>\n"
    )


def expect(condition, message):
    if not condition:
        raise SystemExit(message)


def run_c(root: Path, build: Path, mode: str, extra: list[str]):
    cc = shutil.which("cc") or shutil.which("gcc")
    if cc is None:
        raise SystemExit("cc is required")
    work = build / "docs-science"
    work.mkdir(parents=True, exist_ok=True)
    binary = work / "science_tutorials"
    src = str(build / "src")
    if not binary.is_file():
        subprocess.check_call(
            [
                cc,
                "-I",
                str(root / "include"),
                "-I",
                src,
                str(root / "docs" / "examples" / "science_tutorials.c"),
                "-L" + src,
                "-lseams",
                "-Wl,-rpath," + src,
                "-lm",
                "-o",
                str(binary),
            ]
        )
    env = os.environ.copy()
    env["LD_LIBRARY_PATH"] = src + os.pathsep + env.get("LD_LIBRARY_PATH", "")
    out = subprocess.check_output([str(binary), mode, *extra], cwd=root, env=env, text=True)
    return out


def _language_env(libdir: Path):
    env = os.environ.copy()
    env["SEAMS_LIB"] = str(libdir / "libseams.so")
    env["LD_LIBRARY_PATH"] = str(libdir) + os.pathsep + env.get("LD_LIBRARY_PATH", "")
    return env


def _expect_languages(text, frames):
    for line in (
        "interface hexagonal 352",
        "HHHHHHHHC",
        "nucleus cubic 47",
        "CHHHCHHHHHHHHH",
        "cages sII 512 16",
        "cages sH 51268 4",
    ):
        expect(line in text, text)
    if frames is not None:
        expect("brine t0 hexagonal 733" in text and "brine tend hexagonal 741" in text, text)
    else:
        expect("brine tutorial skipped" in text, text)


def run_julia(root: Path, libdir: Path, frames):
    julia = shutil.which("julia")
    if julia is None:
        print("julia tutorial skipped: julia was not found")
        return
    args = [julia, str(root / "docs" / "examples" / "science_tutorials.jl"), str(root)]
    if frames is not None:
        args.extend(str(path) for path in frames)
    text = subprocess.check_output(args, cwd=root, env=_language_env(libdir), text=True)
    _expect_languages(text, frames)
    print("julia tutorials ok")


def run_r(root: Path, build: Path, libdir: Path, frames):
    rscript = shutil.which("Rscript")
    if rscript is None:
        print("r tutorial skipped: Rscript was not found")
        return
    cc = shutil.which("cc") or shutil.which("gcc")
    if cc is None:
        print("r tutorial skipped: cc was not found")
        return
    work = build / "docs-science"
    work.mkdir(parents=True, exist_ok=True)
    shared = work / "science_r.so"
    subprocess.check_call(
        [
            cc,
            "-shared",
            "-fPIC",
            "-I",
            str(root / "include"),
            "-I",
            str(libdir),
            str(root / "docs" / "examples" / "science_r.c"),
            "-L" + str(libdir),
            "-lseams",
            "-Wl,-rpath," + str(libdir),
            "-o",
            str(shared),
        ]
    )
    args = [rscript, str(root / "docs" / "examples" / "science_tutorials.R"), str(root), str(shared)]
    if frames is not None:
        args.extend(str(path) for path in frames)
    text = subprocess.check_output(args, cwd=root, env=_language_env(libdir), text=True)
    _expect_languages(text, frames)
    print("r tutorials ok")


def main():
    root = Path(sys.argv[1] if len(sys.argv) > 1 else Path(__file__).resolve().parents[2])
    build = Path(sys.argv[2]) if len(sys.argv) > 2 else Path(os.environ.get("SEAMS_BUILD", root / "bbdir"))
    libdir = Path(os.environ.get("SEAMS_LIBDIR", str(build / "src")))
    lib = _lib(libdir)
    traj = root / "input" / "traj"
    figures = root / "docs" / "source" / "science" / "figures"
    figures.mkdir(parents=True, exist_ok=True)

    interface = traj / "interface_f0.lammpstrj"
    pos = read_dump(interface)
    _labels, counts = chill_file(lib, interface)
    kinds = layer_kinds(pos, _labels, np.array([0.0, 0.0, 1.0]))
    stacking = "".join("C" if kind == 0 else "H" for kind in kinds)
    expect(counts[1] == 352 and counts[2] == 434 and counts[3] == 49 and counts[0] == 1, counts)
    expect(stacking == "HHHHHHHHC", stacking)
    expect(abs(cubicity(kinds) - 1.0 / 9.0) < 1e-12, cubicity(kinds))
    print("interface hexagonal 352 water 434 interfacial 49 cubic 1")
    print("interface layers", stacking, "cubicity", f"{cubicity(kinds):.6f}")

    nucleus = traj / "nucleation_f1.lammpstrj"
    npos = read_dump(nucleus)
    nlabels, ncounts = chill_file(lib, nucleus)
    axis_name, axis = strongest_axis(npos, nlabels)
    nkinds = layer_kinds(npos, nlabels, axis)
    nstack = "".join("C" if kind == 0 else "H" for kind in nkinds)
    expect(ncounts[0] == 47 and ncounts[1] == 62 and ncounts[2] == 3589 and ncounts[3] == 240, ncounts)
    expect(axis_name == "x", axis_name)
    expect(nstack == "CHHHCHHHHHHHHH", nstack)
    expect(abs(cubicity(nkinds) - 2.0 / 14.0) < 1e-12, cubicity(nkinds))
    print("nucleus cubic 47 hexagonal 62 water 3589 interfacial 240")
    print("nucleus axis", axis_name, "layers", nstack, "cubicity", f"{cubicity(nkinds):.6f}")
    edge = []
    with interface.open() as handle:
        for line in handle:
            if line.startswith("ITEM: BOX BOUNDS"):
                for _ in range(3):
                    lo, hi = (float(item) for item in handle.readline().split()[:2])
                    edge.append(hi - lo)
                break
    render_interface_v12(pos, _labels, np.asarray(edge), figures)

    sii = signature(lib, traj / "genice_sII.lammpstrj", "512")
    sii_large = signature(lib, traj / "genice_sII.lammpstrj", "51264")
    sh_small = signature(lib, traj / "genice_sH.lammpstrj", "512")
    sh_medium = signature(lib, traj / "genice_sH.lammpstrj", "sh")
    sh_large = signature(lib, traj / "genice_sH.lammpstrj", "51268")
    expect((sii, sii_large, sh_small, sh_medium, sh_large) == (16, 8, 12, 8, 4), "cages")
    print("cages sII 512 16")
    print("cages sII 51264 8")
    print("cages sH 512 12")
    print("cages sH sh 8")
    print("cages sH 51268 4")
    write_figure(
        "cages",
        "Hydrate cages on the GenIce hosts",
        [("sII 5^12", 16), ("sII 5^12 6^4", 8), ("sH 5^12", 12), ("sH 4^3 5^6 6^3", 8), ("sH 5^12 6^8", 4)],
        figures,
    )

    cache = root / "docs" / ".cache" / "brine-front-230K" / BRINE_REV
    frames = fetch_brine(cache)
    if frames is None:
        _bars("Brine frames were not fetched", [("cached", 0), ("fetched", 0)], figures / "brine.svg")
        print(SKIP)
    else:
        early = brine_frame(lib, frames[0])
        late = brine_frame(lib, frames[1])
        expect(early["agree"] == early["oxygen"] == 1496, early)
        expect(late["agree"] == 1496, late)
        expect(early["histogram"]["hexagonal"] == 733 and early["histogram"]["water"] == 619, early)
        expect(early["histogram"]["cubic"] == 1 and early["histogram"]["interfacial"] == 91, early)
        expect(late["histogram"]["hexagonal"] == 741 and late["histogram"]["water"] == 621, late)
        expect(late["histogram"]["cubic"] == 0 and late["histogram"]["interfacial"] == 89, late)
        expect(early["found"] == 1 and late["found"] == 1, "front")
        expect(early["defined"] == 1 and late["defined"] == 1, "defined")
        expect(early["c_ice"] == 0.0 and late["c_ice"] == 0.0, "c_ice")
        expect(early["k"] == 0.0 and late["k"] == 0.0, "k")
        expect(abs(early["thickness"] - 25.5618) < 5e-3, early["thickness"])
        expect(abs(late["thickness"] - 25.3740) < 5e-3, late["thickness"])
        print("brine t0 hexagonal 733 water 619 interfacial 91 cubic 1 agree 1496")
        print(f"brine t0 thickness {early['thickness']:.4f} k {early['k']:.1f}")
        print("brine tend hexagonal 741 water 621 interfacial 89 cubic 0 agree 1496")
        print(f"brine tend thickness {late['thickness']:.4f} k {late['k']:.1f}")
        write_figure(
            "brine",
            "230 K brine front, 20 NaCl pairs, hexagonal waters",
            [("0 ns", 733), ("5 ns", 741)],
            figures,
        )

    if build.is_dir() and (build / "src" / "libseams.so").is_file():
        c_interface = run_c(root, build, "interface", [])
        expect("interface hexagonal 352" in c_interface, c_interface)
        expect("interface layers HHHHHHHHC" in c_interface, c_interface)
        expect("nucleus cubic 47" in c_interface and "CHHHCHHHHHHHHH" in c_interface, c_interface)
        c_cages = run_c(root, build, "cages", [])
        for line in (
            "cages sII 512 16",
            "cages sII 51264 8",
            "cages sH 512 12",
            "cages sH sh 8",
            "cages sH 51268 4",
        ):
            expect(line in c_cages, c_cages)
        if frames is not None:
            c_brine = run_c(root, build, "brine", [str(frames[0]), str(frames[1])])
            expect("brine t0 hexagonal 733" in c_brine, c_brine)
            expect("brine tend hexagonal 741" in c_brine, c_brine)
            expect("k 0.0" in c_brine, c_brine)
        print("c tutorials ok")
    run_julia(root, libdir, frames)
    run_r(root, build, libdir, frames)
    print("science tutorials ok")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
