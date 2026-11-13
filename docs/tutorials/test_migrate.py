#!/usr/bin/env python3
# Copyright (c) d-SEAMS developers
"""Run each 1.x mapping on one input and check the 2.14 sides agree."""

import os
import subprocess
import sys
from ctypes import (
    CDLL,
    POINTER,
    Structure,
    byref,
    c_char_p,
    c_double,
    c_int,
    c_void_p,
)
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
from migrate_map import PAIRS

LABELS = (
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


class Prism(Structure):
    _fields_ = (("size", c_int), ("perfect", c_int), ("imperfect", c_int))


def fail(message):
    print(message, file=sys.stderr)
    raise SystemExit(1)


def run_cli(seams, args):
    env = os.environ.copy()
    env["NO_COLOR"] = "1"
    done = subprocess.run(
        [seams, *args],
        check=False,
        capture_output=True,
        text=True,
        env=env,
    )
    if done.returncode != 0:
        fail("seams {} failed\n{}".format(" ".join(args), done.stderr[-400:]))
    return done.stdout


def counts_from_cli(text):
    found = {}
    for line in text.splitlines():
        parts = line.split()
        if parts and parts[0] == "nop":
            for index in range(2, len(parts) - 1, 2):
                if parts[index + 1].isdigit() or (
                    parts[index + 1].startswith("-") and parts[index + 1][1:].isdigit()
                ):
                    found[parts[index]] = int(parts[index + 1])
    return found


def load(libdir, langdir):
    os.environ["LD_LIBRARY_PATH"] = os.pathsep.join(
        [str(libdir), str(langdir), os.environ.get("LD_LIBRARY_PATH", "")]
    )
    api = CDLL(str(Path(libdir) / "libseams.so"))
    lang = CDLL(str(Path(langdir) / "libseams_lang.so"))
    lang.seams_frame_read.restype = c_void_p
    lang.seams_frame_read.argtypes = [c_char_p, c_int, c_int]
    lang.seams_frame_free.argtypes = [c_void_p]
    lang.seams_frame_nop.restype = c_int
    lang.seams_frame_nop.argtypes = [c_void_p]
    lang.seams_frame_positions.restype = POINTER(c_double)
    lang.seams_frame_positions.argtypes = [c_void_p]
    lang.seams_frame_types.restype = POINTER(c_int)
    lang.seams_frame_types.argtypes = [c_void_p]
    lang.seams_frame_box.argtypes = [c_void_p, POINTER(c_double)]
    for name in ("seams_frame_chill_plus", "seams_frame_chill"):
        fn = getattr(lang, name)
        fn.restype = c_int
        fn.argtypes = [c_void_p, c_double, c_int, c_char_p, POINTER(c_int), POINTER(c_int)]
    lang.seams_frame_rings.restype = c_int
    lang.seams_frame_rings.argtypes = [
        c_void_p, c_double, c_int, c_char_p, c_int, POINTER(c_int), c_int,
    ]
    lang.seams_frame_steinhardt_graph.restype = c_int
    lang.seams_frame_steinhardt_graph.argtypes = [
        c_void_p, c_double, c_int, c_int, c_char_p, POINTER(c_double), POINTER(c_double),
    ]
    lang.seams_frame_cages.restype = c_int
    lang.seams_frame_cages.argtypes = [
        c_void_p, c_double, c_int, c_int, c_char_p, c_int,
        POINTER(c_int), POINTER(c_int), POINTER(c_int), POINTER(c_int), POINTER(c_int),
    ]
    lang.seams_frame_signature_edge.restype = c_int
    lang.seams_frame_signature_edge.argtypes = [
        c_void_p, c_double, c_int, c_char_p, c_char_p, c_char_p, POINTER(c_int), POINTER(c_int),
    ]
    lang.seams_frame_neighbours.restype = c_int
    lang.seams_frame_neighbours.argtypes = [c_void_p, c_double, c_int, c_char_p]
    lang.seams_frame_neighbour_pairs.restype = c_int
    lang.seams_frame_neighbour_pairs.argtypes = [c_void_p]
    lang.seams_frame_rdf_nbins.restype = c_int
    lang.seams_frame_rdf_nbins.argtypes = [c_double, c_int]
    lang.seams_frame_rdf.restype = c_int
    lang.seams_frame_rdf.argtypes = [
        c_void_p, c_double, c_int, c_int, c_int,
        POINTER(c_double), POINTER(c_double), POINTER(c_int), POINTER(c_double),
    ]
    lang.seams_frame_prisms.restype = c_int
    lang.seams_frame_prisms.argtypes = [
        c_void_p, c_double, c_int, c_int, c_char_p, POINTER(Prism), POINTER(c_int),
    ]
    lang.seams_frame_monolayer.restype = c_int
    lang.seams_frame_monolayer.argtypes = [
        c_void_p, c_double, c_int, c_int, c_double, c_char_p, POINTER(c_int), c_int,
    ]
    api.seams_system_create.restype = c_void_p
    api.seams_system_destroy.argtypes = [c_void_p]
    api.seams_system_set_frame.restype = c_int
    api.seams_system_set_frame.argtypes = [
        c_void_p, POINTER(c_double), c_int, POINTER(c_double), POINTER(c_int), POINTER(c_int),
    ]
    api.seams_neighbours_build.restype = c_int
    api.seams_neighbours_build.argtypes = [c_void_p, c_double, c_int, c_int, POINTER(c_void_p)]
    api.seams_neighbours_destroy.argtypes = [c_void_p]
    api.seams_neighbours_csr.restype = c_int
    api.seams_neighbours_csr.argtypes = [
        c_void_p, POINTER(POINTER(c_int)), POINTER(c_int), POINTER(POINTER(c_int)), POINTER(c_int),
    ]
    return api, lang


def open_frame(lang, path, type_filter):
    frame = lang.seams_frame_read(str(path).encode(), 1, type_filter)
    if not frame:
        fail("read failed {}".format(path))
    return frame


def chill_counts(lang, frame, which):
    nop = lang.seams_frame_nop(frame)
    labels = (c_int * nop)()
    counts = (c_int * 9)()
    fn = lang.seams_frame_chill_plus if which == "plus" else lang.seams_frame_chill
    if fn(frame, 3.5, 1, b"cutoff", labels, counts) != 0:
        fail("chill call failed")
    return {LABELS[index]: counts[index] for index in range(9) if counts[index]}


def csr_pairs(api, lang, frame, graph_name, graph_id):
    if lang.seams_frame_neighbours(frame, 3.5, 1, graph_name) != 0:
        fail("frame neighbours failed")
    directed = lang.seams_frame_neighbour_pairs(frame)
    nop = lang.seams_frame_nop(frame)
    xyz = lang.seams_frame_positions(frame)
    types = lang.seams_frame_types(frame)
    box = (c_double * 9)()
    lang.seams_frame_box(frame, box)
    system = api.seams_system_create()
    if api.seams_system_set_frame(system, xyz, nop, box, types, None) != 0:
        fail("set_frame failed")
    handle = c_void_p()
    if api.seams_neighbours_build(system, 3.5, 0, graph_id, byref(handle)) != 0:
        fail("neighbours_build failed")
    offsets = POINTER(c_int)()
    cols = POINTER(c_int)()
    n_off = c_int()
    nnz = c_int()
    if api.seams_neighbours_csr(handle, byref(offsets), byref(n_off), byref(cols), byref(nnz)) != 0:
        fail("csr failed")
    api.seams_neighbours_destroy(handle)
    api.seams_system_destroy(system)
    if directed != nnz.value:
        fail("neighbour pairs {} and {}".format(directed, nnz.value))
    return directed


def expect_line(text, needle):
    if needle not in text:
        fail("missing {!r} in\n{}".format(needle, text))


def main():
    root, seams, libdir, langdir = sys.argv[1:5]
    root = Path(root)
    page = root / "docs" / "orgmode" / "howto" / "migrate-1x.org"
    page_text = page.read_text()
    for pair in PAIRS:
        for name in (pair["old"], *pair["also"], pair["cli"], pair["c_api"]):
            if name not in page_text:
                fail("migration page omits {}".format(name))
    api, lang = load(libdir, langdir)
    cubic = root / "input" / "traj" / "mW_cubic.lammpstrj"
    tip = root / "input" / "traj" / "exampleTraj.lammpstrj"
    guest = root / "input" / "traj" / "genice_sI.lammpstrj"
    ran = set()

    def once(kind, fn):
        if kind in ran:
            return
        fn()
        ran.add(kind)
        print(kind)

    once("read", lambda: check_read(seams, lang, cubic))
    once("chill_plus", lambda: check_chill(seams, lang, cubic, "plus", "chill-plus"))
    once("chill", lambda: check_chill(seams, lang, cubic, "plain", "chill"))
    once("rings", lambda: check_rings(seams, lang, cubic))
    once("steinhardt", lambda: check_stein(seams, lang, cubic))
    once("cages", lambda: check_cages(seams, lang, cubic))
    once("signature", lambda: check_signature(seams, lang, guest))
    once("neighbors", lambda: check_graph(api, lang, cubic, b"cutoff", 0))
    once("knn", lambda: check_graph(api, lang, cubic, b"mutual", 1))
    once("rdf", lambda: check_rdf(seams, lang, cubic))
    once("prisms", lambda: check_prisms(seams, lang, cubic))
    once("monolayer", lambda: check_monolayer(lang, tip))
    once("hbonds", lambda: expect_line(
        run_cli(seams, ["hbonds", str(tip), "--type", "2", "--htype", "1"]), "hbonds 460"))
    once("cn", lambda: expect_line(
        run_cli(seams, ["cn", str(tip), "--pair-types", "1,2", "--cutoff", "4.5"]), "cn 9.144"))
    once("density", lambda: expect_line(
        run_cli(seams, ["density-z", str(tip), "--type", "2", "--bins", "4"]), "0.000887889"))
    once("fingerprint", lambda: expect_line(
        run_cli(seams, ["fingerprint", str(tip), "--type", "2"]), "classes 89"))
    once("pairs", lambda: expect_line(
        run_cli(seams, ["pairs", str(tip), "--site", "1=cationHead,2=anion"]), "count 0"))
    once("domains", lambda: expect_line(
        run_cli(seams, ["domains", str(tip), "--site", "1=cationHead,2=anion", "--subset", "polar"]),
        "largest 750"))
    missing = {pair["kind"] for pair in PAIRS} - ran
    if missing:
        fail("unrun mappings {}".format(sorted(missing)))
    print("migrate {} mappings".format(len(PAIRS)))


def check_read(seams, lang, path):
    text = run_cli(seams, ["read", str(path), "--type", "1"])
    frame = open_frame(lang, path, 1)
    try:
        nop = lang.seams_frame_nop(frame)
    finally:
        lang.seams_frame_free(frame)
    line = text.splitlines()[0].split()
    if line[0] != "nop" or int(line[1]) != nop:
        fail("read {} vs {}".format(line, nop))


def check_chill(seams, lang, path, which, command):
    args = [command, str(path), "--type", "1", "--cutoff", "3.5", "--graph", "cutoff"]
    if command == "chill-plus":
        args.extend(["--rule", "published"])
    text = run_cli(seams, args)
    frame = open_frame(lang, path, 1)
    try:
        got = chill_counts(lang, frame, "plus" if which == "plus" else "plain")
    finally:
        lang.seams_frame_free(frame)
    cli = counts_from_cli(text)
    if cli != got:
        fail("{} cli {} c {}".format(command, cli, got))


def check_rings(seams, lang, path):
    text = run_cli(seams, ["rings", str(path), "--type", "1", "--cutoff", "3.5", "--max", "6"])
    hist = {}
    for line in text.splitlines():
        parts = line.split()
        if len(parts) == 2 and parts[0].isdigit():
            hist[int(parts[0])] = int(parts[1])
    frame = open_frame(lang, path, 1)
    try:
        counts = (c_int * 7)()
        if lang.seams_frame_rings(frame, 3.5, 1, b"cutoff", 6, counts, 7) != 0:
            fail("rings call")
    finally:
        lang.seams_frame_free(frame)
    got = {index: counts[index] for index in range(7) if counts[index]}
    if hist != got:
        fail("rings cli {} c {}".format(hist, got))


def check_stein(seams, lang, path):
    text = run_cli(seams, ["steinhardt", str(path), "--type", "1", "--cutoff", "3.5", "--l", "6"])
    q6 = q6bar = None
    for line in text.splitlines():
        if " q6 " in line:
            q6 = float(line.split(" q6 ", 1)[1].split()[0])
        if " q6bar " in line:
            q6bar = float(line.split(" q6bar ", 1)[1].split()[0])
    frame = open_frame(lang, path, 1)
    try:
        nop = lang.seams_frame_nop(frame)
        local = (c_double * nop)()
        bar = (c_double * nop)()
        if lang.seams_frame_steinhardt_graph(frame, 3.5, 1, 6, b"cutoff", local, bar) != 0:
            fail("steinhardt call")
    finally:
        lang.seams_frame_free(frame)
    mean = sum(local) / nop
    mean_bar = sum(bar) / nop
    if abs(mean - q6) > 5e-5 or abs(mean_bar - q6bar) > 5e-5:
        fail("q6 cli {} {} c {} {}".format(q6, q6bar, mean, mean_bar))


def check_cages(seams, lang, path):
    text = run_cli(seams, ["cages", str(path), "--type", "1", "--cutoff", "3.5", "--graph", "seeded"])
    cli = counts_from_cli(text)
    frame = open_frame(lang, path, 1)
    try:
        n_hex = c_int()
        n_cubic = c_int()
        n_water = c_int()
        n_both = c_int()
        if lang.seams_frame_cages(
            frame, 3.5, 1, 4, b"seeded", 0, None, byref(n_hex), byref(n_cubic), byref(n_water), byref(n_both)
        ) != 0:
            fail("cages call")
    finally:
        lang.seams_frame_free(frame)
    got = {
        "hexagonal": n_hex.value,
        "cubic": n_cubic.value,
        "water": n_water.value,
        "both": n_both.value,
    }
    for name, value in got.items():
        if cli.get(name, 0) != value:
            fail("cages {} cli {} c {}".format(name, cli, got))


def check_signature(seams, lang, path):
    text = run_cli(
        seams,
        ["cages", str(path), "--type", "1", "--cutoff", "3.5", "--graph", "cutoff",
         "--signature", "512", "--cage-edge", "none"],
    )
    cages = atoms = None
    for line in text.splitlines():
        if " cages " in line:
            cages = int(line.split(" cages ", 1)[1].split()[0])
        if " atoms " in line:
            atoms = int(line.split(" atoms ", 1)[1].split()[0])
    frame = open_frame(lang, path, 1)
    try:
        n_cages = c_int()
        n_atoms = c_int()
        if lang.seams_frame_signature_edge(
            frame, 3.5, 1, b"cutoff", b"512", b"none", byref(n_cages), byref(n_atoms)
        ) != 0:
            fail("signature call")
    finally:
        lang.seams_frame_free(frame)
    if cages != n_cages.value or atoms != n_atoms.value:
        fail("signature cli {} {} c {} {}".format(cages, atoms, n_cages.value, n_atoms.value))


def check_graph(api, lang, path, name, graph_id):
    frame = open_frame(lang, path, 1)
    try:
        pairs = csr_pairs(api, lang, frame, name, graph_id)
    finally:
        lang.seams_frame_free(frame)
    if pairs <= 0:
        fail("no pairs for {}".format(name))


def check_rdf(seams, lang, path):
    text = run_cli(
        seams, ["rdf", str(path), "--type", "1", "--pair-types", "1,1", "--cutoff", "6", "--bins", "4"]
    )
    rows = []
    for line in text.splitlines():
        parts = line.split()
        if len(parts) == 3 and parts[0][0].isdigit():
            rows.append((float(parts[1]), int(parts[2])))
    frame = open_frame(lang, path, 1)
    try:
        width = lang.seams_frame_rdf_nbins(6.0, 4)
        radius = (c_double * width)()
        gram = (c_double * width)()
        count = (c_int * width)()
        volume = c_double()
        if lang.seams_frame_rdf(frame, 6.0, 4, 1, 1, radius, gram, count, byref(volume)) != 0:
            fail("rdf call")
    finally:
        lang.seams_frame_free(frame)
    if len(rows) != width:
        fail("rdf rows {} width {}".format(len(rows), width))
    for index, (gee, number) in enumerate(rows):
        if abs(gee - gram[index]) > 1e-4 or number != count[index]:
            fail("rdf bin {} cli {} c {} {}".format(index, rows[index], gram[index], count[index]))


def check_prisms(seams, lang, path):
    text = run_cli(seams, ["prisms", str(path), "--type", "1", "--cutoff", "3.5"])
    cli = {}
    for line in text.splitlines():
        parts = line.split()
        if len(parts) == 6 and parts[0] == "size":
            cli[int(parts[1])] = (int(parts[3]), int(parts[5]))
    frame = open_frame(lang, path, 1)
    try:
        rows = (Prism * 4)()
        n_rows = c_int(4)
        if lang.seams_frame_prisms(frame, 3.5, 1, 6, None, rows, byref(n_rows)) != 0:
            fail("prisms call")
    finally:
        lang.seams_frame_free(frame)
    got = {rows[index].size: (rows[index].perfect, rows[index].imperfect) for index in range(n_rows.value)}
    if cli != got:
        fail("prisms cli {} c {}".format(cli, got))


def check_monolayer(lang, path):
    frame = open_frame(lang, path, 2)
    try:
        first = (c_int * 8)()
        second = (c_int * 8)()
        if lang.seams_frame_monolayer(frame, 3.5, 2, 6, 1.0, None, first, 8) != 0:
            fail("monolayer first")
        if lang.seams_frame_monolayer(frame, 3.5, 2, 6, 1.0, None, second, 8) != 0:
            fail("monolayer second")
    finally:
        lang.seams_frame_free(frame)
    if [first[index] for index in range(8)] != [second[index] for index in range(8)]:
        fail("monolayer calls differ")


if __name__ == "__main__":
    main()
