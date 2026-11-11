"""CHILL+ on the mixed TIP4P dump agrees across the CLI, Python, R, and Julia.

usage: chill_plus_fronts.py SEAMS LIBSEAMS_LANG [--require-all]

Each front runs ``chill_plus`` on frame 1 of ``input/traj/exampleTraj.lammpstrj``
(oxygen, LAMMPS type 2, cutoff 3.5 A) with its default graph, with
``graph="cutoff"``, and with ``graph="mutual"``. The default must be the
cutoff list. Every installed front returns all three graphs. Histograms
must agree, and per-atom labels must agree across Python, R, and Julia.
A front that is not installed is skipped (exit 77 if any is missing)
unless ``--require-all`` is given.
"""

import inspect
import json
import os
import shutil
import subprocess
import sys
from pathlib import Path

ROOT = Path(os.environ.get("DSEAMS_ROOT", Path(__file__).resolve().parents[1]))
TIP = ROOT / "input" / "traj" / "exampleTraj.lammpstrj"
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
EXPECTED = {
    "cutoff": {"interClathrate": 9, "water": 241},
    "mutual": {"clathrate": 1, "interClathrate": 7, "water": 242},
}
GRAPHS = ("default", "cutoff", "mutual")


class Missing(Exception):
    pass


def histogram(labels):
    counts = {}
    for label in labels:
        name = NAMES[label]
        counts[name] = counts.get(name, 0) + 1
    return counts


def cli(seams):
    out = {}
    for graph in GRAPHS:
        args = [seams, "chill-plus", str(TIP), "--cutoff", "3.5", "--type", "2"]
        if graph != "default":
            args += ["--graph", graph]
        proc = subprocess.run(args, capture_output=True, text=True, check=False)
        if proc.returncode != 0:
            raise RuntimeError(f"seams {graph}: {proc.stdout}{proc.stderr}")
        tokens = proc.stdout.split()
        counts = {}
        for name, value in zip(tokens[::2], tokens[1::2]):
            if name in NAMES and int(value) > 0:
                counts[name] = int(value)
        out[graph] = {"counts": counts, "labels": None}
    return out


def _accepts_graph(method):
    return "graph" in inspect.signature(method).parameters


def _classify(frame, graph):
    """One signature on every front: chill_plus(frame, graph)."""
    if _accepts_graph(frame.chill_plus):
        if graph == "default":
            frame.chill_plus()
        else:
            frame.chill_plus(graph=graph)
        return
    if graph in ("default", "cutoff"):
        frame.chill_plus()
        return
    if graph != "mutual":
        raise RuntimeError(f"unknown graph {graph}")
    # The installed method classifies the cutoff neighbour list.
    # The mutual list is the four nearest neighbours with candidates
    # inside cutoff + 1.5, then the same CHILL+ classifier.
    from pydseams import yoda

    nlist = yoda.kNearestNeighbourList(
        frame.cloud, 4, frame.cutoff + 1.5, frame.atom_type, True
    )
    yoda.getCorrelPlus(frame.cloud, nlist, False)
    yoda.getIceTypePlusNoPrint(frame.cloud, nlist, False)


def python(lib):
    """Python front of this tree: ctypes against libseams_lang.

    pydseamslib on PyPI stops at 2.11 and does not carry the 2015 rule.
    """
    import ctypes

    so = ctypes.CDLL(lib)
    so.seams_frame_read.restype = ctypes.c_void_p
    so.seams_frame_read.argtypes = [ctypes.c_char_p, ctypes.c_int, ctypes.c_int]
    so.seams_frame_free.argtypes = [ctypes.c_void_p]
    so.seams_frame_nop.restype = ctypes.c_int
    so.seams_frame_nop.argtypes = [ctypes.c_void_p]
    so.seams_frame_chill_plus_rule.restype = ctypes.c_int
    so.seams_frame_chill_plus_rule.argtypes = [
        ctypes.c_void_p, ctypes.c_double, ctypes.c_int,
        ctypes.c_char_p, ctypes.c_char_p,
        ctypes.POINTER(ctypes.c_int), ctypes.POINTER(ctypes.c_int),
    ]
    so.seams_last_error.restype = ctypes.c_char_p
    frame = so.seams_frame_read(str(TIP).encode(), 1, 2)
    if not frame:
        err = so.seams_last_error()
        raise RuntimeError(err.decode() if err else "seams_frame_read failed")
    try:
        n = so.seams_frame_nop(frame)
        out = {}
        for graph in GRAPHS:
            labels = (ctypes.c_int * n)()
            counts = (ctypes.c_int * 9)()
            rule = b"chill+"
            graph_arg = b"cutoff" if graph == "default" else graph.encode()
            code = so.seams_frame_chill_plus_rule(
                frame, 3.5, 0, graph_arg, rule, labels, counts
            )
            if code != 0:
                err = so.seams_last_error()
                raise RuntimeError(err.decode() if err else "chill_plus failed")
            got = [labels[i] for i in range(n)]
            out[graph] = {"counts": histogram(got), "labels": got}
        return out
    finally:
        so.seams_frame_free(frame)


def r_front(env):
    if shutil.which("Rscript") is None:
        raise Missing("Rscript is not installed")
    script = f"""
suppressPackageStartupMessages(library(dseams))
frame <- read_frame("{TIP}", type = 2L)
one <- function(ice) as.integer(ice$label)
out <- list(default = one(chill_plus(frame, cutoff = 3.5)),
            cutoff = one(chill_plus(frame, cutoff = 3.5, graph = "cutoff")),
            mutual = one(chill_plus(frame, cutoff = 3.5, graph = "mutual")))
cat(sprintf('{{"default": [%s], "cutoff": [%s], "mutual": [%s]}}\\n',
            paste(out$default, collapse = ","), paste(out$cutoff, collapse = ","),
            paste(out$mutual, collapse = ",")))
"""
    proc = subprocess.run(
        ["Rscript", "-e", script], capture_output=True, text=True, env=env, check=False
    )
    if proc.returncode != 0:
        if "there is no package called" in proc.stderr:
            raise Missing("the dseams R package is not installed")
        raise RuntimeError(f"R: {proc.stderr}")
    data = json.loads(proc.stdout.strip().splitlines()[-1])
    return {g: {"counts": histogram(data[g]), "labels": data[g]} for g in GRAPHS}


def julia_front(env):
    julia = shutil.which("julia")
    if julia is None:
        raise Missing("julia is not installed")
    project = ROOT / "bindings" / "julia" / "Dseams"
    script = f"""
using Dseams
frame = read_frame("{TIP}"; type = 2)
one(ice) = join(Int.(ice.labels), ",")
println("{{\\"default\\": [", one(chill_plus(frame; cutoff = 3.5)),
        "], \\"cutoff\\": [", one(chill_plus(frame; cutoff = 3.5, graph = "cutoff")),
        "], \\"mutual\\": [", one(chill_plus(frame; cutoff = 3.5, graph = "mutual")), "]}}")
"""
    proc = subprocess.run(
        [julia, f"--project={project}", "-e", script],
        capture_output=True,
        text=True,
        env=env,
        check=False,
    )
    if proc.returncode != 0:
        raise RuntimeError(f"Julia: {proc.stdout}{proc.stderr}")
    data = json.loads(proc.stdout.strip().splitlines()[-1])
    return {g: {"counts": histogram(data[g]), "labels": data[g]} for g in GRAPHS}


def main(argv):
    if len(argv) < 3:
        print(__doc__, file=sys.stderr)
        return 2
    seams, lib = argv[1], str(Path(argv[2]).resolve())
    require_all = "--require-all" in argv[3:]
    env = os.environ.copy()
    env["SEAMS_LANG_LIB"] = lib
    lib_dir = str(Path(lib).parent)
    src_dir = str(Path(lib).resolve().parents[1] / "src")
    env["LD_LIBRARY_PATH"] = os.pathsep.join(
        p for p in (lib_dir, src_dir, env.get("LD_LIBRARY_PATH", "")) if p
    )
    env["PATH"] = os.pathsep.join(
        p for p in (os.path.expanduser("~/.juliaup/bin"), env.get("PATH", "")) if p
    )

    results = {"cli": cli(seams)}
    missing = []
    for name, run in (("python", lambda: python(lib)), ("r", lambda: r_front(env)),
                      ("julia", lambda: julia_front(env))):
        try:
            results[name] = run()
        except Missing as exc:
            missing.append(f"{name}: {exc}")

    failures = []
    for front, by_graph in results.items():
        for graph in GRAPHS:
            if graph not in by_graph:
                failures.append(f"{front} has no {graph} result")
                continue
            want = EXPECTED["cutoff" if graph == "default" else graph]
            got = by_graph[graph]["counts"]
            print(f"{front:6s} {graph:7s} {json.dumps(got, sort_keys=True)}")
            if got != want:
                failures.append(f"{front} {graph}: {got} != {want}")
        if "default" in by_graph and "cutoff" in by_graph:
            labels = by_graph["default"]["labels"]
            if labels is not None and labels != by_graph["cutoff"]["labels"]:
                failures.append(f"{front}: default labels are not the cutoff labels")
    for graph in GRAPHS:
        labelled = [
            front
            for front, by_graph in results.items()
            if graph in by_graph and by_graph[graph]["labels"] is not None
        ]
        for front in labelled[1:]:
            if results[front][graph]["labels"] != results[labelled[0]][graph]["labels"]:
                failures.append(
                    f"{front} {graph}: labels differ from {labelled[0]}"
                )

    for line in missing:
        print(f"skipped {line}")
    for line in failures:
        print(f"FAIL {line}", file=sys.stderr)
    if failures:
        return 1
    if missing:
        return 1 if require_all else 77
    return 0


if __name__ == "__main__":
    sys.exit(main(sys.argv))
