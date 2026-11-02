"""CHILL+ on the mixed TIP4P dump agrees across the CLI, Python, R, and Julia.

usage: chill_plus_fronts.py SEAMS LIBSEAMS_LANG [--require-all]

Each front runs ``chill_plus`` on frame 1 of ``input/traj/exampleTraj.lammpstrj``
(oxygen, LAMMPS type 2, cutoff 3.5 A) with its default graph, with
``graph="cutoff"``, and with ``graph="mutual"``. The default must be the
cutoff list. Histograms must agree across all four fronts and per-atom labels
across Python, R, and Julia. A front that is not installed is skipped (exit 77
if any is missing) unless ``--require-all`` is given.
"""

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
    "cutoff": {"interClathrate": 12, "water": 238},
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


def python():
    try:
        import pydseams
    except ImportError as exc:
        raise Missing("pydseamslib is not installed") from exc
    frame = pydseams.read(str(TIP), atom_type=2, cutoff=3.5)
    out = {}
    for graph in GRAPHS:
        try:
            counts = frame.chill_plus() if graph == "default" else frame.chill_plus(graph=graph)
        except TypeError as exc:
            raise RuntimeError("pydseamslib has no chill_plus(graph=...)") from exc
        labels = [NAMES.index(pt.iceType.name) for pt in frame.cloud.pts]
        out[graph] = {"counts": {k: v for k, v in counts.items() if v}, "labels": labels}
    return out


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
    env["LD_LIBRARY_PATH"] = os.pathsep.join(
        p for p in (str(Path(lib).parent), env.get("LD_LIBRARY_PATH", "")) if p
    )

    results = {"cli": cli(seams)}
    missing = []
    for name, run in (("python", python), ("r", lambda: r_front(env)),
                      ("julia", lambda: julia_front(env))):
        try:
            results[name] = run()
        except Missing as exc:
            missing.append(f"{name}: {exc}")

    failures = []
    for front, by_graph in results.items():
        for graph in GRAPHS:
            want = EXPECTED["cutoff" if graph == "default" else graph]
            got = by_graph[graph]["counts"]
            print(f"{front:6s} {graph:7s} {json.dumps(got, sort_keys=True)}")
            if got != want:
                failures.append(f"{front} {graph}: {got} != {want}")
        labels = by_graph["default"]["labels"]
        if labels is not None and labels != by_graph["cutoff"]["labels"]:
            failures.append(f"{front}: default labels are not the cutoff labels")
    labelled = [f for f in results if results[f]["default"]["labels"] is not None]
    for graph in GRAPHS:
        for front in labelled[1:]:
            if results[front][graph]["labels"] != results[labelled[0]][graph]["labels"]:
                failures.append(f"{front} {graph}: labels differ from {labelled[0]}")

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
