#!/usr/bin/env python3
"""Regenerate the 230 K ice/brine front frames from the campaign deck.

The deck is the one that wrote results/brine/T230_m{0,20,40}_r0 in the
d-SEAMS reproducibility package: GenIce ice Ih, TIP4P geometry, replication
4 4 6 (1536 waters), the upper half marked liquid, Madrid 2019 ions drawn
with replica r0 (numpy seed 0). GenIce's own default seed is 1000; the
deck does not override it. The LAMMPS input is deck/in.brine.lmp: 100 ps
melt at 400 K with the ice held, 100 ps quench to 230 K, then NPT along z
at 1 bar with the campaign thermostat (Tdamp 200 fs, Pdamp 2000 fs) and the
velocity seed "{replica}3" (03 for r0).

Dumps and restarts stay out of git. The production deck writes
step-0.restart, step-830000.restart, step-1670000.restart and
step-2500000.restart. frames reads those restarts only: in.resume.lmp
with run 0 dumps the coordinates, and d-SEAMS writes the CHILL+ label.
The extended-XYZ files land in out/, which is gitignored.

    python3 paper/data/brine-frames/regenerate.py prepare --work work/brine
    python3 paper/data/brine-frames/regenerate.py run --work work/brine --engine omp --threads 2
    python3 paper/data/brine-frames/regenerate.py resume --work work/brine --lmp lmp --threads 2
    python3 paper/data/brine-frames/regenerate.py frames --restarts data/restarts --lmp lmp --seams seams --lib libyodaLib.so
"""

from __future__ import annotations

import argparse
import ctypes
import math
import os
import shutil
import subprocess
import sys
from pathlib import Path

HERE = Path(__file__).resolve().parent
DECK = HERE / "deck"
OUT = HERE / "out"
PAIRS = (0, 20, 40)
# t0, one third, two thirds, and 5 ns. 830000 and 1670000 are the frames
# choose_frames used to pick from a dump saved every 10000 steps.
FRAME_STEPS = (("t0", 0), ("t1", 830_000), ("t2", 1_670_000), ("tend", 2_500_000))
TEMPERATURE = 230.0
REPLICA = 0
PRESSURE = 1.0
MELT_STEPS = 50000
PROD_STEPS = 2500000  # 5 ns at 2 fs; the wall clock may stop it earlier
DUMP_EVERY = 10000  # 20 ps
# The LAMMPS timer caps each run command. Melt and quench finish in minutes.
# 5 h 40 min of production plus that melt is the 6 hour wall.
TIMEOUT = "5:40:00"
SYMBOL = {1: "O", 2: "H", 3: "Na", 4: "Cl"}
CHILL = (
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
# Campaign table, 230 K, replica r0, from results/brine/summary.txt.
# nmax is the last largest cage cluster; ionliq is the last ions-in-liquid count.
CAMPAIGN = {
    0: {"nmax": 1115, "ionliq": 0, "chillice": 880},
    20: {"nmax": 997, "ionliq": 25, "chillice": 785},
    40: {"nmax": 859, "ionliq": 63, "chillice": 694},
}


def velocity_seed(replica: int) -> str:
    """The deck passes -var seed {replica}3, so replica 0 is 03."""
    return f"{replica}3"


def which(name: str) -> str:
    found = shutil.which(name)
    if not found:
        raise SystemExit(f"{name} not on PATH")
    return found


def prepare(work: Path, genice: str) -> None:
    work.mkdir(parents=True, exist_ok=True)
    script = DECK / "brine_system.py"
    for pairs in PAIRS:
        tag = f"T230_m{pairs}_r0"
        dest = work / tag
        dest.mkdir(parents=True, exist_ok=True)
        subprocess.run(
            [
                sys.executable,
                str(script),
                "--genice",
                genice,
                "--rep",
                "4",
                "4",
                "6",
                "--pairs",
                str(pairs),
                "--seed",
                str(REPLICA),
                "--data",
                str(dest / "brine.data"),
                "--groups",
                str(dest / "groups.lmp"),
            ],
            check=True,
        )


def run_cmd(engine: str, threads: int, cores: str, lmp: str, mpirun: str) -> list[str]:
    if engine == "omp":
        cmd = ["taskset", "-c", cores, lmp, "-sf", "omp", "-pk", f"omp {threads}"]
        env_threads = str(threads)
    elif engine == "mpi":
        cmd = [
            mpirun,
            "--bind-to",
            "core",
            "--map-by",
            "core",
            "-np",
            str(threads),
            "--cpu-set",
            cores,
            lmp,
        ]
        env_threads = "1"
    else:
        raise SystemExit(f"unknown engine {engine}")
    return cmd, env_threads


def launch(work: Path, engine: str, threads: int, lmp: str, mpirun: str) -> list[subprocess.Popen]:
    procs = []
    # Eight cores, three loadings. Two cores each; the last pair stays idle
    # so the three runs do not steal from each other.
    core_sets = ("0,1", "2,3", "4,5")
    for pairs, cores in zip(PAIRS, core_sets):
        tag = f"T230_m{pairs}_r0"
        dest = work / tag
        if not (dest / "brine.data").is_file():
            raise SystemExit(f"missing {dest / 'brine.data'}; run prepare")
        cmd, env_threads = run_cmd(engine, threads, cores, lmp, mpirun)
        env = os.environ.copy()
        env["OMP_NUM_THREADS"] = env_threads
        env["OMP_PROC_BIND"] = "close"
        full = cmd + [
            "-log",
            "log.lammps",
            "-in",
            str(DECK / "in.brine.lmp"),
            "-var",
            "data",
            "brine.data",
            "-var",
            "groups",
            "groups.lmp",
            "-var",
            "T",
            f"{TEMPERATURE:.2f}",
            "-var",
            "P",
            f"{PRESSURE:.1f}",
            "-var",
            "seed",
            velocity_seed(REPLICA),
            "-var",
            "steps_melt",
            str(MELT_STEPS),
            "-var",
            "steps",
            str(PROD_STEPS),
            "-var",
            "dumpevery",
            str(DUMP_EVERY),
            "-var",
            "dump",
            "traj.lammpstrj",
            "-var",
            "timeout",
            TIMEOUT,
        ]
        log = open(dest / "driver.log", "w")
        print("start", tag, "cores", cores, " ".join(full), flush=True)
        procs.append(subprocess.Popen(full, cwd=dest, env=env, stdout=log, stderr=subprocess.STDOUT))
    return procs


def parse_dump(path: Path):
    """Yield (step, box_origin, box_length, atoms) for each frame.

    atoms is a list of (id, type, x, y, z). box lengths are xhi-xlo.
    """
    lines = path.read_text().splitlines()
    i = 0
    frames = []
    while i < len(lines):
        if lines[i] != "ITEM: TIMESTEP":
            i += 1
            continue
        step = int(lines[i + 1])
        # NUMBER OF ATOMS
        n = int(lines[i + 3])
        # BOX BOUNDS
        bounds = []
        row = i + 5
        # skip the ITEM: BOX BOUNDS line, which is i+4
        for k in range(3):
            parts = lines[row + k].split()
            lo, hi = float(parts[0]), float(parts[1])
            bounds.append((lo, hi))
        atom_header = row + 3
        if not lines[atom_header].startswith("ITEM: ATOMS"):
            raise SystemExit(f"unexpected dump header at {path}:{atom_header}")
        atoms = []
        for line in lines[atom_header + 1 : atom_header + 1 + n]:
            parts = line.split()
            atoms.append((int(parts[0]), int(parts[1]), float(parts[2]), float(parts[3]), float(parts[4])))
        origin = [bounds[k][0] for k in range(3)]
        length = [bounds[k][1] - bounds[k][0] for k in range(3)]
        frames.append((step, origin, length, atoms))
        i = atom_header + 1 + n
    return frames


def load_chill(lib_path: str):
    lib = ctypes.CDLL(lib_path)
    # box is ox, oy, oz, lx, ly, lz. Labels are for the oxygen coordinates only.
    lib.seams_chill_plus_graph.argtypes = [
        ctypes.POINTER(ctypes.c_double),
        ctypes.c_int,
        ctypes.POINTER(ctypes.c_double),
        ctypes.c_double,
        ctypes.c_char_p,
        ctypes.POINTER(ctypes.c_int),
    ]
    lib.seams_chill_plus_graph.restype = ctypes.c_int
    return lib


def chill_labels(lib, atoms, origin, length):
    oxygen_at = [i for i, atom in enumerate(atoms) if atom[1] == 1]
    n = len(oxygen_at)
    xyz = (ctypes.c_double * (3 * n))()
    for k, i in enumerate(oxygen_at):
        _aid, _typ, x, y, z = atoms[i]
        xyz[3 * k] = x
        xyz[3 * k + 1] = y
        xyz[3 * k + 2] = z
    box = (ctypes.c_double * 6)(origin[0], origin[1], origin[2], length[0], length[1], length[2])
    labels_o = (ctypes.c_int * n)()
    rc = lib.seams_chill_plus_graph(xyz, n, box, 3.5, b"cutoff", labels_o)
    if rc != 0:
        raise SystemExit(f"CHILL+ failed ({rc})")
    labels = [6] * len(atoms)
    for k, i in enumerate(oxygen_at):
        labels[i] = int(labels_o[k])
    return labels


def cli_chill_counts(seams: str, dump: Path, frame: int | None = None) -> dict[str, int]:
    cmd = [seams]
    if frame is not None:
        cmd.extend(["--frame", str(frame)])
    cmd.extend(["chill-plus", str(dump), "--type", "1", "--cutoff", "3.5"])
    proc = subprocess.run(
        cmd,
        check=False,
        capture_output=True,
        text=True,
    )
    if proc.returncode != 0:
        raise SystemExit(proc.stdout + proc.stderr)
    counts = {}
    parts = proc.stdout.split()
    for i, word in enumerate(parts):
        if word in CHILL and i + 1 < len(parts) and parts[i + 1].isdigit():
            counts[word] = int(parts[i + 1])
    if not counts:
        raise SystemExit(f"chill-plus returned no labels: {proc.stdout}")
    return counts


def write_xyz(path: Path, step, origin, length, atoms, labels):
    # 2 fs timestep. Step 0 is the released front.
    t_ps = step * 0.002
    lattice = (
        f"{length[0]:.6f} 0 0 0 {length[1]:.6f} 0 0 0 {length[2]:.6f}"
    )
    lines = [
        str(len(atoms)),
        (
            f"t_ps={t_ps:.3f} step={step} "
            f'Lattice="{lattice}" origin="{origin[0]:.6f} {origin[1]:.6f} {origin[2]:.6f}" '
            "Properties=species:S:1:pos:R:3:chill:S:1"
        ),
    ]
    order = sorted(range(len(atoms)), key=lambda i: atoms[i][0])
    for i in order:
        _aid, typ, x, y, z = atoms[i]
        name = SYMBOL[typ]
        if typ == 1:
            code = labels[i]
            label = CHILL[code] if 0 <= code < len(CHILL) else str(code)
        else:
            label = "-"
        lines.append(f"{name} {x:.6f} {y:.6f} {z:.6f} {label}")
    path.write_text("\n".join(lines) + "\n")


def mic(a, b, length):
    acc = 0.0
    for k in range(3):
        d = b[k] - a[k]
        d -= length[k] * round(d / length[k])
        acc += d * d
    return math.sqrt(acc)


def largest_cluster(points, length, flags, k=4, candidate=5.5):
    """Largest component of flagged points on the union of k nearest neighbours.

    Neighbours are the k nearest oxygens inside `candidate` angstrom, the
    same graph the campaign uses for the cage cluster (K=4, candidate 5.5).
    An edge joins two flagged points when either lists the other.
    """
    n = len(points)
    if n == 0 or not any(flags):
        return 0
    neigh = [[] for _ in range(n)]
    for i in range(n):
        best = []
        for j in range(n):
            if i == j:
                continue
            d = mic(points[i], points[j], length)
            if d <= candidate:
                best.append((d, j))
        best.sort()
        neigh[i] = [j for _d, j in best[:k]]
    parent = list(range(n))

    def find(a):
        while parent[a] != a:
            parent[a] = parent[parent[a]]
            a = parent[a]
        return a

    def union(a, b):
        a, b = find(a), find(b)
        if a != b:
            parent[b] = a

    for i, row in enumerate(neigh):
        for j in row:
            if j > i and flags[i] and flags[j]:
                union(i, j)
    counts = {}
    for i, flag in enumerate(flags):
        if flag:
            root = find(i)
            counts[root] = counts.get(root, 0) + 1
    return max(counts.values()) if counts else 0


def ions_in_liquid(oxygens, ice_flags, ions, length, cutoff=3.5):
    """Campaign rule: no cage oxygen in the 3.5 A shell means liquid."""
    n_ice = n_front = n_liq = 0
    cut2 = cutoff * cutoff
    for ion in ions:
        labelled = 0
        in_shell = 0
        for oxy, flag in zip(oxygens, ice_flags):
            acc = 0.0
            for k in range(3):
                d = oxy[k] - ion[k]
                d -= length[k] * round(d / length[k])
                acc += d * d
            if acc < cut2:
                in_shell += 1
                labelled += 1 if flag else 0
        if in_shell > 0 and labelled == in_shell:
            n_ice += 1
        elif labelled > 0:
            n_front += 1
        else:
            n_liq += 1
    return n_ice, n_front, n_liq


def write_one_frame_dump(path: Path, step, origin, length, atoms):
    lines = [
        "ITEM: TIMESTEP",
        str(step),
        "ITEM: NUMBER OF ATOMS",
        str(len(atoms)),
        "ITEM: BOX BOUNDS pp pp pp",
        f"{origin[0]:.8f} {origin[0] + length[0]:.8f}",
        f"{origin[1]:.8f} {origin[1] + length[1]:.8f}",
        f"{origin[2]:.8f} {origin[2] + length[2]:.8f}",
        "ITEM: ATOMS id mol type x y z",
    ]
    # One molecule per oxygen so the water selector keeps every type-1 atom
    # with its hydrogens. Ions are their own molecules.
    mol = 0
    # atoms are not grouped; assign mol by runs of O,H,H as written by the
    # data file. Fall back to one molecule per atom if the pattern breaks.
    i = 0
    ordered = sorted(atoms, key=lambda a: a[0])
    mol_of = {}
    while i < len(ordered):
        mol += 1
        _aid, typ, *_rest = ordered[i]
        if typ == 1 and i + 2 < len(ordered) and ordered[i + 1][1] == 2 and ordered[i + 2][1] == 2:
            for j in range(3):
                mol_of[ordered[i + j][0]] = mol
            i += 3
        else:
            mol_of[ordered[i][0]] = mol
            i += 1
    for aid, typ, x, y, z in ordered:
        lines.append(f"{aid} {mol_of[aid]} {typ} {x:.6f} {y:.6f} {z:.6f}")
    path.write_text("\n".join(lines) + "\n")


def cage_flags(seams: str, dump: Path, per_atom: Path) -> dict[int, int]:
    if per_atom.exists():
        per_atom.unlink()
    proc = subprocess.run(
        [
            seams,
            "cages",
            str(dump),
            "--graph",
            "seeded",
            "--complete",
            "--type",
            "1",
            "--cutoff",
            "3.5",
            "--per-atom",
            str(per_atom),
        ],
        check=False,
        capture_output=True,
        text=True,
    )
    if proc.returncode != 0:
        raise SystemExit(proc.stdout + proc.stderr)
    # per-atom dump column "cage": 0 water, 1 hexagonal, 2 cubic, 3 both
    flags = {}
    lines = per_atom.read_text().splitlines()
    i = 0
    while i < len(lines):
        if lines[i].startswith("ITEM: ATOMS"):
            n_line = lines[i - 5] if False else None
            break
        i += 1
    # Walk atoms after the header.
    header = next(k for k, line in enumerate(lines) if line.startswith("ITEM: ATOMS"))
    n = int(lines[header - 3]) if lines[header - 4] == "ITEM: NUMBER OF ATOMS" else None
    # BOX is three lines; NUMBER OF ATOMS value is two lines before BOX item.
    # Simpler: every atom line is numeric.
    for line in lines[header + 1 :]:
        if line.startswith("ITEM:"):
            break
        parts = line.split()
        if len(parts) < 6:
            continue
        flags[int(parts[0])] = int(parts[5])
    if n is not None and len(flags) != n:
        raise SystemExit(f"cage dump has {len(flags)} atoms, expected {n}")
    text = proc.stdout
    return flags, text


def performance(log: Path) -> str:
    last = ""
    if not log.is_file():
        return "no log"
    for line in log.read_text(errors="replace").splitlines():
        if line.startswith("Performance:") or line.startswith("Loop time of"):
            last = line
    return last


def final_step(log: Path) -> int | None:
    """Absolute step from the last thermo row. A restart keeps that step."""
    found = None
    if not log.is_file():
        return None
    for line in log.read_text(errors="replace").splitlines():
        parts = line.split()
        if len(parts) < 4 or not parts[0].isdigit():
            continue
        try:
            step = int(parts[0])
            float(parts[1])
        except ValueError:
            continue
        found = step
    return found


def load_series(work: Path, tag: str):
    """Original dump plus every continuation, one frame per step."""
    dest = work / tag
    paths = []
    original = dest / "traj.lammpstrj"
    if original.is_file():
        paths.append(original)
    paths.extend(sorted(dest.glob("cont*.lammpstrj")))
    if not paths:
        raise SystemExit(f"missing dumps for {tag}")
    merged = {}
    for path in paths:
        for frame in parse_dump(path):
            merged[frame[0]] = frame
    return [merged[step] for step in sorted(merged)]


def dump_restart(lmp: str, restart: Path, dump: Path) -> None:
    """Dump coordinates from a restart with the resume deck and run 0.

    The restart already holds the positions, the box and the step. Run 0
    does not integrate. The scratch directory keeps the deck's restart.a
    and end.restart away from the dataset.
    """
    # One directory per labelled frame, so concurrent loadings do not share it.
    scratch = dump.parent
    scratch.mkdir(parents=True, exist_ok=True)
    dump.unlink(missing_ok=True)
    cmd = [
        lmp,
        "-log",
        "log.lammps",
        "-in",
        str(DECK / "in.resume.lmp"),
        "-var",
        "restart",
        str(restart.resolve()),
        "-var",
        "T",
        f"{TEMPERATURE:.2f}",
        "-var",
        "P",
        f"{PRESSURE:.1f}",
        "-var",
        "steps",
        "0",
        "-var",
        "dumpevery",
        "1",
        "-var",
        "dump",
        str(dump.resolve()),
    ]
    proc = subprocess.run(cmd, cwd=scratch, capture_output=True, text=True)
    if proc.returncode != 0:
        raise SystemExit(proc.stdout + "\n" + proc.stderr)


def write_labelled_xyz(dump: Path, out_xyz: Path, publish: Path | None, seams: str, lib) -> int:
    """CHILL+ the dumped frame and write the extended XYZ."""
    parsed = parse_dump(dump)
    if not parsed:
        raise SystemExit(f"no frame in {dump}")
    step, origin, length, atoms = parsed[-1]
    labels = chill_labels(lib, atoms, origin, length)
    library_counts = {}
    for i, atom in enumerate(atoms):
        if atom[1] != 1:
            continue
        library_counts[CHILL[labels[i]]] = library_counts.get(CHILL[labels[i]], 0) + 1
    cli_counts = cli_chill_counts(seams, dump)
    if library_counts != cli_counts:
        raise SystemExit(f"{out_xyz.name}: library {library_counts} cli {cli_counts}")
    out_xyz.parent.mkdir(parents=True, exist_ok=True)
    write_xyz(out_xyz, step, origin, length, atoms, labels)
    if publish is not None:
        publish.parent.mkdir(parents=True, exist_ok=True)
        shutil.copyfile(out_xyz, publish)
    print(f"wrote {out_xyz.name} atoms={len(atoms)} step={step}", flush=True)
    return step


def frame_from_restart(
    restart: Path,
    out_xyz: Path,
    publish: Path | None,
    lmp: str,
    seams: str,
    lib_path: str,
    expect_step: int | None = None,
) -> None:
    lib = load_chill(lib_path)
    dump = out_xyz.parent / ".work" / out_xyz.stem / "frame.lammpstrj"
    dump_restart(lmp, restart, dump)
    step = write_labelled_xyz(dump, out_xyz, publish, seams, lib)
    if expect_step is not None and step != expect_step:
        raise SystemExit(f"{restart.name} dumped step {step}, expected {expect_step}")


def frames(restarts: Path, out: Path, frames_dir: Path | None, lmp: str, seams: str, lib_path: str) -> None:
    """Write the twelve labelled frames from step restarts only."""
    out.mkdir(parents=True, exist_ok=True)
    for pairs in PAIRS:
        tag = f"T230_m{pairs}_r0"
        for name, step in FRAME_STEPS:
            restart = restarts / tag / f"step-{step}.restart"
            if not restart.is_file():
                raise SystemExit(f"missing restart {restart}")
            publish = None if frames_dir is None else frames_dir / f"{tag}_{name}.xyz"
            frame_from_restart(
                restart,
                out / f"{tag}_{name}.xyz",
                publish,
                lmp,
                seams,
                lib_path,
                expect_step=step,
            )


CHUNK_STEPS = 125000  # 250 ps at 2 fs
UNTIL_STEP = 2500000  # 5 ns


def newest_restart(dest: Path) -> Path:
    """The completed chunk is end.restart. Periodic files are the in-chunk copies."""
    end = dest / "end.restart"
    if end.is_file():
        return end
    periodic = [path for path in (dest / "restart.a", dest / "restart.b") if path.is_file()]
    if not periodic:
        raise SystemExit(f"no restart in {dest}")
    return max(periodic, key=lambda path: path.stat().st_mtime)


def read_logged_step(dest: Path) -> int | None:
    text = dest / "step.txt"
    if text.is_file():
        return int(text.read_text().strip())
    return final_step(dest / "log.lammps")


def loop_steps(log: Path) -> int | None:
    found = None
    if not log.is_file():
        return None
    for line in log.read_text(errors="replace").splitlines():
        if not line.startswith("Loop time of"):
            continue
        parts = line.split()
        if "for" in parts:
            found = int(parts[parts.index("for") + 1])
    return found


def note_progress(dest: Path, started: int) -> int:
    # Thermo rows carry the absolute step, including after read_restart.
    thermo = final_step(dest / "log.lammps")
    if thermo is not None and thermo >= started:
        step = thermo
    else:
        done = loop_steps(dest / "log.lammps")
        if done is None:
            raise SystemExit(f"no step after the run in {dest}")
        step = started + done
    (dest / "step.txt").write_text(f"{step}\n")
    return step


def resume_chunk(dest: Path, lmp: str, threads: int, cores: str, steps: int) -> int:
    """Advance one loading by `steps`. Returns the LAMMPS exit code."""
    source = newest_restart(dest)
    current = dest / "current.restart"
    shutil.copyfile(source, current)
    env = os.environ.copy()
    env["OMP_NUM_THREADS"] = str(threads)
    env["OMP_PROC_BIND"] = "close"
    cmd = [
        "taskset",
        "-c",
        cores,
        lmp,
        "-sf",
        "omp",
        "-pk",
        f"omp {threads}",
        "-log",
        "log.lammps",
        "-in",
        str(DECK / "in.resume.lmp"),
        "-var",
        "restart",
        "current.restart",
        "-var",
        "T",
        f"{TEMPERATURE:.2f}",
        "-var",
        "P",
        f"{PRESSURE:.1f}",
        "-var",
        "steps",
        str(steps),
        "-var",
        "dumpevery",
        str(DUMP_EVERY),
        "-var",
        "dump",
        "cont.lammpstrj",
    ]
    log = open(dest / "resume.log", "a")
    print("resume", dest.name, "cores", cores, "steps", steps, flush=True)
    proc = subprocess.run(cmd, cwd=dest, env=env, stdout=log, stderr=subprocess.STDOUT)
    log.close()
    return proc.returncode


def core_plan(layout: str, threads: int) -> list[tuple[int, str]]:
    """(pairs, cpu set) for one wave of work."""
    if layout == "serial":
        last = threads - 1
        return [(pairs, f"0-{last}") for pairs in PAIRS]
    if layout == "concurrent":
        sets = ("0,1", "2,3", "4,5")
        if threads != 2:
            raise SystemExit("concurrent layout uses 2 threads on each of three loadings")
        return list(zip(PAIRS, sets))
    raise SystemExit(f"unknown layout {layout}")


def resume(work: Path, lmp: str, threads: int, layout: str, store: Path | None) -> None:
    """Run 250 ps chunks until every loading has reached 5 ns."""
    pending = True
    while pending:
        pending = False
        wave = []
        for pairs, cores in core_plan(layout, threads):
            dest = work / f"T230_m{pairs}_r0"
            step = read_logged_step(dest)
            if step is None:
                raise SystemExit(f"no thermo step in {dest}")
            if step >= UNTIL_STEP:
                print(f"{dest.name} already at step {step}", flush=True)
                continue
            pending = True
            take = min(CHUNK_STEPS, UNTIL_STEP - step)
            if layout == "serial":
                code = resume_chunk(dest, lmp, threads, cores, take)
                if code != 0:
                    raise SystemExit(f"{dest.name} failed ({code})")
                note_progress(dest, step)
                publish(dest, store)
            else:
                wave.append((dest, cores, take, step))
        if layout == "concurrent" and wave:
            procs = []
            for dest, cores, take, step in wave:
                procs.append((start_chunk(dest, lmp, threads, cores, take), step))
            codes = [proc.wait() for (proc, _dest), _step in procs]
            print("codes", codes, flush=True)
            if any(codes):
                raise SystemExit(f"LAMMPS failed: {codes}")
            for (proc, dest), step in procs:
                note_progress(dest, step)
                publish(dest, store)
            pending = any(
                (read_logged_step(work / f"T230_m{pairs}_r0") or 0) < UNTIL_STEP for pairs in PAIRS
            )


def start_chunk(dest: Path, lmp: str, threads: int, cores: str, steps: int):
    source = newest_restart(dest)
    current = dest / "current.restart"
    shutil.copyfile(source, current)
    env = os.environ.copy()
    env["OMP_NUM_THREADS"] = str(threads)
    env["OMP_PROC_BIND"] = "close"
    cmd = [
        "taskset",
        "-c",
        cores,
        lmp,
        "-sf",
        "omp",
        "-pk",
        f"omp {threads}",
        "-log",
        "log.lammps",
        "-in",
        str(DECK / "in.resume.lmp"),
        "-var",
        "restart",
        "current.restart",
        "-var",
        "T",
        f"{TEMPERATURE:.2f}",
        "-var",
        "P",
        f"{PRESSURE:.1f}",
        "-var",
        "steps",
        str(steps),
        "-var",
        "dumpevery",
        str(DUMP_EVERY),
        "-var",
        "dump",
        "cont.lammpstrj",
    ]
    log = open(dest / "resume.log", "a")
    print("resume", dest.name, "cores", cores, "steps", steps, flush=True)
    proc = subprocess.Popen(cmd, cwd=dest, env=env, stdout=log, stderr=subprocess.STDOUT)
    return proc, dest


def publish(dest: Path, store: Path | None) -> None:
    if store is None:
        return
    target = store / dest.name
    target.mkdir(parents=True, exist_ok=True)
    step = read_logged_step(dest)
    for name in ("end.restart", "restart.a", "restart.b"):
        src = dest / name
        if src.is_file():
            shutil.copyfile(src, target / name)
    if step is not None and (dest / "end.restart").is_file():
        shutil.copyfile(dest / "end.restart", target / f"step-{step}.restart")
    print(f"checkpoint {dest.name} step {step}", flush=True)


def series(work: Path, seams: str, lib_path: str, out: Path) -> None:
    """One row per saved frame: CHILL+ counts, the cage cluster, and the ions."""
    lib = load_chill(lib_path)
    header = (
        "t_ns,pairs,step,hexagonal,cubic,interfacial,water,clathrate,interClathrate,"
        "chill_ice,ion_liquid,ion_front,ion_ice,nmax"
    )
    rows = [header]
    for pairs in PAIRS:
        tag = f"T230_m{pairs}_r0"
        for step, origin, length, atoms in load_series(work, tag):
            labels = chill_labels(lib, atoms, origin, length)
            counts = {name: 0 for name in CHILL}
            for i, atom in enumerate(atoms):
                if atom[1] == 1:
                    counts[CHILL[labels[i]]] += 1
            one = work / tag / "series-frame.lammpstrj"
            write_one_frame_dump(one, step, origin, length, atoms)
            flags, _text = cage_flags(seams, one, work / tag / "series-cage.lammpstrj")
            oxy, ice, ions = [], [], []
            for aid, typ, x, y, z in atoms:
                if typ == 1:
                    oxy.append((x, y, z))
                    ice.append(flags.get(aid, 0) != 0)
                elif typ in (3, 4):
                    ions.append((x, y, z))
            nmax = largest_cluster(oxy, length, ice, k=4)
            n_ice, n_front, n_liq = ions_in_liquid(oxy, ice, ions, length)
            chill_ice = counts["hexagonal"] + counts["cubic"] + counts["reHex"] + counts["reCubic"]
            rows.append(
                f"{step * 2e-6:.4f},{pairs},{step},{counts['hexagonal']},{counts['cubic']},"
                f"{counts['interfacial']},{counts['water']},{counts['clathrate']},"
                f"{counts['interClathrate']},{chill_ice},{n_liq},{n_front},{n_ice},{nmax}"
            )
            print(rows[-1], flush=True)
    out.write_text("\n".join(rows) + "\n")
    print(f"wrote {out}", flush=True)


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("stage", choices=("prepare", "run", "resume", "frame", "frames", "series", "all"))
    parser.add_argument("--work", type=Path)
    parser.add_argument("--restarts", type=Path)
    parser.add_argument("--restart", type=Path)
    parser.add_argument("--out", type=Path)
    parser.add_argument("--publish", type=Path, help="copy one labelled frame to this path")
    parser.add_argument("--frames-dir", type=Path, help="copy the twelve labelled frames into this directory")
    parser.add_argument("--step", type=int)
    parser.add_argument("--genice", default="")
    parser.add_argument("--lmp", default="")
    parser.add_argument("--mpirun", default="")
    parser.add_argument("--seams", default="")
    parser.add_argument("--lib", default="")
    parser.add_argument("--engine", choices=("omp", "mpi"), default="omp")
    parser.add_argument("--threads", type=int, default=2)
    parser.add_argument("--layout", choices=("concurrent", "serial"), default="concurrent")
    parser.add_argument("--store", type=Path, default=None)
    args = parser.parse_args()
    needs_work = args.stage in ("prepare", "run", "resume", "series", "all")
    if needs_work and args.work is None:
        raise SystemExit("pass --work")
    if args.stage in ("prepare", "all"):
        genice = args.genice or which("genice2")
        prepare(args.work, genice)
    if args.stage in ("run", "all"):
        lmp = args.lmp or which("lmp")
        mpirun = args.mpirun or shutil.which("mpirun") or ""
        if args.engine == "mpi" and not mpirun:
            raise SystemExit("mpirun not found")
        procs = launch(args.work, args.engine, args.threads, lmp, mpirun)
        codes = [p.wait() for p in procs]
        print("codes", codes, flush=True)
        if any(codes):
            raise SystemExit(f"LAMMPS failed: {codes}")
    if args.stage == "resume":
        lmp = args.lmp or which("lmp")
        resume(args.work, lmp, args.threads, args.layout, args.store)
    if args.stage == "frame":
        if args.restart is None or args.out is None:
            raise SystemExit("frame needs --restart and --out")
        lmp = args.lmp or which("lmp")
        seams = args.seams or which("seams")
        if not args.lib:
            raise SystemExit("pass --lib path to libyodaLib.so")
        frame_from_restart(args.restart, args.out, args.publish, lmp, seams, args.lib, args.step)
    if args.stage in ("frames", "all"):
        if args.restarts is None:
            raise SystemExit("frames reads restarts; pass --restarts")
        lmp = args.lmp or which("lmp")
        seams = args.seams or which("seams")
        if not args.lib:
            raise SystemExit("pass --lib path to libyodaLib.so")
        frames(args.restarts, args.out or OUT, args.frames_dir, lmp, seams, args.lib)
    if args.stage == "series":
        seams = args.seams or which("seams")
        if not args.lib:
            raise SystemExit("pass --lib path to libyodaLib.so")
        series(args.work, seams, args.lib, HERE / "front.csv")


if __name__ == "__main__":
    main()
