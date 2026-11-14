# Copyright (c) d-SEAMS developers
"""compute dseams `type I` and `mol water` on water with ions, against the CLI.

usage: e2e_select.py BUILD [--ranks N]

A TIP4P/Ice ice Ic slab (O type 1, H type 2) with its top half jittered,
plus sodium ions (type 3, one molecule each) in the empty half of the box.
`type 1`, `mol water` and a group of the oxygens must give each oxygen the
label of `seams chill-plus --type 1` and `--mol water` on the same frame;
hydrogens and ions read 6. A group of several types with neither keyword,
a type that is not in the group, and `mol water` without molecule IDs are
errors. LMP names the lmp executable (default: lmp on PATH).
"""

import argparse
import os
import random
import shutil
import subprocess
import sys
import tempfile
from pathlib import Path

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE.parent))
import e2e_common as e2e  # noqa: E402
import tip4pice_slab  # noqa: E402

UNCLASSIFIED = 6
IONS = 6

DECK = """units real
atom_style {style}
boundary p p p
read_data system.data
mass * 1.0
pair_style zero 4.0
pair_coeff * *
plugin load {plugin}
group oxygen type 1
{computes}
dump d all custom 1 frame.lammpstrj id mol type x y z {columns}
dump_modify d sort id format float %.10g
run 0
"""


def system(seed=7):
    mols, box = tip4pice_slab.slab(4, 4, 2)
    rng = random.Random(seed)
    top = box[2] / 4
    atoms = []
    for m, (o, h1, h2, _) in enumerate(mols, 1):
        shift = (0.0, 0.0, 0.0)
        if o[2] > top:
            shift = tuple(rng.uniform(-0.6, 0.6) for _ in range(3))
        for t, x in ((1, o), (2, h1), (2, h2)):
            atoms.append((m, t, tuple(c + s for c, s in zip(x, shift))))
    for k in range(IONS):
        pos = (rng.uniform(0, box[0]), rng.uniform(0, box[1]), rng.uniform(0.6, 0.9) * box[2])
        atoms.append((len(mols) + 1 + k, 3, pos))
    return atoms, box


def write_data(path, atoms, box, style):
    with open(path, "w") as fh:
        fh.write("water oxygens, hydrogens and ions\n\n%d atoms\n4 atom types\n\n" % len(atoms))
        for lo, hi in zip(("x", "y", "z"), box):
            fh.write("0 %.6f %slo %shi\n" % (hi, lo, lo))
        fh.write("\nAtoms # %s\n\n" % style)
        for i, (mol, t, x) in enumerate(atoms, 1):
            if style == "full":
                fh.write("%d %d %d 0.0 %.6f %.6f %.6f\n" % (i, mol, t, *x))
            else:
                fh.write("%d %d %.6f %.6f %.6f\n" % (i, t, *x))


def run(launch, work, plugin, computes, columns, style="full"):
    write_data(work / "system.data", *system(), style)
    (work / "in.select.lmp").write_text(DECK.format(
        style=style, plugin=plugin, computes="\n".join(computes), columns=columns))
    for stale in ("frame.lammpstrj", "log.lammps"):
        (work / stale).unlink(missing_ok=True)
    done = subprocess.run(launch + ["-in", "in.select.lmp", "-log", "log.lammps"], cwd=work,
                          env=dict(os.environ, OMP_NUM_THREADS="1"),
                          capture_output=True, text=True)
    log = (work / "log.lammps").read_text() if (work / "log.lammps").is_file() else ""
    return done.returncode, done.stdout + done.stderr + log


def expect_error(launch, work, plugin, compute, message, style="full"):
    rc, out = run(launch, work, plugin, [compute], "c_bad", style)
    if rc == 0 or message not in out:
        sys.stderr.write(out[-3000:])
        sys.exit(f"`{compute}` ({style}) did not fail with: {message}")


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("build")
    ap.add_argument("--ranks", type=int, default=1)
    args = ap.parse_args()

    lmp = e2e.which("lmp", "LMP")
    if lmp is None:
        e2e.skip("LAMMPS (lmp or $LMP) was not found")
    plugin = e2e.find_build_file(args.build, "dseamsplugin.so")
    if plugin is None:
        e2e.skip("dseamsplugin.so was not built (LAMMPS class headers missing)")
    cli = e2e.seams_cli(args.build)
    launch = [lmp]
    if args.ranks > 1:
        mpirun = shutil.which("mpirun")
        if mpirun is None:
            e2e.skip("mpirun was not found")
        launch = [mpirun, "-np", str(args.ranks), lmp]

    with tempfile.TemporaryDirectory(prefix="dseams-lammps-select-") as tmp:
        work = Path(tmp)
        rc, out = run(launch, work, plugin, [
            "compute sel all dseams type 1",
            "compute mol all dseams mol water",
            "compute ox oxygen dseams",
        ], "c_sel c_mol c_ox")
        if rc != 0:
            if "Unrecognized atom style" in out or "Invalid atom style" in out:
                e2e.skip("this LAMMPS has no atom_style full (MOLECULE package)")
            sys.stderr.write(out[-4000:])
            sys.exit("lmp failed")
        dump = work / "frame.lammpstrj"
        types = e2e.read_dump(dump, "type")[0][1]
        oxygens = {i for i, t in types.items() if t == 1}
        columns = {c: e2e.read_dump(dump, c)[0][1] for c in ("c_sel", "c_mol", "c_ox")}
        bad = []
        for flags, column in (
            (("--type", "1"), "c_sel"),
            (("--mol", "water"), "c_mol"),
            (("--type", "1"), "c_ox"),
        ):
            offline = e2e.offline_labels(cli, dump, 1, 3.5, work, extra=flags)[0][1]
            if set(offline) != oxygens:
                sys.exit(f"seams chill-plus {' '.join(flags)} labelled atoms other than the oxygens")
            mine = columns[column]
            for i in sorted(types):
                want = offline[i] if i in oxygens else (-1 if column == "c_ox" else UNCLASSIFIED)
                if mine[i] != want:
                    bad.append((column, i, types[i], mine[i], want))
        labels = sorted({columns["c_sel"][i] for i in oxygens})
        print(f"{len(oxygens)} oxygens, {len(types) - len(oxygens)} other atoms, "
              f"oxygen labels {labels}")
        if bad:
            for column, i, t, mine, want in bad[:20]:
                print(f"{column} id {i} type {t}: {mine}, want {want}", file=sys.stderr)
            sys.exit(f"{len(bad)} labels differ")
        if len(labels) < 2:
            sys.exit("the frame has a single oxygen label; the check is too weak")

        expect_error(launch, work, plugin, "compute bad all dseams",
                     "mixed species need a water type or molecule")
        expect_error(launch, work, plugin, "compute bad all dseams type 4",
                     "atom type 4 is not in the group")
        expect_error(launch, work, plugin, "compute bad all dseams mol water",
                     "mol water needs an atom style with molecule IDs", style="atomic")
    print("compute dseams type and mol water match the CLI")


if __name__ == "__main__":
    main()
