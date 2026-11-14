# Copyright (c) d-SEAMS developers
"""TIP4P/Ice ice Ic slab and liquid in GROMACS, classified by d-SEAMS.

usage: e2e_tip4pice.py BUILD [--plumed] [--steps S] [--record FILE]

Needs gmx (or $GMX). Without --plumed, gmx-seams runs its analyzeFrame
hook on every frame of mdrun's full-precision .trr; GROMACS will not read
a trajectory from a pipe, so the hook runs after mdrun. With --plumed,
mdrun -plumed loads libSeamsPlumed.so and SEAMS classifies the oxygens
inside the MD loop; gmx-seams then also classifies the same steps and the
two in-situ paths must agree. Either way the labels are checked against
the CLI on the dumped frames, and the dumped oxygens against the .trr.
"""

import argparse
import os
import re
import shutil
import subprocess
import sys
import tempfile
from pathlib import Path

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE.parent))
import e2e_common as e2e  # noqa: E402
import tip4pice_slab  # noqa: E402

PLUMED_DECK = """LOAD FILE={plugin}
seams: SEAMS ATOMS=1-{natoms}:4 CUTOFF=3.5 DUMP=insitu.lammpstrj
PRINT ARG=seams.q6,seams.ice STRIDE={every} FILE=colvar
"""


def run(cmd, cwd, env, what):
    done = subprocess.run(cmd, cwd=cwd, env=env, capture_output=True, text=True)
    if done.returncode != 0:
        sys.stderr.write(done.stdout[-3000:] + done.stderr[-3000:])
        sys.exit(f"{what} failed")
    return done.stdout + done.stderr


def data_dir(gmx, env):
    text = run([gmx, "-quiet", "--version"], None, env, "gmx --version")
    found = re.search(r"Data prefix:\s*(\S+)", text)
    if found is None:
        sys.exit("gmx --version did not report its data prefix")
    return Path(found.group(1)) / "share" / "gromacs" / "top"


def stack(ice, liquid, dz, out):
    """Write ice with the liquid shifted up by dz nm, in the ice box."""
    ice_lines = Path(ice).read_text().splitlines()
    liq_lines = Path(liquid).read_text().splitlines()
    atoms = ice_lines[2:-1]
    for line in liq_lines[2:-1]:
        z = float(line[36:44]) + dz
        atoms.append(line[:36] + "%8.3f" % z + line[44:])
    with open(out, "w") as fh:
        fh.write("TIP4P/Ice ice Ic slab and liquid\n%d\n" % len(atoms))
        for k, line in enumerate(atoms):
            res, atom = k // 4 + 1, k + 1
            fh.write("%5d%s%5d%s\n" % (res % 100000, line[5:15], atom % 100000, line[20:44]))
        fh.write(ice_lines[-1] + "\n")
    return len(atoms)


def trr_oxygens(gmx, trr, env):
    """Oxygen x of every .trr frame, in Angstrom, keyed by 1-based atom id."""
    text = run([gmx, "-quiet", "dump", "-f", str(trr)], None, env, "gmx dump")
    frames = []
    current = None
    for line in text.splitlines():
        line = line.strip()
        if line.startswith("x (") and line.endswith("):"):
            current = {}
            frames.append(current)
        elif current is not None and line.startswith("x["):
            idx = int(line[2:line.index("]")])
            if idx % 4 == 0:
                vals = line[line.index("{") + 1:line.index("}")].split(",")
                current[idx + 1] = tuple(float(v) * 10.0 for v in vals)
        elif line.startswith(("v (", "f (")):
            current = None
    return frames


def dump_positions(path):
    """x, y, z of every dump frame keyed by atom id."""
    frames = []
    lines = iter(Path(path).read_text().splitlines())
    n = 0
    for line in lines:
        if line.startswith("ITEM: NUMBER OF ATOMS"):
            n = int(next(lines))
        elif line.startswith("ITEM: ATOMS"):
            cols = line.split()[2:]
            at = [cols.index(c) for c in ("id", "x", "y", "z")]
            frame = {}
            for _ in range(n):
                parts = next(lines).split()
                frame[int(parts[at[0]])] = tuple(float(parts[c]) for c in at[1:])
            frames.append(frame)
    return frames


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("build")
    ap.add_argument("--plumed", action="store_true")
    ap.add_argument("--cells", type=int, nargs=3, default=(6, 6, 3))
    ap.add_argument("--steps", type=int, default=200)
    ap.add_argument("--every", type=int, default=10)
    ap.add_argument("--threads", type=int, default=2)
    ap.add_argument("--record", default="")
    args = ap.parse_args()

    gmx = e2e.which("gmx", "GMX")
    if gmx is None:
        e2e.skip("GROMACS (gmx or $GMX) was not found")
    hook = e2e.find_build_file(args.build, "gmx-seams")
    if hook is None:
        e2e.skip("gmx-seams was not built (set gromacs_legacy_include for GROMACS 2023+)")
    env = dict(os.environ, GMX_MAXBACKUP="-1", GMX_SUPPRESS_DUMP="1",
               OMP_NUM_THREADS=str(args.threads), PLUMED_NUM_THREADS="1")
    plugin = None
    if args.plumed:
        kernel = e2e.plumed_kernel()
        if kernel is None:
            e2e.skip("no libplumedKernel.so in PLUMED_KERNEL or beside plumed on PATH")
        env["PLUMED_KERNEL"] = kernel
        mdrun_help = subprocess.run([gmx, "-quiet", "mdrun", "-h"], env=env,
                                    capture_output=True, text=True).stdout
        if "-plumed" not in mdrun_help:
            e2e.skip("this gmx mdrun is not patched with PLUMED")
        plugin = e2e.find_build_file(args.build, "libSeamsPlumed.so")
        if plugin is None:
            e2e.skip("libSeamsPlumed.so was not built")
        # GROMACS 2024 sorts atoms under single-rank domain decomposition and
        # the PLUMED 2.9 patch then hands PLUMED about 5% of the oxygens with
        # a neighbour's coordinates. The .trr comparison below catches it.
        env["GMX_DD_SINGLE_RANK"] = "0"
    cli = e2e.seams_cli(args.build)
    solvent = data_dir(gmx, env) / "tip4p.gro"

    with tempfile.TemporaryDirectory(prefix="dseams-gromacs-") as tmp:
        work = Path(tmp)
        mols, box = tip4pice_slab.slab(*args.cells)
        tip4pice_slab.write_gro(work / "ice.gro", mols, box)
        for name in ("tip4pice.top", "tip4pice.mdp", "tip4pice-em.mdp"):
            shutil.copy(HERE / name, work)
        # Solvating the slab itself puts waters into the lattice
        # interstitials, so the liquid is a separate box stacked on top.
        half = box[2] / 20
        run([gmx, "-quiet", "solvate", "-cs", str(solvent), "-box", str(box[0] / 10),
             str(box[1] / 10), str(half - 0.2), "-o", "liquid.gro"], work, env, "gmx solvate")
        natoms = stack(work / "ice.gro", work / "liquid.gro", half + 0.1, work / "system.gro")
        top = (work / "tip4pice.top").read_text() + f"SOL {natoms // 4}\n"
        (work / "topol.top").write_text(top)
        mdp = (work / "tip4pice.mdp").read_text()
        mdp = re.sub(r"(?m)^nsteps.*$", f"nsteps = {args.steps}", mdp)
        mdp = re.sub(r"(?m)^nstxout .*$", f"nstxout = {args.every}", mdp)
        (work / "md.mdp").write_text(mdp)
        run([gmx, "-quiet", "grompp", "-f", "tip4pice-em.mdp", "-c", "system.gro", "-p", "topol.top",
             "-o", "em.tpr"], work, env, "grompp (minimisation)")
        run([gmx, "-quiet", "mdrun", "-s", "em.tpr", "-deffnm", "em", "-ntmpi", "1"], work, env,
            "mdrun (minimisation)")
        run([gmx, "-quiet", "grompp", "-f", "md.mdp", "-c", "em.gro", "-p", "topol.top",
             "-o", "md.tpr"], work, env, "grompp")
        mdrun = [gmx, "-quiet", "mdrun", "-s", "md.tpr", "-deffnm", "md", "-ntmpi", "1"]
        if plugin is not None:
            (work / "plumed.dat").write_text(
                PLUMED_DECK.format(plugin=plugin, natoms=natoms, every=args.every))
            mdrun += ["-plumed", "plumed.dat"]
        log = run(mdrun, work, env, "mdrun")
        version = re.search(r"GROMACS version:\s*(\S+)", (work / "md.log").read_text())

        hook_dump = "hook.lammpstrj" if plugin is not None else "insitu.lammpstrj"
        out = run([str(hook), "-quiet", "-s", "md.tpr", "-f", "md.trr", "-select", "name OW",
                   "-dump", hook_dump], work, env, "gmx-seams")
        if plugin is not None:
            # The PLUMED patch logs into md.log; GROMACS 2025 and later, which
            # link PLUMED themselves, write PLUMED.OUT.
            plog = (work / "md.log").read_text() + log
            if (work / "PLUMED.OUT").is_file():
                plog += (work / "PLUMED.OUT").read_text()
            timing = re.search(r"SEAMS calls (\d+) seconds (\S+) per_call (\S+)", plog)
            if timing is None:
                sys.exit("the SEAMS action did not report its timing")
            per_frame = float(timing.group(3))
        else:
            timing = re.search(r"overhead_s (\S+)", out)
            if timing is None:
                sys.exit("gmx-seams did not report its timing")
            per_frame = float(timing.group(1))

        insitu = e2e.read_dump(work / "insitu.lammpstrj", "label")
        offline = e2e.offline_labels(cli, work / "insitu.lammpstrj", len(insitu), 3.5, work,
                                     extra=("--type", "1"))
        frames, atom_frames, bad = e2e.compare(insitu, offline)

        trr = trr_oxygens(gmx, work / "md.trr", env)
        dumped = dump_positions(work / "insitu.lammpstrj")
        if len(trr) != len(dumped):
            sys.exit(f"{len(dumped)} dumped frames but {len(trr)} .trr frames")
        # mdrun hands PLUMED wrapped coordinates; the box is fixed and orthorhombic.
        drift = max(abs((a - b) - side * round((a - b) / side))
                    for t, d in zip(trr, dumped) for k in d
                    for a, b, side in zip(t[k], d[k], box))

        record = {
            "engine": "gromacs-plumed" if plugin is not None else "gromacs",
            "version": version.group(1) if version else "unknown",
            "system": f"TIP4P/Ice ice Ic slab and liquid, {natoms // 4} molecules",
            "steps": args.steps,
            "frames": frames,
            "atom_frames": atom_frames,
            "labels": e2e.census(insitu),
            "insitu_s_per_frame": per_frame,
            "dump_vs_trr_max_abs_A": drift,
        }
        if plugin is not None:
            hooked = e2e.read_dump(work / hook_dump, "label")
            if [s for s, _ in hooked] != [s for s, _ in insitu]:
                sys.exit("mdrun -plumed and gmx-seams saw different steps")
            _, _, cross = e2e.compare(insitu, hooked)
            record["native_dump_mismatches"] = len(cross)
        e2e.finish(args.record, record, bad)
        # gmx dump prints six significant figures of nm.
        if drift > 5e-4:
            sys.exit(f"dumped oxygens differ from the .trr by {drift} A")


if __name__ == "__main__":
    main()
