"""Exercise each seams command and the errors a user can hit.

usage: cli_surface.py SEAMS
"""

import os
import subprocess
import sys
import tempfile


def run(seams, args, env=None):
    return subprocess.run([seams, *args], capture_output=True, text=True, env=env)


def _write_interface(path):
    source = "input/traj/interface_f0.lammpstrj"
    with open(source, encoding="ascii") as src, open(path, "w", encoding="ascii") as dst:
        dst.write(src.read())


def _write_gro(path):
    with open(path, "w", encoding="ascii") as handle:
        handle.write(
            "frame one\n    1\n"
            "    1SOL     OW    1   0.000   0.000   0.000\n"
            "   2.00000   2.00000   2.00000\n"
            "frame two\n    1\n"
            "    1SOL     OW    1   0.100   0.100   0.000\n"
            "   2.00000   2.00000   2.00000\n"
        )


def _write_pdb(path):
    with open(path, "w", encoding="ascii") as handle:
        handle.write(
            "ATOM      1  OW  SOL     1       0.000   0.000   0.000"
            "  1.00  0.00           O\nEND\n"
        )


def _write_run(path):
    # Frame 1 is a molecule, frame 2 is empty, frame 3 is a molecule.
    # A strict read must keep the empty frame's failure after frame 3.
    chunks = []
    for step, n in ((0, 1), (1, 0), (2, 1)):
        body = "" if n == 0 else "1 1 0 0 0\n"
        chunks.append(
            f"ITEM: TIMESTEP\n{step}\nITEM: NUMBER OF ATOMS\n{n}\n"
            "ITEM: BOX BOUNDS pp pp pp\n0 10\n0 10\n0 10\n"
            "ITEM: ATOMS id type x y z\n" + body
        )
    with open(path, "w", encoding="ascii") as handle:
        handle.write("".join(chunks))


def _write_gap(path):
    with open(path, "w", encoding="ascii") as handle:
        handle.write(
            "ITEM: TIMESTEP\n0\nITEM: NUMBER OF ATOMS\n1\n"
            "ITEM: BOX BOUNDS pp pp pp\n0 10\n0 10\n0 10\n"
            "ITEM: ATOMS id type x y z\n1 1 0 0 0\n"
            "ITEM: TIMESTEP\n1\nITEM: NUMBER OF ATOMS\n0\n"
            "ITEM: BOX BOUNDS pp pp pp\n0 10\n0 10\n0 10\n"
            "ITEM: ATOMS id type x y z\n"
        )


def _write_three_types(path):
    with open(path, "w", encoding="ascii") as handle:
        handle.write(
            "ITEM: TIMESTEP\n0\nITEM: NUMBER OF ATOMS\n3\n"
            "ITEM: BOX BOUNDS pp pp pp\n0 20\n0 20\n0 20\n"
            "ITEM: ATOMS id type x y z\n"
            "1 1 0.0 0.0 0.0\n"
            "2 2 3.0 0.0 0.0\n"
            "3 3 0.0 3.0 0.0\n"
        )


def _write_ions(path):
    with open(path, "w", encoding="ascii") as handle:
        handle.write(
            "ITEM: TIMESTEP\n0\nITEM: NUMBER OF ATOMS\n2\n"
            "ITEM: BOX BOUNDS pp pp pp\n0 20\n0 20\n0 20\n"
            "ITEM: ATOMS id type x y z\n"
            "1 1 0.0 0.0 0.0\n"
            "2 2 2.0 0.0 0.0\n"
        )


def main():
    seams = sys.argv[1]
    ice = "input/traj/ice1h.lammpstrj"
    aliased = run(seams, ["rdf", ice, "--types", "1,1", "--cutoff", "3.5"])
    paired = run(seams, ["rdf", ice, "--pair-types", "1,1", "--cutoff", "3.5"])
    if (aliased.returncode != 0 or paired.returncode != 0
            or aliased.stdout != paired.stdout):
        sys.stderr.write(
            "types alias output\n" + aliased.stdout + paired.stdout
            + aliased.stderr + paired.stderr)
        return 1
    if "--types is deprecated; use --pair-types" not in aliased.stderr:
        sys.stderr.write("types alias warning\n" + aliased.stderr)
        return 1
    if "deprecated" in paired.stderr:
        sys.stderr.write("pair-types warning\n" + paired.stderr)
        return 1
    if os.environ.get("SEAMS_HAS_CHEMFILES", "0") != "1":
        for args in (
            ["read", "missing-cli.pdb"],
            ["read", "missing-cli.gro"],
            ["--frame", "1", "--last", "2", "read", "missing-cli.pdb"],
            ["--frame", "1", "--last", "2", "read", "missing-cli.dcd"],
        ):
            probed = run(seams, args)
            if (probed.returncode != 2
                    or "chemfiles is not in this build" not in probed.stderr):
                sys.stderr.write(
                    "missing chemfiles reader " + " ".join(args) + " "
                    + str(probed.returncode) + "\n" + probed.stderr)
                return 1
        # The rest of this surface reads pdb, gro, and dcd.
        return 77
    tip = "input/traj/exampleTraj.lammpstrj"
    cubic = "input/traj/mW_cubic.lammpstrj"
    guest = "input/traj/genice_sI.lammpstrj"
    work = tempfile.mkdtemp(prefix="dseams-cli-")
    three = os.path.join(work, "three.lammpstrj")
    _write_three_types(three)
    empty = os.path.join(work, "empty.lammpstrj")
    with open(empty, "w", encoding="ascii") as handle:
        handle.write(
            "ITEM: TIMESTEP\n0\nITEM: NUMBER OF ATOMS\n0\n"
            "ITEM: BOX BOUNDS pp pp pp\n0.0 10.0\n0.0 10.0\n0.0 10.0\n"
            "ITEM: ATOMS id type x y z\n"
        )
    per_atom = os.path.join(work, "per-atom.lammpstrj")
    classes = os.path.join(work, "classes.lammpstrj")
    cages_atom = os.path.join(work, "cages.lammpstrj")
    gap = os.path.join(work, "gap.lammpstrj")
    _write_gap(gap)
    run_path = os.path.join(work, "run.lammpstrj")
    _write_run(run_path)
    nodot = os.path.join(work, "nodot")
    with open(nodot, "w", encoding="ascii") as handle:
        handle.write(
            "ITEM: TIMESTEP\n0\nITEM: NUMBER OF ATOMS\n1\n"
            "ITEM: BOX BOUNDS pp pp pp\n0 10\n0 10\n0 10\n"
            "ITEM: ATOMS id type x y z\n1 1 0 0 0\n"
        )
    xyz = "input/xyz/configuration.xyz"
    ok = [
        ["read", tip],
        ["--frame", "1", "read", xyz],
        ["read", nodot],
        ["--frame", "2", "read", xyz],
        ["--frame", "1", "--last", "1", "read", xyz],
        ["--format", "json", "read", tip],
        ["--frame", "1", "--last", "2", "--jobs", "2", "read", cubic],
        ["chill", tip, "--cutoff", "3.5", "--type", "2"],
        ["chill-plus", tip, "--cutoff", "3.5", "--type", "2"],
        ["ovito-chill+", tip, "--cutoff", "3.5", "--type", "2"],
        ["chill_plus", ice, "--type", "1"],
        ["chill-plus", ice, "--rule", "examples/chill-plus.lua", "--type", "1"],
        ["chill-plus", tip, "--graph", "mutual", "--type", "2"],
        ["chill-plus", tip, "--layers", "--type", "2"],
        ["cages", tip, "--type", "2"],
        ["rings", tip, "--type", "2"],
        ["prisms", ice, "--type", "1"],
        ["cluster", tip, "--type", "2"],
        ["cages", tip, "--graph", "cutoff", "--cutoff", "3.5", "--type", "2"],
        ["cages", ice, "--graph", "mutual", "--type", "1", "--signature", "512"],
        ["fingerprint", tip, "--type", "2", "--hops", "1"],
        ["fingerprint", tip, "--type", "2", "--by-atom"],
        ["f4", tip, "--type", "2", "--htype", "1"],
        ["steinhardt", tip, "--type", "2", "--cutoff", "3.5"],
        ["rdf", tip, "--type", "2", "--cutoff", "6", "--bins", "8"],
        ["cn", tip, "--pair-types", "2,2", "--cutoff", "3.5"],
        ["hbonds", tip, "--type", "2", "--htype", "1"],
        ["density-z", tip, "--type", "2", "--axis", "x"],
        ["density-z", tip, "--axis", "y", "--bins", "4"],
        ["cages", ice, "--type", "1", "--signature", "512", "--incomplete"],
        ["fingerprint", tip, "--type", "2", "--by-atom"],
        ["ions", tip, "--ion-types", "1,2", "--site", "1=cationHead,2=anion", "--type", "2"],
        ["pairs", tip, "--site", "1=cationHead,2=anion"],
        ["domains", tip, "--site", "1=polar,2=apolar", "--subset", "apolar"],
        ["hbonds", tip, "--type", "2", "--htype", "1", "--donors"],
        ["--tpp", "2", "--block", "32", "--resident", "0.5", "--cell", "4", "read", tip],
        ["chill-plus", ice, "--type", "1", "--tum-layers"],
        ["steinhardt", ice, "--type", "1"],
        ["cages", ice, "--type", "1", "--signature", "hc", "--cage-edge", "none"],
        ["--print-config"],
        ["--features"],
        ["--version"],
        ["--help"],
        ["--format", "json", "chill-plus", ice, "--type", "1"],
        ["--frame", "1", "--last", "2", "--jobs", "1", "--rebuild", "read", cubic],
        ["--strict-input", "read", tip],
        ["--axis", "X", "density-z", tip, "--type", "2"],
        ["--axis", "0", "density-z", tip],
        ["--axis", "1", "density-z", tip],
        ["--axis", "z", "density-z", tip, "--type", "0", "--bins", "4"],
        ["--water-types", "2", "chill", tip, "--type", "2"],
        ["--water-types", "1,2", "chill-plus", tip],
        ["rdf", tip, "--pair-types", "2,1", "--cutoff", "6", "--bins", "4"],
        ["cn", tip, "--pair-types", "1,2", "--cutoff", "3.5", "--bins", "4"],
        ["cn", tip, "--ions", "--site", "1=cationHead,2=anion", "--cutoff", "3.5"],
        ["f4", empty, "--type", "1"],
        ["steinhardt", empty, "--type", "1"],
        ["cages", empty, "--type", "1"],
        ["cages", empty, "--signature", "512", "--type", "1"],
        ["fingerprint", empty, "--type", "1"],
        ["chill", empty, "--type", "1"],
        ["rdf", empty, "--type", "1"],
        ["cages", guest, "--type", "1", "--signature", "512+51268", "--incomplete"],
        ["cages", guest, "--type", "1", "--signature", "512", "--guest-types", "1",
         "--guest-radius", "4"],
        ["cages", guest, "--type", "1", "--signature", "512", "--guest-types", "1",
         "--inside"],
        ["cages", guest, "--type", "1", "--signature", "512", "--shape",
         "--cage-edge", "2.0,4.0"],
        ["cages", "input/traj/nucleation_f1.lammpstrj", "--type", "1", "--per-atom", cages_atom],
        ["cages", ice, "--type", "1", "--per-atom", cages_atom],
        ["cages", ice, "--type", "1", "--graph", "cutoff", "--water-types", "1"],
        ["fingerprint", ice, "--type", "1", "--hops", "1", "--colour-types"],
        ["fingerprint", ice, "--type", "1", "--hops", "1", "--per-atom", classes],
        ["fingerprint", ice, "--type", "1", "--graph", "mutual", "--hops", "1"],
        ["ions", tip, "--ion-types", "1", "--ion-cutoff", "3.0",
         "--site", "1=cationHead,2=anion", "--per-atom", per_atom],
        ["domains", tip, "--site", "1=polar,2=apolar", "--subset", "polar"],
        ["domains", tip, "--site", "1=polar,2=apolar", "--subset", "tail"],
        ["hbonds", tip, "--type", "2", "--htype", "1", "--hdist", "2.5",
         "--hangle", "30"],
        ["read", "input/xyz/configuration.xyz"],
        ["-k", "4", "cages", ice, "--type", "1", "--graph", "mutual"],
        ["--water-types", "1", "cages", ice, "--type", "1"],
        ["--water-types", "1", "cages", ice, "--type", "1", "--graph", "mutual",
         "--signature", "512"],
        ["--frame", "1", "--last", "2", "read", "input/xyz/configuration.xyz"],
        ["f4", tip, "--type", "0", "--htype", "1"],
        ["--frame", "1", "--last", "2", "--jobs", "1", "cages", cubic, "--type", "1"],
        ["--frame", "1", "--last", "1", "read", "input/xyz/configuration.xyz"],
        ["hbonds", empty, "--type", "1", "--htype", "1"],
        ["ions", empty, "--ion-types", "1", "--site", "1=cationHead,2=anion"],
        ["steinhardt", tip, "--type", "2"],
    ]
    failed = False
    published = run(seams, ["chill-plus", tip, "--cutoff", "3.5", "--type", "2"])
    if published.returncode != 0 or "interClathrate 9" not in published.stdout or "water 241" not in published.stdout:
        sys.stderr.write("published counts\n" + published.stdout + published.stderr)
        return 1
    mutual = run(seams, ["chill-plus", tip, "--graph", "mutual", "--cutoff", "3.5", "--type", "2"])
    if mutual.returncode != 0 or "interClathrate 7" not in mutual.stdout or "water 242" not in mutual.stdout:
        sys.stderr.write("mutual counts\n" + mutual.stdout + mutual.stderr)
        return 1
    checked = [
        (["ovito-chill+", ice, "--type", "1"], "hexagonal"),
        (["ovito-chill+", ice, "--type", "1", "--rule", "chill+"], "hexagonal"),
        (["ovito_chill_plus", ice, "--type", "1"], "hexagonal"),
        (["--frame", "1", "--last", "2", "--type", "1", "read", cubic], "nop"),
        (["--frame", "1", "--last", "2", "--type", "0", "read", cubic], "nop"),
        (["--frame", "1", "--last", "2", "--type", "0", "density-z", tip, "--bins", "2"], "rho"),
        (["cages", empty, "--type", "1"], "both"),
        (["cages", empty, "--type", "1"], "water"),
        (["domains", tip, "--site", "1=polar,2=apolar"], "subset"),
        (["--water-types", "1", "cages", ice, "--type", "1", "--signature", "512", "--graph", "cutoff"], "signature"),
        (["cages", ice, "--type", "1", "--signature", "512", "--graph", "cutoff"], "signature"),
        (["--family", "ionicLiquid", "read", tip, "--type", "2"], "nop"),
        (["cages", cubic, "--type", "1", "--graph", "cutoff"], "cubic 4096"),
        (["cages", "input/traj/genice_sI.lammpstrj", "--type", "1", "--signature", "512",
          "--guest-types", "1", "--guest-radius", "2"], "occ 0:2"),
    ]
    hinted = run(seams, ["--type", "-1", "read", three])
    if hinted.returncode != 0 or "more than two" not in hinted.stderr:
        sys.stderr.write("three types\n" + hinted.stdout + hinted.stderr)
        return 1
    for args, needle in checked:
        result = run(seams, args)
        if result.returncode != 0 or needle not in result.stdout:
            sys.stderr.write("miss " + " ".join(args) + "\n" + result.stdout + result.stderr)
            return 1
    for args in ok:
        result = run(seams, args)
        if result.returncode != 0:
            sys.stderr.write("fail " + " ".join(args) + "\n" + result.stderr)
            failed = True
    bad = [
        ["nope", tip],
        ["chill-plus", tip, "--rule", "nope"],
        ["chill", tip, "--graph", "banana"],
        ["--format", "yaml", "read", tip],
        ["rdf", tip, "--bins", "-1"],
        ["rdf", tip, "--pair-types", "2,x"],
        ["density-z", tip, "--axis", "q"],
        ["domains", tip, "--subset", "nope", "--site", "1=tail"],
        ["--strict-input", "read", "missing-file.lammpstrj"],
        ["--strict-input", "--frame", "1", "--last", "2", "read", gap],
        ["--strict-input", "--frame", "1", "--last", "3", "read", run_path],
        ["fingerprint", tip, "--type", "2", "--graph", "seedex"],
        ["fingerprint", tip, "--type", "2", "--graph", "seede"],
        ["fingerprint", tip, "--type", "2", "--graph", "seededX"],
        ["chill", tip, "--type", "2", "--graph", "seedex"],
        ["chill", tip, "--type", "2", "--graph", "seede"],
        ["chill", ice, "--type", "1", "--graph", "seededX"],
        ["domains", tip],
        ["pairs", tip],
        ["ions", tip],
        ["--cage-edge", "nope", "cages", ice, "--type", "1"],
        ["--family", "not-a-family", "read", tip],
        ["--family", "moltenSalt", "chill", tip, "--type", "2"],
        ["--site", "9=nope", "pairs", tip],
        ["cn", tip, "--ions"],
        [],
        ["read"],
        ["--frame", "zz", "read", tip],
        ["--axis", "   ", "density-z", tip],
        ["--axis", "yy", "density-z", tip],
        ["--axis", "yY", "density-z", tip],
        ["--axis", "zz", "density-z", tip],
        ["--axis", "11", "density-z", tip],
        ["--site", "nocolon", "pairs", tip],
        ["--site", "abc=cationHead", "pairs", tip],
        ["--per-atom", os.path.join(work, "no-such-dir", "out.lammpstrj"),
         "fingerprint", ice, "--type", "1", "--hops", "1"],
        ["chill", ice, "--type", "1", "--graph", "knn-union"],
        ["--pair-types", ",", "rdf", tip],
        ["--pair-types", "2,", "cn", tip],
        ["--site", "1=", "pairs", tip],
        ["--site", "=anion", "pairs", tip],
    ]
    ok.extend([
        ["--axis", "Z", "density-z", tip, "--type", "2"],
        ["--axis", "2", "density-z", tip],
        ["--axis", "Y", "density-z", tip],
        ["ovito_chill_plus", ice, "--type", "1"],
        ["ovito-chill+", ice, "--type", "1", "--rule", "chill+"],
        ["--site", "family=waterIce,1=waterO,2=waterH", "pairs", tip],
        ["--site", "1=unspecified,2=donorH,3=acceptor,4=solvent", "pairs", tip],
        ["--water-types", "1,,2", "chill", tip, "--type", "2"],
        ["cages", "input/traj/clathrate-thf.lammpstrj", "--type", "1", "--guest-types", "2"],
    ])
    emitted = run(seams, ["--emit-library", "hex", "fingerprint", ice, "--type", "1", "--hops", "1"])
    if emitted.returncode != 0 or not emitted.stdout.strip():
        sys.stderr.write("emit library failed\n" + emitted.stderr)
        return 1
    lib_a = os.path.join(work, "hex.txt")
    lib_b = os.path.join(work, "deep.txt")
    with open(lib_a, "w", encoding="ascii") as handle:
        handle.write(emitted.stdout)
    deeper = run(seams, ["--emit-library", "deep", "fingerprint", ice, "--type", "1", "--hops", "2"])
    if deeper.returncode != 0:
        sys.stderr.write("emit deep library failed\n" + deeper.stderr)
        return 1
    with open(lib_b, "w", encoding="ascii") as handle:
        handle.write(deeper.stdout)
    matched = [
        ["--library", lib_a, "fingerprint", ice, "--type", "1", "--hops", "1"],
        ["--library", lib_a, "--per-atom", per_atom, "fingerprint", ice, "--type", "1", "--hops", "1"],
        ["--library", lib_a + "," + lib_b, "fingerprint", ice, "--type", "1"],
        ["--library", lib_a + ",,", "fingerprint", ice, "--type", "1", "--hops", "1"],
    ]
    for args in matched:
        result = run(seams, args)
        if result.returncode != 0:
            sys.stderr.write("fail " + " ".join(args) + "\n" + result.stderr)
            return 1
    miss = os.path.join(work, "miss.txt")
    with open(lib_a, encoding="ascii") as handle:
        header = handle.readline()
    with open(miss, "w", encoding="ascii") as handle:
        handle.write(header)
        handle.write("deadbeef nomatch\n")
    unmatched = run(seams, ["--library", miss, "--per-atom", per_atom,
                            "fingerprint", ice, "--type", "1", "--hops", "1"])
    if unmatched.returncode != 0 or "unmatched" not in unmatched.stdout:
        sys.stderr.write("unmatched library\n" + unmatched.stdout + unmatched.stderr)
        return 1
    blank = os.path.join(work, "blank.txt")
    with open(blank, "w", encoding="ascii") as handle:
        handle.write(header.replace(" hops 1 ", " hops 3 "))
    deep_match = run(seams, ["--library", lib_a + "," + blank, "--per-atom", per_atom,
                             "fingerprint", ice, "--type", "1"])
    if deep_match.returncode != 0 or "depth" not in deep_match.stdout:
        sys.stderr.write("depth library\n" + deep_match.stdout + deep_match.stderr)
        return 1
    mixed_depth = run(seams, ["--library", lib_a + "," + blank,
                                "fingerprint", tip, "--type", "2"])
    if mixed_depth.returncode != 0 or "depth" not in mixed_depth.stdout:
        sys.stderr.write("mixed depth\n" + mixed_depth.stdout + mixed_depth.stderr)
        return 1
    ion_is_water = run(seams, ["ions", tip, "--type", "1", "--ion-types", "1,2",
                               "--site", "1=cationHead,2=anion"])
    if ion_is_water.returncode != 0:
        sys.stderr.write("ion type is the water type\n" + ion_is_water.stderr)
        return 1
    env = os.environ.copy()
    env["CHILL_PLUS_RULE"] = "ovito-chill+"
    ruled = run(seams, ["chill-plus", ice, "--type", "1"], env)
    if ruled.returncode != 0:
        sys.stderr.write("env rule failed\n" + ruled.stderr)
        return 1
    for args in (
        ["read", "missing-cli.pdb"],
        ["read", "missing-cli.gro"],
        ["--frame", "1", "--last", "2", "read", "missing-cli.pdb"],
        ["--frame", "1", "--last", "2", "read", "missing-cli.dcd"],
    ):
        probed = run(seams, args)
        if probed.returncode != 0:
            sys.stderr.write("aborted " + " ".join(args) + " " + str(probed.returncode) + "\n")
            return 1
    off = os.environ.copy()
    off["SEAMS_OFFLOAD"] = "1"
    reached = run(seams, ["cages", ice, "--type", "1", "--graph", "cutoff"], off)
    if reached.returncode != 0 or "nop" not in reached.stdout:
        sys.stderr.write("offload cages aborted\n" + reached.stdout)
        return 1
    weird = run(seams, ["--emit-library", 'say "hi"\\\t', "fingerprint", ice,
                        "--type", "1", "--hops", "1"])
    if weird.returncode == 0 and weird.stdout.strip():
        weird_path = os.path.join(work, "weird.txt")
        with open(weird_path, "w", encoding="utf-8") as handle:
            handle.write(weird.stdout)
        escaped = run(seams, ["--format", "json", "--library", weird_path,
                              "fingerprint", ice, "--type", "1", "--hops", "1"])
        if escaped.returncode != 0:
            sys.stderr.write("json escape aborted " + str(escaped.returncode) + "\n")
            return 1
    missing_lib = run(seams, ["--library", os.path.join(work, "absent.txt"),
                              "fingerprint", ice, "--type", "1"])
    if missing_lib.returncode == 0:
        sys.stderr.write("missing library was accepted\n")
        failed = True
    label = "say\"hi\\\t\r\x01"
    weird_json = run(seams, ["--format", "json", "--emit-library", label,
                             "fingerprint", ice, "--type", "1", "--hops", "1"])
    if weird_json.returncode != 0 or "\\u00" not in weird_json.stdout:
        sys.stderr.write("json escape did not print the control character\n")
        failed = True
    interface = os.path.join(work, "interface.lammpstrj")
    ions = os.path.join(work, "ions.lammpstrj")
    gro = os.path.join(work, "two.gro")
    pdb = os.path.join(work, "one.pdb")
    gap = os.path.join(work, "gap.lammpstrj")
    _write_interface(interface)
    _write_ions(ions)
    _write_gro(gro)
    _write_pdb(pdb)
    _write_gap(gap)
    bare = "input/baretraj"
    with open(ice, encoding="ascii") as src, open(bare, "w", encoding="ascii") as dst:
        dst.write(src.read())
    for args in (
        ["read", "input/traj/agi_water_tiny.lammpstrj"],
        ["chill-plus", interface, "--type", "1"],
        ["chill", interface, "--type", "1"],
        ["--water-types", "1", "chill-plus", interface, "--type", "1", "--graph", "mutual"],
        ["cn", ions, "--ions", "--site", "1=cationHead,2=anion", "--cutoff", "6"],
        ["pairs", ions, "--site", "1=cationHead,2=anion"],
        ["pairs", ions, "--site", "1=unspecified"],
        ["pairs", ions, "--site", "1=donorH"],
        ["pairs", ions, "--site", "1=acceptor"],
        ["pairs", ions, "--site", "1=waterO"],
        ["pairs", ions, "--site", "1=waterH"],
        ["pairs", ions, "--site", "1=solvent"],
        ["--frame", "1", "--last", "2", "read", "missing-range.gro"],
        ["--frame", "1", "--last", "2", "read", "missing-range.dcd"],
        ["--frame", "5", "--last", "6", "read", "input/xyz/configuration.xyz"],
        ["read", "input/baretraj"],
        ["--frame", "1", "--last", "2", "read", "input/baretraj"],
        ["--type", "0", "read", "missing-type.pdb"],
        ["--type", "2", "read", "missing-type.pdb"],
        ["--type", "-1", "read", "missing-type.gro"],
        ["f4", ions, "--type", "0", "--htype", "1"],
        ["rdf", ions, "--type", "0"],
        ["cn", ions, "--type", "0"],
        ["cages", ice, "--type", "1", "--signature", "512", "--guest-types", "9", "--guest-radius", "0.1"],
        ["ions", ions, "--ion-types", "1", "--site", "1=cationHead,2=anion", "--type", "2"],
        ["cages", interface, "--type", "1", "--signature", "512", "--graph", "cutoff",
         "--water-types", "1"],
        ["cages", ice, "--type", "1", "--guest-types", "1", "--shape"],
        ["read", bare],
        ["--frame", "1", "--last", "1", "read", "input/xyz/configuration.xyz"],
        ["chill", ice, "--type", "1", "--graph", "seeded"],
        ["chill-plus", empty, "--type", "1"],
        ["f4", ice, "--type", "1", "--htype", "1"],
        ["f4", ice, "--type", "0", "--htype", "9"],
        ["rdf", tip, "--pair-types", " 2 , 2 "],
        ["--water-types", ",,", "chill", tip, "--type", "2"],
        ["--water-types", "1", "ions", tip, "--ion-types", "1",
         "--site", "1=cationHead,2=anion"],
        ["--site", " , 1=cationHead", "pairs", tip],
        ["--site", "family=ionicLiquid,1=cationHead,2=anion", "pairs", tip],
        ["--site", "1=tail,2=polar,3=apolar", "pairs", "input/traj/agi_water_tiny.lammpstrj"],
        ["cages", ice, "--type", "1", "--graph", "knn-union", "--signature", "512"],
        ["--format", "json", "chill", empty, "--type", "1"],
        ["--jobs", "1", "--frame", "1", "--last", "2", "read", cubic],
        ["density-z", tip, "--type", "2", "--bins", "0"],
        ["rdf", tip, "--bins", "0", "--type", "2"],
        ["cn", empty, "--type", "1"],
        ["hbonds", ice, "--type", "1", "--htype", "1", "--hdist", "0.1"],
        ["read", "missing-edge.pdb"],
        ["read", "missing-edge.gro"],
        ["read", "missing-edge.dcd"],
        ["--frame", "1", "--last", "2", "read", "missing-edge.pdb"],
        ["read", gro],
        ["read", pdb],
        ["--frame", "1", "--last", "2", "read", gro],
        ["--frame", "1", "--last", "3", "read", gro],
        ["--frame", "1", "--last", "2", "read", pdb],
        ["--type", "0", "read", pdb],
        ["--type", "2", "read", pdb],
        ["--type", "-1", "read", gro],
        ["--format", "json", "--frame", "1", "--last", "2", "read", gro],
        ["--jobs", "2", "--frame", "1", "--last", "2", "read", "input/xyz/configuration.xyz"],
        ["rdf", empty, "--type", "0"],
        ["rdf", empty, "--type", "-1"],
        ["rdf", "missing-volume.lammpstrj", "--type", "0", "--bins", "2", "--cutoff", "1"],
        ["density-z", pdb, "--bins", "0"],
        ["density-z", pdb, "--type", "0", "--bins", "0", "--axis", "x"],
        ["--water-types", "1,2", "f4", tip, "--type", "0", "--htype", "1"],
        ["--water-types", "1", "f4", ice, "--type", "0", "--htype", "1"],
        ["cn", ions, "--ions", "--site", "9=cationHead", "--cutoff", "6"],
        ["cn", ions, "--ions", "--site", "1=cationHead,2=anion", "--cutoff", "0.01"],
        ["ions", ions, "--ion-types", "1", "--type", "1", "--site", "1=cationHead,2=anion"],
        ["ions", ice, "--ion-types", "9", "--type", "1"],
        ["--graph", "", "fingerprint", ice, "--type", "1", "--hops", "1"],
        ["--graph", "seeded", "fingerprint", ice, "--type", "1", "--hops", "1"],
        ["--graph", "", "chill", ice, "--type", "1"],
        ["cages", ice, "--type", "1", "--signature", "512", "--graph", "mutual"],
        ["--guest-types", "1,,2", "cages", ice, "--type", "1", "--signature", "512"],
        ["--ion-types", "1,,2", "ions", tip, "--site", "1=cationHead,2=anion"],
        ["cages", ice, "--type", "1", "--graph", "knn", "--signature", "512"],
        ["fingerprint", ice, "--type", "1", "--graph", "knn", "--hops", "1"],
        ["fingerprint", ice, "--type", "1", "--graph", "knn-union", "--hops", "1"],
        ["--water-types", "1", "cages", tip, "--type", "2", "--per-atom", cages_atom],
        ["ions", ice, "--ion-types", "1", "--type", "1", "--per-atom", per_atom],
        ["cages", ice, "--type", "1", "--signature", "512", "--guest-types", "9",
         "--guest-radius", "1", "--shape"],
        ["--frame", "1", "--last", "2", "chill", gro, "--type", "0"],
        ["hbonds", pdb, "--type", "0", "--htype", "1", "--donors"],
        ["density-z", empty, "--type", "0", "--bins", "0", "--axis", "z"],
        ["cn", pdb, "--type", "0", "--bins", "0"],
        ["rdf", pdb, "--type", "-1", "--bins", "2"],
        ["steinhardt", pdb, "--type", "0"],
        ["f4", pdb, "--type", "-1", "--htype", "8"],
        ["cages", pdb, "--type", "0", "--signature", "512", "--incomplete", "--shape"],
        ["domains", ions, "--site", "1=polar,2=apolar", "--subset", "polar", "--cutoff", "6"],
        ["--axis", "x", "density-z", tip, "--type", "2", "--bins", "0"],
        ["--axis", "y", "density-z", tip, "--type", "2", "--bins", "0"],
        ["--axis", "z", "density-z", tip, "--type", "2", "--bins", "0"],
    ):
        result = run(seams, args)
        if result.returncode != 0:
            sys.stderr.write("aborted " + " ".join(args) + " " + str(result.returncode) + "\n")
            failed = True
    for args in (
        ["read", "missing-edge.con"],
        ["--pair-types", "nope", "rdf", tip],
        ["--site", "zz=cationHead", "pairs", tip],
        ["--cage-edge", "9", "cages", ice, "--type", "1"],
        ["--family", "nope", "read", tip],
        ["--frame", "1", "--last", "2", "read", "missing-range.con"],
        ["--pair-types", "2,zz", "rdf", ice, "--type", "1"],
        ["rdf", ions, "--pair-types", "2,zz"],
        ["rdf", tip, "--pair-types", "2,2,2"],
        ["cn", tip, "--pair-types", ","],
        ["cages", ice, "--type", "1", "--signature", "+"],
        ["--strict-input", "--frame", "1", "--last", "2", "read", gap],
        ["--strict-input", "--frame", "1", "--last", "3", "read", run_path],
        ["fingerprint", tip, "--type", "2", "--graph", "seedex"],
        ["fingerprint", tip, "--type", "2", "--graph", "seede"],
        ["fingerprint", tip, "--type", "2", "--graph", "seededX"],
        ["chill", tip, "--type", "2", "--graph", "seedex"],
        ["chill", tip, "--type", "2", "--graph", "seede"],
        ["chill", ice, "--type", "1", "--graph", "seededX"],
        ["--water-types", "9", "chill", tip, "--type", "2"],
        ["--water-types", "9", "cages", ice, "--type", "1"],
        ["--mol", "nope", "chill", tip, "--type", "2"],
        ["--mol", "0-1", "chill", tip],
    ):
        result = run(seams, args)
        if result.returncode != 2:
            sys.stderr.write("bad exit " + " ".join(args) + " " + str(result.returncode) + "\n")
            failed = True
    for args in bad:
        result = run(seams, args)
        if result.returncode == 0:
            sys.stderr.write("accepted " + " ".join(args) + "\n")
            failed = True
    if os.path.exists("input/baretraj"):
        os.remove("input/baretraj")
    fault_lib = os.path.join(os.path.dirname(os.path.abspath(seams)), "bindings", "libseams_test_fault.so")
    # The fault allocator replaces operator new. An ASan runtime has to
    # stay first in the preload list, so this block stays on the ordinary build.
    if os.path.exists(fault_lib) and "ASAN_OPTIONS" not in os.environ:
        def fault_env(base=None):
            env = (base or os.environ).copy()
            prior = env.get("LD_PRELOAD", "")
            # An instrumented seams needs libasan ahead of the fault allocator.
            env["LD_PRELOAD"] = prior + ":" + fault_lib if prior else fault_lib
            return env

        count_env = fault_env()
        count_env["SEAMS_TEST_KIND"] = "count"
        counted = subprocess.run([seams, "read", ice, "--type", "1"], capture_output=True, text=True, env=count_env)
        serial = 0
        for line in counted.stderr.splitlines():
            if line.startswith("fault-serial "):
                serial = int(line.split()[1])
        if serial < 50 or counted.returncode != 0:
            sys.stderr.write("fault count\n" + counted.stderr)
            return 1

        def serial_of(args):
            env = fault_env()
            env["SEAMS_TEST_KIND"] = "count"
            proc = subprocess.run([seams, *args], capture_output=True, text=True, env=env)
            found = 0
            for line in proc.stderr.splitlines():
                if line.startswith("fault-serial "):
                    found = int(line.split()[1])
            return found, proc.returncode

        def throw_at(args, at):
            env = fault_env()
            env["SEAMS_TEST_KIND"] = "bad_alloc"
            env["SEAMS_TEST_AT"] = str(at)
            return subprocess.run([seams, *args], capture_output=True, text=True, env=env)

        long_token = "n" * 64
        loud_sets = (
            ["rdf", tip, "--pair-types", long_token + ",1"],
            ["pairs", tip, "--site", long_token + "=cationHead"],
        )
        for loud_args in loud_sets:
            env = fault_env()
            env["SEAMS_TEST_KIND"] = "count"
            env["SEAMS_TRACE_SIZE"] = "1"
            traced = subprocess.run([seams, *loud_args], capture_output=True, text=True, env=env)
            sized = []
            for line in traced.stderr.splitlines():
                parts = line.split()
                if len(parts) == 3 and parts[0] == "size-serial" and parts[2] == "65":
                    sized.append(int(parts[1]))
            if not sized:
                sys.stderr.write("fault trace\n" + traced.stderr[-400:])
                return 1
            hit = False
            for at in sized:
                probed = throw_at(loud_args, at)
                if probed.returncode == 2 and "bad_alloc" in probed.stderr:
                    hit = True
            if not hit:
                sys.stderr.write("fault mismatch " + " ".join(loud_args) + "\n")
                return 1
        for kind, needle in (("int", "unknown error"), ("bad_alloc", "bad_alloc")):
            env = fault_env()
            env["SEAMS_TEST_KIND"] = kind
            env["SEAMS_TEST_AT"] = str(serial - 20)
            probed = subprocess.run([seams, "read", ice, "--type", "1"], capture_output=True, text=True, env=env)
            if probed.returncode != 2 or needle not in probed.stderr:
                sys.stderr.write("fault " + kind + " " + str(probed.returncode) + "\n" + probed.stderr)
                return 1
    by_atom = run(seams, ["fingerprint", ice, "--type", "1", "--by-atom", "--hops", "1"])
    if by_atom.returncode != 0 or "atom " not in by_atom.stdout or "coord " not in by_atom.stdout:
        sys.stderr.write("by-atom fingerprint\n" + by_atom.stdout + by_atom.stderr)
        return 1
    class_dump = os.path.join(work, "classes.lammpstrj")
    classed = run(seams, ["fingerprint", ice, "--type", "1", "--hops", "1", "--per-atom", class_dump])
    if classed.returncode != 0 or not os.path.exists(class_dump):
        sys.stderr.write("class dump\n" + classed.stdout + classed.stderr)
        return 1
    cage_dump = os.path.join(work, "cages-atom.lammpstrj")
    cage_atoms = run(seams, ["cages", ice, "--type", "1", "--per-atom", cage_dump])
    if cage_atoms.returncode != 0 or "hexagonal" not in cage_atoms.stdout:
        sys.stderr.write("cage per-atom\n" + cage_atoms.stdout + cage_atoms.stderr)
        return 1
    signed = run(seams, ["cages", ice, "--type", "1", "--signature", "512",
                         "--guest-types", "1", "--guest-radius", "4"])
    if signed.returncode != 0 or "occ " not in signed.stdout:
        sys.stderr.write("guest occupancy\n" + signed.stdout + signed.stderr)
        return 1
    bare_ions = run(seams, ["ions", tip, "--type", "2", "--ion-types", "1",
                            "--site", "1=cationHead,2=anion"])
    if bare_ions.returncode != 0 or "ions " not in bare_ions.stdout:
        sys.stderr.write("ions without water types\n" + bare_ions.stdout + bare_ions.stderr)
        return 1
    typed_ions = run(seams, ["--water-types", "2", "ions", tip, "--type", "2",
                             "--ion-types", "1", "--site", "1=cationHead,2=anion"])
    if typed_ions.returncode != 0 or "ions " not in typed_ions.stdout:
        sys.stderr.write("ions with water types\n" + typed_ions.stdout + typed_ions.stderr)
        return 1
    ranged = run(seams, ["--frame", "1", "--last", "2", "read", cubic])
    if ranged.returncode != 0:
        sys.stderr.write("two-frame read\n" + ranged.stderr)
        return 1
    shaped = run(seams, ["cages", ice, "--type", "1", "--signature", "hc", "--shape"])
    if shaped.returncode != 0 or "off" not in shaped.stdout:
        sys.stderr.write("cage shape\n" + shaped.stdout + shaped.stderr)
        return 1
    seeded = run(seams, ["--graph", "seeded", "chill-plus", ice, "--type", "1"])
    plain = run(seams, ["chill-plus", ice, "--type", "1"])
    if seeded.returncode != 0 or plain.returncode != 0 or seeded.stdout != plain.stdout:
        sys.stderr.write("seeded chill-plus\n" + seeded.stdout + plain.stdout)
        return 1
    axis_y = run(seams, ["--axis", "Y", "density-z", ice, "--type", "1", "--bins", "3"])
    axis_x = run(seams, ["--axis", "x", "density-z", ice, "--type", "1", "--bins", "3"])
    if axis_y.returncode != 0 or "# z rho" not in axis_y.stdout or axis_y.stdout == axis_x.stdout:
        sys.stderr.write("axis Y\n" + axis_y.stdout + axis_x.stdout)
        return 1
    framed = run(seams, ["--frame", "zz", "read", ice])
    if framed.returncode != 2 or "not a number" not in framed.stderr:
        sys.stderr.write("frame zz\n" + framed.stderr)
        return 1
    bad_types = run(seams, ["--pair-types", "2,zz", "rdf", ice, "--type", "1"])
    if bad_types.returncode != 2 or "bad --pair-types" not in bad_types.stderr:
        sys.stderr.write("bad types\n" + bad_types.stderr)
        return 1
    bad_site = run(seams, ["--site", "zz=cationHead", "pairs", ice])
    if bad_site.returncode != 2 or "not an integer" not in bad_site.stderr:
        sys.stderr.write("bad site\n" + bad_site.stderr)
        return 1
    return 1 if failed else 0


if __name__ == "__main__":
    raise SystemExit(main())
