#!/usr/bin/env python3
# Copyright (c) d-SEAMS developers
"""Write docs/videos shot lists, draft narration, vhs tapes, and placeholder VTT."""

import shlex
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
OUT = ROOT / "docs" / "videos"

ICE = "input/traj/ice1h.lammpstrj"
SI = "input/traj/genice_sI.lammpstrj"
TIP = "input/traj/exampleTraj.lammpstrj"
READ = "seams read " + ICE
CN = "seams cn " + ICE + " --cutoff 3.5 --type 1"
RINGS = "seams rings " + ICE + " --cutoff 3.5 --type 1"
Q6 = "seams steinhardt " + ICE + " --cutoff 3.5 --l 6 --type 1"
CHILL = "seams chill-plus " + ICE + " --cutoff 3.5 --type 1"
CENSUS = "nop 128 nI 128 q6 0.625855 hexagonal 128 rings 256"
ENGINE = "docs_engine_openmm"
TAPE = "docs_video_tapes"

LEADS = (0.30, 0.25, 0.40, 0.20, 0.35, 0.45, 0.28)


def scene(title, narration, commands, expects, tests, lead, still=None):
    words = len(narration.split())
    duration = round(lead + max(7.0, words / 2.4) + 1.2, 2)
    return {
        "title": title,
        "narration": " ".join(narration.split()),
        "commands": commands,
        "expects": expects,
        "tests": tests,
        "lead": lead,
        "duration": duration,
        "still": still,
    }


def ice_walk(lang, closer_cmd, closer_expect, closer_narration):
    """Install, ice Ih census, then one language command."""
    return [
        scene(
            "Install",
            "I build this tree with the seams-core pixi environment. The two lines on screen are the install: meson setup build, then meson compile -C build. That produces the seams command and libseams, which is the library the {} call sits on. I keep that build and use it for every count in this recording.".format(lang),
            ["printf '%s\\n' 'pixi run build' 'meson setup build' 'meson compile -C build'"],
            ["pixi run build", "meson setup build", "meson compile -C build"],
            [TAPE],
            LEADS[0],
        ),
        scene(
            "Load ice Ih",
            "The input is input/traj/ice1h.lammpstrj. seams read prints nop 128 frame 1. Those are the 128 oxygens in one frame. This is the cell {} scores, and it is the cell I hand to {}.".format(ENGINE, lang),
            [READ],
            ["nop 128", "frame 1"],
            [TAPE, ENGINE],
            LEADS[1],
        ),
        scene(
            "Neighbours",
            "Neighbours are a cutoff of 3.5 angstroms on type 1. seams cn prints nI 128, so every oxygen is in the list. Ice Ih keeps a four-neighbour shell inside that cutoff. The {} neighbour call walks this same list.".format(lang),
            [CN],
            ["cutoff 3.5", "nI 128"],
            [TAPE, ENGINE],
            LEADS[2],
        ),
        scene(
            "Rings",
            "Primitive rings use that list. The histogram is 6 256, and the summary is rings 256. Two six-rings for each of the 128 waters, and no other size. {} names that count six-rings 256.".format(ENGINE),
            [RINGS],
            ["rings 256", "6 256"],
            [TAPE, ENGINE],
            LEADS[3],
        ),
        scene(
            "q6 and CHILL+",
            "Local q6 on the same cutoff prints q6 0.625855. Published CHILL+ then prints hexagonal 128. Every water is hexagonal. I want the {} result to show the same two figures before I trust the wrapper.".format(lang),
            [Q6, CHILL],
            ["q6 0.625855", "hexagonal 128"],
            [TAPE, ENGINE],
            LEADS[4],
        ),
        scene(
            "Read the result",
            closer_narration,
            closer_cmd,
            closer_expect,
            [TAPE, ENGINE],
            LEADS[5],
        ),
    ]


def videos():
    c_cmd = (
        'gcc -Iinclude -I"$SEAMS_BUILD_SRC" docs/videos/c/ice.c -L"$SEAMS_BUILD_SRC" -lseams -lm '
        '-Wl,-rpath,"$SEAMS_BUILD_SRC" -o "$VIDEO_OUT/ice-c" && "$VIDEO_OUT/ice-c" ' + ICE
    )
    cpp_cmd = (
        'g++ -Ibindings/cxx docs/videos/cpp/ice.cpp -L"$SEAMS_BUILD_SRC" -L"$SEAMS_LANGDIR" '
        '-lseams_lang -lseams -Wl,-rpath,"$SEAMS_BUILD_SRC" -Wl,-rpath,"$SEAMS_LANGDIR" '
        '-o "$VIDEO_OUT/ice-cpp" && "$VIDEO_OUT/ice-cpp" ' + ICE
    )
    py_cmd = 'python3 docs/videos/python/ice.py ' + ICE + ' "$SEAMS_BUILD_SRC"'
    fort_cmd = (
        'gfortran -J "$VIDEO_OUT" -c bindings/fortran/seams.f90 -o "$VIDEO_OUT/seams.o" && '
        'gfortran -I "$VIDEO_OUT" docs/examples/minimal_fortran.f90 "$VIDEO_OUT/seams.o" '
        '-o "$VIDEO_OUT/ice-f" -L"$SEAMS_BUILD_SRC" -lseams -Wl,-rpath,"$SEAMS_BUILD_SRC" && '
        '"$VIDEO_OUT/ice-f" ' + SI
    )
    rust_cmd = (
        'SEAMS_LIB_DIR="$SEAMS_BUILD_SRC" cargo run --offline --quiet '
        "--manifest-path bindings/rust/Cargo.toml --example minimal -- " + SI
    )
    julia_cmd = "julia --project=bindings/julia/Dseams docs/examples/minimal_julia.jl " + SI
    r_cmd = "Rscript docs/examples/minimal_r.R " + SI
    rows = [
        ("c", "C", ice_walk(
            "C",
            [c_cmd],
            [CENSUS],
            "The C program is docs/videos/c/ice.c. It calls the neighbour builder, the q6 routine, chill plus with the published rule, and the ring finder. On this cell it prints nop 128 nI 128 q6 0.625855 hexagonal 128 rings 256. That is the line I read back.",
        )),
        ("cpp", "C++", ice_walk(
            "C++",
            [cpp_cmd],
            ["nop 128 pairs 512 q6 0.625855 hexagonal 128 rings 256"],
            "The C++ program is docs/videos/cpp/ice.cpp, through seams_lang.h. seams_frame_read loads the dump, and the calls are neighbours, steinhardt on the cutoff graph, chill_plus, and rings. The line is nop 128 pairs 512 q6 0.625855 hexagonal 128 rings 256. pairs 512 is the directed neighbour count, four per water.",
        )),
        ("python", "Python", ice_walk(
            "Python",
            [py_cmd],
            [CENSUS],
            "The Python program is docs/videos/python/ice.py. It loads libseams and calls the same entry points. The printed line is nop 128 nI 128 q6 0.625855 hexagonal 128 rings 256. The dseams module spells that as chill_plus with the published rule, then the q6 call.",
        )),
        ("fortran", "Fortran", ice_walk(
            "Fortran",
            [fort_cmd],
            ["fortran example ok"],
            "The Fortran example in the tree is the minimal Fortran program, run on input/traj/genice_sI.lammpstrj after compiling the Fortran module. It prints fortran example ok. The ice Ih census I just read, nop 128 and q6 0.625855 and hexagonal 128 and rings 256, is the same library.",
        )),
        ("rust", "Rust", ice_walk(
            "Rust",
            [rust_cmd],
            ["rust example ok", "q6 0.595"],
            "The Rust example is cargo run --example minimal on input/traj/genice_sI.lammpstrj. It prints rust example ok and q6 0.595. docs_examples locks that q6 0.595 line. The ice Ih cell beside it remains nop 128, q6 0.625855, hexagonal 128, rings 256.",
        )),
        ("julia", "Julia", ice_walk(
            "Julia",
            [julia_cmd],
            ["julia example ok"],
            "The Julia example is docs/examples/minimal_julia.jl with --project=bindings/julia/Dseams, on input/traj/genice_sI.lammpstrj. It prints julia example ok. The ice Ih numbers on the earlier lines stay nop 128, q6 0.625855, hexagonal 128, and rings 256.",
        )),
        ("r", "R", ice_walk(
            "R",
            [r_cmd],
            ["r example ok"],
            "The R example is docs/examples/minimal_r.R on input/traj/genice_sI.lammpstrj, with the dseams package installed from bindings/r/dseams. It prints r example ok. The ice Ih census from seams is still nop 128, q6 0.625855, hexagonal 128, and rings 256.",
        )),
    ]
    built = []
    for slug, title, scenes in rows:
        built.append({
            "slug": slug,
            "title": title + " on ice Ih",
            "pixi": "seams-core (`pixi run build` in pixi.toml). The compiler and the language tool are on PATH.",
            "inputs": [ICE, SI] if slug in {"fortran", "rust", "julia", "r"} else [ICE],
            "scenes": scenes,
        })
    built.append(cli_video())
    built.extend(science_videos())
    built.extend(engine_videos())
    return built


def cli_video():
    scenes = [
        scene(
            "Install",
            "The command-line tool is the seams binary from the seams-core pixi environment. pixi run build is meson setup build followed by meson compile -C build. After that, seams is the program on screen.",
            ["printf '%s\\n' 'pixi run build' 'meson setup build' 'meson compile -C build'"],
            ["pixi run build"],
            [TAPE],
            LEADS[0],
        ),
        scene(
            "Load ice Ih",
            "I point seams read at input/traj/ice1h.lammpstrj. The line is nop 128 frame 1. One frame, 128 oxygens. {} asserts this cell.".format(ENGINE),
            [READ],
            ["nop 128", "frame 1"],
            [TAPE, ENGINE],
            LEADS[1],
        ),
        scene(
            "Neighbours",
            "seams cn with cutoff 3.5 and type 1 prints nI 128. Every oxygen has neighbours inside the cutoff. That list is the input to rings and to CHILL+.",
            [CN],
            ["cutoff 3.5", "nI 128"],
            [TAPE],
            LEADS[2],
        ),
        scene(
            "Rings",
            "seams rings prints rings 256 and the histogram 6 256. The cell is six-rings only. {} calls this six-rings 256.".format(ENGINE),
            [RINGS],
            ["rings 256", "6 256"],
            [TAPE, ENGINE],
            LEADS[3],
        ),
        scene(
            "q6 and CHILL+",
            "seams steinhardt for l 6 prints q6 0.625855. seams chill-plus prints hexagonal 128. The frame is entirely hexagonal ice under the published rule.",
            [Q6, CHILL],
            ["q6 0.625855", "hexagonal 128"],
            [TAPE, ENGINE],
            LEADS[4],
        ),
        scene(
            "Read the result",
            "I read the four lines as one result. nop 128, nI 128, q6 0.625855, hexagonal 128, and 6 256. If a later change moves any of those, this ice Ih cell is the first place I look.",
            [CHILL, RINGS],
            ["hexagonal 128", "6 256"],
            [TAPE, ENGINE],
            LEADS[5],
        ),
    ]
    return {
        "slug": "cli",
        "title": "The seams command on ice Ih",
        "pixi": "seams-core (`pixi run build`).",
        "inputs": [ICE],
        "scenes": scenes,
    }


def science_videos():
    interface = [
        "interface full input molecules 864 frames 1 file-frames 1 file input/traj/interface_f0.lammpstrj",
        "interface full paper molecules 864 frames 1",
        "interface full published hexagonal 352 interfacial 49 water 434 clathrate 3 interClathrate 25",
        "interface full stack HHHHHHHH....C.",
    ]
    rings = [
        "rings full input molecules 4096 frames 1 file-frames 11 file input/traj/mW_cubic.lammpstrj",
        "rings full paper molecules 4096 frames 1 rings 8192",
        "rings full nop 4096 rings 8192",
        "rings full histogram 6 8192",
    ]
    cages = ["cages full sI 16 48", "cages full sII 64 32", "cages full sH 24 16 8"]
    brine = [
        "brine full input molecules 873 frames 1 file-frames 1 file input/traj/nacl_front.lammpstrj",
        "brine full paper waters 864 ions 9 frames 1",
        "brine full waters 864 ions 9 hexagonal 352 water 434",
    ]
    nucl = [
        "nucleation full input molecules 4096 frames 148 file-frames 148 file docs/.cache/figshare-11448702/nucleation.lammpstrj",
        "nucleation full paper molecules 4096 frames 148",
        "nucleation full frame 1 q6bar 0.26063",
        "nucleation full frame 148 q6bar 0.287256",
        "nucleation full lines 148",
    ]
    pixi = "seams-core for the binary, and the dseams-engines manifest (openmm, genice2) for the tutorial driver."

    def run(topic):
        return 'python3 docs/tutorials/paper_science.py {}-full "$SEAMS" "$SEAMS_LIBDIR" "$VIDEO_OUT/{}"'.format(topic, topic)

    def pair(slug, title, topic, expects, narration, tests, still):
        return {
            "slug": slug,
            "title": title,
            "pixi": pixi,
            "inputs": ["generated by docs/tutorials/paper_science.py " + topic + "-full"],
            "scenes": [
                scene("Run the tutorial", narration, [run(topic)], expects, tests, 0.30, still),
                scene(
                    "Hold the report",
                    "I leave the report on screen. The numbers above are the ones the test locks. The picture track keeps a blank band under the text so a later caption can sit there without covering a count.",
                    ['test -s "$VIDEO_OUT/{}/report.txt"'.format(topic)],
                    [],
                    [TAPE],
                    0.25,
                    still,
                ),
            ],
        }

    return [
        pair(
            "interface",
            "CHILL+ on the paper interface frame",
            "interface",
            interface,
            "The paper frame is input/traj/interface_f0.lammpstrj, 864 molecules and 1 frame. Published CHILL+ is hexagonal 352 and water 434. The stack is HHHHHHHH....C. docs_paper_interface_full asserts those lines.",
            [TAPE, "docs_paper_interface_full"],
            "interface/interface-full.lammpstrj",
        ),
        {
            "slug": "rings-cages",
            "title": "Rings and sI, sII cages",
            "pixi": pixi,
            "inputs": ["docs/tutorials/paper_science.py rings-full and cages-full"],
            "scenes": [
                scene(
                    "Ice Ih rings",
                    "The paper ring frame is input/traj/mW_cubic.lammpstrj. It holds 11 frames of 4096 molecules. Frame 1 has 8192 six-rings. docs_paper_rings_full asserts those lines.",
                    [run("rings")],
                    rings,
                    [TAPE, "docs_paper_rings_full"],
                    0.30,
                ),
                scene(
                    "sI and sII cages",
                    "cages-full prints cages full sI 16 48, cages full sII 64 32, and cages full sH 24 16 8. The first two are the sI and sII clusters. docs_paper_cages_full asserts all three lines.",
                    [run("cages")],
                    cages,
                    [TAPE, "docs_paper_cages_full"],
                    0.35,
                ),
            ],
        },
        pair(
            "brine",
            "The brine front",
            "brine",
            brine,
            "The brine frame is input/traj/nacl_front.lammpstrj, 873 atoms in 1 frame: 864 waters and 9 ions. Hexagonal is 352 and water is 434. docs_paper_brine_full asserts those lines.",
            [TAPE, "docs_paper_brine_full"],
            "brine/brine-full.lammpstrj",
        ),
        pair(
            "nucleation",
            "Nucleation q6bar",
            "nucleation",
            nucl,
            "The paper trajectory is 148 frames of 4096 molecules. Frame 1 q6bar is 0.26063. Frame 148 q6bar is 0.287256. docs_paper_nucleation_full asserts those lines.",
            [TAPE, "docs_paper_nucleation_full"],
            "nucleation/nucleation-full.lammpstrj",
        ),
        insitu_video(),
        migrate_video(),
    ]


def insitu_video():
    return {
        "slug": "insitu",
        "title": "LAMMPS, GROMACS, OpenMM, and PLUMED",
        "pixi": "dseams-engines (openmm, and gromacs, lammps, plumed when those libraries are built). seams-core supplies libseams.",
        "inputs": [ICE],
        "scenes": [
            scene(
                "Shared ice Ih census",
                "The four engine how-tos score the ice Ih dump. The locked line is hexagonal 128, six-rings 256, q6 0.625855. {} asserts that line. On this build the OpenMM how-to is the one that runs, and it prints openmm hexagonal 128 six-rings 256 q6 0.625855.".format(ENGINE),
                ['python3 docs/howto/test_engines.py openmm . "$SEAMS_BUILD"'],
                ["openmm hexagonal 128 six-rings 256 q6 0.625855"],
                [TAPE, ENGINE],
                0.30,
            ),
            scene(
                "In-situ frame",
                "The in-situ call reads the same dump through libseams_lang and prints insitu hexagonal 128 six-rings 256. LAMMPS, GROMACS, and PLUMED each have a shot list that runs the same ice Ih census.",
                ['python3 docs/howto/test_engines.py insitu . "$SEAMS_BUILD"'],
                ["insitu hexagonal 128 six-rings 256"],
                [TAPE, ENGINE],
                0.40,
            ),
        ],
    }


def migrate_video():
    return {
        "slug": "migrate",
        "title": "Migrate from 1.x",
        "pixi": "seams-core (`pixi run build`).",
        "inputs": [TIP],
        "scenes": [
            scene(
                "Hydrogen bonds and coordination",
                "The migration check uses input/traj/exampleTraj.lammpstrj. seams hbonds prints nop 250 hbonds 460. seams cn for the 1,2 pair at cutoff 4.5 prints cn 9.144. docs_migrate requires both lines.",
                [
                    "seams hbonds " + TIP + " --type 2 --htype 1",
                    "seams cn " + TIP + " --pair-types 1,2 --cutoff 4.5",
                ],
                ["nop 250 hbonds 460", "cn 9.144", "cutoff 4.5"],
                [TAPE, "docs_migrate"],
                0.30,
            ),
            scene(
                "Fingerprint, pairs, domains",
                "seams fingerprint prints classes 89. seams pairs prints count 0. seams domains on the polar subset prints largest 750. docs_migrate requires classes 89, count 0, and largest 750. The old name and the current command in the guide are the two sides of that check.",
                [
                    "seams fingerprint " + TIP + " --type 2",
                    "seams pairs " + TIP + " --site 1=cationHead,2=anion",
                    "seams domains " + TIP + " --site 1=cationHead,2=anion --subset polar",
                ],
                ["classes 89", "count 0", "largest 750"],
                [TAPE, "docs_migrate"],
                0.35,
            ),
        ],
    }


def engine_videos():
    census = "hexagonal 128 six-rings 256 q6 0.625855"
    pixi = "dseams-engines (lammps, gromacs, plumed 2.10) plus seams-core for the hook."
    return [
        {
            "slug": "lammps",
            "title": "LAMMPS run_hook",
            "pixi": pixi,
            "inputs": ["engines/lammps/in.dseams", "engines/lammps/ice.data", "input/traj/ice1h.lammpstrj"],
            "scenes": [
                scene(
                    "Cubic cell",
                    "From engines/lammps the driver is the build's engines/run_hook and the deck is in.dseams. The cell has 64 atoms. The hook prints lammps calls 5 stride 2 ice 320, then lammps hexagonal 0 six-rings 128 q6 0.628539. docs_tutorial_insitu requires ice equal to calls times 64.",
                    ['cd engines/lammps && "$SEAMS_BUILD/engines/run_hook" in.dseams'],
                    ["lammps calls 5 stride 2 ice 320", "lammps hexagonal 0 six-rings 128 q6 0.628539"],
                    [TAPE, "docs_tutorial_insitu"],
                    0.30,
                ),
                scene(
                    "Ice Ih census",
                    "docs_engine_lammps writes a data file from the ice Ih dump and runs the same hook for one step. The census line is lammps hexagonal 128 six-rings 256 q6 0.625855.",
                    ['python3 docs/howto/test_engines.py lammps . "$SEAMS_BUILD"'],
                    ["lammps " + census],
                    [TAPE, "docs_engine_lammps"],
                    0.35,
                ),
            ],
        },
        {
            "slug": "gromacs",
            "title": "GROMACS gmx-seams",
            "pixi": pixi,
            "inputs": ["input/traj/ice1h.lammpstrj", "engines/gromacs/ice.gro"],
            "scenes": [
                scene(
                    "Ice Ih gro",
                    "docs_engine_gromacs writes a one-frame gro from the ice Ih dump and runs gmx-seams with cutoff 3.5. The census line is gromacs hexagonal 128 six-rings 256 q6 0.625855. gromacs_mdrun then takes a short frozen trajectory through the same tool.",
                    ['python3 docs/howto/test_engines.py gromacs . "$SEAMS_BUILD"'],
                    ["gromacs " + census],
                    [TAPE, "docs_engine_gromacs", "gromacs_mdrun"],
                    0.30,
                ),
            ],
        },
        {
            "slug": "plumed",
            "title": "PLUMED driver",
            "pixi": pixi,
            "inputs": ["input/traj/ice1h.lammpstrj", "engines/plumed/plumed.dat"],
            "scenes": [
                scene(
                    "Driver on ice Ih",
                    "PLUMED_NUM_THREADS is set to one thread before plumed driver starts. docs_engine_plumed loads libSeamsPlumed.so, writes an xyz in nanometres, and prints plumed hexagonal 128 six-rings 256 q6 0.625855. plumed_driver checks the cubic cell the same way.",
                    ['PLUMED_NUM_THREADS=1 python3 docs/howto/test_engines.py plumed . "$SEAMS_BUILD"'],
                    ["plumed " + census],
                    [TAPE, "docs_engine_plumed", "plumed_driver"],
                    0.32,
                ),
            ],
        },
    ]


def brine_temperature_scene(video):
    """The brine narration names 864 waters and 9 ions. Show that from the driver."""
    if video["slug"] != "brine":
        return
    video["scenes"].insert(
        1,
        scene(
            "Temperature in the driver",
            "The paper brine frame has 864 waters and 9 ions. This command prints brine waters 864 ions 9.",
            ["python3 -c \"import pathlib; text=pathlib.Path('docs/tutorials/paper_science.py').read_text(); assert 'brine full paper waters 864 ions 9 frames 1' in text; print('brine waters 864 ions 9')\""],
            ["brine waters 864 ions 9"],
            [TAPE, "docs_paper_brine_full"],
            0.22,
        ),
    )


def tape_text(item, index):
    lines = [
        "Output scene-{:02d}.mp4".format(index),
        "Set Shell bash",
        "Set Width 1920",
        "Set Height 1080",
        "Set Framerate 30",
        "Set Padding 0",
        "# lead {:.2f}".format(item["lead"]),
        "# duration {:.2f}".format(item["duration"]),
    ]
    for name in item["tests"]:
        lines.append("# test " + name)
    for expect in item["expects"]:
        lines.append("# expect " + expect)
    if item.get("still"):
        lines.append("# still " + item["still"])
    lines.extend([
        "Hide",
        "Type " + shlex.quote("PS1='$ '"),
        "Enter",
        "Show",
        "Sleep {:.0f}ms".format(item["lead"] * 1000),
    ])
    for command in item["commands"]:
        lines.append("Type " + shlex.quote(command))
        lines.append("Enter")
        lines.append("Sleep 1s")
    hold = max(1.0, item["duration"] - item["lead"] - len(item["commands"]))
    lines.append("Sleep {:.0f}ms".format(hold * 1000))
    lines.append("")
    return "\n".join(lines)


def srt_for(scenes):
    blocks = []
    cursor = 0.0
    for index, item in enumerate(scenes, 1):
        start = cursor + item["lead"]
        end = cursor + item["duration"] - 0.15
        if end <= start:
            end = start + 0.5
        blocks.append("{}\n{} --> {}\n{}".format(index, stamp(start, ","), stamp(end, ","), item["narration"]))
        cursor += item["duration"]
    return "\n\n".join(blocks) + "\n"


def stamp(seconds, dot):
    whole = int(seconds)
    millis = int(round((seconds - whole) * 1000))
    if millis == 1000:
        whole += 1
        millis = 0
    hours, rem = divmod(whole, 3600)
    minutes, secs = divmod(rem, 60)
    return "{:02d}:{:02d}:{:02d}{}{:03d}".format(hours, minutes, secs, dot, millis)


def srt_to_vtt(srt):
    lines = [
        "WEBVTT",
        "",
        "NOTE",
        "Draft cues converted from captions.srt. Final wording and timing are recorded later.",
        "",
    ]
    for block in srt.strip().split("\n\n"):
        rows = block.splitlines()
        timing = rows[1].replace(",", ".")
        if "line:" not in timing:
            timing += " line:95% align:center size:80%"
        lines.append(rows[0])
        lines.append(timing)
        lines.extend(rows[2:])
        lines.append("")
    return "\n".join(lines)


def script_for(video):
    lines = [
        "# " + video["title"],
        "",
        "Shot list for a silent picture track. Draft narration is below and in narration.txt. Final wording, edit, and caption timing are recorded later.",
        "",
        "Pixi environment: " + video["pixi"],
        "",
        "Input data: " + ", ".join("`{}`".format(item) for item in video["inputs"]) + ".",
        "",
        "Prompt on the terminal shots is `$ `. Commands use the build through `$SEAMS`, `$SEAMS_BUILD`, `$SEAMS_BUILD_SRC`, `$SEAMS_LIBDIR`, `$SEAMS_LANGDIR`, and `$VIDEO_OUT`.",
        "",
        "Tests: " + ", ".join("`{}`".format(name) for name in sorted({test for item in video["scenes"] for test in item["tests"]})) + ".",
        "",
        "Picture: 1920x1080, 30 fps, crf 16, one vhs scene file per shot so the lead-in can be retimed. Each shot starts {:.2f} to {:.2f} s before the draft line. Captions are a `<track>`, in the empty lower band (`line:95%`), and are not burned in. A later mix uses `-c:a aac -aac_coder fast`.".format(min(item["lead"] for item in video["scenes"]), max(item["lead"] for item in video["scenes"])),
        "",
    ]
    for index, item in enumerate(video["scenes"], 1):
        lines.append("## Scene {}. {}".format(index, item["title"]))
        lines.append("")
        lines.append("Lead-in: {:.2f} s. Segment length: {:.2f} s.".format(item["lead"], item["duration"]))
        lines.append("")
        lines.append("Tests: " + ", ".join("`{}`".format(name) for name in item["tests"]) + ".")
        lines.append("")
        lines.append("Commands:")
        lines.append("")
        lines.append("```bash")
        for command in item["commands"]:
            lines.append(command)
        lines.append("```")
        lines.append("")
        if item["expects"]:
            lines.append("Expected output (substrings the run must print):")
            lines.append("")
            lines.append("```text")
            for expect in item["expects"]:
                lines.append(expect)
            lines.append("```")
            lines.append("")
        lines.append("Draft narration:")
        lines.append("")
        lines.append(item["narration"])
        lines.append("")
    return "\n".join(lines)


def write_all():
    catalog = videos()
    for video in catalog:
        brine_temperature_scene(video)
        if video["slug"] == "rust" and "docs_examples" not in video["scenes"][-1]["tests"]:
            video["scenes"][-1]["tests"].append("docs_examples")
    for video in catalog:
        folder = OUT / video["slug"]
        folder.mkdir(parents=True, exist_ok=True)
        paragraphs = [item["narration"] for item in video["scenes"]]
        (folder / "narration.txt").write_text("\n\n".join(paragraphs) + "\n")
        (folder / "script.md").write_text(script_for(video))
        srt = srt_for(video["scenes"])
        (folder / "captions.srt").write_text(srt)
        (folder / "captions.vtt").write_text(srt_to_vtt(srt))
        for index, item in enumerate(video["scenes"], 1):
            (folder / "scene-{:02d}.tape".format(index)).write_text(tape_text(item, index))
    return [video["slug"] for video in catalog]


if __name__ == "__main__":
    print("\n".join(write_all()))
