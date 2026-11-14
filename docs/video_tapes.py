#!/usr/bin/env python3
# Copyright (c) d-SEAMS developers
"""Run every docs/videos tape headless and check the spoken numbers."""

import os
import re
import shutil
import subprocess
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))
from tape_parse import parse_tape  # noqa: E402

NUMBER = re.compile(r"(?<![\d.])(\d+(?:\.\d+)?)(?![\d.])")
def leaked(text):
    bad = ("/" + "tmp", "/" + "home", "/" + "agent", "cursor" + "/")
    return any(item in text for item in bad)


def fail(message):
    print(message, file=sys.stderr)
    raise SystemExit(1)


def tokens(text):
    return NUMBER.findall(text)


def present(number, text):
    return re.search(r"(?<![\d.])" + re.escape(number) + r"(?![\d.])", text) is not None


def stable(text):
    """Drop run-to-run timers. The census lines stay exact."""
    text = re.sub(r"overhead_s \S+", "overhead_s TIME", text)
    text = re.sub(r"Loop time of \S+", "Loop time of TIME", text)
    text = re.sub(r"Performance:.*", "Performance: TIME", text)
    text = re.sub(r"[\d.]+% CPU use.*", "CPU TIME", text)
    text = re.sub(r"Total wall time:.*", "Total wall time: TIME", text)
    text = re.sub(r"read_data CPU = \S+ seconds", "read_data CPU = TIME", text)
    text = re.sub(r"(?m)^((?:Pair|Neigh|Comm|Output|Modify|Other)\s+\|).*$", r"\1 TIME", text)
    # LAMMPS notes an unset OMP_NUM_THREADS; the thread count line follows.
    text = re.sub(r"(?m)^OMP_NUM_THREADS environment is not set\..*\n", "", text)
    # gmx grompp draws a quote and a seed, and mdrun's timer row moves.
    text = re.sub(r"Setting the LD random seed to \S+", "Setting the LD random seed to SEED", text)
    text = re.sub(r"GROMACS reminds you:.*", "GROMACS reminds you: QUOTE", text)
    text = re.sub(r"(?m)^[ \t]*Time:.*", "       Time: TIME", text)
    return text


def clean(text):
    kept = []
    for line in text.splitlines():
        if leaked(line):
            continue
        kept.append(line)
    return "\n".join(kept).strip()


def lua_is_54(program):
    done = subprocess.run([program, "-v"], capture_output=True, text=True)
    return "Lua 5.4" in done.stdout + done.stderr


def prepare(root, seams, libdir, langdir, build):
    root = Path(root)
    build = Path(build)
    out = build / "video-runs"
    out.mkdir(parents=True, exist_ok=True)
    shim = build / "video-bin"
    shim.mkdir(parents=True, exist_ok=True)
    links = [("seams", seams), ("python3", sys.executable)]
    # The lua tapes type lua5.4, the Debian name. conda-forge and Homebrew
    # install the same interpreter as lua.
    lua = shutil.which("lua")
    if not shutil.which("lua5.4") and lua and lua_is_54(lua):
        links.append(("lua5.4", lua))
    for name, target in links:
        link = shim / name
        if link.is_symlink() or link.exists():
            link.unlink()
        link.symlink_to(target)
    env = os.environ.copy()
    env["PATH"] = str(shim) + os.pathsep + env.get("PATH", "")
    env["NO_COLOR"] = "1"
    env["SEAMS"] = str(seams)
    env["SEAMS_BUILD"] = str(build)
    env["SEAMS_BUILD_SRC"] = str(libdir)
    env["SEAMS_LIBDIR"] = str(libdir)
    env["SEAMS_LANGDIR"] = str(langdir)
    env["VIDEO_OUT"] = str(out)
    env["DSEAMS_ROOT"] = str(root)
    env["SEAMS_LANG_LIB"] = str(Path(langdir) / "libseams_lang.so")
    env["SEAMS_LIB_DIR"] = str(libdir)
    env["SEAMS_INCLUDE"] = str(root / "bindings" / "cxx")
    env["R_LIBS_USER"] = str(build / "docs-r-lib")
    env["LD_LIBRARY_PATH"] = os.pathsep.join(
        [
            str(libdir),
            str(langdir),
            str(build / "chemfiles-prefix" / "lib"),
            str(build / "subprojects" / "readcon-core"),
            str(build / "subprojects" / "readcon-db"),
            env.get("LD_LIBRARY_PATH", ""),
        ]
    )
    if not env.get("FENNEL_LUA"):
        for candidate in ("/usr/local/share/fennel.lua", "/usr/share/lua/5.4/fennel.lua"):
            if Path(candidate).is_file():
                env["FENNEL_LUA"] = candidate
                break
    return root, env


class Absent(Exception):
    """A step exited 77: its engine or input frame is not on this machine."""


def run_command(command, env, cwd):
    done = subprocess.run(
        command,
        shell=True,
        cwd=str(cwd),
        env=env,
        capture_output=True,
        text=True,
    )
    if done.returncode != 0:
        detail = (done.stderr or done.stdout or "")
        if done.returncode == 77:
            lines = (done.stdout + done.stderr).strip().splitlines()
            raise Absent("{}: {}".format(command, lines[-1] if lines else "exit 77"))
        if done.returncode == 127 or "cannot open shared object file" in detail or "No such file or directory" in detail and "lmp" in command:
            print("skipped: {}".format(detail.strip().splitlines()[-1] if detail.strip() else command), file=sys.stderr)
            raise SystemExit(77)
        fail("exit {} for {}\n{}".format(done.returncode, command, detail[-800:]))
    return done.stdout


def scenes_of(folder):
    tapes = sorted(folder.glob("scene-*.tape"))
    if not tapes:
        fail("no tapes in {}".format(folder.name))
    parsed = []
    for tape in tapes:
        try:
            item = parse_tape(tape.read_text())
        except ValueError as exc:
            fail("{}: {}".format(tape.name, exc))
        item["name"] = tape.name
        parsed.append(item)
    return parsed


def check_folder(folder, env, root, record, absent):
    parsed = scenes_of(folder)
    narration = (folder / "narration.txt").read_text().strip()
    paragraphs = [part.strip() for part in re.split(r"\n\s*\n", narration) if part.strip()]
    if len(paragraphs) != len(parsed):
        fail("{} narration has {} paragraphs for {} scenes".format(folder.name, len(paragraphs), len(parsed)))
    script = (folder / "script.md").read_text()
    if leaked(script) or leaked(narration):
        fail("{} contains a local path".format(folder.name))
    srt = (folder / "captions.srt").read_text()
    vtt = (folder / "captions.vtt").read_text()
    if "WEBVTT" not in vtt or "line:95%" not in vtt:
        fail("{} captions.vtt is missing the track cue".format(folder.name))
    if srt.replace(",", ".") not in vtt.replace(" line:95% align:center size:80%", ""):
        # The VTT is the SRT with comma timestamps rewritten and cue settings added.
        body = []
        for block in srt.strip().split("\n\n"):
            rows = block.splitlines()
            rows[1] = rows[1].replace(",", ".")
            body.append("\n".join(rows))
        for block in body:
            timing = block.splitlines()[1]
            if timing not in vtt:
                fail("{} vtt was not converted from srt ({})".format(folder.name, timing))
    recorded = []
    for item in parsed:
        if not 0.2 <= item["lead"] <= 0.45:
            fail("{} {} lead {:.2f} is outside 0.20-0.45".format(folder.name, item["name"], item["lead"]))
        if not item["commands"]:
            fail("{} {} has no command".format(folder.name, item["name"]))
        for name in item["tests"]:
            if name not in script:
                fail("{} script.md does not name {}".format(folder.name, name))
        stdout = []
        try:
            for command in item["commands"]:
                if leaked(command):
                    fail("{} command contains a local path".format(folder.name))
                stdout.append(run_command(command, env, root))
        except Absent as exc:
            absent.append("{} {}: {}".format(folder.name, item["name"], exc))
            continue
        text = "\n".join(stdout)
        for expect in item["expects"]:
            if expect not in text:
                fail("{} {} missing {!r}\n{}".format(folder.name, item["name"], expect, text[-500:]))
        recorded.append("## {}\n{}".format(item["name"], clean(text)))
    path = folder / "recorded.txt"
    have = path.read_text() if path.is_file() else ""
    if len(recorded) < len(parsed):
        # A skipped scene leaves the transcript partial: each scene that ran
        # still matches its own section, the narration, which quotes the
        # skipped scene too, is not checked, and --record does not rewrite it.
        sections = {}
        for block in re.split(r"(?m)^(?=## )", have):
            if block.startswith("## "):
                sections[block.splitlines()[0]] = block.rstrip() + "\n"
        for block in recorded:
            head = block.splitlines()[0]
            if stable(sections.get(head, "")) != stable(block.rstrip() + "\n"):
                fail("{} recorded.txt does not match the tape output for {}".format(folder.name, head[3:]))
        print(folder.name, "scenes", len(recorded), "of", len(parsed), "ran")
        return
    got = "\n".join(recorded).rstrip() + "\n"
    if record:
        path.write_text(got)
        have = got
    if stable(have) != stable(got):
        fail("{} recorded.txt does not match the tape output".format(folder.name))
    missing = [number for number in tokens(narration) if not present(number, path.read_text())]
    if missing:
        fail("{} narration numbers missing from recorded.txt: {}".format(folder.name, ", ".join(missing)))
    print(folder.name, "scenes", len(parsed), "recorded ok")


def check_embeds(root):
    """A claim in embed.txt is an /$/embed/ iframe. No claim renders nothing."""
    from video_page import CLAIMS, players_rst, read_claim, slugs

    videos = root / "docs" / "videos"
    if (videos / "lua-fennel").exists():
        fail("lua-fennel should be the lua slug")
    page = (root / "docs" / "orgmode" / "videos" / "index.org").read_text()
    if "huggingface.co/datasets/d-SEAMS/docs-videos" in page:
        fail("videos page still points at the old dataset")
    html = players_rst(videos)
    if "huggingface.co" in html or "video coming soon" in html:
        fail("video player html still has a placeholder")
    if "odysee.com" in html.replace("odysee.com/$/embed/", ""):
        fail("video player html has an Odysee URL that is not /$/embed/")
    for slug in slugs(videos):
        folder = videos / slug
        if not (folder / "captions.vtt").is_file():
            fail("{} is missing captions.vtt".format(slug))
        claim = read_claim(folder)
        expected = CLAIMS.get(slug)
        if expected is None:
            if claim is not None:
                fail("{} has a claim that was not supplied".format(slug))
        else:
            name, claim_id = expected.split(":", 1)
            if claim != (name, claim_id):
                fail("{} embed.txt is not {}".format(slug, expected))
        if claim is None:
            if "docs/videos/{}/captions.vtt".format(slug) in html:
                fail("{} has no claim and should render nothing".format(slug))
        else:
            src = "https://odysee.com/$/embed/{}:{}".format(claim[0], claim[1])
            if src not in html:
                fail("{} iframe is missing".format(slug))
            if "docs/videos/{}/captions.vtt".format(slug) not in html:
                fail("{} transcript link is missing".format(slug))
    print("video embeds ok")


def main():
    record = "--record" in sys.argv
    args = [item for item in sys.argv[1:] if item != "--record"]
    if len(args) != 5:
        fail("usage: video_tapes.py [--record] ROOT SEAMS LIBDIR LANGDIR BUILD")
    root, env = prepare(*args)
    from video_page import CUTS, NOTES
    check_embeds(root)
    videos = root / "docs" / "videos"
    names = sorted(path.name for path in videos.iterdir() if path.is_dir() and (path / "script.md").is_file())
    if set(names) != set(CUTS) | set(NOTES):
        fail("video folders are {}, cuts plus notes are {}".format(
            ", ".join(names), ", ".join(sorted(set(CUTS) | set(NOTES)))))
    if len(CUTS) != 19:
        fail("expected 19 cuts, found {}".format(len(CUTS)))
    folders = [videos / name for name in CUTS]
    absent = []
    for folder in folders:
        check_folder(folder, env, root, record, absent)
    if absent:
        for line in absent:
            print("skipped:", line)
        raise SystemExit(77)


if __name__ == "__main__":
    main()
