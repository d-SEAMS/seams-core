#!/usr/bin/env python3
# Copyright (c) d-SEAMS developers
"""Run every docs/videos tape headless and check the spoken numbers."""

import os
import re
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


def clean(text):
    kept = []
    for line in text.splitlines():
        if leaked(line):
            continue
        kept.append(line)
    return "\n".join(kept).strip()


def prepare(root, seams, libdir, langdir, build):
    root = Path(root)
    build = Path(build)
    out = build / "video-runs"
    out.mkdir(parents=True, exist_ok=True)
    shim = build / "video-bin"
    shim.mkdir(parents=True, exist_ok=True)
    for name, target in (("seams", seams), ("python3", sys.executable)):
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
        [str(libdir), str(langdir), env.get("LD_LIBRARY_PATH", "")]
    )
    if not env.get("FENNEL_LUA"):
        for candidate in ("/usr/local/share/fennel.lua", "/usr/share/lua/5.4/fennel.lua"):
            if Path(candidate).is_file():
                env["FENNEL_LUA"] = candidate
                break
    return root, env


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


def check_folder(folder, env, root, record):
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
    combined = []
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
        for command in item["commands"]:
            if leaked(command):
                fail("{} command contains a local path".format(folder.name))
            stdout.append(run_command(command, env, root))
        text = "\n".join(stdout)
        for expect in item["expects"]:
            if expect not in text:
                fail("{} {} missing {!r}\n{}".format(folder.name, item["name"], expect, text[-500:]))
        combined.append(text)
        recorded.append("## {}\n{}".format(item["name"], clean(text)))
    blob = "\n".join(combined)
    missing = [number for number in tokens(narration) if not present(number, blob)]
    if missing:
        fail("{} narration numbers missing from the run: {}".format(folder.name, ", ".join(missing)))
    if record:
        (folder / "recorded.txt").write_text("\n".join(recorded).rstrip() + "\n")
    print(folder.name, "scenes", len(parsed), "ok")


def check_embeds(root):
    """Odysee claims come from embed.txt. Missing claims stay a coming-soon note."""
    from video_page import players_rst, read_claim, slugs

    videos = root / "docs" / "videos"
    if (videos / "lua-fennel").exists():
        fail("lua-fennel should be the lua slug")
    page = (root / "docs" / "orgmode" / "videos" / "index.org").read_text()
    if "huggingface.co/datasets/d-SEAMS/docs-videos" in page:
        fail("videos page still points at the old dataset")
    html = players_rst(videos)
    if "huggingface.co" in html:
        fail("video player html still names huggingface")
    for slug in slugs(videos):
        folder = videos / slug
        if not (folder / "embed.txt").is_file():
            fail("{} is missing embed.txt".format(slug))
        if not (folder / "captions.vtt").is_file():
            fail("{} is missing captions.vtt".format(slug))
        claim = read_claim(folder)
        if "docs/videos/{}/captions.vtt".format(slug) not in html:
            fail("{} transcript link is missing".format(slug))
        if claim is None:
            if "video coming soon" not in html:
                fail("coming soon note is missing")
        else:
            src = "https://odysee.com/$/embed/{}:{}".format(claim[0], claim[1])
            if src not in html:
                fail("{} iframe is missing".format(slug))
    print("video embeds ok")


def main():
    record = "--record" in sys.argv
    args = [item for item in sys.argv[1:] if item != "--record"]
    if len(args) != 5:
        fail("usage: video_tapes.py [--record] ROOT SEAMS LIBDIR LANGDIR BUILD")
    root, env = prepare(*args)
    check_embeds(root)
    folders = sorted(path for path in (root / "docs" / "videos").iterdir() if (path / "script.md").is_file())
    if len(folders) < 14:
        fail("expected 14 videos, found {}".format(len(folders)))
    for folder in folders:
        check_folder(folder, env, root, record)


if __name__ == "__main__":
    main()
