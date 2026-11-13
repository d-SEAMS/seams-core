#!/usr/bin/env python3
# Copyright (c) d-SEAMS developers
"""Render silent 1920x1080 picture tracks from the vhs scene tapes.

vhs can play the same .tape files. This renderer is the headless path:
it draws the prompt and the recorded command output, leaves the lower
band empty for a later caption track, and encodes with ffmpeg at crf 16.
"""

import os
import shutil
import subprocess
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))
from tape_parse import parse_tape  # noqa: E402

WIDTH = 1920
HEIGHT = 1080
BAND = 140
FPS = 30


def font(size):
    from PIL import ImageFont

    for path in (
        "/usr/share/fonts/truetype/dejavu/DejaVuSansMono.ttf",
        "/usr/share/fonts/truetype/dejavu/DejaVuSans.ttf",
    ):
        if os.path.isfile(path):
            return ImageFont.truetype(path, size)
    return ImageFont.load_default()


def wrap(text, limit=92):
    rows = []
    for line in text.splitlines() or [""]:
        while len(line) > limit:
            rows.append(line[:limit])
            line = line[limit:]
        rows.append(line)
    return rows[:28]


def frame(prompt_lines, body_lines):
    from PIL import Image, ImageDraw

    image = Image.new("RGB", (WIDTH, HEIGHT), (17, 20, 24))
    draw = ImageDraw.Draw(image)
    mono = font(28)
    y = 36
    for line in prompt_lines + body_lines:
        if y > HEIGHT - BAND - 40:
            break
        draw.text((48, y), line, fill=(236, 232, 224), font=mono)
        y += 36
    return image


def scene_frames(item, recorded):
    prompts = ["$ " + command for command in item["commands"]]
    return [
        ("lead", frame(["$ "], []), item["lead"]),
        ("command", frame(prompts, []), 1.2),
        ("output", frame(prompts, wrap(recorded)), max(1.0, item["duration"] - item["lead"] - 1.2)),
    ]


def encode(images, dest):
    work = dest.parent / (dest.stem + "-frames")
    if work.exists():
        shutil.rmtree(work)
    work.mkdir(parents=True)
    manifest = work / "list.txt"
    rows = []
    for index, (image, seconds) in enumerate(images):
        path = work / "f{:02d}.png".format(index)
        image.save(path)
        rows.append("file '{}'".format(path.name))
        rows.append("duration {:.3f}".format(seconds))
    rows.append("file '{}'".format(Path(rows[-2].split("'")[1]).name))
    manifest.write_text("\n".join(rows) + "\n")
    subprocess.check_call(
        [
            "ffmpeg", "-y", "-f", "concat", "-safe", "0", "-i", str(manifest),
            "-vf", "fps={}".format(FPS),
            "-c:v", "libx264", "-pix_fmt", "yuv420p", "-crf", "16",
            "-movflags", "+faststart",
            str(dest),
        ],
        cwd=str(work),
        stdout=subprocess.DEVNULL,
        stderr=subprocess.DEVNULL,
    )
    shutil.rmtree(work)


def recorded_blocks(path, count):
    if not path.is_file():
        return [""] * count
    parts = []
    current = []
    for line in path.read_text().splitlines():
        if line.startswith("## scene-"):
            if current:
                parts.append("\n".join(current).strip())
            current = []
            continue
        current.append(line)
    if current:
        parts.append("\n".join(current).strip())
    while len(parts) < count:
        parts.append("")
    return parts[:count]


def contact(frames, dest):
    from PIL import Image, ImageDraw

    columns = min(3, len(frames))
    rows = (len(frames) + columns - 1) // columns
    cell_w, cell_h = 640, 360
    sheet = Image.new("RGB", (columns * cell_w, rows * (cell_h + 36)), (255, 255, 255))
    draw = ImageDraw.Draw(sheet)
    label = font(22)
    for index, image in enumerate(frames):
        thumb = image.resize((cell_w, cell_h))
        col, row = index % columns, index // columns
        x, y = col * cell_w, row * (cell_h + 36)
        draw.text((x + 8, y + 4), "scene {:02d}".format(index + 1), fill=(20, 20, 20), font=label)
        sheet.paste(thumb, (x, y + 36))
    sheet.save(dest)


def render_slug(folder, dest):
    tapes = sorted(folder.glob("scene-*.tape"))
    parsed = [parse_tape(tape.read_text()) for tape in tapes]
    blocks = recorded_blocks(folder / "recorded.txt", len(parsed))
    dest.mkdir(parents=True, exist_ok=True)
    outputs = []
    sheets = []
    for index, (item, text) in enumerate(zip(parsed, blocks), 1):
        pieces = scene_frames(item, text)
        images = [(image, seconds) for _name, image, seconds in pieces]
        scene_path = dest / "scene-{:02d}.mp4".format(index)
        encode(images, scene_path)
        outputs.append(scene_path)
        sheets.append(images[-1][0])
    concat = dest / "list.txt"
    concat.write_text("".join("file '{}'\n".format(path.name) for path in outputs))
    final = dest / (folder.name + ".mp4")
    subprocess.check_call(
        [
            "ffmpeg", "-y", "-f", "concat", "-safe", "0", "-i", str(concat),
            "-c", "copy", str(final),
        ],
        cwd=str(dest),
        stdout=subprocess.DEVNULL,
        stderr=subprocess.DEVNULL,
    )
    concat.unlink()
    contact(sheets, dest / "contact.png")
    for name in ("narration.txt", "script.md", "captions.vtt", "recorded.txt"):
        source = folder / name
        if source.is_file():
            shutil.copy(source, dest / name)
    return final


def main():
    root = Path(__file__).resolve().parents[1]
    dest_root = Path(sys.argv[1]) if len(sys.argv) > 1 else root / "build" / "videos"
    folders = sorted(path for path in (root / "docs" / "videos").iterdir() if (path / "script.md").is_file())
    for folder in folders:
        final = render_slug(folder, dest_root / folder.name)
        print(final)


if __name__ == "__main__":
    main()
