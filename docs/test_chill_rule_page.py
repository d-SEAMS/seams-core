#!/usr/bin/env python3
# Copyright (c) d-SEAMS developers
"""Run the 18-atom CHILL+ fixture in both modes and check the page."""

import os
import subprocess
import sys
from ctypes import CDLL, POINTER, c_char_p, c_double, c_int, c_void_p
from pathlib import Path

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


def fail(message):
    print(message, file=sys.stderr)
    raise SystemExit(1)


def read_fixture(path):
    rows = path.read_text().split()
    count = int(rows[0])
    box = [float(rows[1]), float(rows[2]), float(rows[3])]
    coords = []
    cursor = 4
    for _ in range(count):
        coords.append((float(rows[cursor]), float(rows[cursor + 1]), float(rows[cursor + 2])))
        cursor += 3
    focus = int(rows[cursor])
    return box, coords, focus


def write_dump(path, box, coords):
    lines = [
        "ITEM: TIMESTEP",
        "0",
        "ITEM: NUMBER OF ATOMS",
        str(len(coords)),
        "ITEM: BOX BOUNDS pp pp pp",
        "0.0 {:.8f}".format(box[0]),
        "0.0 {:.8f}".format(box[1]),
        "0.0 {:.8f}".format(box[2]),
        "ITEM: ATOMS id type x y z",
    ]
    for index, (x, y, z) in enumerate(coords, 1):
        lines.append("{} 1 {:.6f} {:.6f} {:.6f}".format(index, x, y, z))
    path.write_text("\n".join(lines) + "\n")


def cli_counts(seams, dump, rule):
    env = os.environ.copy()
    env["NO_COLOR"] = "1"
    if rule == "ovito-chill+":
        command = ["ovito-chill+", str(dump), "--type", "1", "--cutoff", "3.5"]
    else:
        command = ["chill-plus", str(dump), "--type", "1", "--cutoff", "3.5", "--rule", rule]
    done = subprocess.run([seams, *command], check=False, capture_output=True, text=True, env=env)
    if done.returncode != 0:
        fail(done.stderr[-400:])
    found = {}
    text = done.stdout
    parts = text.split()
    if not parts or parts[0] != "nop":
        fail(done.stdout)
    for index in range(2, len(parts) - 1, 2):
        if parts[index + 1].lstrip("-").isdigit():
            found[parts[index]] = int(parts[index + 1])
    return found


def c_labels(lang, box, coords, rule):
    xyz = (c_double * (len(coords) * 3))()
    for index, (x, y, z) in enumerate(coords):
        xyz[3 * index] = x
        xyz[3 * index + 1] = y
        xyz[3 * index + 2] = z
    box9 = (c_double * 9)(0, 0, 0, box[0], box[1], box[2], 0, 0, 0)
    types = (c_int * len(coords))(*([1] * len(coords)))
    lang.seams_frame_from_arrays.restype = c_void_p
    lang.seams_frame_from_arrays.argtypes = [
        POINTER(c_double), c_int, POINTER(c_double), POINTER(c_int), POINTER(c_int),
    ]
    lang.seams_frame_free.argtypes = [c_void_p]
    lang.seams_frame_chill_plus_rule.restype = c_int
    lang.seams_frame_chill_plus_rule.argtypes = [
        c_void_p, c_double, c_int, c_char_p, c_char_p, POINTER(c_int), POINTER(c_int),
    ]
    frame = lang.seams_frame_from_arrays(xyz, len(coords), box9, types, None)
    if not frame:
        fail("frame create failed")
    labels = (c_int * len(coords))()
    counts = (c_int * 9)()
    code = lang.seams_frame_chill_plus_rule(
        frame, 3.5, 1, b"cutoff", rule.encode(), labels, counts
    )
    lang.seams_frame_free(frame)
    if code != 0:
        fail("rule {} failed".format(rule))
    return [int(labels[index]) for index in range(len(coords))]


def main():
    root = Path(sys.argv[1])
    seams = sys.argv[2]
    libdir = Path(sys.argv[3])
    langdir = Path(sys.argv[4])
    build = Path(sys.argv[5])
    page = (root / "docs" / "orgmode" / "explanation" / "chill-plus-rules.org").read_text()
    howto = (root / "docs" / "orgmode" / "howto" / "chill-plus.org").read_text()
    migration = (root / "docs" / "orgmode" / "howto" / "migrate-1x.org").read_text()
    if "chill-plus-rules.org" not in howto:
        fail("CHILL+ reference does not link the rule page")
    if "chill-plus-rules.org" not in migration:
        fail("migration guide does not link the rule page")
    for token in ("250", "235", "15", "4096", "all 4096", "SEAMS_CHILL_PUBLISHED", "SEAMS_CHILL_OVITO", "dseams.PUBLISHED", "dseams.OVITO"):
        if token not in page:
            fail("page missing {}".format(token))
    box, coords, focus = read_fixture(root / "input" / "chill-def" / "liquid-four-nearest.txt")
    out = build / "chill-rule"
    out.mkdir(parents=True, exist_ok=True)
    dump = out / "liquid-four-nearest.lammpstrj"
    write_dump(dump, box, coords)
    published_counts = cli_counts(seams, dump, "chill+")
    ovito_counts = cli_counts(seams, dump, "ovito-chill+")
    os.environ["LD_LIBRARY_PATH"] = os.pathsep.join(
        [str(libdir), str(langdir), os.environ.get("LD_LIBRARY_PATH", "")]
    )
    lang = CDLL(str(Path(langdir) / "libseams_lang.so"))
    published = c_labels(lang, box, coords, "chill+")
    ovito = c_labels(lang, box, coords, "ovito-chill+")
    if published[focus] == ovito[focus]:
        fail("focus atom {} did not change label".format(focus))
    line = "fixture focus {} published {} ovito {}".format(
        focus, NAMES[published[focus]], NAMES[ovito[focus]]
    )
    if line not in page:
        fail("page missing {}".format(line))
    changed = [index for index, (left, right) in enumerate(zip(published, ovito)) if left != right]
    if focus not in changed:
        fail("focus atom is not among the changed labels")
    if published_counts.get(NAMES[published[focus]], 0) < 1:
        fail("CLI published counts omit {}".format(NAMES[published[focus]]))
    if ovito_counts.get(NAMES[ovito[focus]], 0) < 1:
        fail("CLI ovito counts omit {}".format(NAMES[ovito[focus]]))
    print(line)
    print("changed", len(changed), "cli", published_counts, ovito_counts)


if __name__ == "__main__":
    main()
