#!/usr/bin/env python3
# Copyright (c) d-SEAMS developers
"""The molecule stills use the Ruhi palette, then one still is drawn."""

import importlib.util
import os
import subprocess
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
SCRIPT = ROOT / "docs" / "scripts" / "render_molecules.py"

# The deck colors. They come from scripts/chill_palette.py.
sys.path.insert(0, str(ROOT / "scripts"))
from chill_palette import COLOR


def _hex(color):
    return "#{:02X}{:02X}{:02X}".format(*(int(round(channel * 255)) for channel in color))


EXPECTED = {name: _hex(COLOR[name]) for name in (
    "hexagonal", "cubic", "interfacial", "clathrate", "interClathrate", "water",
)}


def fail(message):
    print(message, file=sys.stderr)
    raise SystemExit(1)


def load():
    # The import puts the host OpenGL folders first on LD_LIBRARY_PATH.
    # The renders below must start from the suite's path, as a direct run does.
    saved = os.environ.get("LD_LIBRARY_PATH")
    spec = importlib.util.spec_from_file_location("render_molecules", SCRIPT)
    module = importlib.util.module_from_spec(spec)
    try:
        spec.loader.exec_module(module)
    finally:
        if saved is None:
            os.environ.pop("LD_LIBRARY_PATH", None)
        else:
            os.environ["LD_LIBRARY_PATH"] = saved
    return module


def check(module):
    for name, color in EXPECTED.items():
        if name not in module.PALETTE:
            fail("palette has no {}".format(name))
        if module.PALETTE[name].upper() != color.upper():
            fail("{} is {}, expected {}".format(name, module.PALETTE[name], color))
        style = module.atom_style(name)
        if style[0] != module.palette_rgb(name):
            fail("{} color does not come from the palette".format(name))
    if module.PALETTE["reCubic"].upper() != EXPECTED["cubic"]:
        fail("reCubic drifted from cubic")
    if module.PALETTE["reHex"].upper() != EXPECTED["hexagonal"]:
        fail("reHex drifted from hexagonal")
    for code, name in (
        (0, "cubic"),
        (1, "hexagonal"),
        (2, "water"),
        (3, "interfacial"),
        (4, "clathrate"),
        (5, "interClathrate"),
    ):
        if module.CHILL_NAME[code] != name:
            fail("label {} is {}".format(code, module.CHILL_NAME[code]))
        if module.atom_style(module.CHILL_NAME[code])[0] != module.palette_rgb(name):
            fail("label {} is not {}".format(name, EXPECTED[name]))
    water = module.atom_style("water")
    crystal = module.atom_style("hexagonal")
    if water[1] >= crystal[1]:
        fail("water is not smaller than hexagonal ice")
    if water[2] < 0.5:
        fail("water is not semi-transparent")
    if crystal[2] != 0.0:
        fail("hexagonal ice is not opaque")
    for name in ("cubic", "interfacial", "clathrate", "interClathrate"):
        if module.atom_style(name)[2] != 0.0:
            fail("{} is not opaque".format(name))
    if "mixed" not in module.PALETTE:
        fail("palette has no mixed")
    if module.PALETTE["mixed"].upper() == EXPECTED["interfacial"]:
        fail("mixed uses the interfacial yellow")
    if module.PALETTE["mixed"].upper() != "#7570B3":
        fail("mixed is {}".format(module.PALETTE["mixed"]))
    if module.delta_e(module.palette_rgb("mixed"), module.palette_rgb("interfacial")) < 20:
        fail("mixed is too close to interfacial")


def check_cache(module):
    """A still is reused only for the same inputs, and round-trips intact."""
    import tempfile

    import numpy as np
    from PIL import Image

    payload = {
        "positions": np.arange(12, dtype=np.float64).reshape(4, 3),
        "samples": np.int32(64),
        "panel": np.array("ice-ih"),
    }
    version = "5.3.0"
    first = module.cycles_digest(payload, version)
    if module.cycles_digest(dict(payload), version) != first:
        fail("the cache digest is not stable for the same inputs")
    moved = dict(payload, positions=payload["positions"] + 1e-9)
    if module.cycles_digest(moved, version) == first:
        fail("a moved atom kept the cached still")
    if module.cycles_digest(dict(payload, samples=np.int32(65)), version) == first:
        fail("a new sample count kept the cached still")
    if module.cycles_digest(payload, "5.4.0") == first:
        fail("another Blender release kept the cached still")
    saved = os.environ.get("MOLECULES_CACHE")
    remote = os.environ.get("MOLECULES_REMOTE")
    os.environ["MOLECULES_REMOTE"] = "0"
    with tempfile.TemporaryDirectory() as folder:
        os.environ["MOLECULES_CACHE"] = folder
        try:
            if module.cycles_cached(first) is not None:
                fail("an empty cache returned a still")
            image = Image.new("RGBA", (8, 6), (10, 20, 30, 255))
            module.cycles_store(first, image, 12.5)
            hit = module.cycles_cached(first)
            if hit is None:
                fail("a stored still was not found")
            if list(hit[0].getdata()) != list(image.getdata()) or hit[1] != 12.5:
                fail("the cached still or its world size changed")
            os.unlink(os.path.join(folder, first + ".world"))
            if module.cycles_cached(first) is not None:
                fail("a still without its world size was reused")
        finally:
            if saved is None:
                os.environ.pop("MOLECULES_CACHE", None)
            else:
                os.environ["MOLECULES_CACHE"] = saved
            if remote is None:
                os.environ.pop("MOLECULES_REMOTE", None)
            else:
                os.environ["MOLECULES_REMOTE"] = remote


def main():
    module = load()
    check(module)
    check_cache(module)
    if module.cycles_python() is None:
        print("no Python here imports bpy; pixi install -e cycles provides one")
        return 77
    color = subprocess.run(
        [sys.executable, str(SCRIPT), "--check-color", "--engine", "cycles"],
        check=False,
    )
    if color.returncode != 0:
        return color.returncode
    args = list(sys.argv[1:])
    if "--engine" not in args:
        args = ["--engine", "cycles", *args]
    done = subprocess.run([sys.executable, str(SCRIPT), *args], check=False)
    if done.returncode != 0:
        return done.returncode
    if "--out" in args:
        out = Path(args[args.index("--out") + 1])
        pngs = sorted(out.glob("*/overlay.png")) if out.is_dir() else []
        if pngs:
            frame = subprocess.run(
                [sys.executable, str(ROOT / "docs" / "scripts" / "check_frame.py"), *map(str, pngs)],
                check=False,
            )
            return frame.returncode
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
