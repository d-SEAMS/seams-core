#!/usr/bin/env python3
# Copyright (c) d-SEAMS developers
"""The molecule stills use the Ruhi palette, then one still is drawn."""

import importlib.util
import subprocess
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
SCRIPT = ROOT / "docs" / "scripts" / "render_molecules.py"

# The deck colours. They come from scripts/chill_palette.py.
sys.path.insert(0, str(ROOT.parent / "scripts"))
from chill_palette import COLOUR


def _hex(colour):
    return "#{:02X}{:02X}{:02X}".format(*(int(round(channel * 255)) for channel in colour))


EXPECTED = {name: _hex(COLOUR[name]) for name in (
    "hexagonal", "cubic", "interfacial", "clathrate", "interClathrate", "water",
)}


def fail(message):
    print(message, file=sys.stderr)
    raise SystemExit(1)


def load():
    spec = importlib.util.spec_from_file_location("render_molecules", SCRIPT)
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


def check(module):
    for name, colour in EXPECTED.items():
        if name not in module.PALETTE:
            fail("palette has no {}".format(name))
        if module.PALETTE[name].upper() != colour.upper():
            fail("{} is {}, expected {}".format(name, module.PALETTE[name], colour))
        style = module.atom_style(name)
        if style[0] != module.palette_rgb(name):
            fail("{} colour does not come from the palette".format(name))
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


def main():
    check(load())
    colour = subprocess.run(
        [sys.executable, str(SCRIPT), "--check-colour", "--engine", "cycles"],
        check=False,
    )
    if colour.returncode != 0:
        return colour.returncode
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
                [sys.executable, str(ROOT / "scripts" / "check_frame.py"), *map(str, pngs)],
                check=False,
            )
            return frame.returncode
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
