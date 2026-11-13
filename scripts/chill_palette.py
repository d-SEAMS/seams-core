"""One colour for each role, shared by the docs and the paper figures.

A colour is an sRGB triple in 0..1. reCubic is cubic and reHex is
hexagonal, so those two names may repeat a triple. Every other name
has its own triple. Water is not sodium, and chloride is not clathrate.
"""

from __future__ import annotations

CUBIC = (0.000, 0.447, 0.698)
HEXAGONAL = (0.902, 0.624, 0.000)
# science_render.py yellow. Brine, the interface, and the cages use it.
INTERFACIAL = (0.941, 0.894, 0.259)
# Liquid oxygen. Not the sodium violet.
WATER = (0.780, 0.780, 0.780)
# One clathrate green for the cages and the brine spheres. Not chloride.
CLATHRATE = (0.133, 0.533, 0.200)
# Wine. Distinct from coral and from red.
INTER_CLATHRATE = (0.533, 0.133, 0.333)
SODIUM = (0.37, 0.24, 0.60)
CHLORIDE = (0.12, 0.12, 0.12)
# Ruhi coral and the oxygen red. Neither is the interfacial-clathrate wine.
CORAL = (255 / 255, 101 / 255, 93 / 255)
RED = (255 / 255, 13 / 255, 13 / 255)
UNCLASSIFIED = (0.45, 0.45, 0.45)

CHILL = {
    0: ("cubic Ic", CUBIC),
    1: ("hexagonal Ih", HEXAGONAL),
    3: ("interfacial", INTERFACIAL),
    2: ("liquid", WATER),
    4: ("clathrate", CLATHRATE),
    5: ("interfacial clathrate", INTER_CLATHRATE),
}

COLOUR = {
    "cubic": CUBIC,
    "hexagonal": HEXAGONAL,
    "interfacial": INTERFACIAL,
    "water": WATER,
    "clathrate": CLATHRATE,
    "interClathrate": INTER_CLATHRATE,
    "unclassified": UNCLASSIFIED,
    "reCubic": CUBIC,
    "reHex": HEXAGONAL,
    "sodium": SODIUM,
    "chloride": CHLORIDE,
    "coral": CORAL,
    "red": RED,
}

LEGEND = (
    ("hexagonal", HEXAGONAL),
    ("cubic", CUBIC),
    ("interfacial", INTERFACIAL),
    ("liquid", WATER),
    ("clathrate", CLATHRATE),
    ("interf. clathrate", INTER_CLATHRATE),
    ("Na+", SODIUM),
    ("Cl-", CHLORIDE),
)

# Same role, so the same triple is allowed.
ALIASES = {"reCubic": "cubic", "reHex": "hexagonal"}


def assert_distinct(colours: dict | None = None) -> None:
    """Fail when two roles share a triple."""
    chosen = COLOUR if colours is None else colours
    seen: dict[tuple[float, float, float], str] = {}
    for name, colour in chosen.items():
        key = tuple(round(float(channel), 5) for channel in colour)
        role = ALIASES.get(name, name)
        if key in seen and seen[key] != role:
            raise AssertionError(f"{name} and {seen[key]} share {key}")
        seen.setdefault(key, role)
    if chosen["water"] == chosen["sodium"]:
        raise AssertionError("water and sodium share a colour")
    if chosen["chloride"] == chosen["clathrate"]:
        raise AssertionError("chloride and clathrate share a colour")
    distinct = {
        tuple(round(float(channel), 5) for channel in chosen[name])
        for name in ("coral", "interClathrate", "red")
    }
    if len(distinct) != 3:
        raise AssertionError("coral, interfacial clathrate, and red are not distinct")
