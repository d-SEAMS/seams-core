"""One colour for each role, shared by the docs and the paper figures.

A colour is an sRGB triple in 0..1. reCubic is cubic and reHex is
hexagonal, so those two names may repeat a triple. Every other name
has its own triple. Water is not sodium, and chloride is not clathrate.
"""

from __future__ import annotations

import math

CUBIC = (0.000, 0.447, 0.698)
HEXAGONAL = (0.902, 0.624, 0.000)
# science_render.py yellow. Brine, the interface, and the cages use it.
INTERFACIAL = (0.941, 0.894, 0.259)
# Liquid oxygen. Not the sodium violet.
WATER = (0.780, 0.780, 0.780)
# One clathrate green for the cages and the brine spheres. Not chloride.
CLATHRATE = (0.133, 0.533, 0.200)
# Pink, not the sodium violet. CIEDE2000 from sodium is about 31.
INTER_CLATHRATE = (0.90, 0.35, 0.62)
SODIUM = (0.37, 0.24, 0.60)
CHLORIDE = (0.12, 0.12, 0.12)
# Ruhi coral, and a deep red far enough from that coral to read as another role.
CORAL = (255 / 255, 101 / 255, 93 / 255)
RED = (0.55, 0.02, 0.06)
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
    ("interclathrate", INTER_CLATHRATE),
    ("Na+", SODIUM),
    ("Cl-", CHLORIDE),
)

# Same role, so the same triple is allowed.
ALIASES = {"reCubic": "cubic", "reHex": "hexagonal"}


def _lab(colour: tuple[float, float, float]) -> tuple[float, float, float]:
    def channel(value: float) -> float:
        return value / 12.92 if value <= 0.04045 else ((value + 0.055) / 1.055) ** 2.4

    linear = tuple(channel(float(value)) for value in colour)
    matrix = (
        (0.4124564, 0.3575761, 0.1804375),
        (0.2126729, 0.7151522, 0.0721750),
        (0.0193339, 0.1191920, 0.9503041),
    )
    xyz = [sum(row[i] * linear[i] for i in range(3)) for row in matrix]
    white = (0.95047, 1.0, 1.08883)

    def pivot(value: float) -> float:
        return value ** (1.0 / 3.0) if value > 0.008856 else 7.787 * value + 16.0 / 116.0

    mapped = [pivot(xyz[i] / white[i]) for i in range(3)]
    return (116.0 * mapped[1] - 16.0, 500.0 * (mapped[0] - mapped[1]), 200.0 * (mapped[1] - mapped[2]))


def ciede2000(left: tuple[float, float, float], right: tuple[float, float, float]) -> float:
    """CIEDE2000 distance between two sRGB triples in 0..1."""
    lightness_1, green_red_1, blue_yellow_1 = _lab(left)
    lightness_2, green_red_2, blue_yellow_2 = _lab(right)
    chroma_1 = math.hypot(green_red_1, blue_yellow_1)
    chroma_2 = math.hypot(green_red_2, blue_yellow_2)
    chroma_bar = (chroma_1 + chroma_2) / 2.0
    weight = 0.5 * (1.0 - math.sqrt(chroma_bar**7 / (chroma_bar**7 + 25.0**7)))
    green_red_1p = green_red_1 * (1.0 + weight)
    green_red_2p = green_red_2 * (1.0 + weight)
    chroma_1p = math.hypot(green_red_1p, blue_yellow_1)
    chroma_2p = math.hypot(green_red_2p, blue_yellow_2)
    hue_1 = math.atan2(blue_yellow_1, green_red_1p) % (2.0 * math.pi)
    hue_2 = math.atan2(blue_yellow_2, green_red_2p) % (2.0 * math.pi)
    delta_l = lightness_2 - lightness_1
    delta_c = chroma_2p - chroma_1p
    delta_h = hue_2 - hue_1
    if chroma_1p * chroma_2p == 0.0:
        delta_h = 0.0
    elif delta_h > math.pi:
        delta_h -= 2.0 * math.pi
    elif delta_h < -math.pi:
        delta_h += 2.0 * math.pi
    delta_H = 2.0 * math.sqrt(chroma_1p * chroma_2p) * math.sin(delta_h / 2.0)
    lightness_bar = (lightness_1 + lightness_2) / 2.0
    chroma_bar_p = (chroma_1p + chroma_2p) / 2.0
    if chroma_1p * chroma_2p == 0.0:
        hue_bar = hue_1 + hue_2
    else:
        hue_sum = hue_1 + hue_2
        if abs(hue_1 - hue_2) <= math.pi:
            hue_bar = hue_sum / 2.0
        elif hue_sum < 2.0 * math.pi:
            hue_bar = (hue_sum + 2.0 * math.pi) / 2.0
        else:
            hue_bar = (hue_sum - 2.0 * math.pi) / 2.0
    turn = (
        1.0
        - 0.17 * math.cos(hue_bar - math.radians(30.0))
        + 0.24 * math.cos(2.0 * hue_bar)
        + 0.32 * math.cos(3.0 * hue_bar + math.radians(6.0))
        - 0.20 * math.cos(4.0 * hue_bar - math.radians(63.0))
    )
    rotation = math.radians(30.0) * math.exp(-((math.degrees(hue_bar) - 275.0) / 25.0) ** 2)
    chroma_weight = 2.0 * math.sqrt(chroma_bar_p**7 / (chroma_bar_p**7 + 25.0**7))
    scale_l = 1.0 + (0.015 * (lightness_bar - 50.0) ** 2) / math.sqrt(20.0 + (lightness_bar - 50.0) ** 2)
    scale_c = 1.0 + 0.045 * chroma_bar_p
    scale_h = 1.0 + 0.015 * chroma_bar_p * turn
    cross = -math.sin(2.0 * rotation) * chroma_weight
    return math.sqrt(
        (delta_l / scale_l) ** 2
        + (delta_c / scale_c) ** 2
        + (delta_H / scale_h) ** 2
        + cross * (delta_c / scale_c) * (delta_H / scale_h)
    )


def closest_roles(colours: dict | None = None) -> tuple[float, str, str]:
    """Smallest CIEDE2000 between two different roles."""
    chosen = COLOUR if colours is None else colours
    names = list(chosen)
    best = (float("inf"), "", "")
    for index, left in enumerate(names):
        left_role = ALIASES.get(left, left)
        for right in names[index + 1 :]:
            if ALIASES.get(right, right) == left_role:
                continue
            distance = ciede2000(chosen[left], chosen[right])
            if distance < best[0]:
                best = (distance, left, right)
    return best


def assert_distinct(colours: dict | None = None, minimum: float = 20.0) -> float:
    """Fail when two roles share a triple or sit closer than ``minimum``."""
    chosen = COLOUR if colours is None else colours
    seen: dict[tuple[float, float, float], str] = {}
    for name, colour in chosen.items():
        key = tuple(round(float(channel), 5) for channel in colour)
        role = ALIASES.get(name, name)
        if key in seen and seen[key] != role:
            raise AssertionError(f"{name} and {seen[key]} share {key}")
        seen.setdefault(key, role)
    distance, left, right = closest_roles(chosen)
    if distance < minimum:
        raise AssertionError(
            f"{left} and {right} are CIEDE2000 {distance:.2f}, under {minimum:.0f}"
        )
    return distance
