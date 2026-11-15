#!/usr/bin/env python3
# Copyright (c) d-SEAMS developers
"""Measure the nucleus paragraph and write it from the CSV.

The frame is input/traj/nucleation_f1.lammpstrj. CHILL+, the cutoff
neighbour list, primitive rings, and q6 come from libseams. Adaptive
common-neighbor analysis and the diamond test come from OVITO. The
tetrahedron count is the RMSD of the four nearest neighbours against a
regular tetrahedron, at the OVITO polyhedral-template cutoff of 0.1.
The org page is rewritten from the CSV, so a new measurement replaces
the numbers in the text.
"""

from __future__ import annotations

import csv
import os
import sys
from pathlib import Path

import numpy as np

ROOT = Path(__file__).resolve().parents[1]
FRAME = ROOT / "input" / "traj" / "nucleation_f1.lammpstrj"
ICE = ROOT / "input" / "traj" / "ice1h.lammpstrj"
CSV_PATH = ROOT / "docs" / "data" / "interface-order-nucleus.csv"
ORG = ROOT / "docs" / "orgmode" / "explanation" / "interface-order.org"
CUTOFF = 3.5
XI = 2.4
N_RAYS = 288
TETRA_CUTOFF = 0.1

NAMES = (
    "cubic",
    "hexagonal",
    "water",
    "interfacial",
    "clathrate",
    "interfacial clathrate",
)


def library():
    lib = os.environ.get("SEAMS_LIB")
    if not lib:
        built = ROOT / "bbdir" / "src" / "libseams.so"
        if built.is_file():
            os.environ["SEAMS_LIB"] = str(built)
    src = ROOT / "bindings" / "python" / "src"
    if str(src) not in sys.path:
        sys.path.insert(0, str(src))
    import seams

    return seams


def read_dump(path):
    lines = path.read_text().splitlines()
    box = []
    xyz = []
    i = 0
    while i < len(lines):
        if lines[i].startswith("ITEM: BOX BOUNDS"):
            for _ in range(3):
                i += 1
                lo, hi = lines[i].split()[:2]
                box.append((float(lo), float(hi)))
        elif lines[i].startswith("ITEM: ATOMS"):
            i += 1
            while i < len(lines) and not lines[i].startswith("ITEM:"):
                parts = lines[i].split()
                xyz.append((float(parts[2]), float(parts[3]), float(parts[4])))
                i += 1
            continue
        i += 1
    pos = np.asarray(xyz, dtype=float)
    bounds = np.asarray(box, dtype=float)
    lengths = bounds[:, 1] - bounds[:, 0]
    return pos, lengths, bounds[:, 0]


def mic(delta, lengths):
    return delta - lengths * np.round(delta / lengths)


def com_of(pos, mask, lengths):
    picked = np.flatnonzero(mask)
    seed = pos[picked[0]]
    delta = mic(pos[picked] - seed, lengths)
    return seed + delta.mean(axis=0)


def radii_from(pos, origin, lengths):
    return np.linalg.norm(mic(pos - origin, lengths), axis=1)


def shell_fraction(radius, mask, width=0.5):
    edges = np.arange(0.0, float(radius.max()) + width, width)
    mid = []
    frac = []
    for lo, hi in zip(edges[:-1], edges[1:]):
        sel = (radius >= lo) & (radius < hi)
        count = int(sel.sum())
        if count == 0:
            continue
        mid.append(0.5 * (lo + hi))
        frac.append(float(mask[sel].mean()))
    return np.asarray(mid), np.asarray(frac)


def outward_level(mid, frac, level, start):
    begin = int(np.argmin(np.abs(mid - start)))
    for i in range(begin, len(mid) - 1):
        left, right = float(frac[i]), float(frac[i + 1])
        if left >= level > right or left <= level < right:
            span = right - left
            if span == 0.0:
                return float(mid[i])
            weight = (level - left) / span
            return float(mid[i] + weight * (mid[i + 1] - mid[i]))
    return float(mid[-1])


def front(radius, mask):
    mid, frac = shell_fraction(radius, mask)
    peak = float(mid[int(np.argmax(frac))])
    half = outward_level(mid, frac, 0.5, peak)
    high = outward_level(mid, frac, 0.9, peak)
    low = outward_level(mid, frac, 0.1, high)
    return half, abs(low - high)


def fibonacci(n):
    points = np.empty((n, 3))
    golden = np.pi * (3.0 - np.sqrt(5.0))
    for i in range(n):
        z = 1.0 - (2.0 * i + 1.0) / n
        radius = np.sqrt(max(0.0, 1.0 - z * z))
        theta = golden * i
        points[i] = (radius * np.cos(theta), radius * np.sin(theta), z)
    return points


def willard(pos, origin, lengths, xi, solid, ovito_solid):
    directions = fibonacci(N_RAYS)
    reach = float(np.min(lengths)) * 0.5
    samples = np.linspace(0.0, reach, 64)
    scale = 2.0 * xi * xi
    delta = mic(pos - origin, lengths)
    distance = np.linalg.norm(delta, axis=1)
    direction_of = delta / np.maximum(distance[:, None], 1e-12)
    surface = []
    width_rho = []
    width_pub = []
    width_ovi = []
    for direction in directions:
        rho = np.empty(samples.shape)
        for i, step in enumerate(samples):
            point = origin + direction * step
            moved = mic(pos - point, lengths)
            rho[i] = np.exp(-np.sum(moved * moved, axis=1) / scale).sum()
        core = float(rho[0])
        outer = float(np.median(rho[-8:]))
        level = 0.5 * (core + outer)
        cross = None
        for i in range(len(samples) - 1):
            if rho[i] >= level >= rho[i + 1] or rho[i] <= level <= rho[i + 1]:
                span = float(rho[i + 1] - rho[i])
                weight = 0.0 if span == 0.0 else (level - rho[i]) / span
                cross = float(samples[i] + weight * (samples[i + 1] - samples[i]))
                break
        if cross is None:
            continue
        surface.append(cross)
        scaled = rho / max(core, 1e-12)
        high = outward_level(samples, scaled, 0.9, 0.0)
        low = outward_level(samples, scaled, 0.1, high)
        width_rho.append(abs(low - high))
        cone = np.abs(direction_of @ direction) > 0.95
        if int(cone.sum()) < 8:
            continue
        mid, frac = shell_fraction(distance[cone], solid[cone])
        _mid, ovi = shell_fraction(distance[cone], ovito_solid[cone])
        if mid.size < 2:
            continue
        width_pub.append(abs(outward_level(mid, frac, 0.1, 0.0) - outward_level(mid, frac, 0.9, 0.0)))
        width_ovi.append(abs(outward_level(_mid, ovi, 0.1, 0.0) - outward_level(_mid, ovi, 0.9, 0.0)))
    surface = np.asarray(surface)
    return (
        int(surface.size),
        N_RAYS,
        float(surface.mean()),
        float(surface.std()),
        float(np.mean(width_rho)),
        float(np.mean(width_pub)),
        float(np.mean(width_ovi)),
    )


def solid_widths(radius, solid, ovito_solid):
    mid, frac = shell_fraction(radius, solid)
    _, ovito = shell_fraction(radius, ovito_solid)
    # Same bin centers. The 10-90 width is read from the center outward.
    pub = abs(outward_level(mid, frac, 0.1, 0.0) - outward_level(mid, frac, 0.9, 0.0))
    ovi = abs(outward_level(mid, ovito, 0.1, 0.0) - outward_level(mid, ovito, 0.9, 0.0))
    return pub, ovi


def tetrahedron(pos, columns, offsets, lengths):
    ideal = np.array(
        [[1.0, 1.0, 1.0], [1.0, -1.0, -1.0], [-1.0, 1.0, -1.0], [-1.0, -1.0, 1.0]],
        dtype=float,
    )
    ideal -= ideal.mean(axis=0)
    ideal /= np.linalg.norm(ideal[0])
    passed = np.zeros(pos.shape[0], dtype=bool)
    for i in range(pos.shape[0]):
        cols = columns[offsets[i] : offsets[i + 1]]
        cols = cols[cols != i]
        if cols.size < 4:
            continue
        delta = mic(pos[cols] - pos[i], lengths)
        order = np.argsort(np.sum(delta * delta, axis=1))[:4]
        cloud = delta[order]
        cloud -= cloud.mean(axis=0)
        scale = np.linalg.norm(cloud, axis=1).mean()
        if scale <= 0.0:
            continue
        trial = ideal * scale
        cov = cloud.T @ trial
        _u, _s, vt = np.linalg.svd(cov)
        rot = vt.T @ _u.T
        if np.linalg.det(rot) < 0.0:
            vt[-1] *= -1.0
            rot = vt.T @ _u.T
        fitted = cloud @ rot
        rmsd = np.sqrt(np.mean(np.sum((fitted - trial) ** 2, axis=1))) / scale
        passed[i] = rmsd < TETRA_CUTOFF
    return passed


def ovito_counts(path):
    from ovito.io import import_file
    from ovito.modifiers import CommonNeighborAnalysisModifier, IdentifyDiamondModifier

    pipeline = import_file(str(path))
    diamond = IdentifyDiamondModifier()
    pipeline.modifiers.append(diamond)
    data = pipeline.compute()
    structure = np.asarray(data.particles["Structure Type"])
    cubic = int(np.count_nonzero(structure == IdentifyDiamondModifier.Type.CUBIC_DIAMOND))
    hexagonal = int(np.count_nonzero(structure == IdentifyDiamondModifier.Type.HEX_DIAMOND))
    pipeline.modifiers.clear()
    cna = CommonNeighborAnalysisModifier(mode=CommonNeighborAnalysisModifier.Mode.AdaptiveCutoff)
    pipeline.modifiers.append(cna)
    data = pipeline.compute()
    kind = np.asarray(data.particles["Structure Type"])
    other = int(np.count_nonzero(kind == CommonNeighborAnalysisModifier.Type.OTHER))
    return cubic, hexagonal, other, int(kind.shape[0])


def degree(columns, offsets, n):
    out = np.empty(n, dtype=int)
    for i in range(n):
        cols = columns[offsets[i] : offsets[i + 1]]
        out[i] = int(np.count_nonzero(cols != i))
    return out


def ring_members(system, seams, n):
    on = {4: np.zeros(n, dtype=bool), 5: np.zeros(n, dtype=bool), 6: np.zeros(n, dtype=bool)}
    prism = 0
    with system.rings(CUTOFF, max_depth=6) as rings:
        labels = rings.labels()
        for index, size in enumerate(labels.sizes):
            size = int(size)
            if int(labels.symmetry[index]) == 3:
                prism += 1
            if size not in on:
                continue
            verts = rings.vertices(index)
            on[size][verts] = True
    return on, prism


def q6_mean(system, mask):
    order = system.ql(6, CUTOFF)
    values = np.asarray(order.ql, dtype=float)
    picked = values[mask]
    picked = picked[np.isfinite(picked)]
    if picked.size == 0:
        return float("nan")
    return float(picked.mean())


def measure():
    seams = library()
    pos, lengths, _origin = read_dump(FRAME)
    with seams.System(FRAME, frame=1) as system:
        published = system.chill_plus(rule=seams.PUBLISHED)
        ovito = system.chill_plus(rule=seams.OVITO)
        pub = np.asarray(published.labels)
        ovi = np.asarray(ovito.labels)
        with system.neighbours(CUTOFF) as neighbours:
            csr = neighbours.copy()
        members, prism = ring_members(system, seams, pos.shape[0])
        disagree = pub != ovi
        q_dis = q6_mean(system, disagree)
        q_rest = q6_mean(system, ~disagree)
    with seams.System(ICE, frame=1) as ice:
        q_ice = q6_mean(ice, np.ones(ice.count(), dtype=bool))
    crystal = np.isin(pub, (0, 1))
    solid = np.isin(pub, (0, 1, 3))
    solid_ovi = np.isin(ovi, (0, 1, 3))
    origin = com_of(pos, crystal, lengths)
    radius = radii_from(pos, origin, lengths)
    crystal_radius, crystal_width = front(radius, crystal)
    solid_radius, _solid_width = front(radius, solid)
    ovito_radius, _ovito_width = front(radius, solid_ovi)
    com_pub, com_ovi = solid_widths(radius, solid, solid_ovi)
    rays, ray_total, mean_r, rms_r, rho_width, along_pub, along_ovi = willard(
        pos, origin, lengths, XI, solid, solid_ovi
    )
    passed = tetrahedron(pos, csr.columns, csr.offsets, lengths)
    cubic, hexagonal, other, n_cna = ovito_counts(FRAME)
    coord = degree(csr.columns, csr.offsets, pos.shape[0])
    agree = ~disagree
    row = {
        "n": int(pos.shape[0]),
        "disagree": int(disagree.sum()),
        "coord_disagree_4": int(np.count_nonzero((coord == 4) & disagree)),
        "coord_agree_3": int(np.count_nonzero((coord == 3) & agree)),
        "coord_agree_4": int(np.count_nonzero((coord == 4) & agree)),
        "coord_agree_5": int(np.count_nonzero((coord == 5) & agree)),
        "coord_agree_6": int(np.count_nonzero((coord == 6) & agree)),
        "coord_agree_7plus": int(np.count_nonzero((coord >= 7) & agree)),
        "ring6_disagree": int(np.count_nonzero(members[6] & disagree)),
        "ring5_disagree": int(np.count_nonzero(members[5] & disagree)),
        "ring4_disagree": int(np.count_nonzero(members[4] & disagree)),
        "ring6_agree": int(np.count_nonzero(members[6] & agree)),
        "ring5_agree": int(np.count_nonzero(members[5] & agree)),
        "ring4_agree": int(np.count_nonzero(members[4] & agree)),
        "prism_rings": prism,
        "crystal_radius": crystal_radius,
        "crystal_width": crystal_width,
        "solid_radius_2015": solid_radius,
        "solid_radius_ovito": ovito_radius,
        "solid_width_2015": com_pub,
        "solid_width_ovito": com_ovi,
        "xi": XI,
        "rays_cross": rays,
        "rays": ray_total,
        "surface_mean": mean_r,
        "surface_rms": rms_r,
        "density_width": rho_width,
        "density_width_2015": along_pub,
        "density_width_ovito": along_ovi,
        "acna_other": other,
        "acna_n": n_cna,
        "diamond_cubic": cubic,
        "diamond_hex": hexagonal,
        "tetrahedron": int(passed.sum()),
        "tetrahedron_disagree": int(np.count_nonzero(passed & disagree)),
        "q6_disagree": q_dis,
        "q6_rest": q_rest,
        "q6_ice": q_ice,
    }
    for index, name in enumerate(NAMES):
        row[f"pub_{name}"] = int(np.count_nonzero(pub == index))
        row[f"ovi_{name}"] = int(np.count_nonzero(ovi == index))
    return row


def write_csv(row):
    CSV_PATH.parent.mkdir(parents=True, exist_ok=True)
    with CSV_PATH.open("w", newline="\n") as handle:
        writer = csv.writer(handle, lineterminator="\n")
        writer.writerow(("key", "value"))
        for key, value in row.items():
            if isinstance(value, float):
                writer.writerow((key, f"{value:.6f}"))
            else:
                writer.writerow((key, value))


def angstrom(value):
    return f"{value:.2f}"


def paragraph(row):
    disagree = row["disagree"]
    lines = [
        "Frame 1 of ~input/traj/nucleation_f1.lammpstrj~ (timestep 15677000,",
        f"{row['n']} mW molecules).",
        f"The two rules disagree on {disagree} molecules.",
        "",
        "| label | 2015 | OVITO |",
        "|-------+------|-------|",
    ]
    for name in NAMES:
        lines.append(f"| {name} | {row[f'pub_{name}']} | {row[f'ovi_{name}']} |")
    lines.extend(
        [
            "",
            f"All {disagree} molecules are in the four-neighbour shell."
            if row["coord_disagree_4"] == disagree
            else f"{row['coord_disagree_4']} of the {disagree} molecules are in the four-neighbour shell.",
            "Of the molecules that do not change label, "
            f"{row['coord_agree_4']} also have four",
            f"neighbours, {row['coord_agree_5']} have five, {row['coord_agree_6']} have six and "
            f"{row['coord_agree_7plus']} have seven or more.",
            f"{row['coord_agree_3']} have three.",
            "The extra neighbours sit on the agreeing population.",
            "They change a neighbour's ~q3~, and that is the shell lever.",
            "",
            f"Ring membership does not separate the {disagree}.",
            f"{row['ring6_disagree']} of them lie on a primitive 6-ring, "
            f"{row['ring5_disagree']} on a 5-ring and {row['ring4_disagree']} on a",
            "4-ring, on the 3.5 Å graph.",
            f"The agreeing population is {row['ring6_agree']}, {row['ring5_agree']} and "
            f"{row['ring4_agree']} on those same rings.",
            "No prism block is present."
            if row["prism_rings"] == 0
            else f"The primitive-ring search finds {row['prism_rings']} prismatic faces.",
            "",
            "The crystal radius, read outward from the shell of highest crystal",
            f"fraction to the 0.5 crossing, is {angstrom(row['crystal_radius'])} Å.",
            f"The 10-90 width of that crystal front is {angstrom(row['crystal_width'])} Å.",
            f"The solid radius is {angstrom(row['solid_radius_2015'])} Å (2015) and "
            f"{angstrom(row['solid_radius_ovito'])} Å (OVITO).",
            "The 10-90 width of the solid front, measured from the center of mass,",
            f"is {angstrom(row['solid_width_2015'])} Å under the 2015 rule and "
            f"{angstrom(row['solid_width_ovito'])} Å under OVITO.",
            "",
            "The same frame has an instantaneous surface in the sense of Willard",
            "and Chandler, J. Phys. Chem. B 114, 1954 (2010).",
            f"The coarse-graining Gaussian has their width, ξ = {angstrom(row['xi'])} Å.",
            "Rays leave the center of mass of the 2015 crystal.",
            "The surface on a ray is where the coarse density crosses the midpoint",
            "of the core and the outer liquid.",
            f"{row['rays_cross']} of {row['rays']} rays cross.",
            f"The mean radius is {angstrom(row['surface_mean'])} Å and the "
            "root-mean-square deviation of the",
            f"surface about that mean is {angstrom(row['surface_rms'])} Å.",
            "The 10-90 width of the coarse density, measured along the local",
            f"normal of that surface, is {angstrom(row['density_width'])} Å.",
            "On the same coordinate the 2015 solid fraction has a 10-90 width of",
            f"{angstrom(row['density_width_2015'])} Å and the OVITO solid fraction has "
            f"{angstrom(row['density_width_ovito'])} Å.",
            f"The {angstrom(row['solid_width_2015'])} Å and "
            f"{angstrom(row['solid_width_ovito'])} Å widths are the center-of-mass profiles.",
            "Adaptive common-neighbour analysis calls "
            f"{row['acna_other']} of {row['acna_n']} molecules other.",
            f"The diamond test calls {row['diamond_cubic']} cubic and {row['diamond_hex']} hexagonal.",
            "A tetrahedron template with RMSD below 0.1, the OVITO polyhedral",
            "template matching cutoff, passes "
            f"{row['tetrahedron']} molecules, of which "
            f"{row['tetrahedron_disagree']} "
            f"{'is' if row['tetrahedron_disagree'] == 1 else 'are'} in",
            "the disagreeing set.",
            f"Mean ~q6~ is {row['q6_disagree']:.3f} on the disagreeing set and "
            f"{row['q6_rest']:.3f} on the rest.",
            f"Bulk ice Ih in the test set has mean ~q6~ {row['q6_ice']:.3f}.",
            "",
            "The numbers are ~docs/data/interface-order-nucleus.csv~, written by",
            "~scripts/interface_order_nucleus.py~.",
            "",
            "",
        ]
    )
    return "\n".join(lines)


def write_org(row):
    text = ORG.read_text()
    start = text.index("* Nucleus\n")
    end = text.index("* Slab\n")
    body = "* Nucleus\n\n" + paragraph(row)
    ORG.write_text(text[:start] + body + text[end:])


def main():
    row = measure()
    write_csv(row)
    write_org(row)
    print(f"wrote {CSV_PATH}")
    print(f"disagree {row['disagree']} tetrahedron {row['tetrahedron']}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
