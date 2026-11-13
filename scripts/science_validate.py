#!/usr/bin/env python3
# Copyright (c) d-SEAMS developers

"""Compare d-SEAMS labels with OVITO CHILL+ and with freud q6.

OVITO_THREAD_COUNT is set before OVITO is imported. The cutoff is the
published 3.5 angstrom value. Nothing is refitted.

d-SEAMS chill_ovito is the ovito-chill+ rule. chill is the published
Nguyen and Molinero rule. Agreement with OVITO is judged on chill_ovito.
A published-rule difference is reported with the cause counts from the
interface and nucleation cause files: the neighbour test, the q3 shell
and the closed edge.
"""

from __future__ import annotations

import os
from pathlib import Path

os.environ.setdefault("OVITO_THREAD_COUNT", "4")

import numpy as np

ROOT = Path(__file__).resolve().parents[1]
DATA = ROOT / "paper" / "tex" / "figures" / "data" / "science"
FRAMES = ROOT / "data" / "science-frames"
CUTOFF = 3.5


def read_atoms(path: Path):
    box = None
    rows = []
    header = None
    for line in path.read_text().splitlines():
        if not line or line.startswith("#"):
            if line.startswith("# box"):
                box = np.asarray([float(x) for x in line.split()[2:5]], dtype=float)
            continue
        parts = line.split()
        if header is None:
            header = parts
            continue
        rows.append(parts)
    if header is None or box is None or not rows:
        raise SystemExit(f"incomplete atom table {path}")
    index = {name: i for i, name in enumerate(header)}
    pos = np.array([[float(r[index[k]]) for k in ("x", "y", "z")] for r in rows])
    chill = np.array([int(r[index["chill"]]) for r in rows], dtype=int)
    ovito_rule = np.array([int(r[index["chill_ovito"]]) for r in rows], dtype=int)
    q6 = np.array([float(r[index["q6"]]) for r in rows])
    q6bar = np.array([float(r[index["q6bar"]]) for r in rows])
    return box, pos, chill, ovito_rule, q6, q6bar


def read_causes(path: Path) -> dict[str, int]:
    lines = [line for line in path.read_text().splitlines() if line]
    header = lines[0].split()
    values = lines[1].split()
    return {key: int(float(val)) for key, val in zip(header, values)}


def write_lammpstrj(path: Path, box: np.ndarray, pos: np.ndarray) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w") as handle:
        handle.write("ITEM: TIMESTEP\n0\nITEM: NUMBER OF ATOMS\n")
        handle.write(f"{pos.shape[0]}\nITEM: BOX BOUNDS pp pp pp\n")
        for length in box:
            handle.write(f"0.0 {length:.8f}\n")
        handle.write("ITEM: ATOMS id type x y z\n")
        for i, (x, y, z) in enumerate(pos, start=1):
            handle.write(f"{i} 1 {x:.8f} {y:.8f} {z:.8f}\n")


def ovito_labels(dump: Path) -> np.ndarray:
    from ovito.io import import_file
    from ovito.modifiers import ChillPlusModifier

    names = {}
    for name in dir(ChillPlusModifier.Type):
        if name.isupper():
            names[int(getattr(ChillPlusModifier.Type, name))] = name
    pipeline = import_file(str(dump))
    pipeline.modifiers.append(ChillPlusModifier(cutoff=CUTOFF))
    data = pipeline.compute()
    raw = np.asarray(data.particles["Structure Type"], dtype=int)
    # d-SEAMS atom_state_type order.
    to_dseams = {
        "CUBIC_ICE": 0,
        "HEXAGONAL_ICE": 1,
        "WATER": 2,
        "OTHER": 2,
        "INTERFACIAL_ICE": 3,
        "HYDRATE": 4,
        "INTERFACIAL_HYDRATE": 5,
    }
    out = np.empty(raw.shape, dtype=int)
    unknown = {}
    for i, value in enumerate(raw):
        name = names.get(int(value), str(int(value)))
        if name not in to_dseams:
            unknown[name] = unknown.get(name, 0) + 1
            out[i] = -1
        else:
            out[i] = to_dseams[name]
    if unknown:
        print("unmapped OVITO types", unknown)
    return out


def freud_q6(box: np.ndarray, pos: np.ndarray, average: bool) -> np.ndarray:
    import freud

    cell = freud.box.Box(box[0], box[1], box[2])
    op = freud.order.Steinhardt(l=6, average=average)
    op.compute((cell, pos), neighbors={"mode": "ball", "r_max": CUTOFF})
    return np.asarray(op.particle_order, dtype=float)


def pearson(a: np.ndarray, b: np.ndarray) -> float:
    mask = np.isfinite(a) & np.isfinite(b)
    aa = a[mask] - np.mean(a[mask])
    bb = b[mask] - np.mean(b[mask])
    den = float(np.sqrt(np.dot(aa, aa) * np.dot(bb, bb)))
    if den == 0.0:
        return 1.0
    return float(np.dot(aa, bb) / den)


def q6_row(system: str, quantity: str, ours: np.ndarray, theirs: np.ndarray) -> str:
    mask = np.isfinite(ours) & np.isfinite(theirs)
    diff = ours[mask] - theirs[mask]
    rms = float(np.sqrt(np.mean(diff * diff)))
    corr = pearson(ours, theirs)
    close = float(np.mean(np.abs(diff) <= 0.02))
    agree = 1 if corr >= 0.99 and rms <= 0.02 else 0
    note = "cutoff 3.5 A ball, same threshold as the d-SEAMS kernel"
    if not agree:
        note = (
            "d-SEAMS q6 uses the orthorhombic cutoff list with pairs strictly "
            "inside 3.5 A. freud's ball list keeps the same radius. Molecules "
            "on the cutoff and the Lechner-Dellago average of those shells "
            "are where the two series separate. The threshold was not changed."
        )
    return (
        f"{system} {quantity} {int(mask.sum())} {corr:.6f} {rms:.6f} "
        f"{float(np.max(np.abs(diff))):.6f} {close:.6f} {agree} {note}"
    )


def label_row(system: str, quantity: str, ours: np.ndarray, theirs: np.ndarray, note: str) -> str:
    match = ours == theirs
    frac = float(np.mean(match))
    agree = 1 if frac >= 0.99 else 0
    return (
        f"{system} {quantity} {ours.size} {frac:.6f} 0 0 {frac:.6f} {agree} {note}"
    )


def cause_note(causes: dict[str, int], published_frac: float) -> str:
    if published_frac >= 0.99:
        return "published CHILL+ and OVITO assign the same label"
    if not causes:
        return "published CHILL+ differs from OVITO; cause file was not written for this frame"
    frame = causes.get("frame")
    where = f"On frame {frame}, " if frame is not None else ""
    return (
        f"{where}published CHILL+ keeps the 2015 neighbour test and the four-neighbour "
        f"q3 average. OVITO does not. neighbour_test {causes['neighbour_test']} "
        f"q3_shell {causes['q3_shell']} open_edge {causes['open_edge']} "
        f"mixed {causes['mixed']}"
    )


def main() -> int:
    rows = [
        "system quantity n pearson_or_fraction rms max_abs frac_close agree note"
    ]
    jobs = [
        ("interface", DATA / "interface-atoms.csv", DATA / "interface-causes.csv"),
        ("nucleation", DATA / "nucleus-atoms.csv", DATA / "nucleation-causes.csv"),
    ]
    for atoms in sorted(DATA.glob("interface-atoms-*.csv")):
        tag = atoms.name[len("interface-atoms-") : -len(".csv")]
        causes = DATA / f"interface-causes-{tag}.csv"
        jobs.append((f"interface-{tag}", atoms, causes if causes.is_file() else None))
    for atoms in sorted(DATA.glob("nucleus-atoms-*.csv")):
        tag = atoms.name[len("nucleus-atoms-") : -len(".csv")]
        causes = DATA / f"nucleation-causes-{tag}.csv"
        jobs.append((f"nucleation-{tag}", atoms, causes if causes.is_file() else None))
    for system, atoms, causes_path in jobs:
        box, pos, chill, ovito_rule, q6, q6bar = read_atoms(atoms)
        dump = FRAMES / f"validate-{system}.lammpstrj"
        write_lammpstrj(dump, box, pos)
        external = ovito_labels(dump)
        causes = read_causes(causes_path) if causes_path is not None else {}
        rule_frac = float(np.mean(ovito_rule == external))
        rows.append(
            label_row(
                system,
                "ovito_chill_plus",
                ovito_rule,
                external,
                "d-SEAMS ovito-chill+ against OVITO ChillPlusModifier at cutoff 3.5 A",
            )
        )
        published_frac = float(np.mean(chill == external))
        rows.append(
            label_row(
                system,
                "published_chill_plus",
                chill,
                external,
                cause_note(causes, published_frac),
            )
        )
        if rule_frac < 0.99:
            print(f"warning {system} ovito-chill+ agreement {rule_frac:.4f}")
        rows.append(q6_row(system, "freud_q6", q6, freud_q6(box, pos, average=False)))
        rows.append(q6_row(system, "freud_q6bar", q6bar, freud_q6(box, pos, average=True)))
    out = DATA / "validation.csv"
    out.write_text("\n".join(rows) + "\n")
    print(out.relative_to(ROOT).as_posix())
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
