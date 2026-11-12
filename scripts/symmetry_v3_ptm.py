#!/usr/bin/env python3
# Copyright (c) d-SEAMS developers
# SPDX-License-Identifier: MIT
"""Join PTM labels onto report_thermal output and time PTM on ice Ic."""

from __future__ import annotations

import os
import sys
import time

import numpy as np
from ovito.data import DataCollection, Particles, SimulationCell
from ovito.modifiers import PolyhedralTemplateMatchingModifier
from ovito.pipeline import Pipeline, StaticSource

ICE = {0: "Ih", 1: "Ic", 2: "Isd", -1: "other"}
PTM = {6: "Ic", 7: "Ih"}


def read_frame(path):
    lines = open(path, encoding="utf-8").read().splitlines()
    name = lines[0].split()[1]
    temperature = float(lines[1].split()[1])
    role = lines[2].split()[1]
    n = int(lines[3].split()[1])
    box = np.array([[float(x) for x in lines[5 + i].split()] for i in range(3)])
    xyz = np.zeros((n, 3))
    for i, line in enumerate(lines[9 : 9 + n]):
        parts = line.split()
        xyz[i] = [float(parts[0]), float(parts[1]), float(parts[2])]
    return name, temperature, role, box, xyz


def read_labels(path):
    kernel, chill, geom, rmsd, cosines = [], [], [], [], []
    for line in open(path, encoding="utf-8"):
        if line.startswith("index"):
            continue
        parts = line.split()
        kernel.append(int(parts[1]))
        chill.append(int(parts[2]))
        geom.append(int(parts[3]))
        rmsd.append([float(parts[4]), float(parts[5]), float(parts[6])])
        if len(parts) >= 13:
            cosines.append([float(parts[9]), float(parts[10]), float(parts[11]), float(parts[12])])
        else:
            cosines.append([-2, -2, -2, -2])
    return (
        np.array(kernel),
        np.array(chill),
        np.array(geom),
        np.array(rmsd),
        np.array(cosines),
    )


def shell_call(rmsd, cosines):
    """Second-shell call. A bond-graph residual below 0.45 A names the template.
    Otherwise four cosines that all sit outside [0.65, 0.85] name the eclipsed count.
    """
    usable = [(i, rmsd[i]) for i in range(3) if rmsd[i] >= 0.0]
    if usable and min(v for _, v in usable) < 0.45:
        return min(usable, key=lambda item: item[1])[0]
    known = [c for c in cosines if c > -1.0]
    if len(known) < 4 or any(0.65 <= c <= 0.85 for c in known):
        return None
    n_ecl = sum(c > 0.75 for c in known)
    return {0: 1, 1: 0, 2: 2}.get(n_ecl, -1)


def ptm_labels(xyz, box, cutoff):
    data = DataCollection()
    particles = Particles()
    particles.create_property("Position", data=np.ascontiguousarray(xyz))
    cell = SimulationCell(pbc=(True, True, True))
    mat = np.zeros((3, 4))
    mat[:, :3] = box.T
    cell[...] = mat
    data.objects.append(particles)
    data.objects.append(cell)
    pipe = Pipeline(source=StaticSource(data=data))
    mod = PolyhedralTemplateMatchingModifier(rmsd_cutoff=cutoff)
    for structure in mod.structures:
        structure.enabled = True
    pipe.modifiers.append(mod)
    out = pipe.compute()
    return np.asarray(out.particles["Structure Type"])


def family(code):
    return PTM.get(int(code), "other")


def matrix(kernel, ptm):
    kinds = ["Ic", "Ih", "Isd", "other"]
    print("agreement kernel\\\\ptm", " ".join(f"{k:>8}" for k in kinds))
    for row_name, row_code in (("Ic", 1), ("Ih", 0), ("Isd", 2), ("other", -1)):
        counts = []
        for col in kinds:
            n = 0
            for k, p in zip(kernel, ptm):
                if ICE[int(k)] == row_name and family(p) == col:
                    n += 1
            counts.append(n)
        print(f"{row_name:>8}", " ".join(f"{c:8d}" for c in counts))


def verdict(kernel, chill_family, geom, ptm_code):
    """geom is 0 Ih, 1 Ic, 2 Isd, -1 not a shell, -2 ambiguous.

    A shell whose best template RMSD is at least 0.55 A is not bulk ice.
    The kernel's unlabelled value matches that. A crystal label does not.
    """
    ptm_ice = {6: 1, 7: 0}.get(int(ptm_code), -1)
    if geom == -2:
        return "ambiguous"
    if geom < 0:
        geom = -1
    votes = {"kernel": kernel == geom, "chill": chill_family == geom, "ptm": ptm_ice == geom}
    return ",".join(name for name, ok in votes.items() if ok) or "none"


def chill_family(label):
    if label in (0, 7):
        return 1
    if label in (1, 8):
        return 0
    return -1


def report(frame_path):
    name, temperature, role, box, xyz = read_frame(frame_path)
    kernel, chill, geom, rmsd, cosines = read_labels(frame_path + ".labels")
    print(f"\n== {name} T={temperature} role={role} n={len(xyz)}")
    base = ptm_labels(xyz, box, 0.1)
    matrix(kernel, base)
    who = {}
    contested = 0
    for i in range(len(xyz)):
        cf = chill_family(int(chill[i]))
        pf = {6: 1, 7: 0}.get(int(base[i]), -1)
        if not (kernel[i] == cf == pf):
            contested += 1
            key = verdict(int(kernel[i]), cf, int(geom[i]), int(base[i]))
            who[key] = who.get(key, 0) + 1
    print(f"contested {contested}")
    for key in sorted(who):
        print(f"  shell_says {key} {who[key]}")
    calls = {}
    locked = 0
    for i in range(len(xyz)):
        call = shell_call(rmsd[i], cosines[i])
        if call is None:
            key = "borderline"
        else:
            locked += 1
            cf = chill_family(int(chill[i]))
            pf = {6: 1, 7: 0}.get(int(base[i]), -1)
            votes = []
            if int(kernel[i]) == call:
                votes.append("kernel")
            if cf == call:
                votes.append("chill")
            if pf == call:
                votes.append("ptm")
            key = ",".join(votes) or "none"
        calls[key] = calls.get(key, 0) + 1
    print(f"second_shell locked {locked} of {len(xyz)}")
    for key in sorted(calls):
        print(f"  call {key} {calls[key]}")
    print("rmsd_cutoff changed_from_0.1")
    for cutoff in (0.05, 0.1, 0.15, 0.2, 0.3):
        labels_c = base if cutoff == 0.1 else ptm_labels(xyz, box, cutoff)
        changed = int(np.count_nonzero(labels_c != base))
        cubic = int(np.count_nonzero(labels_c == 6))
        hexa = int(np.count_nonzero(labels_c == 7))
        print(f"  {cutoff:.2f} changed {changed} cubic {cubic} hex {hexa}")


def time_ptm():
    lattice = 4.0 * 2.75 / np.sqrt(3.0)
    basis = np.array(
        [
            [0, 0, 0],
            [0.5, 0.5, 0],
            [0.5, 0, 0.5],
            [0, 0.5, 0.5],
            [0.25, 0.25, 0.25],
            [0.75, 0.75, 0.25],
            [0.75, 0.25, 0.75],
            [0.25, 0.75, 0.75],
        ]
    )

    def cell(rep):
        axes = [np.arange(rep) for _ in range(3)]
        ijk = np.stack(np.meshgrid(*axes, indexing="ij"), axis=-1).reshape(-1, 3)
        pts = (ijk[:, None, :] + basis[None, :, :]).reshape(-1, 3) * lattice
        box = np.eye(3) * (rep * lattice)
        return pts, box

    print("\n== ptm structure id, not a cage name")
    for rep in (3,):
        pts, box = cell(rep)
        ptm_labels(pts, box, 0.1)
        samples = []
        for _ in range(6):
            t0 = time.perf_counter()
            ptm_labels(pts, box, 0.1)
            samples.append((time.perf_counter() - t0) * 1000.0)
        ordered = sorted(samples)
        mid = 0.5 * (ordered[2] + ordered[3])
        print(f"ptm-ice molecules={len(pts)} median_ms={mid} raw {' '.join(f'{v:.6f}' for v in samples)}")


def main():
    os.environ.setdefault("OMP_NUM_THREADS", "1")
    for path in sys.argv[1:]:
        if path.endswith(".frame"):
            report(path)
    time_ptm()


if __name__ == "__main__":
    main()
