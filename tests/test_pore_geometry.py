#!/usr/bin/env python3
# Copyright (c) d-SEAMS developers
"""Analytic oracles for triclinic density, pore shells, and host-guest contacts."""

from __future__ import annotations

import ctypes
import math
import sys
from pathlib import Path

import pore_api


def fail(message: str) -> None:
    raise SystemExit(message)


def near(got: float, want: float, rel: float = 1e-9, abs_tol: float = 0.0) -> None:
    scale = max(abs(want), abs(got), 1.0)
    if abs(got - want) > max(abs_tol, rel * scale):
        fail(f"{got} != {want}")


def flat(profile, density: float) -> None:
    if not profile["number_density"]:
        fail("empty profile")
    total = 0.0
    for rho, mass, volume in zip(
        profile["number_density"], profile["mass_density"], profile["volume"]
    ):
        near(rho, density)
        near(mass, density)
        if not volume > 0.0:
            fail("non-positive bin volume")
        total += rho * volume
    near(total, density * sum(profile["volume"]), abs_tol=1e-9)


def place(box, sx, sy, sz):
    ox, oy, oz, lx, ly, lz, xy, xz, yz = box
    return (
        ox + lx * sx + xy * sy + xz * sz,
        oy + ly * sy + yz * sz,
        oz + lz * sz,
    )


def uniform(lib, box, axis: int, n: int = 8):
    xyz = []
    types = []
    for i in range(n):
        frac = [0.5, 0.5, 0.5]
        frac[axis] = (i + 0.5) / n
        xyz.extend(place(box, *frac))
        types.append(1)
    handle = lib.system(xyz, box, types)
    try:
        return lib.fractional(handle, n, n, axis)
    finally:
        lib.close(handle)


def main() -> None:
    if len(sys.argv) < 2:
        fail("usage: test_pore_geometry.py libseams")
    lib = pore_api.PoreLib(Path(sys.argv[1]))
    ortho = (0.0, 0.0, 0.0, 10.0, 8.0, 12.0, 0.0, 0.0, 0.0)
    triclinic = (0.0, 0.0, 0.0, 10.0, 8.0, 12.0, 2.0, 1.0, -1.5)
    volume = 10.0 * 8.0 * 12.0
    for box in (ortho, triclinic):
        for axis in (0, 1, 2):
            profile = uniform(lib, box, axis)
            flat(profile, 8.0 / volume)
            near(profile["volume"][0], volume / 8.0, abs_tol=1e-12)
            near(profile["coordinate"][0], 0.5 / 8.0, abs_tol=1e-12)

    slab = (0.0, 0.0, 0.0, 4.0, 5.0, 10.0, 0.0, 0.0, 0.0)
    xyz = []
    for i in range(5):
        xyz.extend((2.0, 2.5, (i + 0.5) * 2.0))
    handle = lib.system(xyz, slab, [1] * 5)
    try:
        profile = lib.fractional(handle, 5, 5, 2)
    finally:
        lib.close(handle)
    flat(profile, 5.0 / 200.0)
    for value in profile["volume"]:
        near(value, 40.0, abs_tol=1e-12)

    cube = (0.0, 0.0, 0.0, 6.0, 6.0, 6.0, 0.0, 0.0, 0.0)
    points = []
    points.extend(place(cube, 0.2, 0.2, 0.5))
    points.extend(place(cube, 0.8, 0.8, 0.4))
    handle = lib.system(points, cube, [1, 1])
    try:
        directed = lib.fractional(handle, 2, 2, 0, direction=(1.0, 1.0, 0.0))
    finally:
        lib.close(handle)
    flat(directed, 2.0 / 216.0)
    near(directed["volume"][0], 0.5 * 216.0)
    near(directed["volume"][1], 0.5 * 216.0)

    tilted = []
    for i in range(8):
        tilted.extend(place(triclinic, 0.4, 0.6, (i + 0.5) / 8.0))
    handle = lib.system(tilted, triclinic, [1] * 8)
    try:
        along = lib.fractional(handle, 8, 8, 2, direction=(0.0, 0.0, 1.0))
        frac = lib.fractional(handle, 8, 8, 2)
    finally:
        lib.close(handle)
    for left, right in zip(along["number_density"], frac["number_density"]):
        near(left, right)
    for left, right in zip(along["volume"], frac["volume"]):
        near(left, right)

    lx, ly, lz = 10.0, 10.0, 20.0
    radius = 2.5
    n_shell = 12
    shell_xyz = []
    for i in range(n_shell):
        theta = 2.0 * math.pi * (i + 0.5) / n_shell
        z = (i + 0.5) / n_shell * lz
        shell_xyz.extend((5.0 + radius * math.cos(theta), 5.0 + radius * math.sin(theta), z))
    shell_box = (0.0, 0.0, 0.0, lx, ly, lz, 0.0, 0.0, 0.0)
    handle = lib.system(shell_xyz, shell_box, [1] * n_shell)
    try:
        profiles = lib.profiles(handle, 6, 4, n_shell, 6.0, 2, (0.5, 0.5, 0.0))
    finally:
        lib.close(handle)
    shell = math.pi * (9.0 - 4.0) * lz
    near(profiles["radial"]["volume"][2], shell)
    near(profiles["radial"]["number_density"][2], n_shell / shell)
    for i, rho in enumerate(profiles["radial"]["number_density"]):
        if i != 2:
            near(rho, 0.0, abs_tol=0.0)
    sector = lx * ly * lz / 4.0
    for value in profiles["azimuthal"]["volume"]:
        near(value, sector, rel=1e-8)
    flat(profiles["azimuthal"], n_shell / (lx * ly * lz))
    flat(profiles["axial"], n_shell / (lx * ly * lz))

    long_box = (0.0, 0.0, 0.0, 5.0, 6.0, 20.0, 0.0, 0.0, 0.0)
    handle = lib.system([1.0, 1.0, 1.0], long_box, [8])
    try:
        frame = lib.frame(handle, framework=8, axis=-1)
    finally:
        lib.close(handle)
    if frame["axis_index"] != 2:
        fail(f"auto axis {frame['axis_index']}")

    walls = (0.0, 0.0, 0.0, 10.0, 10.0, 10.0, 0.0, 0.0, 0.0)
    wall_xyz = [7.0, 5.0, 1.0, 3.0, 5.0, 2.0, 5.0, 7.0, 3.0, 5.0, 3.0, 4.0]
    handle = lib.system(wall_xyz, walls, [4, 4, 4, 4])
    try:
        centre = lib.frame(handle, framework=4, axis=2)
    finally:
        lib.close(handle)
    near(centre["centre_frac"][0], 0.5, abs_tol=1e-8)
    near(centre["centre_frac"][1], 0.5, abs_tol=1e-8)

    guest_box = (0.0, 0.0, 0.0, 10.0, 10.0, 10.0, 0.0, 0.0, 0.0)
    guest_xyz = [1.0, 1.0, 1.0, 2.0, 1.0, 1.0, 1.0, 3.0, 1.0, 4.0, 1.0, 1.0]
    handle = lib.system(guest_xyz, guest_box, [1, 3, 4, 3])
    try:
        made = lib.host_guest(handle, 1, 0, 2.5)
        status = lib.lib.seams_host_guest(
            handle, 0, 0, 1.0, None, None, None, None, 0, None, None, 0, None, None, 0
        )
    finally:
        lib.close(handle)
    if made["nearest_index"] != [1]:
        fail(f"nearest {made['nearest_index']}")
    if made["nearest_type"] != [3]:
        fail(f"nearest type {made['nearest_type']}")
    near(made["distance"][0], 1.0, abs_tol=1e-12)
    if made["framework_types"] != [3, 4]:
        fail(f"types {made['framework_types']}")
    if made["contacts"] != [1, 1]:
        fail(f"contacts {made['contacts']}")
    if status != pore_api.NULL:
        fail(f"null counts status {status}")

    handle = lib.system([0.0, 0.0, 0.0], ortho, [1])
    try:
        bins = lib.status_fractional(handle, 0, 2)
        axis = lib.status_fractional(handle, 4, 3)
        zero = lib.status_fractional(handle, 4, 2, direction=(0.0, 0.0, 0.0))
        mass = lib.status_fractional(handle, 4, 2, masses=(-1.0,))
        axis_out = (ctypes.c_double * 3)()
        centre_out = (ctypes.c_double * 3)()
        frac_out = (ctypes.c_double * 3)()
        index = ctypes.c_int()
        diagonal = lib.lib.seams_pore_frame(
            handle,
            0,
            -1,
            (ctypes.c_double * 3)(1.0, 1.0, 0.0),
            3,
            None,
            0,
            axis_out,
            centre_out,
            frac_out,
            ctypes.byref(index),
        )
    finally:
        lib.close(handle)
    if bins != pore_api.RANGE or axis != pore_api.RANGE:
        fail(f"bin status {bins} axis status {axis}")
    if zero != pore_api.ARG or mass != pore_api.ARG:
        fail(f"direction status {zero} mass status {mass}")
    if diagonal != pore_api.ARG:
        fail(f"diagonal pore axis status {diagonal}")
    print("pore geometry python oracles passed")


if __name__ == "__main__":
    main()
