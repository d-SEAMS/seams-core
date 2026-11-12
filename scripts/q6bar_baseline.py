# Copyright (c) d-SEAMS developers

"""Pure-Python Lechner-Dellago q6bar on a 3.5 A shell.

The spherical harmonic is the same closed form as ylm::evalUnit for
l = 6: phi = atan2(dx, dy), Condon-Shortley norms, and the average
includes the central molecule. No threshold is fitted. The cutoff is
3.5 A, the shell used by interface_order q6.

    python3 scripts/q6bar_baseline.py frame.lammpstrj
"""

from __future__ import annotations

import math
import os
import sys
from multiprocessing import get_context
from pathlib import Path

CUTOFF = 3.5
# 0.03125 * sqrt(13/pi) and the other six Condon-Shortley norms for l = 6.
N6 = (
    0.03125 * math.sqrt(13.0 / math.pi),
    0.0625 * math.sqrt(273.0 / (2.0 * math.pi)),
    0.015625 * math.sqrt(1365.0 / math.pi),
    0.03125 * math.sqrt(1365.0 / math.pi),
    0.09375 * math.sqrt(91.0 / (2.0 * math.pi)),
    0.09375 * math.sqrt(1001.0 / math.pi),
    0.015625 * math.sqrt(3003.0 / math.pi),
)


def legendre6(abs_m: int, zpow) -> float:
    if abs_m == 0:
        return 231.0 * zpow[6] - 315.0 * zpow[4] + 105.0 * zpow[2] - 5.0
    if abs_m == 1:
        return 33.0 * zpow[5] - 30.0 * zpow[3] + 5.0 * zpow[1]
    if abs_m == 2:
        return 33.0 * zpow[4] - 18.0 * zpow[2] + 1.0
    if abs_m == 3:
        return 11.0 * zpow[3] - 3.0 * zpow[1]
    if abs_m == 4:
        return 11.0 * zpow[2] - 1.0
    if abs_m == 5:
        return zpow[1]
    return 1.0


def add_ylm6(dx: float, dy: float, dz: float, acc_re: list[float], acc_im: list[float], row: int) -> bool:
    r2 = dx * dx + dy * dy + dz * dz
    if r2 <= 0.0:
        return False
    inv = 1.0 / math.sqrt(r2)
    x = dx * inv
    y = dy * inv
    z = dz * inv
    z1 = z
    z2 = z1 * z
    z3 = z2 * z
    z4 = z3 * z
    z5 = z4 * z
    z6 = z5 * z
    zpow = (1.0, z1, z2, z3, z4, z5, z6)
    # m = 0..6 of (y + i x)^m. m = 0 is stored once; +m and -m differ for m > 0.
    wre = 1.0
    wim = 0.0
    for m in range(7):
        amp = N6[m] * legendre6(m, zpow)
        ineg = 6 - m
        acc_re[row + ineg] += amp * wre
        acc_im[row + ineg] += -amp * wim
        if m != 0:
            ipos = 6 + m
            signed = -amp if (m & 1) else amp
            acc_re[row + ipos] += signed * wre
            acc_im[row + ipos] += signed * wim
        nxt_re = wre * y - wim * x
        nxt_im = wre * x + wim * y
        wre = nxt_re
        wim = nxt_im
    return True


def q6bar_values(xs: list[float], ys: list[float], zs: list[float], lengths: tuple[float, float, float], cutoff: float = CUTOFF) -> list[float]:
    """Lechner-Dellago q6bar of every molecule. Orthorhombic minimum image."""
    n = len(xs)
    if n == 0:
        return []
    lx, ly, lz = lengths
    hx = 0.5 * lx
    hy = 0.5 * ly
    hz = 0.5 * lz
    cut2 = cutoff * cutoff
    nx = max(1, int(lx / cutoff))
    ny = max(1, int(ly / cutoff))
    nz = max(1, int(lz / cutoff))
    ox = min(xs)
    oy = min(ys)
    oz = min(zs)
    heads: list[list[int]] = [[] for _ in range(nx * ny * nz)]
    cell_of = [0] * n
    for i in range(n):
        ix = int((xs[i] - ox) / cutoff)
        iy = int((ys[i] - oy) / cutoff)
        iz = int((zs[i] - oz) / cutoff)
        if ix < 0:
            ix = 0
        elif ix >= nx:
            ix = nx - 1
        if iy < 0:
            iy = 0
        elif iy >= ny:
            iy = ny - 1
        if iz < 0:
            iz = 0
        elif iz >= nz:
            iz = nz - 1
        slot = (ix * ny + iy) * nz + iz
        cell_of[i] = slot
        heads[slot].append(i)
    neigh: list[list[int]] = [[] for _ in range(n)]
    for ix in range(nx):
        for iy in range(ny):
            for iz in range(nz):
                slot = (ix * ny + iy) * nz + iz
                here = heads[slot]
                if not here:
                    continue
                seen = set()
                for dxc in (-1, 0, 1):
                    jx = ix + dxc
                    if jx < 0:
                        jx += nx
                    elif jx >= nx:
                        jx -= nx
                    for dyc in (-1, 0, 1):
                        jy = iy + dyc
                        if jy < 0:
                            jy += ny
                        elif jy >= ny:
                            jy -= ny
                        for dzc in (-1, 0, 1):
                            jz = iz + dzc
                            if jz < 0:
                                jz += nz
                            elif jz >= nz:
                                jz -= nz
                            slot_j = (jx * ny + jy) * nz + jz
                            if slot_j in seen:
                                continue
                            seen.add(slot_j)
                            other = heads[slot_j]
                            if not other:
                                continue
                            for i in here:
                                xi = xs[i]
                                yi = ys[i]
                                zi = zs[i]
                                bucket = neigh[i]
                                for j in other:
                                    if j <= i:
                                        continue
                                    dx = xs[j] - xi
                                    dy = ys[j] - yi
                                    dz = zs[j] - zi
                                    if dx > hx:
                                        dx -= lx
                                    elif dx < -hx:
                                        dx += lx
                                    if dy > hy:
                                        dy -= ly
                                    elif dy < -hy:
                                        dy += ly
                                    if dz > hz:
                                        dz -= lz
                                    elif dz < -hz:
                                        dz += lz
                                    if dx * dx + dy * dy + dz * dz < cut2:
                                        bucket.append(j)
                                        neigh[j].append(i)
    width = 13
    qre = [0.0] * (n * width)
    qim = [0.0] * (n * width)
    nused = [0] * n
    for i in range(n):
        row = i * width
        xi = xs[i]
        yi = ys[i]
        zi = zs[i]
        for j in neigh[i]:
            dx = xs[j] - xi
            dy = ys[j] - yi
            dz = zs[j] - zi
            if dx > hx:
                dx -= lx
            elif dx < -hx:
                dx += lx
            if dy > hy:
                dy -= ly
            elif dy < -hy:
                dy += ly
            if dz > hz:
                dz -= lz
            elif dz < -hz:
                dz += lz
            if add_ylm6(dx, dy, dz, qre, qim, row):
                nused[i] += 1
    for i in range(n):
        count = nused[i]
        if count == 0:
            continue
        inv = 1.0 / count
        row = i * width
        for m in range(width):
            qre[row + m] *= inv
            qim[row + m] *= inv
    prefactor = 4.0 * math.pi / 13.0
    out = [0.0] * n
    for i in range(n):
        row = i * width
        acc_re = [0.0] * width
        acc_im = [0.0] * width
        for m in range(width):
            acc_re[m] = qre[row + m]
            acc_im[m] = qim[row + m]
        ncontrib = 1
        for j in neigh[i]:
            jrow = j * width
            for m in range(width):
                acc_re[m] += qre[jrow + m]
                acc_im[m] += qim[jrow + m]
            ncontrib += 1
        inv = 1.0 / ncontrib
        total = 0.0
        for m in range(width):
            re = acc_re[m] * inv
            im = acc_im[m] * inv
            total += re * re + im * im
        out[i] = math.sqrt(prefactor * total)
    return out


def q6bar_checksum(xs, ys, zs, lengths, cutoff: float = CUTOFF) -> float:
    values = q6bar_values(xs, ys, zs, lengths, cutoff)
    return sum(values)


def read_frames(path: Path, limit: int | None = None) -> list[tuple[list[float], list[float], list[float], tuple[float, float, float]]]:
    frames = []
    xs: list[float] = []
    ys: list[float] = []
    zs: list[float] = []
    box: list[tuple[float, float]] = []
    mode = None
    for line in path.read_text().splitlines():
        if line.startswith("ITEM: TIMESTEP"):
            if xs:
                lengths = tuple(hi - lo for lo, hi in box[:3])
                frames.append((xs, ys, zs, lengths))  # type: ignore[arg-type]
                if limit is not None and len(frames) >= limit:
                    return frames
            xs, ys, zs, box = [], [], [], []
            mode = None
            continue
        if line.startswith("ITEM: BOX"):
            mode = "box"
            continue
        if line.startswith("ITEM: ATOMS"):
            mode = "atoms"
            continue
        if line.startswith("ITEM:"):
            mode = None
            continue
        if mode == "box":
            parts = line.split()
            box.append((float(parts[0]), float(parts[1])))
        elif mode == "atoms":
            parts = line.split()
            if len(parts) < 5:
                continue
            xs.append(float(parts[2]))
            ys.append(float(parts[3]))
            zs.append(float(parts[4]))
    if xs:
        lengths = tuple(hi - lo for lo, hi in box[:3])
        frames.append((xs, ys, zs, lengths))  # type: ignore[arg-type]
    return frames


def data_dir() -> Path:
    return Path(os.environ.get("SI4_DATA", "si4-data"))


_WORK: list = []


def _frame_sum(index: int) -> float:
    xs, ys, zs, lengths = _WORK[index]
    return q6bar_checksum(xs, ys, zs, lengths)


def checksum_frames(frames, workers: int) -> float:
    global _WORK
    _WORK = frames
    if workers <= 1 or len(frames) <= 1:
        return sum(_frame_sum(i) for i in range(len(frames)))
    context = get_context("fork")
    with context.Pool(processes=workers) as pool:
        parts = pool.map(_frame_sum, range(len(frames)))
    return sum(parts)


def main() -> None:
    if len(sys.argv) < 2:
        raise SystemExit("usage: q6bar_baseline.py <dump> [nframes]")
    path = Path(sys.argv[1])
    limit = int(sys.argv[2]) if len(sys.argv) > 2 else None
    frames = read_frames(path, limit)
    workers = int(os.environ.get("Q6_WORKERS", "1"))
    total = checksum_frames(frames, workers)
    print(f"frames {len(frames)} checksum {total:.6f}")


if __name__ == "__main__":
    main()
