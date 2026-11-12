# Copyright (c) d-SEAMS developers

"""Speeds of the NpT basal and secondary-prism slabs.

The speed of one front is L / (2 tau). L is the liquid length in the
first frame, the box length along the face normal minus the 2015 solid
thickness. tau is the first time the solid fraction of that rule
reaches 0.95. This is the direct-coexistence definition of Espinosa,
Navarro, Sanz, Valeriani and Vega, J. Chem. Phys. 145, 211922 (2016),
equation 2. A slope of the solid thickness, divided by two, is kept
beside it. The mean of the four seeds is a percentile bootstrap,
10000 draws, Hyndman-Fan type 7, generator seed 270.

    python3 scripts/npt_rates.py --bin bbdir-eessi/tests/interface_order
"""

from __future__ import annotations

import argparse
import os
import subprocess
import time
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
FACES = ("basal", "prism2")
AXIS = {"basal": 1, "prism2": 2}
TEMPERATURES = (240, 250, 260, 270)
SEEDS = (1, 2, 3, 4)
DUMP_PS = 10.0
N_DRAW = 10000
BOOT_SEED = 270


def quantile7(sorted_vals: list[float], q: float) -> float:
    if len(sorted_vals) == 1:
        return sorted_vals[0]
    pos = q * (len(sorted_vals) - 1)
    lo = int(pos)
    frac = pos - lo
    if lo + 1 >= len(sorted_vals):
        return sorted_vals[-1]
    return sorted_vals[lo] * (1.0 - frac) + sorted_vals[lo + 1] * frac


def xorshift(state: int) -> tuple[int, int]:
    state &= (1 << 64) - 1
    state ^= (state << 13) & ((1 << 64) - 1)
    state ^= (state >> 7)
    state ^= (state << 17) & ((1 << 64) - 1)
    return state, state


def bootstrap_mean(samples: list[float]) -> tuple[float, float, float]:
    mean = sum(samples) / len(samples)
    state = BOOT_SEED
    means = []
    n = len(samples)
    for _ in range(N_DRAW):
        acc = 0.0
        for _i in range(n):
            state, drawn = xorshift(state)
            acc += samples[drawn % n]
        means.append(acc / n)
    means.sort()
    return mean, quantile7(means, 0.025), quantile7(means, 0.975)


def first_box(path: Path) -> tuple[float, float, float]:
    lengths = []
    with path.open() as handle:
        for line in handle:
            if line.startswith("ITEM: BOX"):
                for _ in range(3):
                    lo, hi = next(handle).split()[:2]
                    lengths.append(float(hi) - float(lo))
                break
    return lengths[0], lengths[1], lengths[2]


def parse_growth(text: str) -> list[dict]:
    rows = []
    for line in text.splitlines():
        if not line.startswith("frame ") or line.startswith("frame t"):
            continue
        p = line.split()
        rows.append(
            {
                "t": float(p[2]),
                "solid_pub": float(p[7]),
                "solid_ovi": float(p[8]),
                "thick_pub": float(p[9]),
                "thick_ovi": float(p[10]),
            }
        )
    return rows


def slope(xs: list[float], ys: list[float]) -> float:
    n = len(xs)
    if n < 2:
        return float("nan")
    mx = sum(xs) / n
    my = sum(ys) / n
    num = sum((x - mx) * (y - my) for x, y in zip(xs, ys))
    den = sum((x - mx) ** 2 for x in xs)
    if den <= 0.0:
        return float("nan")
    return num / den


def front_speed(rows: list[dict], key_solid: str, key_thick: str, liquid: float) -> tuple[float, float, float]:
    tau = float("nan")
    for row in rows:
        if row[key_solid] >= 0.95:
            tau = row["t"]
            break
    direct = liquid / (2.0 * tau) if tau > 0.0 else float("nan")
    xs = []
    ys = []
    for row in rows:
        if row[key_thick] <= 0.0:
            continue
        if not (0.55 <= row[key_solid] <= 0.95):
            continue
        xs.append(row["t"])
        ys.append(row[key_thick])
    fitted = slope(xs, ys) / 2.0
    return direct, fitted, tau


def run_growth(binary: Path, dump: Path) -> str:
    env = os.environ.copy()
    lib = binary.resolve().parents[1] / "src"
    env["LD_LIBRARY_PATH"] = str(lib) + os.pathsep + env.get("LD_LIBRARY_PATH", "")
    env["OMP_NUM_THREADS"] = "1"
    done = subprocess.run(
        [str(binary), "growth", str(dump), str(DUMP_PS)],
        cwd=ROOT,
        env=env,
        check=True,
        capture_output=True,
        text=True,
    )
    return done.stdout


def fmt(value: float) -> str:
    if value != value:
        return "nan"
    return f"{value:.6f}"


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--bin", type=Path, default=ROOT / "bbdir-eessi" / "tests" / "interface_order")
    parser.add_argument("--traj", type=Path, default=ROOT / "data" / "interface" / "npt")
    parser.add_argument("--out", type=Path, default=ROOT / "docs" / "orgmode" / "explanation" / "npt-growth.tsv")
    args = parser.parse_args()
    cores = os.cpu_count() or 1
    print(f"cores {cores}")
    t0 = time.perf_counter()
    lines = [
        "face\ttemperature\tseed\tnframes\tliquid_A\t"
        "u_pub_A_per_ps\tu_ovi_A_per_ps\t"
        "slope_pub_A_per_ps\tslope_ovi_A_per_ps\t"
        "tau_pub_ps\ttau_ovi_ps\tseconds"
    ]
    block: dict[tuple[str, int], list[tuple[float, float, float, float]]] = {}
    for face in FACES:
        for temperature in TEMPERATURES:
            block[(face, temperature)] = []
            for seed in SEEDS:
                dump = args.traj / f"{face}-t{temperature}-s{seed}.lammpstrj"
                if not dump.exists():
                    raise SystemExit(f"missing {dump}")
                started = time.perf_counter()
                rows = parse_growth(run_growth(args.bin, dump))
                box = first_box(dump)
                liquid = box[AXIS[face]] - rows[0]["thick_pub"]
                up, sp, tp = front_speed(rows, "solid_pub", "thick_pub", liquid)
                uo, so, to = front_speed(rows, "solid_ovi", "thick_ovi", liquid)
                elapsed = time.perf_counter() - started
                block[(face, temperature)].append((up, uo, sp, so))
                lines.append(
                    f"{face}\t{temperature}\t{seed}\t{len(rows)}\t{liquid:.4f}\t"
                    f"{fmt(up)}\t{fmt(uo)}\t{fmt(sp)}\t{fmt(so)}\t{fmt(tp)}\t{fmt(to)}\t{elapsed:.2f}"
                )
                print(
                    f"{face} T={temperature} seed={seed} u_pub={fmt(up)} u_ovi={fmt(uo)} "
                    f"tau={fmt(tp)} {elapsed:.1f}s",
                    flush=True,
                )
    lines.append(f"# cores {cores}")
    lines.append(f"# analysis_seconds {time.perf_counter() - t0:.2f}")
    lines.append("# bootstrap n_draw 10000 seed 270 quantile hyndman-fan-7")
    lines.append("# u is L/(2 tau) with tau at solid fraction 0.95, in A/ps")
    for face in FACES:
        for temperature in TEMPERATURES:
            rows = block[(face, temperature)]
            usable = [row for row in rows if row[0] == row[0] and row[1] == row[1]]
            if len(usable) < 2:
                lines.append(f"# summary {face} T {temperature} frozen {len(usable)}")
                print(f"{face} T {temperature} frozen {len(usable)}")
                continue
            pub = [row[0] for row in usable]
            ovi = [row[1] for row in usable]
            pm, plo, phi = bootstrap_mean(pub)
            om, olo, ohi = bootstrap_mean(ovi)
            diff = bootstrap_mean([a - b for a, b in zip(pub, ovi)])
            excludes = diff[1] > 0.0 or diff[2] < 0.0
            lines.append(
                f"# summary {face} T {temperature} n {len(usable)} "
                f"u_pub {pm:.6f} lo {plo:.6f} hi {phi:.6f} "
                f"u_ovi {om:.6f} lo {olo:.6f} hi {ohi:.6f} "
                f"diff {diff[0]:.6f} lo {diff[1]:.6f} hi {diff[2]:.6f} "
                f"excludes_zero {int(excludes)} "
                f"u_pub_A_per_ns {pm * 1000:.2f} u_ovi_A_per_ns {om * 1000:.2f}"
            )
            print(
                f"{face} T {temperature} pub {pm * 1000:.1f} [{plo * 1000:.1f}, {phi * 1000:.1f}] A/ns "
                f"ovi {om * 1000:.1f} excludes {excludes}"
            )
    args.out.parent.mkdir(parents=True, exist_ok=True)
    args.out.write_text("\n".join(lines) + "\n")
    print(f"wrote {args.out}")


if __name__ == "__main__":
    main()
