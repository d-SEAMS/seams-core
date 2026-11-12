# Copyright (c) d-SEAMS developers

"""Growth rates of the sixteen mW slabs, with bootstrap intervals.

Reads data/interface/campaign/tTTTT-sS.lammpstrj. Each file is one
temperature and one velocity seed. interface_order growth supplies the
solid thickness of both CHILL+ rules. The slope is ordinary least
squares on the frames whose 2015 solid fraction is between 0.15 and
0.90 and whose thickness is still defined. The mean over the four
seeds is a percentile bootstrap, 10000 draws, Hyndman-Fan type 7,
seed 270. A paired difference excludes zero when its 2.5 and 97.5
percentiles lie on the same side of zero.

    python3 scripts/growth_rates.py --bin bbdir-eessi/tests/interface_order
"""

from __future__ import annotations

import argparse
import os
import subprocess
import time
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
TEMPERATURES = (240, 250, 260, 270, 272, 273)
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


def bootstrap_mean(samples: list[float], n_draw: int = N_DRAW, seed: int = BOOT_SEED) -> tuple[float, float, float]:
    mean = sum(samples) / len(samples)
    state = seed or 1
    means = []
    n = len(samples)
    for _ in range(n_draw):
        acc = 0.0
        for _i in range(n):
            state, drawn = xorshift(state)
            acc += samples[drawn % n]
        means.append(acc / n)
    means.sort()
    return mean, quantile7(means, 0.025), quantile7(means, 0.975)


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


def parse_growth(text: str) -> list[dict]:
    rows = []
    for line in text.splitlines():
        if not line.startswith("frame ") or line.startswith("frame t"):
            continue
        p = line.split()
        rows.append(
            {
                "t": float(p[2]),
                "dis": int(p[4]),
                "face_pub": float(p[5]),
                "face_ovi": float(p[6]),
                "solid_pub": float(p[7]),
                "solid_ovi": float(p[8]),
                "thick_pub": float(p[9]),
                "thick_ovi": float(p[10]),
                "crystal_pub": float(p[11]),
                "crystal_ovi": float(p[12]),
            }
        )
    return rows


def rate_of(rows: list[dict], key: str) -> tuple[float, int]:
    xs = []
    ys = []
    for row in rows:
        if row["thick_pub"] <= 0.0 or row["thick_ovi"] <= 0.0:
            continue
        if not (0.15 < row["solid_pub"] < 0.90):
            continue
        xs.append(row["t"])
        ys.append(row[key])
    return slope(xs, ys), len(xs)


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


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--bin", type=Path, default=ROOT / "bbdir-eessi" / "tests" / "interface_order")
    parser.add_argument("--traj", type=Path, default=ROOT / "data" / "interface" / "campaign")
    parser.add_argument("--out", type=Path, default=ROOT / "docs" / "orgmode" / "explanation" / "growth-campaign.tsv")
    args = parser.parse_args()
    cores = os.cpu_count() or 1
    print(f"cores {cores}")
    t0 = time.perf_counter()
    lines = [
        "temperature\tseed\tnframes\trate_solid_pub_A_per_ps\trate_solid_ovi_A_per_ps\t"
        "rate_crystal_pub_A_per_ps\trate_crystal_ovi_A_per_ps\tn_fit\tseconds"
    ]
    per_t: dict[int, list[tuple[float, float, float, float]]] = {t: [] for t in TEMPERATURES}
    for temperature in TEMPERATURES:
        for seed in SEEDS:
            dump = args.traj / f"t{temperature}-s{seed}.lammpstrj"
            if not dump.exists():
                raise SystemExit(f"missing {dump}")
            started = time.perf_counter()
            rows = parse_growth(run_growth(args.bin, dump))
            elapsed = time.perf_counter() - started
            rsp, nfit = rate_of(rows, "thick_pub")
            rso, _ = rate_of(rows, "thick_ovi")
            rcp, _ = rate_of(rows, "crystal_pub")
            rco, _ = rate_of(rows, "crystal_ovi")
            per_t[temperature].append((rsp, rso, rcp, rco))
            lines.append(
                f"{temperature}\t{seed}\t{len(rows)}\t{rsp:.8f}\t{rso:.8f}\t{rcp:.8f}\t{rco:.8f}\t{nfit}\t{elapsed:.2f}"
            )
            print(f"t={temperature} seed={seed} nfit={nfit} pub={rsp:.5f} ovi={rso:.5f} {elapsed:.1f}s", flush=True)
    lines.append(f"# cores {cores}")
    lines.append(f"# analysis_seconds {time.perf_counter() - t0:.2f}")
    lines.append("# bootstrap n_draw 10000 seed 270 quantile hyndman-fan-7")
    for temperature in TEMPERATURES:
        block = per_t[temperature]
        pub = [row[0] for row in block]
        ovi = [row[1] for row in block]
        pm, plo, phi = bootstrap_mean(pub)
        om, olo, ohi = bootstrap_mean(ovi)
        diff = bootstrap_mean([a - b for a, b in zip(pub, ovi)])
        excludes = diff[1] > 0.0 or diff[2] < 0.0
        lines.append(
            f"# summary T {temperature} solid_pub_mean {pm:.8f} lo {plo:.8f} hi {phi:.8f} "
            f"solid_ovi_mean {om:.8f} lo {olo:.8f} hi {ohi:.8f} "
            f"diff_mean {diff[0]:.8f} lo {diff[1]:.8f} hi {diff[2]:.8f} "
            f"excludes_zero {int(excludes)}"
        )
        print(
            f"T {temperature} pub {pm:.5f} [{plo:.5f}, {phi:.5f}] "
            f"ovi {om:.5f} [{olo:.5f}, {ohi:.5f}] "
            f"diff {diff[0]:.5f} [{diff[1]:.5f}, {diff[2]:.5f}] sig {excludes}"
        )
    args.out.parent.mkdir(parents=True, exist_ok=True)
    args.out.write_text("\n".join(lines) + "\n")
    print(f"wrote {args.out}")


if __name__ == "__main__":
    main()
