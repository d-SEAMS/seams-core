# Copyright (c) d-SEAMS developers

"""Front speeds under several readings of the same trajectories.

Espinosa, Navarro, Sanz, Valeriani and Vega, J. Chem. Phys. 145,
211922 (2016), equation 2, page 211922-4:

    u = L_growth / (2 tau)

L_growth is about half the box, the liquid that freezes, and the 2
is the pair of interfaces. tau is the time until the slab has frozen.

This script also reports the slope of the solid count,

    u = (dN/dt) * Lz / (2 N)

and the slope of the solid thickness divided by two. Both equal
equation 2 when the freeze is linear and the liquid is half the box.
Lengths are angstroms and times are picoseconds, so a speed in
angstroms per nanosecond is 1000 times a speed in angstroms per
picosecond. One angstrom per nanosecond is 0.1 m/s.

    python3 scripts/espinosa_rates.py --bin bbdir-eessi/tests/interface_order
"""

from __future__ import annotations

import argparse
import os
import subprocess
import time
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
PUBLISHED = {
    240: 490.0,
    250: 375.0,
    260: 266.0,
}


def quantile7(sorted_vals: list[float], q: float) -> float:
    if len(sorted_vals) == 1:
        return sorted_vals[0]
    pos = q * (len(sorted_vals) - 1)
    lo = int(pos)
    frac = pos - lo
    if lo + 1 >= len(sorted_vals):
        return sorted_vals[-1]
    return sorted_vals[lo] * (1.0 - frac) + sorted_vals[lo + 1] * frac


def xorshift(state: int) -> int:
    state &= (1 << 64) - 1
    state ^= (state << 13) & ((1 << 64) - 1)
    state ^= state >> 7
    state ^= (state << 17) & ((1 << 64) - 1)
    return state


def bootstrap_mean(samples: list[float], n_draw: int = 10000, seed: int = 270):
    mean = sum(samples) / len(samples)
    state = seed
    means = []
    n = len(samples)
    for _ in range(n_draw):
        acc = 0.0
        for _i in range(n):
            state = xorshift(state)
            acc += samples[state % n]
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


def first_box(path: Path) -> tuple[float, float, float]:
    lengths = []
    with path.open() as handle:
        for line in handle:
            if line.startswith("ITEM: BOX"):
                for _ in range(3):
                    lo, hi = next(handle).split()[:2]
                    lengths.append(float(hi) - float(lo))
                return lengths[0], lengths[1], lengths[2]
    raise SystemExit(f"no box in {path}")


def parse_growth(text: str) -> list[dict]:
    rows = []
    for line in text.splitlines():
        if not line.startswith("frame ") or line.startswith("frame t"):
            continue
        p = line.split()
        rows.append(
            {
                "t": float(p[2]),
                "n": int(p[3]),
                "solid": float(p[7]),
                "thick": float(p[9]),
                "crystal": float(p[11]),
            }
        )
    return rows


def run_growth(binary: Path, dump: Path, dt: float) -> str:
    env = os.environ.copy()
    lib = binary.resolve().parents[1] / "src"
    env["LD_LIBRARY_PATH"] = str(lib) + os.pathsep + env.get("LD_LIBRARY_PATH", "")
    env["OMP_NUM_THREADS"] = "1"
    done = subprocess.run(
        [str(binary), "growth", str(dump), str(dt)],
        cwd=ROOT,
        env=env,
        check=True,
        capture_output=True,
        text=True,
    )
    return done.stdout


def first_time(rows: list[dict], level: float) -> float:
    for row in rows:
        if row["solid"] >= level:
            return row["t"]
    return float("nan")


def measure(rows: list[dict], lz: float) -> dict[str, float]:
    n = rows[0]["n"]
    liquid = lz - rows[0]["thick"]
    out = {"liquid": liquid, "lz": lz, "n": float(n)}
    for level, name in ((0.95, "u95"), (0.99, "u99")):
        tau = first_time(rows, level)
        out[name] = liquid / (2.0 * tau) if tau == tau and tau > 0.0 else float("nan")
        out[name + "_box"] = lz / (2.0 * tau) if tau == tau and tau > 0.0 else float("nan")
        out[name + "_box_notwo"] = lz / tau if tau == tau and tau > 0.0 else float("nan")
        out["tau" + name[1:]] = tau
    xs = []
    ns = []
    hs = []
    hc = []
    for row in rows:
        if not (0.55 <= row["solid"] <= 0.90):
            continue
        if row["thick"] <= 0.0:
            continue
        xs.append(row["t"])
        ns.append(row["solid"] * n)
        hs.append(row["thick"])
        if row["crystal"] > 0.0:
            hc.append((row["t"], row["crystal"]))
    dndt = slope(xs, ns)
    dhdt = slope(xs, hs)
    out["u_count"] = dndt * lz / (2.0 * n)
    out["u_thick"] = dhdt / 2.0
    if len(hc) >= 2:
        out["u_crystal"] = slope([p[0] for p in hc], [p[1] for p in hc]) / 2.0
    else:
        out["u_crystal"] = float("nan")
    out["nfit"] = float(len(xs))
    return out


def angstrom_per_ns(angstrom_per_ps: float) -> float:
    return angstrom_per_ps * 1000.0


def self_test() -> None:
    # 60 A of liquid, two fronts, frozen at 200 ps: 0.15 A/ps = 150 A/ns.
    assert abs(60.0 / (2.0 * 200.0) * 1000.0 - 150.0) < 1.0e-9
    # The count slope of a linear freeze of half the molecules.
    # dN/dt = (N/2) / tau, u = dN/dt * Lz / (2 N) = Lz / (4 tau).
    # With liquid = Lz/2 this is liquid / (2 tau).
    n = 4000
    lz = 90.0
    tau = 200.0
    dndt = (n / 2.0) / tau
    u = dndt * lz / (2.0 * n)
    assert abs(u - (lz / 2.0) / (2.0 * tau)) < 1.0e-12
    # 1 A/ns = 0.1 m/s. 490 A/ns is 49 m/s, not 4.9 and not 490.
    assert abs(490.0 * 0.1 - 49.0) < 1.0e-12


def fmt(value: float) -> str:
    if value != value:
        return "nan"
    return f"{angstrom_per_ns(value):.3f}"


def main() -> None:
    self_test()
    parser = argparse.ArgumentParser()
    parser.add_argument("--bin", type=Path, default=ROOT / "bbdir-eessi" / "tests" / "interface_order")
    parser.add_argument("--new", type=Path, default=Path(os.environ.get("SI4_DATA", "si4-data")))
    parser.add_argument("--old", type=Path, default=ROOT / "data" / "interface" / "npt")
    parser.add_argument(
        "--out",
        type=Path,
        default=ROOT / "docs" / "orgmode" / "explanation" / "espinosa-gap.tsv",
    )
    args = parser.parse_args()
    jobs = []
    for temperature in (240, 260):
        for seed in (1, 2, 3, 4):
            new = args.new / f"prism2-t{temperature}-s{seed}.lammpstrj"
            old = args.old / f"prism2-t{temperature}-s{seed}.lammpstrj"
            if new.exists():
                jobs.append(("4000-2ps", temperature, seed, new, 5.0))
            if old.exists():
                jobs.append(("8064-0.2ps", temperature, seed, old, 10.0))
    if not jobs:
        raise SystemExit("no trajectories")
    lines = [
        "protocol\ttemperature\tseed\tliquid_A\tlz_A\tnfit\t"
        "u95_A_per_ns\tu99_A_per_ns\t"
        "u_count_A_per_ns\tu_thick_A_per_ns\tu_crystal_A_per_ns\t"
        "u95_fullbox_over_2tau\tu95_fullbox_over_tau\t"
        "printed_A_per_ns\tprinted_over_10"
    ]
    grouped: dict[tuple[str, int], list[dict[str, float]]] = {}
    t0 = time.perf_counter()
    for protocol, temperature, seed, path, dt in jobs:
        rows = parse_growth(run_growth(args.bin, path, dt))
        _lx, _ly, lz = first_box(path)
        measured = measure(rows, lz)
        grouped.setdefault((protocol, temperature), []).append(measured)
        lines.append(
            "\t".join(
                [
                    protocol,
                    str(temperature),
                    str(seed),
                    f"{measured['liquid']:.3f}",
                    f"{measured['lz']:.3f}",
                    f"{int(measured['nfit'])}",
                    fmt(measured["u95"]),
                    fmt(measured["u99"]),
                    fmt(measured["u_count"]),
                    fmt(measured["u_thick"]),
                    fmt(measured["u_crystal"]),
                    fmt(measured["u95_box"]),
                    fmt(measured["u95_box_notwo"]),
                    f"{PUBLISHED[temperature]:.1f}",
                    f"{PUBLISHED[temperature] / 10.0:.1f}",
                ]
            )
        )
        print(
            f"{protocol} T={temperature} s={seed} "
            f"u95={fmt(measured['u95'])} count={fmt(measured['u_count'])} "
            f"thick={fmt(measured['u_thick'])}",
            flush=True,
        )
    lines.append(f"# cores {os.cpu_count() or 1}")
    lines.append(f"# analysis_seconds {time.perf_counter() - t0:.2f}")
    lines.append("# speeds are A/ns. u95 = L_liquid/(2 tau) at solid fraction 0.95.")
    lines.append("# u_count = (dN/dt) * Lz / (2 N) on solid fraction 0.55 to 0.90.")
    lines.append("# printed values are Espinosa et al. Table I, page 211922-4.")
    for key, block in grouped.items():
        protocol, temperature = key
        for name in ("u95", "u_count", "u_thick", "u_crystal"):
            samples = [row[name] for row in block if row[name] == row[name]]
            if len(samples) < 2:
                continue
            mean, lo, hi = bootstrap_mean(samples)
            lines.append(
                f"# summary {protocol} T {temperature} {name} "
                f"n {len(samples)} mean {angstrom_per_ns(mean):.3f} "
                f"lo {angstrom_per_ns(lo):.3f} hi {angstrom_per_ns(hi):.3f} "
                f"over_printed {angstrom_per_ns(mean) / PUBLISHED[temperature]:.3f} "
                f"over_corrected {angstrom_per_ns(mean) / (PUBLISHED[temperature] / 10.0):.3f}"
            )
            print(
                f"{protocol} T {temperature} {name} "
                f"{angstrom_per_ns(mean):.1f} [{angstrom_per_ns(lo):.1f}, {angstrom_per_ns(hi):.1f}]"
            )
    args.out.parent.mkdir(parents=True, exist_ok=True)
    args.out.write_text("\n".join(lines) + "\n")
    if args.out.stat().st_size > 50000:
        raise SystemExit(f"{args.out} exceeds 50 KB")
    print(f"wrote {args.out} ({args.out.stat().st_size} bytes)")


if __name__ == "__main__":
    main()
