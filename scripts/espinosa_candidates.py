# Copyright (c) d-SEAMS developers

"""Candidate readings of the mW secondary-prism speed.

Espinosa, Navarro, Sanz, Valeriani and Vega, J. Chem. Phys. 145,
211922 (2016), equation 2:

    u = L / (2 tau)

L is the liquid length, about half the box, and the 2 is the pair of
interfaces. tau is the time until that liquid has frozen. This script
also reports the other readings of L and of tau, the Lechner-Dellago
q6 count on the 3.5 A shell, the potential-energy clock, the three
faces, and the 5 fs step.

Speeds in the table are angstroms per nanosecond. One angstrom per
nanosecond is 0.1 m/s and 10 cm/s. The mW length sigma is 2.3925 A.

    python3 scripts/espinosa_candidates.py --bin bbdir/tests/interface_order
"""

from __future__ import annotations

import argparse
import os
import subprocess
import time
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
SIGMA_A = 2.3925
PUBLISHED = {240: 490.0, 260: 266.0}
# Shi, Fulford, Salvalaglio and Molteni, J. Chem. Phys. 162, 054714
# (2025), Figure 3. Marker and error-bar caps, read from the figure.
SHI = {
    ("prism1", 240): (98.29, 85.20, 111.39),
    ("prism1", 260): (62.14, 51.55, 72.74),
    ("basal", 240): (88.07, 76.41, 99.73),
    ("basal", 260): (39.54, 27.16, 51.91),
}
LEVELS = (0.70, 0.80, 0.90, 0.95, 0.99)


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


def angstrom_per_ns(angstrom_per_ps: float) -> float:
    return angstrom_per_ps * 1000.0


def run_cmd(binary: Path, args: list[str]) -> str:
    env = os.environ.copy()
    lib = binary.resolve().parents[1] / "src"
    env["LD_LIBRARY_PATH"] = str(lib) + os.pathsep + env.get("LD_LIBRARY_PATH", "")
    env["OMP_NUM_THREADS"] = "1"
    done = subprocess.run(
        [str(binary), *args],
        cwd=ROOT,
        env=env,
        check=True,
        capture_output=True,
        text=True,
    )
    return done.stdout


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


def first_time(rows: list[dict], key: str, level: float) -> float:
    for row in rows:
        if row[key] >= level:
            return row["t"]
    return float("nan")


def window_slope(rows: list[dict], key: str, lo: float, hi: float) -> tuple[float, int]:
    xs = []
    ys = []
    for row in rows:
        if lo <= row["solid"] <= hi and row[key] > 0.0:
            xs.append(row["t"])
            ys.append(row[key])
    return slope(xs, ys), len(xs)


def chill_speeds(rows: list[dict], lz: float) -> dict[str, float]:
    n = float(rows[0]["n"])
    liquid = lz - rows[0]["thick"]
    slab = rows[0]["thick"]
    out = {"liquid": liquid, "slab": slab, "lz": lz, "n": n, "solid0": rows[0]["solid"]}
    for level in LEVELS:
        tau = first_time(rows, "solid", level)
        tag = f"{int(level * 100)}"
        out[f"tau{tag}"] = tau
        out[f"u_liq2_{tag}"] = liquid / (2.0 * tau) if tau == tau and tau > 0.0 else float("nan")
        out[f"u_box2_{tag}"] = lz / (2.0 * tau) if tau == tau and tau > 0.0 else float("nan")
        out[f"u_box1_{tag}"] = lz / tau if tau == tau and tau > 0.0 else float("nan")
        out[f"u_slab2_{tag}"] = slab / (2.0 * tau) if tau == tau and tau > 0.0 else float("nan")
        # One interface travels liquid/2. Dividing by two again double counts.
        out[f"u_one2_{tag}"] = (liquid / 2.0) / (2.0 * tau) if tau == tau and tau > 0.0 else float("nan")
        out[f"u_one1_{tag}"] = (liquid / 2.0) / tau if tau == tau and tau > 0.0 else float("nan")
    peak = max(row["solid"] for row in rows)
    tau_full = first_time(rows, "solid", min(0.995, peak)) if peak >= 0.99 else float("nan")
    out["peak"] = peak
    out["tau_full"] = tau_full
    out["u_full"] = liquid / (2.0 * tau_full) if tau_full == tau_full and tau_full > 0.0 else float("nan")
    dndt, nfit = window_slope(rows, "solid", 0.55, 0.90)
    # solid is a fraction, so d(fraction)/dt * n is dN/dt.
    dhdt, _ = window_slope(rows, "thick", 0.55, 0.90)
    dcdt, _ = window_slope(rows, "crystal", 0.55, 0.90)
    out["nfit"] = float(nfit)
    out["u_count"] = dndt * n * lz / (2.0 * n) if dndt == dndt else float("nan")
    out["u_thick"] = dhdt / 2.0 if dhdt == dhdt else float("nan")
    out["u_crystal"] = dcdt / 2.0 if dcdt == dcdt else float("nan")
    return out


def parse_q6(text: str) -> tuple[list[float], list[dict], list[int], list[int]]:
    thresholds: list[float] = []
    rows = []
    ice: list[int] = []
    liq: list[int] = []
    for line in text.splitlines():
        if line.startswith("thresholds "):
            thresholds = [float(x) for x in line.split()[1:]]
        elif line.startswith("hist_ice "):
            ice = [int(x) for x in line.split()[1:]]
        elif line.startswith("hist_liq "):
            liq = [int(x) for x in line.split()[1:]]
        elif line.startswith("frame "):
            p = line.split()
            bar_at = p.index("bar")
            local_at = p.index("local")
            rows.append(
                {
                    "t": float(p[2]),
                    "n": int(p[3]),
                    "bar": [int(x) for x in p[bar_at + 1 : local_at]],
                    "local": [int(x) for x in p[local_at + 1 :]],
                }
            )
    return thresholds, rows, ice, liq


def mislabel_rates(thresholds: list[float], ice: list[int], liq: list[int], t: float) -> tuple[float, float]:
    n_ice = sum(ice)
    n_liq = sum(liq)
    # Bins are 0.02 wide. Bin k holds q in [0.02k, 0.02k+0.02).
    k = int(t / 0.02)
    ice_below = sum(ice[:k]) / n_ice
    liq_above = sum(liq[k:]) / n_liq
    return ice_below, liq_above


def mislabel_threshold(thresholds: list[float], ice: list[int], liq: list[int]) -> tuple[float, float]:
    """Threshold where the two mislabel rates are equal.

    Espinosa et al., page 211922-3: the cut is the point at which the
    fraction of bulk ice labeled liquid equals the fraction of bulk
    liquid labeled ice, and that fraction is about half a percent.
    When a gap separates the two histograms, that point is the low edge
    of the gap, not a cut inside it.
    """
    n_ice = sum(ice)
    n_liq = sum(liq)
    if n_ice == 0 or n_liq == 0:
        return float("nan"), float("nan")
    chosen = None
    chosen_gap = float("nan")
    for t in thresholds:
        ice_below, liq_above = mislabel_rates(thresholds, ice, liq, t)
        if ice_below <= 0.005 and liq_above <= 0.005:
            chosen = t
            chosen_gap = abs(ice_below - liq_above)
            break
    if chosen is None:
        best_gap = 1.0e9
        chosen = thresholds[0]
        for t in thresholds:
            ice_below, liq_above = mislabel_rates(thresholds, ice, liq, t)
            gap = abs(ice_below - liq_above)
            if gap < best_gap:
                best_gap = gap
                chosen = t
                chosen_gap = gap
    return chosen, chosen_gap


def fraction_series(rows: list[dict], thresholds: list[float], level: float, kind: str) -> list[dict]:
    if not thresholds:
        return []
    # nearest grid point
    k = min(range(len(thresholds)), key=lambda i: abs(thresholds[i] - level))
    out = []
    for row in rows:
        counts = row[kind]
        frac = counts[k] / row["n"] if row["n"] else float("nan")
        out.append({"t": row["t"], "n": row["n"], "solid": frac, "thick": 0.0, "crystal": 0.0})
    return out


def q6_speed(rows: list[dict], lz: float) -> float:
    if len(rows) < 2:
        return float("nan")
    n = float(rows[0]["n"])
    xs = []
    ys = []
    for row in rows:
        if 0.55 <= row["solid"] <= 0.90:
            xs.append(row["t"])
            ys.append(row["solid"] * n)
    dndt = slope(xs, ys)
    if dndt != dndt:
        return float("nan")
    return dndt * lz / (2.0 * n)


def parse_energy(path: Path, step_fs: float) -> list[tuple[float, float, float]]:
    """Production (time_ps, pe, long-edge) after reset_timestep 0."""
    blocks: list[list[tuple[float, float, float]]] = []
    current: list[tuple[float, float, float]] = []
    for line in path.read_text(errors="replace").splitlines():
        parts = line.split()
        if len(parts) < 7:
            if current and line.startswith("Loop time"):
                blocks.append(current)
                current = []
            continue
        try:
            step = int(parts[0])
            pe = float(parts[2])
            edges = (float(parts[4]), float(parts[5]), float(parts[6]))
        except ValueError:
            continue
        if step == 0 and current:
            blocks.append(current)
            current = []
        current.append((float(step), pe, max(edges)))
    if current:
        blocks.append(current)
    if not blocks:
        return []
    prod = max(blocks, key=len)
    return [(step * step_fs * 1.0e-3, pe, edge) for step, pe, edge in prod]


def energy_speed(series: list[tuple[float, float, float]]) -> dict[str, float]:
    out = {"u_energy": float("nan"), "tau_energy": float("nan"), "drop": float("nan"), "lz0": float("nan")}
    if len(series) < 5:
        return out
    p0 = sum(p for _, p, _ in series[:3]) / 3.0
    p1 = sum(p for _, p, _ in series[-5:]) / 5.0
    drop = p0 - p1
    out["drop"] = drop
    out["lz0"] = series[0][2]
    if drop <= 1.0:
        return out
    target = p0 - 0.95 * drop
    tau = float("nan")
    for t, pe, _edge in series:
        if pe <= target:
            tau = t
            break
    out["tau_energy"] = tau
    # Half the long edge is the liquid that freezes, and two fronts share it.
    if tau == tau and tau > 0.0:
        out["u_energy"] = (series[0][2] / 2.0) / (2.0 * tau)
    return out


def self_test() -> None:
    # Half the box is liquid: L_box/tau is four times L_liquid/(2 tau).
    liquid = 45.0
    lz = 90.0
    tau = 200.0
    correct = liquid / (2.0 * tau)
    four = lz / tau
    assert abs(four / correct - 4.0) < 1.0e-12
    # 1 A/ns = 0.1 m/s = 10 cm/s. sigma = 2.3925 A.
    assert abs(116.0 * 0.1 - 11.6) < 1.0e-12
    assert abs(116.0 * 10.0 - 1160.0) < 1.0e-12
    assert abs(116.0 / SIGMA_A - 48.485) < 1.0e-2
    # Table/10 times sigma is not a length conversion of a correct angstrom speed.
    assert abs((490.0 / 10.0) * SIGMA_A - 117.23) < 1.0e-2


def fmt(value: float) -> str:
    if value != value:
        return "nan"
    return f"{angstrom_per_ns(value):.3f}"


def discover(data: Path) -> list[tuple[str, int, int, float, Path, Path]]:
    found = []
    for path in sorted(data.glob("*.lammpstrj")):
        name = path.name[: -len(".lammpstrj")]
        log = path.with_suffix(".log")
        parts = name.split("-")
        face = parts[0]
        step = 10.0
        temperature = None
        seed = None
        for part in parts[1:]:
            if part.startswith("dt"):
                step = float(part[2:])
            elif part.startswith("t"):
                temperature = int(part[1:])
            elif part.startswith("s"):
                seed = int(part[1:])
        if temperature is None or seed is None:
            continue
        found.append((face, temperature, seed, step, path, log))
    return found


def main() -> None:
    self_test()
    parser = argparse.ArgumentParser()
    parser.add_argument("--bin", type=Path, default=ROOT / "bbdir" / "tests" / "interface_order")
    parser.add_argument("--data", type=Path, default=Path(os.environ.get("SI4_DATA", "si4-data")))
    parser.add_argument(
        "--out",
        type=Path,
        default=ROOT / "docs" / "orgmode" / "explanation" / "espinosa-gap.tsv",
    )
    args = parser.parse_args()
    jobs = discover(args.data)
    if not jobs:
        raise SystemExit(f"no trajectories in {args.data}")
    header = [
        "face",
        "temperature",
        "seed",
        "step_fs",
        "liquid_A",
        "slab_A",
        "lz_A",
        "solid0",
        "peak",
        "nfit",
        "u_eq2_095",
        "u_eq2_099",
        "u_eq2_full",
        "u_count",
        "u_thick",
        "u_crystal",
        "u_box_over_2tau",
        "u_box_over_tau",
        "u_slab_over_2tau",
        "u_one_over_2tau",
        "u_tau70",
        "u_tau80",
        "u_tau90",
        "u_q6_mislabel",
        "q6_threshold",
        "q6_gap",
        "u_q6_0p358",
        "u_q6_local_0p5",
        "u_energy",
        "energy_drop_kcal",
        "u_m_per_s",
        "u_cm_per_s",
        "u_sigma_per_ns",
        "printed",
        "printed_over_10",
        "printed_over_4",
        "shi_same_face",
    ]
    lines = ["\t".join(header)]
    grouped: dict[tuple[str, int, float], list[dict[str, float]]] = {}
    t0 = time.perf_counter()
    for face, temperature, seed, step, path, log in jobs:
        rows = parse_growth(run_cmd(args.bin, ["growth", str(path), "5.0"]))
        if not rows:
            print(f"empty {path.name}", flush=True)
            continue
        lz = rows and None
        # box length from the q6 line, which follows the long axis
        q6_text = run_cmd(args.bin, ["q6", str(path), "5.0"])
        thresholds, q6_rows, ice, liq = parse_q6(q6_text)
        # growth thickness is along the long axis; lz is that length at frame 1
        # read it from the first q6 frame via a side channel: re-parse
        lz = float("nan")
        for line in q6_text.splitlines():
            if line.startswith("frame 1 "):
                lz = float(line.split()[5])
                break
        measured = chill_speeds(rows, lz)
        q_t, q_gap = mislabel_threshold(thresholds, ice, liq)
        q_rows = fraction_series(q6_rows, thresholds, q_t, "bar")
        q_fixed = fraction_series(q6_rows, thresholds, 0.358, "bar")
        # nearest grid point to 0.50 for the local, un-averaged q6
        q_local = fraction_series(q6_rows, thresholds, 0.50, "local")
        measured["u_q6"] = q6_speed(q_rows, lz)
        measured["u_q6_358"] = q6_speed(q_fixed, lz)
        measured["u_q6_local"] = q6_speed(q_local, lz)
        measured["q6_t"] = q_t
        measured["q6_gap"] = q_gap
        if log.exists():
            series = parse_energy(log, step)
            measured.update(energy_speed(series))
        else:
            measured.update({"u_energy": float("nan"), "drop": float("nan")})
        count_ns = angstrom_per_ns(measured["u_count"])
        shi = SHI.get((face, temperature))
        shi_s = f"{shi[0]:.2f}" if shi else "nan"
        printed = PUBLISHED.get(temperature, float("nan"))
        fields = [
            face,
            str(temperature),
            str(seed),
            f"{step:.0f}",
            f"{measured['liquid']:.3f}",
            f"{measured['slab']:.3f}",
            f"{measured['lz']:.3f}",
            f"{measured['solid0']:.4f}",
            f"{measured['peak']:.4f}",
            f"{int(measured['nfit'])}",
            fmt(measured["u_liq2_95"]),
            fmt(measured["u_liq2_99"]),
            fmt(measured["u_full"]),
            fmt(measured["u_count"]),
            fmt(measured["u_thick"]),
            fmt(measured["u_crystal"]),
            fmt(measured["u_box2_95"]),
            fmt(measured["u_box1_95"]),
            fmt(measured["u_slab2_95"]),
            fmt(measured["u_one2_95"]),
            fmt(measured["u_liq2_70"]),
            fmt(measured["u_liq2_80"]),
            fmt(measured["u_liq2_90"]),
            fmt(measured["u_q6"]),
            f"{q_t:.3f}" if q_t == q_t else "nan",
            f"{q_gap:.4f}" if q_gap == q_gap else "nan",
            fmt(measured["u_q6_358"]),
            fmt(measured["u_q6_local"]),
            fmt(measured["u_energy"]),
            f"{measured['drop']:.1f}" if measured["drop"] == measured["drop"] else "nan",
            f"{count_ns * 0.1:.3f}",
            f"{count_ns * 10.0:.1f}",
            f"{count_ns / SIGMA_A:.3f}",
            f"{printed:.1f}" if printed == printed else "nan",
            f"{printed / 10.0:.1f}" if printed == printed else "nan",
            f"{printed / 4.0:.1f}" if printed == printed else "nan",
            shi_s,
        ]
        if len(fields) != len(header):
            raise SystemExit(f"width {len(fields)} != {len(header)}")
        lines.append("\t".join(fields))
        if abs(step - 10.0) < 1.0e-6:
            grouped.setdefault((face, temperature, step), []).append(measured)
        print(
            f"{face} T={temperature} s={seed} dt={step:.0f} "
            f"eq2={fmt(measured['u_liq2_95'])} count={fmt(measured['u_count'])} "
            f"q6={fmt(measured['u_q6'])} energy={fmt(measured['u_energy'])} "
            f"box/tau={fmt(measured['u_box1_95'])}",
            flush=True,
        )
    lines.append(f"# cores {os.cpu_count() or 1}")
    lines.append(f"# analysis_seconds {time.perf_counter() - t0:.2f}")
    lines.append("# speeds except the three unit columns are A/ns.")
    lines.append("# u_eq2 is L_liquid/(2 tau). u_count is (dN/dt) Lz/(2 N) on solid fraction 0.55 to 0.90.")
    lines.append("# u_box_over_tau is Lz/tau at solid fraction 0.95, the factor-of-four misreading.")
    lines.append("# u_one_over_2tau divides by two a second time, after taking one interface's travel.")
    lines.append("# u_q6_mislabel uses Lechner-Dellago q6bar, cutoff 3.5 A.")
    lines.append("# The threshold is the low edge of the gap where both mislabel rates are at most 0.5 percent.")
    lines.append("# u_q6_0p358 uses the constant 0.358 quoted for the earlier seeding work.")
    lines.append("# u_energy is (Lz/2)/(2 tau) with tau when the potential has fallen 95 percent of the way to its plateau.")
    lines.append("# 1 A/ns = 0.1 m/s = 10 cm/s. u_sigma_per_ns divides the count slope by sigma = 2.3925 A.")
    lines.append("# printed values are Espinosa et al. Table I. shi_same_face is Figure 3 of Shi et al., J. Chem. Phys. 162, 054714 (2025).")
    for key, block in grouped.items():
        face, temperature, step = key
        for name in ("u_liq2_95", "u_count", "u_thick", "u_q6", "u_q6_358", "u_energy", "u_box1_95", "u_full"):
            samples = [row[name] for row in block if row.get(name) == row.get(name)]
            if len(samples) < 1:
                continue
            if len(samples) == 1:
                mean = samples[0]
                lo = hi = mean
            else:
                mean, lo, hi = bootstrap_mean(samples)
            printed = PUBLISHED.get(temperature)
            over_p = angstrom_per_ns(mean) / printed if printed else float("nan")
            over_c = angstrom_per_ns(mean) / (printed / 10.0) if printed else float("nan")
            over_4 = angstrom_per_ns(mean) / (printed / 4.0) if printed else float("nan")
            lines.append(
                f"# summary {face} dt {step:.0f} T {temperature} {name} "
                f"n {len(samples)} mean {angstrom_per_ns(mean):.3f} "
                f"lo {angstrom_per_ns(lo):.3f} hi {angstrom_per_ns(hi):.3f} "
                f"over_printed {over_p:.3f} over_corrected {over_c:.3f} over_printed_div4 {over_4:.3f}"
            )
    args.out.parent.mkdir(parents=True, exist_ok=True)
    args.out.write_text("\n".join(lines) + "\n")
    if args.out.stat().st_size > 80000:
        raise SystemExit(f"{args.out} exceeds 80 KB")
    print(f"wrote {args.out} ({args.out.stat().st_size} bytes)")


if __name__ == "__main__":
    main()
