# Copyright (c) d-SEAMS developers

"""Count-slope speeds and a Wilson-Frenkel curve for the growth campaign.

The estimator is the one in espinosa_candidates.py: the slope of the
2015 CHILL+ solid count on solid fraction 0.55 to 0.90,

    u = (dN/dt) * Lz / (2 N),

and the Lechner-Dellago q6bar count on the same window. CHILL and
CHILL+ thresholds are not changed. The q6 cut is the mislabel point
already used for Table I, or the constant 0.358. Seeds are summarised
by the same percentile bootstrap.

The curve on each face is

    u(T) = A * exp(-B/T) * (1 - exp(-C * (Tm - T) / T)),

with Tm = 274.6 K from Espinosa et al. Table II. A, B and C are free.
The maximum is read off that curve.

Speeds in the tables are angstroms per nanosecond.

    python3 scripts/growth_curve_rates.py --bin bbdir/tests/interface_order
"""

from __future__ import annotations

import argparse
import math
import os
import time
from concurrent.futures import ProcessPoolExecutor
from pathlib import Path

import espinosa_candidates as cand

ROOT = Path(__file__).resolve().parents[1]
TM = 274.6
FACES = ("basal", "prism1", "prism2")
TEMPERATURES = (220, 230, 240, 250, 260, 265)
# Espinosa et al., J. Chem. Phys. 145, 211922 (2016), Table I, as printed.
# The one-front reading of that column is the printed number divided by 4.
TABLE_I = {240: 490.0, 250: 375.0, 260: 266.0}
# Table II: ln(u / (A/ns)) = C1 + C2*dT + C3*dT^2, dT = Tm - T.
# This fit follows the printed column. Dividing by 4 is the one-front reading.
TABLE_II = (4.1711, 0.096266, -0.001148)
MIN_FIT = 8


def table_ii(temperature: float) -> float:
    c1, c2, c3 = TABLE_II
    delta = TM - temperature
    return math.exp(c1 + c2 * delta + c3 * delta * delta)


def table_ii_maximum() -> tuple[float, float]:
    """Temperature where the printed quadratic in ln(u) is largest."""
    _c1, c2, c3 = TABLE_II
    delta = -c2 / (2.0 * c3)
    return TM - delta, table_ii(TM - delta)


def wilson_frenkel(temperature: float, amp: float, barrier: float, drive: float) -> float:
    if temperature <= 0.0 or temperature >= TM:
        return 0.0
    delta = TM - temperature
    return amp * math.exp(-barrier / temperature) * (1.0 - math.exp(-drive * delta / temperature))


def curve_maximum(amp: float, barrier: float, drive: float) -> tuple[float, float]:
    best_t = TM
    best_u = -1.0
    # 0.05 K is finer than the seed-to-seed scatter of the maximum.
    t = 200.0
    while t < TM:
        speed = wilson_frenkel(t, amp, barrier, drive)
        if speed > best_u:
            best_u = speed
            best_t = t
        t += 0.05
    return best_t, best_u


def _shape(temperature: float, barrier: float, drive: float) -> float:
    delta = TM - temperature
    return math.exp(-barrier / temperature) * (1.0 - math.exp(-drive * delta / temperature))


def fit_curve(temps: list[float], speeds: list[float]) -> tuple[float, float, float, float]:
    """Return A, B, C and the sum of squared residuals, in A/ns and kelvin.

    A is the linear least-squares scale of a fixed (B, C) shape. B and C
    are read from a grid, then from a finer grid around the best cell.
    """
    if len(temps) < 3:
        return float("nan"), float("nan"), float("nan"), float("nan")

    def solve(barrier: float, drive: float) -> tuple[float, float] | None:
        if barrier <= 0.0 or drive <= 0.0:
            return None
        weights = [_shape(temperature, barrier, drive) for temperature in temps]
        den = sum(weight * weight for weight in weights)
        if den <= 0.0:
            return None
        amp = sum(speed * weight for speed, weight in zip(speeds, weights)) / den
        if amp <= 0.0:
            return None
        rss = sum((amp * weight - speed) ** 2 for speed, weight in zip(speeds, weights))
        return amp, rss

    def search(b0: float, b1: float, db: float, c0: float, c1: float, dc: float, best):
        barrier = b0
        while barrier <= b1 + 1.0e-9:
            drive = c0
            while drive <= c1 + 1.0e-9:
                got = solve(barrier, drive)
                if got is not None and (best is None or got[1] < best[3]):
                    best = (got[0], barrier, drive, got[1])
                drive += dc
            barrier += db
        return best

    best = search(200.0, 8000.0, 200.0, 0.2, 8.0, 0.2, None)
    if best is None:
        return float("nan"), float("nan"), float("nan"), float("nan")
    _amp, barrier, drive, _rss = best
    best = search(barrier - 200.0, barrier + 200.0, 20.0, drive - 0.2, drive + 0.2, 0.02, best)
    _amp, barrier, drive, _rss = best
    best = search(max(50.0, barrier - 20.0), barrier + 20.0, 5.0, max(0.05, drive - 0.02), drive + 0.02, 0.005, best)
    return best


def self_test() -> None:
    cand.self_test()
    # The printed quadratic peaks once, below Tm, and /4 does not move it.
    t_max, u_max = table_ii_maximum()
    assert 230.0 < t_max < 236.0
    assert abs(table_ii(240.0) - 458.3) < 0.5
    assert abs((u_max / 4.0) / (table_ii(t_max) / 4.0) - 1.0) < 1.0e-12
    # A Wilson-Frenkel curve with a barrier vanishes at Tm and has one maximum.
    peak_t, peak_u = curve_maximum(300.0, 2000.0, 2.0)
    assert 200.0 < peak_t < TM
    assert peak_u > wilson_frenkel(220.0, 300.0, 2000.0, 2.0) or peak_t < 230.0
    assert wilson_frenkel(TM, 300.0, 2000.0, 2.0) == 0.0
    temps = [220.0, 230.0, 240.0, 250.0, 260.0, 265.0]
    speeds = [wilson_frenkel(t, 280.0, 1600.0, 2.2) for t in temps]
    amp, barrier, drive, loss = fit_curve(temps, speeds)
    assert loss < 1.0
    got_t, got_u = curve_maximum(amp, barrier, drive)
    true_t, true_u = curve_maximum(280.0, 1600.0, 2.2)
    assert abs(got_t - true_t) < 1.0
    assert abs(got_u - true_u) < 2.0


def score_one(task: tuple[str, int, int, float, str, str, str]) -> dict:
    face, temperature, seed, step, dump, log, binary = task
    rows = cand.parse_growth(cand.run_cmd(Path(binary), ["growth", dump, "5.0"]))
    q6_text = cand.run_cmd(Path(binary), ["q6", dump, "5.0"])
    thresholds, q6_rows, ice, liq = cand.parse_q6(q6_text)
    lz = float("nan")
    for line in q6_text.splitlines():
        if line.startswith("frame 1 "):
            lz = float(line.split()[5])
            break
    measured = cand.chill_speeds(rows, lz)
    q_t, q_gap = cand.mislabel_threshold(thresholds, ice, liq)
    q_rows = cand.fraction_series(q6_rows, thresholds, q_t, "bar")
    q_fixed = cand.fraction_series(q6_rows, thresholds, 0.358, "bar")
    measured["u_q6"] = cand.q6_speed(q_rows, lz)
    measured["u_q6_358"] = cand.q6_speed(q_fixed, lz)
    measured["q6_t"] = q_t
    measured["q6_gap"] = q_gap
    measured["face"] = face
    measured["temperature"] = temperature
    measured["seed"] = seed
    measured["step"] = step
    measured["log"] = log
    return measured


def mean_block(samples: list[float]) -> tuple[float, float, float]:
    if not samples:
        return float("nan"), float("nan"), float("nan")
    if len(samples) == 1:
        return samples[0], samples[0], samples[0]
    return cand.bootstrap_mean(samples)


def ns(value: float) -> float:
    return cand.angstrom_per_ns(value)


def fmt(value: float, digits: int = 3) -> str:
    if value != value:
        return "nan"
    return f"{value:.{digits}f}"


def main() -> None:
    self_test()
    parser = argparse.ArgumentParser()
    parser.add_argument("--bin", type=Path, default=ROOT / "bbdir" / "tests" / "interface_order")
    parser.add_argument("--data", type=Path, default=Path(os.environ.get("SI4_DATA", "si4-data")))
    parser.add_argument("--out", type=Path, default=ROOT / "docs" / "orgmode" / "explanation" / "growth-curve.tsv")
    parser.add_argument("--fit", type=Path, default=ROOT / "docs" / "orgmode" / "explanation" / "growth-curve-fit.csv")
    parser.add_argument(
        "--compare",
        type=Path,
        default=ROOT / "docs" / "orgmode" / "explanation" / "growth-curve-compare.csv",
    )
    args = parser.parse_args()
    jobs = []
    for face, temperature, seed, step, path, log in cand.discover(args.data):
        if face not in FACES or temperature not in TEMPERATURES:
            continue
        if abs(step - 10.0) > 1.0e-6:
            continue
        jobs.append((face, temperature, seed, step, str(path), str(log), str(args.bin)))
    if not jobs:
        raise SystemExit(f"no growth-curve trajectories in {args.data}")
    cores = os.cpu_count() or 1
    t0 = time.perf_counter()
    workers = min(cores, len(jobs))
    with ProcessPoolExecutor(max_workers=workers) as pool:
        measured = list(pool.map(score_one, jobs))
    elapsed = time.perf_counter() - t0
    header = [
        "face",
        "temperature",
        "seed",
        "nfit",
        "peak",
        "liquid_A",
        "lz_A",
        "u_count",
        "u_q6",
        "u_q6_0p358",
        "q6_threshold",
        "used",
    ]
    lines = ["\t".join(header)]
    by_face: dict[str, dict[int, list[dict]]] = {face: {} for face in FACES}
    for row in sorted(measured, key=lambda item: (item["face"], item["temperature"], item["seed"])):
        used = int(row["nfit"]) >= MIN_FIT and row["u_count"] == row["u_count"] and row["peak"] >= 0.90
        fields = [
            row["face"],
            str(row["temperature"]),
            str(row["seed"]),
            str(int(row["nfit"])),
            f"{row['peak']:.4f}",
            f"{row['liquid']:.3f}",
            f"{row['lz']:.3f}",
            fmt(ns(row["u_count"])),
            fmt(ns(row["u_q6"])),
            fmt(ns(row["u_q6_358"])),
            f"{row['q6_t']:.3f}" if row["q6_t"] == row["q6_t"] else "nan",
            "yes" if used else "no",
        ]
        lines.append("\t".join(fields))
        if used:
            by_face[row["face"]].setdefault(int(row["temperature"]), []).append(row)
        print(
            f"{row['face']} T={row['temperature']} s={row['seed']} "
            f"count={fmt(ns(row['u_count']))} q6={fmt(ns(row['u_q6']))} "
            f"nfit={int(row['nfit'])} peak={row['peak']:.3f} used={used}",
            flush=True,
        )
    have = {(row["face"], int(row["temperature"]), int(row["seed"])) for row in measured}
    for face in FACES:
        for temperature in TEMPERATURES:
            for seed in (1, 2, 3, 4):
                if (face, temperature, seed) not in have:
                    lines.append(f"# missing {face} T {temperature} seed {seed}")
    lines.append(f"# cores {cores}")
    lines.append(f"# analysis_seconds {elapsed:.2f}")
    lines.append("# speeds are A/ns. u_count is the 2015 CHILL+ solid-count slope.")
    lines.append("# u_q6 is Lechner-Dellago q6bar at the mislabel cut. u_q6_0p358 uses 0.358.")
    lines.append("# used=yes requires at least 8 frames in solid fraction 0.55 to 0.90 and a peak of 0.90.")
    lines.append("# The bootstrap is 10000 draws, Hyndman-Fan type 7, generator seed 270.")

    summary_header = [
        "face",
        "temperature",
        "n",
        "u_count",
        "u_count_lo",
        "u_count_hi",
        "u_q6",
        "u_q6_lo",
        "u_q6_hi",
        "u_q6_0p358",
        "u_q6_0p358_lo",
        "u_q6_0p358_hi",
    ]
    summaries: list[str] = [",".join(summary_header)]
    fit_rows = [
        "face,n_temperatures,A_A_per_ns,B_K,C,rss,T_max_K,u_max,T_max_lo,T_max_hi,u_max_lo,u_max_hi,T_meas_K,u_meas"
    ]
    compare = [
        "face,temperature,n,u_count,u_count_lo,u_count_hi,"
        "table1_printed,table1_over_4,table2_printed,table2_over_4,"
        "shi,shi_lo,shi_hi"
    ]
    for face in FACES:
        temps: list[float] = []
        means: list[float] = []
        blocks: list[list[float]] = []
        for temperature in TEMPERATURES:
            block = by_face[face].get(temperature, [])
            if not block:
                continue
            count_samples = [ns(row["u_count"]) for row in block]
            q6_samples = [ns(row["u_q6"]) for row in block if row["u_q6"] == row["u_q6"]]
            q6_fixed = [ns(row["u_q6_358"]) for row in block if row["u_q6_358"] == row["u_q6_358"]]
            mean, lo, hi = mean_block(count_samples)
            q_mean, q_lo, q_hi = mean_block(q6_samples)
            f_mean, f_lo, f_hi = mean_block(q6_fixed)
            summaries.append(
                ",".join(
                    [
                        face,
                        str(temperature),
                        str(len(count_samples)),
                        fmt(mean),
                        fmt(lo),
                        fmt(hi),
                        fmt(q_mean),
                        fmt(q_lo),
                        fmt(q_hi),
                        fmt(f_mean),
                        fmt(f_lo),
                        fmt(f_hi),
                    ]
                )
            )
            lines.append(
                f"# summary {face} T {temperature} n {len(count_samples)} "
                f"u_count {fmt(mean)} lo {fmt(lo)} hi {fmt(hi)} "
                f"u_q6 {fmt(q_mean)} lo {fmt(q_lo)} hi {fmt(q_hi)} "
                f"u_q6_0p358 {fmt(f_mean)} lo {fmt(f_lo)} hi {fmt(f_hi)}"
            )
            temps.append(float(temperature))
            means.append(mean)
            blocks.append(count_samples)
            printed = TABLE_I.get(temperature, float("nan")) if face == "prism2" else float("nan")
            fitted = table_ii(temperature) if face == "prism2" else float("nan")
            shi = cand.SHI.get((face, temperature))
            compare.append(
                ",".join(
                    [
                        face,
                        str(temperature),
                        str(len(count_samples)),
                        fmt(mean),
                        fmt(lo),
                        fmt(hi),
                        fmt(printed, 1),
                        fmt(printed / 4.0, 1) if printed == printed else "nan",
                        fmt(fitted, 1),
                        fmt(fitted / 4.0, 1) if fitted == fitted else "nan",
                        fmt(shi[0]) if shi else "nan",
                        fmt(shi[1]) if shi else "nan",
                        fmt(shi[2]) if shi else "nan",
                    ]
                )
            )
        amp, barrier, drive, rss = fit_curve(temps, means)
        if amp != amp:
            fit_rows.append(f"{face},{len(temps)},nan,nan,nan,nan,nan,nan,nan,nan,nan,nan,nan,nan")
            continue
        meas_i = max(range(len(means)), key=lambda i: means[i])
        meas_t = temps[meas_i]
        meas_u = means[meas_i]
        peak_t, peak_u = curve_maximum(amp, barrier, drive)
        # Resample seeds inside each temperature, refit, and read the maximum again.
        state = 270
        peak_ts: list[float] = []
        peak_us: list[float] = []
        n_draw = 1000
        for _ in range(n_draw):
            sample_means = []
            ok = True
            for block in blocks:
                acc = 0.0
                n = len(block)
                for _i in range(n):
                    state = cand.xorshift(state)
                    acc += block[state % n]
                sample_means.append(acc / n)
            s_amp, s_bar, s_drv, s_rss = fit_curve(temps, sample_means)
            if s_amp != s_amp or s_rss != s_rss:
                ok = False
            if not ok:
                continue
            s_t, s_u = curve_maximum(s_amp, s_bar, s_drv)
            peak_ts.append(s_t)
            peak_us.append(s_u)
        if peak_ts:
            peak_ts.sort()
            peak_us.sort()
            t_lo = cand.quantile7(peak_ts, 0.025)
            t_hi = cand.quantile7(peak_ts, 0.975)
            u_lo = cand.quantile7(peak_us, 0.025)
            u_hi = cand.quantile7(peak_us, 0.975)
        else:
            t_lo = t_hi = peak_t
            u_lo = u_hi = peak_u
        # The reported maximum is the fit to the seed means. The interval is the bootstrap cloud.
        fit_rows.append(
            ",".join(
                [
                    face,
                    str(len(temps)),
                    fmt(amp, 3),
                    fmt(barrier, 1),
                    fmt(drive, 3),
                    fmt(rss, 2),
                    fmt(peak_t, 2),
                    fmt(peak_u, 3),
                    fmt(t_lo, 2),
                    fmt(t_hi, 2),
                    fmt(u_lo, 3),
                    fmt(u_hi, 3),
                    fmt(meas_t, 0),
                    fmt(meas_u, 3),
                ]
            )
        )
        lines.append(
            f"# wilson {face} nT {len(temps)} A {fmt(amp, 3)} B {fmt(barrier, 1)} C {fmt(drive, 3)} "
            f"rss {fmt(rss, 2)} T_max {fmt(peak_t, 2)} u_max {fmt(peak_u, 3)} "
            f"T_lo {fmt(t_lo, 2)} T_hi {fmt(t_hi, 2)} u_lo {fmt(u_lo, 3)} u_hi {fmt(u_hi, 3)} "
            f"boot {len(peak_ts)}"
        )
    ii_t, ii_u = table_ii_maximum()
    lines.append(
        f"# table2_maximum T {ii_t:.2f} printed {ii_u:.1f} over4 {ii_u / 4.0:.1f}"
    )
    lines.append("# table1_over_4 and table2_over_4 are the one-front reading of the printed mW column.")
    lines.append("# shi is Figure 3 of Shi et al., J. Chem. Phys. 162, 054714 (2025), marker and error-bar caps.")
    args.out.parent.mkdir(parents=True, exist_ok=True)
    # The per-seed rows and the temperature means share growth-curve.tsv.
    # The means are repeated as growth-curve-summary.csv so a sheet can read them.
    summary_path = args.out.with_name("growth-curve-summary.csv")
    args.out.write_text("\n".join(lines) + "\n")
    summary_path.write_text("\n".join(summaries) + "\n")
    args.fit.write_text("\n".join(fit_rows) + "\n")
    args.compare.write_text("\n".join(compare) + "\n")
    for path in (args.out, summary_path, args.fit, args.compare):
        if path.stat().st_size > 80000:
            raise SystemExit(f"{path} exceeds 80 KB")
        print(f"wrote {path.name} ({path.stat().st_size} bytes)")


if __name__ == "__main__":
    main()
