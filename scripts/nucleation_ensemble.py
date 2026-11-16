#!/usr/bin/env python3
# Copyright (c) d-SEAMS developers

"""Bootstrap the ring-share lead across seeded liquid quenches.

Each seed directory holds frames.tsv and summary.tsv from
nucleation-pathway. The lead is the half-rise lag already printed in
the summary: ice half-rise minus six-ring-share half-rise, in
picoseconds. A positive lag means the ring share moves first.

The mean lag gets a percentile bootstrap interval. The interval
excludes zero when both ends sit on the same side of zero.

The post-critical stacking fraction is the cubic share inside the
largest cluster at the half-rise of that cluster. It is reported only
when the run starts as a liquid and that label's own size series
climbs. The barrier is the size at the foot of the final strict climb.
A series that never climbs does not enter the basin. It is not an
MFPT critical size.

The per-frame largest-cluster columns are copied to one table so a
later mean-first-passage fit can reuse them.
"""

from __future__ import annotations

import argparse
import math
import random
from pathlib import Path

# A starting frame is a liquid when the published core cluster is at
# most this size, the published crystal population is at most this
# count, and the six-ring share is still on the liquid side of the
# figshare window (which opened at 0.552 with a 59-molecule cluster).
LIQUID_CORE_MAX = 5
LIQUID_CRYSTAL_MAX = 20
LIQUID_SHARE_MAX = 0.50

N_BOOT = 10000
BOOT_SEED = 20261111

LEAD_NAMES = (
    "chill_core",
    "chill_ext",
    "ovito_core",
    "ovito_ext",
    "q6bar",
    "symmetry",
)

# Requested comparisons. CHILL+ is the published rule, both masks.
PRIMARY = ("chill_core", "chill_ext", "q6bar", "symmetry")

CLUSTER_COLUMNS = (
    "n_chill_core",
    "c_chill_core",
    "h_chill_core",
    "n_chill_ext",
    "c_chill_ext",
    "h_chill_ext",
    "n_ovito_core",
    "c_ovito_core",
    "h_ovito_core",
    "n_ovito_ext",
    "c_ovito_ext",
    "h_ovito_ext",
    "n_q6",
    "c_q6",
    "h_q6",
    "n_sym",
    "c_sym",
    "h_sym",
    "ring5",
    "ring6",
    "chill_cubic",
    "chill_hex",
    "chill_inter",
    "q6bar_ice",
    "sym_cubic",
    "sym_hex",
)


def load_tsv(path: Path) -> list[dict[str, str]]:
    lines = path.read_text().splitlines()
    if not lines:
        return []
    header = lines[0].split("\t")
    rows = []
    for line in lines[1:]:
        if not line.strip():
            continue
        parts = line.split("\t")
        rows.append(dict(zip(header, parts)))
    return rows


def parse_summary(path: Path) -> dict[str, dict[str, str]]:
    out: dict[str, dict[str, str]] = {}
    for line in path.read_text().splitlines():
        parts = line.split("\t")
        if len(parts) < 2:
            continue
        kind, name = parts[0], parts[1]
        fields = {parts[i]: parts[i + 1] for i in range(2, len(parts) - 1, 2)}
        out[f"{kind}:{name}"] = fields
    return out


def ring_share(row: dict[str, str]) -> float:
    n5 = int(row["ring5"])
    n6 = int(row["ring6"])
    den = n5 + n6
    return n6 / den if den else 0.0


def is_liquid(row: dict[str, str]) -> bool:
    core = int(row["n_chill_core"])
    crystal = int(row["chill_cubic"]) + int(row["chill_hex"])
    return (
        core <= LIQUID_CORE_MAX
        and crystal <= LIQUID_CRYSTAL_MAX
        and ring_share(row) < LIQUID_SHARE_MAX
    )


def percentile(sorted_values: list[float], p: float) -> float:
    if not sorted_values:
        return float("nan")
    if len(sorted_values) == 1:
        return sorted_values[0]
    k = (len(sorted_values) - 1) * (p / 100.0)
    lo = math.floor(k)
    hi = math.ceil(k)
    if lo == hi:
        return sorted_values[lo]
    w = k - lo
    return sorted_values[lo] * (1.0 - w) + sorted_values[hi] * w


def bootstrap_mean(values: list[float], n_boot: int, rng_seed: int) -> dict[str, float]:
    rng = random.Random(rng_seed)
    n = len(values)
    if n == 0:
        return {
            "n": 0,
            "mean": float("nan"),
            "lo": float("nan"),
            "hi": float("nan"),
            "excludes_zero": 0,
        }
    means = []
    for _ in range(n_boot):
        draw = 0.0
        for _i in range(n):
            draw += values[rng.randrange(n)]
        means.append(draw / n)
    means.sort()
    lo = percentile(means, 2.5)
    hi = percentile(means, 97.5)
    excludes = (lo > 0.0 and hi > 0.0) or (lo < 0.0 and hi < 0.0)
    return {
        "n": n,
        "mean": sum(values) / n,
        "lo": lo,
        "hi": hi,
        "excludes_zero": int(excludes),
    }


def grow_barrier(series: list[float]) -> float | None:
    """Size at the foot of the final strict climb, or None if it never climbs."""
    if not series:
        return None
    imax = max(range(len(series)), key=lambda i: series[i])
    foot = imax
    while foot > 0 and series[foot - 1] < series[foot]:
        foot -= 1
    if not series[imax] > series[foot]:
        return None
    return series[foot]


SIZE_COLUMN = {
    "chill_core": "n_chill_core",
    "chill_ext": "n_chill_ext",
    "ovito_core": "n_ovito_core",
    "ovito_ext": "n_ovito_ext",
    "q6bar": "n_q6",
    "symmetry": "n_sym",
}


def fnum(text: str) -> float:
    try:
        return float(text)
    except (TypeError, ValueError):
        return float("nan")


def collect(seed_dirs: list[Path]) -> tuple[list[dict], list[str], list[dict]]:
    leads = []
    cluster_lines = []
    stacking = []
    header = (
        "seed\tframe\ttimestep\ttime_ps\tring_share\t"
        + "\t".join(CLUSTER_COLUMNS)
    )
    cluster_lines.append(header)
    for seed_dir in seed_dirs:
        seed = seed_dir.name
        if seed.startswith("seed-"):
            seed = seed[len("seed-"):]
        frames = load_tsv(seed_dir / "frames.tsv")
        summary = parse_summary(seed_dir / "summary.tsv")
        if not frames:
            continue
        ps = 1.0
        for line in (seed_dir / "summary.tsv").read_text().splitlines():
            parts = line.split("\t")
            if parts[:1] == ["ps_per_frame"] and len(parts) == 2:
                ps = float(parts[1])
        t0 = int(frames[0]["frame"])
        liquid = is_liquid(frames[0])
        max_core = max(int(row["n_chill_core"]) for row in frames)
        final_core = int(frames[-1]["n_chill_core"])
        barriers = {
            name: grow_barrier([float(row[column]) for row in frames])
            for name, column in SIZE_COLUMN.items()
        }
        # The core occupies its basin when its own series climbs and the
        # last frame is still at or above that barrier.
        core_barrier = barriers["chill_core"]
        grown = core_barrier is not None and final_core >= core_barrier
        rises = {
            "chill_core": max(int(r["chill_cubic"]) + int(r["chill_hex"]) for r in frames)
            - (int(frames[0]["chill_cubic"]) + int(frames[0]["chill_hex"])),
            "chill_ext": max(
                int(r["chill_cubic"]) + int(r["chill_hex"]) + int(r["chill_inter"]) for r in frames
            )
            - (
                int(frames[0]["chill_cubic"])
                + int(frames[0]["chill_hex"])
                + int(frames[0]["chill_inter"])
            ),
            "ovito_core": max(int(r["ovito_cubic"]) + int(r["ovito_hex"]) for r in frames)
            - (int(frames[0]["ovito_cubic"]) + int(frames[0]["ovito_hex"])),
            "ovito_ext": max(
                int(r["ovito_cubic"]) + int(r["ovito_hex"]) + int(r["ovito_inter"]) for r in frames
            )
            - (
                int(frames[0]["ovito_cubic"])
                + int(frames[0]["ovito_hex"])
                + int(frames[0]["ovito_inter"])
            ),
            "q6bar": max(int(r["q6bar_ice"]) for r in frames) - int(frames[0]["q6bar_ice"]),
            "symmetry": max(int(r["sym_cubic"]) + int(r["sym_hex"]) for r in frames)
            - (int(frames[0]["sym_cubic"]) + int(frames[0]["sym_hex"])),
        }
        share_rise = max(ring_share(r) for r in frames) - ring_share(frames[0])
        for row in frames:
            time_ps = (int(row["frame"]) - t0) * ps
            values = [seed, row["frame"], row["timestep"], f"{time_ps:.4f}", f"{ring_share(row):.6f}"]
            values.extend(row[col] for col in CLUSTER_COLUMNS)
            cluster_lines.append("\t".join(values))
        one = {
            "seed": seed,
            "liquid": int(liquid),
            "grown": int(grown),
            "frame1_core": int(frames[0]["n_chill_core"]),
            "frame1_crystal": int(frames[0]["chill_cubic"]) + int(frames[0]["chill_hex"]),
            "frame1_share": ring_share(frames[0]),
            "max_core": max_core,
            "final_core": final_core,
            "frames": len(frames),
            "share_rise": share_rise,
        }
        # Population half-rises are kept for a liquid whose ring share
        # and that ice label both move. Diamond needs a larger rise than
        # a handful of labels, which otherwise makes a one-molecule lag.
        minimum = {
            "chill_core": 20,
            "chill_ext": 20,
            "ovito_core": 20,
            "ovito_ext": 20,
            "q6bar": 20,
            "symmetry": 20,
        }
        for name in LEAD_NAMES:
            fields = summary.get(f"ring_lead:{name}", {})
            ok = (
                fields.get("ok") == "1"
                and liquid
                and share_rise >= 0.02
                and rises[name] >= minimum[name]
            )
            one[f"{name}_ok"] = int(ok)
            one[f"{name}_ps"] = fnum(fields.get("half_ps")) if ok else float("nan")
            one[f"{name}_frames"] = fnum(fields.get("half_frames")) if ok else float("nan")
        cross_ps = float("nan")
        cross_stack = float("nan")
        cross_size = float("nan")
        if core_barrier is not None:
            for row in frames:
                if float(row["n_chill_core"]) >= core_barrier:
                    cross_ps = (int(row["frame"]) - t0) * ps
                    cubic = int(row["c_chill_core"])
                    hexagonal = int(row["h_chill_core"])
                    den = cubic + hexagonal
                    cross_stack = cubic / den if den else float("nan")
                    cross_size = int(row["n_chill_core"])
                    break
        for name in ("chill_core", "chill_ext", "ovito_core", "ovito_ext", "q6bar", "symmetry"):
            fields = summary.get(f"nucleus:{name}", {})
            barrier = barriers[name]
            half_ok = fields.get("half_ok") == "1" and liquid and barrier is not None
            stacking.append(
                {
                    "seed": seed,
                    "definition": name,
                    "liquid": int(liquid),
                    "grown": int(grown and liquid),
                    "half_ok": int(half_ok),
                    "half_size": fnum(fields.get("half_size")) if half_ok else float("nan"),
                    "half_stack": fnum(fields.get("half_stack")) if half_ok else float("nan"),
                    "half_ps": fnum(fields.get("half_ps")) if half_ok else float("nan"),
                    "max_size": fnum(fields.get("max_size")) if liquid else float("nan"),
                    "max_stack": fnum(fields.get("max_stack")) if liquid and half_ok else float("nan"),
                    "cross50_ps": cross_ps if name == "chill_core" else float("nan"),
                    "cross50_size": cross_size if name == "chill_core" else float("nan"),
                    "cross50_stack": cross_stack if name == "chill_core" else float("nan"),
                }
            )
        leads.append(one)
    return leads, cluster_lines, stacking


def write_leads(path: Path, leads: list[dict]) -> None:
    cols = [
        "seed",
        "liquid",
        "grown",
        "frame1_core",
        "frame1_crystal",
        "frame1_share",
        "max_core",
        "final_core",
        "frames",
        "share_rise",
    ]
    for name in LEAD_NAMES:
        cols.extend((f"{name}_ok", f"{name}_frames", f"{name}_ps"))
    lines = ["\t".join(cols)]
    for row in leads:
        cells = []
        for col in cols:
            value = row[col]
            if isinstance(value, float):
                cells.append(f"{value:.6f}" if math.isfinite(value) else "nan")
            else:
                cells.append(str(value))
        lines.append("\t".join(cells))
    path.write_text("\n".join(lines) + "\n")


def write_bootstrap(path: Path, leads: list[dict]) -> list[dict]:
    rows = []
    lines = [
        "signal\tn\tmean_ps\tmedian_ps\tmin_ps\tmax_ps\tci95_lo\tci95_hi\texcludes_zero\tn_positive\tn_negative\tn_zero"
    ]
    for name in PRIMARY:
        values = [
            row[f"{name}_ps"]
            for row in leads
            if row[f"{name}_ok"] and math.isfinite(row[f"{name}_ps"])
        ]
        stat = bootstrap_mean(values, N_BOOT, BOOT_SEED)
        n_pos = sum(1 for v in values if v > 0)
        n_neg = sum(1 for v in values if v < 0)
        n_zero = sum(1 for v in values if v == 0)
        ordered = sorted(values)
        median = percentile(ordered, 50) if ordered else float("nan")
        vmin = ordered[0] if ordered else float("nan")
        vmax = ordered[-1] if ordered else float("nan")
        lines.append(
            f"{name}\t{int(stat['n'])}\t{stat['mean']:.6f}\t{median:.6f}\t{vmin:.6f}\t{vmax:.6f}\t"
            f"{stat['lo']:.6f}\t{stat['hi']:.6f}\t{int(stat['excludes_zero'])}\t{n_pos}\t{n_neg}\t{n_zero}"
        )
        rows.append({"signal": name, **stat, "n_positive": n_pos})
    path.write_text("\n".join(lines) + "\n")
    return rows


def write_stacking(path: Path, stacking: list[dict]) -> None:
    cols = [
        "seed",
        "definition",
        "liquid",
        "grown",
        "half_ok",
        "half_size",
        "half_stack",
        "half_ps",
        "max_size",
        "max_stack",
        "cross50_ps",
        "cross50_size",
        "cross50_stack",
    ]
    lines = ["\t".join(cols)]
    for row in stacking:
        cells = []
        for col in cols:
            value = row[col]
            if isinstance(value, float):
                cells.append(f"{value:.6f}" if math.isfinite(value) else "nan")
            else:
                cells.append(str(value))
        lines.append("\t".join(cells))
    path.write_text("\n".join(lines) + "\n")


def plot_leads(path: Path, leads: list[dict], boot: list[dict]) -> None:
    import matplotlib

    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    fig, ax = plt.subplots(figsize=(7.2, 4.2))
    labels = []
    for i, name in enumerate(PRIMARY):
        values = [
            row[f"{name}_ps"]
            for row in leads
            if row[f"{name}_ok"] and math.isfinite(row[f"{name}_ps"])
        ]
        xs = [i + (k - (len(values) - 1) / 2) * 0.03 for k in range(len(values))]
        ax.scatter(xs, values, s=18, color="black", zorder=3)
        stat = next(row for row in boot if row["signal"] == name)
        ax.plot([i, i], [stat["lo"], stat["hi"]], color="C0", lw=2, zorder=2)
        ax.plot(i, stat["mean"], marker="D", color="C0", zorder=4)
        labels.append(name.replace("_", " "))
    ax.axhline(0.0, color="0.5", lw=0.8)
    ax.set_xticks(range(len(PRIMARY)))
    ax.set_xticklabels(labels)
    ax.set_ylabel("ring-share lead / ps")
    ax.set_title("Ice half-rise minus six-ring-share half-rise")
    fig.tight_layout()
    path.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(path, dpi=140)


def self_test() -> None:
    same = bootstrap_mean([10.0] * 8, 2000, 1)
    if not (same["lo"] > 9.9 and same["hi"] < 10.1 and same["excludes_zero"] == 1):
        raise SystemExit(f"constant sample failed: {same}")
    split = bootstrap_mean([-5.0, 5.0] * 10, 4000, 2)
    if split["excludes_zero"] != 0:
        raise SystemExit(f"balanced sample should include zero: {split}")
    again = bootstrap_mean([1.0, 2.0, 3.0, 4.0], 1000, 7)
    once = bootstrap_mean([1.0, 2.0, 3.0, 4.0], 1000, 7)
    if again["lo"] != once["lo"] or again["hi"] != once["hi"]:
        raise SystemExit("bootstrap is not repeatable")
    print("self-test ok")


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("seed_dirs", nargs="*", type=Path)
    parser.add_argument("--out-dir", type=Path, required=False)
    parser.add_argument("--figure", type=Path, default=None)
    parser.add_argument("--self-test", action="store_true")
    args = parser.parse_args()
    if args.self_test:
        self_test()
        if not args.seed_dirs:
            return
    if not args.out_dir:
        raise SystemExit("--out-dir is required unless --self-test is used alone")
    args.out_dir.mkdir(parents=True, exist_ok=True)
    leads, cluster_lines, stacking = collect(args.seed_dirs)
    write_leads(args.out_dir / "leads.tsv", leads)
    (args.out_dir / "clusters.tsv").write_text("\n".join(cluster_lines) + "\n")
    boot = write_bootstrap(args.out_dir / "bootstrap.tsv", leads)
    write_stacking(args.out_dir / "stacking.tsv", stacking)
    if args.figure:
        plot_leads(args.figure, leads, boot)
    print(args.out_dir / "bootstrap.tsv")


if __name__ == "__main__":
    main()
