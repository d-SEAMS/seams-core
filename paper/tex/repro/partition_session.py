#!/usr/bin/env python3
"""One session of the block cut and the linkcell cell-grid cut.

Strong scaling uses 110592 and 884736 molecules at 1, 2, 4 and 8 ranks.
Weak scaling keeps 110592 molecules on each rank: 24x24x24, 24x24x48,
24x48x48 and 48x48x48 cells. Interface slabs are ice/water at 110592
and 884736 molecules and hydrate/liquid at 105984 and 883200. One
thread per rank. Both cuts run in one process: one warm-up of each,
then twelve block/cell pairs. Ranks are bound to cores.

The host is one Intel Xeon socket in a KVM guest: 8 cores, 1 thread
per core, AVX-512, 320 MiB of L3, and MemTotal 49342816 kB. The guest
reports 2400 MHz and one NUMA node.

    python3 paper/tex/repro/partition_session.py --crystal

Pass no flag to include the interface slabs as well. --crystal is the
uniform-crystal session: strong scaling at both sizes and weak scaling.

The benchmark line is:

    OMP_NUM_THREADS=1 RAYON_NUM_THREADS=1 OMP_WAIT_POLICY=PASSIVE \\
    OMP_DYNAMIC=FALSE OMP_PROC_BIND=TRUE OMP_PLACES=cores \\
    OMPI_MCA_btl=^openib OMPI_MCA_pml=ob1 \\
    mpirun --bind-to core --map-by core:PE=1 --rank-by core -n N \\
        build-mpi/tests/rank_bench 48 12 rings
"""

from __future__ import annotations

import csv
import os
import subprocess
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[3]
DATA = ROOT / "paper" / "tex" / "figures" / "data"
OUT = ROOT / "paper" / "tex" / "figures"
CSV_PATH = DATA / "partition-session.csv"
TABLE_PATH = DATA / "tab-partition.tex"
PROSE_PATH = DATA / "partition-prose.tex"
BENCH = ROOT / "build-mpi" / "tests" / "rank_bench"
BUILD_SO = ROOT / "build-mpi" / "src" / "libyodaLib.so"
LIBS = ":".join(
    [
        str(ROOT / "build-mpi" / "src"),
        str(ROOT / "build-mpi" / "subprojects" / "readcon-core"),
        str(ROOT / "build-mpi" / "subprojects" / "readcon-db"),
        str(ROOT / "build-mpi" / "subprojects" / "linkcell"),
        str(ROOT / "build-mpi" / "subprojects" / "vesin"),
        str(ROOT / "build-mpi" / "subprojects" / "minimage"),
        "/usr/lib/x86_64-linux-gnu/openmpi/lib",
        "/cvmfs/software.eessi.io/versions/2026.06/software/linux/x86_64/intel/"
        "sapphirerapids/software/FlexiBLAS/3.5.0-GCC-15.2.0/lib",
    ]
)
OPS = ("neigh", "chill", "q6", "rings")
OP_LABEL = {"neigh": "neighbours", "chill": "CHILL+", "q6": "$q_6$", "rings": "rings"}
STRONG = ((24, 24, 24, 110592), (48, 48, 48, 884736))
WEAK = ((1, 24, 24, 24), (2, 24, 24, 48), (4, 24, 48, 48), (8, 48, 48, 48))
REPS = 12
COLOURS = {"neigh": "#0072B2", "chill": "#D55E00", "q6": "#009E73", "rings": "#CC79A7"}


def quantile_linear(ordered, q):
    if len(ordered) == 1:
        return ordered[0]
    pos = q * (len(ordered) - 1)
    lo = int(pos)
    frac = pos - lo
    if lo + 1 >= len(ordered):
        return ordered[-1]
    return ordered[lo] * (1.0 - frac) + ordered[lo + 1] * frac


def summarize(values):
    ordered = sorted(values)
    mid = len(ordered) // 2
    if len(ordered) % 2 == 0:
        median = 0.5 * (ordered[mid - 1] + ordered[mid])
    else:
        median = ordered[mid]
    q1 = quantile_linear(ordered, 0.25)
    q3 = quantile_linear(ordered, 0.75)
    return median, q1, q3, q3 - q1


def run_bench(argv, ranks, label, threads=1, reps=REPS):
    env = os.environ.copy()
    env.pop("LD_PRELOAD", None)
    env.pop("SEAMS_RING_PROFILE", None)
    env.pop("SEAMS_FORCE_CUT", None)
    env["LD_LIBRARY_PATH"] = LIBS
    env["OMP_NUM_THREADS"] = str(threads)
    env["RAYON_NUM_THREADS"] = "1"
    env["OMP_WAIT_POLICY"] = "PASSIVE"
    env["OMP_DYNAMIC"] = "FALSE"
    env["OMP_PROC_BIND"] = "TRUE"
    env["OMP_PLACES"] = "cores"
    env["OMPI_MCA_btl"] = "^openib"
    env["OMPI_MCA_pml"] = "ob1"
    env["SEAMS_CUT_LOG"] = "1"
    cmd = [
        "mpirun",
        "--bind-to",
        "core",
        "--map-by",
        "core:PE=1" if threads == 1 else f"socket:PE={threads}",
        "--rank-by",
        "core",
        "-n",
        str(ranks),
        str(BENCH),
        *argv,
    ]
    proc = subprocess.run(cmd, env=env, capture_output=True, text=True, cwd=ROOT / "input")
    with (DATA / "partition-raw.txt").open("a") as raw:
        raw.write(f"# cores {os.cpu_count()} ranks {ranks} threads {threads} {label}\n")
        raw.write(proc.stdout)
        if proc.stderr:
            raw.write(proc.stderr)
    if proc.returncode != 0:
        sys.stderr.write(proc.stderr)
        raise SystemExit(f"benchmark failed ({proc.returncode}): {' '.join(cmd)}")
    got = {
        "block": {"walls": [], "halos": [], "imbalances": []},
        "cell": {"walls": [], "halos": [], "imbalances": []},
    }
    ratios = []
    rules = []
    for line in proc.stdout.splitlines():
        parts = line.split()
        if len(parts) == 6 and parts[0] in got:
            got[parts[0]]["walls"].append(float(parts[1]))
            got[parts[0]]["halos"].append(float(parts[4]))
            got[parts[0]]["imbalances"].append(float(parts[5]))
    for line in proc.stderr.splitlines():
        if line.startswith("cut "):
            parts = line.split()
            if "ratio" in parts:
                ratios.append(float(parts[parts.index("ratio") + 1]))
            if "rule" in parts:
                rules.append(parts[parts.index("rule") + 1])
    for kind in ("block", "cell"):
        if len(got[kind]["walls"]) < reps:
            sys.stderr.write(proc.stderr)
            raise SystemExit(f"expected {reps} {kind} samples, got {len(got[kind]['walls'])}")
    got["ratio"] = max(ratios) if ratios else 1.0
    got["rule"] = rules[-1] if rules else ""
    got["threads"] = threads
    return got


def row_of(kind, scale, molecules, ranks, op, pack, ratio):
    median, q1, q3, iqr = summarize(pack["walls"])
    halo, _, _, _ = summarize(pack["halos"])
    imbalance, _, _, _ = summarize(pack["imbalances"])
    return {
        "cut": kind,
        "scale": scale,
        "molecules": molecules,
        "ranks": ranks,
        "operation": op,
        "median_ms": f"{median:.6f}",
        "q1_ms": f"{q1:.6f}",
        "q3_ms": f"{q3:.6f}",
        "iqr_ms": f"{iqr:.6f}",
        "halo_fraction": f"{halo:.6f}",
        "imbalance": f"{imbalance:.6f}",
        "block_ratio": f"{ratio:.6f}",
        "rule": "",
        "threads": "1",
        "speedup": "",
        "efficiency": "",
    }


def fill_efficiency(rows):
    base = {}
    for row in rows:
        if int(row["ranks"]) == 1:
            base[(row["cut"], row["scale"], row["operation"], int(row["molecules"]))] = float(
                row["median_ms"]
            )
    for row in rows:
        atoms = int(row["molecules"])
        med = float(row["median_ms"])
        n = int(row["ranks"])
        threads = int(row.get("threads") or 1)
        workers = n * threads
        if row["scale"] == "hybrid":
            origin = base[(row["cut"], "strong", row["operation"], atoms)]
            row["speedup"] = f"{origin / med:.6f}"
            row["efficiency"] = f"{origin / (med * workers):.6f}"
        elif row["scale"] in ("strong", "interface", "nucleation"):
            origin = base[(row["cut"], row["scale"], row["operation"], atoms)]
            row["speedup"] = f"{origin / med:.6f}"
            row["efficiency"] = f"{origin / (med * workers):.6f}"
        else:
            key = (row["cut"], row["scale"], row["operation"], 110592 if row["scale"] == "weak" else atoms)
            origin = base[key]
            row["speedup"] = f"{origin / med:.6f}"
            row["efficiency"] = f"{origin / med:.6f}"
        row["molecules"] = atoms


def add_pair(rows, scale, molecules, ranks, op, got):
    for kind in ("block", "cell"):
        row = row_of(kind, scale, molecules, ranks, op, got[kind], got["ratio"])
        row["rule"] = got.get("rule", "")
        row["threads"] = str(got.get("threads", 1))
        rows.append(row)


def measure():
    DATA.mkdir(parents=True, exist_ok=True)
    (DATA / "partition-raw.txt").write_text("")
    rows = []
    for mx, my, mz, molecules in STRONG:
        for ranks in (1, 2, 4, 8):
            for op in OPS:
                print(f"strong {molecules} n={ranks} {op}", flush=True)
                got = run_bench([str(mx), str(my), str(mz), str(REPS), op], ranks, f"strong {molecules} {op}")
                add_pair(rows, "strong", molecules, ranks, op, got)
                dump(rows)
    for ranks, mx, my, mz in WEAK:
        molecules = mx * my * mz * 8
        for op in OPS:
            print(f"weak {molecules} n={ranks} {op}", flush=True)
            got = run_bench([str(mx), str(my), str(mz), str(REPS), op], ranks, f"weak {molecules} {op}")
            add_pair(rows, "weak", molecules, ranks, op, got)
            dump(rows)
    if "--crystal" not in sys.argv:
        for kind, which, molecules in (
            ("ice", 110, 110592),
            ("ice", 885, 884736),
            ("hydrate", 110, 105984),
            ("hydrate", 885, 883200),
        ):
            for ranks in (1, 2, 4, 8):
                for op in OPS:
                    print(f"interface {kind} {molecules} n={ranks} {op}", flush=True)
                    got = run_bench(
                        ["slab", kind, str(which), str(REPS), op],
                        ranks,
                        f"interface {kind} {molecules} {op}",
                    )
                    add_pair(rows, "interface", molecules, ranks, op, got)
                    dump(rows)
    fill_efficiency(rows)
    dump(rows)
    return rows


def dump(rows):
    DATA.mkdir(parents=True, exist_ok=True)
    fields = [
        "cut",
        "scale",
        "molecules",
        "ranks",
        "operation",
        "median_ms",
        "q1_ms",
        "q3_ms",
        "iqr_ms",
        "speedup",
        "efficiency",
        "halo_fraction",
        "imbalance",
        "block_ratio",
        "rule",
        "threads",
    ]
    with CSV_PATH.open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields)
        writer.writeheader()
        writer.writerows(rows)


def pick(rows, scale, molecules, op, cut, ranks):
    for row in rows:
        if (
            row["scale"] == scale
            and int(row["molecules"]) == molecules
            and row["operation"] == op
            and row["cut"] == cut
            and int(row["ranks"]) == ranks
        ):
            return row
    raise KeyError((scale, molecules, op, cut, ranks))


def write_table(rows):
    lines = [
        "\\begin{table}[tbp]",
        "\\centering",
        "\\caption{\\textbf{Strong scaling of the block cut and the cell-grid cut.}",
        "Block and cell alternate in one loop, twelve times, ranks bound to cores.",
        "Efficiency is the one-rank time divided by the eight-rank time and by eight.",
        "IQR is the type-7 interquartile range of the eight-rank walls.",
        "A gap smaller than the two ranges is not a difference.}",
        "\\label{tab:partition-strong}",
        "\\begin{tabular}{llrrrrrrr}",
        "\\toprule",
        "Molecules & Operation & Block 8 & Block IQR & Cell 8 & Cell IQR & Cell eff. & Imbalance & Halo \\\\",
        "\\midrule",
    ]
    last = None
    for molecules in (110592, 884736):
        for op in OPS:
            block = pick(rows, "strong", molecules, op, "block", 8)
            cell = pick(rows, "strong", molecules, op, "cell", 8)
            shown = str(molecules) if molecules != last else ""
            last = molecules
            lines.append(
                f"{shown} & {OP_LABEL[op]} & {float(block['median_ms']):.3f} & {float(block['iqr_ms']):.3f} "
                f"& {float(cell['median_ms']):.3f} & {float(cell['iqr_ms']):.3f} "
                f"& {float(cell['efficiency']):.3f} & {float(cell['imbalance']):.3f} "
                f"& {float(cell['halo_fraction']):.3f} \\\\"
            )
    lines.extend(["\\bottomrule", "\\end{tabular}", "\\end{table}", ""])
    lines.extend(
        [
            "\\begin{table}[tbp]",
            "\\centering",
            "\\caption{\\textbf{Weak scaling at 110592 molecules per rank.}",
            "Block and cell alternate in one loop. IQR is the type-7 range of the eight-rank walls.}",
            "\\label{tab:partition-weak}",
            "\\begin{tabular}{lrrrrrr}",
            "\\toprule",
            "Operation & Block 8 & Block IQR & Cell 8 & Cell IQR & Cell eff. & Imbalance \\\\",
            "\\midrule",
        ]
    )
    for op in OPS:
        block = pick(rows, "weak", 884736, op, "block", 8)
        cell = pick(rows, "weak", 884736, op, "cell", 8)
        lines.append(
            f"{OP_LABEL[op]} & {float(block['median_ms']):.3f} & {float(block['iqr_ms']):.3f} "
            f"& {float(cell['median_ms']):.3f} & {float(cell['iqr_ms']):.3f} "
            f"& {float(cell['efficiency']):.3f} & {float(cell['imbalance']):.3f} \\\\"
        )
    lines.extend(["\\bottomrule", "\\end{tabular}", "\\end{table}", ""])
    if not any(row["scale"] == "interface" for row in rows):
        TABLE_PATH.write_text("\n".join(lines))
        return
    lines.extend(
        [
            "\\begin{table}[tbp]",
            "\\centering",
            "\\caption{\\textbf{Interface slabs, eight ranks, block and cell in one loop.}",
            "The block-load ratio is the measured work on the block ranks.",
            "At or above 1.025 the automatic rule selects the cell-grid cut.}",
            "\\label{tab:partition-interface}",
            "\\begin{tabular}{llrrrrrr}",
            "\\toprule",
            "Frame & Operation & Block 8 & Block IQR & Cell 8 & Cell IQR & Ratio & Imbalance \\\\",
            "\\midrule",
        ]
    )
    last = None
    for molecules, frame in (
        (110592, "ice/water"),
        (884736, "ice/water"),
        (105984, "hydrate/liquid"),
        (883200, "hydrate/liquid"),
    ):
        for op in OPS:
            block = pick(rows, "interface", molecules, op, "block", 8)
            cell = pick(rows, "interface", molecules, op, "cell", 8)
            mark = f"{frame} {molecules}"
            shown = mark if mark != last else ""
            last = mark
            lines.append(
                f"{shown} & {OP_LABEL[op]} & {float(block['median_ms']):.3f} & {float(block['iqr_ms']):.3f} "
                f"& {float(cell['median_ms']):.3f} & {float(cell['iqr_ms']):.3f} "
                f"& {float(cell['block_ratio']):.3f} & {float(cell['imbalance']):.3f} \\\\"
            )
    lines.extend(["\\bottomrule", "\\end{tabular}", "\\end{table}", ""])
    TABLE_PATH.write_text("\n".join(lines))


def ranges_overlap(a, b):
    return not (float(a["q3_ms"]) < float(b["q1_ms"]) or float(b["q3_ms"]) < float(a["q1_ms"]))


def write_prose(rows):
    bits = []
    for molecules in (110592, 884736):
        for op in OPS:
            cell = pick(rows, "strong", molecules, op, "cell", 8)
            block = pick(rows, "strong", molecules, op, "block", 8)
            ratio = float(cell["block_ratio"])
            choice = "the cell-grid cut" if ratio >= 1.025 else "the block cut"
            verdict = "The ranges overlap." if ranges_overlap(block, cell) else (
                "The ranges do not overlap."
            )
            bits.append(
                f"At {molecules} molecules the eight-rank {OP_LABEL[op]} wall is "
                f"{float(block['median_ms']):.3f} ms on the block cut "
                f"(IQR {float(block['iqr_ms']):.3f}, imbalance {float(block['imbalance']):.3f}) "
                f"and {float(cell['median_ms']):.3f} ms on the cell-grid cut "
                f"(IQR {float(cell['iqr_ms']):.3f}, imbalance {float(cell['imbalance']):.3f}, "
                f"halo {float(cell['halo_fraction']):.3f}). "
                f"The block-load ratio is {ratio:.3f}, so the automatic rule selects {choice}. "
                f"{verdict}"
            )
    interface_frames = (
        (
            (110592, "ice/water"),
            (884736, "ice/water"),
            (105984, "hydrate/liquid"),
            (883200, "hydrate/liquid"),
        )
        if any(row["scale"] == "interface" for row in rows)
        else ()
    )
    for molecules, frame in interface_frames:
        for op in OPS:
            cell = pick(rows, "interface", molecules, op, "cell", 8)
            block = pick(rows, "interface", molecules, op, "block", 8)
            ratio = float(cell["block_ratio"])
            choice = "the cell-grid cut" if ratio >= 1.025 else "the block cut"
            verdict = "The ranges overlap." if ranges_overlap(block, cell) else (
                "The ranges do not overlap."
            )
            bits.append(
                f"On the {frame} slab of {molecules} molecules the block-load ratio for "
                f"{OP_LABEL[op]} is {ratio:.3f}, so the automatic rule selects {choice}. "
                f"Eight ranks: block {float(block['median_ms']):.3f} ms "
                f"(IQR {float(block['iqr_ms']):.3f}) and cell {float(cell['median_ms']):.3f} ms "
                f"(IQR {float(cell['iqr_ms']):.3f}, imbalance {float(cell['imbalance']):.3f}). "
                f"{verdict}"
            )
    for op in OPS:
        cell = pick(rows, "weak", 884736, op, "cell", 8)
        block = pick(rows, "weak", 884736, op, "block", 8)
        bits.append(
            f"Weak scaling at eight ranks, {OP_LABEL[op]}: the block cut is "
            f"{float(block['median_ms']):.3f} ms (IQR {float(block['iqr_ms']):.3f}, "
            f"efficiency {float(block['efficiency']):.3f}, imbalance {float(block['imbalance']):.3f}) "
            f"and the cell-grid cut is {float(cell['median_ms']):.3f} ms "
            f"(IQR {float(cell['iqr_ms']):.3f}, efficiency {float(cell['efficiency']):.3f}, "
            f"imbalance {float(cell['imbalance']):.3f}). "
            f"{'The ranges overlap.' if ranges_overlap(block, cell) else 'The ranges do not overlap.'}"
        )
    PROSE_PATH.write_text("\n\n".join(bits) + "\n")


def write_figure(rows):
    import matplotlib

    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    fig, axes = plt.subplots(1, 2, figsize=(7.2, 3.4), constrained_layout=True)
    for ax, molecules, title in (
        (axes[0], 110592, "110592 molecules"),
        (axes[1], 884736, "884736 molecules"),
    ):
        ax.plot([1, 8], [1, 8], color="#000000", lw=0.8, ls="--", label="ideal")
        for op in OPS:
            for cut, style, marker in (("cell", "-", "o"), ("block", "--", "s")):
                series = [
                    row
                    for row in rows
                    if row["scale"] == "strong"
                    and int(row["molecules"]) == molecules
                    and row["operation"] == op
                    and row["cut"] == cut
                ]
                series.sort(key=lambda row: int(row["ranks"]))
                label = None
                if molecules == 110592:
                    label = f"{OP_LABEL[op].replace('$', '')} {cut}"
                ax.plot(
                    [int(row["ranks"]) for row in series],
                    [float(row["speedup"]) for row in series],
                    marker=marker,
                    ms=3.5,
                    ls=style,
                    color=COLOURS[op],
                    lw=1.2,
                    label=label,
                )
        ax.set_title(title)
        ax.set_xlabel("Ranks (1 thread each)")
        ax.set_ylabel("Speedup vs 1 rank")
        ax.set_xticks([1, 2, 4, 8])
        ax.set_ylim(0, 9)
        ax.spines["top"].set_visible(False)
        ax.spines["right"].set_visible(False)
    axes[0].legend(frameon=False, fontsize=7, ncol=2)
    fig.savefig(OUT / "partition-scale.pdf")
    fig.savefig(OUT / "partition-scale.svg")
    plt.close(fig)
    if not any(row["scale"] == "interface" for row in rows):
        return

    frames = (
        (110592, "ice/water\n110592"),
        (884736, "ice/water\n884736"),
        (105984, "hydrate\n105984"),
        (883200, "hydrate\n883200"),
    )
    fig, axes = plt.subplots(1, 4, figsize=(7.2, 2.6), sharey=False, constrained_layout=True)
    for ax, (molecules, title) in zip(axes, frames):
        for i, op in enumerate(OPS):
            block = pick(rows, "interface", molecules, op, "block", 8)
            cell = pick(rows, "interface", molecules, op, "cell", 8)
            x0 = i * 2
            def whisker(row):
                med = float(row["median_ms"])
                return [[med - float(row["q1_ms"])], [float(row["q3_ms"]) - med]]

            ax.bar(
                x0,
                float(block["median_ms"]),
                width=0.8,
                color=COLOURS[op],
                alpha=0.35,
                yerr=whisker(block),
                error_kw={"ecolor": "#333333", "lw": 0.6, "capsize": 2},
            )
            ax.bar(
                x0 + 0.8,
                float(cell["median_ms"]),
                width=0.8,
                color=COLOURS[op],
                yerr=whisker(cell),
                error_kw={"ecolor": "#333333", "lw": 0.6, "capsize": 2},
            )
        ax.set_title(title, fontsize=8)
        ax.set_xticks([i * 2 + 0.4 for i in range(len(OPS))])
        ax.set_xticklabels([OP_LABEL[op].replace("$", "") for op in OPS], fontsize=6, rotation=30)
        ax.spines["top"].set_visible(False)
        ax.spines["right"].set_visible(False)
    axes[0].set_ylabel("Eight-rank wall (ms)")
    fig.savefig(OUT / "partition-interface.pdf")
    fig.savefig(OUT / "partition-interface.svg")
    plt.close(fig)


def measure_v16():
    """884736-molecule crystal, pure MPI, and the same eight cores as hybrid teams."""
    DATA.mkdir(parents=True, exist_ok=True)
    (DATA / "partition-raw.txt").write_text("")
    rows = []
    molecules = 884736
    for ranks in (1, 2, 4, 8):
        for op in OPS:
            print(f"v16 {molecules} n={ranks} t=1 {op}", flush=True)
            got = run_bench(["48", "48", "48", str(REPS), op], ranks, f"v16 {op}", threads=1)
            add_pair(rows, "strong", molecules, ranks, op, got)
            dump(rows)
    for ranks, threads in ((4, 2), (2, 4)):
        for op in ("chill", "q6"):
            print(f"v16 hybrid n={ranks} t={threads} {op}", flush=True)
            got = run_bench(
                ["48", "48", "48", str(REPS), op],
                ranks,
                f"v16 hybrid {op}",
                threads=threads,
            )
            add_pair(rows, "hybrid", molecules, ranks, op, got)
            dump(rows)
    fill_efficiency(rows)
    dump(rows)
    return rows


def write_v16(rows):
    lines = [
        "\\begin{table}[tbp]",
        "\\centering",
        "\\caption{\\textbf{884736 molecules, eight cores.} "
        "Pure MPI is one thread on each rank. Hybrid splits the same cores into ranks and threads. "
        "Efficiency is the one-rank, one-thread time divided by the wall and by eight. "
        "The rule column is the cut the cost rule selects.}",
        "\\label{tab:partition-v16}",
        "\\begin{tabular}{llrrrrrr}",
        "\\toprule",
        "Layout & Operation & Block & IQR & Cell & IQR & Cell eff. & Rule \\\\",
        "\\midrule",
    ]
    layouts = [(8, 1, "strong", "8$\\times$1")]
    layouts += [(4, 2, "hybrid", "4$\\times$2"), (2, 4, "hybrid", "2$\\times$4")]
    for ranks, threads, scale, label in layouts:
        ops = OPS if scale == "strong" else ("chill", "q6")
        first = True
        for op in ops:
            block = pick(rows, scale, 884736, op, "block", ranks)
            cell = pick(rows, scale, 884736, op, "cell", ranks)
            if int(block.get("threads") or 1) != threads:
                continue
            shown = label if first else ""
            first = False
            lines.append(
                f"{shown} & {OP_LABEL[op]} & {float(block['median_ms']):.3f} & {float(block['iqr_ms']):.3f} "
                f"& {float(cell['median_ms']):.3f} & {float(cell['iqr_ms']):.3f} "
                f"& {float(cell['efficiency']):.3f} & {cell.get('rule') or ''} \\\\"
            )
    lines.extend(["\\bottomrule", "\\end{tabular}", "\\end{table}", ""])
    path = DATA / "tab-v16.tex"
    path.write_text("\n".join(lines))
    print(f"wrote {path}")


def measure_v17():
    """884736 crystal, hybrid layouts, the bond-trimmed rings, and the local frames."""
    DATA.mkdir(parents=True, exist_ok=True)
    (DATA / "partition-raw.txt").write_text("")
    rows = []
    reps = 8
    molecules = 884736
    for ranks in (1, 8):
        for op in OPS:
            print(f"v17 {molecules} n={ranks} t=1 {op}", flush=True)
            got = run_bench(
                ["48", "48", "48", str(reps), op], ranks, f"v17 {op}", threads=1, reps=reps
            )
            add_pair(rows, "strong", molecules, ranks, op, got)
            dump(rows)
    for ranks, threads in ((4, 2), (2, 4)):
        for op in ("chill", "q6"):
            print(f"v17 hybrid n={ranks} t={threads} {op}", flush=True)
            got = run_bench(
                ["48", "48", "48", str(reps), op],
                ranks,
                f"v17 hybrid {op}",
                threads=threads,
                reps=reps,
            )
            add_pair(rows, "hybrid", molecules, ranks, op, got)
            dump(rows)
    frames = (
        (ROOT / "input" / "traj" / "interface_f0.lammpstrj", 864, "interface"),
        (ROOT / "input" / "traj" / "nucleation_f1.lammpstrj", 4096, "nucleation"),
    )
    for path, count, name in frames:
        for ranks in (1, 8):
            for op in OPS:
                print(f"v17 {name} {count} n={ranks} {op}", flush=True)
                got = run_bench(
                    ["file", str(path), str(reps), op],
                    ranks,
                    f"v17 {name} {op}",
                    threads=1,
                    reps=reps,
                )
                add_pair(rows, name, count, ranks, op, got)
                dump(rows)
    fill_efficiency(rows)
    dump(rows)
    return rows


def write_v17(rows):
    lines = [
        "\\begin{table}[tbp]",
        "\\centering",
        "\\caption{\\textbf{884736 molecules on eight cores.} "
        "Pure MPI is one thread on each rank. A repeated ring cut reuses the cutoff rows and the trimmed bond ball. "
        "A moved frame rebuilds the rows of its own share. On eight thermal frames the distinct-frame ring efficiency is 0.794. "
        "CHILL+ and q6 on the same frames are 0.905 and 0.903. "
        "Hybrid uses the same eight cores and is slower: the eight-rank window is 117645 atoms, about 2.8 MiB, next to a 2 MiB L2, "
        "while the two-rank window is 451584 atoms and 340605 cells, about 11 MiB. "
        "The ring skin stays at floor(6/2) cutoffs. "
        "Efficiency divides the one-rank, one-thread time by the wall and by the workers. "
        "Eight repeats, type-7 median and IQR.}",
        "\\label{tab:partition-v17}",
        "\\begin{tabular}{llrrrrrr}",
        "\\toprule",
        "Layout & Operation & Block & IQR & Cell & IQR & Cell eff. & Rule \\\\",
        "\\midrule",
    ]
    layouts = [(8, 1, "strong", "8$\\times$1")]
    layouts += [(4, 2, "hybrid", "4$\\times$2"), (2, 4, "hybrid", "2$\\times$4")]
    for ranks, threads, scale, label in layouts:
        ops = OPS if scale == "strong" else ("chill", "q6")
        first = True
        for op in ops:
            block = pick(rows, scale, 884736, op, "block", ranks)
            cell = pick(rows, scale, 884736, op, "cell", ranks)
            if int(block.get("threads") or 1) != threads:
                continue
            shown = label if first else ""
            first = False
            lines.append(
                f"{shown} & {OP_LABEL[op]} & {float(block['median_ms']):.3f} & {float(block['iqr_ms']):.3f} "
                f"& {float(cell['median_ms']):.3f} & {float(cell['iqr_ms']):.3f} "
                f"& {float(cell['efficiency']):.3f} & {cell.get('rule') or ''} \\\\"
            )
    lines.extend(["\\bottomrule", "\\end{tabular}", "\\end{table}", ""])
    path = DATA / "tab-v17.tex"
    path.write_text("\n".join(lines))
    print(f"wrote {path}")


def main():
    for path in (BENCH, BUILD_SO):
        if not path.exists():
            raise SystemExit(f"missing {path}")
    if "--v16" in sys.argv:
        rows = measure_v16()
        write_v16(rows)
        return
    if "--v17" in sys.argv:
        rows = measure_v17()
        write_v17(rows)
        return
    if "--plot-only" in sys.argv:
        with CSV_PATH.open() as handle:
            rows = list(csv.DictReader(handle))
    else:
        rows = measure()
    write_table(rows)
    write_prose(rows)
    write_figure(rows)
    print(f"wrote {TABLE_PATH}")
    print(f"wrote {OUT / 'partition-scale.pdf'}")


if __name__ == "__main__":
    main()
