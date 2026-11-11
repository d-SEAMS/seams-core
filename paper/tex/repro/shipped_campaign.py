#!/usr/bin/env python3
"""Release campaign: rings, neighbours+rings, shipped cages, RSS, and the 1e6 frame.

Writes paper/tex/figures/data/shipped-campaign.json, shipped_macros.tex and
tab-shipped.tex. The paper reads those files and does not repeat the medians.
"""

from __future__ import annotations

import json
import os
import re
import resource
import subprocess
import sys
import time
from pathlib import Path

_HERE = Path(__file__).resolve()
ROOT = Path(os.environ["SEAMS_ROOT"]) if "SEAMS_ROOT" in os.environ else _HERE.parents[3]
sys.path.insert(0, str(ROOT / "tests"))
from bench_h2h_rings import (  # noqa: E402
    read_lammps_frame,
    run_ovito_worker,
    save_frame_npz,
    type7,
)

OUT = Path(os.environ.get("SEAMS_OUT", ROOT / "paper/tex/figures/data"))
H2H = ROOT / "tests/bench_h2h_rings.py"
BASE = Path(os.environ.get("SEAMS_BASE", "seams-base/bbdir-release"))
TIP = Path(os.environ.get("SEAMS_TIP", str(ROOT / "bbdir-release")))
CORE = os.environ.get("SEAMS_CORE", "0")
ORTH = Path(os.environ.get(
    "SEAMS_ORTH", ROOT / "input/traj/orth-1e5-ring-scale.lammpstrj"))
MILLION = Path(os.environ.get("SEAMS_1E6", "v37-frames/ice-1e6.lammpstrj"))
N_ORTH = 12
N_BIG = 6
STAGES = ("io_ms", "neigh_ms", "knn_ms", "knn_pair_ms", "rings_ms", "affil_ms", "seeded_ms")


def cell(values: list[float]) -> dict[str, float]:
    med = type7(values, 0.5)
    q1 = type7(values, 0.25)
    q3 = type7(values, 0.75)
    return {"n": len(values), "median": med, "q1": q1, "q3": q3, "iqr": q3 - q1}


def fmt(stat: dict[str, float], digits: int = 3) -> str:
    return (
        f"{stat['median']:.{digits}f} "
        f"({stat['q1']:.{digits}f}, {stat['q3']:.{digits}f}, {stat['iqr']:.{digits}f})"
    )


def host_name(raw: str) -> str:
    """Drop registered-mark spellings. lscpu prints Intel(R) Xeon(R)."""
    text = raw.replace("(R)", "").replace("(r)", "").replace("®", "")
    return " ".join(text.split())


LOGS = Path(os.environ.get("SEAMS_LOGS", ROOT / "paper/tex/repro/logs"))


def machine_banner(command: str) -> str:
    host = subprocess.check_output(["hostname", "-f"], text=True).strip()
    nproc = subprocess.check_output(["nproc"], text=True).strip()
    nproc_all = subprocess.check_output(["nproc", "--all"], text=True).strip()
    affinity = subprocess.check_output(
        ["taskset", "-pc", str(os.getpid())], text=True
    ).strip()
    lscpu = subprocess.check_output(["lscpu"], text=True).rstrip()
    return "\n".join([
        f"HOSTNAME {host}",
        f"nproc {nproc}",
        f"nproc_all {nproc_all}",
        f"affinity {affinity}",
        f"pin taskset -c {CORE}",
        f"COMMAND {command}",
        "lscpu <<LSC",
        lscpu,
        "LSC",
        "",
    ])


def write_log(name: str, command: str, body: str) -> None:
    LOGS.mkdir(parents=True, exist_ok=True)
    text = machine_banner(command) + body
    if "(R)" in text.split("lscpu <<LSC", 1)[-1].split("LSC", 1)[-1]:
        raise SystemExit(f"(R) leaked into {name}")
    (LOGS / name).write_text(text if text.endswith("\n") else text + "\n")


def seams_env(arm: str) -> dict[str, str]:
    env = os.environ.copy()
    env["OMP_NUM_THREADS"] = "1"
    env["OMP_DYNAMIC"] = "FALSE"
    env.pop("LD_LIBRARY_PATH", None)
    if arm == "base":
        return env
    # The tip binary's library sits beside the build.
    env["LD_LIBRARY_PATH"] = str(TIP / "src")
    # readcon stays on the tip runpath when LD_LIBRARY_PATH is only the lib dir
    # and the binary is the build-tree seams. Keep the binary in the build tree.
    return env


def run_stages(binary: Path, traj: Path) -> dict[str, float]:
    env = os.environ.copy()
    env["OMP_NUM_THREADS"] = "1"
    env["OMP_DYNAMIC"] = "FALSE"
    env.pop("LD_LIBRARY_PATH", None)
    out = subprocess.check_output(
        ["taskset", "-c", CORE, str(binary), str(traj), "1", "1", "1"],
        text=True,
        env=env,
    )
    row = {}
    for key in STAGES:
        match = re.search(rf"^{key} ([0-9.eE+-]+)$", out, re.M)
        if not match:
            raise RuntimeError(out)
        row[key] = float(match.group(1))
    row["neigh_rings_ms"] = row["neigh_ms"] + row["rings_ms"]
    return row


def run_cages(arm: str, mode: str, traj: Path, cap_base: bool) -> dict[str, float]:
    binary = (BASE if arm == "base" else TIP) / "seams"
    env = seams_env(arm)
    if arm == "tip":
        # Run the build-tree binary so its runpath finds readcon. The copied
        # library path is unnecessary when TIP is the build directory.
        env.pop("LD_LIBRARY_PATH", None)
    args = ["cages"] if mode == "seeded" else ["cages", "--graph", "cutoff"]
    cmd = ["/usr/bin/time", "-v", "taskset", "-c", CORE, str(binary), *args, str(traj)]

    def preexec() -> None:
        if cap_base and arm == "base":
            cap = 12 * 1024 * 1024 * 1024
            resource.setrlimit(resource.RLIMIT_AS, (cap, cap))

    t0 = time.perf_counter()
    proc = subprocess.run(
        cmd,
        text=True,
        env=env,
        stdout=subprocess.PIPE,
        stderr=subprocess.PIPE,
        preexec_fn=preexec if (cap_base and arm == "base") else None,
    )
    row = {
        "wall_ms": (time.perf_counter() - t0) * 1000.0,
        "rc": proc.returncode,
    }
    rss = re.search(r"Maximum resident set size \(kbytes\): (\d+)", proc.stderr)
    if rss:
        row["rss_kb"] = float(rss.group(1))
    tally = next((ln for ln in proc.stdout.splitlines() if ln.startswith("nop ")), "")
    row["tally"] = tally
    if proc.returncode != 0:
        text = (proc.stdout + "\n" + proc.stderr)
        row["oom"] = "std::bad_alloc" in text
    return row


def run_ovito(npz: Path, stage: str) -> dict[str, float]:
    out = run_ovito_worker(
        script_path=str(H2H),
        frame_npz=str(npz),
        cutoff=3.5,
        max_ring_size=6,
        ovito_threads=1,
        stage=stage,
    )
    return {"stage_ms": float(out["ms"]), "rings": float(out["count"])}


def sample(fn, n: int) -> list[dict]:
    fn()  # warmup
    return [fn() for _ in range(n)]


def tex_escape(text: str) -> str:
    return text


def summary_line(name: str, stat: dict | None, digits: int) -> str:
    if stat is None:
        return f"SUMMARY {name} n=0\n"
    return f"SUMMARY {name} n={stat['n']} {fmt(stat, digits)}\n"


def series_log(name: str, command: str, samples: list[str], summaries: list[str]) -> None:
    body = "".join(f"SAMPLE {line}\n" for line in samples)
    body += "".join(summaries)
    write_log(name, command, body)


def emit_logs(record, b_rows, t_rows, ov_rings, ov_bonds, cages, big_cages,
              ov_big_rings, ov_big_bonds) -> None:
    orth = str(ORTH)
    big = str(MILLION)
    series_log(
        "orth-stages-base.log",
        f"taskset -c {CORE} {BASE}/tests/bench_stages {orth} 1 1 1",
        [f"rings_ms {r['rings_ms']:.6f} neigh_rings_ms {r['neigh_rings_ms']:.6f}" for r in b_rows],
        [
            summary_line("ShipRingsBase", record["orth"]["rings_base"], 3),
            summary_line("ShipNeighBase", record["orth"]["neigh_rings_base"], 3),
        ],
    )
    series_log(
        "orth-stages-tip.log",
        f"taskset -c {CORE} {TIP}/tests/bench_stages {orth} 1 1 1",
        [f"rings_ms {r['rings_ms']:.6f} neigh_rings_ms {r['neigh_rings_ms']:.6f}" for r in t_rows],
        [
            summary_line("ShipRingsTip", record["orth"]["rings_tip"], 3),
            summary_line("ShipNeighTip", record["orth"]["neigh_rings_tip"], 3),
        ],
    )
    series_log(
        "orth-ovito-rings.log",
        "ovito FindRingsModifier in process, bonds built before the timer",
        [f"stage_ms {r['stage_ms']:.6f} rings {int(r['rings'])}" for r in ov_rings],
        [summary_line("ShipRingsOvito", record["orth"]["rings_ovito"], 3)],
    )
    series_log(
        "orth-ovito-bonds.log",
        "ovito CreateBondsModifier and FindRingsModifier inside one timer",
        [f"stage_ms {r['stage_ms']:.6f} rings {int(r['rings'])}" for r in ov_bonds],
        [summary_line("ShipNeighOvito", record["orth"]["neigh_rings_ovito"], 3)],
    )
    for arm, mode, med, rss in (
        ("base", "cutoff", "ShipCutBase", "ShipRssCutBase"),
        ("tip", "cutoff", "ShipCutTip", "ShipRssCutTip"),
        ("base", "seeded", "ShipSeedBase", "ShipRssSeedBase"),
        ("tip", "seeded", "ShipSeedTip", "ShipRssSeedTip"),
    ):
        rows = cages[f"{arm}_{mode}"]
        graph = "" if mode == "seeded" else " --graph cutoff"
        series_log(
            f"orth-cages-{arm}-{mode}.log",
            f"taskset -c {CORE} seams cages{graph} {orth}",
            [
                f"wall_ms {r['wall_ms']:.6f} rss_kb {r.get('rss_kb', float('nan')):.0f} {r.get('tally', '')}"
                for r in rows
            ],
            [
                summary_line(med, record["orth"][{"cutoff": "cutoff", "seeded": "seeded"}[mode] + f"_{arm}"], 3),
                summary_line(rss, record["orth"][f"rss_{mode}_{arm}"], 3),
            ],
        )
    for arm, mode in (("base", "cutoff"), ("base", "seeded"), ("tip", "cutoff"), ("tip", "seeded")):
        rows = big_cages[f"{arm}_{mode}"]
        graph = "" if mode == "seeded" else " --graph cutoff"
        extra = []
        if arm == "tip":
            key = "cutoff_tip" if mode == "cutoff" else "seeded_tip"
            extra.append(summary_line(
                "ShipBigCutTip" if mode == "cutoff" else "ShipBigSeedTip",
                record["ice1e6"][key], 3))
            if mode == "cutoff":
                extra.append(summary_line("ShipBigRssCutTip", record["ice1e6"]["rss_cutoff_tip"], 3))
        series_log(
            f"ice1e6-cages-{arm}-{mode}.log",
            f"taskset -c {CORE} seams cages{graph} {big}",
            [
                f"wall_ms {r['wall_ms']:.6f} rss_kb {r.get('rss_kb', float('nan')):.0f} "
                f"rc {r.get('rc', 0)} oom {int(bool(r.get('oom')))} {r.get('tally', '')}"
                for r in rows
            ],
            extra,
        )
    series_log(
        "ice1e6-ovito-rings.log",
        "ovito FindRingsModifier in process on the large frame",
        [f"stage_ms {r['stage_ms']:.6f}" for r in ov_big_rings],
        [summary_line("ShipBigRingsOvito", record["ice1e6"]["rings_ovito"], 3)],
    )
    series_log(
        "ice1e6-ovito-bonds.log",
        "ovito CreateBondsModifier and FindRingsModifier on the large frame",
        [f"stage_ms {r['stage_ms']:.6f}" for r in ov_big_bonds],
        [summary_line("ShipBigBondsOvito", record["ice1e6"]["bonds_ovito"], 3)],
    )
    write_log(
        "frame-sizes.log",
        "molecule counts and the dense neighbour matrix",
        "\n".join([
            f"NOP {record['orth_nop']}",
            f"BIG_NOP {record['big_nop']}",
            f"MATRIX_GIB {record['matrix_gib_1e6']:.3f}",
            f"BASE_1E6_OOM {int(bool(record['base_1e6_oom']))}",
            "HOST_MODEL Intel Xeon Processor",
            "",
        ]),
    )


def write_outputs(record: dict) -> None:
    record["host_model"] = host_name(record["host_model"])
    if "(R)" in record["host_model"] or record["host_model"] != "Intel Xeon Processor":
        raise SystemExit(f"host name {record['host_model']!r}")
    OUT.mkdir(parents=True, exist_ok=True)
    (OUT / "shipped-campaign.json").write_text(json.dumps(record, indent=2) + "\n")
    macros = ["% Generated by paper/tex/repro/shipped_campaign.py. Do not edit."]
    rows = []

    def add(name: str, stat: dict | None, unit: str) -> None:
        if stat is None:
            macros.append(f"\\newcommand{{\\{name}}}{{--}}")
            return
        macros.append(f"\\newcommand{{\\{name}}}{{{stat['median']:.3f}}}")
        macros.append(f"\\newcommand{{\\{name}Q}}{{{fmt(stat, 3)}}}")

    pairs = record["orth"]
    add("ShipRingsBase", pairs["rings_base"], "ms")
    add("ShipRingsTip", pairs["rings_tip"], "ms")
    add("ShipRingsOvito", pairs["rings_ovito"], "ms")
    add("ShipNeighBase", pairs["neigh_rings_base"], "ms")
    add("ShipNeighTip", pairs["neigh_rings_tip"], "ms")
    add("ShipNeighOvito", pairs["neigh_rings_ovito"], "ms")
    add("ShipCutBase", pairs["cutoff_base"], "ms")
    add("ShipCutTip", pairs["cutoff_tip"], "ms")
    add("ShipSeedBase", pairs["seeded_base"], "ms")
    add("ShipSeedTip", pairs["seeded_tip"], "ms")
    add("ShipRssCutBase", pairs["rss_cutoff_base"], "kB")
    add("ShipRssCutTip", pairs["rss_cutoff_tip"], "kB")
    add("ShipRssSeedBase", pairs["rss_seeded_base"], "kB")
    add("ShipRssSeedTip", pairs["rss_seeded_tip"], "kB")
    big = record["ice1e6"]
    add("ShipBigCutTip", big["cutoff_tip"], "ms")
    add("ShipBigSeedTip", big["seeded_tip"], "ms")
    add("ShipBigRingsOvito", big["rings_ovito"], "ms")
    add("ShipBigBondsOvito", big["bonds_ovito"], "ms")
    add("ShipBigRssCutTip", big["rss_cutoff_tip"], "kB")
    macros.append(f"\\newcommand{{\\ShipHost}}{{{record['host_model']}}}")
    macros.append(f"\\newcommand{{\\ShipNop}}{{{record['orth_nop']}}}")
    macros.append(f"\\newcommand{{\\ShipBigNop}}{{{record['big_nop']}}}")
    macros.append(
        f"\\newcommand{{\\ShipMatrixGib}}{{{record['matrix_gib_1e6']:.3f}}}"
    )
    (OUT / "shipped_macros.tex").write_text("\n".join(macros) + "\n")

    def row(label: str, *stats: dict | None) -> str:
        cells = []
        for stat in stats:
            cells.append("--" if stat is None else fmt(stat))
        return label + " & " + " & ".join(cells) + r" \\"

    lines = [
        row("Rings", pairs["rings_base"], pairs["rings_tip"], pairs["rings_ovito"]),
        row("Neighbours and rings", pairs["neigh_rings_base"], pairs["neigh_rings_tip"], pairs["neigh_rings_ovito"]),
        row(r"\texttt{seams cages --graph cutoff}", pairs["cutoff_base"], pairs["cutoff_tip"], None),
        row(r"\texttt{seams cages}", pairs["seeded_base"], pairs["seeded_tip"], None),
        row("Cutoff resident set (kB)", pairs["rss_cutoff_base"], pairs["rss_cutoff_tip"], None),
        row("Seeded resident set (kB)", pairs["rss_seeded_base"], pairs["rss_seeded_tip"], None),
    ]
    (OUT / "tab-shipped.tex").write_text("\n".join(lines) + "\n")
    print("wrote", OUT / "shipped-campaign.json")


def main() -> None:
    host = subprocess.check_output(["lscpu"], text=True)
    model = next(ln.split(":", 1)[1].strip() for ln in host.splitlines() if ln.startswith("Model name"))
    orth = read_lammps_frame(str(ORTH), 1)
    npz = Path("v40-orth.npz")
    save_frame_npz(orth, str(npz))
    big = read_lammps_frame(str(MILLION), 1)
    npz_big = Path("v40-1e6.npz")
    save_frame_npz(big, str(npz_big))

    base_stages = Path(BASE / "tests/bench_stages")
    tip_stages = Path(TIP / "tests/bench_stages")
    print("sampling stages", flush=True)
    b_rows = sample(lambda: run_stages(base_stages, ORTH), N_ORTH)
    t_rows = sample(lambda: run_stages(tip_stages, ORTH), N_ORTH)
    print("sampling ovito orth", flush=True)
    ov_rings = sample(lambda: run_ovito(npz, "rings"), N_ORTH)
    ov_bonds = sample(lambda: run_ovito(npz, "bonds"), N_ORTH)
    print("sampling cages orth", flush=True)
    cages = {}
    for arm in ("base", "tip"):
        for mode in ("cutoff", "seeded"):
            print(arm, mode, flush=True)
            cages[f"{arm}_{mode}"] = sample(lambda arm=arm, mode=mode: run_cages(arm, mode, ORTH, False), N_ORTH)
    print("sampling 1e6", flush=True)
    big_cages = {}
    for arm in ("base", "tip"):
        for mode in ("cutoff", "seeded"):
            print("1e6", arm, mode, flush=True)
            big_cages[f"{arm}_{mode}"] = sample(
                lambda arm=arm, mode=mode: run_cages(arm, mode, MILLION, True), N_BIG
            )
    ov_big_rings = sample(lambda: run_ovito(npz_big, "rings"), N_BIG)
    ov_big_bonds = sample(lambda: run_ovito(npz_big, "bonds"), N_BIG)

    def ok(rows: list[dict], key: str) -> dict | None:
        vals = [row[key] for row in rows if row.get("rc", 0) == 0 and key in row]
        if len(vals) < 6:
            return None
        return cell(vals)

    orth_rec = {
        "rings_base": cell([r["rings_ms"] for r in b_rows]),
        "rings_tip": cell([r["rings_ms"] for r in t_rows]),
        "rings_ovito": cell([r["stage_ms"] for r in ov_rings]),
        "neigh_rings_base": cell([r["neigh_rings_ms"] for r in b_rows]),
        "neigh_rings_tip": cell([r["neigh_rings_ms"] for r in t_rows]),
        "neigh_rings_ovito": cell([r["stage_ms"] for r in ov_bonds]),
        "cutoff_base": ok(cages["base_cutoff"], "wall_ms"),
        "cutoff_tip": ok(cages["tip_cutoff"], "wall_ms"),
        "seeded_base": ok(cages["base_seeded"], "wall_ms"),
        "seeded_tip": ok(cages["tip_seeded"], "wall_ms"),
        "rss_cutoff_base": ok(cages["base_cutoff"], "rss_kb"),
        "rss_cutoff_tip": ok(cages["tip_cutoff"], "rss_kb"),
        "rss_seeded_base": ok(cages["base_seeded"], "rss_kb"),
        "rss_seeded_tip": ok(cages["tip_seeded"], "rss_kb"),
    }
    base_oom = all(row.get("oom") for row in big_cages["base_cutoff"] + big_cages["base_seeded"])
    record = {
        "host_model": model,
        "orth_nop": int(orth.positions.shape[0]),
        "big_nop": int(big.positions.shape[0]),
        "base_1e6_oom": base_oom,
        "matrix_gib_110592": (110592 * ((110592 + 63) >> 6) * 8) / 1024**3,
        "matrix_gib_1e6": None,
        "orth": orth_rec,
        "ice1e6": {
            "cutoff_tip": ok(big_cages["tip_cutoff"], "wall_ms"),
            "seeded_tip": ok(big_cages["tip_seeded"], "wall_ms"),
            "rss_cutoff_tip": ok(big_cages["tip_cutoff"], "rss_kb"),
            "rings_ovito": cell([r["stage_ms"] for r in ov_big_rings]),
            "bonds_ovito": cell([r["stage_ms"] for r in ov_big_bonds]),
            "base_cutoff_abort": [row["wall_ms"] for row in big_cages["base_cutoff"]],
            "base_seeded_abort": [row["wall_ms"] for row in big_cages["base_seeded"]],
        },
    }
    nbig = record["big_nop"]
    record["matrix_gib_1e6"] = (nbig * ((nbig + 63) >> 6) * 8) / 1024**3
    emit_logs(record, b_rows, t_rows, ov_rings, ov_bonds, cages, big_cages,
              ov_big_rings, ov_big_bonds)
    write_outputs(record)
    for key, stat in orth_rec.items():
        if stat:
            print(f"{key} {fmt(stat)}")
    print("base_1e6_oom", base_oom, "matrix_GiB", record["matrix_gib_1e6"])


if __name__ == "__main__":
    main()
