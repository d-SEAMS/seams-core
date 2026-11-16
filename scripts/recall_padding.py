#!/usr/bin/env python3
"""False positives, interface layers, and classify cost of open-face padding.

Bulk mW liquid is a 4096-molecule hold at 230 K and at 300 K. The frame
written as the thermostat switches is left out, leaving ten frames at the
hold temperature. Brine frames are the oxygen rows of the pinned
d-SEAMS/brine-front-230K set. The interface lane is the oxygen dump the
TIP4P/Ice basal rule already writes. CHILL+ is the published rule.

A false positive is a non-negative shell label on an atom CHILL+ does not
call cubic or hexagonal. Timing is one process: eight replicates, padding
off then on, at 4096 and at the 3x3x3 tile of that liquid (110592).
"""

from __future__ import annotations

import argparse
import csv
import os
import subprocess
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
HEADER = [
    "record",
    "lane",
    "temperature_K",
    "frame",
    "layer",
    "z_lo_A",
    "molecules",
    "chill_ice",
    "dseams_off",
    "dseams_on",
    "fp_off",
    "fp_on",
    "kind_changed",
    "agree_off",
    "agree_on",
    "padding",
    "rep",
    "seconds",
    "cores_online",
    "cores_affinity",
    "cores_used",
    "source",
]


def child_env(binary: Path, threads: str) -> dict[str, str]:
    env = os.environ.copy()
    env["OMP_NUM_THREADS"] = threads
    env["OMP_DYNAMIC"] = "FALSE"
    origin = binary.resolve().parent
    libs = [
        origin.parent / "src",
        origin.parent / "chemfiles-prefix" / "lib",
        origin.parent / "subprojects" / "readcon-core",
        origin.parent / "subprojects" / "readcon-db",
        Path("/lib/x86_64-linux-gnu"),
    ]
    extra = env.get("SEAMS_RUNTIME_LIBDIR", "")
    if extra:
        libs.insert(0, Path(extra))
    parts = [str(path) for path in libs if path.is_dir()]
    old = env.get("LD_LIBRARY_PATH", "")
    if old:
        parts.append(old)
    if parts:
        env["LD_LIBRARY_PATH"] = ":".join(parts)
    return env


def run_probe(binary: Path, args: list[str], threads: str) -> list[str]:
    done = subprocess.run(
        [str(binary), *args],
        cwd=ROOT,
        env=child_env(binary, threads),
        text=True,
        capture_output=True,
    )
    if done.returncode != 0:
        sys.stderr.write(done.stderr)
        raise SystemExit(f"padding probe exited {done.returncode}")
    return [line for line in done.stdout.splitlines() if line.strip()]


def blank() -> dict[str, str]:
    return {name: "" for name in HEADER}


def liquid_row(fields: list[str], source: str) -> dict[str, str]:
    row = blank()
    row.update(
        {
            "record": "liquid",
            "lane": fields[1],
            "temperature_K": fields[2],
            "frame": fields[3],
            "layer": fields[4],
            "z_lo_A": fields[5],
            "molecules": fields[6],
            "chill_ice": fields[7],
            "dseams_off": fields[8],
            "dseams_on": fields[9],
            "fp_off": fields[10],
            "fp_on": fields[11],
            "kind_changed": fields[12],
            "agree_off": fields[13],
            "agree_on": fields[14],
            "source": source,
        }
    )
    return row


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--probe", type=Path, required=True)
    parser.add_argument("--out", type=Path, required=True)
    parser.add_argument("--reps", type=int, default=8)
    parser.add_argument("--ctest", default="not-run")
    parser.add_argument("--pytest", default="not-run")
    args = parser.parse_args()
    if args.reps < 8:
        raise SystemExit("timing needs at least 8 replicates")
    if not args.probe.is_file():
        raise SystemExit(f"missing {args.probe}")

    liquid = {
        "230": ROOT / "data" / "science-frames" / "padding-liquid-230" / "mw-4096" / "quench.lammpstrj",
        "300": ROOT / "data" / "science-frames" / "padding-liquid-300" / "mw-4096" / "quench.lammpstrj",
    }
    interface = ROOT / "data" / "science-frames" / "tip4p-basal-230K.lammpstrj"
    brine = ROOT / "data" / "science-frames" / "brine-hf"
    for path in (*liquid.values(), interface):
        if not path.is_file():
            raise SystemExit(f"missing {path.name}")
    xyz_files = sorted(brine.glob("T230_*.xyz"))
    if len(xyz_files) < 8:
        raise SystemExit("brine frame set is short")

    rows: list[dict[str, str]] = []
    liquid_source = (
        "mW bulk liquid, 4096 molecules, 2000-step 350 K melt, "
        "100 ps hold, seed 20261114, ten frames after the thermostat switch"
    )
    for temperature, path in liquid.items():
        lines = run_probe(args.probe, ["liquid", str(path), temperature, "2", "11"], "1")
        if len(lines) != 10:
            raise SystemExit(f"{temperature} K returned {len(lines)} frames")
        for line in lines:
            fields = line.split(",")
            if fields[0] != "liquid" or fields[6] != "4096":
                raise SystemExit(f"unexpected liquid row {line}")
            rows.append(liquid_row(fields, liquid_source))

    interface_source = "TIP4P/Ice basal oxygen dump at 230 K, ten frames"
    layer_lines = run_probe(
        args.probe,
        ["layers", str(interface), "tip4p-basal", "230", "8", "80", "8", "2"],
        "1",
    )
    for line in layer_lines:
        fields = line.split(",")
        row = liquid_row(fields, interface_source)
        row["record"] = fields[0]
        row["lane"] = fields[1]
        rows.append(row)

    for path in xyz_files:
        tag = path.stem
        lines = run_probe(args.probe, ["xyz", str(path), tag, "230", "2"], "1")
        for line in lines:
            fields = line.split(",")
            row = liquid_row(fields, "d-SEAMS/brine-front-230K oxygen rows")
            row["record"] = fields[0]
            row["lane"] = fields[1]
            rows.append(row)

    team = str(len(os.sched_getaffinity(0)))
    timed = run_probe(
        args.probe,
        ["time", str(liquid["300"]), "2", str(args.reps)],
        team,
    )
    cores = [line.split(",") for line in timed if line.startswith("cores,")]
    if len(cores) != 1 or len(cores[0]) != 4:
        raise SystemExit("timing did not record a core count")
    online, affinity, used = cores[0][1:]
    if used == "":
        used = affinity
    samples = [line.split(",") for line in timed if line.startswith("time,")]
    if len(samples) != args.reps * 4:
        raise SystemExit(f"timing returned {len(samples)} samples")
    for fields in samples:
        molecules = fields[1]
        source = (
            "300 K bulk-liquid frame"
            if molecules == "4096"
            else "3x3x3 tile of the 300 K bulk-liquid frame"
        )
        row = blank()
        row.update(
            {
                "record": "time",
                "lane": "classify",
                "molecules": molecules,
                "padding": fields[2],
                "rep": fields[3],
                "seconds": fields[4],
                "cores_online": online,
                "cores_affinity": affinity,
                "cores_used": used,
                "source": source,
            }
        )
        rows.append(row)

    def summarize(record: str, lane: str, temperature: str, source: str) -> None:
        chosen = [
            row
            for row in rows
            if row["record"] == record and row["lane"] == lane and row["temperature_K"] == temperature
        ]
        if not chosen:
            return
        total = blank()
        total["record"] = "summary"
        total["lane"] = lane
        total["temperature_K"] = temperature
        total["frame"] = str(len(chosen))
        total["source"] = source
        for key in (
            "molecules",
            "chill_ice",
            "dseams_off",
            "dseams_on",
            "fp_off",
            "fp_on",
            "kind_changed",
            "agree_off",
            "agree_on",
        ):
            total[key] = str(sum(int(row[key]) for row in chosen))
        rows.append(total)

    for temperature in ("230", "300"):
        summarize("liquid", "mw-bulk", temperature, "sum of the ten hold frames")
    summarize("frame", "tip4p-basal", "230", "sum of the ten basal frames")
    for lane in sorted({row["lane"] for row in rows if row["record"] == "frame" and row["lane"] != "tip4p-basal"}):
        summarize("frame", lane, "230", "one pinned brine frame")

    fp_off = sum(int(row["fp_off"]) for row in rows if row["record"] == "liquid")
    fp_on = sum(int(row["fp_on"]) for row in rows if row["record"] == "liquid")
    gate = blank()
    gate["record"] = "gate"
    gate["lane"] = "mw-bulk"
    gate["fp_off"] = str(fp_off)
    gate["fp_on"] = str(fp_on)
    if fp_on > fp_off:
        gate["source"] = (
            f"padding adds {fp_on - fp_off} false positives on the bulk liquid; "
            "open-face padding is opt-in"
        )
    else:
        gate["source"] = (
            "padding adds no false positive on the bulk liquid; "
            "the open-face completion stays on"
        )
    rows.append(gate)

    for name, status in (("ctest", args.ctest), ("pytest", args.pytest)):
        suite = blank()
        suite["record"] = "suite"
        suite["lane"] = name
        suite["source"] = status
        rows.append(suite)

    args.out.parent.mkdir(parents=True, exist_ok=True)
    with args.out.open("w", newline="") as handle:
        writer = csv.DictWriter(handle, HEADER, lineterminator="\n")
        writer.writeheader()
        writer.writerows(rows)
    print(f"wrote {args.out}")
    print(gate["source"])


if __name__ == "__main__":
    main()
