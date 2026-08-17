#!/usr/bin/env python3
"""Assemble the reproducibility-campaign outputs into one manifest.

Usage: aggregate.py <results_dir>

Reads the text outputs the Snakemake DAG produced and emits the manifest
JSON on stdout. Missing files become nulls rather than errors so a partial
run still yields an inspectable object; the DAG itself decides completeness.
"""
import json
import pathlib
import re
import sys


def read(path):
    try:
        return path.read_text()
    except OSError:
        return ""


def parse_scaling(text):
    rows = {}
    for line in text.splitlines():
        m = re.match(
            r"\s*(\d+)\s+([\d.]+)\s+([\d.]+)\s+([\d.]+)\s+([\d.]+|skipped)",
            line,
        )
        if m:
            rows[int(m.group(1))] = {
                "neighListO_ms": float(m.group(2)),
                "byIndex_ms": float(m.group(3)),
                "getCorrelPlus_ms": float(m.group(4)),
                "ringNetwork_ms": (
                    None if m.group(5) == "skipped" else float(m.group(5))
                ),
            }
    return rows


def parse_kv(text):
    out = {}
    for line in text.splitlines():
        m = re.match(r"([A-Za-z/_+-]+)\s+([-\d.]+)\s*$", line.strip())
        if m:
            out[m.group(1)] = float(m.group(2))
    return out


def parse_pipeline(text):
    rows = {}
    for line in text.splitlines():
        m = re.match(
            r"\s*(\d+)\s+([\d.]+)\s+([\d.]+)\s+([\d.]+)\s+([\d.]+)\s+([\d.]+)",
            line,
        )
        if m:
            rows[int(m.group(1))] = {
                "neigh_ms": float(m.group(2)),
                "knn_ms": float(m.group(3)),
                "knn_pair_ms": float(m.group(4)),
                "rings_ms": float(m.group(5)),
                "affil_ms": float(m.group(6)),
            }
    return rows


def parse_incremental_file(text):
    out = {}
    for line in text.splitlines():
        m = re.match(r"(\S+)\s+(.*)$", line.strip())
        if not m:
            continue
        key, val = m.group(1), m.group(2)
        if key in ("nAtoms", "frames", "sites", "recomp", "balls"):
            out[key] = val
        elif key == "identical":
            out[key] = val
        else:
            try:
                out[key] = float(val)
            except ValueError:
                out[key] = val
    return out


def parse_overhead(text):
    rows = {}
    for line in text.splitlines():
        m = re.match(r"\s*(\d+)\s+([\d.]+)\s+([\d.]+)\s+([\d.]+)", line)
        if m:
            rows[int(m.group(1))] = {
                "vesin_ms": float(m.group(2)),
                "neighListO_ms": float(m.group(3)),
            }
    return rows


def parse_strong(text):
    for line in text.splitlines():
        m = re.match(
            r"\s*(\d+)\s+\d+\s+\d+\s+\d+\s+([\d.]+)\s+([\d.]+)"
            r"(?:\s+([\d.]+)\s+([\d.]+))?",
            line,
        )
        if m:
            row = {"neigh_ms": float(m.group(2)), "ql_ms": float(m.group(3))}
            if m.group(4) is not None:
                row["index_ms"] = float(m.group(4))
                row["rings_ms"] = float(m.group(5))
            return row
    return None


def parse_trajectory(text):
    frames = []
    summary = None
    for line in text.splitlines():
        if line.startswith("# steady"):
            m = re.search(
                r"ringFull_ms ([\d.]+) ringInc_ms ([\d.]+) "
                r"affilBatch_ms ([\d.]+) affilInc_ms ([\d.]+) allEqual (\w+)",
                line,
            )
            if m:
                summary = {
                    "ringFull_ms": float(m.group(1)),
                    "ringInc_ms": float(m.group(2)),
                    "affilBatch_ms": float(m.group(3)),
                    "affilInc_ms": float(m.group(4)),
                    "allEqual": m.group(5) == "yes",
                }
            continue
        if line.startswith("#"):
            continue
        parts = line.split()
        if len(parts) == 11:
            frames.append(
                {
                    "frame": int(parts[0]),
                    "nRings": int(parts[2]),
                    "ringFull_ms": float(parts[3]),
                    "ringInc_ms": float(parts[4]),
                    "ringRecomputed": int(parts[5]),
                    "affilBatch_ms": float(parts[6]),
                    "affilInc_ms": float(parts[7]),
                    "affilReclassified": int(parts[8]),
                    "ringsEqual": parts[9] == "yes",
                    "affilEqual": parts[10] == "yes",
                }
            )
    return {"frames": frames, "steady": summary}


def parse_ql(text):
    tools = {}
    current = None
    for line in text.splitlines():
        m = re.match(r"tool=(\w+) lattice=(\w+)", line)
        if m:
            current = tools.setdefault(m.group(1), {}).setdefault(
                m.group(2), {}
            )
            continue
        if current is None:
            continue
        for key, value in re.findall(r"(ql\d|vor_q\d)=([\d.]+)", line):
            current[key] = float(value)
    return tools


def load_json(path):
    try:
        return json.loads(path.read_text())
    except (OSError, ValueError):
        return None


def digest_figshare_demos(obj):
    if not obj:
        return None
    return {
        key: {k: v for k, v in entry.items() if k != "outputs"}
        for key, entry in obj.items()
    }


def digest_figshare_incremental(obj):
    if not obj:
        return None
    rows = obj.get("per_frame", [])
    return {
        "trajectory": obj.get("trajectory"),
        "doi": obj.get("doi"),
        "frames": obj.get("frames"),
        "batch_ring_seconds": obj.get("batch_ring_seconds"),
        "incremental_ring_seconds": obj.get("incremental_ring_seconds"),
        "first_frame": rows[0] if rows else None,
        "last_frame": rows[-1] if rows else None,
        "peak_seeded_ddc_atoms": max(
            (r["seeded_ddc_atoms"] for r in rows), default=None
        ),
    }


def main():
    out_dir = pathlib.Path(sys.argv[1])
    conditions = {}
    for line in read(out_dir / "conditions.txt").splitlines():
        if ":" in line:
            key, _, value = line.partition(":")
            conditions[key.strip()] = value.strip()

    strong = {}
    for f in sorted(out_dir.glob("tip-strong-t*.txt")):
        threads = int(re.search(r"t(\d+)", f.name).group(1))
        strong[threads] = parse_strong(read(f))

    ql = parse_ql(read(out_dir / "ql-python.txt"))
    ql.update(parse_ql(read(out_dir / "ql-dseams.txt")))

    manifest = {
        "conditions": conditions,
        "identity_gate": "Fail:              0"
        in read(out_dir / "tip-test.log"),
        "tip": {
            "scaling": parse_scaling(read(out_dir / "tip-scaling.txt")),
            "cages": parse_kv(read(out_dir / "tip-cages.txt")),
            "overhead": parse_overhead(read(out_dir / "tip-overhead.txt")),
            "strong": strong,
        },
        "base": {
            "scaling": parse_scaling(read(out_dir / "base-scaling.txt")),
            "cages": parse_kv(read(out_dir / "base-cages.txt")),
        },
        "pipeline": parse_pipeline(read(out_dir / "tip-pipeline.txt")),
        "incremental": {
            int(p.stem.split("-")[-1]): parse_incremental_file(read(p))
            for p in sorted(out_dir.glob("tip-incremental-*.txt"))
        },
        "stages": {
            "cubic": parse_kv(read(out_dir / "tip-stages-cubic.txt")),
            "nucleation": parse_kv(read(out_dir / "tip-stages-nucleation.txt")),
        },
        "trajectory_incremental": parse_trajectory(
            read(out_dir / "trajectory-incremental.txt")
        ),
        "ql_compare": ql,
        "figshare": {
            "demos": digest_figshare_demos(
                load_json(
                    out_dir / "figshare-demos" / "figshare-demos.json"
                )
            ),
            "incremental": digest_figshare_incremental(
                load_json(out_dir / "figshare-incremental.json")
            ),
        },
    }
    json.dump(manifest, sys.stdout, indent=2, sort_keys=True)
    print()


if __name__ == "__main__":
    main()
