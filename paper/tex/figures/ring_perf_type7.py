#!/usr/bin/env python3
"""Type 7 medians of the one-thread ring samples and the neighbour samples.

The samples are the committed JSON beside this script. The stored
median in that JSON is not the value printed here. The ring record
must set ovito_threads to 1.
"""

import json
import math
from pathlib import Path


def quantile_type7(values, p):
    xs = sorted(float(v) for v in values)
    n = len(xs)
    if n < 6:
        raise ValueError(f"need n >= 6, got {n}")
    h = (n - 1) * p + 1
    j = int(math.floor(h))
    g = h - j
    xj = xs[j - 1]
    if j >= n:
        return xj
    return xj + g * (xs[j] - xj)


def fmt(x):
    return f"{x:.3f}"


def cell(samples):
    med = quantile_type7(samples, 0.5)
    q1 = quantile_type7(samples, 0.25)
    q3 = quantile_type7(samples, 0.75)
    return med, q1, q3


def tex_cell(samples):
    med, q1, q3 = cell(samples)
    return rf"\num{{{fmt(med)}}} (\iqr{{{fmt(q1)}}}{{{fmt(q3)}}})"


def main():
    root = Path(__file__).resolve().parent / "data" / "ring-perf"
    ring = json.loads((root / "ring-head2head.json").read_text())
    neigh = json.loads((root / "neigh-head2head.json").read_text())
    ovito_threads = int(ring["thread_policy"]["ovito_threads"])
    dseams_threads = [int(value) for value in ring["thread_policy"]["dseams_threads_reported"]]
    if ovito_threads != 1 or any(value != 1 for value in dseams_threads):
        raise SystemExit(
            f"ring record is not one thread: ovito {ovito_threads} dseams {dseams_threads}"
        )
    machine = ring["machine"]
    lines = []
    lines.append(
        "hardware uname -m "
        + machine["uname_m"]
        + " nproc "
        + str(machine["nproc"])
        + " CPU "
        + machine["cpu_model"]
    )
    lines.append(
        "protocol repeats "
        + str(ring["protocol"]["repeats"])
        + " median type 7 quartile type 7"
    )
    lines.append("git_sha " + ring["git_sha"])
    lines.append(
        "ring ovito_threads "
        + str(ovito_threads)
        + " dseams_threads "
        + " ".join(str(value) for value in dseams_threads)
    )
    lines.append(
        "neigh ovito_threads " + str(neigh["thread_policy"]["ovito_threads"])
    )

    ring_rows = []
    for name, body in ring["datasets"].items():
        stage = body["ring_stage_ms"]
        dse = tex_cell(stage["dseams_omp1"]["samples_ms"])
        ovi = tex_cell(stage["ovito"]["samples_ms"])
        med, q1, q3 = cell(stage["dseams_omp1"]["samples_ms"])
        omed, oq1, oq3 = cell(stage["ovito"]["samples_ms"])
        lines.append(
            f"ring {name} n {body['n_atoms']} dseams_omp1 "
            f"median {fmt(med)} q1 {fmt(q1)} q3 {fmt(q3)} "
            f"ovito median {fmt(omed)} q1 {fmt(oq1)} q3 {fmt(oq3)}"
        )
        ring_rows.append(
            f"{name} & {body['n_atoms']} & {dse} & {ovi} \\\\"
        )
        lines.append(
            f"rings {name} total {body['rings_total']} size6 {body['rings_size6']}"
        )

    tools = (
        ("dseams_omp1", "d-SEAMS"),
        ("freud_aabb_cutoff", "freud AABB"),
        ("freud_linkcell_cutoff", "freud link cell"),
        ("ovito_cutoff", "OVITO"),
        ("scipy_ckdtree_cutoff", "scipy"),
    )
    projected = neigh["linkcell_policy"]["projections"]
    neigh_rows = []
    for name, body in neigh["datasets"].items():
        stage = body["cutoff_stage"]
        bits = [name, str(body["n_atoms"])]
        for key, _label in tools:
            if key in stage and "samples_ms" in stage[key]:
                med, q1, q3 = cell(stage[key]["samples_ms"])
                lines.append(
                    f"cutoff {name} {key} median {fmt(med)} "
                    f"q1 {fmt(q1)} q3 {fmt(q3)}"
                )
                bits.append(tex_cell(stage[key]["samples_ms"]))
            else:
                proj_key = f"{name}/cutoff"
                if key == "freud_linkcell_cutoff" and proj_key in projected:
                    pred = projected[proj_key]["projected_sample_ms"]
                    lines.append(
                        f"cutoff {name} {key} projected {fmt(pred)}"
                    )
                    bits.append(
                        rf"\num{{{fmt(pred)}}}\textsuperscript{{p}}"
                    )
                else:
                    lines.append(f"cutoff {name} {key} not measured")
                    bits.append("")
        # bits layout is awkward because tex_cell already has two columns.
        # Rebuild the row from parts that are either 'num & iqr' or ' & '.
        neigh_rows.append(bits)

    out_dir = Path(__file__).resolve().parent / "data"
    ring_tex = out_dir / "tab-ring-perf.tex"
    ring_body = ["% snapper:off"]
    for row in ring_rows:
        ring_body.append(row)
    ring_body.append("% snapper:on")
    ring_tex.write_text("\n".join(ring_body) + "\n")

    # Flatten neighbour rows. Each tool contributes two column fields
    # already joined by & inside tex_cell, or a single blank pair.
    neigh_tex = out_dir / "tab-neigh-perf.tex"
    body = ["% snapper:off"]
    for bits in neigh_rows:
        # bits[0] name, bits[1] n, then one entry per tool
        fields = [bits[0], bits[1]]
        for entry in bits[2:]:
            fields.append(entry)
        body.append(" & ".join(fields) + r" \\")
    body.append("% snapper:on")
    neigh_tex.write_text("\n".join(body) + "\n")

    log = out_dir / "ring-perf-type7.log"
    log.write_text("\n".join(lines) + "\n")
    print("\n".join(lines))


if __name__ == "__main__":
    main()
