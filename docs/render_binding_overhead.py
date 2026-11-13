#!/usr/bin/env python3
# Copyright (c) d-SEAMS developers
"""Draw the binding-overhead page from benchmarks/binding-overhead.csv.

The PNG is written at docs build time and is not stored in git.
"""

from __future__ import annotations

import csv
import sys
from collections import defaultdict
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
CSV = ROOT / "benchmarks" / "binding-overhead.csv"
PNG = ROOT / "docs" / "source" / "_static" / "binding-overhead.png"
MD = ROOT / "docs" / "explanation" / "binding-overhead.md"
RST = ROOT / "docs" / "source" / "explanation" / "binding-overhead.rst"
ORG = ROOT / "docs" / "orgmode" / "explanation" / "binding-overhead.org"

ORDER = ("c", "cxx", "rust", "fortran", "python", "julia", "r", "lua", "freud", "ovito")
OPS = ("read", "knn", "steinhardt", "chill", "rings")
OP_NAME = {
    "read": "read",
    "knn": "knn",
    "steinhardt": "q6",
    "chill": "CHILL",
    "rings": "rings",
}


def load():
    rows = []
    with CSV.open() as handle:
        for row in csv.DictReader(handle):
            row["molecules"] = int(row["molecules"])
            row["threads"] = int(row["threads"])
            row["median"] = float(row["median_ms"])
            row["iqr"] = float(row["iqr_ms"])
            rows.append(row)
    return rows


def lookup(rows, binding, op, molecules, threads):
    for row in rows:
        if (row["binding"], row["operation"], row["molecules"], row["threads"]) == (
            binding, op, molecules, threads,
        ):
            return row
    return None


def ratio(rows, binding, op, molecules, threads):
    got = lookup(rows, binding, op, molecules, threads)
    base = lookup(rows, "c", op, molecules, threads)
    if got is None or base is None or base["median"] <= 0:
        return None
    return got["median"] / base["median"]


def cell(row):
    if row is None:
        return ""
    return "{:.3f} ({:.3f})".format(row["median"], row["iqr"])


def table_lines(rows):
    threads = sorted({row["threads"] for row in rows})
    sizes = sorted({row["molecules"] for row in rows})
    lines = []
    for op in OPS:
        lines.append("")
        lines.append("### {}".format(OP_NAME[op]))
        lines.append("")
        header = ["binding"]
        for molecules in sizes:
            for width in threads:
                header.append("{} mol, {} thr".format(molecules, width))
        lines.append("| " + " | ".join(header) + " |")
        lines.append("| " + " | ".join("---" for _ in header) + " |")
        for binding in ORDER:
            values = [binding]
            any_value = False
            for molecules in sizes:
                for width in threads:
                    row = lookup(rows, binding, op, molecules, width)
                    if row is not None:
                        any_value = True
                    values.append(cell(row))
            if any_value:
                lines.append("| " + " | ".join(values) + " |")
    return "\n".join(lines)


def _python_rings(rows):
    value = ratio(rows, "python", "rings", 4096, 1)
    if value is None:
        return ""
    return " Python's ctypes ring call is {:.2f} times C on that cell.".format(value)


def prose(rows):
    model = rows[0]["cpu_model"]
    cores = rows[0]["cpu_cores"]
    sizes = sorted({row["molecules"] for row in rows})
    largest = max(sizes)
    wrappers = ("cxx", "rust", "fortran", "python", "julia", "r", "lua")
    handle_ops = ("read", "knn", "rings")
    bits = []
    for op in handle_ops:
        for width in (1, int(cores)):
            ratios = []
            for binding in wrappers:
                value = ratio(rows, binding, op, largest, width)
                if value is not None:
                    ratios.append((value, binding, width))
            if ratios:
                worst = max(ratios)
                bits.append((worst[0], worst[1], op, worst[2]))
    handle_peak = max(bits) if bits else None
    lua_bits = []
    for op in ("chill", "steinhardt"):
        for width in (1, int(cores)):
            value = ratio(rows, "lua", op, largest, width)
            if value is not None:
                lua_bits.append((value, op, width))
    paragraphs = []
    paragraphs.append(
        "The run used {} cores of {}. No host name is recorded. "
        "One process per binding stayed open for both thread widths. "
        "The serial phase set OMP_NUM_THREADS and PLUMED_NUM_THREADS to 1. "
        "The threaded phase set both to {} and called seams_set_num_threads. "
        "Each cell below is the median of 8 repeats after one warmup, "
        "in milliseconds, with the interquartile range in parentheses.".format(
            cores, model, cores
        )
    )
    if handle_peak:
        value, binding, op, width = handle_peak
        paragraphs.append(
            "Read, the mutual 4-neighbour list, and the primitive-ring search "
            "return a handle. On the {}-molecule cell the slowest of those "
            "wrappers was {} {} at {} thread{}, {:.2f} times the direct C call. "
            "The file parser and the search run inside libseams. The wrapper "
            "passes a pointer and, for rings and neighbours, destroys the handle.".format(
                largest, binding, OP_NAME[op], width, "" if width == 1 else "s", value
            )
        )
    if lua_bits:
        value, op, width = max(lua_bits)
        paragraphs.append(
            "CHILL and q6 write one number per molecule, and the timed q6 call "
            "also writes q6bar. C, C++, Fortran, Rust, Julia, R, and the ctypes "
            "module fill one contiguous buffer and hand its address to libseams. "
            "Lua then stores every value in a table. On the {}-molecule cell, "
            "Lua {} at {} thread{} is {:.2f} times the C median. That gap is the "
            "table write.".format(
                largest, OP_NAME[op], width, "" if width == 1 else "s", value
            )
        )
    paragraphs.append(
        "Two copies that did not belong in the wrappers are gone. "
        "Lua's CSR export walked the neighbour lists into a temporary buffer "
        "and then into the table. It now reads the CSR that the handle already "
        "owns. Julia's orthorhombic frame setter copied a 3-by-n Float64 matrix "
        "that was already the layout seams_system_set_ortho copies. It now passes "
        "that pointer. Rust, C++, and Julia allocate the q6bar buffer only when "
        "the caller asks for the Lechner-Dellago average. The timed call asks "
        "for it, so both buffers are required."
    )
    guest = []
    for binding, label in (("freud", "freud"), ("ovito", "OVITO")):
        present = sorted({
            OP_NAME[row["operation"]] for row in rows if row["binding"] == binding
        })
        if present:
            guest.append("{} times {}".format(label, ", ".join(present)))
    if guest:
        paragraphs.append(
            "freud and OVITO are separate programs, timed in the same session "
            "on the same oxygen coordinates. {} "
            "freud's neighbour query is the one-way 4-nearest list inside "
            "5 angstroms, the candidate shell the mutual list uses, and its q6 "
            "call includes the averaged parameter. "
            "OVITO reads the dump, builds a 4-neighbour finder, runs the "
            "Steinhardt modifier with the averaged parameters, and runs "
            "FindRings after creating bonds at 3.5 angstroms. "
            "OVITO has no original CHILL modifier, and freud has no reader, "
            "no CHILL, and no primitive-ring search, so those cells are empty. "
            "seams_rings_find also classifies prisms. That walk is about one "
            "second at 4096 molecules, where every binding was repeated. "
            "The C search was repeated at 32768 molecules as well. "
            "A 110592-molecule call did not return within 180 seconds, "
            "so that cell is empty.{} The other languages return the "
            "ring handle and stay with the C time.".format(
                " ".join(item + "." for item in (
                    sentence[0].upper() + sentence[1:] for sentence in guest
                )),
                _python_rings(rows),
            )
        )
    paragraphs.append(
        "At 8 threads the seams neighbour, q6, and CHILL calls take longer "
        "than the one-thread calls, on every language at once. freud and "
        "OVITO move the other way and get faster. The bindings are following "
        "the library, not adding a cost of their own."
    )
    return "\n\n".join(paragraphs)


def draw(rows):
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    sizes = sorted({row["molecules"] for row in rows})
    widths = sorted({row["threads"] for row in rows})
    bindings = [name for name in ORDER if any(row["binding"] == name for row in rows)]
    fig, axes = plt.subplots(len(widths), len(OPS), figsize=(14, 6), sharex=True, squeeze=False)
    for row_i, width in enumerate(widths):
        for col, op in enumerate(OPS):
            ax = axes[row_i][col]
            for binding in bindings:
                xs = []
                ys = []
                for molecules in sizes:
                    value = ratio(rows, binding, op, molecules, width)
                    if value is None:
                        continue
                    xs.append(molecules)
                    ys.append(value)
                if xs:
                    ax.plot(xs, ys, marker="o", label=binding)
            ax.axhline(1.0, color="0.5", linewidth=0.8)
            ax.set_yscale("log")
            ax.set_title("{} , {} thread{}".format(OP_NAME[op], width, "" if width == 1 else "s"))
            if row_i == len(widths) - 1:
                ax.set_xlabel("molecules")
            if col == 0:
                ax.set_ylabel("median / C")
    handles, labels = axes[0][0].get_legend_handles_labels()
    fig.legend(handles, labels, loc="upper center", ncol=min(len(labels), 8))
    fig.tight_layout(rect=(0, 0, 1, 0.92))
    PNG.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(PNG, dpi=120)
    plt.close(fig)


def write_pages(rows):
    body = prose(rows)
    tables = table_lines(rows)
    MD.parent.mkdir(parents=True, exist_ok=True)
    MD.write_text(
        "# Binding overhead\n\n"
        "Ice Ih at 4096, 32768, and 110592 molecules. "
        "GenIce builds a missing cell. The oxygen sites are kept and scaled "
        "from nanometres to angstroms.\n\n"
        + body + "\n" + tables + "\n\n"
        "The figure is drawn when the docs are built. "
        "It is the median ratio to the direct C call.\n",
        encoding="utf-8",
    )
    RST.parent.mkdir(parents=True, exist_ok=True)
    RST.write_text(
        "Binding overhead\n"
        "================\n\n"
        "Ice Ih at 4096, 32768, and 110592 molecules. "
        "GenIce builds a missing cell. The oxygen sites are kept and scaled "
        "from nanometres to angstroms.\n\n"
        + body + "\n\n"
        ".. figure:: ../_static/binding-overhead.png\n"
        "   :alt: Median time of each binding divided by the direct C call\n\n"
        "   Median time divided by the direct C call. "
        "The vertical axis is logarithmic so a wrapper near 1 and a separate "
        "analysis sit on the same axes.\n\n"
        "The table is milliseconds, median and interquartile range.\n\n"
        ".. code-block:: text\n\n"
        + "\n".join("   " + line if line else "" for line in tables.splitlines())
        + "\n\n"
        "``benchmarks/binding-overhead.csv`` is the summary. "
        "The sample log is written beside the generated cells and is not stored in git. "
        "``python3 scripts/bench_binding_overhead.py`` repeats the session. "
        "The Snakemake rule ``binding_overhead`` calls that script.\n",
        encoding="utf-8",
    )
    ORG.parent.mkdir(parents=True, exist_ok=True)
    ORG.write_text(
        "#+TITLE: Binding overhead\n\n"
        "Ice Ih at 4096, 32768, and 110592 molecules.\n\n"
        + body + "\n\n"
        "The Sphinx page is docs/source/explanation/binding-overhead.rst. "
        "The figure is drawn when the docs are built.\n",
        encoding="utf-8",
    )


def main(argv=None):
    del argv
    if not CSV.is_file():
        print("missing", CSV, file=sys.stderr)
        return 1
    rows = load()
    if not rows:
        print("empty csv", file=sys.stderr)
        return 1
    draw(rows)
    write_pages(rows)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
