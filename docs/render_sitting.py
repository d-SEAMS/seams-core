#!/usr/bin/env python3
# Copyright (c) d-SEAMS developers
"""Build the rc3 sitting tables from benchmarks/rc3-sitting.csv.

The median of an even sample is the mean of the two central values.
The interquartile range is the type-7 quartile gap, the same linear
position q*(n-1) that scripts/h2h_protocol.py uses. --check recomputes
both from samples_ms and refuses a row that does not match the stored
summary. The Sphinx configuration calls write_tables before it reads
the book, so the numbers are not copied by hand.
"""

import csv
import pathlib
import sys

ROOT = pathlib.Path(__file__).resolve().parents[1]
CSV = ROOT / "benchmarks" / "rc3-sitting.csv"
TRACE = ROOT / "benchmarks" / "trace-threads.csv"
POOL = ROOT / "benchmarks" / "pool-audit.csv"
LARGE = ROOT / "benchmarks" / "paper-cells" / "sitting-large.csv"
CUTOFF = ROOT / "benchmarks" / "cutoff-rings.csv"
OVERVIEW_OUT = ROOT / "docs" / "source" / "explanation" / "performance-overview-tables.rst"
LARGE_OUT = ROOT / "docs" / "source" / "explanation" / "sitting-large-tables.rst"
CUTOFF_OUT = ROOT / "docs" / "source" / "explanation" / "cutoff-rings-tables.rst"
OUT = ROOT / "docs" / "source" / "explanation" / "performance-tables.rst"
TRACE_OUT = ROOT / "docs" / "source" / "explanation" / "trace-threads-tables.rst"
POOL_OUT = ROOT / "docs" / "source" / "explanation" / "pool-audit-tables.rst"
POOL_KERNELS = (
    "rings",
    "q4-local",
    "q6-local",
    "q4-cutoff",
    "q6-cutoff",
    "chill",
    "symmetry",
    "hop",
)
POOL_SIZES = ("4096", "32768", "110592", "884736")
H2H_OUT = ROOT / "docs" / "source" / "howto" / "literature-speed.rst"
# Guest codes in the rc3 reproduction how-to. Those PLUMED columns were
# recorded with PLUMED_NUM_THREADS unset, so both are one thread. The
# 884736 sitting sets the variable and is not in this list.
H2H_FALL = (
    ("freud", "q4"),
    ("freud", "q6"),
    ("OVITO", "CHILL+"),
    ("OVITO", "PTM"),
    ("OVITO", "rings"),
    ("mdapy", "CHILL+"),
    ("mdapy", "q4"),
    ("mdapy", "q6"),
    ("mdapy", "rings"),
    ("vesin", "neighbours"),
    ("rings-statistics-matter", "rings"),
)
H2H_FLAT = (("PLUMED", "q4"), ("PLUMED", "q6"))
H2H_SIZES = ("110592", "884736")

ORDER = ("neighbours", "q4", "q6", "CHILL+", "rings", "cages", "PTM", "king")


def quantile_linear(sorted_vals, q):
    if len(sorted_vals) == 1:
        return sorted_vals[0]
    pos = q * (len(sorted_vals) - 1)
    lo = int(pos)
    frac = pos - lo
    if lo + 1 >= len(sorted_vals):
        return sorted_vals[-1]
    return sorted_vals[lo] * (1.0 - frac) + sorted_vals[lo + 1] * frac


def summarize(samples):
    ordered = sorted(float(v) for v in samples)
    count = len(ordered)
    mid = count // 2
    if count % 2 == 0:
        median = 0.5 * (ordered[mid - 1] + ordered[mid])
    else:
        median = ordered[mid]
    return count, median, quantile_linear(ordered, 0.75) - quantile_linear(ordered, 0.25)


def load(path):
    with path.open(newline="") as handle:
        return list(csv.DictReader(handle))


def check(rows, tolerance):
    errors = []
    for index, row in enumerate(rows, start=2):
        samples = [float(part) for part in row["samples_ms"].split()]
        count, median, iqr = summarize(samples)
        if int(row["n"]) != count:
            errors.append("line {}: n {} != {} samples".format(index, row["n"], count))
        # samples_ms is rounded, so the recomputed summary can move by
        # one stored unit from the value written at full precision.
        if abs(median - float(row["median_ms"])) > tolerance:
            errors.append("line {}: median {} != {:.6f}".format(index, row["median_ms"], median))
        if abs(iqr - float(row["iqr_ms"])) > tolerance:
            errors.append("line {}: iqr {} != {:.6f}".format(index, row["iqr_ms"], iqr))
    return errors


def cell(text):
    return text.replace("\n", " ")


def write_trace(rows, path):
    grouped = {}
    for row in rows:
        grouped.setdefault(row["crystal"], []).append(row)
    lines = [
        ".. list-table::",
        "   :header-rows: 0",
        "",
        "   * - These tables are written by ``docs/render_sitting.py`` from",
        "       ``benchmarks/trace-threads.csv``. Edit the CSV, not this file.",
        "",
    ]
    header = ("cell", "molecules", "threads", "cores", "n", "median_ms", "iqr_ms")
    for name in sorted(grouped):
        block = sorted(grouped[name], key=lambda row: (int(row["molecules"]), int(row["threads"])))
        title = "signature cages, {}".format(name)
        lines.append(title)
        lines.append("^" * len(title))
        lines.append("")
        lines.append(".. csv-table::")
        lines.append("   :header: " + ", ".join(header))
        lines.append("")
        for row in block:
            lines.append("   " + ", ".join(cell(row[key]) for key in header))
        lines.append("")
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text("\n".join(lines) + "\n")


def check_pool(rows):
    errors = check(rows, 1.5e-4)
    if len(rows) != 128:
        errors.append("pool-audit has {} rows, expected 128".format(len(rows)))
    widths = sorted({int(row["threads"]) for row in rows})
    if widths != [1, 2, 4, 8]:
        errors.append("pool threads are {}, expected 1 2 4 8".format(widths))
    kernels = sorted({row["kernel"] for row in rows})
    if kernels != sorted(POOL_KERNELS):
        errors.append("pool kernels are {}".format(kernels))
    by_size = {}
    fingerprints = {}
    for row in rows:
        by_size[row["molecules"]] = by_size.get(row["molecules"], 0) + 1
        key = (row["kernel"], row["molecules"])
        fingerprints.setdefault(key, set()).add(row["fingerprint"])
    for size in POOL_SIZES:
        if by_size.get(size) != 32:
            errors.append("molecules {} has {} rows, expected 32".format(size, by_size.get(size, 0)))
    for key, found in sorted(fingerprints.items()):
        if len(found) != 1:
            errors.append("{} fingerprint changes with threads".format(key))
    return errors


def write_pool(rows, path):
    grouped = {}
    order = {name: index for index, name in enumerate(POOL_KERNELS)}
    for row in rows:
        grouped.setdefault(row["molecules"], []).append(row)
    lines = [
        ".. list-table::",
        "   :header-rows: 0",
        "",
        "   * - Written by ``docs/render_sitting.py`` from ``benchmarks/pool-audit.csv``.",
        "       Each size is one 32-row table: eight kernels at 1, 2, 4 and 8 threads.",
        "       ``median_ms`` and ``iqr_ms`` are recomputed from ``samples_ms``.",
        "",
    ]
    header = ("kernel", "threads", "cores", "n", "median_ms", "iqr_ms")
    for molecules in POOL_SIZES:
        block = sorted(
            grouped.get(molecules, []),
            key=lambda row: (order.get(row["kernel"], 99), int(row["threads"])),
        )
        title = "ice Ic, {} molecules".format(molecules)
        lines.append(title)
        lines.append("^" * len(title))
        lines.append("")
        lines.append(".. csv-table::")
        lines.append("   :header: " + ", ".join(header))
        lines.append("")
        for row in block:
            lines.append("   " + ", ".join(cell(row[key]) for key in header))
        lines.append("")
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text("\n".join(lines) + "\n")


def h2h_index(rows):
    found = {}
    for row in rows:
        found[(row["code"], row["operation"], row["molecules"], row["threads"])] = float(row["median_ms"])
    return found


def h2h_pairs(found):
    pairs = []
    for code, operation in H2H_FALL + H2H_FLAT:
        for molecules in H2H_SIZES:
            left = found.get((code, operation, molecules, "1"))
            right = found.get((code, operation, molecules, "8"))
            if left is None or right is None:
                continue
            pairs.append((code, operation, molecules, left, right, left / right))
    return pairs


def check_h2h(rows):
    errors = []
    found = h2h_index(rows)
    seen = {(code, operation, molecules) for code, operation, molecules, _a, _b, _r in h2h_pairs(found)}
    for code, operation in H2H_FALL + H2H_FLAT:
        for molecules in H2H_SIZES:
            if (code, operation, molecules) not in seen:
                errors.append("missing {} {} {}".format(code, operation, molecules))
    for code, operation, molecules, left, right, ratio in h2h_pairs(found):
        if (code, operation) in H2H_FLAT:
            if not 0.95 <= ratio <= 1.05:
                errors.append(
                    "{} {} {} columns differ ({:.3f}); both were recorded with PLUMED_NUM_THREADS unset".format(
                        code, operation, molecules, ratio
                    )
                )
        elif not right < left:
            errors.append("{} {} {} did not fall ({} -> {})".format(code, operation, molecules, left, right))
    return errors


def write_h2h(rows, path):
    lines = [
        ".. list-table::",
        "   :header-rows: 0",
        "",
        "   * - Written by ``docs/render_sitting.py`` from ``benchmarks/rc3-sitting.csv``.",
        "       ``t1`` and ``t8`` are the medians in the threads=1 and threads=8 columns.",
        "       ``ratio`` is ``t1 / t8``.",
        "       A PLUMED row was measured on one thread because ``PLUMED_NUM_THREADS`` was unset.",
        "       Both of its columns are that measurement.",
        "",
        "One thread against eight",
        "^^^^^^^^^^^^^^^^^^^^^^^^",
        "",
        ".. csv-table::",
        "   :header: code, operation, molecules, t1_ms, t8_ms, ratio, note",
        "",
    ]
    for code, operation, molecules, left, right, ratio in h2h_pairs(h2h_index(rows)):
        note = ""
        if code == "PLUMED":
            note = "one thread; PLUMED_NUM_THREADS unset"
        lines.append(
            "   {}, {}, {}, {:.3f}, {:.3f}, {:.2f}, {}".format(
                code, operation, molecules, left, right, ratio, note
            )
        )
    lines.append("")
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text("\n".join(lines) + "\n")


def check_large(rows):
    errors = check(rows, 1.5e-4)
    if len(rows) != 104:
        errors.append("sitting-large has {} rows, expected 104".format(len(rows)))
    widths = sorted({int(row["threads"]) for row in rows})
    if widths != [1, 2, 4, 8]:
        errors.append("sitting-large threads are {}".format(widths))
    sizes = sorted({row["molecules"] for row in rows})
    if sizes != ["884736"]:
        errors.append("sitting-large molecules are {}".format(sizes))
    # The detail column is the agreement check on one shared ice Ic frame.
    for row in rows:
        detail = row.get("detail", "")
        code, operation = row["code"], row["operation"]
        if code in ("d-SEAMS", "rings-statistics-matter") and operation == "rings":
            if "agree_six=" not in detail:
                errors.append("rings {} threads {} missing six-ring agreement".format(code, row["threads"]))
        if operation == "q6" and code in ("freud", "mdapy"):
            if "agree_freud_max_abs=" not in detail or "agree_mdapy_max_abs=" not in detail:
                errors.append("q6 {} threads {} missing the q6 agreement".format(code, row["threads"]))
        if operation == "q6" and code in ("d-SEAMS", "pyscal3", "PLUMED"):
            if "agree_pyscal3_max_abs=" not in detail or "agree_PLUMED_max_abs=" not in detail:
                errors.append("q6 {} threads {} missing the pyscal3/PLUMED agreement".format(code, row["threads"]))
            if code == "PLUMED" and "PLUMED_NUM_THREADS=" not in detail:
                errors.append("PLUMED q6 threads {} did not set PLUMED_NUM_THREADS".format(row["threads"]))
        if operation == "neighbours" and code in ("vesin", "freud"):
            if "agree=d-SEAMS/vesin/freud" not in detail:
                errors.append("neighbours {} threads {} missing the pair agreement".format(code, row["threads"]))
        if operation == "neighbours" and code in ("d-SEAMS", "pyscal3"):
            if "agree=d-SEAMS,pyscal3" not in detail:
                errors.append("neighbours {} threads {} missing the pyscal3 agreement".format(code, row["threads"]))
        if operation == "coordination" and code == "OVITO":
            if "coordination_sum=" not in detail or "dseams_pairs=" not in detail:
                errors.append("coordination threads {} missing the pair sum".format(row["threads"]))
        if operation == "chill+" and code in ("d-SEAMS", "OVITO"):
            if "cubic_dseams=" not in detail or "cubic_ovito=" not in detail:
                errors.append("chill+ {} threads {} missing the cubic agreement".format(code, row["threads"]))
        if (code, operation) in (("d-SEAMS", "classifier"), ("OVITO", "ptm-ice")):
            if "agree_ice_ic_cubic_diamond=" not in detail:
                errors.append("{} {} threads {} missing the classifier agreement".format(code, operation, row["threads"]))
    return errors


def write_large(rows, path):
    grouped = {}
    for row in rows:
        grouped.setdefault(row["operation"], []).append(row)
    order = (
        "rings",
        "q4",
        "q6",
        "q4+q6",
        "neighbours",
        "coordination",
        "chill+",
        "classifier",
        "ptm-default",
        "ptm-ice",
    )
    names = [name for name in order if name in grouped]
    names.extend(sorted(set(grouped) - set(names)))
    lines = [
        ".. list-table::",
        "   :header-rows: 0",
        "",
        "   * - Written by ``docs/render_sitting.py`` from ``benchmarks/paper-cells/sitting-large.csv``.",
        "       Ice Ic at 884736 molecules. Each code is timed at 1, 2, 4 and 8 threads.",
        "       ``median_ms`` and ``iqr_ms`` are recomputed from ``samples_ms``.",
        "",
    ]
    header = ("code", "version", "operation", "threads", "n", "median_ms", "iqr_ms")
    for name in names:
        block = sorted(grouped[name], key=lambda row: (row["code"], int(row["threads"])))
        title = "884736 molecules, {}".format(name)
        lines.append(title)
        lines.append("^" * len(title))
        lines.append("")
        lines.append(".. csv-table::")
        lines.append("   :header: " + ", ".join(header))
        lines.append("")
        for row in block:
            lines.append("   " + ", ".join(cell(row[key]) for key in header))
        lines.append("")
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text("\n".join(lines) + "\n")


def check_cutoff(rows):
    errors = []
    timed = [row for row in rows if row.get("samples_ms")]
    errors.extend(check(timed, 1.5e-3))
    if len(rows) != 32:
        errors.append("cutoff rings has {} rows, expected 32".format(len(rows)))
    widths = sorted({int(row["threads"]) for row in rows})
    if widths != [1, 2, 4, 8]:
        errors.append("cutoff rings threads are {}".format(widths))
    codes = sorted({row["code"] for row in rows})
    if codes != ["cell-list", "linkcell"]:
        errors.append("cutoff rings codes are {}".format(codes))
    return errors


def write_cutoff(rows, path):
    order = {"total": 0, "link": 1, "convert": 2, "search": 3}
    block = sorted(rows, key=lambda row: (order.get(row["operation"], 9), row["code"], int(row["threads"])))
    lines = [
        ".. list-table::",
        "   :header-rows: 0",
        "",
        "   * - Written by ``docs/render_sitting.py`` from ``benchmarks/cutoff-rings.csv``.",
        "       Ice Ic at 884736 molecules. ``linkcell`` is the 1ce3be3b builder.",
        "       ``cell-list`` is the orthorhombic list inside the ring call.",
        "       ``total`` medians are recomputed from ``samples_ms``.",
        "",
        "Full cutoff ring call",
        "^^^^^^^^^^^^^^^^^^^^^^",
        "",
        ".. csv-table::",
        "   :header: code, operation, threads, n, median_ms, iqr_ms",
        "",
    ]
    header = ("code", "operation", "threads", "n", "median_ms", "iqr_ms")
    for row in block:
        lines.append("   " + ", ".join(cell(row[key]) for key in header))
    lines.append("")
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text("\n".join(lines) + "\n")


def ratio_series(rows, fields, unset_plumed=False):
    grouped = {}
    for row in rows:
        key = tuple(row[field] for field in fields)
        grouped.setdefault(key, {})[row["threads"]] = float(row["median_ms"])
    series = []
    for key in sorted(grouped):
        widths = grouped[key]
        if "1" not in widths or "8" not in widths or widths["8"] <= 0:
            continue
        ratio = widths["1"] / widths["8"]
        if unset_plumed and "PLUMED" in key:
            label = "one thread; PLUMED_NUM_THREADS unset"
        elif 0.95 <= ratio <= 1.05:
            label = "flat"
        elif ratio > 1:
            label = "falls"
        else:
            label = "rises"
        series.append(key + (widths["1"], widths["8"], ratio, label))
    return series


def append_ratio_table(lines, title, rows, fields, unset_plumed=False):
    lines.append(title)
    lines.append("^" * len(title))
    lines.append("")
    header = list(fields) + ["t1_ms", "t8_ms", "ratio", "with_threads"]
    lines.append(".. csv-table::")
    lines.append("   :header: " + ", ".join(header))
    lines.append("")
    for item in ratio_series(rows, fields, unset_plumed):
        text = [str(part) for part in item[:-4]]
        left, right, ratio, label = item[-4:]
        text.extend(["{:.3f}".format(left), "{:.3f}".format(right), "{:.3f}".format(ratio), label])
        lines.append("   " + ", ".join(text))
    lines.append("")


def write_overview(books, path):
    lines = [
        ".. list-table::",
        "   :header-rows: 0",
        "",
        "   * - Written by ``docs/render_sitting.py`` from the four CSVs.",
        "       ``ratio`` is the threads=1 median divided by the threads=8 median.",
        "       ``flat`` means that ratio is between 0.95 and 1.05.",
        "       A PLUMED row in the rc3 table says ``one thread; PLUMED_NUM_THREADS unset``.",
        "       Both of its columns were measured on one thread because that variable was unset.",
        "       The large-cell PLUMED row set ``PLUMED_NUM_THREADS``, so it is not that label.",
        "",
        "What each file holds",
        "^^^^^^^^^^^^^^^^^^^^",
        "",
        ".. csv-table::",
        "   :header: file, rows, threads, molecules",
        "",
    ]
    for name, rows in books:
        threads = " ".join(str(width) for width in sorted({int(row["threads"]) for row in rows}))
        molecules = " ".join(sorted({row["molecules"] for row in rows}, key=int))
        lines.append("   {}, {}, {}, {}".format(name, len(rows), threads, molecules))
    lines.append("")
    append_ratio_table(
        lines,
        "rc3 sitting, one thread against eight",
        books[0][1],
        ("code", "operation", "molecules"),
        unset_plumed=True,
    )
    append_ratio_table(
        lines,
        "884736-molecule sitting, one thread against eight",
        books[1][1],
        ("code", "operation"),
    )
    append_ratio_table(
        lines,
        "signature search, one thread against eight",
        books[2][1],
        ("crystal", "cell", "molecules"),
    )
    append_ratio_table(
        lines,
        "pool-sized kernels, one thread against eight",
        books[3][1],
        ("kernel", "molecules"),
    )
    lines.append("Pool audit at 4096 molecules, four threads against eight")
    lines.append("^" * len("Pool audit at 4096 molecules, four threads against eight"))
    lines.append("")
    lines.append(".. csv-table::")
    lines.append("   :header: kernel, t4_ms, t8_ms, t8_over_t4")
    lines.append("")
    four = {}
    eight = {}
    for row in books[3][1]:
        if row["molecules"] != "4096":
            continue
        if row["threads"] == "4":
            four[row["kernel"]] = float(row["median_ms"])
        elif row["threads"] == "8":
            eight[row["kernel"]] = float(row["median_ms"])
    for kernel in POOL_KERNELS:
        if kernel not in four or kernel not in eight or eight[kernel] <= 0:
            continue
        lines.append(
            "   {}, {:.3f}, {:.3f}, {:.3f}".format(
                kernel, four[kernel], eight[kernel], eight[kernel] / four[kernel]
            )
        )
    lines.append("")
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text("\n".join(lines) + "\n")


def write_tables(rows, path):
    grouped = {}
    for row in rows:
        grouped.setdefault(row["operation"], []).append(row)
    names = [name for name in ORDER if name in grouped]
    names.extend(sorted(set(grouped) - set(names)))
    lines = [
        ".. list-table::",
        "   :header-rows: 0",
        "",
        "   * - These tables are written by ``docs/render_sitting.py`` from",
        "       ``benchmarks/rc3-sitting.csv``. Edit the CSV, not this file.",
        "",
    ]
    header = ("code", "version", "molecules", "threads", "cores", "n", "median_ms", "iqr_ms")
    for name in names:
        block = sorted(
            grouped[name],
            key=lambda row: (int(row["molecules"]), int(row["threads"]), row["code"], row["version"]),
        )
        lines.append(name)
        lines.append("^" * len(name))
        lines.append("")
        lines.append(".. csv-table::")
        lines.append("   :header: " + ", ".join(header))
        lines.append("")
        for row in block:
            lines.append("   " + ", ".join(cell(row[key]) for key in header))
        lines.append("")
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text("\n".join(lines) + "\n")


def main(argv):
    if "--check-h2h" in argv:
        rows = load(CSV)
        errors = check_h2h(rows)
        if errors:
            print("\n".join(errors), file=sys.stderr)
            return 1
        print("h2h thread scaling ok", len(h2h_pairs(h2h_index(rows))), "rows")
        return 0
    if "--check-cutoff" in argv:
        rows = load(CUTOFF)
        errors = check_cutoff(rows)
        if errors:
            print("\n".join(errors), file=sys.stderr)
            return 1
        print("cutoff rings csv ok", len(rows), "rows")
        return 0
    if "--check-large" in argv:
        rows = load(LARGE)
        errors = check_large(rows)
        if errors:
            print("\n".join(errors), file=sys.stderr)
            return 1
        print("sitting-large csv ok", len(rows), "rows")
        return 0
    if "--check-pool" in argv:
        rows = load(POOL)
        errors = check_pool(rows)
        if errors:
            print("\n".join(errors), file=sys.stderr)
            return 1
        print("pool-audit csv ok", len(rows), "rows")
        return 0
    if "--check-trace" in argv:
        rows = load(TRACE)
        errors = check(rows, 1.5e-4)
        widths = sorted({int(row["threads"]) for row in rows})
        if widths != [1, 2, 4, 8]:
            errors.append("threads are {}, expected 1 2 4 8".format(widths))
        if errors:
            print("\n".join(errors), file=sys.stderr)
            return 1
        print("trace-threads csv ok", len(rows), "rows")
        return 0
    rows = load(CSV)
    errors = check(rows, 1.5e-3)
    if errors:
        print("\n".join(errors), file=sys.stderr)
        return 1
    if "--check" in argv:
        print("sitting csv ok", len(rows), "rows")
        return 0
    write_tables(rows, OUT)
    trace_rows = load(TRACE)
    trace_errors = check(trace_rows, 1.5e-4)
    if trace_errors:
        print("\n".join(trace_errors), file=sys.stderr)
        return 1
    write_trace(trace_rows, TRACE_OUT)
    h2h_errors = check_h2h(rows)
    if h2h_errors:
        print("\n".join(h2h_errors), file=sys.stderr)
        return 1
    write_h2h(rows, H2H_OUT)
    pool_rows = load(POOL)
    pool_errors = check_pool(pool_rows)
    if pool_errors:
        print("\n".join(pool_errors), file=sys.stderr)
        return 1
    write_pool(pool_rows, POOL_OUT)
    large_rows = load(LARGE)
    large_errors = check_large(large_rows)
    if large_errors:
        print("\n".join(large_errors), file=sys.stderr)
        return 1
    write_large(large_rows, LARGE_OUT)
    cutoff_rows = load(CUTOFF)
    cutoff_errors = check_cutoff(cutoff_rows)
    if cutoff_errors:
        print("\n".join(cutoff_errors), file=sys.stderr)
        return 1
    write_cutoff(cutoff_rows, CUTOFF_OUT)
    write_overview(
        (
            ("rc3-sitting.csv", rows),
            ("sitting-large.csv", large_rows),
            ("trace-threads.csv", trace_rows),
            ("pool-audit.csv", pool_rows),
        ),
        OVERVIEW_OUT,
    )
    print("wrote", OUT.relative_to(ROOT), len(rows), "rows")
    print("wrote", TRACE_OUT.relative_to(ROOT), len(trace_rows), "rows")
    return 0


if __name__ == "__main__":
    raise SystemExit(main(sys.argv[1:]))
