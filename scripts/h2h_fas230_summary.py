#!/usr/bin/env python3
# Copyright (c) d-SEAMS developers
"""Median and 95% bootstrap interval for the fastatomstruct 2.3.0 sitting.

The arm order is one shuffle of the campaign arms from seed 230220,
using Python 3.11+ Random.shuffle (randbelow). Even rounds run that
order. Odd rounds run it reversed. A cell is one operation and one
ice Ic size. Twenty rounds, no replacement of a sample.

The 95% interval is the 2.5 and 97.5 percentiles of 10000 bootstrap
medians. The percentile indexes are int(0.025 * draws) and
int(0.975 * draws) - 1, the same indexes as bench/measure.py. The
reported median is the median of the original samples. An even count
uses the mean of the two central values.

The ratio is time(fastatomstruct 2.3.0) / time(fastatomstruct 2.2.0),
paired by round. The A/A band is the same interval for two 2.2.0
arms. The ratio is faster or slower only when its interval lies
entirely outside that band. A ring-count or q6 disagreement is a
correctness finding, not a timing call.
"""

from __future__ import annotations

import argparse
import csv
import random
import statistics
import sys
from collections import defaultdict

SEED = 230220
ROUNDS = 20
BOOT_DRAWS = 10000
BOOT_SEED = 230220
Q6_ATOL = 1e-6
OPERATIONS = ("q6", "rings")
MOLECULES = (1000, 10648, 110592)
REPETITIONS = {1000: 5, 10648: 11, 110592: 24}
ARMS_AB = ("seams", "fas220", "fas230")
ARMS_AA = ("fas220a", "fas220b")
CAMPAIGNS = (("ab", ARMS_AB), ("aa", ARMS_AA))

SUMMARY_FIELDS = (
    "row_kind",
    "campaign",
    "operation",
    "molecules",
    "code",
    "version",
    "arm",
    "n",
    "median_ms",
    "ci_lo_ms",
    "ci_hi_ms",
    "median_ratio",
    "ratio_lo",
    "ratio_hi",
    "aa_median_ratio",
    "aa_lo",
    "aa_hi",
    "call",
    "agree",
    "detail",
    "boot_draws",
    "boot_seed",
    "code_revision",
    "host",
    "build",
)


def base_order(arms, seed):
    order = list(arms)
    random.Random(seed).shuffle(order)
    return tuple(order)


def round_order(base, round_index):
    if round_index % 2 == 0:
        return tuple(base)
    return tuple(reversed(base))


def iter_jobs():
    """q6 before rings, smaller cells first, the version sitting before A/A."""
    for operation in OPERATIONS:
        for molecules in MOLECULES:
            for campaign, arms in CAMPAIGNS:
                base = base_order(arms, SEED)
                for rnd in range(ROUNDS):
                    order = round_order(base, rnd)
                    yield {
                        "campaign": campaign,
                        "operation": operation,
                        "molecules": molecules,
                        "round": rnd,
                        "arm_order": order,
                    }


def boot_seed(*parts):
    text = ",".join(str(part) for part in parts)
    hash_value = 2166136261
    for char in text:
        hash_value ^= ord(char)
        hash_value = (hash_value * 16777619) & 0xFFFFFFFF
    return BOOT_SEED + (hash_value % 1000000)


def bootstrap_median_ci(values, draws, seed):
    if len(values) < 2:
        raise ValueError(f"need at least two samples, got {len(values)}")
    if draws < 100:
        raise ValueError("bootstrap draws")
    rng = random.Random(seed)
    count = len(values)
    medians = []
    for _ in range(draws):
        sample = [values[rng.randrange(count)] for _ in range(count)]
        medians.append(statistics.median(sample))
    medians.sort()
    lo_i = int(0.025 * draws)
    hi_i = int(0.975 * draws) - 1
    if lo_i < 0 or hi_i >= draws or lo_i > hi_i:
        raise ValueError("bootstrap percentile index")
    return statistics.median(values), medians[lo_i], medians[hi_i]


def outside(lo, hi, band_lo, band_hi):
    if hi < band_lo:
        return "faster"
    if lo > band_hi:
        return "slower"
    return "within"


def _num(text):
    if text is None or text == "":
        return None
    return float(text)


def _cell_key(row):
    return (row["campaign"], row["operation"], int(row["molecules"]))


def _check_row(row, expected):
    if row["seed"] != str(SEED):
        raise SystemExit("sample seed does not match the recorded seed")
    if int(row["warmup"]) != 0:
        raise SystemExit("a sample was pre-warmed")
    if row["periodic"] != "1":
        raise SystemExit("a sample was not periodic")
    if row["cpus"] != "0" or row["threads"] != "1":
        raise SystemExit("a sample was not pinned to one core and one thread")
    if int(row["round"]) != expected["round"] or int(row["position"]) != expected["position"]:
        raise SystemExit("sample position does not match the schedule")
    if row["arm"] != expected["arm"]:
        raise SystemExit("sample arm does not match the schedule")
    if row["arm_order"] != ",".join(expected["arm_order"]):
        raise SystemExit("sample order does not match the schedule")
    molecules = int(row["molecules"])
    if int(row["repetitions"]) != REPETITIONS[molecules]:
        raise SystemExit("repetitions do not match the ice Ic cell")
    if not row["code_revision"] or not row["host"] or not row["build"]:
        raise SystemExit("sample is missing code_revision, host, or build")


def _index(records):
    found = {}
    for row in records:
        key = (
            row["campaign"],
            row["operation"],
            int(row["molecules"]),
            int(row["round"]),
            row["arm"],
        )
        if key in found:
            raise SystemExit("duplicate sample")
        found[key] = row
    return found


def _expected_slots():
    slots = {}
    for job in iter_jobs():
        for position, arm in enumerate(job["arm_order"]):
            key = (job["campaign"], job["operation"], job["molecules"], job["round"], arm)
            slots[key] = {**job, "position": position, "arm": arm}
    return slots


def _census(row):
    return tuple(int(row[name]) for name in ("ring3", "ring4", "ring5", "ring6"))


def _q6_detail(rows):
    worst = max(float(row["max_abs"]) for row in rows)
    return worst, f"max_abs {worst:.6e}"


def _ring_detail(rows):
    by_arm = {}
    stable = True
    for row in rows:
        census = _census(row)
        previous = by_arm.get(row["arm"])
        if previous is None:
            by_arm[row["arm"]] = census
        elif previous != census:
            stable = False
    text = "; ".join(
        f"{arm} {' '.join(str(v) for v in counts)}" for arm, counts in sorted(by_arm.items())
    )
    same = stable and len(set(by_arm.values())) == 1
    return same, text


def _agree_cell(rows, operation):
    flags = {row["agree"] for row in rows}
    if flags != {"1"} and flags != {"0"} and flags != {"0", "1"}:
        raise SystemExit("agree column is not 0 or 1")
    if operation == "rings":
        same, detail = _ring_detail(rows)
        return same and flags == {"1"}, detail
    worst, detail = _q6_detail(rows)
    return worst <= Q6_ATOL and flags == {"1"}, detail


def _identity(rows):
    revision = {row["code_revision"] for row in rows}
    host = {row["host"] for row in rows}
    if len(revision) != 1 or len(host) != 1:
        raise SystemExit("code_revision or host changes inside a cell")
    return revision.pop(), host.pop()


def _times(rows, arm):
    chosen = [row for row in rows if row["arm"] == arm]
    chosen.sort(key=lambda row: int(row["round"]))
    if [int(row["round"]) for row in chosen] != list(range(ROUNDS)):
        raise SystemExit("rounds are not 0..19")
    return [float(row["time_ms"]) for row in chosen]


def _version(rows, arm):
    values = {row["version"] for row in rows if row["arm"] == arm}
    codes = {row["code"] for row in rows if row["arm"] == arm}
    builds = {row["build"] for row in rows if row["arm"] == arm}
    if len(values) != 1 or len(codes) != 1 or len(builds) != 1:
        raise SystemExit("version, code, or build changes inside an arm")
    return codes.pop(), values.pop(), builds.pop()


def _blank_ratio():
    return {
        "median_ratio": "",
        "ratio_lo": "",
        "ratio_hi": "",
        "aa_median_ratio": "",
        "aa_lo": "",
        "aa_hi": "",
        "call": "",
    }


def _time_row(campaign, operation, molecules, arm, rows, agree, detail, revision, host):
    code, version, build = _version(rows, arm)
    samples = _times(rows, arm)
    seed = boot_seed("time", campaign, operation, molecules, arm)
    median, lo, hi = bootstrap_median_ci(samples, BOOT_DRAWS, seed)
    row = {
        "row_kind": "time",
        "campaign": campaign,
        "operation": operation,
        "molecules": str(molecules),
        "code": code,
        "version": version,
        "arm": arm,
        "n": str(len(samples)),
        "median_ms": f"{median:.6f}",
        "ci_lo_ms": f"{lo:.6f}",
        "ci_hi_ms": f"{hi:.6f}",
        "agree": "1" if agree else "0",
        "detail": detail,
        "boot_draws": str(BOOT_DRAWS),
        "boot_seed": str(seed),
        "code_revision": revision,
        "host": host,
        "build": build,
    }
    row.update(_blank_ratio())
    return row


def _ratio_values(rows, left, right):
    paired = []
    for rnd in range(ROUNDS):
        got = {row["arm"]: float(row["time_ms"]) for row in rows if int(row["round"]) == rnd}
        if left not in got or right not in got:
            raise SystemExit("a round is missing an arm")
        if got[left] <= 0.0 or got[right] <= 0.0:
            raise SystemExit("non-positive sample")
        paired.append(got[right] / got[left])
    return paired


def assemble(records, full):
    slots = _expected_slots()
    found = _index(records)
    if full:
        missing = [key for key in slots if key not in found]
        extra = [key for key in found if key not in slots]
        if missing or extra:
            raise SystemExit(
                f"sample grid is incomplete: missing {len(missing)} extra {len(extra)}"
            )
    else:
        extra = [key for key in found if key not in slots]
        if extra:
            raise SystemExit("sample is outside the schedule")
    for key, row in found.items():
        _check_row(row, slots[key])

    grouped = defaultdict(list)
    for key, row in found.items():
        grouped[(key[0], key[1], key[2])].append(row)

    if full:
        for campaign, _arms in CAMPAIGNS:
            for operation in OPERATIONS:
                for molecules in MOLECULES:
                    if (campaign, operation, molecules) not in grouped:
                        raise SystemExit("missing cell")

    summaries = []
    ratio_inputs = {}
    for (campaign, operation, molecules), rows in sorted(
        grouped.items(), key=lambda item: (item[0][1], item[0][2], item[0][0])
    ):
        agree, detail = _agree_cell(rows, operation)
        revision, host = _identity(rows)
        arms = ARMS_AB if campaign == "ab" else ARMS_AA
        for arm in arms:
            summaries.append(
                _time_row(
                    campaign, operation, molecules, arm, rows, agree, detail, revision, host
                )
            )
        bag = ratio_inputs.setdefault((operation, molecules), {})
        if campaign == "ab":
            bag["ab"] = _ratio_values(rows, "fas220", "fas230")
            bag["agree"] = agree
            bag["detail"] = detail
            bag["revision"] = revision
            bag["host"] = host
            bag["build"] = _version(rows, "fas230")[2]
        else:
            bag["aa"] = _ratio_values(rows, "fas220a", "fas220b")
            bag["aa_agree"] = agree

    for operation in OPERATIONS:
        for molecules in MOLECULES:
            bag = ratio_inputs.get((operation, molecules))
            if bag is None:
                continue
            if "ab" not in bag or "aa" not in bag:
                if full:
                    raise SystemExit("ratio is missing a campaign")
                continue
            seed = boot_seed("ratio", "ab", operation, molecules)
            aa_seed = boot_seed("ratio", "aa", operation, molecules)
            median, lo, hi = bootstrap_median_ci(bag["ab"], BOOT_DRAWS, seed)
            aa_median, aa_lo, aa_hi = bootstrap_median_ci(bag["aa"], BOOT_DRAWS, aa_seed)
            agree = bag["agree"] and bag["aa_agree"]
            call = "correctness" if not agree else outside(lo, hi, aa_lo, aa_hi)
            summaries.append({
                "row_kind": "ratio",
                "campaign": "ab",
                "operation": operation,
                "molecules": str(molecules),
                "code": "fastatomstruct",
                "version": "2.3.0/2.2.0",
                "arm": "fas230/fas220",
                "n": str(ROUNDS),
                "median_ms": "",
                "ci_lo_ms": "",
                "ci_hi_ms": "",
                "median_ratio": f"{median:.6f}",
                "ratio_lo": f"{lo:.6f}",
                "ratio_hi": f"{hi:.6f}",
                "aa_median_ratio": f"{aa_median:.6f}",
                "aa_lo": f"{aa_lo:.6f}",
                "aa_hi": f"{aa_hi:.6f}",
                "call": call,
                "agree": "1" if agree else "0",
                "detail": bag["detail"],
                "boot_draws": str(BOOT_DRAWS),
                "boot_seed": str(seed),
                "code_revision": bag["revision"],
                "host": bag["host"],
                "build": bag["build"],
            })
    return summaries


def write_summary(records, path, full):
    rows = assemble(records, full)
    with path.open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=SUMMARY_FIELDS, lineterminator="\n")
        writer.writeheader()
        writer.writerows(rows)
    return rows


def _fake_row(job, position, arm, time_ms, agree="1", max_abs="0"):
    return {
        "campaign": job["campaign"],
        "seed": str(SEED),
        "round": str(job["round"]),
        "position": str(position),
        "arm_order": ",".join(job["arm_order"]),
        "arm": arm,
        "code": "d-SEAMS" if arm == "seams" else "fastatomstruct",
        "version": "2.14.0" if arm == "seams" else ("2.3.0" if arm == "fas230" else "2.2.0"),
        "upstream_revision": "rev",
        "operation": job["operation"],
        "molecules": str(job["molecules"]),
        "repetitions": str(REPETITIONS[job["molecules"]]),
        "time_ms": f"{time_ms:.6f}",
        "ring3": "0",
        "ring4": "0",
        "ring5": "0",
        "ring6": str(2 * job["molecules"]),
        "q6_mean": "0.6285393576845176",
        "q6_min": "0.6285393576766468",
        "q6_max": "0.6285393576900145",
        "q6_sha256": "abc",
        "agree": agree,
        "max_abs": max_abs,
        "cpus": "0",
        "threads": "1",
        "periodic": "1",
        "warmup": "0",
        "code_revision": "461b129842bc9153ff8c4559e64555c73388b874",
        "host": "Intel Xeon Processor family 6 model 207 stepping 2",
        "build": "test",
    }


def self_check():
    order = base_order(ARMS_AB, SEED)
    if sorted(order) != sorted(ARMS_AB) or len(set(order)) != 3:
        raise SystemExit("shuffle dropped an arm")
    if base_order(ARMS_AB, SEED) != order:
        raise SystemExit("shuffle is not repeatable")
    if round_order(order, 0) != order or round_order(order, 1) != tuple(reversed(order)):
        raise SystemExit("round flip is wrong")
    if round_order(order, 2) != order:
        raise SystemExit("even rounds do not restore the shuffled order")
    aa = base_order(ARMS_AA, SEED)
    if sorted(aa) != sorted(ARMS_AA):
        raise SystemExit("A/A shuffle dropped an arm")
    constant, lo, hi = bootstrap_median_ci([4.0] * 20, 2000, 3)
    if constant != 4.0 or lo != 4.0 or hi != 4.0:
        raise SystemExit("constant bootstrap interval moved")
    again = bootstrap_median_ci([1.0, 2.0, 3.0, 4.0], 1000, 7)
    if again != bootstrap_median_ci([1.0, 2.0, 3.0, 4.0], 1000, 7):
        raise SystemExit("bootstrap is not repeatable")
    if statistics.median([1.0, 2.0, 3.0, 4.0]) != 2.5:
        raise SystemExit("even median is not the mean of the two central samples")
    if outside(0.5, 0.6, 0.9, 1.1) != "faster":
        raise SystemExit("faster band")
    if outside(1.2, 1.4, 0.9, 1.1) != "slower":
        raise SystemExit("slower band")
    if outside(0.95, 1.05, 0.9, 1.1) != "within":
        raise SystemExit("overlapping band")
    if outside(0.9, 1.0, 0.9, 1.1) != "within":
        raise SystemExit("a touching interval is not entirely outside the band")

    times = {"seams": 10.0, "fas220": 100.0, "fas230": 50.0, "fas220a": 100.0, "fas220b": 100.0}
    records = []
    for job in iter_jobs():
        if job["molecules"] != 1000 or job["operation"] != "q6":
            continue
        for position, arm in enumerate(job["arm_order"]):
            records.append(_fake_row(job, position, arm, times[arm]))
    rows = assemble(records, full=False)
    ratio = [row for row in rows if row["row_kind"] == "ratio"]
    if len(ratio) != 1 or ratio[0]["call"] != "faster" or ratio[0]["median_ratio"] != "0.500000":
        raise SystemExit(f"half-speed ratio was not faster: {ratio}")
    if ratio[0]["code_revision"] == "" or ratio[0]["host"] == "":
        raise SystemExit("summary row dropped code_revision or host")
    broken = []
    for row in records:
        copied = dict(row)
        copied["agree"] = "0"
        copied["max_abs"] = "1.0e-03"
        broken.append(copied)
    bad = [row for row in assemble(broken, full=False) if row["row_kind"] == "ratio"]
    if len(bad) != 1 or bad[0]["call"] != "correctness":
        raise SystemExit("a q6 mismatch was reported as a timing call")
    print("fas230 summary self-check ok")


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--samples", default="")
    parser.add_argument("--summary", default="")
    parser.add_argument("--self-check", action="store_true")
    args = parser.parse_args()
    if args.self_check:
        self_check()
        return 0
    if not args.samples or not args.summary:
        print("usage: h2h_fas230_summary.py --samples S --summary O", file=sys.stderr)
        return 2
    from pathlib import Path
    records = list(csv.DictReader(Path(args.samples).open(newline="")))
    write_summary(records, Path(args.summary), full=True)
    print(f"wrote {args.summary}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
