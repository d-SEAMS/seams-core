# Copyright (c) d-SEAMS developers
"""The head-to-head campaign harness keeps its protocol.

Checks the statistics, the interleaved schedule, the append-only record and
the provenance on every row, the path scrubber, the output reduction the
agreement check compares, the asv cell bookkeeping of the lanes, and the
six-round floor. No tool under test is imported and no sample is timed.
"""

import json
import math
import os
import re
import subprocess
import sys
import tempfile
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "scripts"))

try:
    import h2h_campaign as hc
except ImportError as exc:
    print(f"skip: {exc}")
    raise SystemExit(77)

RESULT_COLUMNS = [
    "result", "params", "version", "started_at", "duration", "stats_ci_99_a",
    "stats_ci_99_b", "stats_q_25", "stats_q_75", "stats_number", "stats_repeat",
    "samples", "profile",
]


def sample(workload, system, arm, r, wall):
    return {
        "type": "sample", "status": "ok", "workload": workload, "system": system,
        "arm": arm, "round": r, "wall_ms": wall, "cpu_ms": wall,
        "rss_peak_kb": 1000, "rss_before_kb": 900, "molecules": 512, "os_threads": 1,
    }


def check_statistics():
    s = hc.stats([4.0, 1.0, 3.0, 2.0])
    assert s["n"] == 4 and s["median"] == 2.5, s
    assert abs(s["q1"] - 1.75) < 1e-12 and abs(s["q3"] - 3.25) < 1e-12, s
    assert hc.stats([7.0])["median"] == 7.0
    assert hc.stats([])["n"] == 0
    band = hc.ratio_band([1.1, 0.9, 1.0])
    assert (band["min"], band["median"], band["max"]) == (0.9, 1.0, 1.1), band
    assert hc.ratio_band([])["n"] == 0


def check_verdicts():
    seams = [10.0, 10.2, 9.9, 10.1, 10.0, 10.05]
    twin = [10.1, 10.0, 10.0, 10.2, 9.9, 10.1]
    arms = {
        "seams": seams,
        "seams_aa": twin,
        "slow": [20.0 + 0.1 * i for i in range(6)],
        "fast": [5.0 + 0.1 * i for i in range(6)],
        "level": [10.0, 10.3, 9.8, 10.1, 10.15, 9.95],
    }
    rows = [sample("cutoff", "x", arm, r, wall)
            for arm, walls in arms.items() for r, wall in enumerate(walls)]
    got = {c["arm"]: c.get("verdict") for c in hc.cell_summary(rows)}
    assert got["slow"] == "seams faster", got
    assert got["fast"] == "competitor faster", got
    assert got["level"] == "within A/A band", got
    assert got["seams"] is None and got["seams_aa"] is None, got
    text = hc.render_summary(rows)
    assert "seams faster" in text and "competitor faster" in text


def check_schedule():
    cells = [("cutoff", "a", ("seams", "seams_aa", "x")), ("knn", "a", ("seams", "seams_aa"))]
    plan = list(hc.schedule(cells, 4))
    rounds = {r: [(w, s, a) for rr, _p, _o, w, s, a in plan if rr == r] for r in range(4)}
    assert rounds[0] == rounds[2], rounds
    assert rounds[1] == list(reversed(rounds[0])), rounds
    assert rounds[3] == list(reversed(rounds[2])), rounds
    for r, seq in rounds.items():
        assert len(seq) == len(set(seq)) == 5, (r, seq)
    orders = {(rr % 2, o) for rr, _p, o, *_ in plan}
    assert orders == {(0, "fwd"), (1, "rev")}, orders


def check_record():
    prov = hc.provenance("run-a", "2026-01-02")
    for key in ("host", "cpu_model", "lscpu", "kernel", "code_revision", "date", "run_id",
                "schema"):
        assert prov[key], key
    assert re.fullmatch(r"bench-[0-9a-f]{12}", prov["host"]), prov["host"]
    with tempfile.TemporaryDirectory(prefix="dseams-h2h-") as work:
        path = Path(work) / "record.jsonl"
        hc.Record(path, prov).append({"type": "first", "value": 1})
        before = path.read_text(encoding="utf-8")
        later = hc.Record(path, prov)
        later.append({"type": "second", "value": 2})
        assert path.read_text(encoding="utf-8").startswith(before)
        rows = hc.read_rows(path)
        assert [r["type"] for r in rows] == ["first", "second"]
        assert all(all(r[k] == v for k, v in prov.items()) for r in rows)
        try:
            later.append({"type": "bad", "value": float("nan")})
        except ValueError:
            pass
        else:
            raise AssertionError("a NaN row was accepted")
        assert len(hc.read_rows(path)) == 2
        assert not hc.has_plan(path, "run-a")
        later.append({"type": "plan"})
        assert hc.has_plan(path, "run-a") and not hc.has_plan(path, "run-b")


def check_scrub():
    tmp = tempfile.gettempdir()
    home = os.path.expanduser("~")
    text = f"failed in {ROOT}/src/a.cpp, {tmp}/work/frame.npz, {home}/cache/x"
    clean = hc.scrub(text)
    assert str(ROOT) not in clean and "<repo>/src/a.cpp" in clean, clean
    assert f"{tmp}/" not in clean and f"{home}/cache" not in clean, clean
    assert clean.startswith("failed in"), clean


def check_canonical():
    pairs = hc.canonical("cutoff", ("pairs", [0, 1, 2], [1, 2, 0]), 3, "freud")[0]
    rows = hc.canonical("cutoff", ("rows", [[1, 2], [0, 2], [0, 1]]), 3, "mdapy")[0]
    csr = hc.canonical("cutoff", ("csr", [0, 2, 4, 6], [1, 2, 0, 2, 0, 1]), 3, "seams")[0]
    for form in (pairs, rows, csr):
        assert form["n_pairs"] == 3 and form["degree_min"] == 2, form
    assert pairs["pairs_digest"] == rows["pairs_digest"] == csr["pairs_digest"]
    assert pairs["degree_digest"] == rows["degree_digest"] == csr["degree_digest"]
    other = hc.canonical("cutoff", ("pairs", [0, 1], [1, 2]), 3, "freud")[0]
    assert other["pairs_digest"] != pairs["pairs_digest"]
    hist = ("hist", {3: 2, 4: 1, 7: 5})
    seams = hc.canonical("rings", hist, 8, "seams")[0]
    ovito = hc.canonical("rings", hist, 8, "ovito")[0]
    assert seams["hist"] == {"3": 2, "4": 1, "5": 0, "6": 0} and seams["total"] == 3, seams
    assert ovito["hist"] == {"4": 1, "5": 0, "6": 0} and ovito["total"] == 1, ovito


def check_arms():
    assert hc.applicable("seams", "cutoff", 110592)
    assert not hc.applicable("vesin", "knn", 512)
    assert hc.applicable("cycless", "rings", 512)
    assert not hc.applicable("cycless", "rings", 4096)
    assert hc.REFERENCE in hc.ARMS["cutoff"] and hc.TWIN in hc.ARMS["rings"]
    for workload, arms in hc.ARMS.items():
        assert arms[:2] == (hc.REFERENCE, hc.TWIN), (workload, arms)
    assert hc.MIN_ROUNDS >= 6 and hc.ROUNDS >= hc.MIN_ROUNDS


def check_lane_cells():
    found = {
        "mod.C.time_x": {
            "type": "time", "param_names": ["operation", "threads"],
            "params": [["'a'", "'b'"], ["1", "2", "4"]], "version": "v1",
        },
        "mod.D.peakmem_y": {"type": "peakmemory", "param_names": [], "params": [],
                            "version": "v2"},
    }
    lane = hc.Lane("scaled", "scaling", "mod.C.time_x", cpus="threads")
    cells = hc.lane_cells(lane, found)
    assert len(cells) == 6
    assert cells[0]["label"] == "mod.C.time_x('a', 1)" and cells[0]["index"] == 0
    assert cells[5]["label"] == "mod.C.time_x('b', 4)" and cells[5]["index"] == 5
    assert cells[4]["params"] == {"operation": "'b'", "threads": "2"}
    pattern = hc.cell_regex(cells[1])
    assert re.search(pattern, cells[1]["label"])
    assert not re.search(pattern, cells[1]["label"] + "x")
    assert not re.search(pattern, cells[2]["label"])
    only = hc.Lane("one", "scaling", "mod.C.time_x", cpus="threads",
                   only=(("operation", ("'b'",)),))
    assert [c["index"] for c in hc.lane_cells(only, found)] == [3, 4, 5]
    peak = hc.lane_cells(hc.Lane("peak", "kernel", "mod.D.peakmem_y"), found)
    assert len(peak) == 1 and peak[0]["type"] == "peakmem" and peak[0]["label"] == "mod.D.peakmem_y"
    eight = list(range(8))
    assert hc.cell_cpus(lane, cells[0], "3", eight) == "0"
    assert hc.cell_cpus(lane, cells[2], "3", eight) == "0-3"
    assert hc.cell_cpus(lane, cells[2], "3", [0, 1]) is None
    assert hc.cell_cpus(hc.Lane("o", "kernel", "b"), cells[0], "3", eight) == "3"
    assert hc.cell_cpus(hc.Lane("a", "kernel", "b", cpus="all"), cells[0], "3", eight) == "0-7"
    assert hc.cpu_list([2, 4, 5]) == "2,4,5" and hc.cpu_list([6]) == "6"
    assert not hc.h2h_cell_runs({"workload": "'knn'", "system": "'ice512'", "arm": "'vesin'"})
    assert hc.h2h_cell_runs({"workload": "'cutoff'", "system": "'ice512'", "arm": "'vesin'"})
    assert hc.h2h_cell_runs({"system": "'ice512'", "arm": "'rsm'"})
    assert not hc.h2h_cell_runs({"system": "'ice4096'", "arm": "'cycless'"})
    try:
        hc.select_lanes("no-such-lane")
    except SystemExit:
        pass
    else:
        raise AssertionError("an unknown lane was accepted")
    names = {lane.name for lane in hc.select_lanes("kernel,rings,scaling")}
    assert "omp-scaling" in names and "h2h-rings-time" not in names
    assert {lane.name for lane in hc.select_lanes("h2h")} >= {"h2h-neighbours-time"}


def write_results(directory, benchmark, entry):
    host = directory / "h"
    host.mkdir(parents=True, exist_ok=True)
    (host / "abc-existing-python.json").write_text(
        json.dumps({"result_columns": RESULT_COLUMNS, "results": {benchmark: entry}}),
        encoding="utf-8")


def check_harvest():
    with tempfile.TemporaryDirectory(prefix="dseams-h2h-") as work:
        results = Path(work)
        flat = {"benchmark": "mod.C.time_x", "index": 0, "params": {}}
        write_results(results, "mod.C.time_x",
                      [[0.001], [], "v", 1, 0.1, [0], [0], [0.001], [0.001], [1], [1],
                       [[0.001]], None])
        assert hc.harvest(results, "h", flat) == {"value": 0.001, "samples": 1}
        assert hc.harvest(results, "h", {**flat, "benchmark": "mod.other"}) is None
        hc.clear_results(results, "h")
        assert hc.harvest(results, "h", flat) is None
        grid = {"benchmark": "mod.C.time_x", "params": {"threads": "4"}}
        write_results(results, "mod.C.time_x",
                      [[None, None, 0.5], [["1", "2", "4"]], "v", 1, 0.1, [], [], [], [],
                       [], [], [None, None, [0.5]], None])
        assert hc.harvest(results, "h", {**grid, "index": 2}) == {"value": 0.5, "samples": 1}
        assert hc.harvest(results, "h", {**grid, "index": 0}) is None
        peak = {"benchmark": "mod.D.peakmem_y", "index": 0, "params": {}}
        write_results(results, "mod.D.peakmem_y",
                      [[1234567], [], "w", 1, 0.1, None, None, None, None, None, None,
                       None, None])
        assert hc.harvest(results, "h", peak) == {"value": 1234567.0, "samples": 1}
        write_results(results, "mod.D.peakmem_y",
                      [[float("nan")], [], "w", 1, 0.1, None, None, None, None, None, None,
                       None, None])
        assert hc.harvest(results, "h", peak) is None


def check_lane_summary():
    def asv(lane, arm, r, value, label="mod.C.time_x(4)", **extra):
        return {"type": "asv", "lane": lane, "group": "scaling", "cell": label, "kind": "time",
                "params": {"threads": label[label.index("(") + 1:-1]}, "arm": arm,
                "round": r, "status": "ok", "value": value, **extra}
    rows = []
    for r in range(6):
        for threads, base in (("1", 0.08), ("4", 0.02)):
            label = f"mod.C.time_x({threads})"
            rows.append(asv("scaled", "A", r, base, label))
            rows.append(asv("scaled", "A2", r, base * 1.01, label))
    stats = {c["cell"]: c for c in hc.lane_stat_rows(rows)}
    one = stats["mod.C.time_x(1)"]
    assert one["A"]["n"] == 6 and math.isclose(one["A"]["median"], 80.0), one
    assert math.isclose(one["aa"]["median"], 1.01), one
    text = hc.render_lanes(rows)
    assert "speed-up over one thread" in text and "4.00" in text, text


def check_cli_floor_and_plan():
    script = ROOT / "scripts" / "h2h_campaign.py"
    with tempfile.TemporaryDirectory(prefix="dseams-h2h-") as work:
        record = str(Path(work) / "record.jsonl")
        short = subprocess.run(
            [sys.executable, str(script), "run", "--head-to-head", "rings", "--rounds", "3",
             "--dry-run", "--record", record],
            capture_output=True, text=True, check=False)
        assert short.returncode != 0 and "at least 6 rounds" in short.stderr, short.stderr
        plan = subprocess.run(
            [sys.executable, str(script), "run", "--head-to-head", "neighbours", "--systems",
             "ice512", "sII", "--dry-run", "--record", record],
            capture_output=True, text=True, check=False)
        assert plan.returncode == 0, plan.stderr
        expected = 2 * hc.ROUNDS * (len(hc.ARMS["cutoff"]) + len(hc.ARMS["knn"]))
        assert f"samples {expected}" in plan.stdout, plan.stdout
        assert not Path(record).exists(), "a dry run wrote a record"


def main() -> int:
    checks = [(name, fn) for name, fn in sorted(globals().items())
              if name.startswith("check_") and callable(fn)]
    for name, fn in checks:
        fn()
        print(f"{name} ok")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
