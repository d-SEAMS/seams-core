# Copyright (c) d-SEAMS developers
"""The head-to-head campaign harness keeps its protocol.

Checks the statistics, the interleaved schedule, the append-only record and
the provenance on every row, the path scrubber, the output reduction the
agreement check compares, the asv cell bookkeeping and filtered parameter
grid of the lanes, the twelve-round floor, required CLI inputs, thread-count
consistency, and the ancestry/tree/report integrity of any committed records.
No tool under test is imported and no sample is timed.
"""

import json
import getpass
import ipaddress
import math
import os
import platform
import re
import subprocess
import sys
import tempfile
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "scripts"))
sys.path.insert(0, str(ROOT / "benchmarks"))

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

ALLOWED_EMAILS: frozenset[str] = frozenset()
ALLOWED_PUBLIC_HOSTS = frozenset({
    "dseams.github.io", "github.com", "pypi.org", "readthedocs.io", "zenodo.org",
})
FILE_SUFFIXES = frozenset({
    "cpp", "jsonl", "lammpstrj", "md", "npy", "py", "so",
})


def sample(workload, system, arm, r, wall, os_threads=1, core="0"):
    return {
        "type": "sample", "status": "ok", "workload": workload, "system": system,
        "arm": arm, "tool": hc.TOOL_OF.get(arm, arm), "round": r, "wall_ms": wall,
        "cpu_ms": wall, "rss_peak_kb": 1000, "rss_before_kb": 900, "molecules": 512,
        "os_threads": os_threads, "single_os_thread": os_threads == 1, "core": core,
        "calls_per_sample": 1, "warmup_calls": 0,
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
    rows.extend({
        "type": "agreement", "workload": "cutoff", "system": "x", "arm": arm,
        "reference": "seams", "equal": True, "matches_brute_force": True,
    } for arm in arms)
    got = {c["arm"]: c.get("verdict") for c in hc.cell_summary(rows)}
    assert got["slow"] == "seams faster", got
    assert got["fast"] == "competitor faster", got
    assert got["level"] == "within A/A band", got
    assert got["seams"] is None and got["seams_aa"] is None, got
    text = hc.render_summary(rows)
    assert "seams faster" in text and "competitor faster" in text


def check_disagreement_blocks_speed_verdicts():
    rows = [
        sample("knn", "x", "seams", r, 10.0) for r in range(6)
    ] + [
        sample("knn", "x", "seams_aa", r, 10.0) for r in range(6)
    ] + [
        sample("knn", "x", "wrong", r, 20.0) for r in range(6)
    ] + [
        {
            "type": "agreement", "workload": "knn", "system": "x",
            "arm": "seams", "reference": "seams", "equal": True,
            "matches_brute_force": True,
        },
        {
            "type": "agreement", "workload": "knn", "system": "x",
            "arm": "wrong", "reference": "seams", "equal": False,
            "matches_brute_force": False,
        },
    ]
    wrong = next(c for c in hc.cell_summary(rows) if c["arm"] == "wrong")
    assert wrong["verdict"].startswith("not comparable:"), wrong
    assert "faster" not in wrong["verdict"] and "slower" not in wrong["verdict"], wrong
    text = hc.render_summary(rows)
    line = next(line for line in text.splitlines() if "| wrong |" in line)
    assert "| - |" in line and "not comparable:" in line, line


def check_disagreement_classification():
    record = ROOT / "data" / "h2h-campaign" / "2026-11-18-r2.jsonl"
    rows = hc.read_rows(record)
    found = hc.disagreement_classification(rows)
    assert len(found) == 15, found
    assert sum(x["workload"] == "knn" for x in found) == 8, found
    assert sum(x["workload"] == "rings" for x in found) == 7, found
    assert {x["reason_class"] for x in found} <= {
        "competitor bug", "different definition/cutoff convention",
        "tolerance", "d-SEAMS bug",
    }, found
    stacking = [x for x in found if x["workload"] == "knn" and x["system"] == "stacking"]
    assert len(stacking) == 3, stacking
    assert all(x["matching_side"] == "neither (brute force)" for x in stacking), stacking
    assert all(x["reason_class"] == "tolerance" for x in stacking), stacking


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
    prov = hc.provenance("run-a", "2026-01-02", "0")
    for key in ("host", "core", "cpu_model", "lscpu", "kernel", "code_revision", "date",
                "run_id", "schema", "engine_tree"):
        assert prov[key], key
    assert re.fullmatch(r"bench-[0-9a-f]{12}", prov["host"]), prov["host"]
    assert re.fullmatch(r"[0-9a-f]{40}", prov["engine_tree"]), prov["engine_tree"]
    head_tree = subprocess.check_output(
        ["git", "rev-parse", f"{prov['code_revision']}^{{tree}}"],
        cwd=ROOT, text=True).strip()
    assert prov["engine_tree"] == head_tree, (prov["engine_tree"], head_tree)
    assert "base_revision" in prov
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


def check_machine_capture_without_lscpu():
    real_which = hc.shutil.which
    real_run = hc.subprocess.run
    try:
        hc.shutil.which = lambda _name: None
        summary = hc.machine_summary()
        assert summary["cpu_model"] and summary["lscpu"] and summary["lscpu_sha256"], summary
        assert "unavailable (" in summary["lscpu"], summary

        hc.shutil.which = lambda _name: "lscpu"

        def missing_lscpu(command, *args, **kwargs):
            if command == ["lscpu"]:
                raise FileNotFoundError("missing")
            return real_run(command, *args, **kwargs)

        hc.subprocess.run = missing_lscpu
        summary = hc.machine_summary()
        assert summary["cpu_model"] and "unavailable (" in summary["lscpu"], summary
    finally:
        hc.shutil.which = real_which
        hc.subprocess.run = real_run


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
    assert hc.MIN_ROUNDS >= 12 and hc.ROUNDS >= hc.MIN_ROUNDS
    for name, ncell in hc.ICE_CELLS.items():
        assert hc.SYSTEM_MOLECULES[name] == 8 * (ncell ** 3), name
    for name in (*hc.MIXED, *hc.DUMPS):
        pos, _ = hc.build_frame(name)
        assert hc.SYSTEM_MOLECULES[name] == int(pos.shape[0]), name
    for key in ("OMP_NUM_THREADS", "MKL_NUM_THREADS", "OPENBLAS_NUM_THREADS",
                "NUMBA_NUM_THREADS", "TBB_NUM_THREADS", "RAYON_NUM_THREADS"):
        assert hc.PINNED_ENV[key] == "1", key


def check_campaign_h2h_grid():
    import campaign_h2h as ch

    assert ch.NeighbourLists.param_names == ["cell"]
    assert ch.RingStatistics.param_names == ["cell"]
    neigh_cells = ch.NeighbourLists.params[0]
    assert not any(w == "knn" and a == "vesin" for w, _s, a in neigh_cells)
    assert any(w == "cutoff" and a == "vesin" for w, _s, a in neigh_cells)
    ring_cells = ch.RingStatistics.params[0]
    assert not any(a == "cycless" and hc.SYSTEM_MOLECULES[s] > 1000 for s, a in ring_cells)
    assert any(a == "cycless" and s == "ice512" for s, a in ring_cells)


def check_lane_cells():
    found = {
        "mod.C.time_x": {
            "type": "time", "param_names": ["operation", "threads"],
            "params": [["'a'", "'b'"], ["1", "2", "4"]], "version": "v1",
        },
        "mod.D.peakmem_y": {"type": "peakmemory", "param_names": [], "params": [],
                            "version": "v2"},
        "campaign_h2h.NeighbourLists.time_call": {
            "type": "time", "param_names": ["cell"],
            "params": [[
                "('cutoff', 'ice512', 'vesin')",
                "('cutoff', 'ice4096', 'vesin')",
                "('knn', 'ice512', 'seams')",
            ]],
            "version": "v3",
        },
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
    h2h_lane = hc.Lane("h2h-n", "h2h", "campaign_h2h.NeighbourLists.time_call", twin=False,
                       only=(("system", ("ice512",)),))
    h2h_cells = hc.lane_cells(h2h_lane, found)
    assert [c["index"] for c in h2h_cells] == [0, 2], h2h_cells
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
    assert not hc.h2h_cell_runs({"cell": "('knn', 'ice512', 'vesin')"})
    assert hc.h2h_cell_runs({"cell": "('cutoff', 'ice512', 'vesin')"})
    assert hc.h2h_cell_runs({"system": "'ice512'", "arm": "'rsm'"})
    assert not hc.h2h_cell_runs({"cell": "('ice4096', 'cycless')"})
    by_name = {item.name: item for item in hc.LANES}
    assert by_name["omp-session"].calls_per_sample == 32
    assert by_name["omp-session"].warmup_calls == 0
    assert by_name["cutoff-ice-neighbours"].calls_per_sample == 1
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
    def asv(lane, arm, r, value, label="mod.C.time_x(4)", calls=1, **extra):
        return {"type": "asv", "lane": lane, "group": "scaling", "cell": label, "kind": "time",
                "params": {"threads": label[label.index("(") + 1:-1]}, "arm": arm,
                "round": r, "status": "ok", "value": value, "calls_per_sample": calls, **extra}
    rows = []
    for r in range(6):
        for threads, base in (("1", 0.08), ("4", 0.02)):
            label = f"mod.C.time_x({threads})"
            rows.append(asv("scaled", "A", r, base, label, calls=32))
            rows.append(asv("scaled", "A2", r, base * 1.01, label, calls=32))
    stats = {c["cell"]: c for c in hc.lane_stat_rows(rows)}
    one = stats["mod.C.time_x(1)"]
    assert one["A"]["n"] == 6 and math.isclose(one["A"]["median"], 80.0), one
    assert math.isclose(one["aa"]["median"], 1.01), one
    assert one["calls_per_sample"] == 32
    text = hc.render_lanes(rows)
    assert "speed-up over one thread" in text and "4.00" in text, text
    assert "| scaled | threads=1 | 32 | 6 |" in text, text


def check_absolute_paths():
    class Args:
        record = "data/h2h-campaign/x.jsonl"
        raw = "bbdir-h2h/x/raw"
        work = ""
        marker = "bbdir-h2h/x/lanes.done"

    args = Args()
    hc.absolute_paths(args)
    assert Path(args.record).is_absolute() and args.record.endswith("data/h2h-campaign/x.jsonl")
    assert Path(args.raw).is_absolute() and Path(args.marker).is_absolute()
    assert args.work == "", "an empty option stays empty so that the default applies"


def check_lane_counters():
    rows = [
        {"type": "asv", "run_id": "a", "group": "kernel", "cell": "k(1)", "status": "ok"},
        {"type": "asv", "run_id": "a", "group": "kernel", "cell": "k(2)", "status": "failed"},
        {"type": "asv", "run_id": "a", "group": "rings", "cell": "r", "status": "timeout"},
        {"type": "asv", "run_id": "a", "group": "kernel", "cell": "k(3)", "status": "skipped"},
        {"type": "asv", "run_id": "a", "group": "inventory", "cell": "i", "status": "failed"},
        {"type": "asv", "run_id": "b", "group": "kernel", "cell": "k(1)", "status": "ok"},
    ]
    with tempfile.TemporaryDirectory(prefix="dseams-h2h-") as work:
        path = Path(work) / "record.jsonl"
        path.write_text("".join(json.dumps(r) + "\n" for r in rows), encoding="utf-8")
        assert hc.lane_samples(path, "a") == 1 and hc.lane_samples(path, "b") == 1
        broken = sorted(r["cell"] for r in hc.lane_failures(path, "a"))
        assert broken == ["k(2)", "r"], "the inventory and skipped cells are not failures"
        assert hc.lane_failures(path, "b") == []


def check_attempts():
    def plan(run_id, lanes):
        return {"type": "plan", "kind": "lanes", "run_id": run_id, "lanes": lanes}

    def asv(run_id, status, r=0, **extra):
        return {"type": "asv", "run_id": run_id, "lane": "scaled", "group": "scaling",
                "cell": "mod.C.time_x(1)", "kind": "time", "params": {}, "arm": "A",
                "round": r, "status": status, "code_revision": "0123456789abcdef", **extra}

    rows = [plan("first", ["scaled"]), plan("other", ["different"]), plan("second", ["scaled"])]
    rows += [asv("first", "failed", detail="no information stored about the machine"),
             asv("first", "skipped", r=1, detail="earlier round: failed"),
             asv("other", "failed"),
             asv("second", "ok", value=0.001), asv("second", "ok", r=1, value=0.002)]
    replaced = hc.replaced_runs(rows)
    assert replaced == {"first": "second"}, replaced
    current = [r for r in rows if r.get("run_id") not in replaced]
    earlier = [r for r in rows if r.get("run_id") in replaced]
    stats = hc.lane_stat_rows(current)
    assert [c["cell"] for c in stats if c["lane"] == "scaled"] == ["mod.C.time_x(1)"]
    assert stats[0]["A"]["n"] == 2, "the replaced attempt must not leak into the cell"
    text = hc.render_attempts(earlier, replaced)
    assert "| first | second | 01234567 |" in text and "asv/failed 1" in text, text
    assert "no information stored" in text, text

    later_empty = rows + [plan("third", ["scaled"]), asv("third", "failed")]
    assert hc.replaced_runs(later_empty)["third"] == "second", (
        "a later run that recorded no sample does not replace one that did")


def check_header_and_inventory():
    def asv(run_id, revision, group, threads, value):
        return {"type": "asv", "run_id": run_id, "lane": "scaled", "group": group,
                "cell": f"mod.C.time_x({threads})", "kind": "time", "arm": "A", "round": 0,
                "params": {"threads": str(threads)}, "status": "ok", "value": value,
                "core": "0", "calls_per_sample": 1,
                "code_revision": revision, "base_revision": "c" * 40, "engine_tree": "d" * 40,
                "host": "bench-test", "kernel": "6.0", "lscpu": "cpu"}

    rows = [asv("lanes", "a" * 40, "scaling", 1, 0.08), asv("lanes", "a" * 40, "scaling", 8, 0.02),
            asv("inventory", "b" * 40, "inventory", 1, 0.5),
            asv("inventory", "b" * 40, "inventory", 8, 0.1)]
    text = hc.render_summary(rows)
    assert f"- lanes: code_revision {'a' * 40}; base_revision {'c' * 40}; " in text, text
    assert f"- inventory: code_revision {'b' * 40}; " in text and f"engine_tree {'d' * 40}" in text
    assert "; core 0" in text
    assert f"- lanes: code_revision {'a' * 40}" in hc.render_summary(rows, "lanes")
    assert "b" * 40 not in hc.render_summary(rows, "lanes")
    speed = hc.render_lanes(rows)
    assert speed.count("| scaled | - | 1.00 |") == 1, "the inventory is not a speed-up row"


def check_thread_report_consistency():
    rows = [
        sample("cutoff", "ice512", "seams", 0, 1.0, os_threads=1),
        sample("cutoff", "ice512", "seams_aa", 0, 1.0, os_threads=1),
        sample("cutoff", "ice512", "vesin", 0, 1.5, os_threads=8),
        sample("cutoff", "ice512", "mdapy", 0, 1.1, os_threads=3),
    ]
    ts = hc.thread_summary(rows)
    assert ts["total"] == 4 and ts["single"] == 2 and ts["multi"] == 2, ts
    assert ts["by_tool"] == {"mdapy": [3], "seams": [1], "vesin": [8]}, ts
    text = hc.render_summary(rows)
    assert "2/4 good samples had 1 OS thread; 2/4 good samples had >1 OS thread." in text, text
    assert "mdapy (3), vesin (8)" in text, text
    assert "every sample is single-threaded" not in text.lower()


def check_failure_detail():
    skipped = "n/a [100.00%]\n asv: skipped: NotImplementedError('plumed is not importable')\n"
    assert hc.failure_detail(skipped) == "NotImplementedError('plumed is not importable')"
    compile_log = "\n".join(
        ["cc -c x.cpp", "x.cpp:9:3: error: no matching function for call to 'knn'",
         "note: candidate expects 2 arguments", "ninja: build stopped", "make: failed",
         "asv: benchmark failed (exit status 1)"])
    found = hc.failure_detail(compile_log)
    assert found.splitlines()[0].endswith("no matching function for call to 'knn'"), found
    assert found.endswith("asv: benchmark failed (exit status 1)"), found
    assert hc.failure_detail(f"error in {ROOT}/src/a.cpp").startswith("error in <repo>/src")


def check_cli_floor_and_plan():
    script = ROOT / "scripts" / "h2h_campaign.py"
    with tempfile.TemporaryDirectory(prefix="dseams-h2h-") as work:
        record = str(Path(work) / "record.jsonl")
        raw = str(Path(work) / "raw")
        wdir = str(Path(work) / "work")
        missing_core = subprocess.run(
            [sys.executable, str(script), "run", "--head-to-head", "rings",
             "--stamp", "2026-11-18-test", "--dry-run",
             "--record", record, "--raw", raw, "--work", wdir],
            capture_output=True, text=True, check=False)
        assert missing_core.returncode != 0, "missing --core must fail"
        missing_raw = subprocess.run(
            [sys.executable, str(script), "run", "--head-to-head", "rings", "--core", "0",
             "--stamp", "2026-11-18-test", "--dry-run",
             "--record", record, "--work", wdir],
            capture_output=True, text=True, check=False)
        assert missing_raw.returncode != 0, "missing --raw must fail"
        short = subprocess.run(
            [sys.executable, str(script), "run", "--head-to-head", "rings", "--rounds", "8",
             "--core", "0", "--stamp", "2026-11-18-test", "--dry-run",
             "--record", record, "--raw", raw, "--work", wdir],
            capture_output=True, text=True, check=False)
        assert short.returncode != 0 and "at least 12 rounds" in short.stderr, short.stderr
        plan = subprocess.run(
            [sys.executable, str(script), "run", "--head-to-head", "neighbours", "--systems",
             "ice512", "sII", "--core", "0", "--stamp", "2026-11-18-test", "--dry-run",
             "--record", record, "--raw", raw, "--work", wdir],
            capture_output=True, text=True, check=False)
        assert plan.returncode == 0, plan.stderr
        expected = 2 * hc.ROUNDS * (len(hc.ARMS["cutoff"]) + len(hc.ARMS["knn"]))
        assert f"samples {expected}" in plan.stdout, plan.stdout
        assert not Path(record).exists(), "a dry run wrote a record"


def campaign_leaks(text: str, current_user: str = "", node: str = "") -> list[str]:
    leaks = []
    machine_prefixes = tuple(
        p for p in (str(ROOT), os.path.expanduser("~"), tempfile.gettempdir())
        if p and p != "/"
    )
    for prefix in machine_prefixes:
        if prefix in text:
            leaks.append("current-machine path")
    user_path = re.compile(
        r"(?:/(?:home|Users)/[A-Za-z0-9_.-]+(?:/|$)|"
        r"[A-Za-z]:\\+Users\\+[A-Za-z0-9_.-]+(?:\\+|$))"
    )
    if user_path.search(text):
        leaks.append("user path")
    if current_user and len(current_user) > 2 and re.search(
            rf"(?<![A-Za-z0-9_.-]){re.escape(current_user)}(?![A-Za-z0-9_.-])", text):
        leaks.append("current username")
    if node and len(node) > 3 and node in text:
        leaks.append("current hostname")

    token_parts = (
        "ghp_", "github_pat_", "gho_", "hf_", "sk-", "AKIA", "xoxb-", "xoxp-",
    )
    if any(part in text for part in token_parts):
        leaks.append("token")
    assignment = re.compile(
        r"(?im)\b[A-Za-z_][A-Za-z0-9_.-]*\s*=\s*"
        r"(?:[0-9a-f]{32,}|[A-Za-z0-9+/]{40,}={0,2})\b"
    )
    if assignment.search(text):
        leaks.append("secret assignment")

    email_re = re.compile(r"\b[A-Za-z0-9._%+-]+@[A-Za-z0-9.-]+\.[A-Za-z]{2,}\b")
    if any(email.lower() not in ALLOWED_EMAILS for email in email_re.findall(text)):
        leaks.append("email")

    ip_re = re.compile(r"(?<![\d.])(?:\d{1,3}\.){3}\d{1,3}(?![\d.])")
    for candidate in ip_re.findall(text):
        try:
            address = ipaddress.ip_address(candidate)
        except ValueError:
            continue
        if address.is_private or address.is_loopback or address.is_link_local:
            leaks.append("private IP")
            break

    host_re = re.compile(
        r"(?<![@A-Za-z0-9_-])"
        r"(?:[a-z0-9](?:[a-z0-9-]{0,61}[a-z0-9])?\.)+[a-z]{2,63}"
        r"(?![A-Za-z0-9_-])"
    )
    for candidate in host_re.findall(text):
        host = candidate.lower()
        if host.rsplit(".", 1)[-1] in FILE_SUFFIXES:
            continue
        if host not in ALLOWED_PUBLIC_HOSTS:
            leaks.append("hostname")
            break
    return leaks


def check_campaign_leak_checker_negative_fixtures():
    fixtures = {
        "current username": "owner=" + getpass.getuser(),
        "user path": "/" + "home" + "/sample/private.txt",
        "windows user path": "C:" + "\\Users" + "\\sample\\private.txt",
        "token": "token=" + "ghp" + "_" + "a" * 24,
        "secret assignment": "credential=" + "a" * 48,
        "email": "contact=" + "person" + "@" + "example.invalid",
        "private IP": "address=" + ".".join(("10", "2", "3", "4")),
        "hostname": "machine=" + "build" + "." + "internal",
    }
    for kind, fixture in fixtures.items():
        assert campaign_leaks(
            fixture, current_user=getpass.getuser(), node=platform.node()
        ), (kind, fixture)


def check_campaign_records():
    data_dir = ROOT / "data" / "h2h-campaign"
    if not data_dir.is_dir():
        return
    records = sorted(data_dir.glob("*.jsonl"))
    if not records:
        return
    head = subprocess.check_output(["git", "rev-parse", "HEAD"], cwd=ROOT, text=True).strip()
    tree_cache: dict[str, str] = {}
    for record_path in records:
        rows = hc.read_rows(record_path)
        assert rows, f"{record_path} is empty"
        for idx, row in enumerate(rows, 1):
            assert row.get("code_dirty") is False, (record_path.name, idx, row.get("code_dirty"))
            rev = row.get("code_revision", "")
            assert rev, (record_path.name, idx)
            if rev not in tree_cache:
                anc = subprocess.run(
                    ["git", "merge-base", "--is-ancestor", rev, head],
                    cwd=ROOT, check=False)
                assert anc.returncode == 0, f"{record_path.name}:{idx} {rev} is not an ancestor of {head}"
                tree_cache[rev] = subprocess.check_output(
                    ["git", "rev-parse", f"{rev}^{{tree}}"], cwd=ROOT, text=True).strip()
            assert row.get("engine_tree") == tree_cache[rev], (
                record_path.name, idx, row.get("engine_tree"), tree_cache[rev])
            assert row.get("core") is not None and str(row["core"]).isdigit(), (record_path.name, idx)
            if row.get("type") in ("sample", "asv"):
                assert "calls_per_sample" in row and int(row["calls_per_sample"]) >= 1, (
                    record_path.name, idx)
                if row.get("type") == "sample" and row.get("status") == "ok":
                    if row.get("tool") == "seams":
                        assert row.get("os_threads") == 1, (record_path.name, idx, row)
            if row.get("type") == "plan" and (
                "workloads" in row
                or any(name in {lane.name for lane in hc.LANES}
                       for name in (row.get("lanes") or []))
            ):
                assert int(row.get("rounds", 0)) >= 12, (record_path.name, idx, row)
        assert not hc.replaced_runs(rows), f"{record_path.name} contains replaced rescue runs"
        md_path = record_path.with_name(
            "2026-11-18-r3.md"
            if record_path.name == "2026-11-18-r2.jsonl"
            else record_path.with_suffix(".md").name
        )
        assert md_path.is_file(), f"missing summary {md_path}"
        md_text = md_path.read_text(encoding="utf-8")
        expected_md = hc.render_summary(rows) + hc.render_lanes(rows)
        assert md_text == expected_md, f"{md_path.name} does not match {record_path.name}"
        ts = hc.thread_summary(rows)
        assert (f"{ts['single']}/{ts['total']} good samples had 1 OS thread; "
                f"{ts['multi']}/{ts['total']} good samples had >1 OS thread.") in md_text
    campaign_files = sorted(data_dir.glob("*.jsonl")) + sorted(data_dir.glob("*.md"))
    for path in campaign_files:
        leaks = campaign_leaks(
            path.read_text(encoding="utf-8"),
            current_user=getpass.getuser(),
            node=platform.node(),
        )
        assert not leaks, f"{path.name} contains sensitive data: {sorted(set(leaks))}"


def main() -> int:
    checks = [(name, fn) for name, fn in sorted(globals().items())
              if name.startswith("check_") and callable(fn)]
    for name, fn in checks:
        fn()
        print(f"{name} ok")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
