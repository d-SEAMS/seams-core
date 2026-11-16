"""Committor driver on a real 512-molecule mW cell.

Short LAMMPS runs exercise every stage the 4096-molecule shards use:
the climb with its restart, the climb table, both kinds of extension,
basin A from the climb's own plateau beside the fresh-liquid reference,
and shots that stop at the first basin they reach. Skipped when LMP is unset. NUCLEATION_PATHWAY, or
bbdir/nucleation-pathway, is the census.
"""

from __future__ import annotations

import csv
import os
import subprocess
import sys
import time
from dataclasses import replace
from pathlib import Path

import pytest

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "scripts"))
import committor_shoot as cs  # noqa: E402
import lmp_launch  # noqa: E402

SEED = 20261111
TINY = cs.Plan(
    n_cell=8,
    melt_steps=1000,
    liquid_steps=3000,
    dump_steps=1000,
    shot_steps=3000,
    shots=2,
    extend_steps=2000,
    limit_steps=8000,
    equil_ps=0.0,
)

pytestmark = pytest.mark.skipif(not os.environ.get("LMP"), reason="LMP is unset")


@pytest.fixture(scope="module")
def tools() -> tuple[list[str], str]:
    pathway = os.environ.get("NUCLEATION_PATHWAY") or str(ROOT / "bbdir" / "nucleation-pathway")
    if not Path(pathway).is_file():
        pytest.skip(f"no census at {pathway}")
    return lmp_launch.launcher(), pathway


@pytest.fixture(scope="module")
def climb(tmp_path_factory, tools) -> Path:
    _lmp, pathway = tools
    frames = tmp_path_factory.mktemp("climb")
    subprocess.run(
        [sys.executable, str(ROOT / "input" / "nucleation" / "regenerate.py"),
         "--n", "8", "--seeds", str(SEED), "--melt", "1000", "--quench", "4000",
         "--dump-every", "1000", "--temperature", "180", "--seams", pathway,
         "--dest", str(frames)],
        check=True, cwd=ROOT,
    )
    return frames


def climb_table(path: Path) -> dict[str, str]:
    rows = cs.read_climb_outcomes(path)
    assert list(rows) == [SEED]
    return rows[SEED]


def test_climb_writes_restart_and_three_series(climb):
    directory = cs.seed_dir(climb, SEED, TINY)
    assert (directory / "quench.restart").is_file()
    rows = cs.read_frames(directory / "frames.tsv")
    assert len(rows) == 5
    assert cs.has_series(rows)
    for row in rows:
        assert 0 <= int(row["n_chill_core"]) <= int(row["n_chill_ext"]) <= 512
        assert 0 <= int(row["n_q6solid"]) <= 512


def test_climb_table_and_both_extensions(climb, tools, tmp_path):
    lmp, pathway = tools
    directory = cs.seed_dir(climb, SEED, TINY)
    table = tmp_path / "climb-0-of-1.csv"
    cs.write_climb(table, str(SEED), "0/1", climb, "test", pathway, plan=TINY)
    header = next(l for l in table.read_text().splitlines() if not l.startswith("#"))
    assert header.split(",") == cs.CLIMB_FIELDS
    row = climb_table(table)
    assert row["outcome"] == "no-nucleation"
    assert row["t_end_ps"] == "40" and row["frames"] == "5" and row["segments"] == "1"
    core = cs.column(cs.read_frames(directory / "frames.tsv"), "n_chill_core")
    assert float(row["n_chill_core_max"]) == max(core)
    assert row["n_q6solid_max"] != "" and row["n_chill_ext_max"] != ""
    times = cs.frame_times(cs.read_frames(directory / "frames.tsv"))
    want = cs.basin_a_row(core, times, TINY.plateau_mark, TINY.plateau_ps, TINY.equil_ps)
    assert {key: row[key] for key in want} == want
    assert row["basin_b_ps"] == "" and row["liquid_2ns_ref"] == ""

    assert cs.extend_seed(SEED, directory, lmp, pathway, 1, TINY) == "restart 2000 steps"
    assert (directory / "extend-1.restart").is_file()
    (directory / "extend-1.restart").unlink()
    result = cs.extend_seed(SEED, directory, lmp, pathway, 1, TINY)
    assert result == "dump-mb 2000 steps"
    assert cs.extend_seed(SEED, directory, lmp, pathway, 1, TINY) == "at the limit"

    rows = cs.read_frames(directory / "frames.tsv")
    steps = [int(r["timestep"]) for r in rows]
    assert steps == list(range(steps[0], steps[0] + 9000, 1000))
    assert [r["frame"] for r in rows] == [str(i) for i in range(1, 10)]
    cs.write_climb(table, str(SEED), "0/1", climb, "test", pathway, plan=TINY)
    row = climb_table(table)
    velocity = cs.derived_seed(SEED, "extend", 2)
    assert row["extension"] == f"restart;dump-mb:{velocity}"
    assert row["segments"] == "3" and row["t_end_ps"] == "80" and row["frames"] == "9"


def test_shots_stop_when_resolved(climb, tools, tmp_path):
    lmp, pathway = tools
    directory = cs.seed_dir(climb, SEED, TINY)
    data = tmp_path / "shot.data"
    cs.shooting_frame(directory, int(cs.read_frames(directory / "frames.tsv")[2]["timestep"]),
                      data)
    long = replace(TINY, shot_steps=100_000)
    started = time.perf_counter()
    quick = cs.run_shot(lmp, tmp_path / "quick", data, 7, 1.0e9, pathway, 1, long)
    assert quick["fate"] == "uncommitted" and quick["t_end_ps"] == "10"
    assert 0.0 < float(quick["wall_s"]) < 120
    assert "Total wall time" not in (tmp_path / "quick" / "log.lammps").read_text()
    assert not list((tmp_path / "quick").glob("frame.*.lammpstrj"))
    assert time.perf_counter() - started < 120

    held = replace(long, basin_b=0.0, basin_b_hold_ps=20.0)
    committed = cs.run_shot(lmp, tmp_path / "held", data, 7, -1.0, pathway, 1, held)
    assert committed["fate"] == "committed" and committed["t_end_ps"] == "30"
    assert "Total wall time" not in (tmp_path / "held" / "log.lammps").read_text()

    never = replace(TINY, basin_b=1.0e9)
    open_window = cs.run_shot(lmp, tmp_path / "open", data, 7, -1.0, pathway, 1, never)
    assert open_window["fate"] == "unresolved" and open_window["t_end_ps"] == "30"
    trace = (tmp_path / "open" / "trace.tsv").read_text().splitlines()
    assert trace[0].split("\t") == ["t_ps", *cs.SERIES]
    assert [line.split("\t")[0] for line in trace[1:]] == ["10", "20", "30"]


def test_shard_end_to_end(climb, tools, tmp_path):
    lmp, pathway = tools
    directory = cs.seed_dir(climb, SEED, TINY)
    rows = cs.read_frames(directory / "frames.tsv")
    core = cs.column(rows, "n_chill_core")
    top = max(core)
    first = core.index(top)
    if top < 1 or first == 0:
        pytest.skip(f"the tiny climb never rose: {core}")
    # Cuts scaled to the tiny cell: its first visit to its own maximum,
    # with no hold, since the tiny climb is only 40 ps long.
    plan = replace(TINY, cut=top - 1, hold_ps=0.0, sizes=(int(top),), plateau_mark=top)
    times = cs.frame_times(rows)
    assert cs.nucleation_onset(core, plan.cut, TINY.hold_ps, times) is None
    onset = cs.nucleation_onset(core, plan.cut, plan.hold_ps, times)
    assert onset == first

    out = tmp_path / "shot-0-of-1.csv"
    cs.write_shots(out, str(SEED), "0/1", climb, "test", lmp, pathway, 1, plan)
    body = [l for l in out.read_text().splitlines() if not l.startswith("#")]
    shots = list(csv.DictReader(body))
    assert [s["shot"] for s in shots] == ["0", "1"]
    assert body[0].split(",") == cs.SHOT_FIELDS
    for shot in shots:
        assert shot["fate"] in ("committed", "uncommitted", "unresolved")
        assert shot["shoot_n"] == str(int(top))
        assert 10 <= float(shot["t_end_ps"]) <= 30
    base = list(csv.DictReader(
        l for l in (tmp_path / "basins-0-of-1.csv").read_text().splitlines()
        if not l.startswith("#")))
    assert base[0]["onset_ps"] == cs.fmt(times[onset])
    assert base[0]["basin_a"] == cs.fmt(cs.core_p95(core[:first]))
    assert base[0]["first_passage_ps"] == cs.fmt(times[first])
    assert base[0]["basin_b"] == cs.fmt(plan.basin_b) and base[0]["liquid_2ns_ref"] != ""
    assert base[0]["liquid_seed"] == str(cs.derived_seed(SEED, "liquid"))
    assert (directory / "liquid-2ns" / "frames.tsv").is_file()
    points = list(csv.DictReader(
        l for l in (tmp_path / "points-0-of-1.csv").read_text().splitlines()
        if not l.startswith("#")))
    assert points[0]["timestep"] == rows[onset]["timestep"]
    assert (directory / "shots" / f"n{int(top)}-step{rows[onset]['timestep']}" / "shot.data").is_file()
    assert not out.with_name(out.name + ".partial").exists()
    assert cs.merge_shot_dir(tmp_path) == [(float(int(top)), s["fate"]) for s in shots]
