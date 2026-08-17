#!/usr/bin/env python3
"""Run the v2 pipeline on the v1 figshare demonstration trajectories.

The d-SEAMS 1.0 paper (10.1021/acs.jcim.0c00031) demonstrated on five
LAMMPS trajectories archived in the figshare project "d-SEAMS Datasets"
(https://figshare.com/projects/d-SEAMS_Datasets/73545). This script
fetches those exact deposits (MD5-verified against figshare's records)
and runs each corresponding analysis through require("dseams") in
yodaStruct.

Subcommands (all paths relative to the repository root):

  fetch <trajdir>                       download + verify the deposits
  demos <yoda_src> <yoda_build> <trajdir> <outdir>
                                        run the five deposits via require("dseams")

The Python-bindings demonstrations live as percent-format notebooks
under repro/notebooks/; jupytext converts them and papermill executes
them via the figshare_notebook Snakemake rule.
"""
import hashlib
import json
import os
import pathlib
import shutil
import subprocess
import sys
import time
import urllib.request

FILES = {
    "nucleation": {
        "name": "nucleation.lammpstrj",
        "url": "https://ndownloader.figshare.com/files/20418309",
        "md5": "c28fc209a2854d2ba49e8f9b55accbe8",
        "doi": "10.6084/m9.figshare.11448702.v1",
    },
    "chillplus": {
        "name": "mW_cubic.lammpstrj",
        "url": "https://ndownloader.figshare.com/files/20418330",
        "md5": "3fb29979d1548e0b3a8c9630e155a3d9",
        "doi": "10.6084/m9.figshare.11448720.v1",
    },
    "nanotube": {
        "name": "dump-240-square.lammpstrj",
        "url": "https://ndownloader.figshare.com/files/20418540",
        "md5": "9a33476aea679b54a8aefca280d06479",
        "doi": "10.6084/m9.figshare.11448768.v1",
    },
    "monolayer": {
        "name": "dump-6-320-310.lammpstrj",
        "url": "https://ndownloader.figshare.com/files/20418408",
        "md5": "09a4a595052f399ac94d583474e5ae89",
        "doi": "10.6084/m9.figshare.11448741.v1",
    },
    "rdf2d": {
        "name": "dump-320.lammpstrj",
        "url": "https://ndownloader.figshare.com/files/20418318",
        "md5": "71eb895b197600dc89ac64d5dc6b15a5",
        "doi": "10.6084/m9.figshare.11448711.v1",
    },
}

def example_lua_root():
    """Directory that contains example_lua/ (yodaStruct, not this tree)."""
    marker = pathlib.Path("example_lua/chillPlus/iceType/vars.lua")
    here = pathlib.Path.cwd()
    env = os.environ.get("YODASTRUCT_ROOT")
    candidates = []
    if env:
        candidates.append(pathlib.Path(env))
    candidates.extend(
        [
            here,
            here.parent / "yodaStruct",
            here.parent / "seams-base-repro",
        ]
    )
    for root in candidates:
        if (root / marker).is_file():
            return root
    raise FileNotFoundError(
        "example_lua not found; set YODASTRUCT_ROOT or clone "
        "https://github.com/d-SEAMS/yodaStruct next to this tree"
    )


# Each demo pairs a figshare trajectory with the 2.x Lua library
# (require("dseams")). frame="last" is the crystallized end of the
# nucleation run.
DEMOS = [
    {"key": "chillPlus", "traj": "chillplus", "frame": 1},
    {"key": "bulkTopologicalCriterion", "traj": "nucleation", "frame": "last"},
    {"key": "iceNanotube", "traj": "nanotube", "frame": 1},
    {"key": "monolayer", "traj": "monolayer", "frame": 1},
    {"key": "rdf2D", "traj": "rdf2d", "frame": 1},
]

SCRIPT = pathlib.Path(__file__).resolve().parent.parent / "lua" / "figshare_demo.lua"


def find_lua():
    for name in ("lua5.4", "lua-5.4", "lua5.3", "lua"):
        path = shutil.which(name)
        if path:
            return path
    raise FileNotFoundError("lua 5.3/5.4 not on PATH; add lua to the repro env")

def md5sum(path):
    h = hashlib.md5()
    with open(path, "rb") as f:
        for chunk in iter(lambda: f.read(1 << 20), b""):
            h.update(chunk)
    return h.hexdigest()


def fetch(trajdir):
    trajdir = pathlib.Path(trajdir)
    trajdir.mkdir(parents=True, exist_ok=True)
    for key, meta in FILES.items():
        dest = trajdir / meta["name"]
        if dest.exists() and md5sum(dest) == meta["md5"]:
            print(f"{key}: present, md5 ok")
            continue
        print(f"{key}: fetching {meta['url']}")
        tmp = dest.with_suffix(".part")
        with urllib.request.urlopen(meta["url"]) as r, open(tmp, "wb") as f:
            while True:
                chunk = r.read(1 << 20)
                if not chunk:
                    break
                f.write(chunk)
        got = md5sum(tmp)
        if got != meta["md5"]:
            tmp.unlink()
            sys.exit(f"{key}: md5 mismatch {got} != {meta['md5']}")
        tmp.rename(dest)
        print(f"{key}: fetched, md5 ok")


def count_frames(path):
    n = 0
    with open(path, "rb") as f:
        for line in f:
            if line.startswith(b"ITEM: TIMESTEP"):
                n += 1
    return n


def summarize_run(outrun):
    """Inventory the files a demo produced; keep the last line of small
    text outputs so the manifest carries the headline numbers."""
    out = {}
    if not outrun.is_dir():
        return out
    for p in sorted(outrun.rglob("*")):
        if not p.is_file():
            continue
        entry = {"bytes": p.stat().st_size}
        if p.suffix in {".txt", ".dat"} and entry["bytes"] < 1 << 16:
            lines = [
                ln for ln in p.read_text(errors="replace").splitlines() if ln.strip()
            ]
            if lines:
                entry["last_line"] = lines[-1]
        out[str(p.relative_to(outrun))] = entry
    return out


def demos(yoda_src, yoda_build, trajdir, outdir):
    outdir = pathlib.Path(outdir)
    yoda_src = pathlib.Path(yoda_src).resolve()
    yoda_build = pathlib.Path(yoda_build).resolve()
    so = yoda_build / "dseams_core.so"
    if not so.is_file():
        sys.exit(f"missing {so}; build yodaStruct first")
    if not SCRIPT.is_file():
        sys.exit(f"missing {SCRIPT}")
    lua = find_lua()
    env = os.environ.copy()
    env["LUA_PATH"] = str(yoda_src / "lua" / "?.lua") + ";;"
    env["LUA_CPATH"] = str(yoda_build / "?.so") + ";;"
    results = {}
    for demo in DEMOS:
        key = demo["key"]
        rundir = outdir / key
        rundir.mkdir(parents=True, exist_ok=True)
        traj = (pathlib.Path(trajdir) / FILES[demo["traj"]]["name"]).resolve()
        frames = count_frames(traj)
        frame = frames if demo["frame"] == "last" else int(demo["frame"])
        env["FIGSHARE_TRAJ"] = str(traj)
        env["FIGSHARE_OUTDIR"] = str(rundir)
        env["FIGSHARE_FRAME"] = str(frame)
        t0 = time.perf_counter()
        proc = subprocess.run(
            [lua, str(SCRIPT), key],
            capture_output=True,
            text=True,
            env=env,
            cwd=str(yoda_src),
        )
        secs = time.perf_counter() - t0
        (rundir / "stdout.log").write_text(proc.stdout)
        (rundir / "stderr.log").write_text(proc.stderr)
        results[key] = {
            "trajectory": FILES[demo["traj"]]["name"],
            "doi": FILES[demo["traj"]]["doi"],
            "frames_in_deposit": frames,
            "frame": frame,
            "returncode": proc.returncode,
            "seconds": round(secs, 2),
            "stdout": proc.stdout.strip(),
            "outputs": summarize_run(rundir),
        }
        status = "ok" if proc.returncode == 0 else f"rc={proc.returncode}"
        print(f"{key}: {status} in {secs:.1f}s", flush=True)
        if proc.returncode != 0:
            print(proc.stderr, file=sys.stderr)
    (outdir / "figshare-demos.json").write_text(
        json.dumps(results, indent=2, sort_keys=True) + "\n"
    )
    failed = [k for k, v in results.items() if v["returncode"] != 0]
    if failed:
        sys.exit(f"failed demos: {failed}")


def main():
    cmd = sys.argv[1] if len(sys.argv) > 1 else ""
    if cmd == "fetch" and len(sys.argv) == 3:
        fetch(sys.argv[2])
    elif cmd == "demos" and len(sys.argv) == 6:
        demos(sys.argv[2], sys.argv[3], sys.argv[4], sys.argv[5])
    else:
        sys.exit(__doc__)


if __name__ == "__main__":
    main()
