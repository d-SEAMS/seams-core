#!/usr/bin/env python3
"""Run the v2 pipeline on the v1 figshare demonstration trajectories.

The d-SEAMS 1.0 paper (10.1021/acs.jcim.0c00031) demonstrated on five
LAMMPS trajectories archived in the figshare project "d-SEAMS Datasets"
(https://figshare.com/projects/d-SEAMS_Datasets/73545). This script
fetches those exact deposits (MD5-verified against figshare's records)
and runs each corresponding published example workflow, unmodified,
with the current yoda binary.

Subcommands (all paths relative to the repository root):

  fetch <trajdir>                       download + verify the deposits
  demos <yoda> <trajdir> <outdir>       run the five example workflows (Lua CLI)

The Python-bindings demonstrations live as percent-format notebooks
under repro/notebooks/; jupytext converts them and papermill executes
them via the figshare_notebook Snakemake rule.
"""
import hashlib
import json
import pathlib
import re
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

# Each demo pairs a figshare trajectory with the example workflow that
# consumed it in the 1.0 paper. The example_lua trees live in
# https://github.com/d-SEAMS/yodaStruct . frame="last" retargets the
# workflow at the final frame of the deposit (the crystallized end of
# the nucleation run); frame=None keeps the example's own frame settings.
DEMOS = [
    {
        "key": "chillPlus",
        "config": "example_lua/chillPlus/config.yml",
        "vars": "example_lua/chillPlus/iceType/vars.lua",
        "traj": "chillplus",
        "frame": None,
    },
    {
        "key": "bulkTopologicalCriterion",
        "config": "example_lua/bulkTopologicalCriterion/config.yml",
        "vars": "example_lua/bulkTopologicalCriterion/iceType/vars.lua",
        "traj": "nucleation",
        "frame": "last",
    },
    {
        "key": "iceNanotube",
        "config": "example_lua/iceNanotube/strictCriterion/config.yml",
        "vars": "example_lua/iceNanotube/strictCriterion/iceType/vars.lua",
        "traj": "nanotube",
        "frame": None,
    },
    {
        "key": "monolayer",
        "config": "example_lua/monolayer/config.yml",
        "vars": "example_lua/monolayer/iceType/vars.lua",
        "traj": "monolayer",
        "frame": None,
    },
    {
        "key": "rdf2D",
        "config": "example_lua/rdf2D-example/config.yml",
        "vars": "example_lua/rdf2D-example/iceType/vars.lua",
        "traj": "rdf2d",
        "frame": None,
    },
]

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


def patch_demo(demo, trajdir, rundir):
    """Write the patched vars.lua and config.yml for one demo; return the
    config path and the frame count of its trajectory."""
    traj = (pathlib.Path(trajdir) / FILES[demo["traj"]]["name"]).resolve()
    frames = count_frames(traj)
    outrun = rundir / "run"
    outrun.mkdir(parents=True, exist_ok=True)

    vtext = pathlib.Path(demo["vars"]).read_text()
    vtext = re.sub(
        r"^outDir\s*=.*$",
        f'outDir="{outrun}/";',
        vtext,
        flags=re.M,
    )
    if demo["frame"] == "last":
        vtext = re.sub(
            r"^targetFrame\s*=.*$", f"targetFrame={frames};", vtext, flags=re.M
        )
        vtext = re.sub(
            r"^finalFrame\s*=.*$", f"finalFrame={frames};", vtext, flags=re.M
        )
    vpath = rundir / "vars.lua"
    vpath.write_text(vtext)

    ctext = pathlib.Path(demo["config"]).read_text()
    ctext = re.sub(
        r'^trajectory:.*$', f'trajectory: "{traj}"', ctext, flags=re.M
    )
    ctext = re.sub(
        r'^variables:.*$', f'variables: "{vpath}"', ctext, flags=re.M
    )
    cpath = rundir / "config.yml"
    cpath.write_text(ctext)
    return cpath, frames


def summarize_run(outrun):
    """Inventory the files a demo produced; keep the last line of small
    text outputs so the manifest carries the headline numbers."""
    out = {}
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


def demos(yoda, trajdir, outdir):
    outdir = pathlib.Path(outdir)
    results = {}
    for demo in DEMOS:
        key = demo["key"]
        rundir = outdir / key
        cpath, frames = patch_demo(demo, trajdir, rundir)
        t0 = time.perf_counter()
        proc = subprocess.run(
            [yoda, "-c", str(cpath)],
            capture_output=True,
            text=True,
        )
        secs = time.perf_counter() - t0
        (rundir / "stdout.log").write_text(proc.stdout)
        (rundir / "stderr.log").write_text(proc.stderr)
        results[key] = {
            "trajectory": FILES[demo["traj"]]["name"],
            "doi": FILES[demo["traj"]]["doi"],
            "frames_in_deposit": frames,
            "returncode": proc.returncode,
            "seconds": round(secs, 2),
            "outputs": summarize_run(rundir / "run"),
        }
        status = "ok" if proc.returncode == 0 else f"rc={proc.returncode}"
        print(f"{key}: {status} in {secs:.1f}s", flush=True)
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
    elif cmd == "demos" and len(sys.argv) == 5:
        demos(sys.argv[2], sys.argv[3], sys.argv[4])
    else:
        sys.exit(__doc__)


if __name__ == "__main__":
    main()
