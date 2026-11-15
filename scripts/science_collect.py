#!/usr/bin/env python3
# Copyright (c) d-SEAMS developers

"""Analyse the finished 230/240/250 K basal dumps and render the end frames.

Expects data/science-frames/tip4p-basal-<T>K.lammpstrj and the matching
meta file. Runs the threaded science-front kernels, writes the speed
table, checks the last frame of each trajectory against OVITO and freud,
and renders the closing slab.
"""

from __future__ import annotations

import os
import subprocess
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
DATA = ROOT / "paper" / "tex" / "figures" / "data" / "science"
FRAMES = ROOT / "data" / "science-frames"
TEMPS = (230, 240, 250)
COMPAT = "/cvmfs/software.eessi.io/versions/2026.06/compat/linux/x86_64/lib64"


def env_for_front() -> dict[str, str]:
    env = os.environ.copy()
    libs = [COMPAT, "/lib/x86_64-linux-gnu"]
    current = [part for part in env.get("LD_LIBRARY_PATH", "").split(":") if part]
    env["LD_LIBRARY_PATH"] = ":".join(libs + [part for part in current if part not in libs])
    env["OMP_NUM_THREADS"] = env.get("SCIENCE_THREADS", "2")
    return env


def meta_value(temp: int, key: str) -> str:
    for line in (DATA / f"interface-meta-{temp}K.csv").read_text().splitlines():
        parts = line.split()
        if parts and parts[0] == key:
            return parts[1]
    raise SystemExit(f"missing {key} in the {temp} K meta file")


def run(cmd: list[str], env: dict[str, str] | None = None) -> None:
    print(" ".join(cmd), flush=True)
    subprocess.run(cmd, check=True, env=env, cwd=ROOT)


def main() -> int:
    if len(sys.argv) != 2:
        raise SystemExit("usage: science_collect.py FIGDIR")
    figdir = Path(sys.argv[1])
    figdir.mkdir(parents=True, exist_ok=True)
    front = os.environ.get("SCIENCE_FRONT", str(ROOT / "bbdir" / "science-front"))
    py = os.environ.get("SCIENCE_PYTHON", sys.executable)
    threads = os.environ.get("SCIENCE_THREADS", "2")
    child = env_for_front()
    for temp in TEMPS:
        dump = FRAMES / f"tip4p-basal-{temp}K.lammpstrj"
        if not dump.is_file():
            raise SystemExit(f"missing {dump}")
        out = DATA / f"interface-{temp}K"
        out.mkdir(parents=True, exist_ok=True)
        run(
            [
                front,
                "interface",
                str(dump),
                "--ps",
                meta_value(temp, "dump_ps"),
                "--bin",
                "1",
                "--out-dir",
                str(out),
                "--threads",
                threads,
            ],
            child,
        )
    run([py, str(ROOT / "scripts" / "science_growth.py")], child)
    for name in (
        "interface-layers",
        "interface-fronts",
        "interface-velocity",
        "interface-atoms",
        "interface-causes",
    ):
        source = DATA / "interface-230K" / f"{name}.csv"
        target = DATA / f"{name}.csv"
        target.write_text(source.read_text())
    # Last frame of each temperature, for OVITO and freud.
    sys.path.insert(0, str(ROOT / "scripts"))
    import science_render as render

    for temp in TEMPS:
        dump = FRAMES / f"tip4p-basal-{temp}K.lammpstrj"
        nframes = render.count_frames(dump)
        box, pos = render.read_dump_frame(dump, nframes)
        work = FRAMES / "render-work" / f"end-{temp}"
        work.mkdir(parents=True, exist_ok=True)
        one = work / "one.lammpstrj"
        render.write_one_frame(one, box, pos)
        render.run_front(one, work, "interface")
        (DATA / f"interface-atoms-{temp}K-end.csv").write_text((work / "interface-atoms.csv").read_text())
        (DATA / f"interface-causes-{temp}K-end.csv").write_text((work / "interface-causes.csv").read_text())
        print(f"labeled {temp} K end frame {nframes}", flush=True)
    run([py, str(ROOT / "scripts" / "science_validate.py")], child)
    render_env = child.copy()
    render_env["QT_QPA_PLATFORM"] = "ovitoheadless"
    render_env["LIBGL_ALWAYS_SOFTWARE"] = "1"
    render_env["OVITO_THREAD_COUNT"] = "1"
    render_env.setdefault("VK_DRIVER_FILES", "/usr/share/vulkan/icd.d/lvp_icd.json")
    render_env["VK_ICD_FILENAMES"] = render_env["VK_DRIVER_FILES"]
    render_env["LP_NUM_THREADS"] = "1"
    run(
        [
            py,
            str(ROOT / "scripts" / "science_render.py"),
            "--out",
            str(figdir),
            "--interface-dump",
            str(FRAMES / "tip4p-basal-230K.lammpstrj"),
            "--interface-frame",
            "-1",
            "--nucleus-dump",
            str(FRAMES / "nucleation.lammpstrj"),
            "--nucleus-min-frame",
            "2",
            "--nucleus-late-frame",
            "-1",
        ],
        render_env,
    )
    run([py, str(ROOT / "scripts" / "science_figures.py"), "--out", str(figdir)], child)
    print("collected", flush=True)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
