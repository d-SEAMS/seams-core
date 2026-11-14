#!/usr/bin/env python3
"""Run the minimal examples and require each public symbol to name one.

Usage: run_examples.py ROOT BUILD
ROOT is the source tree. BUILD is the meson build tree.
"""

import importlib.util
import os
import re
import shutil
import subprocess
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT / "docs"))

import check_c_symbols as symbols  # noqa: E402
import check_examples  # noqa: E402

FIXTURE = "input/traj/genice_sI.lammpstrj"
EXAMPLES = {
    # The C page cites the tour, which calls every exported function.
    # minimal_c.c is the short quickstart and does not name the macros.
    "c": "docs/examples/c_api_tour.c",
    "cxx": "docs/examples/minimal_cpp.cpp",
    "python": "docs/examples/minimal_py.py",
    "rust": "docs/examples/minimal_rust.rs",
    "fortran": "docs/examples/minimal_fortran.f90",
    "julia": "docs/examples/minimal_julia.jl",
    "r": "docs/examples/minimal_r.R",
}


def library_path(build):
    parts = [str(build / "src"), str(build / "bindings")]
    current = os.environ.get("LD_LIBRARY_PATH")
    if current:
        parts.append(current)
    return os.pathsep.join(parts)


def run(cmd, env=None, cwd=None):
    subprocess.check_call(cmd, env=env, cwd=cwd)


def page_blobs():
    errors = []
    corpus = check_examples.harvest(errors)
    if errors:
        raise SystemExit("\n".join(errors))
    pages = []
    for folder in ("api", "engines", "tutorials"):
        pages.extend(sorted((ROOT / "docs" / folder).glob("*.md")))
    pages.extend(sorted((ROOT / "docs" / "orgmode").rglob("*.org")))
    blobs = []
    for doc in pages:
        blobs.append(doc.read_text() + "\n" + corpus.text(doc))
    return blobs


def require_example_links():
    groups = {
        "c": symbols.c_header_symbols(),
        "cxx": symbols.cxx_symbols(),
        "python": symbols.python_symbols(),
        "rust": symbols.rust_symbols(),
        "fortran": symbols.fortran_symbols(),
        "julia": symbols.julia_symbols(),
        "r": symbols.r_symbols(),
    }
    blobs = page_blobs()
    missing = []
    for lang, names in groups.items():
        needle = EXAMPLES[lang]
        for name in names:
            pattern = re.compile(r"(?<!\w)" + re.escape(name) + r"(?!\w)")
            if not any(needle in blob and pattern.search(blob) for blob in blobs):
                missing.append(lang + ":" + name)
    if missing:
        raise SystemExit("symbols without an example link:\n" + "\n".join(missing))
    print("example links", sum(len(v) for v in groups.values()))


def compile_and_run(build):
    work = build / "docs-examples"
    work.mkdir(parents=True, exist_ok=True)
    fixture = str(ROOT / FIXTURE)
    env = os.environ.copy()
    env["LD_LIBRARY_PATH"] = library_path(build)
    cc = shutil.which("cc") or shutil.which("gcc")
    cxx = shutil.which("c++") or shutil.which("g++")
    if cc is None or cxx is None:
        raise SystemExit("cc and c++ are required")
    c_bin = work / "minimal_c"
    run(
        [
            cc,
            "-I",
            str(ROOT / "include"),
            "-I",
            str(build / "src"),
            str(ROOT / "docs/examples/minimal_c.c"),
            "-L" + str(build / "src"),
            "-lseams",
            "-Wl,-rpath," + str(build / "src"),
            "-o",
            str(c_bin),
        ]
    )
    run([str(c_bin), fixture], env=env)
    cpp_bin = work / "minimal_cpp"
    run(
        [
            cxx,
            "-I",
            str(ROOT / "bindings/cxx"),
            str(ROOT / "docs/examples/minimal_cpp.cpp"),
            "-L" + str(build / "src"),
            "-L" + str(build / "bindings"),
            "-lseams_lang",
            "-lseams",
            "-Wl,-rpath," + str(build / "src"),
            "-Wl,-rpath," + str(build / "bindings"),
            "-o",
            str(cpp_bin),
        ]
    )
    run([str(cpp_bin), fixture], env=env)
    run(
        [sys.executable, str(ROOT / "docs/examples/minimal_py.py"), fixture, str(build / "src")],
        env=env,
    )
    seams = build / "seams"
    if not seams.is_file():
        raise SystemExit("seams CLI was not built")
    run(["sh", str(ROOT / "docs/examples/minimal_cli.sh"), str(seams), fixture], env=env)
    gfortran = shutil.which("gfortran")
    if gfortran is None:
        raise SystemExit("gfortran is not on PATH")
    run(
        [
            gfortran,
            "-J",
            str(work),
            "-c",
            str(ROOT / "bindings/fortran/seams.f90"),
            "-o",
            str(work / "seams.o"),
        ]
    )
    f_bin = work / "minimal_fortran"
    run(
        [
            gfortran,
            "-I",
            str(work),
            str(ROOT / "docs/examples/minimal_fortran.f90"),
            str(work / "seams.o"),
            "-o",
            str(f_bin),
            "-L" + str(build / "src"),
            "-lseams",
            "-Wl,-rpath," + str(build / "src"),
        ]
    )
    run([str(f_bin), fixture], env=env)
    cargo = shutil.which("cargo")
    if cargo is None:
        raise SystemExit("cargo is not on PATH")
    rust_env = env.copy()
    rust_env["SEAMS_LIB_DIR"] = str(build / "src")
    run(
        [
            cargo,
            "run",
            "--offline",
            "--quiet",
            "--manifest-path",
            str(ROOT / "bindings/rust/Cargo.toml"),
            "--example",
            "minimal",
            "--",
            fixture,
        ],
        env=rust_env,
        cwd=ROOT,
    )


def binding_helpers():
    path = ROOT / "docs" / "examples" / "run_binding.py"
    spec = importlib.util.spec_from_file_location("run_binding", path)
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


def run_julia(build, fixture):
    julia = binding_helpers().julia_bin()
    env = os.environ.copy()
    env["DSEAMS_ROOT"] = str(ROOT)
    env["SEAMS_LANG_LIB"] = str(build / "bindings" / "libseams_lang.so")
    env["LD_LIBRARY_PATH"] = library_path(build)
    run(
        [
            julia,
            "--project=" + str(ROOT / "bindings/julia/Dseams"),
            str(ROOT / "docs/examples/minimal_julia.jl"),
            fixture,
        ],
        env=env,
        cwd=ROOT,
    )


def run_r(build, fixture):
    script, env = binding_helpers().install_r_package(build)
    run([script, str(ROOT / "docs/examples/minimal_r.R"), fixture], env=env)


def run_lua(build, fixture):
    lua = shutil.which("lua5.4") or shutil.which("lua")
    if lua is None:
        raise SystemExit("lua is not on PATH")
    modules = list((build / "bindings").rglob("dseams_core.so"))
    if not modules:
        raise SystemExit("dseams_core.so was not built")
    env = os.environ.copy()
    env["LUA_PATH"] = str(ROOT / "bindings" / "lua") + "/?.lua;;"
    env["LUA_CPATH"] = str(modules[0].parent) + "/?.so;;"
    env["LD_LIBRARY_PATH"] = library_path(build)
    env["OMP_NUM_THREADS"] = "1"
    done = subprocess.run(
        [lua, str(ROOT / "docs/examples/minimal_lua.lua"), fixture],
        env=env,
        cwd=ROOT,
        text=True,
        stdout=subprocess.PIPE,
        stderr=subprocess.STDOUT,
    )
    if done.returncode != 0 or "lua clathrate 46" not in done.stdout:
        raise SystemExit("lua example failed\n{}".format(done.stdout[-800:]))


def main():
    if len(sys.argv) != 3:
        raise SystemExit("usage: run_examples.py ROOT BUILD")
    build = Path(sys.argv[2])
    require_example_links()
    fixture = str(ROOT / FIXTURE)
    compile_and_run(build)
    run_julia(build, fixture)
    run_r(build, fixture)
    run_lua(build, fixture)
    print("docs examples ok")


if __name__ == "__main__":
    main()
