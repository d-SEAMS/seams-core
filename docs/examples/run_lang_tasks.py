#!/usr/bin/env python3
"""Run one language example of the ice Ih census.

Exit status 77 means the compiler or the binding is not on this machine.
A successful run prints one line. ``parity`` runs every language, checks
that the successful lines match, and writes the table Sphinx includes.
"""

import os
import shutil
import subprocess
import sys
import tempfile
from pathlib import Path

ROOT = Path(__file__).resolve().parents[2]
EXAMPLES = ROOT / "docs" / "examples"
CANON = "pairs 512 rings 0 0 0 0 0 0 256 chill 0 128 0 0 0 0 0 0 0"

LANGS = ("c", "cpp", "python", "rust", "julia", "lua", "fortran", "r")


class Skip(Exception):
    pass


def library_path(build):
    parts = [
        str(build / "src"),
        str(build / "bindings"),
        str(build / "subprojects" / "readcon-core" / "cargo-target" / "release"),
        str(build / "subprojects" / "readcon-db" / "cargo-target" / "release"),
    ]
    current = os.environ.get("LD_LIBRARY_PATH")
    if current:
        parts.append(current)
    return os.pathsep.join(parts)


def base_env(build):
    env = os.environ.copy()
    env["DSEAMS_ROOT"] = str(ROOT)
    env["LD_LIBRARY_PATH"] = library_path(build)
    env["OMP_NUM_THREADS"] = "1"
    lang = build / "bindings" / "libseams_lang.so"
    if lang.is_file():
        env["SEAMS_LANG_LIB"] = str(lang)
    return env


def run(cmd, env, cwd=ROOT):
    completed = subprocess.run(cmd, env=env, cwd=cwd, text=True, capture_output=True)
    if completed.returncode != 0:
        sys.stderr.write(completed.stderr)
        sys.stderr.write(completed.stdout)
        raise SystemExit(completed.returncode)
    return completed.stdout


def compile_and_run(cmd, exe, env):
    run(cmd, env)
    return run([str(exe)], env)


def lang_c(build, env):
    cc = shutil.which("cc") or shutil.which("gcc")
    if not cc:
        raise Skip("cc is not on PATH")
    work = Path(tempfile.mkdtemp(prefix="seams-lang-c-"))
    exe = work / "lang_tasks_c"
    lib = build / "src"
    text = compile_and_run(
        [
            cc,
            "-I",
            str(ROOT / "include"),
            "-I",
            str(build / "src"),
            str(EXAMPLES / "lang_tasks_c.c"),
            "-o",
            str(exe),
            "-L" + str(lib),
            "-lseams",
            "-Wl,-rpath," + str(lib),
        ],
        exe,
        env,
    )
    sys.stdout.write(text)


def lang_cpp(build, env):
    cxx = shutil.which("c++") or shutil.which("g++")
    if not cxx:
        raise Skip("c++ is not on PATH")
    work = Path(tempfile.mkdtemp(prefix="seams-lang-cpp-"))
    exe = work / "lang_tasks_cpp"
    src = build / "src"
    bindings = build / "bindings"
    text = compile_and_run(
        [
            cxx,
            "-std=c++17",
            "-I",
            str(ROOT / "bindings" / "cxx"),
            str(EXAMPLES / "lang_tasks_cpp.cpp"),
            "-o",
            str(exe),
            "-L" + str(bindings),
            "-L" + str(src),
            "-lseams_lang",
            "-lseams",
            "-Wl,-rpath," + str(bindings),
            "-Wl,-rpath," + str(src),
        ],
        exe,
        env,
    )
    sys.stdout.write(text)


def lang_python(build, env):
    try:
        import numpy  # noqa: F401
    except ImportError as err:
        raise Skip("numpy is not importable") from err
    if "SEAMS_LANG_LIB" not in env:
        raise Skip("libseams_lang.so was not built")
    text = run([sys.executable, str(EXAMPLES / "lang_tasks_python.py")], env)
    sys.stdout.write(text)


def lang_rust(build, env):
    rustc = shutil.which("rustc")
    if not rustc:
        raise Skip("rustc is not on PATH")
    if not (build / "src" / "libseams.so").is_file() and not (build / "src" / "libseams.a").is_file():
        raise Skip("libseams was not built")
    work = Path(tempfile.mkdtemp(prefix="seams-lang-rust-"))
    exe = work / "lang_tasks_rs"
    src = build / "src"
    bindings = build / "bindings"
    text = compile_and_run(
        [
            rustc,
            "--edition",
            "2021",
            "-o",
            str(exe),
            str(EXAMPLES / "lang_tasks.rs"),
            "-L",
            str(bindings),
            "-L",
            str(src),
            "-C",
            "link-arg=-Wl,-rpath," + str(bindings),
            "-C",
            "link-arg=-Wl,-rpath," + str(src),
        ],
        exe,
        env,
    )
    sys.stdout.write(text)


def lang_julia(build, env):
    julia = shutil.which("julia")
    if not julia:
        raise Skip("julia is not on PATH")
    if "SEAMS_LANG_LIB" not in env:
        raise Skip("libseams_lang.so was not built")
    env = dict(env)
    env["JULIA_LOAD_PATH"] = str(ROOT / "bindings" / "julia" / "Dseams") + os.pathsep
    text = run([julia, "--startup-file=no", str(EXAMPLES / "lang_tasks.jl")], env)
    sys.stdout.write(text)


def lua_module(build):
    found = list((build / "bindings" / "lua").glob("*/dseams_core.so"))
    return found[0] if found else None


def lang_lua(build, env):
    lua = shutil.which("lua5.4") or shutil.which("lua")
    if not lua:
        raise Skip("lua is not on PATH")
    module = lua_module(build)
    if module is None:
        raise Skip("the Lua module was not built")
    env = dict(env)
    env["LUA_CPATH"] = str(module.parent / "?.so") + ";;"
    env["LUA_PATH"] = str(ROOT / "bindings" / "lua" / "?.lua") + ";;"
    text = run([lua, str(EXAMPLES / "lang_tasks.lua")], env)
    sys.stdout.write(text)


def lang_fortran(build, env):
    gfortran = shutil.which("gfortran")
    if not gfortran:
        raise Skip("gfortran is not on PATH")
    if "SEAMS_LANG_LIB" not in env:
        raise Skip("libseams_lang.so was not built")
    work = Path(tempfile.mkdtemp(prefix="seams-lang-f90-"))
    exe = work / "lang_tasks_f90"
    src = build / "src"
    bindings = build / "bindings"
    text = compile_and_run(
        [
            gfortran,
            str(EXAMPLES / "lang_tasks.f90"),
            "-o",
            str(exe),
            "-L" + str(bindings),
            "-L" + str(src),
            "-lseams_lang",
            "-lseams",
            "-Wl,-rpath," + str(bindings),
            "-Wl,-rpath," + str(src),
        ],
        exe,
        env,
    )
    sys.stdout.write(text)


def rscript():
    found = shutil.which("Rscript")
    if found:
        return found
    raise Skip("Rscript is not on PATH")


def lang_r(build, env):
    if shutil.which("Rscript") is None:
        raise Skip("Rscript is not on PATH")
    sys.path.insert(0, str(ROOT / "docs" / "examples"))
    import run_binding

    script, renv = run_binding.install_r_package(build)
    merged = dict(env)
    merged.update(renv)
    text = run([script, str(EXAMPLES / "lang_tasks.R")], merged)
    sys.stdout.write(text)


DRIVERS = {
    "c": lang_c,
    "cpp": lang_cpp,
    "python": lang_python,
    "rust": lang_rust,
    "julia": lang_julia,
    "lua": lang_lua,
    "fortran": lang_fortran,
    "r": lang_r,
}


def one(name, build):
    env = base_env(build)
    before = sys.stdout.tell() if sys.stdout.seekable() else None
    try:
        DRIVERS[name](build, env)
    except Skip as err:
        print(err, file=sys.stderr)
        raise SystemExit(77) from err
    if before is None:
        return
    # Compiled drivers print through the child, which this process already
    # forwarded. Nothing else to do.


def capture(name, build):
    env = base_env(build)
    try:
        completed = subprocess.run(
            [sys.executable, str(Path(__file__).resolve()), name, str(build)],
            text=True,
            capture_output=True,
            env=env,
        )
    except OSError as err:
        return 77, "", str(err)
    return completed.returncode, completed.stdout.strip(), completed.stderr.strip()


def write_table(rows, destination):
    skipped = [name for name, code, line, _err in rows if code == 77]
    if skipped:
        note = "Skipped in this environment: " + ", ".join(f"``{name}``" for name in skipped) + "."
    else:
        note = "No binding was skipped in this environment."
    body = [
        note,
        "",
        "The table is written by ``docs/examples/run_lang_tasks.py parity``.",
        "A row that ran prints the same integers as every other row that ran.",
        "",
        ".. list-table::",
        "   :header-rows: 1",
        "   :widths: 12 10 78",
        "",
        "   * - Language",
        "     - Status",
        "     - Census",
    ]
    for name, code, line, _err in rows:
        if code == 0:
            status = "ran"
            census = line
        elif code == 77:
            status = "skipped"
            census = ""
        else:
            status = f"failed ({code})"
            census = ""
        body.append(f"   * - {name}")
        body.append(f"     - {status}")
        body.append(f"     - ``{census}``" if census else "     -")
    destination.write_text("\n".join(body) + "\n", encoding="utf-8")


def parity(build, write):
    rows = []
    lines = []
    for name in LANGS:
        code, line, err = capture(name, build)
        rows.append((name, code, line, err))
        if code == 0:
            lines.append((name, line))
            if line != CANON:
                print(f"{name} printed {line!r}, expected {CANON!r}", file=sys.stderr)
        elif code != 77:
            print(f"{name} exited {code}: {err}", file=sys.stderr)
    if write:
        write_table(rows, ROOT / "docs" / "source" / "lang" / "parity.rst")
    if not lines:
        print("no language ran", file=sys.stderr)
        return 1
    first = lines[0][1]
    for name, line in lines[1:]:
        if line != first:
            print(f"{name} disagrees with {lines[0][0]}", file=sys.stderr)
            return 1
    if any(code not in (0, 77) for _name, code, _line, _err in rows):
        return 1
    if any(line != CANON for _name, line in lines):
        return 1
    skipped = [name for name, code, _line, _err in rows if code == 77]
    print(CANON)
    if skipped:
        print("skipped: " + ", ".join(skipped), file=sys.stderr)
    return 0


def main():
    if len(sys.argv) < 3:
        print(f"usage: {sys.argv[0]} <{'|'.join(LANGS)}|parity> <build>", file=sys.stderr)
        return 2
    name = sys.argv[1]
    build = Path(sys.argv[2]).resolve()
    if name == "parity":
        write = "--write" in sys.argv[3:]
        return parity(build, write)
    if name not in DRIVERS:
        print(f"unknown language {name}", file=sys.stderr)
        return 2
    one(name, build)
    return 0


if __name__ == "__main__":
    sys.exit(main())
