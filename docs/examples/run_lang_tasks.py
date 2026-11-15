#!/usr/bin/env python3
# Copyright (c) d-SEAMS developers
"""Run the C or C++ example of the ice Ih census.

Exit status 77 means the compiler is not on this machine. The Python,
Rust, Julia, Lua, Fortran and R copies of this census run in the binding
repositories.
A successful run prints one line. ``parity`` runs every language, checks
that the successful lines match, and writes the table Sphinx includes.
"""

import os
import subprocess
import sys
import tempfile
from pathlib import Path

ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT / "docs"))
from toolchain import c_compiler, cxx_compiler, link_command  # noqa: E402
EXAMPLES = ROOT / "docs" / "examples"
CANON = "pairs 512 rings 0 0 0 0 0 0 256 chill 0 128 0 0 0 0 0 0 0"

LANGS = ("c", "cpp")


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
    cc = c_compiler()
    work = Path(tempfile.mkdtemp(prefix="seams-lang-c-"))
    exe = work / "lang_tasks_c"
    lib = build / "src"
    text = compile_and_run(
        link_command(
            cc,
            [
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
        ),
        exe,
        env,
    )
    sys.stdout.write(text)


def lang_cpp(build, env):
    cxx = cxx_compiler()
    work = Path(tempfile.mkdtemp(prefix="seams-lang-cpp-"))
    exe = work / "lang_tasks_cpp"
    src = build / "src"
    bindings = build / "bindings"
    text = compile_and_run(
        link_command(
            cxx,
            [
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
        ),
        exe,
        env,
    )
    sys.stdout.write(text)


DRIVERS = {
    "c": lang_c,
    "cpp": lang_cpp,
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
        note = "No language was skipped in this environment."
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
