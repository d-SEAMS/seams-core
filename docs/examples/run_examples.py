#!/usr/bin/env python3
# Copyright (c) d-SEAMS developers
"""Run the C, C++ and CLI minimal examples and require each public symbol to name one.

Usage: run_examples.py ROOT BUILD
ROOT is the source tree. BUILD is the meson build tree. The Python, Rust,
Fortran, Julia, R and Lua minimal examples live in the binding
repositories.
"""

import os
import re
import subprocess
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT / "docs"))
from toolchain import c_compiler, cxx_compiler, link_command  # noqa: E402

import check_c_symbols as symbols  # noqa: E402
import check_examples  # noqa: E402

FIXTURE = "input/traj/genice_sI.lammpstrj"
EXAMPLES = {
    # The C page cites the tour, which calls every exported function.
    # minimal_c.c is the short quickstart and does not name the macros.
    "c": "docs/examples/c_api_tour.c",
    "cxx": "docs/examples/minimal_cpp.cpp",
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
    cc = c_compiler()
    cxx = cxx_compiler()
    c_bin = work / "minimal_c"
    run(
        link_command(
            cc,
            [
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
            ],
        )
    )
    run([str(c_bin), fixture], env=env)
    cpp_bin = work / "minimal_cpp"
    run(
        link_command(
            cxx,
            [
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
            ],
        )
    )
    run([str(cpp_bin), fixture], env=env)
    seams = build / "seams"
    if not seams.is_file():
        raise SystemExit("seams CLI was not built")
    run(["sh", str(ROOT / "docs/examples/minimal_cli.sh"), str(seams), fixture], env=env)


def main():
    if len(sys.argv) != 3:
        raise SystemExit("usage: run_examples.py ROOT BUILD")
    build = Path(sys.argv[2])
    require_example_links()
    compile_and_run(build)
    print("docs examples ok")


if __name__ == "__main__":
    main()
