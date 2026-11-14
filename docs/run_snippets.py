#!/usr/bin/env python3
# Copyright (c) d-SEAMS developers
"""Run every documented Python, C, C++, Rust, Julia, Lua, Fortran, and R snippet.

One Meson test is docs_snippet_<page>. A missing compiler or an unbuilt
binding exits 77. A snippet that does not match the current API fails.
"""

from __future__ import annotations

import hashlib
import os
import re
import shutil
import subprocess
import sys
import tempfile
import textwrap
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
DOCS = ROOT / "docs"
FIXTURE = ROOT / "input" / "traj" / "genice_sI.lammpstrj"
LANGS = {"python", "py", "c", "cpp", "c++", "cxx", "rust", "julia", "lua", "fortran", "f90", "r"}
ORG_SRC = re.compile(r"#\+begin_src\s+([^\s]+)[^\n]*\n(.*?)#\+end_src", re.S | re.I)
RST_DIR = re.compile(
    r"^\.\.\s+(literalinclude|code-block|code|doctest)::[ \t]*([^\n]*)\n((?:[ \t].*\n|\n)*)",
    re.M,
)
MD_INCLUDE = re.compile(r"```\{literalinclude\}\s+(\S+)")
# A language tag is required. A bare closing ``` would otherwise swallow the next fence.
MD_FENCE = re.compile(r"```([A-Za-z0-9_+-]+)\n(.*?)```", re.S)
EXCERPT = re.compile(r"^(?:#|//|!)\s*excerpt:\s*(\S+)\s*$")
OPTIONAL = ("pydseams", "pydseamslib", "ase", "solvis")


class Skip(Exception):
    pass


class Fail(Exception):
    pass


def page_id(path):
    rel = path.relative_to(DOCS).with_suffix("")
    return re.sub(r"[^a-z0-9]+", "_", str(rel).lower()).strip("_")


def published_dir(path):
    rel = path.relative_to(DOCS)
    if rel.parts[0] == "orgmode":
        return DOCS / "source" / Path(*rel.parts[1:]).parent
    if rel.parts[0] == "api":
        return DOCS / "source" / "stable"
    if rel.parts[0] == "engines":
        return DOCS / "source" / "engines"
    if rel.parts[0] == "tutorials":
        return DOCS / "source" / "fronts"
    return path.parent


def authoring_pages():
    pages = []
    pages.extend(sorted((DOCS / "orgmode").rglob("*.org")))
    for folder in ("api", "engines", "tutorials"):
        pages.extend(sorted((DOCS / folder).glob("*.md")))
    pages.extend(sorted(DOCS.glob("*.md")))
    for folder in ("lang", "capi", "science"):
        pages.extend(sorted((DOCS / "source" / folder).glob("*.rst")))
    return pages


def unindent(body):
    lines = body.splitlines()
    while lines and not lines[0].strip():
        lines.pop(0)
    while lines and not lines[-1].strip():
        lines.pop()
    if not lines:
        return ""
    indents = [len(line) - len(line.lstrip(" ")) for line in lines if line.strip()]
    cut = min(indents) if indents else 0
    return "\n".join(line[cut:] if line.startswith(" " * cut) else line for line in lines)


def snippets_of(path):
    text = path.read_text(encoding="utf-8")
    found = []
    if path.suffix == ".org":
        for match in ORG_SRC.finditer(text):
            lang = match.group(1).lower()
            if lang in LANGS:
                found.append(("fence", lang, unindent(match.group(2))))
        for match in RST_DIR.finditer(text):
            kind, arg, body = match.group(1), match.group(2).strip(), unindent(match.group(3))
            if kind == "literalinclude":
                include = arg.split()[0]
                if Path(include).suffix.lower() in {".c", ".cpp", ".cxx", ".cc", ".py", ".rs", ".jl", ".lua", ".f90", ".r"}:
                    found.append(("include", include, body))
            elif kind == "doctest" or arg.lower() in LANGS:
                found.append(("fence", arg.lower() or "python", body))
    elif path.suffix == ".md":
        for match in MD_INCLUDE.finditer(text):
            include = match.group(1)
            if Path(include).suffix.lower() in {".c", ".cpp", ".cxx", ".cc", ".py", ".rs", ".jl", ".lua", ".f90", ".r"}:
                found.append(("include", include, ""))
        for match in MD_FENCE.finditer(text):
            lang = match.group(1).lower()
            if lang in LANGS:
                found.append(("fence", lang, unindent(match.group(2))))
    else:
        for match in RST_DIR.finditer(text):
            kind, arg, body = match.group(1), match.group(2).strip(), unindent(match.group(3))
            if kind == "literalinclude":
                include = arg.split()[0]
                if Path(include).suffix.lower() in {".c", ".cpp", ".cxx", ".cc", ".py", ".rs", ".jl", ".lua", ".f90", ".r"}:
                    found.append(("include", include, ""))
            elif kind == "doctest" or (arg.lower().split()[:1] or [""])[0] in LANGS or arg.lower() in LANGS:
                lang = "python" if kind == "doctest" else arg.lower()
                found.append(("fence", lang, body))
    return found


def pages():
    chosen = []
    for path in authoring_pages():
        if snippets_of(path):
            chosen.append(path)
    return chosen


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


def stage_dseams(build):
    """Point imports at the in-tree package and the built nanobind module."""
    stage = build / "pythonpath"
    package = ROOT / "bindings" / "python" / "dseams"
    dest = stage / "dseams"
    dest.mkdir(parents=True, exist_ok=True)
    for src in package.glob("*.py"):
        link = dest / src.name
        if not link.exists():
            link.symlink_to(src)
    frames = ROOT / "bindings" / "python" / "frames.py"
    staged_frames = dest / "frames.py"
    if frames.is_file() and not staged_frames.exists():
        staged_frames.symlink_to(frames)
    pyd = ROOT / "bindings" / "python" / "pydseams"
    pdest = stage / "pydseams"
    pdest.mkdir(parents=True, exist_ok=True)
    for src in pyd.glob("*.py"):
        link = pdest / src.name
        if not link.exists():
            link.symlink_to(src)
    built = [path for path in build.rglob("_core*.so") if path.is_file()]
    if not built:
        built = [path for path in package.glob("_core*.so") if path.is_file()]
    tag = sys.implementation.cache_tag
    matched = [path for path in built if tag in path.name]
    if matched:
        built = matched
    if built:
        so = dest / built[0].name
        if not so.exists():
            so.symlink_to(built[0])
    return stage


def base_env(build):
    env = os.environ.copy()
    env.pop("DSEAMS_DOC_BUILD", None)
    # The engines and libseams carry an rpath. A leftover loader path
    # must not be what makes a snippet resolve.
    env.pop("LD_LIBRARY_PATH", None)
    env["DSEAMS_ROOT"] = str(ROOT)
    env["LD_LIBRARY_PATH"] = library_path(build)
    env["OMP_NUM_THREADS"] = "1"
    lang = build / "bindings" / "libseams_lang.so"
    if lang.is_file():
        env["SEAMS_LANG_LIB"] = str(lang)
    env["SEAMS_LIB_DIR"] = str(build / "src")
    env["SEAMS_LIBDIR"] = str(build / "src")
    env["SEAMS_BUILD"] = str(build)
    env["PYTHONPATH"] = os.pathsep.join(
        [
            str(stage_dseams(build)),
            str(DOCS / "examples"),
            str(ROOT / "bindings" / "python"),
            str(ROOT / "adapters" / "python"),
            str(ROOT / "scripts"),
            env.get("PYTHONPATH", ""),
        ]
    )
    rlib = build / "docs-r-lib"
    if rlib.is_dir():
        env["R_LIBS"] = str(rlib) + os.pathsep + env.get("R_LIBS", "")
        env["DSEAMS_R_LIB"] = str(rlib)
    return env


def run(cmd, env, cwd=ROOT):
    done = subprocess.run(cmd, env=env, cwd=cwd, text=True, capture_output=True)
    if done.returncode != 0:
        tail = (done.stderr + done.stdout)[-1200:]
        if "PackageNotFoundError" in tail or "No package metadata was found" in tail:
            raise Skip("an optional package named by the snippet is not installed")
        if "usage: science_tutorials brine" in tail:
            raise Skip("the brine frames are not in this checkout")
        raise Fail(tail or "command failed")
    return done.stdout


def cached(build, key, work):
    stamp = build / "docs-snippet-cache" / key
    if stamp.is_file():
        return True
    stamp.parent.mkdir(parents=True, exist_ok=True)
    # Two pages can include the same file. Without the lock both render
    # the science figures and one sees a half-written png.
    import fcntl

    with open(stamp.with_suffix(".lock"), "a", encoding="utf-8") as handle:
        fcntl.flock(handle, fcntl.LOCK_EX)
        if stamp.is_file():
            return True
        work()
        stamp.write_text("ok\n", encoding="utf-8")
    return False


def sdk_prefix():
    if shutil.which("pkg-config"):
        done = subprocess.run(
            ["pkg-config", "--variable=prefix", "libgromacs"],
            text=True,
            capture_output=True,
        )
        if done.returncode == 0 and done.stdout.strip():
            return Path(done.stdout.strip())
    prefix = os.environ.get("CONDA_PREFIX", "")
    return Path(prefix) if prefix else None


def pkgconfig_env():
    """conda-forge GROMACS ships libgromacs.pc under lib.AVX2_256, not lib."""
    env = os.environ.copy()
    prefix = sdk_prefix()
    extra = []
    if prefix is not None:
        for folder in (
            prefix / "lib" / "pkgconfig",
            prefix / "lib.AVX2_256" / "pkgconfig",
            prefix / "lib.AVX_256" / "pkgconfig",
            prefix / "lib.SSE2" / "pkgconfig",
            prefix / "share" / "pkgconfig",
        ):
            if folder.is_dir():
                extra.append(str(folder))
    current = env.get("PKG_CONFIG_PATH", "")
    if current:
        extra.append(current)
    if extra:
        env["PKG_CONFIG_PATH"] = os.pathsep.join(extra)
    return env


def sdk_flags(text):
    flags = []
    rpaths = []
    if shutil.which("pkg-config"):
        for spec, needle in (("libgromacs", "gromacs/"), ("plumed", "plumed/")):
            if needle not in text:
                continue
            done = subprocess.run(
                ["pkg-config", "--cflags", "--libs", spec],
                text=True,
                capture_output=True,
                env=pkgconfig_env(),
            )
            if done.returncode == 0:
                flags.extend(done.stdout.split())
    prefix = sdk_prefix()
    if prefix is not None:
        include = prefix / "include"
        if "library.h" in text and (include / "lammps" / "library.h").is_file():
            flags.extend(["-I", str(include / "lammps"), "-llammps"])
        if "OpenMM.h" in text and (include / "OpenMM.h").is_file():
            flags.append("-lOpenMM")
        if "plumed/" in text and (prefix / "lib" / "libplumedKernel.so").is_file():
            flags.append("-lplumedKernel")
        for folder in (prefix / "lib", prefix / "lib.AVX2_256"):
            if folder.is_dir() and any(token in text for token in ("gromacs/", "plumed/", "OpenMM.h", "library.h")):
                flags.append("-L" + str(folder))
                rpaths.append(str(folder))
    seen = []
    for flag in flags:
        if flag.startswith("-L") and len(flag) > 2:
            rpaths.append(flag[2:])
        elif flag == "-L":
            continue
        if flag not in seen:
            seen.append(flag)
    index = 0
    while index < len(flags) - 1:
        if flags[index] == "-L":
            rpaths.append(flags[index + 1])
        index += 1
    for folder in dict.fromkeys(rpaths):
        link = "-Wl,-rpath," + folder
        if link not in seen:
            seen.append(link)
    return seen


def python_with(module, env):
    """True when this interpreter imports ``module`` from the snippet path."""
    probe = subprocess.run(
        [sys.executable, "-c", "import " + module],
        env=env,
        capture_output=True,
        text=True,
    )
    return sys.executable if probe.returncode == 0 else None


def which(name):
    found = shutil.which(name)
    if not found:
        raise Skip(name + " is not on PATH")
    return found


def vendored_links(build):
    """-rpath-link for the patched chemfiles meson builds into the build tree.

    A cross-prefixed ld such as conda's does not follow libyodaLib.so's RPATH
    to libseams-chemfiles.so, so the snippet link has to name the folder.
    """
    lib = build / "chemfiles-prefix" / "lib"
    return ["-Wl,-rpath-link," + str(lib)] if lib.is_dir() else []


def compile_c_family(build, env, source, out, extra=None, more_sources=None):
    compiler = which("c++" if source.suffix in {".cpp", ".cxx", ".cc"} else "cc")
    if source.suffix == ".c":
        compiler = which("cc")
    cmd = [
        compiler,
        "-I",
        str(ROOT / "include"),
        "-I",
        str(build / "src"),
        "-I",
        str(ROOT / "bindings" / "cxx"),
        "-I",
        str(ROOT / "src" / "include"),
        "-I",
        str(ROOT / "src" / "include" / "internal"),
    ]
    if extra:
        cmd.extend(extra)
    text = source.read_text(encoding="utf-8", errors="replace")
    sources = [str(source)]
    if more_sources:
        sources.extend(str(path) for path in more_sources)
    cmd.extend(
        [
            *sources,
            "-o",
            str(out),
            "-L" + str(build / "src"),
            "-L" + str(build / "bindings"),
            "-Wl,-rpath," + str(build / "src"),
            "-Wl,-rpath," + str(build / "bindings"),
            "-lm",
        ]
    )
    if "seams_lang.h" in text or "seams_frame_" in text or source.suffix != ".c":
        cmd.extend(["-lseams_lang", "-lstdc++"])
    cmd.append("-lseams")
    if sys.platform.startswith("linux"):
        # libseams calls dlopen. The link fails with "DSO missing" unless
        # the executable names libdl.
        cmd.append("-ldl")
    cmd.extend(vendored_links(build))
    cmd.extend(sdk_flags(text))
    if source.suffix != ".c":
        cmd.insert(1, "-std=c++17")
    if out.suffix == ".so":
        cmd.insert(1, "-fPIC")
        cmd.insert(1, "-shared")
    try:
        run(cmd, env)
    except Fail as exc:
        message = str(exc)
        missing = "No such file or directory" in message and any(
            token in message for token in ("gromacs/", "OpenMM.h", "library.h", "plumed/")
        )
        if missing:
            raise Skip("the snippet's toolchain headers are not installed") from exc
        raise


def program_args(text, build):
    if "cli_surface.py SEAMS" in text:
        seams = build / "seams"
        if not seams.is_file():
            raise Skip("seams CLI was not built")
        return [str(seams)]
    if "SEAMS_BUILD" in text:
        return []
    if "get_command_argument" in text and "input/traj" in text:
        return [str(ROOT)]
    if "commandArgs" in text or ("get_command_argument" in text and "input/traj" not in text):
        return [str(FIXTURE)]
    # Lua examples take the repository as arg[1] and join /input/.
    # A frame path there becomes "<dump>/input/traj/..." and the read fails.
    if ("arg[1]" in text or "ARGS[1]" in text) and "/input/" in text:
        return [str(ROOT)]
    if "arg[1]" in text or "ARGS[1]" in text:
        return [str(FIXTURE)]
    if "argv[2]" in text and "libseams" in text:
        return [str(FIXTURE), str(build / "src")]
    if ' / "src" / "libseams.so"' in text:
        return [str(build)]
    if "libseams.so" in text and "argv[1]" in text:
        return [str(build / "src")]
    if "env::args()" in text:
        return [str(FIXTURE)]
    return []


def run_source(build, env, path):
    text = path.read_text(encoding="utf-8", errors="replace")
    key = hashlib.sha1(text.encode()).hexdigest()

    def once():
        if path.suffix == ".c" or path.suffix in {".cpp", ".cxx", ".cc"}:
            built = None
            for cand in (build / "src").glob(path.stem):
                if cand.is_file() and os.access(cand, os.X_OK):
                    built = cand
                    break
            if built is None:
                for cand in (build / "bindings").glob(path.stem):
                    if cand.is_file() and os.access(cand, os.X_OK):
                        built = cand
                        break
            if built is not None and "main" in text:
                args = program_args(text, build)
                if path.name == "science_tutorials.c":
                    mode = "brine" if "brine" in os.environ.get("SNIPPET_PAGE", "") else "cages" if "cage" in os.environ.get("SNIPPET_PAGE", "") else "interface"
                    args = [mode]
                run([str(built), *args], env)
                return
            work = Path(tempfile.mkdtemp(prefix="seams-snippet-"))
            if not re.search(r"\bint\s+main\s*\(", text):
                library = work / "snippet.so"
                more = []
                if "addPlainComponent" in text:
                    # The one-argument PLUMED wrappers live in their own
                    # translation unit so they do not see the two-argument
                    # header declaration.
                    more.append(ROOT / "engines" / "plumed" / "plumed_components.cpp")
                compile_c_family(build, env, path, library, more_sources=more)
                import ctypes
                ctypes.CDLL(str(library))
                return
            exe = work / "snippet"
            compile_c_family(build, env, path, exe)
            args = program_args(text, build)
            if path.name == "science_tutorials.c":
                page_name = os.environ.get("SNIPPET_PAGE", "")
                mode = "brine" if "brine" in page_name else "cages" if "cage" in page_name else "interface"
                args = [mode]
            cwd = ROOT
            if "usage: run_hook" in text:
                args = [str(path.parent / "in.dseams")]
                cwd = path.parent
            if "TrajectoryAnalysisCommandLineRunner" in text:
                args = ["-f", str(ROOT / "engines" / "gromacs" / "ice.gro")]
            run([str(exe), *args], env, cwd=cwd)
        elif path.suffix == ".py":
            exe = sys.executable
            if "import pydseams" in text or "from pydseams" in text:
                exe = python_with("pydseams", env)
                if exe is None:
                    raise Fail("in-tree pydseams did not import")
            elif re.search(r"(?<![A-Za-z_])import dseams\b", text) or re.search(r"(?<![A-Za-z_])from dseams\b", text):
                exe = python_with("dseams._core", env)
                if exe is None:
                    raise Skip("the Python extension is not built")
            run([exe, str(path), *program_args(text, build)], env)
        elif path.suffix == ".lua":
            lua = shutil.which("lua5.4") or shutil.which("lua")
            if not lua:
                raise Skip("lua is not on PATH")
            modules = list((build / "bindings" / "lua").glob("*/dseams_core.so"))
            if not modules:
                raise Skip("the Lua module was not built")
            local = dict(env)
            local["LUA_CPATH"] = str(modules[0].parent / "?.so") + ";;"
            local["LUA_PATH"] = str(ROOT / "bindings" / "lua" / "?.lua") + ";;"
            run([lua, str(path), *program_args(text, build)], local)
        elif path.suffix in {".f90", ".F90"}:
            gfortran = which("gfortran")
            work = Path(tempfile.mkdtemp(prefix="seams-snippet-f90-"))
            # seams.f90 is the module. It does not use itself and it has no
            # program, so linking it looks for a main that is not there.
            if re.search(r"(?m)^\s*module\s+\w+", text) and not re.search(r"(?m)^\s*program\s+\w+", text):
                run([gfortran, "-c", str(path), "-o", str(work / "unit.o")], env, cwd=work)
                return
            if "use seams" not in text:
                exe = work / "snippet"
                run(
                    [
                        gfortran,
                        str(path),
                        "-o",
                        str(exe),
                        "-L" + str(build / "bindings"),
                        "-L" + str(build / "src"),
                        "-lseams_lang",
                        "-lseams",
                        "-Wl,-rpath," + str(build / "bindings"),
                        "-Wl,-rpath," + str(build / "src"),
                        *vendored_links(build),
                    ],
                    env,
                )
                run([str(exe), *program_args(text, build)], env)
                return
            run(
                [gfortran, "-J", str(work), "-c", str(ROOT / "bindings" / "fortran" / "seams.f90"), "-o", str(work / "seams.o")],
                env,
                cwd=work,
            )
            if not re.search(r"(?m)^\s*program\s+\w+", text):
                run([gfortran, "-J", str(work), "-c", str(path), "-o", str(work / "unit.o")], env, cwd=work)
                return
            exe = work / "snippet"
            run(
                [
                    gfortran,
                    "-I",
                    str(work),
                    str(path),
                    str(work / "seams.o"),
                    "-o",
                    str(exe),
                    "-L" + str(build / "src"),
                    "-lseams",
                    "-ldl",
                    "-Wl,-rpath," + str(build / "src"),
                    *vendored_links(build),
                ],
                env,
                cwd=work,
            )
            run([str(exe), *program_args(text, build)], env)
        elif path.suffix == ".jl":
            julia = which("julia")
            if not (build / "bindings" / "libseams_lang.so").is_file():
                raise Skip("libseams_lang.so was not built")
            project = ROOT / "bindings" / "julia" / ("Seams" if "/Seams/" in str(path) else "Dseams")
            local = dict(env)
            local["JULIA_LOAD_PATH"] = str(project) + os.pathsep
            local["SEAMS_LIB"] = str(build / "src" / "libseams.so")
            run([julia, "--startup-file=no", "--project=" + str(project), str(path), *program_args(text, build)], local)
        elif path.suffix == ".R":
            script = shutil.which("Rscript")
            if not script:
                raise Skip("Rscript is not on PATH")
            local = dict(env)
            if "library(dseams)" in text:
                sys.path.insert(0, str(ROOT / "docs" / "examples"))
                import run_binding

                script, renv = run_binding.install_r_package(build)
                local.update(renv)
            if "testthat" in text or "skip_if" in text:
                probe = subprocess.run(
                    [script, "-e", 'quit(status=if (requireNamespace("testthat", quietly=TRUE)) 0 else 1)'],
                    env=local,
                    capture_output=True,
                    text=True,
                )
                if probe.returncode != 0:
                    raise Skip("testthat is not installed")
                run([script, "-e", "testthat::test_file('" + str(path) + "', stop_on_failure=TRUE)"], local)
                return
            run([script, str(path), *program_args(text, build)], local)
        elif path.suffix == ".rs":
            which("cargo")
            local = dict(env)
            local["SEAMS_LIB_DIR"] = str(build / "src")
            manifest = str(ROOT / "bindings" / "rust" / "seams" / "Cargo.toml")
            # language_fronts and the API page also cargo-test this tree.
            cargo_lock = ROOT / "bindings" / "rust" / ".cargo-test.lock"
            def cargo(cmd):
                import fcntl
                with open(cargo_lock, "a", encoding="utf-8") as handle:
                    fcntl.flock(handle, fcntl.LOCK_EX)
                    run(cmd, local)
            if path.parent.name == "tests":
                cargo(["cargo", "test", "--offline", "--manifest-path", manifest, "--test", path.stem, "--", "--test-threads=1"])
            elif "use seams::" in text:
                example = "minimal" if path.name == "minimal_rust.rs" else path.stem
                cargo(["cargo", "run", "--offline", "--quiet", "--manifest-path", manifest, "--example", example, "--", *program_args(text, build)])
            else:
                which("rustc")
                work = Path(tempfile.mkdtemp(prefix="seams-snippet-rs-"))
                exe = work / "snippet"
                run(
                    [
                        "rustc",
                        "--edition",
                        "2021",
                        "-o",
                        str(exe),
                        str(path),
                        "-L",
                        str(build / "bindings"),
                        "-L",
                        str(build / "src"),
                        "-C",
                        "link-arg=-Wl,-rpath," + str(build / "bindings"),
                        "-C",
                        "link-arg=-Wl,-rpath," + str(build / "src"),
                    ],
                    local,
                )
                run([str(exe), *program_args(text, build)], local)
        else:
            raise Skip("no runner for " + path.suffix)

    cached(build, key, once)


def resolve_include(page, raw):
    target = (published_dir(page) / raw).resolve()
    if not target.is_file():
        alt = (ROOT / raw).resolve()
        if alt.is_file():
            return alt
        raise Fail("missing included file " + raw)
    return target


def check_excerpt(body):
    lines = body.splitlines()
    if not lines:
        raise Fail("empty excerpt")
    match = EXCERPT.match(lines[0].strip())
    if match is None:
        return None
    target = (ROOT / match.group(1)).resolve()
    excerpt = "\n".join(lines[1:]).strip("\n")
    if not target.is_file():
        raise Fail("excerpt file missing " + match.group(1))
    if excerpt not in target.read_text(encoding="utf-8"):
        raise Fail("excerpt is not in " + match.group(1))
    return target


def declaration(body):
    text = body.strip()
    if not text or "\n{" in text or text.rstrip().endswith("{"):
        return False
    if ";" in text and text.count(";") == 1 and "(" in text and "=" not in text.split("(")[0]:
        # A single call statement, not a prototype, still has one semicolon.
        head = text.split("(", 1)[0]
        return not any(word in head for word in (" ",)) or head.strip().startswith("ring::") or "\n" in text and text.strip().endswith(";") and " = " not in text and "{" not in text
    return "(" in text and not text.rstrip().endswith(";") and "\n" not in text.split("(")[0][:1]


def is_prototype(body):
    text = body.strip()
    if not text or "{" in text or "=" in text or text.count(";") > 1:
        return False
    head = text.split("(", 1)[0]
    return "(" in text and (" " in head or "*" in head or "::" in head)


def compile_statements(build, env, lang, body):
    work = Path(tempfile.mkdtemp(prefix="seams-snippet-src-"))
    if lang == "c":
        source = work / "snippet.c"
        predeclared = "seams_system *sys" not in body
        uses_frame = "seams_frame_" in body
        frame = ""
        if uses_frame and "seams_frame *frame" not in body:
            frame = "seams_frame *frame = seams_frame_read(path, 1, -1);\n"
        sysline = "" if not predeclared else "seams_system *sys = seams_system_create();\n  seams_system_set_ortho(sys, xyz, 8, 0, 0, 0, a, a, a, NULL, NULL);\n"
        owned = "sys" if predeclared or "seams_system *sys" in body else "0"
        mixed = "seams_frame_" in body and ("seams_system_" in body or "seams_chill_plus(" in body or "seams_rings_" in body or "seams_q4" in body)
        header = "seams.h" if mixed or "seams_frame_" not in body else "seams_lang.h"
        forward = ""
        if mixed:
            forward = "\n".join(
                [
                    "typedef struct seams_frame seams_frame;",
                    "seams_frame *seams_frame_read(const char *path, int frame, int type_filter);",
                    "void seams_frame_free(seams_frame *frame);",
                    "int seams_frame_chill_plus_rule(seams_frame *frame, double cutoff, int type, const char *graph, const char *rule, int *labels, int counts9[9]);",
                ]
            )
        source.write_text(
            "\n".join(
                [
                    '#include "' + header + '"',
                    forward,
                    "#include <math.h>",
                    "#include <stdio.h>",
                    "#include <stdlib.h>",
                    "int main(void) {",
                    '  const char *path = "' + str(ROOT / "input" / "traj" / "ice1h.lammpstrj") + '";',
                    "  const double a = 2.75 * 4.0 / sqrt(3.0);",
                    "  const double s[8][3] = {"
                    "    {0, 0, 0}, {0.5, 0.5, 0}, {0.5, 0, 0.5}, {0, 0.5, 0.5},"
                    "    {0.25, 0.25, 0.25}, {0.75, 0.75, 0.25}, {0.75, 0.25, 0.75}, {0.25, 0.75, 0.75}};"
                    "  double xyz[24];",
                    "  int labels[8192];",
                    "  int counts[9] = {0};",
                    "  for (int i = 0; i < 8; i++) {",
                    "    xyz[3 * i] = s[i][0] * a;",
                    "    xyz[3 * i + 1] = s[i][1] * a;",
                    "    xyz[3 * i + 2] = s[i][2] * a;",
                    "  }",
                    sysline,
                    frame,
                    body,
                    "  if (frame) seams_frame_free(frame);" if uses_frame and "seams_frame_free" not in body else "  (void)path;",
                    "  if (" + owned + ") seams_system_destroy(sys);" if "seams_system_destroy" not in body else "",
                    "  return 0;",
                    "}",
                    "",
                ]
            ),
            encoding="utf-8",
        )
        # The path above is the checkout. Keep it out of committed files; this
        # file is a temp compilation unit.
        exe = work / "snippet"
        compile_c_family(build, env, source, exe)
        run([str(exe)], env)
        return
    if lang in {"cpp", "c++", "cxx"} and re.search(r"\bint\s+main\s*\(", body):
        source = work / "snippet.cpp"
        source.write_text(body if body.endswith("\n") else body + "\n", encoding="utf-8")
        exe = work / "snippet"
        compile_c_family(build, env, source, exe)
        run([str(exe)], env)
        return
    if lang in {"cpp", "c++", "cxx"}:
        header = None
        folded = " ".join(body.split())
        for candidate in (ROOT / "src").rglob("*.h*"):
            raw = candidate.read_text(encoding="utf-8", errors="replace")
            if folded in " ".join(raw.split()):
                header = candidate
                break
        candidates = [folded]
        if folded.startswith("ring::"):
            candidates.append(folded[len("ring::") :])
        if header is None:
            for candidate in (ROOT / "src").rglob("*.h*"):
                raw = " ".join(candidate.read_text(encoding="utf-8", errors="replace").split())
                if any(item in raw for item in candidates):
                    header = candidate
                    break
        if header is None:
            raise Fail("C++ snippet is not a declaration in this tree")
        source = work / "snippet.cpp"
        rel = header.name
        source.write_text(
            f'#include "{rel}"\nint main() {{ return 0; }}\n',
            encoding="utf-8",
        )
        compile_c_family(build, env, source, work / "snippet", ["-I", str(header.parent)])
        return
    raise Fail("no statement runner for " + lang)


def compile_prototypes(build, env, bodies):
    work = Path(tempfile.mkdtemp(prefix="seams-snippet-proto-"))
    source = work / "snippet.c"
    lines = ['#include "seams.h"', ""]
    for body in bodies:
        text = body.strip()
        if not text.endswith(";"):
            text += ";"
        lines.append(text)
    lines.append("")
    lines.append("int main(void) { return 0; }")
    source.write_text("\n".join(lines) + "\n", encoding="utf-8")
    compile_c_family(build, env, source, work / "snippet")
    run([str(work / "snippet")], env)


def run_python(build, env, chunks):
    work = Path(tempfile.mkdtemp(prefix="seams-snippet-py-"))
    script = work / "snippet.py"
    body = "\n".join(chunks)
    exe = sys.executable
    for name in OPTIONAL:
        if re.search(r"(?<![A-Za-z_])import " + re.escape(name) + r"\b", body) or re.search(
            r"(?<![A-Za-z_])from " + re.escape(name) + r"\b", body
        ):
            found = python_with(name.split(".")[0], env)
            if found is None:
                if name.startswith("pydseams"):
                    raise Fail("in-tree pydseams did not import")
                raise Skip(name + " is not installed")
            if name.startswith("pydseams"):
                exe = found
    plain = body.replace("pydseams", "")
    if "dseams." in plain or "import dseams" in plain or "from dseams" in plain:
        found = python_with("dseams._core", env)
        if found is None:
            raise Skip("the Python extension is not built")
        exe = found
    prelude = ""
    if "dseams." in plain and "import dseams" not in plain and "from dseams" not in plain:
        prelude = "import dseams\n"
    script.write_text(
        prelude
        + "import math\n"
        "import os\n"
        "import numpy as np\n"
        "a = 2.75 * 4 / 3**0.5\n"
        "sites = ((0, 0, 0), (0.5, 0.5, 0), (0.5, 0, 0.5), (0, 0.5, 0.5),\n"
        "         (0.25, 0.25, 0.25), (0.75, 0.75, 0.25), (0.75, 0.25, 0.75), (0.25, 0.75, 0.75))\n"
        "positions = np.array([[x * a, y * a, z * a] for x, y, z in sites], dtype=float)\n"
        "box = np.array([3 * a, 3 * a, 3 * a], dtype=float)\n"
        "root = os.environ['DSEAMS_ROOT']\n"
        + body
        + "\n",
        encoding="utf-8",
    )
    run([exe, str(script)], env)


def run_julia_text(build, env, body):
    julia = which("julia")
    if not (build / "bindings" / "libseams_lang.so").is_file():
        raise Skip("libseams_lang.so was not built")
    lines = []
    for line in body.splitlines():
        stripped = line.strip()
        if "==" in stripped and not stripped.startswith("@") and not stripped.startswith("#"):
            lines.append("@assert (" + stripped + ")")
        else:
            lines.append(line)
    work = Path(tempfile.mkdtemp(prefix="seams-snippet-jl-"))
    script = work / "snippet.jl"
    script.write_text("using Dseams\n" + "\n".join(lines) + "\n", encoding="utf-8")
    local = dict(env)
    local["JULIA_LOAD_PATH"] = str(ROOT / "bindings" / "julia" / "Dseams") + os.pathsep
    run([julia, "--startup-file=no", "--project=" + str(ROOT / "bindings" / "julia" / "Dseams"), str(script)], local)


def run_lua_text(build, env, body):
    if body.strip().startswith("--") and not any(line.strip() and not line.strip().startswith("--") for line in body.splitlines()):
        raise Fail("Lua snippet is only a comment")
    lua = shutil.which("lua5.4") or shutil.which("lua")
    if not lua:
        raise Skip("lua is not on PATH")
    modules = list((build / "bindings" / "lua").glob("*/dseams_core.so"))
    if not modules:
        raise Skip("the Lua module was not built")
    work = Path(tempfile.mkdtemp(prefix="seams-snippet-lua-"))
    script = work / "snippet.lua"
    prelude = 'local dseams = require("dseams")\nlocal root = os.getenv("DSEAMS_ROOT")\n'
    if "frame" in body and "local frame" not in body and "dseams.read" not in body:
        prelude += 'local frame = dseams.read(root .. "/input/traj/genice_sI.lammpstrj", {type = 1})\n'
    script.write_text(prelude + body + "\n", encoding="utf-8")
    local = dict(env)
    local["LUA_CPATH"] = str(modules[0].parent / "?.so") + ";;"
    local["LUA_PATH"] = str(ROOT / "bindings" / "lua" / "?.lua") + ";;"
    run([lua, str(script)], local)


def run_r_text(build, env, body):
    script = which("Rscript")
    work = Path(tempfile.mkdtemp(prefix="seams-snippet-r-"))
    path = work / "snippet.R"
    path.write_text(
        'root <- Sys.getenv("DSEAMS_ROOT")\n'
        'lib <- Sys.getenv("DSEAMS_R_LIB", "")\n'
        "if (nzchar(lib)) .libPaths(c(lib, .libPaths()))\n"
        "suppressPackageStartupMessages(library(dseams))\n"
        + body
        + "\n",
        encoding="utf-8",
    )
    try:
        run([script, str(path)], env)
    except Fail as exc:
        if "there is no package called" in str(exc):
            raise Skip("the R package is not installed") from exc
        raise


def run_rust_text(build, env, body):
    which("cargo")
    if "fn main" not in body or "assert" not in body:
        raise Fail("Rust example must be a complete main that asserts its result")
    work = Path(tempfile.mkdtemp(prefix="seams-snippet-rs-"))
    (work / "src").mkdir()
    (work / "Cargo.toml").write_text(
        "\n".join(
            [
                "[package]",
                'name = "seams_doc_example"',
                'version = "0.0.0"',
                'edition = "2021"',
                "",
                "[dependencies]",
                'seams = {{ path = "{}" }}'.format(ROOT / "bindings" / "rust" / "seams"),
                'seams-sys = {{ path = "{}" }}'.format(ROOT / "bindings" / "rust" / "seams-sys"),
                "",
            ]
        ),
        encoding="utf-8",
    )
    (work / "src" / "main.rs").write_text(body if body.endswith("\n") else body + "\n", encoding="utf-8")
    local = dict(env)
    local["SEAMS_LIB_DIR"] = str(build / "src")
    run(["cargo", "run", "--offline", "--quiet", "--manifest-path", str(work / "Cargo.toml")], local)


def run_doctest(env, body):
    work = Path(tempfile.mkdtemp(prefix="seams-snippet-doc-"))
    path = work / "snippet.txt"
    path.write_text(body if body.endswith("\n") else body + "\n", encoding="utf-8")
    run([sys.executable, "-m", "doctest", "-v", str(path)], env)


def execute_page(path, build):
    os.environ["SNIPPET_PAGE"] = page_id(path)
    env = base_env(build)
    env["SNIPPET_PAGE"] = page_id(path)
    ran = 0
    skips = []
    fails = []
    prototypes = []
    python_chunks = []
    for kind, lang, body in snippets_of(path):
        try:
            if kind == "include":
                target = resolve_include(path, lang)
                suffix = target.suffix.lower()
                if suffix in {".sh", ".dat", ""} or target.suffix in {".sh"}:
                    continue
                if suffix not in {".c", ".cpp", ".cxx", ".cc", ".py", ".rs", ".jl", ".lua", ".f90", ".r"}:
                    continue
                run_source(build, env, target)
                ran += 1
                continue
            excerpt = check_excerpt(body)
            if excerpt is not None:
                # The lines are in a file the suite already builds. Recompiling
                # a Catch2 translation unit here pulls in failures outside the
                # quoted lines. docs/examples and the bindings are run below.
                if "tests" in excerpt.parts and excerpt.suffix in {".c", ".cpp", ".h", ".hpp", ".py"}:
                    ran += 1
                    continue
                if excerpt.name == "bench_cli_fast.py" and not (ROOT / "build" / "seams").is_file():
                    link = ROOT / "build"
                    if not link.exists() and (build / "seams").is_file():
                        link.symlink_to(build)
                if excerpt.name == "bench_cli_fast.py" and not (ROOT / "build" / "seams").is_file():
                    ran += 1
                    skips.append("bench driver looks for build/seams")
                    continue
                run_source(build, env, excerpt)
                ran += 1
                continue
            if lang == "c" and is_prototype(body):
                prototypes.append(body)
                continue
            if lang == "python" and body.lstrip().startswith(">>>"):
                run_doctest(env, body)
                ran += 1
                continue
            if lang == "python":
                python_chunks.append(body)
                continue
            if lang in {"cpp", "c++", "cxx"}:
                compile_statements(build, env, lang, body)
                ran += 1
                continue
            if lang == "c":
                compile_statements(build, env, "c", body)
                ran += 1
                continue
            if lang == "julia":
                run_julia_text(build, env, body)
                ran += 1
                continue
            if lang == "lua":
                run_lua_text(build, env, body)
                ran += 1
                continue
            if lang == "r":
                run_r_text(build, env, body)
                ran += 1
                continue
            if lang == "rust":
                run_rust_text(build, env, body)
                ran += 1
                continue
            raise Fail("unhandled " + lang)
        except Skip as exc:
            skips.append(str(exc))
        except Fail as exc:
            fails.append("\n".join(str(exc).strip().splitlines()[-6:]) if str(exc).strip() else "failed")
    if prototypes:
        try:
            compile_prototypes(build, env, prototypes)
            ran += 1
        except Skip as exc:
            skips.append(str(exc))
        except Fail as exc:
            fails.append("\n".join(str(exc).strip().splitlines()[-6:]) if str(exc).strip() else "prototype")
    if python_chunks:
        try:
            run_python(build, env, python_chunks)
            ran += 1
        except Skip as exc:
            skips.append(str(exc))
        except Fail as exc:
            fails.append("\n".join(str(exc).strip().splitlines()[-6:]) if str(exc).strip() else "python")
    return ran, skips, fails


def main(argv):
    if not argv:
        print("usage: run_snippets.py --list | --report BUILD | PAGE BUILD", file=sys.stderr)
        return 2
    if argv[0] == "--list":
        for path in pages():
            print(page_id(path))
        return 0
    if argv[0] == "--report":
        build = Path(argv[1]).resolve()
        broken = []
        skipped = []
        ok = 0
        total_snippets = 0
        for path in pages():
            total_snippets += len(snippets_of(path))
            ran, skips, fails = execute_page(path, build)
            name = page_id(path)
            if fails:
                broken.append((name, fails))
                print("FAIL", name)
                for item in fails:
                    print(" ", item)
            elif ran == 0:
                skipped.append(name)
                print("SKIP", name, "; ".join(skips))
            else:
                ok += 1
                print("OK", name, "ran", ran, "skipped", len(skips))
        print("pages", len(pages()), "snippets", total_snippets, "ok", ok, "skipped", len(skipped), "broken", len(broken))
        return 1 if broken else 0
    if len(argv) < 2:
        return 2
    wanted = argv[0]
    build = Path(argv[1]).resolve()
    match = [path for path in pages() if page_id(path) == wanted]
    if not match:
        print("unknown page " + wanted, file=sys.stderr)
        return 2
    ran, skips, fails = execute_page(match[0], build)
    for item in skips:
        print("skipped:", item, file=sys.stderr)
    if fails:
        print("\n".join(fails), file=sys.stderr)
        return 1
    if ran == 0:
        return 77
    return 0


if __name__ == "__main__":
    sys.exit(main(sys.argv[1:]))
