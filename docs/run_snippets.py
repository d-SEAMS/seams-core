#!/usr/bin/env python3
# Copyright (c) d-SEAMS developers
"""Run every documented Python, C and C++ snippet.

One Meson test is docs_snippet_<page>. A missing compiler exits 77. A
snippet that does not match the current API fails.

Python snippets import ``dseams`` and ``pydseams`` from the installed
pydseamslib package, built from d-SEAMS/PydSEAMSlib. Nothing in this tree
is put on the import path for them. When the package is missing the page
fails and names the pip line that installs it.

Rust, Julia, Lua, Fortran and R examples live in the binding repositories.
A page here that still carries one fails with the repository it belongs in.
"""

from __future__ import annotations

import hashlib
import os
import re
import shlex
import shutil
import subprocess
import sys
import tempfile
import textwrap
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "scripts"))
from shared_library import built_library  # noqa: E402
from project_version import project_version  # noqa: E402
DOCS = ROOT / "docs"
FIXTURE = ROOT / "input" / "traj" / "genice_sI.lammpstrj"
LANGS = {"python", "py", "c", "cpp", "c++", "cxx", "rust", "julia", "lua", "fortran", "f90", "r"}
# Fences and includes in these languages are found so that a page carrying
# one fails; the examples belong to the binding repositories.
MOVED = {
    "rust": "https://github.com/d-SEAMS/dseams-rs",
    "julia": "https://github.com/d-SEAMS/Dseams.jl",
    "lua": "https://github.com/d-SEAMS/yodaStruct",
    "fortran": "https://github.com/d-SEAMS/dseams-fortran",
    "f90": "https://github.com/d-SEAMS/dseams-fortran",
    "r": "https://github.com/d-SEAMS/dseams-r",
}
MOVED_SUFFIX = {".rs": "rust", ".jl": "julia", ".lua": "lua", ".f90": "fortran", ".r": "r"}
ORG_SRC = re.compile(r"#\+begin_src\s+([^\s]+)[^\n]*\n(.*?)#\+end_src", re.S | re.I)
RST_DIR = re.compile(
    r"^\.\.\s+(literalinclude|code-block|code|doctest)::[ \t]*([^\n]*)\n((?:[ \t].*\n|\n)*)",
    re.M,
)
MD_INCLUDE = re.compile(r"```\{literalinclude\}\s+(\S+)")
# A language tag is required. A bare closing ``` would otherwise swallow the next fence.
MD_FENCE = re.compile(r"```([A-Za-z0-9_+-]+)\n(.*?)```", re.S)
EXCERPT = re.compile(r"^(?:#|//|!)\s*excerpt:\s*(\S+)\s*$")
OPTIONAL = ("ase", "solvis")
BINDING = ("dseams", "pydseams")


def install_hint():
    return "install it with: pip install pydseamslib=={}".format(project_version())


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
        str(build / "chemfiles-prefix" / "lib"),
        str(build / "subprojects" / "readcon-core"),
        str(build / "subprojects" / "readcon-db"),
        str(build / "subprojects" / "readcon-core" / "cargo-target" / "release"),
        str(build / "subprojects" / "readcon-db" / "cargo-target" / "release"),
    ]
    current = os.environ.get("LD_LIBRARY_PATH")
    if current:
        parts.append(current)
    return os.pathsep.join(parts)


def base_env(build):
    env = os.environ.copy()
    env.pop("DSEAMS_DOC_BUILD", None)
    # The engines and libseams carry an rpath. A leftover loader path
    # must not be what makes a snippet resolve.
    env.pop("LD_LIBRARY_PATH", None)
    env["DSEAMS_ROOT"] = str(ROOT)
    env["LD_LIBRARY_PATH"] = library_path(build)
    env["OMP_NUM_THREADS"] = "1"
    lang = built_library(build / "bindings", "seams_lang", "SEAMS_LANG_LIB")
    if lang is not None:
        env["SEAMS_LANG_LIB"] = str(lang)
    env["SEAMS_LIB_DIR"] = str(build / "src")
    env["SEAMS_LIBDIR"] = str(build / "src")
    lib = built_library(build / "src", "seams", "SEAMS_LIB")
    if lib is not None:
        env["SEAMS_LIB"] = str(lib)
    env["SEAMS_BUILD"] = str(build)
    oracle = build / "bindings" / "binding_oracle"
    if oracle.is_file():
        # The binding test suites that pages run compare against it.
        env["SEAMS_CPP_ORACLE"] = str(oracle)
        env["SEAMS_SOURCE_ROOT"] = str(ROOT)
    env["PYTHONPATH"] = os.pathsep.join(
        [
            str(DOCS / "examples"),
            str(ROOT / "adapters" / "python"),
            str(ROOT / "scripts"),
            env.get("PYTHONPATH", ""),
        ]
    )
    return env


def tsan_host(env):
    """Environment for a binary we just compiled against the sanitized libseams.

    That binary is not itself instrumented, so the thread runtime is not
    in its link. Preloading it resolves __tsan_*. An executable the build
    already linked with -fsanitize=thread must not use this.
    """
    runtime = os.environ.get("TSAN_RUNTIME", "")
    if not runtime or not os.path.isfile(runtime):
        return env
    out = dict(env)
    prior = out.get("LD_PRELOAD", "")
    if runtime not in prior.split(":"):
        out["LD_PRELOAD"] = runtime if not prior else runtime + ":" + prior
    return out


def run(cmd, env, cwd=ROOT):
    # A nested interpreter does not go through the asan-path wrapper, so
    # it does not map the thread runtime unless this environment does.
    if Path(cmd[0]).name.startswith("python"):
        env = tsan_host(env)
    done = subprocess.run(cmd, env=env, cwd=cwd, text=True, capture_output=True)
    if done.returncode != 0:
        tail = (done.stderr + done.stdout)[-1200:]
        if done.returncode == 77:
            raise Skip(tail.strip().splitlines()[-1] if tail.strip() else "skipped")
        if "PackageNotFoundError" in tail or "No package metadata was found" in tail:
            raise Skip("an optional package named by the snippet is not installed")
        if "No module named 'openmm'" in tail or 'No module named "openmm"' in tail:
            raise Skip("openmm is not installed")
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


def compiler_includes():
    folders = []
    for name in ("CPATH", "CPLUS_INCLUDE_PATH", "C_INCLUDE_PATH"):
        folders.extend(Path(part) for part in os.environ.get(name, "").split(os.pathsep) if part)
    return list(dict.fromkeys(folders))


def library_file(name):
    """The shared library the linker would take for ``name``, or None."""
    compiler = os.environ.get("CC", "").strip()
    if compiler:
        found = subprocess.run(
            [*shlex.split(compiler), "-print-file-name=" + name],
            text=True,
            capture_output=True,
        ).stdout.strip()
        if os.path.isabs(found) and os.path.isfile(found):
            return found
    # clang links from LIBRARY_PATH but -print-file-name does not search it.
    for part in os.environ.get("LIBRARY_PATH", "").split(os.pathsep):
        if part and os.path.isfile(os.path.join(part, name)):
            return os.path.join(part, name)
    return None


def sdk_flags(text):
    flags = []
    rpaths = []
    lammps_pc = False
    if shutil.which("pkg-config"):
        # A LAMMPS CMake install names its pkg-config file liblammps.pc.
        specs = (("libgromacs", "gromacs/"), ("plumed", "plumed/"), ("liblammps", "library.h"), ("lammps", "library.h"))
        for spec, needle in specs:
            if needle not in text or (needle == "library.h" and lammps_pc):
                continue
            done = subprocess.run(
                ["pkg-config", "--cflags", "--libs", spec],
                text=True,
                capture_output=True,
                env=pkgconfig_env(),
            )
            if done.returncode == 0:
                flags.extend(done.stdout.split())
                if needle == "library.h":
                    lammps_pc = True
                    where = {}
                    for name in ("includedir", "libdir"):
                        where[name] = subprocess.run(
                            ["pkg-config", "--variable=" + name, spec],
                            text=True,
                            capture_output=True,
                            env=pkgconfig_env(),
                        ).stdout.strip()
                    # Snippets include "library.h", which is under includedir/lammps.
                    if where["includedir"]:
                        flags.extend(["-I", str(Path(where["includedir"]) / "lammps")])
                    # The environment's lib may hold another liblammps.so.0, so
                    # this one comes first on the runpath.
                    if where["libdir"]:
                        rpaths.append(where["libdir"])
    prefix = sdk_prefix()
    if prefix is not None:
        include = prefix / "include"
        if "library.h" in text and not lammps_pc and (include / "lammps" / "library.h").is_file():
            flags.extend(["-I", str(include / "lammps"), "-llammps"])
            if library_file("libmpi.so"):
                flags.append("-lmpi")
        if "OpenMM.h" in text and (include / "OpenMM.h").is_file():
            flags.append("-lOpenMM")
        if "plumed/" in text and (prefix / "lib" / "libplumedKernel.so").is_file():
            flags.append("-lplumedKernel")
        for folder in (prefix / "lib", prefix / "lib.AVX2_256"):
            if folder.is_dir() and any(token in text for token in ("gromacs/", "plumed/", "OpenMM.h", "library.h")):
                flags.append("-L" + str(folder))
                rpaths.append(str(folder))
    if "plumed/" in text and "-lplumedKernel" not in flags:
        # The kernel symbols, such as PLMD::OpenMP, are in libplumedKernel.
        kernel = library_file("libplumedKernel.so")
        if kernel:
            flags.extend(["-L" + os.path.dirname(kernel), "-lplumedKernel"])
    if "OpenMM.h" in text and "-lOpenMM" not in flags:
        # Module systems put OpenMM on the compiler search paths, as
        # engines/meson.build finds it, and leave LD_LIBRARY_PATH alone.
        for include in compiler_includes():
            lib = include.parent / "lib"
            if (include / "OpenMM.h").is_file() and (lib / "libOpenMM.so").is_file():
                flags.extend(["-L" + str(lib), "-lOpenMM"])
                break
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
    # pkg-config leaves out -L for folders already on LIBRARY_PATH, so the
    # library links and then has no runpath.
    for flag in flags:
        if not flag.startswith("-l"):
            continue
        found = library_file("lib{}.so".format(flag[2:]))
        if found:
            rpaths.append(os.path.dirname(os.path.realpath(found)))
    if "-fopenmp" in flags:
        # clang's -fopenmp links libomp from its own folder, which the loader does not search.
        for runtime in ("libomp.so", "libgomp.so"):
            found = library_file(runtime)
            if found:
                rpaths.append(os.path.dirname(os.path.realpath(found)))
                break
    for folder in dict.fromkeys(rpaths):
        link = "-Wl,-rpath," + folder
        if link not in seen:
            seen.append(link)
    return seen


def python_with(module, env):
    """True when this interpreter imports ``module`` from the snippet path."""
    probe = subprocess.run(
        [sys.executable, "-c", "import " + module],
        env=tsan_host(env),
        capture_output=True,
        text=True,
    )
    return sys.executable if probe.returncode == 0 else None


def imports_of(text, names):
    """The names from ``names`` that ``text`` imports.

    >>> imports_of("import dseams\\nfrom pydseams import read", BINDING)
    ['dseams', 'pydseams']
    >>> imports_of("from dseams_insitu.engine import InSitu", BINDING)
    []
    """
    found = []
    for name in names:
        pattern = r"(?<![A-Za-z_.])(?:import|from)\s+" + re.escape(name) + r"(?![A-Za-z0-9_])"
        if re.search(pattern, text):
            found.append(name)
    return found


def binding_python(text, env):
    """The interpreter for a snippet, failing when it needs pydseamslib and it is absent."""
    for name in imports_of(text, BINDING):
        if python_with(name, env) is None:
            raise Fail("{} is not importable; {}".format(name, install_hint()))
    return sys.executable


def which(name):
    found = shutil.which(name)
    if not found:
        raise Skip(name + " is not on PATH")
    return found


def vendored_links(build):
    """Link path for the patched chemfiles meson builds into the build tree.

    A cross-prefixed ld such as conda's does not follow libyodaLib.so's RPATH
    to libseams-chemfiles.so, so the snippet link has to name the folder.
    -rpath-link is an ELF option. Apple ld rejects it and takes -rpath.
    """
    lib = build / "chemfiles-prefix" / "lib"
    if not lib.is_dir():
        return []
    if sys.platform == "darwin":
        return ["-L" + str(lib), "-Wl,-rpath," + str(lib)]
    return ["-Wl,-rpath-link," + str(lib)]


def compile_c_family(build, env, source, out, extra=None, more_sources=None):
    from toolchain import c_compiler, cxx_compiler, link_command

    compiler = cxx_compiler() if source.suffix in {".cpp", ".cxx", ".cc"} else c_compiler()
    if source.suffix == ".c":
        compiler = c_compiler()
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
        # The project default; seams.hpp uses std::span.
        cmd.insert(1, "-std=c++20")
    if out.suffix == ".so":
        cmd.insert(1, "-fPIC")
        cmd.insert(1, "-shared")
    cmd = link_command(cmd[0], cmd[1:])
    try:
        run(cmd, env)
    except Fail as exc:
        message = str(exc)
        missing = any(
            token in message for token in ("gromacs/", "OpenMM.h", "library.h", "plumed/")
        ) and any(
            phrase in message
            for phrase in ("No such file or directory", "file not found", "fatal error")
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
    if "<repository root>" in text:
        return [str(ROOT)]
    if 'Path(sys.argv[1]) if len(sys.argv) > 1 else root / "build"\n' in text:
        return [str(build)]
    if "input/traj" in text and any(token in text for token in ("get_command_argument", "arg[1]", "ARGS[1]")):
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


def moved_message(lang, what):
    return "{} is {} code; it belongs in {} and this page links there".format(what, lang, MOVED[lang])


def run_source(build, env, path):
    moved = MOVED_SUFFIX.get(path.suffix.lower())
    if moved:
        raise Fail(moved_message(moved, path.name))
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
            exe = binding_python(text, env)
            run([exe, str(path), *program_args(text, build)], env)
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


C_TYPES = r"(?:const\s+)?(?:unsigned\s+)?(?:int|size_t|double|char|float|seams_\w+)\b"


def c_declared(body):
    """Names a C fragment declares itself.

    >>> sorted(c_declared("size_t n_cages = 0, n_atoms = 0;\\nseams_rings_destroy(rings);"))
    ['n_atoms', 'n_cages']
    >>> sorted(c_declared("int *off = malloc(n * sizeof *off), *cols = NULL;"))
    ['cols', 'off']
    """
    names = set()
    for stmt in re.split(r"[;{}]", body):
        match = re.match(r"\s*" + C_TYPES + r"(.*)", stmt, re.S)
        if match is None:
            continue
        depth = 0
        outside = []
        for ch in match.group(1):
            if ch == "(":
                depth += 1
            elif ch == ")":
                depth -= 1
            elif depth == 0:
                outside.append(ch)
        for part in "".join(outside).split(","):
            name = re.match(r"[\s*]*(?:const\s+)?\**\s*([A-Za-z_]\w*)", part)
            if name:
                names.add(name.group(1))
    return names


# What a header example may use without declaring: an 8-atom cubic ice cell
# in sys, buffers of the atom count, and handles built from sys.
C_BUFFERS = [
    ("n", "  size_t n = 8;"),
    ("xyz", "  double *xyz = doc_cell;"),
    ("labels", "  int labels_buf[8192] = {0};\n  int *labels = labels_buf;"),
    # Ring size is the histogram index. seams_rings_histogram(..., 13)
    # writes counts[0] through counts[12]. CHILL+ uses the first 9.
    ("counts", "  int counts[13] = {0};"),
    ("box9", "  const double box9[9] = {0, 0, 0, a, a, a, 0, 0, 0};"),
    ("frame0", "  const double *frame0 = doc_cell;"),
    ("frames", "  const double *frames[2] = {doc_cell, doc_cell};"),
    ("nframes", "  size_t nframes = 2;"),
] + [(name, "  int {}[8192] = {{0}};".format(name)) for name in ("ice", "per_atom", "kinds", "ids", "types", "sizes", "symmetry")] + [
    (name, "  double {}[8192] = {{0}};".format(name)) for name in ("q4", "q6", "q6bar")
]
C_HANDLES = [
    ("nb", "seams_neighbours", "seams_neighbours_build(sys, 3.5, 0, SEAMS_GRAPH_CUTOFF, &nb)", "seams_neighbours_destroy"),
    ("rings", "seams_rings", "seams_rings_find(sys, 3.5, 0, SEAMS_GRAPH_CUTOFF, 6, &rings)", "seams_rings_destroy"),
    ("cages", "seams_cage_set", "seams_cages_find(sys, 3.5, 0, SEAMS_GRAPH_CUTOFF, 0, 0, &cages)", "seams_cages_destroy"),
    ("fp", "seams_topo", "seams_fingerprint(sys, 3.5, 0, 4, SEAMS_GRAPH_CUTOFF, 2, 0, &fp)", "seams_fingerprint_destroy"),
    ("front", "seams_brine", "seams_brine_analyse(sys, front_ice, 8, 1, 2, 0, 2, &front)", "seams_brine_destroy"),
]


def c_context(body, predeclared):
    """Declarations and cleanup for the names a C fragment uses but does not declare."""
    declared = c_declared(body)

    def wanted(name):
        return name not in declared and re.search(r"\b" + name + r"\b", body) is not None

    context = [line for name, line in C_BUFFERS if wanted(name)]
    cleanup = []
    for name, ctype, make, destroy in C_HANDLES:
        if not predeclared or not wanted(name):
            continue
        if name == "front":
            context.append("  int front_ice[8] = {0};")
        context.append("  {} *{} = NULL;".format(ctype, name))
        make_call = make.split("(", 1)[0] + "("
        if make_call not in body:
            context.append("  (void){};".format(make))
        if destroy + "(" + name not in body:
            cleanup.append("  {}({});".format(destroy, name))
    return context, cleanup


def compile_statements(build, env, lang, body):
    work = Path(tempfile.mkdtemp(prefix="seams-snippet-src-"))
    if lang == "c":
        source = work / "snippet.c"
        uses_frame = "seams_frame_" in body
        uses_sys = "seams_system_" in body or re.search(r"\bsys\b", body) is not None
        declares_sys = "seams_system *sys" in body
        wants_handle = any(
            re.search(r"\b" + name + r"\b", body) for name, _ctype, _make, _destroy in C_HANDLES
        )
        frame = ""
        if uses_frame and "seams_frame *frame" not in body:
            frame = "seams_frame *frame = seams_frame_read(path, 1, -1);\n"
        # A frame-only call does not need a system handle. Injecting one
        # pulls seams.h symbols into a unit that includes only seams_lang.h.
        # A neighbour, ring, cage, or fingerprint handle is built from sys,
        # so that example still gets one.
        needs_sys = uses_sys or wants_handle
        if needs_sys and not declares_sys:
            sysline = (
                "seams_system *sys = seams_system_create();\n"
                "  seams_system_set_ortho(sys, doc_cell, 8, 0, 0, 0, a, a, a, NULL, NULL);\n"
            )
            owned = "sys"
        elif needs_sys:
            sysline = ""
            owned = "sys"
        else:
            sysline = ""
            owned = "0"
        predeclared = owned == "sys" and not declares_sys
        mixed = uses_frame and (
            needs_sys or "seams_chill_plus(" in body or "seams_rings_" in body or "seams_q4" in body
        )
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
        context, cleanup = c_context(body, predeclared)
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
                    "  double doc_cell[24];",
                    "  for (int i = 0; i < 8; i++) {",
                    "    doc_cell[3 * i] = s[i][0] * a;",
                    "    doc_cell[3 * i + 1] = s[i][1] * a;",
                    "    doc_cell[3 * i + 2] = s[i][2] * a;",
                    "  }",
                    sysline,
                    frame,
                    *context,
                    body,
                    *cleanup,
                    "  if (frame) seams_frame_free(frame);" if uses_frame and "seams_frame_free" not in body else "  (void)path;",
                    "  if (sys) seams_system_destroy(sys);" if owned == "sys" and "seams_system_destroy" not in body else "",
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
    if lang in {"cpp", "c++", "cxx"} and re.search(r"\bseams::(System|Error|Status)\b", body):
        # Statements on the seams.hpp wrapper run inside a main.
        source = work / "snippet.cpp"
        source.write_text('#include "seams.hpp"\nint main() {\n' + body + "\n  return 0;\n}\n", encoding="utf-8")
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
    prelude = ""
    plain = body.replace("pydseams", "")
    if re.search(r"(?<![A-Za-z_])dseams\.", plain) and not imports_of(plain, ("dseams",)):
        prelude = "import dseams\n"
    exe = binding_python(prelude + body, env)
    for name in imports_of(body, OPTIONAL):
        if python_with(name, env) is None:
            raise Skip(name + " is not installed")
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


def run_doctest(env, body):
    binding_python(body, env)
    for name in imports_of(body, OPTIONAL):
        if python_with(name, env) is None:
            raise Skip(name + " is not installed")
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
                if suffix in MOVED_SUFFIX:
                    raise Fail(moved_message(MOVED_SUFFIX[suffix], lang))
                run_source(build, env, target)
                ran += 1
                continue
            if lang in MOVED:
                raise Fail(moved_message(lang, "a " + lang + " block"))
            excerpt = check_excerpt(body)
            if excerpt is not None:
                # The lines are in a file the suite already builds. Recompiling
                # a Catch2 translation unit here pulls in failures outside the
                # quoted lines. docs/examples and the bindings are run below.
                if "tests" in excerpt.parts and excerpt.suffix in {".c", ".cpp", ".h", ".hpp", ".py"}:
                    ran += 1
                    continue
                if excerpt.name == "bench_cli_fast.py":
                    # check_excerpt found the quoted line. The driver is a
                    # timing run that rewrites markdown/cli-fast-paths.md,
                    # so a snippet test must not execute it.
                    ran += 1
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
