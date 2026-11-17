# Copyright (c) d-SEAMS developers
"""Compiler the snippet tests must use.

Meson puts the C and C++ compilers that built the library, and the OpenMP
link arguments those compilers need, in the test environment. ``cc`` on
``PATH`` is a different toolchain: in the pixi jobs it is conda's gcc, whose
sysroot cannot link a library built by the job's host compiler.
"""

import os
import shlex
import subprocess


def _need(name):
    value = os.environ.get(name, "").strip()
    if not value:
        raise SystemExit(
            name + " is unset. The test environment must name the compiler meson used."
        )
    return value


def c_compiler():
    return _need("CC")


def cxx_compiler():
    return _need("CXX")


def fortran_compiler():
    value = os.environ.get("FC", "").strip()
    if not value:
        raise SystemExit("FC is unset. The test environment must name meson's Fortran compiler.")
    return value


def gfortran_link():
    """``-lgfortran`` plus the directory ``FC`` actually ships it in.

    clang does not search the GCC library directory, so a bare
    ``-lgfortran`` misses the runtime next to ``gfortran-15``.
    """
    found = subprocess.run(
        [fortran_compiler(), "-print-file-name=libgfortran.so"],
        text=True,
        capture_output=True,
        check=False,
    ).stdout.strip()
    if found.startswith("/") and os.path.isfile(found):
        folder = os.path.dirname(found)
        return ["-L" + folder, "-lgfortran", "-Wl,-rpath," + folder]
    return ["-lgfortran"]


def openmp_link():
    raw = os.environ.get("SEAMS_OPENMP_LINK", "").strip()
    return shlex.split(raw) if raw else []


def link_command(compiler, args):
    """``compiler`` plus ``args``, with the OpenMP link line at the end."""
    extra = [flag for flag in openmp_link() if flag not in args]
    return [compiler, *args, *extra]
