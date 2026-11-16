#!/usr/bin/env python3
"""LAMMPS launcher from the environment.

LMP is a command line, not only a path. It may already carry an MPI
launcher:

    export LMP="mpirun -np 8 lmp"

LMP_NP prefixes "mpirun -np LMP_NP" when LMP does not already start
with mpirun, mpiexec or srun. LMP_NP=all uses every processor of the
machine, the count of nproc --all. Plain nproc is not used: it honours
OMP_NUM_THREADS, which the drivers set to 1, so it would return 1.

    export LMP=lmp
    export LMP_NP=$(nproc --all)

LMP_NP unset, empty or 1 runs LMP as given. LMP_MPIRUN replaces the
mpirun word when the MPI launcher has another name.
"""

from __future__ import annotations

import os
import shlex
import shutil
from pathlib import Path

MPI_LAUNCHERS = ("mpirun", "mpiexec", "srun", "orterun", "prun")


def ranks(text: str | None = None) -> int:
    text = os.environ.get("LMP_NP", "") if text is None else text
    text = text.strip()
    if not text:
        return 1
    if text == "all":
        return os.cpu_count() or 1
    value = int(text)
    if value < 1:
        raise SystemExit(f"LMP_NP must be at least 1, not {text}")
    return value


def launcher(explicit: str = "") -> list[str]:
    """Words that start LAMMPS. Input-script arguments follow them."""
    text = explicit or os.environ.get("LMP", "")
    if not text:
        raise SystemExit("LMP is unset")
    words = shlex.split(text)
    if not words:
        raise SystemExit("LMP is empty")
    head = Path(words[0]).name
    if head not in MPI_LAUNCHERS:
        n = ranks()
        if n > 1:
            words = [os.environ.get("LMP_MPIRUN", "mpirun"), "-np", str(n), *words]
    first = words[0]
    if not Path(first).is_file() and shutil.which(first) is None:
        raise SystemExit(f"LMP not found: {first}")
    return words


def describe(explicit: str = "") -> str:
    """Launcher as recorded in CSV headers. Paths are reduced to names."""
    words = launcher(explicit)
    return " ".join(Path(w).name if "/" in w else w for w in words)


def self_test() -> None:
    saved = {key: os.environ.get(key) for key in ("LMP", "LMP_NP", "LMP_MPIRUN")}
    try:
        os.environ["LMP"] = "sh"
        os.environ.pop("LMP_NP", None)
        if launcher() != ["sh"]:
            raise SystemExit(f"plain launcher {launcher()}")
        os.environ["LMP_NP"] = "4"
        os.environ["LMP_MPIRUN"] = "env"
        if launcher() != ["env", "-np", "4", "sh"]:
            raise SystemExit(f"LMP_NP launcher {launcher()}")
        os.environ["LMP"] = "mpiexec -n 2 lmp"
        words = shlex.split(os.environ["LMP"])
        if Path(words[0]).name not in MPI_LAUNCHERS:
            raise SystemExit("mpiexec was not seen as a launcher")
        if ranks("all") != (os.cpu_count() or 1):
            raise SystemExit("LMP_NP=all is not the processor count")
        os.environ.pop("LMP")
        try:
            launcher()
        except SystemExit as exc:
            if "unset" not in str(exc):
                raise
        else:
            raise SystemExit("an unset LMP was accepted")
    finally:
        for key, value in saved.items():
            if value is None:
                os.environ.pop(key, None)
            else:
                os.environ[key] = value
    print("self-test ok")


if __name__ == "__main__":
    self_test()
