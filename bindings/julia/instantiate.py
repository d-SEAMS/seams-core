# Copyright (c) d-SEAMS developers
"""Install the Dseams project's recorded dependencies, including Preferences.

``using Dseams`` under ``--project=bindings/julia/Dseams`` requires
Preferences. The project names it; the depot does not have it until
``dev/develop_jll.jl`` has instantiated the environment. That script also
points libseams_jll at the in-tree stub and then removes the ``[sources]``
table it wrote into the tracked project file.
"""

import fcntl
import subprocess
from pathlib import Path

ROOT = Path(__file__).resolve().parents[2]
PROJECT = ROOT / "bindings" / "julia" / "Dseams"


def ensure(julia: str, env=None) -> None:
    lock_path = PROJECT / ".instantiate.lock"
    with lock_path.open("a+") as lock:
        fcntl.flock(lock.fileno(), fcntl.LOCK_EX)
        subprocess.run(
            [julia, f"--project={PROJECT}", str(PROJECT / "dev" / "develop_jll.jl")],
            cwd=PROJECT,
            env=env,
            check=True,
        )
