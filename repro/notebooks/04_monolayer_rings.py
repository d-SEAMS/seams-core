# ---
# jupyter:
#   jupytext:
#     formats: py:percent
#     text_representation:
#       extension: .py
#       format_name: percent
#   kernelspec:
#     display_name: Python 3
#     language: python
#     name: python3
# ---

# %% [markdown]
# # Ring statistics in the confined monolayer deposit
#
# The d-SEAMS 1.0 paper classified rings in a quasi-2D confined water
# monolayer (figshare doi:10.6084/m9.figshare.11448741.v1). The deposit
# holds a 9000-atom system; the published analysis slices the confined
# region $x \in [0, 50]$, which selects the ~314 monolayer molecules.
# This notebook reruns that analysis via the Python bindings.

# %%
import tempfile
from pathlib import Path

from pydseamslib import Trajectory

ROOT = next(
    p
    for p in [Path.cwd(), *Path.cwd().parents]
    if (p / "repro" / "figshare").is_dir()
)
TRAJ = ROOT / "repro" / "figshare" / "dump-6-320-310.lammpstrj"

traj = Trajectory(
    TRAJ,
    frame=1,
    atom_type=2,
    cutoff=3.5,
    region=([0, 0, 0], [50, 0, 0]),
)
traj.n_atoms

# %% [markdown]
# Without the slice the same frame holds 3000 oxygens; the region
# restricts the analysis to the confined monolayer.

# %%
assert 250 < traj.n_atoms < 400

# %% [markdown]
# Ring classification on the hydrogen-bonded network, with coverage
# areas relative to the confining sheet ($50 \times 50$ Å$^2$).

# %%
out = tempfile.mkdtemp(prefix="dseams_monolayer_")
counts = traj.monolayer_rings(out, sheet_area=50.0 * 50.0, max_depth=4)
counts

# %% [markdown]
# The monolayer is dominated by four-membered rings, the square-ice
# motif of confined water; triangles are the defects.

# %%
assert counts[4]["count"] > counts[3]["count"]
assert counts[4]["coverage_xy"] > counts[3]["coverage_xy"]
