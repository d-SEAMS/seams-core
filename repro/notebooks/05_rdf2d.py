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
# # In-plane RDF of the monolayer deposit
#
# The d-SEAMS 1.0 paper computed the in-plane radial distribution
# function of confined monolayer water (figshare
# doi:10.6084/m9.figshare.11448711.v1), alongside the same ring
# classification as the monolayer example. This notebook reruns both via
# the Python bindings on the deposited trajectory.

# %%
import tempfile
from pathlib import Path

from pydseamslib import Trajectory

ROOT = next(
    p
    for p in [Path.cwd(), *Path.cwd().parents]
    if (p / "repro" / "figshare").is_dir()
)
TRAJ = ROOT / "repro" / "figshare" / "dump-320.lammpstrj"

traj = Trajectory(
    TRAJ,
    frame=1,
    atom_type=2,
    cutoff=3.5,
    region=([0, 0, 0], [50, 0, 0]),
)
traj.n_atoms

# %% [markdown]
# Ring statistics of the sliced monolayer, as in the published workflow.

# %%
out = tempfile.mkdtemp(prefix="dseams_rdf2d_")
counts = traj.monolayer_rings(out, sheet_area=50.0 * 50.0, max_depth=4)
counts

# %% [markdown]
# The in-plane $g(r)$: for square ice the first peak sits at the O--O
# hydrogen-bond distance (~2.7 Å) and a lattice peak follows near the
# diagonal distance.

# %%
r, g = traj.rdf_2d(out, cutoff=12.0, binwidth=0.05)
first_peak_g, first_peak_r = max(
    (gv, rv) for rv, gv in zip(r, g) if rv < 4.0
)
round(first_peak_r, 3), round(first_peak_g, 3)

# %%
import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt

fig, ax = plt.subplots(figsize=(5.5, 3.4), constrained_layout=True)
ax.plot(r, g)
ax.set_xlabel(r"r (Å)")
ax.set_ylabel("g(r)")
ax.axhline(1.0, lw=0.6, color="k")
fig

# %%
assert 2.4 < first_peak_r < 3.2
assert first_peak_g > 1.5
