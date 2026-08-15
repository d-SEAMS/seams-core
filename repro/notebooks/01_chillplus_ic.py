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
# # CHILL+ on the ice Ic deposit
#
# The d-SEAMS 1.0 paper (doi:10.1021/acs.jcim.0c00031) demonstrated the
# CHILL+ classifier on a 4096-molecule mW cubic-ice trajectory, archived
# on figshare as doi:10.6084/m9.figshare.11448720.v1. This notebook
# reruns that analysis through the 2.0 Python bindings: every step the
# published Lua workflow performs, on the exact deposited file.
#
# mW is a single-site model, so the bonded graph is the cutoff neighbour
# list rather than a hydrogen-bonded network.

# %%
from pathlib import Path

from pydseamslib import Trajectory

ROOT = next(
    p
    for p in [Path.cwd(), *Path.cwd().parents]
    if (p / "repro" / "figshare").is_dir()
)
TRAJ = ROOT / "repro" / "figshare" / "mW_cubic.lammpstrj"

frame = Trajectory(TRAJ, frame=1, atom_type=1, cutoff=3.5, bonded="cutoff")
frame.n_atoms, frame.box

# %% [markdown]
# CHILL+ classifies each molecule from the number of staggered and
# eclipsed bonds among its four nearest neighbours. On a cubic-ice
# lattice every molecule should come out cubic.

# %%
counts = frame.classify_chill_plus()
counts

# %%
assert counts.get("cubic", 0) == frame.n_atoms, counts

# %% [markdown]
# The neighbour-averaged Steinhardt parameter separates the same frame
# from a liquid: for cubic ice $\bar q_6$ sits near 0.51 for every
# molecule (a liquid scatters far below).

# %%
q6 = frame.steinhardt(order_l=6)
mean_q6bar = sum(q6["ql_bar"]) / len(q6["ql_bar"])
round(mean_q6bar, 4)

# %%
assert mean_q6bar > 0.45

# %% [markdown]
# The topological view of the same frame: on the cutoff bonded graph
# every six-ring is DDC-affiliated and none HC-affiliated, which is the
# cage statement of "100% ice Ic".

# %%
affil = frame.cage_affiliation()
n_six = len(affil["six_rings"])
n_ddc = sum(affil["ddc"])
n_hc = sum(affil["hc"])
n_six, n_ddc, n_hc

# %%
assert n_ddc == n_six and n_hc == 0
