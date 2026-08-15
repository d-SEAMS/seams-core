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
# # Prisms in the quasi-1D ice nanotube deposit
#
# The d-SEAMS 1.0 paper identified prismatic blocks in a square ice
# nanotube (figshare doi:10.6084/m9.figshare.11448768.v1). This notebook
# reruns the strict-criterion analysis via the Python bindings: the
# hydrogen-bonded network, the primitive rings, and the prism
# identification with per-size counts and height coverage.

# %%
import tempfile
from pathlib import Path

from pydseamslib import Trajectory

ROOT = next(
    p
    for p in [Path.cwd(), *Path.cwd().parents]
    if (p / "repro" / "figshare").is_dir()
)
TRAJ = ROOT / "repro" / "figshare" / "dump-240-square.lammpstrj"

traj = Trajectory(TRAJ, frame=1, atom_type=2, cutoff=3.5)
traj.n_atoms, traj.box

# %% [markdown]
# The nanotube is hydrogen-bonded water: rings walk the H-bond network,
# not the raw distance cutoff.

# %%
hb_degree = sum(len(r) - 1 for r in traj.hbonds) / len(traj.hbonds)
ring_sizes = {}
for r in traj.rings:
    ring_sizes[len(r)] = ring_sizes.get(len(r), 0) + 1
round(hb_degree, 3), dict(sorted(ring_sizes.items()))

# %% [markdown]
# Prism identification stacks pairs of basal rings connected by lateral
# bonds. `nPrisms.dat` records, per ring size: the total and deformed
# prism counts and the height coverage percentage.

# %%
out = tempfile.mkdtemp(prefix="dseams_nanotube_")
traj.find_prisms(output_dir=out + "/", max_depth=6)
nprisms = (Path(out) / "topoINT" / "nPrisms.dat").read_text().splitlines()
fields = nprisms[-1].split()[1:]
prisms = {}
i = 0
while i + 3 < len(fields):
    size = int(fields[i])
    prisms[size] = {
        "prisms": int(fields[i + 1]),
        "deformed": int(fields[i + 2]),
        "coverage_pct": float(fields[i + 3]),
    }
    i += 4
prisms

# %% [markdown]
# The square (tetragonal) nanotube shows up as four-membered prisms with
# nonzero height coverage; no other prism size appears.

# %%
assert prisms[4]["prisms"] > 0
assert prisms[4]["coverage_pct"] > 0
assert all(v["prisms"] == 0 for s, v in prisms.items() if s != 4)
