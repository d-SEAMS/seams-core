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
# # Cage growth across the nucleation deposit
#
# The d-SEAMS 1.0 paper analysed a 4096-molecule mW nucleation
# trajectory (figshare doi:10.6084/m9.figshare.11448702.v1) with the
# bulk topological criterion at a single crystallized frame. The 2.0
# pipeline classifies every frame: the incremental ring updater and the
# incremental cage-affiliation updater carry answers between frames, and
# the seeded (hysteresis) classifier labels atoms from the mutual and
# union four-nearest graphs.
#
# Every frame is also refereed: the incremental ring set is compared
# against a batch recomputation, so this notebook is simultaneously the
# exactness demonstration on real data.

# %%
import json
import time
from pathlib import Path

from pydseamslib import Trajectory, _core

ROOT = next(
    p
    for p in [Path.cwd(), *Path.cwd().parents]
    if (p / "repro" / "figshare").is_dir()
)
TRAJ = ROOT / "repro" / "figshare" / "nucleation.lammpstrj"
RESULTS = ROOT / "repro" / "results"
RESULTS.mkdir(parents=True, exist_ok=True)

n_frames = sum(
    1 for line in open(TRAJ, "rb") if line.startswith(b"ITEM: TIMESTEP")
)
n_frames

# %% [markdown]
# One `Trajectory` object walks the deposit; `load_frame` keeps the
# incremental updaters alive between frames.

# %%
traj = Trajectory(TRAJ, frame=1, atom_type=1, cutoff=3.5, bonded="cutoff")
rows = []
t_batch = t_incr = 0.0
for f in range(1, n_frames + 1):
    if f > 1:
        traj.load_frame(f)

    t0 = time.perf_counter()
    rings = traj.rings
    t_incr += time.perf_counter() - t0

    t0 = time.perf_counter()
    batch = _core.ringNetwork(traj.bonds_by_index, 6)
    t_batch += time.perf_counter() - t0
    assert sorted(map(sorted, rings)) == sorted(map(sorted, batch)), f

    affil = traj.cage_affiliation()
    seeded = traj.seeded_affiliation()
    rows.append(
        {
            "frame": f,
            "natoms": traj.n_atoms,
            "six_rings": len(affil["six_rings"]),
            "ring_sources_recomputed": traj.rings_recomputed_sources,
            "affiliation_reclassified": affil["reclassified"],
            "hc_rings": sum(affil["hc"]),
            "ddc_rings": sum(affil["ddc"]),
            "seeded_hc_atoms": sum(seeded["hc"]),
            "seeded_ddc_atoms": sum(seeded["ddc"]),
        }
    )
rows[0], rows[-1]

# %% [markdown]
# The crystallizing run grows cage structure monotonically on the large
# scale: DDC- and HC-affiliated rings both climb, and the seeded
# per-atom labels follow. In a diffusive liquid nearly every bond moves
# between stored frames, so the incremental updater recomputes almost
# every source here; its payoff is quasi-static trajectories, and the
# refereed equality above is the point being demonstrated on real data.

# %%
import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt

fig, ax = plt.subplots(1, 2, figsize=(10, 3.6), constrained_layout=True)
frames = [r["frame"] for r in rows]
ax[0].plot(frames, [r["ddc_rings"] for r in rows], label="DDC rings")
ax[0].plot(frames, [r["hc_rings"] for r in rows], label="HC rings")
ax[0].set_xlabel("frame")
ax[0].set_ylabel("affiliated six-rings")
ax[0].legend(frameon=False)
ax[1].plot(frames, [r["seeded_ddc_atoms"] for r in rows], label="DDC atoms")
ax[1].plot(frames, [r["seeded_hc_atoms"] for r in rows], label="HC atoms")
ax[1].set_xlabel("frame")
ax[1].set_ylabel("seeded cage atoms")
ax[1].legend(frameon=False)
fig.savefig(RESULTS / "figshare-nucleation-growth.png", dpi=150)
fig

# %% [markdown]
# The per-frame table feeds the reproducibility manifest.

# %%
manifest = {
    "trajectory": TRAJ.name,
    "doi": "10.6084/m9.figshare.11448702.v1",
    "frames": n_frames,
    "batch_ring_seconds": round(t_batch, 3),
    "incremental_ring_seconds": round(t_incr, 3),
    "per_frame": rows,
}
(RESULTS / "figshare-incremental.json").write_text(
    json.dumps(manifest, indent=2, sort_keys=True) + "\n"
)
{k: manifest[k] for k in ("frames", "batch_ring_seconds", "incremental_ring_seconds")}

# %%
assert rows[-1]["ddc_rings"] > 4 * rows[0]["ddc_rings"]
assert all(r["seeded_ddc_atoms"] <= r["natoms"] for r in rows)
