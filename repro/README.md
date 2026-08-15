# Reproducing the paper numbers

One Snakemake DAG rebuilds both trees, gates on a green test suite, runs
every bench behind the paper's tables, and writes a single manifest. On a
Slurm cluster the heavy rules execute through HyperQueue inside one
exclusive allocation.

## On a Slurm cluster

From a login node with network access, in the repository root:

```sh
repro/elja_submit.sh prep
ELJA_ACCOUNT=<account> repro/elja_submit.sh submit
```

`prep` solves the `repro` pixi environment, downloads the meson subproject
wraps for both trees, creates the baseline worktree at the pinned commit,
records the tip SHA, and fetches the HyperQueue binary. `submit` queues one
exclusive job whose body starts an `hq` server plus a worker and runs the
DAG. The result is `repro/results/paper_manifest.json` next to the raw
per-bench outputs and the conditions record (node, SHAs, load, CPU).

## Anywhere else

Without Slurm or `hq` the same DAG runs the commands directly:

```sh
repro/elja_submit.sh prep       # or perform its steps by hand
SEAMS_NO_HQ=1 pixi run -e repro -- snakemake -s repro/Snakefile --cores all
```

## The v1 figshare demonstrations

The 1.0 paper demonstrated on five LAMMPS trajectories archived in the
figshare project [d-SEAMS
Datasets](https://figshare.com/projects/d-SEAMS_Datasets/73545). The DAG
fetches those exact deposits (MD5-verified against figshare's records,
`repro/scripts/figshare_demos.py`) and demonstrates on them twice over:

- The published Lua workflows run unmodified with the current binary
  (CHILL+ on the Ic lattice, the bulk topological criterion on the
  crystallized end of the nucleation run, the quasi-one-dimensional
  nanotube, the monolayer, and the in-plane 2D RDF).
- The percent-format notebooks under `repro/notebooks/` rerun the same
  analyses through the `pydseamslib` Python bindings. Jupytext converts
  each source to ipynb and papermill executes it; execution is the
  test (each notebook asserts its own headline numbers) and the
  executed `.ipynb` files are the artifacts. The nucleation notebook
  classifies every frame with the incremental ring and affiliation
  updaters, referees the incremental rings against a batch
  recomputation on each frame, runs the seeded classifier, and writes
  `figshare-incremental.json`.

On clusters whose compute nodes are offline the download happens
during `prep`; anywhere else the `figshare_fetch` rule downloads on
first use into `repro/figshare/`.

## What the DAG produces

| Output | Feeds |
|---|---|
| `tip-scaling.txt`, `base-scaling.txt` | the neighbour/CHILL+ scaling table |
| `tip-cages.txt`, `base-cages.txt` | the ring-and-cage pipeline paragraph |
| `tip-overhead.txt` | the vesin-overhead note |
| `tip-strong-t{1,2,4}.txt` | the thread sweep |
| `figshare-demos/figshare-demos.json` | the five v1 Lua workflows on the figshare deposits |
| `notebooks/*.ipynb` | the executed Python-bindings notebooks on the same deposits |
| `figshare-incremental.json` | per-frame incremental + seeded run over the nucleation deposit |
| `trajectory-incremental.txt` | incremental rings and cage affiliation, with the exactness referee (nonzero exit on any inequality) |
| `ql-dseams.txt`, `ql-python.txt` | the freud/pyscal3 descriptor comparison |
| `tip-test.log` | the identity gate: a red suite aborts the DAG |
| `paper_manifest.json` | everything above, parsed into one object |

The baseline commit, sizes, trajectory, and repetition counts live in
`repro/config.yaml`. Timings on a shared or loaded machine are not
comparable to the paper's exclusive-node numbers; the conditions file
records the load so a reviewer can tell.
