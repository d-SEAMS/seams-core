# Reproducibility deposit (Zenodo)

The deposit is the measurement campaign, not the manuscript. Do not put
`rg_main.org`, the CPC PDF, the cover letter, or the highlights in this
payload. Pin the software to a release tag, not a commit.

## Payload

| Path | What |
|---|---|
| `paper_manifest.json` | parsed timings, identity gate, conditions |
| `conditions.txt` | host, job, SHAs, CPU, load |
| `tip-*.txt`, `base-*.txt` | raw bench stdout |
| `tip-incremental-*.txt` | hop-bound updater vs rebuild at each *n* |
| `tip-pipeline.txt` | rings, k-NN, affiliation vs *n* |
| `tip-stages-*.txt` | per-stage times on cubic and nucleation frames |
| `config.yaml` | sizes, reps, baseline SHA |
| `Snakefile` | the DAG that produced the files |
| `figshare-demos/` | the five v1 deposits through `require("dseams")` |
| `figshare-incremental.json` | per-frame incremental rings and seeded labels |

Fill `payload/` from `repro/results/` after `repro/elja_submit.sh run`
exits 0. The driver is `repro/stage_zenodo.sh`. Pin the software
in the record to a release tag (`v2.6.0`), not a commit.

## What this is not

A CPC-library software tarball. That is the tagged release of
seams-core / PydSEAMSlib / yodaStruct. This record is the exclusive-node
campaign that the paper's timing figures read.
