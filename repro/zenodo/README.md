# Reproducibility deposit (Zenodo)

The deposit is the measurement campaign, not the manuscript. Do not put
`rg_main.org`, the CPC PDF, the cover letter, or the highlights in this
payload. Pin the software to a release tag, not a commit.

## Payload

| Path | What |
|---|---|
| `paper_manifest.json` | parsed timings, identity gates, source graph, parity evidence, conditions |
| `source-manifest.json` | exact revisions used by the campaign |
| `workflow-parity.json` | executable CLI/Python/Lua agreement report |
| `conditions.txt` | host, job, SHAs, CPU, load |
| `tip-*.txt`, `base-*.txt` | raw bench stdout |
| `tip-incremental-*.txt` | hop-bound updater vs rebuild at each *n* |
| `tip-pipeline.txt` | rings, k-NN, affiliation vs *n* |
| `tip-stages-*.txt` | per-stage times on cubic and nucleation frames |
| `config.yaml` | sizes, reps, baseline SHA |
| `ecosystem-lock.json` | immutable frontend and linkcell source inputs |
| `Snakefile` | the DAG that produced the files |
| `figshare-demos/` | the five v1 deposits through `require("dseams")` |
| `figshare-incremental.json` | per-frame incremental rings and seeded labels |

Fill `payload/` from `repro/results/` after `repro/elja_submit.sh run`
exits 0. The driver is `repro/stage_zenodo.sh`. Pin each software record
to the paper release tag, verify that each locked revision belongs to that
tag, and use the tags rather than commits as the deposit provenance.

## What this is not

A CPC-library software tarball. That is the tagged release of
seams-core / PydSEAMSlib / yodaStruct. This record is the exclusive-node
campaign that the paper's timing figures read.
