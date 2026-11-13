#!/bin/sh
# Build paper/tex/figures/brine_front.pdf. The PDF is gitignored.
# Frames are read from d-SEAMS/brine-front-230K at hf_revision in config.yaml.
set -eu
snakemake -s paper/data/brine-frames/Snakefile -j1 paper/tex/figures/brine_front.pdf
