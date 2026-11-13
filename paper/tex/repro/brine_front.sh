#!/bin/sh
# Build paper/tex/figures/brine_front.pdf. The PDF is gitignored.
# Restarts are read from the d-SEAMS/brine-frames dataset in config.yaml.
set -eu
snakemake -s paper/data/brine-frames/Snakefile -j1 paper/tex/figures/brine_front.pdf
