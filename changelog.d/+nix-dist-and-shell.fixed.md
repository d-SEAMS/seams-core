The Nix package source includes `scripts/dist_strip_vcs.py`, which Meson reads when it finds Python. The development shell does not include Snakemake, which depends on the broken `stopit` package.
