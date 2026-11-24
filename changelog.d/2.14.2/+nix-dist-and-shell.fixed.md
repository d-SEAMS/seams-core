The Nix package source includes `scripts/dist_strip_vcs.py`, which Meson reads when it finds Python. The development shell omits the workflow runner that depends on the broken `stopit` package.
