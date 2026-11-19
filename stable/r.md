<!-- https://docs.dseams.info/stable/r.md -->

# R {#r}

The R package `dseams` is [d-SEAMS/dseams-r](https://github.com/d-SEAMS/dseams-r). It links `libseams_lang` (`bindings/cxx/seams_lang.h`).

```text
remotes::install_github("d-SEAMS/dseams-r")
```

The package reads a frame through `libseams_lang`. On `input/traj/ice1h.lammpstrj` (128 oxygens, type 1) published CHILL+ at 3.5 angstrom is 128 hexagonal:

```default
nop 128 frame 1 box 15.6457 14.7071 18.073
nop 128 hexagonal 128
```

The first classification is the [R tutorial](../fronts/r.md).
