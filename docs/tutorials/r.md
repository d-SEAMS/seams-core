# R tutorial

d-SEAMS 2.14.0.

```bash
meson setup build -Dwith_tests=true
meson compile -C build
export DSEAMS_ROOT=$PWD
export LD_LIBRARY_PATH=$PWD/build/bindings:$PWD/build/src:$LD_LIBRARY_PATH
R CMD INSTALL bindings/r/dseams
Rscript -e 'library(dseams); print(packageVersion("dseams"))'
```

`dseams` is not on CRAN yet. `R CMD INSTALL` compiles the package's own
copy of the engine, which `tools/vendor.sh` copies out of the checkout
`DSEAMS_ROOT` names; `SEAMS_LIBDIR` and `SEAMS_INCLUDE` link the
`build/bindings` library instead.

The vignette runs on the example dumps in `inst/extdata`, which ship with
the package. With `DSEAMS_ROOT` set, the testthat suite also checks the
cubic dump in `input/traj`. CI runs the same library through:

```bash
DSEAMS_ROOT=$PWD meson test -C build language_fronts --print-errorlogs
```


## Structure still

```{figure} ../tutorials/molecules/r/overlay-web.png
:alt: Molecule still for this front
:width: 480px

Ice Ih colored by published CHILL+, drawn with chemparseplot. `python docs/scripts/render_molecules.py --slug r` writes the 2400 by 1800 PNG.
```

