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

The vignette is evaluated when `DSEAMS_ROOT` is set. CI runs the same
library through:

```bash
DSEAMS_ROOT=$PWD meson test -C build language_fronts --print-errorlogs
```
