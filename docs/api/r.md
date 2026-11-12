# R

d-SEAMS 2.14.0. The R package is `bindings/r/dseams`. It calls
`libseams_lang` through cpp11. `chill_plus`, `rings`, `cages`,
`signature_cages`, `symmetry`, `steinhardt`, and `q4_q6` are the
exported frames.

Build the library, then install the package against it:

```bash
meson setup build -Dwith_tests=true
meson compile -C build
DSEAMS_ROOT=$PWD R CMD INSTALL bindings/r/dseams
```

CI runs the testthat suite when `Rscript` is on `PATH` and
`DSEAMS_ROOT` points at this tree. The meson target is
`seams_lang_fixtures`, which loads the R tests. A missing R skips
that part of the suite.

```bash
DSEAMS_ROOT=$PWD meson test -C build seams_lang_fixtures --print-errorlogs
```
