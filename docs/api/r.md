# R

d-SEAMS 2.14.0. The R package is `bindings/r/dseams`. It calls
`libseams_lang` through cpp11. The exported frames are `read_frame`,
`nframes`, `positions`, `chill_plus`, `cages`, `primitive_rings`,
`signature_cages`, `symmetry`, `steinhardt`, `q4_q6`, and `jump_tau90`.
`print` on a frame dispatches to `print.dseams_frame`.

Build the library, then install the package against it:

```bash
meson setup build -Dwith_tests=true
meson compile -C build
DSEAMS_ROOT=$PWD R CMD INSTALL bindings/r/dseams
```

CI runs the testthat file when `Rscript` is on `PATH` and
`DSEAMS_ROOT` points at this tree. The meson target is `language_fronts`
(`tests/run_language_fronts.sh`). A missing R skips that front.
`docs/examples/minimal_r.R` reads the GenIce sI
cell and calls `steinhardt`, `chill_plus`, and `primitive_rings`.
`meson test docs_examples` runs it.

```{literalinclude} ../../../docs/examples/minimal_r.R
:language: r
```
`meson test docs_api_r` installs the package against `libseams_lang` and
runs the script below.

```{literalinclude} ../../../docs/examples/r_api.R
:language: r
```

`seams_lang_fixtures` is the C++ driver of `seams_lang.h`, not the R package.

```bash
DSEAMS_ROOT=$PWD meson test -C build language_fronts --print-errorlogs
```

```{literalinclude} ../../../bindings/r/dseams/tests/testthat/test-ice.R
:language: r
:lines: 9-36,77-91
```
