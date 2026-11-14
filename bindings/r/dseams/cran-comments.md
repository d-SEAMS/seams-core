## Submission

This is a new submission.

The package compiles a copy of the d-SEAMS engine (C++20, in `src/seams`)
together with the R glue. `tools/vendor.sh` makes that copy from the
seams-core repository when the source package is built; installation does
not use the network. In the copy, `[[nodiscard]]` is removed (the
upstream build silences those warnings with `-Wno-unused-result`) and the
engine's console output and `exit` go through `src/dseams_r_io.cpp` to
`Rprintf`/`REprintf` and an R error.

## R CMD check results

0 errors | 0 warnings | 1 note

* checking CRAN incoming feasibility ... NOTE
  New submission

Checked with `R CMD check --as-cran` on R 4.5.3 (x86_64 Linux, GCC 14,
`-g -O2 -Wall -pedantic`). The R-CMD-check workflow in the repository
also checks Ubuntu (R release and devel), macOS and Windows.

The examples and tests run with at most two threads when
`_R_CHECK_LIMIT_CORES_` is set.
