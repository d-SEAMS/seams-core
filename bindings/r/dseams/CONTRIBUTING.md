# Contributing to dseams

The R package lives in `bindings/r/dseams` of
[seams-core](https://github.com/HaoZeke/seams-core). The repository
[contributing guide](https://github.com/d-SEAMS/seams-core/blob/main/CONTRIBUTING.md)
and the [code of conduct](CODE_OF_CONDUCT.md) apply here.

## Checks

`src/seams` and `inst/extdata` are produced by `tools/vendor.sh` from the
seams-core checkout. They are not stored in git. `configure` runs that
script when they are missing.

```sh
sh tools/vendor.sh
R CMD build .
R CMD check --as-cran dseams_*.tar.gz
```

Use at most two compile jobs (`MAKEFLAGS=-j2`). Tests and examples use at
most two OpenMP threads when `_R_CHECK_LIMIT_CORES_` is set.

Coverage of the R sources:

```r
covr::package_coverage()
```

The C++ engine is covered by the seams-core tests. `covr` reports the R
functions and the small amount of glue that R executes.
