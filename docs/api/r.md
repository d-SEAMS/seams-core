# R

The R binding is the R package `dseams`. It lives in its own repository,
[d-SEAMS/dseams-r](https://github.com/d-SEAMS/dseams-r), with its API reference, tutorials,
tests and release notes. It links the engine through the C API
(`libseams`, `include/seams.h`) and the C ABI `libseams_lang`
(`bindings/cxx/seams_lang.h`) that this repository builds and installs.

The binding documentation is at <https://r.dseams.info/>, and the R walkthrough is [on the video page](https://docs.dseams.info/videos/embeds/r.html).

Install it from a source checkout or a GitHub release source tarball using
that repository's README.

The package is not on CRAN or r-universe yet; those registry routes are
coming. `docs/orgmode/contributing/repositories.org` describes how the repositories fit together.
