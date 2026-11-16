# R

The R binding is the R package `dseams`. It lives in its own repository,
[d-SEAMS/dseams-r](https://github.com/d-SEAMS/dseams-r), with its API reference, tutorials,
tests and release notes. It links the engine through the C API
(`libseams`, `include/seams.h`) and the C ABI `libseams_lang`
(`bindings/cxx/seams_lang.h`) that this repository builds and installs.

The [README](https://github.com/d-SEAMS/dseams-r#user-content-installation) of that repository
gives the install command and a first example.

The binding release with the same version as seams-core is built
against that seams-core tag. `docs/orgmode/contributing/repositories.org`
describes how the repositories fit together.
