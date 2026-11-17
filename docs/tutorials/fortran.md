# Fortran

The Fortran binding is the fpm package `dseams`. It lives in its own repository,
[d-SEAMS/dseams-fortran](https://github.com/d-SEAMS/dseams-fortran), with its API reference, tutorials,
tests and release notes. It links the engine through the C API
(`libseams`, `include/seams.h`) and the C ABI `libseams_lang`
(`bindings/cxx/seams_lang.h`) that this repository builds and installs.

```toml
[dependencies]
dseams = { git = "https://github.com/d-SEAMS/dseams-fortran", tag = "v2.14.0" }
```

The binding release with the same version as seams-core is built
against that seams-core tag. `docs/orgmode/development/repositories.org`
describes how the repositories fit together.
