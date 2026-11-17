# Julia

The Julia binding is the Julia package `Dseams`. It lives in its own repository,
[d-SEAMS/Dseams.jl](https://github.com/d-SEAMS/Dseams.jl), with its API reference, tutorials,
tests and release notes. It links the engine through the C API
(`libseams`, `include/seams.h`) and the C ABI `libseams_lang`
(`bindings/cxx/seams_lang.h`) that this repository builds and installs.

The [README](https://github.com/d-SEAMS/Dseams.jl#user-content-a-local-engine-build) of that repository
gives the install command and a first example.

The binding release with the same version as seams-core is built
against that seams-core tag. `docs/orgmode/development/repositories.org`
describes how the repositories fit together.
