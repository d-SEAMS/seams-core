# Rust

The Rust binding is the crates `dseams-sys` and `dseams`. It lives in its own repository,
[d-SEAMS/dseams-rs](https://github.com/d-SEAMS/dseams-rs), with its API reference, tutorials,
tests and release notes. It links the engine through the C API
(`libseams`, `include/seams.h`) and the C ABI `libseams_lang`
(`bindings/cxx/seams_lang.h`) that this repository builds and installs.

```bash
cargo add dseams
```

The binding release with the same version as seams-core is built
against that seams-core tag. `docs/orgmode/contributing/repositories.org`
describes how the repositories fit together.
