# Rust

The Rust binding is the crates `dseams-sys` and `dseams`. It lives in its own repository,
[d-SEAMS/dseams-rs](https://github.com/d-SEAMS/dseams-rs), with its API reference, tutorials,
tests and release notes. It links the engine through the C API
(`libseams`, `include/seams.h`) and the C ABI `libseams_lang`
(`bindings/cxx/seams_lang.h`) that this repository builds and installs.

The binding documentation is at <https://rust.dseams.info/>, and the Rust walkthrough is [on the video page](https://docs.dseams.info/videos/embeds/rust.html).

```bash
cargo add dseams --git https://github.com/d-SEAMS/dseams-rs
```

This uses the repository directly; a crates.io release is a separate
registry step. `docs/orgmode/contributing/repositories.org` describes how the repositories fit together.
