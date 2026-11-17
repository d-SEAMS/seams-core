# Rust

| crate | role | library |
| --- | --- | --- |
| `dseams-sys` | raw FFI generated from `include/seams.h` | `libseams` |
| `dseams` | safe handles over that FFI | `libseams` |

Both live in [d-SEAMS/dseams-rs](https://github.com/d-SEAMS/dseams-rs). Neither crate is on crates.io. The crates.io name `seams` is an unrelated project. The README gives the dependency line and the `SEAMS_LIBDIR` build.

The first classification is the {doc}`Rust tutorial </fronts/rust>`.
