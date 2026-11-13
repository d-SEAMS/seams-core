# Rust tutorial

`seams-sys` is the FFI crate. `seams` is the safe wrapper. Point `SEAMS_LIB_DIR` at `libseams.so` before building.

`bindings/rust/seams/tests/ice.rs` is the test `cargo test` runs:

```{literalinclude} ../../../bindings/rust/seams/tests/ice.rs
:language: rust
:lines: 29-63
```

The ice example builds a 64-atom cell:

```bash
export SEAMS_LIB_DIR=$PWD/build/src
export LD_LIBRARY_PATH=$PWD/build/src
cargo run --manifest-path bindings/rust/Cargo.toml --example ice
cargo test --manifest-path bindings/rust/Cargo.toml -- --test-threads=1
```

The integration test checks both CHILL+ rules, q6 above 0.3, four neighbours, and a six-ring fraction above one half.
