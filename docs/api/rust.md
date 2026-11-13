# Rust

`seams-sys` is the raw FFI crate. `seams` is the safe wrapper. Build both with
`SEAMS_LIB_DIR` pointing at the directory that contains `libseams.so`.

```rust
let mut sys = seams::System::new()?;
sys.set_ortho(&xyz, lx, ly, lz)?;
let (labels, counts) = sys.chill_plus(3.5, seams_sys::SEAMS_GRAPH_CUTOFF, seams_sys::SEAMS_CHILL_PUBLISHED)?;
let (q4, q6) = sys.q4q6(3.5)?;
let degree = sys.neighbour_degree(3.5)?;
let (six, fraction) = sys.six_ring_fraction(3.5)?;
```

`set_frame` takes the nine-number cell, including a non-zero tilt.
`seams::set_num_threads` and `seams::version` match the C API.
`bindings/rust/seams/tests/ice.rs` is the `cargo test` that calls every
safe wrapper. `meson test docs_api_rust` runs that test.
`examples/ice.rs` prints the cubic count.

```{literalinclude} ../../../bindings/rust/seams/tests/ice.rs
:language: rust
:lines: 25-57
```
