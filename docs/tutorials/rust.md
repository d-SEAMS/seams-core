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

## The rest of the call set

`bindings/rust/seams/tests/cpp_match.rs` is the test `meson test docs_api_rust` runs.
It calls `chill`, `steinhardt` with the Lechner-Dellago average, `rings` and `vertices`, `cages`, `knn`, `rdf`, `running_cn`, `coordination`, `hbonds`, `density_axis`, `fingerprint`, `site_table`, `contact_pairs`, `domains`, `prisms`, `monolayer`, and `signature_cages`, then checks each value against `tests/binding_oracle.cpp`.
On ice Ih and the sI guest frame the test prints:

```
rustcall hexagonal 128
rustcall q6bar 0.562448843702
rustcall vert 6 0
rustcall cages 128 0 0 0
rustcall knn 512
rustcall sheet 256
rustcall sig 2 40
```


## Structure still

```{figure} ../tutorials/molecules/rust/overlay-web.png
:alt: Molecule still for this front
:width: 480px

Ice Ih coloured by published CHILL+, drawn with chemparseplot. `python docs/scripts/render_molecules.py --slug rust` writes the 2400 by 1800 PNG.
```

