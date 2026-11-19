<!-- https://docs.dseams.info/stable/rust.md -->

# Rust {#rust}

| crate | role | library |
| --- | --- | --- |
| dseams-sys | raw FFI generated from include/seams.h | libseams |
| dseams | safe handles over that FFI | libseams |

Both live in [d-SEAMS/dseams-rs](https://github.com/d-SEAMS/dseams-rs).

```bash
cargo add dseams@2.14.0
```

The safe crate calls the same C API as the CLI. On `input/traj/ice1h.lammpstrj` (128 oxygens, type 1) published CHILL+ at 3.5 angstrom is 128 hexagonal:

```default
nop 128 frame 1 box 15.6457 14.7071 18.073
nop 128 hexagonal 128
```

The first classification is the [Rust tutorial](../fronts/rust.md).
