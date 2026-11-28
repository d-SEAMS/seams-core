<!-- https://docs.dseams.info/fronts/rust.md -->

# Rust {#rust}

Classify ice Ih from Rust. The recording uses `input/traj/ice1h.lammpstrj`: `seams read` prints `nop 128`, and published CHILL+ at a 3.5 angstrom cutoff prints `nop 128 hexagonal 128`.

The binding is [d-SEAMS/dseams-rs](https://github.com/d-SEAMS/dseams-rs). `dseams-sys` is the raw FFI for `include/seams.h`. `dseams` is the safe wrapper.

```bash
cargo add dseams
```

The same two engine commands, without the crates, print the ice Ih census:

```bash
seams read input/traj/ice1h.lammpstrj
seams chill-plus input/traj/ice1h.lammpstrj --cutoff 3.5 --type 1
```

```default
nop 128 frame 1 box 15.6457 14.7071 18.073
nop 128 hexagonal 128
```

Crate names are the [Rust API](../stable/rust.md).
