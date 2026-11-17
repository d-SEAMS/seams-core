# Rust

Classify ice Ih from Rust. The recording uses `input/traj/ice1h.lammpstrj`: `seams read` prints `nop 128`, and published CHILL+ at a 3.5 angstrom cutoff prints `nop 128 hexagonal 128`.

```{eval-rst}
.. include:: ../videos/embeds/rust.rst
```

The binding is [d-SEAMS/dseams-rs](https://github.com/d-SEAMS/dseams-rs). `dseams-sys` is the raw FFI for `include/seams.h`. `dseams` is the safe wrapper. Neither crate is on crates.io. The crates.io name `seams` is an unrelated project. The README gives the dependency line and the `SEAMS_LIBDIR` build against this tree.

The same two engine commands, without the crates:

```bash
seams read input/traj/ice1h.lammpstrj
seams chill-plus input/traj/ice1h.lammpstrj --cutoff 3.5 --type 1
```

Crate names are the {doc}`Rust API </stable/rust>`.
