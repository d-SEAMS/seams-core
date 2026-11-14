# dseams-sys

Raw FFI declarations for the [d-SEAMS](https://dseams.info) C API
(`include/seams.h`). Most code wants the safe [`dseams`](https://crates.io/crates/dseams)
crate instead.

The crates.io name `seams` belongs to an unrelated project, so the
packages are `dseams-sys` and `dseams`. The library names are unchanged:

```toml
[dependencies]
dseams-sys = "2.14.0"
```

```rust
use seams_sys::{seams_version_string, SEAMS_OK};
```

## How libseams is found

1. `SEAMS_LIB_DIR` set: link `libseams` from that directory (a meson build
   tree or an install prefix) and add it to the rpath.
2. Feature `system`: find `seams >= 2.14.0` through pkg-config
   (`PKG_CONFIG_PATH=$prefix/lib/pkgconfig`). The installed library keeps
   every reader it was configured with, including `.con` and chemfiles.
3. Otherwise build the bundled copy under `vendor/` with the `cc` crate:
   the C API, the analyses behind it, and Eigen 3.4.0. This needs a
   C++20 compiler with OpenMP (GCC 11 or newer) on x86_64; the shell
   classifier uses AVX2. The bundled build reads LAMMPS dumps and XYZ,
   not `.con` or the chemfiles formats.

In a git checkout without `vendor/`, the bundled build compiles the
sources at the repository root and takes Eigen from `EIGEN3_INCLUDE_DIR`
or pkg-config. `python3 bindings/rust/seams-sys/vendor.py` fills
`vendor/` before `cargo package`.

## License

The crate and d-SEAMS are MIT. The bundled Eigen is MPL-2.0, with
`AlignedBox.h` under BSD-3-Clause and `BFloat16.h` under Apache-2.0; the
texts are in `vendor/eigen/`.
