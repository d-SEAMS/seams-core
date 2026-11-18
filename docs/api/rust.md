# Rust

| crate | role | library |
| --- | --- | --- |
| `dseams-sys` | raw FFI generated from `include/seams.h` | `libseams` |
| `dseams` | safe handles over that FFI | `libseams` |

Both live in [d-SEAMS/dseams-rs](https://github.com/d-SEAMS/dseams-rs).

```bash
cargo add dseams@2.14.0
```

The first classification is the {doc}`Rust tutorial </fronts/rust>`.
