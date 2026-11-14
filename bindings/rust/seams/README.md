# dseams

Safe Rust calls into [d-SEAMS](https://dseams.info): CHILL and CHILL+
ice labels, Steinhardt order parameters, primitive rings, cages, RDFs,
hydrogen bonds and topological fingerprints for molecular simulation
frames.

The crates.io name `seams` belongs to an unrelated project, so the
package is `dseams`. The library is still called `seams`:

```toml
[dependencies]
dseams = "2.14.0"
```

```rust
fn main() -> Result<(), seams::Error> {
    let a = 2.75 * 4.0 / 3.0_f64.sqrt();
    let sites = [
        (0.0, 0.0, 0.0), (0.5, 0.5, 0.0), (0.5, 0.0, 0.5), (0.0, 0.5, 0.5),
        (0.25, 0.25, 0.25), (0.75, 0.75, 0.25), (0.75, 0.25, 0.75), (0.25, 0.75, 0.75),
    ];
    let mut xyz = Vec::new();
    for (x, y, z) in sites {
        xyz.extend([x * a, y * a, z * a]);
    }
    let mut sys = seams::System::new()?;
    sys.set_ortho(&xyz, a, a, a)?;
    let (labels, _counts) = sys.chill_plus(
        3.5,
        seams::ffi::SEAMS_GRAPH_CUTOFF,
        seams::ffi::SEAMS_CHILL_PUBLISHED,
    )?;
    assert_eq!(labels.len(), 8);
    Ok(())
}
```

By default `dseams-sys` compiles a bundled libseams (x86_64, C++20 and
OpenMP). Enable the `system` feature to link an installed libseams
through pkg-config instead; see the `dseams-sys` README.

Licensed under MIT.
