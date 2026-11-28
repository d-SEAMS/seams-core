<!-- https://docs.dseams.info/stable/fortran.md -->

# Fortran {#fortran}

The Fortran binding is the fpm Git dependency `dseams` in [d-SEAMS/dseams-fortran](https://github.com/d-SEAMS/dseams-fortran).

| module | role |
| --- | --- |
| seams | the module programs call |
| seams_c | iso_c_binding image of include/seams.h |

The package links an installed `libseams`. Kinds are `seams_real` (`real64`) and `seams_int` (`int32`). Coordinates have shape `(3, n)`.

```toml
[dependencies]
dseams = { git = "https://github.com/d-SEAMS/dseams-fortran" }
```

Programs call the `seams` module, which links an installed `libseams`. Coordinates have shape `(3, n)`. On `input/traj/ice1h.lammpstrj` (128 oxygens, type 1) published CHILL+ at 3.5 angstrom is 128 hexagonal:

```default
nop 128 frame 1 box 15.6457 14.7071 18.073
nop 128 hexagonal 128
```

The first classification is the [Fortran tutorial](../fronts/fortran.md).
