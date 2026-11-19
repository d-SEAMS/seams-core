# Fortran

The Fortran binding is the fpm Git dependency `dseams` in [d-SEAMS/dseams-fortran](https://github.com/d-SEAMS/dseams-fortran).

| module | role |
| --- | --- |
| `seams` | the module programs call |
| `seams_c` | `iso_c_binding` image of `include/seams.h` |

The package links an installed `libseams`. Kinds are `seams_real` (`real64`) and `seams_int` (`int32`). Coordinates have shape `(3, n)`.

```toml
[dependencies]
dseams = { git = "https://github.com/d-SEAMS/dseams-fortran" }
```

The first classification is the {doc}`Fortran tutorial </fronts/fortran>`.
