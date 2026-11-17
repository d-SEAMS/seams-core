# Fortran

The fpm package `dseams` is [d-SEAMS/dseams-fortran](https://github.com/d-SEAMS/dseams-fortran).

| module | role |
| --- | --- |
| `seams` | the module programs call |
| `seams_c` | `iso_c_binding` image of `include/seams.h` |

The package links an installed `libseams` and does not build the engine. Kinds are `seams_real` (`real64`) and `seams_int` (`int32`). Coordinates have shape `(3, n)`.

The README's fpm dependency names the git tag `v2.14.0`. That tag is not on the repository, so the dependency line does not resolve.

The first classification is the {doc}`Fortran tutorial </fronts/fortran>`.
