# Fortran

`bindings/fortran/seams.f90` is an `iso_c_binding` module. Every procedure is
`bind(C)` to the symbol in `seams.h`.

`seams_set_num_threads` and `seams_get_num_threads` are in the same module.
`bindings/fortran/test_seams.f90` is the program the language-front script
compiles. It calls every `bind(C)` procedure in the module on a 64-atom
cubic cell.

```{literalinclude} ../../../bindings/fortran/test_seams.f90
:language: fortran
:lines: 37-74
```
