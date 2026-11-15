# Fortran

`seams` is the module programs use. `seams_c` is the `iso_c_binding`
image of `include/seams.h`. Kinds are `seams_real` (`real64`) and
`seams_int` (`int32`).

A handle is a derived type. Assignment shares it. The last finalizer,
or `close`, destroys the C object. Pass `stat` and `errmsg` to catch a
failure. Omit both and a failure stops the program with the library
message.

Coordinates have shape `(3, n)`.

```bash
cd bindings/fortran
fpm test --link-flag "-L$SEAMS_LIB_DIR -Wl,-rpath,$SEAMS_LIB_DIR"
ford ford.md
```

`fpm test` runs `test_seams.f90` and `test/test_high.f90`. Meson runs
the same programs as `fortran_capi` and `fortran_high` when a Fortran
compiler is found. The native reference is the FORD build. The
narrative page is `docs/api/fortran.md`.
