# Fortran tutorial

`bindings/fortran/seams.f90` binds every C symbol with `iso_c_binding`. `test_seams.f90` is the program you compile.

```bash
mkdir -p build/fmod
gfortran -Jbuild/fmod -c bindings/fortran/seams.f90 -o build/fmod/seams.o
gfortran -Ibuild/fmod bindings/fortran/test_seams.f90 build/fmod/seams.o \
  -o build/test_seams -Lbuild/src -lseams -Wl,-rpath,build/src
build/test_seams
```

It prints `fortran ok` after the 64-atom cell is entirely cubic under both CHILL+ rules, q6 is above 0.3, and the six-ring fraction is above one half.

The program the script compiles:

```{literalinclude} ../../../bindings/fortran/test_seams.f90
:language: fortran
:lines: 37-74
```
