# Fortran tutorial

`bindings/fortran/src/seams.f90` is the Fortran module. `src/seams_c.f90` is the `iso_c_binding` image of `seams.h`. `test_seams.f90` calls that image. `test/test_high.f90` calls the derived-type API.

```bash
mkdir -p build/fmod
gfortran -Jbuild/fmod -c bindings/fortran/src/seams_kinds.f90 -o build/fmod/seams_kinds.o
gfortran -Jbuild/fmod -Ibuild/fmod -c bindings/fortran/src/seams_c.f90 -o build/fmod/seams_c.o
gfortran -Jbuild/fmod -Ibuild/fmod -c bindings/fortran/src/seams_support.f90 -o build/fmod/seams_support.o
gfortran -Jbuild/fmod -Ibuild/fmod -c bindings/fortran/src/seams.f90 -o build/fmod/seams.o
gfortran -Jbuild/fmod -Ibuild/fmod -c bindings/fortran/src/seams_methods.f90 -o build/fmod/seams_methods.o
gfortran -Ibuild/fmod bindings/fortran/test_seams.f90 build/fmod/seams_c.o \
  -o build/test_seams -Lbuild/src -lseams -Wl,-rpath,build/src
build/test_seams "$PWD"
```

It prints `fortran ok` after the 64-atom cell is entirely cubic under both CHILL+ rules, q6 is above 0.3, and the six-ring fraction is above one half.

The program the script compiles:

```{literalinclude} ../../../bindings/fortran/test_seams.f90
:language: fortran
:lines: 37-74
```
