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

A shorter call from your own program:

```fortran
use seams
sys = seams_system_create()
code = seams_system_set_ortho(sys, xyz, n, 0d0, 0d0, 0d0, lx, ly, lz, &
                              c_null_ptr, c_null_ptr)
code = seams_chill_plus(sys, 3.5d0, 0, SEAMS_GRAPH_CUTOFF, &
                        SEAMS_CHILL_PUBLISHED, labels, counts)
```
