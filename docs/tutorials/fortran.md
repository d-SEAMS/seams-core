# Fortran tutorial

`bindings/fortran/seams.f90` binds every C symbol with `iso_c_binding`. `test_seams.f90` is the program you compile.

```bash
gfortran -J/tmp/fmod -c bindings/fortran/seams.f90 -o /tmp/fmod/seams.o
gfortran -I/tmp/fmod bindings/fortran/test_seams.f90 /tmp/fmod/seams.o \
  -o /tmp/test_seams -Lbuild/src -lseams -Wl,-rpath,build/src
/tmp/test_seams
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
