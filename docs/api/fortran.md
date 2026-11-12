# Fortran

`bindings/fortran/seams.f90` is an `iso_c_binding` module. Every procedure is
`bind(C)` to the symbol in `seams.h`.

```fortran
use seams
sys = seams_system_create()
code = seams_system_set_ortho(sys, xyz, n, 0d0, 0d0, 0d0, lx, ly, lz, &
                              c_null_ptr, c_null_ptr)
code = seams_chill_plus(sys, 3.5d0, 0, SEAMS_GRAPH_CUTOFF, &
                        SEAMS_CHILL_PUBLISHED, labels, counts)
code = seams_q4q6(sys, 3.5d0, 0, q4, q6)
code = seams_rings_find(sys, 3.5d0, 0, SEAMS_GRAPH_CUTOFF, 6, rings)
```

`seams_set_num_threads` and `seams_get_num_threads` are in the same module.
`test_seams.f90` is the example and the test: it labels a cubic cell, checks
q6, and counts six-rings.
