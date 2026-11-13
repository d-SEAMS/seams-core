# Fortran

`bindings/fortran/seams.f90` is an `iso_c_binding` module. Every procedure is
`bind(C)` to the symbol in `seams.h`.

`seams_set_num_threads` and `seams_get_num_threads` are in the same module.
`bindings/fortran/test_seams.f90` is the program the language-front script
compiles. It calls every `bind(C)` procedure in the module on a 64-atom
cubic cell. `meson test docs_api_fortran` compiles that program with
gfortran and runs it.
`docs/examples/minimal_fortran.f90` reads
`input/traj/genice_sI.lammpstrj` and calls the same C entry points.
`meson test docs_examples` runs it.

The module exports `SEAMS_OK`, `SEAMS_CHILL_PUBLISHED`, `SEAMS_CHILL_OVITO`,
`SEAMS_GRAPH_CUTOFF`, `SEAMS_GRAPH_MUTUAL`, `SEAMS_LABEL_CUBIC`,
`seams_version_string`, `seams_last_error`, `seams_set_num_threads`,
`seams_get_num_threads`, `seams_system_create`, `seams_system_destroy`,
`seams_system_set_ortho`, `seams_system_set_frame`, `seams_chill_plus`,
`seams_q4q6`, `seams_neighbours_build`, `seams_neighbours_destroy`,
`seams_neighbours_csr`, `seams_rings_find`, `seams_rings_destroy`,
`seams_rings_histogram`, and `seams_rings_fraction`.

```{literalinclude} ../../../docs/examples/minimal_fortran.f90
:language: fortran
```

```{literalinclude} ../../../bindings/fortran/test_seams.f90
:language: fortran
:lines: 37-74
```
