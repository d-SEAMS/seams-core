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
`meson test docs_api_fortran` compiles that program and runs it with the repository root as its argument.
The same program calls `seams_chill`, `seams_steinhardt`, `seams_cages`, `seams_signature_cages`, `seams_neighbours_knn`, `seams_rdf`, `seams_running_cn`, `seams_coordination`, `seams_hbonds`, `seams_density_axis`, `seams_fingerprint`, `seams_contact_pairs`, `seams_domains`, `seams_prisms`, `seams_monolayer`, `seams_rings_vertices`, `seams_site_table`, and `seams_system_read`.
The lines it prints are:

```
fortran call cubic 64
fortran call q6 0.628539361055
fortran call cages 0 64 0 0
fortran call sig 0 0
fortran call cn 4.063492063492
fortran call hbond 0
fortran call sheet 128
fortran call vert 6 0
fortran call read 128 46
```

`cubic 64` is the published CHILL+ count. `q6` is the mean of `seams_q4q6`.
`cages` is the seeded census. `sig 0 0` is the `512` signature on this cubic cell.
`read 128 46` is ice Ih and the sI guest frame.

The program the script compiles:

```{literalinclude} ../../../bindings/fortran/test_seams.f90
:language: fortran
:lines: 37-74
```


## Structure still

```{figure} ../tutorials/molecules/fortran/overlay-web.png
:alt: Molecule still for this front
:width: 480px

Ice Ih coloured by published CHILL+, drawn with chemparseplot. `python docs/scripts/render_molecules.py --slug fortran` writes the 2400 by 1800 PNG.
```

